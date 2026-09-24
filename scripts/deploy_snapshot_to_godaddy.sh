#!/usr/bin/env bash
set -euo pipefail

# Build or reuse an OpenAntigens snapshot, then publish its public_site/ to the
# configured static host. Credentials and destination are supplied locally.

usage() {
  cat <<'USAGE'
Usage:
  scripts/deploy_snapshot_to_godaddy.sh [options] [-- build-fresh-snapshot args...]

Options:
  --snapshot-name NAME      Snapshot folder name. Default: openantigen-YYYY-MM-DD
  --snapshot-root DIR       Snapshot root. Default: snapshots
  --existing-snapshot DIR   Publish an existing snapshot directory and skip build
  --skip-build              Publish snapshots/<snapshot-name>/public_site without building
  --dry-run                 Print and dry-run rsync; do not build or publish
  -h, --help                Show this help

Environment:
  GODADDY_TARGET            Required rsync destination (set locally):
                            user@host:/home/user/public_html/
                            Use an SSH host alias with key-based auth. The script
                            verifies non-interactive SSH access before any build,
                            compression, or publish step.
  GODADDY_IDENTITY_FILE     Optional private-key path. When set, both SSH preflight
                            and rsync use this key with IdentitiesOnly=yes.
  AGDESIGN2_BIN             Optional CLI path. Default: run the checkout through
                            .venv/bin/python -m agdesign2.cli.
  RSYNC_DELETE              1 mirrors deletions on GoDaddy (via --delete-before,
                            quota-safe; .ftpquota + .well-known/ are always
                            excluded so the mirror can't wipe host-managed files),
                            0 keeps remote extras. Default: 1
  COMPRESS                  1 gzips public_site assets on disk + drops dead weight
                            (structures/, redundant OT JSON) before publishing to
                            fit the host disk quota, 0 skips. Idempotent. Default: 1
  PURGE_SUCURI              1 clears the Sucuri firewall cache after a successful
                            publish so visitors get the fresh build immediately,
                            0 skips. A requested purge must succeed. Default: 1
  SITE_URL                  Public URL used for post-deploy smoke checks.
                            Default: https://openantigens.org
  SUCURI_API_KEY            Sucuri Firewall API key (dashboard -> API -> API Details).
  SUCURI_API_SECRET         Sucuri Firewall API secret. Required (with the key) for
                            the purge; supply via env or SUCURI_ENV_FILE. Never commit.
  SUCURI_ENV_FILE           Optional host-local file sourced for the two Sucuri vars.
                            Default: ~/.config/openantigen/sucuri.env (chmod 600)
  JOBS                      build-fresh-snapshot --jobs value. Default: 4
  ASSET_JOBS                build-fresh-snapshot --asset-jobs value. Default: 4

Examples:
  scripts/deploy_snapshot_to_godaddy.sh --snapshot-name openantigen-2026-06-01
  scripts/deploy_snapshot_to_godaddy.sh --existing-snapshot snapshots/openantigen-2026-06-01
  GODADDY_TARGET=user@host:/home/user/public_html/ scripts/deploy_snapshot_to_godaddy.sh
USAGE
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

SNAPSHOT_ROOT="${SNAPSHOT_ROOT:-snapshots}"
SNAPSHOT_NAME="${SNAPSHOT_NAME:-openantigen-$(date +%Y-%m-%d)}"
EXISTING_SNAPSHOT=""
SKIP_BUILD=0
DRY_RUN=0
BUILD_ARGS=()

while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --snapshot-name)
      if [[ "$#" -lt 2 ]]; then
        echo "ERROR: --snapshot-name requires a value." >&2
        exit 2
      fi
      SNAPSHOT_NAME="$2"
      shift 2
      ;;
    --snapshot-root)
      if [[ "$#" -lt 2 ]]; then
        echo "ERROR: --snapshot-root requires a value." >&2
        exit 2
      fi
      SNAPSHOT_ROOT="$2"
      shift 2
      ;;
    --existing-snapshot)
      if [[ "$#" -lt 2 ]]; then
        echo "ERROR: --existing-snapshot requires a value." >&2
        exit 2
      fi
      EXISTING_SNAPSHOT="$2"
      SKIP_BUILD=1
      shift 2
      ;;
    --skip-build)
      SKIP_BUILD=1
      shift
      ;;
    --dry-run)
      DRY_RUN=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    --)
      shift
      BUILD_ARGS=("$@")
      break
      ;;
    *)
      echo "ERROR: unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

AGDESIGN2_BIN="${AGDESIGN2_BIN:-}"
GODADDY_TARGET="${GODADDY_TARGET:?Set GODADDY_TARGET=user@host:/home/user/public_html/ (keep host and user out of version control)}"
GODADDY_IDENTITY_FILE="${GODADDY_IDENTITY_FILE:-}"
RSYNC_DELETE="${RSYNC_DELETE:-1}"
COMPRESS="${COMPRESS:-1}"
PURGE_SUCURI="${PURGE_SUCURI:-1}"
SITE_URL="${SITE_URL:-https://openantigens.org}"
SUCURI_ENV_FILE="${SUCURI_ENV_FILE:-${HOME}/.config/openantigen/sucuri.env}"
JOBS="${JOBS:-4}"
ASSET_JOBS="${ASSET_JOBS:-4}"

if [[ -z "${GODADDY_TARGET}" ]]; then
  echo "ERROR: GODADDY_TARGET is empty." >&2
  exit 1
fi

if [[ "${GODADDY_TARGET}" != */ ]]; then
  GODADDY_TARGET="${GODADDY_TARGET}/"
fi

if [[ -n "${EXISTING_SNAPSHOT}" ]]; then
  SNAPSHOT_DIR="${EXISTING_SNAPSHOT}"
else
  SNAPSHOT_DIR="${SNAPSHOT_ROOT}/${SNAPSHOT_NAME}"
fi
PUBLIC_SITE="${SNAPSHOT_DIR}/public_site"

require_command() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "ERROR: $1 is required." >&2
    exit 1
  fi
}

run_agdesign2() {
  if [[ -n "${AGDESIGN2_BIN}" ]]; then
    if [[ ! -x "${AGDESIGN2_BIN}" ]]; then
      echo "ERROR: AGDESIGN2_BIN is not executable: ${AGDESIGN2_BIN}" >&2
      exit 1
    fi
    "${AGDESIGN2_BIN}" "$@"
    return
  fi
  if [[ ! -x .venv/bin/python ]]; then
    echo "ERROR: .venv/bin/python is not executable; set AGDESIGN2_BIN." >&2
    exit 1
  fi
  PYTHONPATH="${REPO_ROOT}/src${PYTHONPATH:+:${PYTHONPATH}}" \
    .venv/bin/python -m agdesign2.cli "$@"
}

load_sucuri_credentials() {
  if [[ "${PURGE_SUCURI}" != "1" ]]; then
    return 0
  fi
  if [[ ( -z "${SUCURI_API_KEY:-}" || -z "${SUCURI_API_SECRET:-}" ) && -f "${SUCURI_ENV_FILE}" ]]; then
    # shellcheck disable=SC1090
    source "${SUCURI_ENV_FILE}"
  fi
  if [[ -z "${SUCURI_API_KEY:-}" || -z "${SUCURI_API_SECRET:-}" ]]; then
    echo "ERROR: Sucuri credentials are required when PURGE_SUCURI=1." >&2
    echo "Set SUCURI_API_KEY/SUCURI_API_SECRET, write ${SUCURI_ENV_FILE} (chmod 600), or set PURGE_SUCURI=0." >&2
    exit 1
  fi
  require_command curl
}

preflight_remote() {
  if [[ "${GODADDY_TARGET}" != *:* ]]; then
    echo "ERROR: GODADDY_TARGET must be an rsync SSH target such as alias:public_html/." >&2
    exit 1
  fi
  local ssh_target="${GODADDY_TARGET%%:*}"
  local remote_path="${GODADDY_TARGET#*:}"
  if [[ -z "${ssh_target}" ]]; then
    echo "ERROR: GODADDY_TARGET has no SSH host." >&2
    exit 1
  fi
  if [[ -z "${remote_path}" ]]; then
    echo "ERROR: GODADDY_TARGET has no remote directory." >&2
    exit 1
  fi
  local quoted_remote_path
  printf -v quoted_remote_path '%q' "${remote_path}"
  echo "[openantigen-godaddy] checking SSH access and publish root on ${ssh_target}"
  if ! ssh "${SSH_ARGS[@]}" "${ssh_target}" \
    "test -d ${quoted_remote_path} && test -w ${quoted_remote_path}"; then
    echo "ERROR: non-interactive SSH preflight failed or the publish root is not writable for ${ssh_target}." >&2
    echo "Verify the host key, deployment key, and remote directory before retrying." >&2
    exit 1
  fi
}

validate_public_site() {
  local site_dir="$1"
  if [[ ! -d "${site_dir}" ]]; then
    echo "ERROR: public_site directory not found: ${site_dir}" >&2
    exit 1
  fi
  if [[ ! -s "${site_dir}/index.html" ]]; then
    echo "ERROR: public_site is missing a non-empty index.html: ${site_dir}" >&2
    exit 1
  fi
  if [[ ! -s "${site_dir}/portal_metadata.json" ]]; then
    echo "ERROR: public_site is missing portal_metadata.json: ${site_dir}" >&2
    exit 1
  fi
  local html_count
  html_count="$(find "${site_dir}" -type f -name '*.html' | wc -l | tr -d ' ')"
  if [[ "${html_count}" -lt 1 ]]; then
    echo "ERROR: public_site contains no HTML files: ${site_dir}" >&2
    exit 1
  fi
  if [[ ! -s "${site_dir}/downloads/agdesign2_portal_index.json" ]]; then
    echo "ERROR: public_site is missing the human portal index: ${site_dir}" >&2
    exit 1
  fi
  local validation_python="${REPO_ROOT}/.venv/bin/python"
  if [[ ! -x "${validation_python}" ]]; then
    validation_python="$(command -v python3)"
  fi
  PYTHONPATH="${REPO_ROOT}/src:${REPO_ROOT}" "${validation_python}" -c \
    'from pathlib import Path; from agdesign2.release_validation import validate_release_data, validate_structure_coverage; p=Path(__import__("sys").argv[1]); validate_structure_coverage(p); validate_release_data(p, require_full_catalog=True)' \
    "${site_dir}"
  echo "[openantigen-godaddy] validated public_site=${site_dir} html_files=${html_count}"
}

validate_compressed_assets() {
  local site_dir="$1"
  local legacy_asset
  legacy_asset="$(find "${site_dir}" -type f \( -name '*.js.gz' -o -name '*.css.gz' -o -name '*.svg.gz' \) -print -quit)"
  if [[ -n "${legacy_asset}" ]]; then
    echo "ERROR: compressed asset retains the obsolete .gz suffix: ${legacy_asset}" >&2
    exit 1
  fi
}

purge_sucuri_cache() {
  # Clear the Sucuri firewall cache so visitors get the fresh build immediately.
  # Report pages and assets may remain cached behind the proxy after publish.
  # Needs the Sucuri Firewall API key+secret (dashboard -> API -> API Details),
  # supplied via env (SUCURI_API_KEY / SUCURI_API_SECRET) or a host-local file
  # (SUCURI_ENV_FILE). NEVER commit these.
  if [[ "${PURGE_SUCURI}" != "1" ]]; then
    return 0
  fi
  echo "[openantigen-godaddy] purging Sucuri firewall cache"
  local resp
  if ! resp="$(
    curl --silent --show-error --fail --location --config - <<CURL_CONFIG
url = "https://waf.sucuri.net/api?v2&k=${SUCURI_API_KEY}&s=${SUCURI_API_SECRET}&a=clear_cache"
CURL_CONFIG
  )"; then
    echo "ERROR: Sucuri cache purge request failed." >&2
    return 1
  fi
  echo "[openantigen-godaddy] Sucuri response: ${resp}"
  if printf '%s' "${resp}" | grep -q '"status":1'; then
    echo "[openantigen-godaddy] Sucuri cache cleared"
  else
    echo "ERROR: Sucuri cache purge was not confirmed." >&2
    return 1
  fi
}

verify_live_site() {
  local asset
  local report_script
  local refs
  local citation_doi
  citation_doi="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["citation"]["doi"])' "${PUBLIC_SITE}/portal_metadata.json")"
  refs="$(grep -Eo '(portal\.css|portal-index-data\.js|portal-disease-index\.js|portal-index\.js)\?v=[^"[:space:]]+' "${PUBLIC_SITE}/index.html" | sort -u)"
  if [[ -z "${refs}" ]]; then
    echo "ERROR: could not resolve versioned index assets for live verification." >&2
    return 1
  fi
  echo "[openantigen-godaddy] verifying live index assets"
  while IFS= read -r asset; do
    curl --silent --show-error --fail --location --compressed \
      --output /dev/null "${SITE_URL%/}/${asset}"
  done <<< "${refs}"

  report_script="$(find "${PUBLIC_SITE}/report_scripts" -type f -name '*.js' -print | LC_ALL=C sort | sed -n '1p')"
  if [[ -z "${report_script}" ]]; then
    echo "ERROR: no report script found for live verification." >&2
    return 1
  fi
  asset="${report_script#${PUBLIC_SITE}/}"
  curl --silent --show-error --fail --location --compressed \
    --output /dev/null "${SITE_URL%/}/${asset}"
  local live_dir
  live_dir="$(mktemp -d)"
  trap 'rm -rf "${live_dir}"' RETURN
  mkdir -p "${live_dir}/downloads" "${live_dir}/mouse/downloads"
  for asset in portal_metadata.json downloads/agdesign2_portal_index.json mouse/portal_metadata.json mouse/downloads/agdesign2_portal_index.json; do
    mkdir -p "${live_dir}/$(dirname "${asset}")"
    curl --silent --show-error --fail --location --compressed \
      --output "${live_dir}/${asset}" "${SITE_URL%/}/${asset}"
    cmp "${live_dir}/${asset}" "${PUBLIC_SITE}/${asset}"
  done
  curl --silent --show-error --fail --location --head \
    "${SITE_URL%/}/open_targets_disease_associations.tsv" >/dev/null
  for asset in constructs.html mouse/constructs.html; do
    curl --silent --show-error --fail --location --compressed \
      --output "${live_dir}/page.html" "${SITE_URL%/}/${asset}"
    if ! cmp -s \
      <(sed -n '/id="after-export"/,/<\/section>/p' "${PUBLIC_SITE}/${asset}") \
      <(sed -n '/id="after-export"/,/<\/section>/p' "${live_dir}/page.html"); then
      echo "ERROR: public construct guidance is stale: ${asset}" >&2
      return 1
    fi
  done
  for asset in agent-guide.html llms.txt downloads/openantigens.bib downloads/openantigens.ris \
               mouse/agent-guide.html mouse/llms.txt mouse/downloads/openantigens.bib mouse/downloads/openantigens.ris; do
    curl --silent --show-error --fail --location --compressed \
      --output "${live_dir}/page.html" "${SITE_URL%/}/${asset}"
    if ! grep -Fq "${citation_doi}" "${live_dir}/page.html"; then
      echo "ERROR: public citation document is stale: ${asset}" >&2
      return 1
    fi
  done
  for asset in "$(find "${PUBLIC_SITE}/reports" -type f -name '*.html' | LC_ALL=C sort | sed -n '1p')" \
               "$(find "${PUBLIC_SITE}/mouse/reports" -type f -name '*.html' | LC_ALL=C sort | sed -n '1p')"; do
    asset="${asset#${PUBLIC_SITE}/}"
    curl --silent --show-error --fail --location --compressed \
      --output "${live_dir}/page.html" "${SITE_URL%/}/${asset}"
    if ! grep -Fq "${citation_doi}" "${live_dir}/page.html"; then
      echo "ERROR: public report citation is stale: ${asset}" >&2
      return 1
    fi
  done
  echo "[openantigen-godaddy] live asset verification passed"
}

require_command rsync
require_command ssh
SSH_ARGS=(-o BatchMode=yes -o ConnectTimeout=10 -o StrictHostKeyChecking=yes)
if [[ -n "${GODADDY_IDENTITY_FILE}" ]]; then
  if [[ ! -f "${GODADDY_IDENTITY_FILE}" ]]; then
    echo "ERROR: GODADDY_IDENTITY_FILE does not exist: ${GODADDY_IDENTITY_FILE}" >&2
    exit 1
  fi
  SSH_ARGS+=(-i "${GODADDY_IDENTITY_FILE}" -o IdentitiesOnly=yes)
fi
load_sucuri_credentials
preflight_remote

echo "[openantigen-godaddy] snapshot_dir=${SNAPSHOT_DIR}"
echo "[openantigen-godaddy] godaddy_target=${GODADDY_TARGET}"
echo "[openantigen-godaddy] rsync_delete=${RSYNC_DELETE}"

if [[ "${DRY_RUN}" == "1" ]]; then
  echo "[openantigen-godaddy] dry-run: skipping snapshot build"
elif [[ "${SKIP_BUILD}" == "1" ]]; then
  echo "[openantigen-godaddy] skipping snapshot build"
else
  echo "[openantigen-godaddy] building resumable fresh snapshot"
  run_agdesign2 build-fresh-snapshot \
    --snapshot-root "${SNAPSHOT_ROOT}" \
    --snapshot-name "${SNAPSHOT_NAME}" \
    --jobs "${JOBS}" \
    --asset-jobs "${ASSET_JOBS}" \
    --no-structure-images \
    --resume \
    --verbose \
    "${BUILD_ARGS[@]}"
fi

validate_public_site "${PUBLIC_SITE}"

# Shrink public_site to fit the GoDaddy disk quota: gzip the browser-fetched
# assets on disk (served transparently via the portal .htaccess) and drop dead
# weight (unreferenced structures/, the redundant Open Targets JSON). Idempotent,
# so re-deploying an already-compressed snapshot is a no-op. COMPRESS=0 to skip.
if [[ "${DRY_RUN}" != "1" && "${COMPRESS}" == "1" ]]; then
  echo "[openantigen-godaddy] compressing public_site (gzip JS/CSS/SVG; drop unreferenced structures + redundant OT JSON)"
  run_agdesign2 compress-public-site "${PUBLIC_SITE}" \
    --drop-unreferenced-structures --drop-redundant-downloads --verbose
fi
validate_compressed_assets "${PUBLIC_SITE}"

RSYNC_ARGS=(-az)
RSYNC_RSH_COMMAND="ssh"
for ssh_arg in "${SSH_ARGS[@]}"; do
  printf -v quoted_ssh_arg '%q' "${ssh_arg}"
  RSYNC_RSH_COMMAND+=" ${quoted_ssh_arg}"
done
RSYNC_ARGS+=(-e "${RSYNC_RSH_COMMAND}")
# The target is now the apex root, which holds host-managed files that are not
# part of the site and must survive a --delete mirror: GoDaddy's quota tracker
# and any ACME/cert challenge dir. We never ship these, so exclude them from
# both transfer and deletion.
RSYNC_ARGS+=(--exclude='.ftpquota' --exclude='.well-known/')
if [[ "${RSYNC_DELETE}" == "1" ]]; then
  # --delete-before (not the default delete-after): the host is at its disk
  # quota, so replaced files must be removed *before* the new ones land, or the
  # remote transiently holds both and re-trips quota.
  RSYNC_ARGS+=(--delete-before)
elif [[ "${RSYNC_DELETE}" != "0" ]]; then
  echo "ERROR: RSYNC_DELETE must be 0 or 1." >&2
  exit 1
fi
if [[ "${DRY_RUN}" == "1" ]]; then
  RSYNC_ARGS+=(--dry-run --itemize-changes --stats)
fi

echo "[openantigen-godaddy] publishing ${PUBLIC_SITE}/"
rsync "${RSYNC_ARGS[@]}" "${PUBLIC_SITE}/" "${GODADDY_TARGET}"

if [[ "${DRY_RUN}" == "1" ]]; then
  echo "[openantigen-godaddy] dry-run complete; no remote files were changed"
else
  echo "[openantigen-godaddy] published ${SNAPSHOT_DIR} to ${GODADDY_TARGET}"
  purge_sucuri_cache
  verify_live_site
fi
