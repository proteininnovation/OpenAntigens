#!/usr/bin/env bash
# Regression-panel build and checks for OpenAntigens.
#
#   scripts/regression.sh build <reference_snapshot> [name]   # e.g. snapshots/openantigen-2026-11-01
#   scripts/regression.sh check <panel_snapshot> [baseline_panel_snapshot]
#
# build: the production warm-rebuild command on the panel's targets. It clones cache/ and data/ from
#   the reference production snapshot (APFS clone on macOS, reflink copy elsewhere) and runs
#   `build-fresh-snapshot --resume --reanalyze-reports --ortholog-jobs 6` with --target-panel. Every
#   build step and the publish gate are the production code; the gate's catalog minimums scale to the panel.
# check: (1) panel expectations and, given a baseline, every change against it; (2) the production
#   pre-publish checks, scripts/release_checks.py, on both portals; (3) a Playwright walk-through of
#   every page of both portals. Exits non-zero if any of the three fails.
#
# Environment: PYTHON (default .venv/bin/python), SNAPSHOT_ROOT (default snapshots),
#   PANEL (default data/regression/panel.tsv), PLAYWRIGHT_NODE_PATH (a node_modules dir with playwright).
set -euo pipefail
cd "$(dirname "$0")/.."
PYTHON=${PYTHON:-.venv/bin/python}
ROOT=${SNAPSHOT_ROOT:-snapshots}
PANEL=${PANEL:-data/regression/panel.tsv}

clone() {
  if cp -c -R "$1" "$2" 2>/dev/null; then return; fi
  cp -R --reflink=auto "$1" "$2"
}

case "${1:-}" in
  build)
    reference=${2:?reference production snapshot required}
    name=${3:-regression-$(date +%Y-%m-%d)-$(git rev-parse --short HEAD)}
    snapshot="$ROOT/$name"
    [ -e "$snapshot" ] && { echo "$snapshot already exists" >&2; exit 1; }
    for part in cache data; do [ -d "$reference/$part" ] || { echo "$reference/$part missing" >&2; exit 1; }; done
    mkdir -p "$snapshot"
    clone "$reference/cache" "$snapshot/cache"
    clone "$reference/data" "$snapshot/data"
    PYTHONPATH=src "$PYTHON" -m agdesign2.cli build-fresh-snapshot \
      --target-panel "$PANEL" --snapshot-root "$ROOT" --snapshot-name "$name" \
      --resume --reanalyze-reports --ortholog-jobs 6
    echo "$snapshot"
    ;;
  check)
    snapshot=${2:?panel snapshot required}
    baseline=${3:-}
    status=0
    PYTHONPATH=src "$PYTHON" -m agdesign2.cli check-regression-panel "$snapshot" --panel "$PANEL" \
      ${baseline:+--baseline "$baseline"} || status=1
    PYTHONPATH=src "$PYTHON" scripts/release_checks.py --site "$snapshot/public_site" \
      --human-json "$snapshot/outputs/surfy_batch" --mouse-json "$snapshot/mouse_data/reports" || status=1
    NODE_PATH=${PLAYWRIGHT_NODE_PATH:-${NODE_PATH:-}} node scripts/check_regression_site.cjs \
      "$snapshot/public_site" "$snapshot/browser_check" || status=1
    [ $status -eq 0 ] && echo "regression checks passed: $snapshot" || echo "regression checks FAILED: $snapshot" >&2
    exit $status
    ;;
  *)
    sed -n '2,20p' "$0"
    exit 2
    ;;
esac
