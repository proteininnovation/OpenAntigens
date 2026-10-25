# API and MCP release operations

The managed Node app contains local query metadata and exact summary records.
Immutable target assets live on the existing GoDaddy static site. Both portals
use the same scientific release ID and the same publication controls.

The former self-contained full-catalogue package measured 245,426,866 bytes and
failed the upload budget. The own-static layout keeps those scientific records
and checksums while moving target assets outside the app ZIP. Full-catalogue
current-plus-previous full build passed complete source/static validation. Its
final code-only update measures 12,909,859 bytes compressed and 82,198,385 bytes
expanded, within the unchanged app budgets; its receipt records inherited
scientific validation and current checks separately. Its 19,898 verified static
target assets total 447,159,606 bytes. Local profiles passed both publication
modes and all five HTTP/MCP tools, with about 0.44-second reported app startup
and 348,897,280 bytes peak RSS in the largest-target check. Actual GoDaddy
transfer, headers, limits and proxy gates remain before launch. See the [plan](API_MCP_PLAN.md) and
[host verification summary](API_MCP_HOST_VERIFICATION.md) for the current gates.

## Build and review

Run from the canonical `proteininnovation/OpenAntigens` checkout with the
project Python environment active. Build both portals and finish static
compression and asset removal first. Build the public projection from the
selected frozen reports against the final human root containing `mouse/`:

```sh
PYTHONPATH=src python -m agdesign2.cli build-public-api \
  --human-summary /absolute/path/frozen-human/batch_summary.json \
  --mouse-summary /absolute/path/frozen-mouse/batch_summary.json \
  --public-site /absolute/path/snapshot/public_site \
  --output-dir /absolute/path/api-data \
  --software-commit SCIENTIFIC_PROJECTOR_COMMIT_SHA \
  --snapshot-id FROZEN_SNAPSHOT_ID \
  --previous-bundle /absolute/path/previous-api-data

PYTHONPATH=src python -m agdesign2.cli validate-public-api /absolute/path/api-data
```

Replace the commit and snapshot placeholders with the selected frozen values.
For the first release, omit `--previous-bundle`. The output must be new or empty.
Source reports, original downloads and prior releases remain unchanged.

Human disease evidence contains the published source-order top 100 associations
per target. `summary.disease_coverage` and disease-page `coverage` distinguish
`projected_count` from nullable `source_count`, with selection `portal_top_100`.
The original full website TSV stays outside the app. Its `source_reference`
records a checksum and labels its URL as live; that URL can advance later.
Interactive searches use gene names and accessions, rather than disease IDs.

Build the managed app from that validated projection:

```sh
PYTHONPATH=src python scripts/build_api_service_release.py \
  --bundle-dir /absolute/path/api-data \
  --public-site-dir /absolute/path/snapshot/public_site \
  --output /absolute/path/staging/openantigens-api.zip \
  --dry-run
```

A dry run validates and measures the real temporary app ZIP. It preserves the
source, static site and requested output paths. Repeat without `--dry-run` to
stage verified immutable query assets locally and write the ZIP plus
`openantigens-api.artifact.json`. Existing outputs cause an error.

The builder first verifies the complete logical source, closed file roles,
scientific fingerprints, schemas and all eleven projection export formats. It
then compiles a fresh own-static query runtime. App startup uses only the sealed
local metadata, compact search indexes, exact summary JSONL and per-record
SHA/offset tables. Search selects bounded records locally and fetches no target
assets. Other read operations fetch only the selected target's immutable packs.

The app seal is `service-build.json.data_transport`, format
`openantigens-own-static-query-v1`. It pins the query root and original source
root checksum. The query root pins the release metadata, schema, indexes,
summary records, offsets and target transport index. The code-only
`service_build_id` remains separate from scientific release IDs and data pins.

Target URLs follow this fixed pattern on the existing static origin:

```text
https://openantigens.org/api-runtime/{release_id}/{species}/{target_key}-{part}-{compressed_sha256}.oap
```

The filename checksum makes a physical repack a new URL, including a repack of
the same scientific release. The `.oap` file contains the exact verified gzip
bytes. Shared human and mouse Apache rules serve it as opaque octet-stream data
with no HTTP content encoding and an immutable success response. Compression
preserves its bytes. Actual MIME, cache and encoding behavior through the host
and edge remains a launch check.

Normal packaging adds only catalog-authorized, hash-verified assets under
`public_site/api-runtime/`. Matching existing bytes are reused. A collision with
different bytes or a symlink blocks packaging. All destinations are checked
before additions. Valid unreferenced assets can remain after a later failure;
the builder preserves them for controlled reuse and performs no rollback
cleanup of immutable data. It writes no app ZIP or review record until every
required local asset passes.

The hard ZIP ceiling is less than 100,000,000 bytes. Defaults reserve 20% of that
budget and of a 500,000,000-byte expanded app budget: effective limits are
80,000,000 and 400,000,000 bytes. The review record separates `logical_data_*`,
`physical_data_*`, static asset counts and bytes, and retained export copies.
It also records hashes, remaining headroom and gates needing evidence.
Compressed target data is limited to two parts and 4 MiB per target; decoded
target parts together are limited to 32 MiB. Oversize input blocks the whole
release. Summary records remain bounded to 1 MiB. The seven external projection
bulk formats have a separate 1 GiB staging limit; small schema/citation files
remain bounded to 1 MiB. Neither bound includes the original website disease TSV.

Standard output is compact. The complete app file map and authorized query asset
map stay in the separate artifact JSON. Its `artifact_record` stdout field is
an absolute local path and stays outside public provenance and the deployed ZIP.
A dry run reports that field as null.

Review dependencies from the locked service directory:

```sh
npm ci --ignore-scripts
npm test
npm audit --omit=dev --audit-level=high --json > /absolute/path/staging/npm-audit.json
```

Pass the saved report with `--dependency-audit-report` to record its hash. Review
its exit status, freshness and findings; the package builder records the report
without asserting that the audit or host gates passed.

## Update code against reviewed data

For a runtime-only fix, use the separate explicit repack command. Supply the
prior reviewed record and its ZIP checksum from the known full-validation
receipt through an independent trusted record. Reading a checksum from the
candidate record or computing it from the candidate ZIP does not establish that
independent pin.

```sh
PYTHONPATH=src python scripts/build_api_service_release.py repack \
  --reviewed-record /absolute/path/reviewed/app.artifact.json \
  --expected-zip-sha256 INDEPENDENT_REVIEWED_ZIP_SHA256 \
  --public-site-dir /absolute/path/snapshot/public_site \
  --output /absolute/path/staging/app-code-update.zip \
  --dry-run
```

The repack verifies the entire old archive, extracts only data and the old seal,
and runs the current query schema, retirement/privacy policy and every required
static asset/export/source-manifest checksum. It reuses every data byte, path,
seal and scientific ID, while replacing all old code with current allowlisted
sources, dependencies and licenses. The new record identifies the parent ZIP and
labels scientific validation as inherited, with current checks rerun. Existing
budgets, exclusive outputs and dry-run protections apply.

Use a normal full build for projection, public schema, transport or data changes,
or a data-boundary fix outside the current query validator's coverage. A current
policy failure blocks the repack; it rewrites no data and uses no fallback.
After review, repeat without `--dry-run`, run production installation/audit and
test the new code artifact before authorizing publication.

## Test and configure the app

Use Node 22 in a separate managed app outside the cPanel static document roots.
Reuse the same app ID for preview, publication, updates and rollback. Confirm the
selected plan's app, memory, disk, inode, startup, replica and transfer allowances
and charges before provisioning. The provider FAQ describes preview removal
after 30 days; verify its visibility and expiry before using it for staging.

```text
Build: npm ci --omit=dev --ignore-scripts
Start: npm start
```

Use the host-provided `process.env.PORT`, set `OPENANTIGENS_DATA_DIR=data`, and
configure exact allowed hostnames and origins. Keep the generated app seal with
its query data. Startup performs no origin asset fetch. Local readiness proves
local metadata validation; it does not establish remote target availability.
The runtime rejects corrupt, missing or changed data and follows no redirects.

Keep forwarding disabled until the actual ingress path is verified. The
supported configuration uses `OPENANTIGENS_CLIENT_IDENTITY_VERIFIED=1`, exact IPs
in `OPENANTIGENS_TRUSTED_PROXY_PEERS`, and `OPENANTIGENS_PROXY_HOPS` from 1 to 8.
Verify the trusted peers, header chain and absence of an ingress bypass before
enabling it. Invalid forwarding uses the socket peer for admission and supplies
no verified analytics identity.

Sanitized request records go to standard output, including refusals and aborts.
Use the separate [local traffic report workflow](API_MCP_ANALYTICS.md); its script
stays outside the app ZIP. Plausible export stays off until owner approval and
verification of plan, goals and trusted identity. Its filtered visitor measures
remain separate from the local request counts.

Use explicit known target queries for bounded HTTP and MCP checks:

```sh
python scripts/check_api_service.py \
  --base-url https://api.openantigens.org \
  --release-id EXACT_RELEASE_ID \
  --human-query EGFR --mouse-query Egfr
```

The checker requires positive full-summary matches in both species, exactly five
tools and only four schema/citation references. It checks compact provenance
against the immutable public manifest reference. Default pacing preserves the
20/minute peer rate and burst of eight. Limits are three attempts per read,
twenty attempts total and thirty seconds of cumulative waiting, honoring numeric
`Retry-After`. It follows no redirects and sends no credentials. For local tests,
bind with `OPENANTIGENS_HOST=127.0.0.1` and use the loopback URL.

## Publish and retain immutable data

Obtain approval for the checked artifacts, exact app/domain, plan costs and
rollback path. Recheck canonical `origin`. Local build commands publish nothing.

1. Set `OPENANTIGENS_PUBLICATION_STATE=maintenance` through trusted app controls.
   Verify that HTTP and MCP current discovery close with
   `publication_in_progress`; valid literal reads remain available.
2. Publish the static artifact through the protected deployer, supplying the
   reviewed app record:

   ```sh
   scripts/deploy_snapshot_to_godaddy.sh \
     --existing-snapshot /absolute/path/snapshot \
     --api-artifact-record /absolute/path/staging/openantigens-api.artifact.json \
     --dry-run
   ```

   After review and authorization, repeat without `--dry-run`. The deployer
   verifies the reviewed ZIP and sealed catalog, derives the exact asset/export
   inventory and hashes local bytes. Immutable files transfer additively after
   collision checks; checksum readback must converge before the general mirror.
   Unknown files under `api-runtime/` or `data/releases/` never enter that mirror.
3. Keep both immutable prefixes protected from receiver deletions, including
   `--delete-before`. Old and new app overlap can require three scientific
   releases, plus distinct physical asset sets for a same-release code repack.
   Retain every asset needed by the active and rollback app catalogs. Cleanup
   requires a separate approval after old traffic drains and rollback retention
   is reviewed; the deployer performs no automatic garbage collection.
4. Verify normal static URLs, checksums, opaque `.oap` headers and both portals
   through the actual edge. Update the same managed app ID with the checked ZIP.
   Keep current discovery held while validating literal HTTP and real MCP reads,
   target availability, both species, source counts and representative exports.
5. Run the checker with `--expect-maintenance` and, after its HTTPS checks pass,
   `--activate-discovery-local-site /absolute/path/snapshot/public_site`. This
   updates both local portal metadata files with verified service/MCP URLs.
   Publish those metadata files through the authorized static workflow and
   verify the public copies. Set publication state to `active`, restart, and
   repeat normal discovery/client checks before leaving an external hold.

The static copy updates files in place. Keep a controlled publication window;
separate app and static updates do not form an atomic transaction. A failed
upload, collision, missing asset or mismatch keeps current discovery held.
Pinned reads use their own source namespace and never substitute another release.

For a code fault, restore a compatible previous app with its identical data pins.
For a data mismatch, close discovery, restore matching human/mouse mutable pages
and metadata, then restore the matching app. Immutable historical assets remain
in place. Purge and verify actual mutable caches, recheck both interfaces, and
reopen discovery only after recovery passes. Preserve controlled local rollback
artifacts and enough capacity; remove no user material to create headroom.
