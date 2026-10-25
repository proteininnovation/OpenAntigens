# API and MCP release operations

Build the Node service from a frozen, validated public data bundle. Keep the code
build ID separate from the data release ID. A code update can preserve the same
human and mouse scientific release.

The full-catalogue local trial on October 25, 2026 measured a 245,426,866-byte
ZIP, above the 100,000,000-byte hard ceiling and 80,000,000-byte effective budget.
The current self-contained GoDaddy launch is blocked. A separately sealed local
runtime passed both publication states and all five HTTP/MCP read tools; it
reached readiness in about 60 seconds and peaked at 1,213,923,328 bytes RSS.
Those functional diagnostics do not establish an admissible upload or managed
hosting capacity. The rejected ZIP was withheld from publication.
The commands below retain the existing gate. Keep the complete scientific
contract and select a separately approved hosting target with sufficient
artifact, memory and startup allowances before launch.

The public API projection includes up to 100 human disease associations per
target in the published website's source order. Evidence metadata reports the
projected count, full recorded source count and selection coverage. The complete
projection download covers these selected public records. The website's
original full disease-association TSV remains a separate static download
outside the app bundle. Interactive searches use gene names and accessions;
disease-ID filtering is outside this bounded projection's query contract.
Check `summary.disease_coverage` and disease-page `coverage`: the selection is
`portal_top_100`, with separate `projected_count` and `source_count` values.
Unknown source counts stay null. `source_reference` labels the original TSV
link as live and records its checksum; that URL can advance with the website.

Read the [public-safe host verification summary](API_MCP_HOST_VERIFICATION.md) for endpoint checks and remaining launch gates. Account-specific evidence stays in the ignored local folder; local package and smoke results cover the prepared service only.

## Prepare the local artifact

Run these commands from the canonical `proteininnovation/OpenAntigens` checkout
with the Python package installed. Build both portals and finish static
compression and asset removal first. Then project the frozen human and mouse
reports against that final `public_site` directory, which contains `mouse/`.
Use the frozen scientific/projector commit and snapshot identifier:

```sh
PYTHONPATH=src python -m agdesign2.cli build-public-api \
  --human-summary /absolute/path/frozen-human/batch_summary.json \
  --mouse-summary /absolute/path/frozen-mouse/batch_summary.json \
  --public-site /absolute/path/public_site \
  --output-dir /absolute/path/api-data \
  --software-commit SCIENTIFIC_PROJECTOR_COMMIT_SHA \
  --snapshot-id FROZEN_SNAPSHOT_ID \
  --previous-bundle /absolute/path/previous-api-data

PYTHONPATH=src python -m agdesign2.cli validate-public-api /absolute/path/api-data
```

Replace the commit and snapshot placeholders with the selected frozen values.
For the first release, omit `--previous-bundle`. The output directory must be
new or empty. The builder copies current and retained previous downloads into
the final static artifact and checks their hashes. Preserve those bytes through
static publication; revalidate and repackage if any final asset changes.

Package the matching application and data only after these checks pass:

```sh
PYTHONPATH=src python scripts/build_api_service_release.py \
  --bundle-dir /absolute/path/api-data \
  --public-site-dir /absolute/path/public_site \
  --output /absolute/path/staging/openantigens-api.zip \
  --dry-run
```

Repeat without `--dry-run` to create the ZIP and the separate
`openantigens-api.artifact.json` review record. Existing output files cause an
error. A dry run validates inputs and measures a real temporary ZIP while
preserving the requested output paths and all inputs.
Standard output contains the measurements and gates plus an absolute local
`artifact_record` path. The full file-hash map stays in the separate artifact
JSON; local paths stay outside the deployed ZIP and public provenance. A dry
run reports `artifact_record: null` because it writes no review record.

The ZIP contains the allowlisted Node source, package manifest and lock,
`LICENSE`, `DATA_LICENSE.md`, `service-build.json`, and sealed physical runtime
data. Before packing, the builder validates every logical public shard and all
eleven approved static download roles. Manifest membership alone authorizes no
extra file. The packager independently checks each role, path, media type,
encoding and reference, scans public text, and checks citation exports against
their frozen source.

Runtime packing groups related target JSON shards into bounded gzip containers
that retain each original logical path and exact bytes. Public manifests,
release IDs, schemas and checksums remain unchanged. The app keeps indexes,
manifests and the four small schema/citation references as ordinary files.
The seven website bulk roles (JSONL, TSV, FASTA and complete-record gzip files)
stay in static publication and source staging; their bytes stay outside the
app ZIP. Complete export-byte parity is verified at staging. Runtime metadata
does not establish remote static-download byte parity.
Each of those seven external projection exports has a separate 1 GiB
engineering limit, checked with streamed staging hashes and export parity.
The four small schema/citation files remain limited to 1 MiB each. This external
file limit is separate from query-memory bounds and the physical app budgets;
measure static storage and transfer capacity for the actual release. The
original full disease TSV is outside these eleven projection download roles.

`service-build.json` records the independent code build ID and a `data_transport`
seal with the source root-manifest SHA-256 and transport-index SHA-256. Startup
requires those pins for packed data. Keep the generated seal with the exact
runtime directory; packed mode uses bounded reads and creates no extracted files.
The builder verifies the canonical origin, both species, canonical release IDs,
manifest checksums, every retained file, complete-record downloads, and matching
final static copies. It rejects symlinks, unexpected data files, unsafe paths,
changed bytes, missing static copies, and exceeded local size budgets.

The full local engineering profile was compiled with local projector commit
`cce700c05861442ba93657576ad781fe0cfa7bc3`. Its exact Git tree
`4317716f47968a0924cb6919dd6be73c86307a3d` is available in published commit
`f77344d`. The profile retains its original release ID and provenance. For a
future launch build, use the selected published canonical revision and finalize
both portals and their privacy/guidance pages before projecting the data.

The hard ZIP ceiling is **less than 100,000,000 bytes**. Defaults reserve 20% of a
100,000,000-byte compressed budget and a 500,000,000-byte expanded budget. These
are local engineering limits. Check actual subscription limits, expanded disk,
inodes, build and runtime memory, process count, and startup time before
publication. The artifact records measured physical compressed and expanded sizes,
per-file hashes, remaining local headroom, and host checks that still require
evidence. It also measures the retained static download copies and the minimum
combined app/download storage. Reserve capacity for the rest of both portals,
installed Node dependencies, logs and the previous rollback artifacts separately.
Separate `logical_data_*` and `physical_data_*` measurements distinguish the
validated scientific dataset from its deployment representation.
If the measured file count requires ZIP64, verify the host's extractor supports
that archive before upload.
The packager also enforces the service's startup read limits: a 1 MiB root
manifest, 96 MiB release manifests, 32 MiB per species index, and 64 MiB of
compact indexes across all retained species and releases. Target shards remain
limited to 8 MiB, and evidence, construct, sequence and schema shards to 1 MiB.
The measured full-catalogue file-map lower bound of 71,707,842 bytes exceeds the
earlier 64 MiB manifest bound. The 96 MiB manifest limit and 64 MiB aggregate
index limit are engineering trial bounds. The transport index and each decoded container
are limited to 16 MiB; each compressed container is limited to 17 MiB. Measure
the selected full catalogue's actual ZIP, expanded files, startup and peak
memory before approving a package for the selected host. A package that exceeds
these bounds needs a measured data/runtime design change before publication.

Review dependencies from the locked service directory:

```sh
npm ci --ignore-scripts
npm test
npm audit --omit=dev --audit-level=high --json > /absolute/path/staging/npm-audit.json
```

Pass the saved JSON with `--dependency-audit-report` to record its hash. The
builder labels the report as provided and unverified: review its exit status,
freshness, dependency tree, and findings before publication. Keep test and audit
evidence beside the artifact, outside the deployed ZIP.

Use the bounded Python checker for a selected service and literal release:

```sh
python scripts/check_api_service.py \
  --base-url https://api.openantigens.org \
  --release-id YOUR_64_CHARACTER_RELEASE_ID \
  --human-query EGFR \
  --mouse-query Egfr
```

Choose explicit target queries from the selected release. The checker makes
twelve logical HTTP and MCP reads, checks both species, and compares metadata
and records across both interfaces. Release references contain only the schema
and three citation formats; the service exposes five tools and no bulk download
discovery route or tool.
Provenance includes the source commit and snapshot plus a counted, checksummed
reference to the immutable manifest's input-hash map; the checker verifies that
reference while the full map remains in the stored release manifest.
For a local service, use `http://127.0.0.1:<port>` and set
`OPENANTIGENS_HOST=127.0.0.1` when starting the service. Remote checks require
HTTPS. The checker accepts no credentials and follows no redirects. Add
`--expect-maintenance` to verify that HTTP and MCP current discovery are closed
while literal-release reads remain available. Pinned metadata then omits the
current ID and retained list. Use `--retained-release-id` for a caller-known
previous release: the checker probes its literal metadata and both species
through HTTP and MCP, for eighteen logical reads. Active discovery must match
the expected retained set; maintenance checks use the supplied IDs. Checks pace
requests for the default 20/minute peer limit and burst of eight. They honor
numeric `Retry-After`, with at most three attempts per read, twenty attempts in
total, and thirty seconds of cumulative waiting. Load tests need separate checks.

## Configure the managed Node app

Use a separate managed Node app with a stable app ID and a dedicated hostname,
such as `api.openantigens.org`. Its application directory belongs outside the
cPanel static document roots. Configure Node **22**, use the host-provided
`process.env.PORT`, and select these commands:

```text
Build: npm ci --omit=dev --ignore-scripts
Start: npm start
```

The source uses plain JavaScript, so it needs no transpilation. Set
`OPENANTIGENS_DATA_DIR=data` relative to the extracted app directory. Set the
allowed hostnames and origins explicitly for the selected domain and preview.
Keep account credentials and static-site SSH keys outside the app. Check the
runtime environment and TLS/proxy behavior using the selected host's controls.

The service writes sanitized request records to standard output by default,
including refusals and aborts. Archive the complete intended log interval and
use the separate [local traffic report workflow](API_MCP_ANALYTICS.md). The
report script stays outside the deployed ZIP. Plausible export stays off until
the owner approves it and verifies the plan, goals and trusted client identity.
Its bot filtering and coverage are separate from the local request counts.

Keep forwarded-client identity disabled until the actual ingress path is
verified. The supported configuration uses
`OPENANTIGENS_CLIENT_IDENTITY_VERIFIED=1`, an exact IP allowlist in
`OPENANTIGENS_TRUSTED_PROXY_PEERS`, and `OPENANTIGENS_PROXY_HOPS` from 1 to 8.
Confirm the trusted peers, appended header chain and absence of an ingress
bypass before enabling it. Untrusted or malformed forwarding falls back to the
socket peer for admission and supplies no verified analytics identity.

Confirm the account's available app slot, tier, deployment method, capacity and
cost before creating an app. Reuse the **same app ID** for preview, publication,
updates and rollback. A preview can be public; use a restricted preview only
when the host supplies and verifies that restriction. Preserve the previous
verified ZIP, artifact record, static artifact and configuration on controlled
local storage.
The provider FAQ describes preview removal after 30 days. Confirm the selected
app's preview visibility and expiry, publishing scope and charges before using
it for staging. The artifact leaves that account-specific gate unverified.

## Publish both species and the service

Obtain approval for the concrete artifacts, app/domain, cost, account changes and
rollback path. Recheck canonical `origin` before publication. The local build
commands perform no upload or deployment.

1. Set `OPENANTIGENS_PUBLICATION_STATE=maintenance` through the selected app's
   trusted host configuration and verify the restarted app. This closes
   current-release discovery. Verify
   `503 publication_in_progress`, `Retry-After` and `Cache-Control: no-store`.
   Previously retained literal-release reads remain available when validated.
   Use a verified host maintenance control when the app cannot keep discovery
   closed through activation.
2. Publish the finalized human and mouse static artifacts together through the
   protected static deployer. Include every advertised retained release under
   `/data/releases/{release_id}/`. Purge the actual static cache and verify both
   species and downloads through ordinary public URLs.
3. Update the same Node app ID with the checked ZIP and verify startup,
   readiness and routing. Keep discovery closed while checking the final domain
   through its actual TLS proxy. Run HTTP smoke checks and the Node MCP contract
   tests, then check real MCP clients.
4. Confirm that both portal metadata files and literal API and MCP reads name
   the same release. Compare final manifest hashes, species counts,
   representative records and sequence exports. Activate and publish local
   discovery metadata with the procedure below while current discovery stays
   closed.
   The hold controls dynamic current and retained-release discovery. Literal
   release IDs, static immutable downloads and historical retirement evidence
   remain readable.
5. Set `OPENANTIGENS_PUBLICATION_STATE=active` through the same trusted host
   configuration and verify the restart. Run the normal smoke check and final
   client checks before leaving any external maintenance gate.

During step 4, run the checker against the final HTTPS origin with
`--expect-maintenance` and
`--activate-discovery-local-site /absolute/path/public_site`. The checker first
verifies closed current discovery and both species' literal reads. The explicit
activation option then updates both local `portal_metadata.json` files:
it preserves their release metadata, sets `public_api.state` to `active`, and
adds `service_url` and `mcp_url`. It checks both local release IDs before writing;
a failed smoke leaves both files unchanged. Publish these two metadata files
through the authorized static workflow while the publication gate remains
controlled, then verify the ordinary public URLs before step 5. Shared
agent guidance reads those metadata fields, so both portals gain the verified
connection instructions together.

The static rsync copy updates files in place. It can expose mixed content during
the transfer. Schedule that publication window and retain an external
maintenance gate when the service configuration would reopen discovery on
startup. Treat failed transfers, failed starts and release mismatches as an
incomplete publication; keep discovery closed until recovery is verified.

## Roll back

For a code fault, restore a compatible previous app build with the same data
release ID. For a data or portal mismatch, close discovery, restore the matching
previous human and mouse static artifact, purge and verify the static cache,
then restore the matching app/current marker. Verify both interfaces before
reopening discovery. A literal bad release returns its error; preserve its ID
rather than substituting another release's data.

Retain current and previous validated data only while the measured package and
host budgets permit both. Include their static downloads in each new static
artifact. A previously known retired release has published tombstone metadata
and an archive URL when one exists. The bundle keeps at most 1,000 tombstones;
IDs outside that recorded history return the unknown-release response. Test
recovery from a failed static copy and failed app start before
launch. Keep rollback capacity available; deleting user data or backups requires
specific authorization.
