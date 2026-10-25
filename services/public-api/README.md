# OpenAntigens public service

A Node.js 22 application serves frozen public human and mouse records through the HTTP API and five read-only MCP tools. Both interfaces use the same manifest-backed provider. Scientific calculations stay in the Python release builder.

## Build and run

```sh
cd services/public-api
npm ci --ignore-scripts
npm run build
OPENANTIGENS_DATA_DIR=/absolute/path/to/bundle npm start
```

`npm run build` validates the application imports. JavaScript ESM runs directly; no transpilation or native extension is needed. Production installation uses `npm ci --omit=dev --ignore-scripts`. The managed app entry point is `src/server.js`.

The full October 25 catalogue works in local runtime checks but its measured 245 MB ZIP exceeds GoDaddy's documented managed upload limit below 100 MB. The packager enforces that limit. Public launch is blocked pending a separately approved hosting decision; see [host verification](../../docs/API_MCP_HOST_VERIFICATION.md).

The Python bundle builder must finish and validate both portals before the service starts. `OPENANTIGENS_DATA_DIR` names the directory that contains the root `manifest.json`. Startup rejects duplicate JSON object members, then checks release IDs, schemas, sizes, hashes, closed artifact roles, target-owned references and the complete scientific fingerprint. A failed check leaves `/readyz` unavailable. The service keeps validated release metadata and compact search rows in memory; selected summaries load on demand. The service makes no upstream scientific requests.

Packaged data uses hashed gzip containers whose base64 members preserve the exact original logical JSON bytes. The app-code `service-build.json` supplies `data_transport` pins for the root manifest and transport index. Packed data without those independent pins is refused. Every container/member is checked at startup and bounded request reads recheck integrity, with at most two request-local containers of 16 MiB decoded each. There is no extraction or runtime data writing. This seal assumes that a data-directory writer cannot modify app code; deployment permissions and isolation remain hosting checks.

Unpacked developer bundles retain complete local export and gzip-parity validation. For packed runtimes, the original website's seven bulk-export bodies stay outside the app. Their exact approved metadata remains sealed; the producer and package builder verify their bytes and parity before packaging. Only these seven static bulk roles may be up to 1 GiB each; unpacked validation hashes them as streams. The four schema/citation references remain bounded at 1 MiB. API payload and query-shard limits stay unchanged. Runtime readiness verifies query data and its physical inventory, and makes no remote download-byte verification claim.

| Setting                          | Default                                                 | Purpose                                                              |
| -------------------------------- | ------------------------------------------------------- | -------------------------------------------------------------------- |
| `PORT`                           | `3000`                                                  | Managed host listener port                                           |
| `OPENANTIGENS_HOST`              | `0.0.0.0`                                               | Managed ingress binding; use `127.0.0.1` for a local-only listener   |
| `OPENANTIGENS_DATA_DIR`          | `./data`                                                | Validated public bundle directory                                    |
| `OPENANTIGENS_ALLOWED_HOSTS`     | `api.openantigens.org,localhost,127.0.0.1,[::1]`        | Exact hostnames; excludes ports                                      |
| `OPENANTIGENS_ALLOWED_ORIGINS`   | `https://openantigens.org,https://www.openantigens.org` | Exact browser origins, including scheme and port                     |
| `OPENANTIGENS_PUBLICATION_STATE` | unset                                                   | Use `active` (default) or `maintenance`; invalid values stop startup |

`OPENANTIGENS_ANALYTICS_CLIENT_SECRET` supplies the optional HMAC secret. `OPENANTIGENS_ANALYTICS_PLAUSIBLE=1` requests collector delivery; `OPENANTIGENS_CLIENT_IDENTITY_VERIFIED=1` confirms that the deployment's client identity path has been verified. Both flags default off. Direct mode uses the socket peer and ignores forwarded headers. Verified proxy mode also requires `OPENANTIGENS_TRUSTED_PROXY_PEERS`, a comma-separated exact IP allowlist, and `OPENANTIGENS_PROXY_HOPS`, an integer from 1 to 8. Only a listed socket peer can supply one strictly parsed `X-Forwarded-For` header; the configured hop count selects from its right end. Malformed, duplicate or short chains use only the raw peer's rate bucket and provide no analytics identity. Verify the actual proxy chain, header behavior and direct-app bypass isolation before enabling this mode.

Reverse proxies must preserve the expected Host. Forwarded Host and Origin headers do not override the direct request. Forwarded client identity is used only under the verified exact-proxy configuration above. Limits apply per process and socket peer: 20 requests per minute per client with a burst of 8, 90 globally with a burst of 10, two concurrent admitted requests globally, one per client, and a five-second handler budget. `Retry-After` reports the token refill or concurrency wait in seconds. There is no waiting queue. A shared ingress IP therefore shares its client quota; host launch checks must verify proxy behavior and any edge quota.

API and MCP responses include compact provenance: software/snapshot identity, recorded source versions/dates, notes, and an input-hash count with an immutable manifest URL, byte checksum and JSON pointer. The complete input-hash map stays in that manifest. Storage schemas and release hash validation retain the full map. The small `projection_methods` descriptor declares the recorded disease selection and order before a target is fetched; legacy releases report null when this method was unrecorded.

The runtime writes a sanitized local JSON ledger for every request handled by the app, including rejections, rate limits and aborts. Records contain fixed operation/species categories, status, duration, response bytes and a request ID. Raw URLs, queries, target IDs, user agents, client addresses, annotations, environment values and exceptions stay out of the ledger. Traffic blocked upstream or rejected before the app's HTTP request handler falls outside this coverage. There are no application sessions, cookies or administrative routes.

Plausible delivery is off by default. Explicit opt-in also requires verified client identity; only admitted business lookups and their known business failures are eligible. Rejected spam, health checks and protocol handshakes remain in the local ledger. Delivery records distinguish disabled, skipped, filtered, accepted and failed events. Collector attempts are capped at 100 per process per UTC day, with four in flight and no retry queue. Restarts or replicas can reset or multiply that cap and require an operator quota check. Tests and profiles use no live collector.

Daily HMAC client IDs use the UTC day, verified client IP and a safe authentic user agent. They require a secret of at least 32 characters. They support approximate UTC-day visitor counts, not a count of people. Without those prerequisites the ledger records no client ID. Plausible country/city and visitor reports cover only eligible delivered events and remain subject to its bot filtering.

## HTTP and MCP

`GET /v1/releases` and `GET /v1/releases/current` discover immutable IDs. Data routes require a literal release ID. `/openapi.json` publishes exact arguments, output schemas, errors and limits. `/healthz` reports process liveness; `/readyz` reports bundle readiness. GET/HEAD immutable data supports ETags and 304 responses. Discovery, health, errors and MCP use `no-store`.

`POST /mcp` uses official `@modelcontextprotocol/server` 2.4.0 and `@modelcontextprotocol/node` 2.1.2 with `responseMode: 'json'` and `legacy: 'stateless'`. It exposes `get_release`, `search_targets`, `get_target`, `list_constructs` and `get_sequence`. Current MCP 2026-07-28 uses per-request negotiation with JSON responses. The SDK's legacy 2025-11-25 flow returns finite one-shot SSE responses. That compatibility behavior is covered by tests and needs verification through the production proxy. GET and DELETE do not open MCP streams; there is no session store or resumable event history. JSON-RPC batches are rejected at the bounded request boundary.

Release metadata links only schema and citation references. Static website downloads remain available through the portals; the API and MCP expose interactive lookups rather than catalogue downloads. One oversized result can be recovered through its complete static-record fallback.

Each FASTA record starts with a unique encoded release/target/construct identifier for standard parsers. Its metadata headers retain the release ID, actual identifier and sequence source/accession, sequence basis, exact warnings and evidence, immutable manifest and complete-record URLs, and a sequence record selector. Header values use RFC3986 encoding; warning/evidence lists and the selector use canonical JSON. HTTP FASTA exports also supply metadata Link and attachment headers. The family section projects stored `family_context`; portal evidence views provide full-length family matrices and construct homolog exports. Sequence routes return exact released sequences for the requested species. `records.jsonl.gz` contains the complete public projection.

During `maintenance`, HTTP current discovery and MCP `get_release` without an ID return `publication_in_progress`; literal release reads continue. Pinned metadata omits current/retained discovery lists and returns `discovery_state: publication_in_progress`. HTTP supplies `Retry-After: 60`. Maintenance is a trusted deployment setting, changed outside the public API. The publication hold controls dynamic discovery consistency. Literal release records, including retirement history, stay readable. The root manifest retains at most 1,000 recent release tombstones. Discovery returns their count and the literal-release lookup endpoint; direct retired-ID lookups return the recorded tombstone. Known retired IDs return `410 retired_release` with the retirement ID and recorded archive availability; unknown IDs return 404.

Ordinary results are limited to 256 KiB; sequences to 1 MiB. Overflow errors return the manifest-listed complete `records.jsonl.gz` download. Sequence overflow also supplies the exact target/construct `record_selector` for recovery. Collection defaults are 5 records, maximum 10, including similarity sections. Target search requires at least two trimmed literal query characters. Its fixed window contains at most 20 matches; `refine_query` signals that a narrower query is needed. Wildcards and unfiltered catalogue reads are rejected. Human disease sections retain the first 100 frozen rows in published source order, matching the portal. Their `coverage` reports the original source count, projected count, fixed limit and live source reference with its recorded SHA-256. Mouse coverage is not applicable. The complete public projection includes this declared selection; the original website TSV retains the full source feed. `disease_id` is rejected as a search argument for this partial view. Storage limits differ: each species index 32 MiB, all retained indexes together 64 MiB, release manifests 96 MiB, target shards 8 MiB, and section/construct/sequence shards 1 MiB. The 96 MiB manifest and 64 MiB aggregate index caps are a measurement trial for the full catalogue. Actual full-release startup memory, package size and host capacity must pass before launch. Unknown, duplicate or inappropriate arguments are rejected. Coordinates are 1-based inclusive in the declared sequence basis.

Processing success does not establish experimental validation. Model scores and predictions retain their recorded meaning. Sequence similarity does not establish specificity, epitope identity or cross-reactivity. Empty results do not establish biological absence. Cite the release, source accession and supplied attribution.

## Verification

```sh
npm test
OPENANTIGENS_TEST_DATA_DIR=/absolute/path/to/python-built/bundle npm test
OPENANTIGENS_PACKED_FIXTURE=/absolute/path/to/sealed/runtime/data npm test
npm audit --omit=dev --audit-level=high
```

Tests use `node:test`, official MCP clients in current and legacy modes, and independent JSON-RPC wire requests. They check species parity, exact sequences, filter/cursor behavior, malformed messages, protocol header mismatches, methods/media/CORS, quotas, deadlines, maintenance, overflow downloads, byte changes, symlinks and readiness. The optional environment variable runs the same HTTP and MCP provider checks against the actual Python-generated bundle.

The producer schema snapshot in `test/schema.json` is test data. Runtime schemas always come from the selected hash-verified release. Local checks establish application behavior. Production acceptance still requires the intended agent client, actual host memory/cold-start/concurrency measurements, the final proxy/WAF and the approved publication/rollback procedure.
