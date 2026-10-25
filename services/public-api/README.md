# OpenAntigens public service

A Node.js 22 application serves frozen public human and mouse records through the HTTP API and six read-only MCP tools. Both interfaces use the same manifest-backed provider. Scientific calculations stay in the Python release builder.

## Build and run

```sh
cd services/public-api
npm ci --ignore-scripts
npm run build
OPENANTIGENS_DATA_DIR=/absolute/path/to/bundle npm start
```

`npm run build` validates the application imports. JavaScript ESM runs directly; no transpilation or native extension is needed. Production installation uses `npm ci --omit=dev --ignore-scripts`. The managed app entry point is `src/server.js`.

The Python bundle builder must finish and validate both portals before the service starts. `OPENANTIGENS_DATA_DIR` names the directory that contains the root `manifest.json`. Startup checks release IDs, schemas, file sizes, SHA-256 values, paths, download identities and symlinks. It reconstructs the semantic release hash from bounded shards and checks exact coverage of the complete-record gzip fallback. A failed check leaves `/readyz` unavailable. Startup streams hashes and retains compact indexes; requests read bounded manifest-listed shards. The service makes no upstream scientific requests.

| Setting                          | Default                                                 | Purpose                                                                 |
| -------------------------------- | ------------------------------------------------------- | ----------------------------------------------------------------------- |
| `PORT`                           | `3000`                                                  | Managed host listener port                                              |
| `OPENANTIGENS_HOST`              | `0.0.0.0`                                               | Managed ingress binding; use `127.0.0.1` for a local-only listener      |
| `OPENANTIGENS_DATA_DIR`          | `./data`                                                | Validated public bundle directory                                       |
| `OPENANTIGENS_ALLOWED_HOSTS`     | `api.openantigens.org,localhost,127.0.0.1,[::1]`        | Exact hostnames; excludes ports                                         |
| `OPENANTIGENS_ALLOWED_ORIGINS`   | `https://openantigens.org,https://www.openantigens.org` | Exact browser origins, including scheme and port                        |
| `OPENANTIGENS_PUBLICATION_STATE` | unset                                                   | Set `maintenance` to pause current release discovery during publication |

Reverse proxies must preserve the expected Host. Forwarded Host, Origin and client IP headers do not override the direct request. Limits apply per process and socket peer: 60 requests per minute per client, 300 globally, eight concurrent admitted requests, and a five-second handler budget. There is no waiting queue. A shared ingress IP therefore shares its client quota; host launch checks must verify proxy behavior and any edge quota.

API and MCP responses include compact provenance: software/snapshot identity, recorded source versions/dates, notes, and an input-hash count with an immutable manifest URL, byte checksum and JSON pointer. The complete input-hash map stays in that manifest. Storage schemas and release hash validation retain the full map.

Logs contain fixed startup/readiness categories and the listener port. Requests, query strings, client addresses, annotations, environment values and exceptions are not logged. There are no application sessions, cookies, credentials or administrative routes.

## HTTP and MCP

`GET /v1/releases` and `GET /v1/releases/current` discover immutable IDs. Data routes require a literal release ID. `/openapi.json` publishes exact arguments, output schemas, errors and limits. `/healthz` reports process liveness; `/readyz` reports bundle readiness. GET/HEAD immutable data supports ETags and 304 responses. Discovery, health, errors and MCP use `no-store`.

`POST /mcp` uses official `@modelcontextprotocol/server` 2.4.0 and `@modelcontextprotocol/node` 2.1.2 with `responseMode: 'json'` and `legacy: 'stateless'`. It exposes `get_release`, `search_targets`, `get_target`, `list_constructs`, `get_sequence` and `get_downloads`. Current MCP 2026-07-28 uses per-request negotiation with JSON responses. The SDK's legacy 2025-11-25 flow returns finite one-shot SSE responses. That compatibility behavior is covered by tests and needs verification through the production proxy. GET and DELETE do not open MCP streams; there is no session store or resumable event history. JSON-RPC batches are rejected at the bounded request boundary.

Each FASTA record starts with a unique encoded release/target/construct identifier for standard parsers. Saved FASTA headers retain the release ID, actual identifier and sequence source/accession, sequence basis, exact warnings and evidence, immutable manifest and complete-record URLs, and a sequence record selector. Header values use RFC3986 encoding; warning/evidence lists and the selector use canonical JSON. HTTP FASTA exports also supply metadata Link and attachment headers. The family section projects stored `family_context`; portal evidence views provide full-length family matrices and construct homolog exports. Sequence routes return exact released sequences for the requested species. `records.jsonl.gz` contains the complete public projection.

During `maintenance`, HTTP current discovery and MCP `get_release` without an ID return `publication_in_progress`; literal release reads continue. Pinned metadata omits current/retained discovery lists and returns `discovery_state: publication_in_progress`. HTTP supplies `Retry-After: 60`. Maintenance is a trusted deployment setting, changed outside the public API. The publication hold controls dynamic discovery consistency. Literal release records, including retirement history, stay readable. The root manifest retains at most 1,000 recent release tombstones. Known retired IDs return `410 retired_release` with the retirement ID and recorded archive availability; unknown IDs return 404.

Ordinary results are limited to 256 KiB; sequences to 1 MiB. Overflow errors return the manifest-listed complete `records.jsonl.gz` download. Collection defaults are 25 records, maximum 100; similarity sections have a maximum of 50. Storage limits differ: each species index 32 MiB, all retained indexes together 64 MiB, release manifests 64 MiB, target shards 8 MiB, and section/construct/sequence shards 1 MiB. The 64 MiB manifest and aggregate index caps are a measurement trial for the full catalogue. Actual full-release startup memory, package size and host capacity must pass before launch. Unknown, duplicate or inappropriate arguments are rejected. Coordinates are 1-based inclusive in the declared sequence basis.

Processing success does not establish experimental validation. Model scores and predictions retain their recorded meaning. Sequence similarity does not establish specificity, epitope identity or cross-reactivity. Empty results do not establish biological absence. Cite the release, source accession and supplied attribution.

## Verification

```sh
npm test
OPENANTIGENS_TEST_DATA_DIR=/absolute/path/to/python-built/bundle npm test
npm audit --omit=dev --audit-level=high
```

Tests use `node:test`, official MCP clients in current and legacy modes, and independent JSON-RPC wire requests. They check species parity, exact sequences, filter/cursor behavior, malformed messages, protocol header mismatches, methods/media/CORS, quotas, deadlines, maintenance, overflow downloads, byte changes, symlinks and readiness. The optional environment variable runs the same HTTP and MCP provider checks against the actual Python-generated bundle.

The producer schema snapshot in `test/schema.json` is test data. Runtime schemas always come from the selected hash-verified release. Local checks establish application behavior. Production acceptance still requires the intended agent client, actual host memory/cold-start/concurrency measurements, the final proxy/WAF and the approved publication/rollback procedure.
