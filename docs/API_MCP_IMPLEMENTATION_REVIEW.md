# API/MCP implementation review

The implementation follows `API_MCP_PLAN.md` in the `feat/public-api-mcp` worktree branch. The publication target is a PR in `proteininnovation/OpenAntigens`. Production deployment and account changes require separate authorization.

## Review seats and completion rule

The same independent agents remain assigned throughout the review rounds:

- Software engineer: implementation, callers, contracts, tests, resource bounds, and deployment integration.
- Protein scientist: sequence/identity/coordinate correctness, evidence and warning preservation, useful results, and human/mouse distinctions.
- Security reviewer: data exposure, input/protocol boundaries, file integrity, abuse controls, dependency/package scope, and publication safety.
- Interface reviewer: HTTP/MCP/CLI discoverability, error recovery, documentation, and shared human/mouse agent guidance.

Each round reviews actual code and runnable outputs. The author judges each finding, implements the required corrections, and resumes the same reviewers for a full reread. All seats must pass with no remaining material implementation defects before the PR. Review agreement is distinct from experimental validation or production hosting verification.

## Implementation decisions

- Keep the existing Python scientific pipeline and compile explicit public data from frozen reports and finalized portal indexes. Node serves the validated immutable projection through one provider shared by HTTP and MCP.
- Keep target identifiers separate from sequence accessions/sources and preserve the seven construct categories and raw classification.
- Native target-species exports use stored released sequence values. The machine builder performs no homolog remapping and does not consult ambient caches or fetch scientific data.
- Official MCP server 2.4.0 returns JSON for the current protocol and a finite one-response SSE body for legacy stateless calls. Preserve its supported framing with request/result limits, rather than add a custom framing conversion. Long-lived GET streams and sessions remain outside the service.
- Shared guidance connects to a service only when deployment metadata advertises active URLs. Building a package alone does not advertise a live endpoint.

## Review rounds

### Round 1

All four independent seats reviewed actual code and local outputs. Verdict: fixes required. The author accepted the following material findings; duplicate findings from different seats share a correction.

| Finding | Priority | Decision and required correction |
| --- | --- | --- |
| Protein scientist: final predicted model can disagree with the target reference | P1 | Accept. Validate the actual final published model with the shared sequence-match helper; reject a mismatch rather than trust an index flag. |
| Protein scientist: selected-region cysteine warning reused for the full reference | P2 | Accept. Bind warning text to the exported interval and preserve the recorded evidence/context without claiming an in-range partner lies outside it. |
| All correctness/interface seats: FASTA loses source, basis, warnings and release provenance | P2 | Accept. Make HTTP and bulk writers losslessly preserve the same release/source/basis/warning/evidence metadata and immutable references; test decoding and exact sequence bytes. |
| Protein scientist: positive literature counts have no usable query link | P2 | Accept. Preserve the actual frozen portal query link or use the same documented offline URL construction with clear provenance; perform no new literature retrieval. |
| Security: metadata symlink permits an external file overwrite | P2 | Accept. Preflight both destination metadata files and preserve containment through writes; test both species and unchanged outside-file bytes. |
| Software engineer: explicit-ID metadata bypasses maintenance discovery hold | P2 | Accept. Serve pinned metadata without current/retained discovery during maintenance; keep data reads available and adapt the independent smoke checker. |
| Software engineer: a third build forgets a known retired ID | P2 | Accept. Propagate bounded tombstone metadata, distinguish 410 from unknown 404, and preserve retained artifact bytes. |
| Interface reviewer: health/readiness/OpenAPI responses use incorrect schemas | P2 | Accept. Describe the actual bodies, statuses, and headers in OpenAPI and exercise them against the server. |
| Interface reviewer: OpenAPI links an absent `/agents.md` | P2 | Accept. Link the actual shared agent guide. |
| Interface reviewer: build/validate CLI examples absent | P3 | Accept. Add explicit frozen-input commands to the deployment runbook. |

The scientist also recorded coverage observations about full-length family context and construct-level homolog sequences. V1 exports exact target-species native sequences and the declared public projection; it performs no homolog remapping. These observations do not invalidate the checked exports. Keep the projection's scope explicit in the interface documentation; any later homolog-export extension needs the shared frozen mapping contract and a scientific review.

The software reviewer raised a literal-metadata cache issue while the final author check was fixing it. The reviewer rechecked the current source and closed it: metadata/discovery are `no-store`; immutable data routes keep ETags.

Pre-fix checks: 613 Python unit tests passed; 11 Node 22.23.4 service tests passed with the Python-built dual-species bundle. Independent narrow checks of archived public EGFR, ERBB2 and PTPRC in both species verified all 112 native construct sequences against their declared reference slices. This is selected-record verification, not full-catalogue validation.

### Round 2

All four original seats completed a full rereview. The protein scientist passed the revised source/results: all four scientific R1 findings closed, 112 native exports and six references remained exact, and all six positive literature counts retained their frozen query links. The software and interface seats closed their R1 findings and identified the same additional export edge case. The security seat closed the R1 symlink finding and identified a permission-preservation issue.

| Finding | Priority | Decision and correction |
| --- | --- | --- |
| SEC-IMP-R2-1: atomic metadata replacement changes public files from `0644` to `0600` | P2 | Accept. Preserve the verified destination permissions through update and rollback; test both species and a non-default mode. No production file was changed. |
| SYS-I2-1 / interface: encoded FASTA exceeds its result cap and loses the fallback download | P2 | Accept. Keep the validated release/species/record context while checking the encoded representation; return the complete immutable download and selector on overflow. Test a schema-valid warning-heavy record and Unicode/reserved-character expansion. |

R2 checks included 52 focused Python tests and 16 Node 22 tests. These passing tests did not cover the two newly reproduced edge cases; targeted regressions are required.

### Full-catalogue author check

Aggregate inspection of the archived October 17 public release found 9,949 targets and 78,125 reported constructs. Its declared frozen source reports total 3,745,998,719 bytes. The initial whole-report map and 8 MiB flat-manifest limit were unsuitable for this scale. Source compilation now releases each raw report before reading the next, validator checks stream per target, and portable fingerprint/export checks retain their safeguards. The release-manifest budget is 96 MiB and the combined retained-index budget is 64 MiB; individual index, target and small-shard limits remain 32, 8 and 1 MiB. These are local engineering bounds.

The output provenance is also compact: the response identifies the complete input-hash map in the immutable manifest instead of repeating thousands of hashes in every result. A 10,001-input regression protects the ordinary response budget. Source/data hashes remain complete in the release artifact.

A local full-catalogue build completed on an isolated copy of the archived site. Its 226,326 members listed in the manifest total 2,087,991,808 bytes; the release manifest adds 79,027,981 bytes. Both indexes together occupy 6,238,429 bytes. Exact physical package size, Node startup, memory and retrieval checks follow below when measured. These local measurements do not establish the managed host's actual limits. The current GoDaddy account was inspected read-only; detailed account observations remain in ignored local evidence, while the tracked host record is public-safe.

### Round 3

The software and protein-science seats passed the revised source and local results. The interface and security seats closed earlier findings and identified two further boundaries:

| Finding | Priority | Decision and correction |
| --- | --- | --- |
| SEC-IMP-R3-1: a correctly checksummed, manifest-listed raw file can enter the bundle | P2 | Accept. Manifest membership is insufficient. Enforce the fixed approved download roles, media types, encodings and canonical paths; require exact reachable storage graphs and validate or regenerate untyped text/citation exports. Rehashed private database, disguised text and orphan-shard regressions must fail before packaging. |
| Interface: a valid long retirement history overwhelms release discovery | P2 | Accept. Keep bounded history in storage and literal retired-ID lookups; return a compact retirement count/lookup reference in discovery. |

The public metadata mode and encoded FASTA recovery fixes passed independent reproductions. No hosting or provider settings changed.

### User requirements added before round 4

The user requested interactive API/MCP usage limits and visible traffic accounting. The service now has five tools, explicit short queries, small result windows and shared admission quotas. Bulk analysis remains on the static website. The user chose the existing Plausible property plus complete service reports: sanitized logs cover all application requests, including bots and refusals; optional Plausible events supply eligible visitor geography with explicit filtering, billing and delivery coverage. Provider delivery is off by default. New local reports use private output files and expose their imported window and coverage.

### Round 4

All four original seats completed full rereviews including the new limits, proxy identity, analytics and report interfaces. They closed prior findings and required the following corrections:

| Finding | Priority | Decision and correction |
| --- | --- | --- |
| Scientist/software: exact sequence accession is missing from search | P2 | Accept. Include the recorded sequence accession and canonical isoform ID in shared search text; verify actual producer-to-HTTP/MCP queries for mixed mouse RefSeq provenance. |
| Software: SDK-handled MCP errors are recorded as successes | P2 | Accept. Observe bounded SDK-generated error responses without changing their wire format or recording their payloads. Missing arguments and unknown tools must appear as errors for both protocol versions, including HTTP 200 tool errors. |
| Interface: required search query is optional in OpenAPI | P2 | Accept. Mark `q` explicitly required and test the served GET/HEAD contracts. |
| Interface: report grouping limits lose known outcome counts in HTML | P2 | Accept. Preserve known counters and distinguish missing source dimensions from grouping loss. Unfiltered JSON, CSV and HTML totals must agree; filtered views must disclose incomplete attribution. |
| SEC-IMP-R4-1: report output links can overwrite unrelated files | P2 | Accept. Preflight all outputs and parent directories, refuse existing outputs, publish private files exclusively without following links, and preserve outside/input files on every failure path. |
| SEC-IMP-R4-2: collector redirects can forward client identity elsewhere | P2 | Accept. Reject redirects, perform no retry, and prove with localhost-only tests that no second collector receives identity or event data. |

Independent checks included 49 Node tests and 76 focused Python tests before these new regressions. The root full Python suite passed 655 tests. Only six intended human/mouse guide, plain-text guide and privacy golden hashes changed. Rendered report layout remains unverified because the browser explicitly blocked the file URL and prohibited workarounds; safer offline DOM interaction and source checks cover filtering, sorting and counts.

### Round 5

All four original seats reread the R4 corrections, packed transport, disease coverage and traffic reporting. Science and interface passed. Software and security independently reproduced one shared publication defect:

| Finding | Priority | Decision and correction |
| --- | --- | --- |
| SEC-IMP-R5-1 / SYS-R5-1: duplicate JSON keys hide a value from validation while its raw bytes enter the ZIP | P2 | Accept. Reject duplicate keys at every public JSON boundary, including manifests, schemas, logical shards, JSONL exports, packed members and runtime reads. Test rehashed hidden values, keys that decode to the same text, and nested duplicates. Preserve valid numeric spellings and exact original bytes. |

The disease projection follows the website's first 100 associations in source order, with explicit projected/source counts and original-source checksums. The full website TSV is preserved. The owner explicitly confirmed this scope on October 25, 2026. Target search excludes disease-ID filtering to avoid misleading membership results from this selection.

Physical packing preserves logical paths, checksums and canonical bytes. Related target shards share bounded gzip containers. The app keeps sealed metadata for the seven static bulk roles; staging checks their complete byte parity. Runtime does not claim to verify remote static copies.

### Round 6

All four original seats resumed and passed the corrected source and local behavior. Software and security independently rejected their original exploits at validation, packaging and Node startup, with generic errors and no emitted ZIP. Escaped and nested duplicate variants also failed. Valid numeric bytes, HTTP/MCP framing and scientific outputs stayed unchanged.

The final root Python suite passed 688 tests; the Node 22 suite passed 83 tests. Independent scientific checks retained exact agreement for 112 selected native constructs, six references and six frozen literature links. Fresh human and mouse build paths contained the shared guidance and privacy changes. Offline report interaction and count checks passed; rendered report layout remains unverified because of the browser's explicit file-access block.

Source review converged. Full-catalogue physical packaging and runtime measurements followed as separate gates. No production service, DNS or analytics settings were changed.

### Full-catalogue results

The complete source/static validation and physical transport checks passed. The upload gate then rejected the measured ZIP and published neither a ZIP nor an artifact record. A separate private runtime directory was generated from that already validated immutable source for local functional diagnostics; it is not an admissible upload artifact.

| Measurement | Result |
| --- | --- |
| Targets / constructs | 9,949 / 78,125, across both species |
| Logical files including manifests | 226,329 files / 2,167,020,277 bytes |
| Physical runtime data | 9,963 files / 316,899,235 bytes |
| Immutable static copies | 24 files / 1,148,459,779 bytes |
| Measured ZIP | 245,426,866 bytes: fails 80,000,000-byte local budget and hard limit below 100,000,000 bytes |
| Node 22.23.4 maintenance readiness / peak RSS | 59.901 seconds / 1,154,908,160 bytes |
| Node 22.23.4 active readiness / peak RSS | 60.929 seconds / 1,213,923,328 bytes |
| Functional checks | All five tools through HTTP and MCP in both modes, both species, 18 paced requests per mode, matching sequences and hashes, clean exits |
| Production dependency install and audit | Passed; zero reported vulnerabilities |

These local macOS runs kept default quotas, loopback binding, trusted forwarding disabled and Plausible delivery off. The profile contains one full release. Two full retained releases, production concurrency/latency, the selected host's runtime limits and report visual rendering remain unverified. Local diagnostics establish full-catalogue function; they do not establish public launch readiness.

Lossless Brotli feasibility was checked without changing source or dependencies. Thirty deterministic samples covered the target-container size deciles; the three largest containers were measured separately. Exact-byte round trips passed. Correctly weighted sampling estimated about 189 MB including measured metadata compression and calculated ZIP headers, before small app overhead. This is an estimate, not a complete Brotli archive measurement or a mathematical size bound. It does not support implementing that codec as a fix for the current upload gate.

Keep the full scientific contract and the current strict package limits. The self-contained GoDaddy managed upload path is blocked. A different hosting target or data-hosting design requires the separate scoped decision already defined in the plan; no provider, cost, access or account change was made.

### Round 7

All four original seats checked the full-catalogue evidence and revised deployment conclusions. They confirmed local functionality, unchanged source/data integrity and the blocked launch. The interface seat found three stale plan statements, also accepted by the software seat:

| Finding | Priority | Decision and correction |
| --- | --- | --- |
| Plan still promises empty-query catalogue access and bulk-download discovery | P2 | Accept. Require at least two trimmed literal query characters, describe the fixed 20-match window, and direct bulk analysis to the static website. Runtime and OpenAPI already enforce this policy. |

The opening design now also states the measured hosting blocker. These documentation corrections change no runtime, scientific data, limits or deployment settings. The same four seats resume for closure in round 8.

### Round 8 and conclusion

All four original seats passed the final source and documentation. The three stale contract statements are closed; no material implementation finding remains. Runtime and scientific data were unchanged, so the reviewers reused the demonstrated checks instead of repeating them for documentation edits.

The source is ready for a PR. Public launch remains blocked by the measured GoDaddy upload limit. Report visual QA and the remaining host, capacity, proxy, TLS, retention, analytics and recovery checks retain their stated scope. Review convergence authorizes no deployment or account change.

## Fixed-GoDaddy architecture revision

The owner confirmed that GoDaddy remains the provider and the plan must change. This supersedes the earlier recommendation to choose a different hosting target. The 245 MB ZIP and 1.21 GB startup peak above describe the previous all-data transport; they are not measurements of the revised reader.

### Round 9 plan review

The same four seats reviewed the small managed app plus immutable target assets on the existing GoDaddy static site. The plan keeps full summaries local for zero-origin search and pins target-container hashes for bounded reads. It preserves every scientific record, public schema, release ID, tool, quota and traffic/privacy rule. Full source validation remains in staging.

The prototype with exact local summary lookup files measured 6,470,661 bytes compressed and 41,144,863 bytes expanded. This establishes sizing feasibility only. Implementation, full current-plus-previous packaging, memory, cold/warm reads, failure recovery, retention and actual GoDaddy behavior still require verification.

| Finding | Priority | Decision and correction |
| --- | --- | --- |
| Inherited security wording still prohibits every runtime URL fetch | P2 | Accept. Permit only sealed fixed-origin static data reads and the separately controlled optional collector; continue to prohibit client-directed URLs, uploads and proxying. |
| Inherited isolation fallback suggests another provider | P2 | Accept. Require GoDaddy-supported isolation controls or explicit owner acceptance of the documented residual risk; otherwise keep launch blocked within the fixed provider decision. |

The related prompt-injection and integrity rows now distinguish user-directed fetching from sealed asset reads and local readiness from per-read remote validation. No runtime or host action occurred during this plan review. The same seats resume for round 10 closure before implementation.

### Round 10 plan closure

All four original seats passed the revised fixed-GoDaddy plan for implementation. The inherited fetch and provider wording is closed. The static sender must also restrict transfers to catalog-authorized assets while preserving the existing receiver asset union. No runtime implementation, host change or deployment is implied by plan approval. The same seats will review the resulting code and measured full-data artifacts.

### Round 11 implementation review

The same seats reviewed the implemented catalog, local summary reader, bounded static reader, packaging and additive publication path. Scientific and interface checks passed, including an independent replay of 291 legacy operations and exact agreement for all 118 selected FASTA exports. Source startup/search made zero target-origin reads. The full Python run had only the two intended human/mouse `.htaccess` golden changes; those hashes were reviewed and their three focused tests passed.

| Finding | Priority | Decision and correction |
| --- | --- | --- |
| SYS-I11-1: retirement archive metadata bypasses the public-text screen | P2 | Accept. Apply one shared Python guard to source roots and compact query roots, and the matching code-owned policy to all Node root formats and retired responses. Allow safe public HTTPS permalinks or null; reject private text, encoded markers, credentials, unsafe hosts, query/fragment and controls. Preserve the 1,024-character limit and fail rather than discard a tombstone. |
| SYS-I11-2: originless immutable responses omit `Vary: Origin` | P2 | Accept. Merge `Origin` into the response's cache variation even when no Origin header is present; preserve other tokens and wildcard semantics. Check both species, allowed origins, HEAD, 304 and denial responses. |

The inherited service README provider-choice sentence was also corrected during review. Runtime fixes passed 110 Node tests; producer/compiler checks passed with matching policy declarations. Scientific payloads, public schema bytes and release IDs stayed unchanged.

The normal two-full-release package passed mandatory source/static validation and size gates: 12,907,468 bytes ZIP, 82,192,952 bytes expanded, 35 files. Its 18 query-data files occupy 81,976,791 bytes; 19,898 authorized target assets occupy 447,159,606 bytes on the static side. This R11 archive predates the two fixes and is retained as a validated data baseline, not the final code release. Its ZIP SHA is `819a943fbbe31ca72660b67b5bd87047d901c6263d68972bd56602930e786c82`.

Software and security accepted the narrow code-only repack described in the plan for these fixes: independently pin that full-build baseline, validate the unchanged data with current policy, verify static bytes, and replace every old code file. Scientific validation is inherited explicitly. Default new-data builds retain complete validation. The two local static header files were refreshed separately with before/after evidence; scientific files and original artifacts stayed unchanged.

The same seats resume for round 12 to review the fixes and repack implementation. Final artifact/runtime measurements and live GoDaddy gates remain separate; no deployment or host change occurred.

### Round 12 verification and artifact

The corrected source passed the full 717-test Python suite and all 110 Node tests. Only the two intended human/mouse `.htaccess` golden hashes changed in this revision. The independent scientific replay remains exact for 291 operations and 118 FASTA exports; no scientific payload, public schema bytes or release ID changed.

The final code-only artifact is 12,909,859 bytes ZIP, 82,198,385 bytes expanded and 36 entries. Its SHA-256 is `2b1b4b1ab8d5066ab26614b51067e9fed4d24119fe4962cd87cd00a0f51005c5`. All 18 data entries and their seal are byte-identical to the independently pinned full-build baseline; all 17 code/license files match the current checkout. Production installation and audit passed with zero reported vulnerabilities. The receipt labels scientific validation as inherited and current metadata/static checks as rerun.

| Final local Node 22.23.4 check | Result |
| --- | --- |
| Active / maintenance reported app startup | 436.5 ms / 434.9 ms |
| Active / maintenance peak OS RSS | 256,163,840 / 257,327,104 bytes |
| Largest-target diagnostic peak OS RSS | 348,897,280 bytes |
| HTTP/MCP parity | All five tools, both species and both full releases, passed in both modes |
| API attempts / retries | 87 / 0, with default quotas and expected maintenance responses |
| Target assets at startup / metadata / search | 0 / 0 / 0 |
| Cold / immediate warm reads | One asset for the first target read, zero additional asset reads for warm calls |
| Local fixture reads / denials / failures | 9 / 0 / 0 |

The largest exercised container was 10,710,380 bytes decoded and 710,914 bytes compressed. All application children exited cleanly, the fixture stopped and input/code/data hashes remained unchanged. The first HTTP readiness observations include deliberate request pacing; startup values above are the app's own reported startup duration. Local transport, cache and latency checks do not establish live GoDaddy behavior. Real two-client concurrency was unavailable because a second loopback source address could not bind. External `ps` sampling was blocked; Node sampling and its OS peak counter remained available. Traffic-report visual QA remains unverified under the existing browser restriction.

All four original seats passed round 12 after independently closing the two findings and checking the code-only artifact. No material implementation finding remains. GoDaddy is the fixed provider, the measured app now fits its documented upload limit, and the existing PR can carry this implementation. No app, deployment, DNS, billing, access or live analytics change was made. Live GoDaddy checks require the separately authorized publication procedure.

### Round 13 scientist workflows and adversarial stress

The local stress test independently decoded every current-release target container: 5,328 human targets, 4,621 mouse targets and 78,125 constructs. All exported construct sequences matched their declared native reference intervals, lengths and SHA-256 values. Sequence accession/source, coordinates, warnings and evidence agreed between construct and sequence records. Construct pages had no duplicate identities or missing rows. All seven categories occurred in both species; human disease projections stayed within 100 rows and mouse disease remained inapplicable. These checks establish internal export consistency, not experimental validation or current external-database accuracy.

The reviewed R12 package then passed 19 workflow checks through actual loopback MCP with its production quotas: both retained releases and species, target search, every category/species pair, exact sequence retrieval, complete EGFR construct pagination, summary/evidence scope, 12 invalid tool requests and burst recovery. Unsupported sponsor and direct PDB-ID filters were rejected. Twelve simultaneous calls returned eight successes and four 429 responses with Retry-After; subsequent calls recovered. All 70 calls appeared in the application traffic report, including 17 expected errors and four admission refusals, with no malformed or missing application records in that run. Plausible export remained off.

The same review seats found three corrections:

| Finding | Priority | Correction and verification |
| --- | --- | --- |
| The pinned MCP SDK exposes cancellation at `context.mcpReq.signal`; `context.signal` was undefined | P2 | Pass the actual request signal. Before the fix, two disconnected clients could occupy both target-fetch slots until the independent three-second deadline. New current/legacy regressions require cancellation within one second, empty storage slots and a successful subsequent call. Independent checks observed cancellation and empty slots at 150 ms for both protocols. |
| Node's default 30-second connection-check interval delayed the configured six-second header timeout | P2 | Set the check interval to one second. A raw socket that sends incomplete headers every 500 ms now closes near six seconds; the regression requires a timeout response or connection reset before nine seconds and verifies subsequent health. |
| MCP guidance did not explain the meaning of a PDB construct sequence | P2 | State that exports use native reference intervals at mapped boundaries, and direct users to the PDB entry/study for deposited changes, tags, fusions and partners. Official current/legacy client tests verify this shared guidance. |

Independent malformed-message, prototype-name, body/header-size, HTTP framing, slow-body, pagination/cursor, concurrency, recovery, privacy and request-accounting probes used synthetic fixtures and loopback connections. Node parser failures and incomplete-header timeouts precede the application handler and produce no application ledger record. The analytics guide now states that host/edge security logs are required for this traffic.

Workflow limits remain explicit. PDB category filtering applies to one selected target; exact PDB accession matching requires filtering its paged results. Structure evidence combines PDB and predicted records, so a predicted model can require later pages and may have no downloadable asset. Literature provides recorded counts and query links, rather than paper-level construct provenance. Sponsor relationships and filters are absent from the current schema; adding them requires a later data/schema extension. Search results cannot establish absence of sponsorship.

After the fixes, all 113 Node tests passed with zero skips, including real Python-generated query, packed and overflow fixtures. The runtime-only R13 package is 12,909,969 bytes ZIP, 82,198,601 bytes expanded and 36 entries; its SHA-256 is `6165272a5270ecd0bc1d91364a9742df60c4d905f6cf81c9bab9aa3abdf8d4f3`. All 18 data entries and the data seal remain byte-identical to R12; packaged source matches the revised checkout. Production installation, a fresh dependency audit with zero reported vulnerabilities, and packaged readiness/tool-inventory/release smoke checks passed. Scientific validation is inherited, with the additional full current-release consistency scan above. The Python implementation was unchanged; the earlier 717-test result was not rerun for these JavaScript fixes.

All four original review seats passed after rechecking the fixes. No material finding remained in this bounded local scope. All stress traffic stayed local. Live GoDaddy proxy/TLS/WAF behavior, verified client identity/geography, host/edge log coverage, process/replica quota behavior and extended-duration load remain publication gates. No deployment or host settings changed.

### Round 14 human-facing MCP quick start

Both For AI agents guides now start with connection steps and scientist example questions. The connection card reads each portal's own metadata and shows a copyable HTTPS MCP address only for an active publication. Failed metadata reads and unpublished services retain website/download instructions; copy controls have manual-selection fallbacks. The shared Copy AI prompt describes all five tools, release pinning and one-sequence retrieval. The PDB native-interval caveat and current sponsor limitation remain explicit.

The same four seats reviewed the changes and rechecked corrections. The science seat found that a fixed portal species conflicted with opposite-species examples; the prompt now follows the question's explicit species and uses the portal species as its default. The interface seat found that the initial/no-JavaScript text presumed publication; it now describes publication conditionally and points to the native View prompt control. HTTPS host-case/default-port normalization is accepted while credentials, query strings, fragments, other schemes and unexpected paths stay hidden.

All four seats passed. All 718 Python tests passed, and browser regressions checked actual human/mouse builds at desktop and 390/320-pixel widths, active/unpublished/invalid metadata, normalized addresses, and clipboard success/failure. The golden diff contains exactly five intended artifact types per species: guide HTML/JavaScript, homepage, plain-text guide and CSS. No scientific data, backend service code, host settings or deployment changed. Live client setup and publication checks remain part of the deployment procedure.
