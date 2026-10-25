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

The disease projection follows the website's first 100 associations in source order, with explicit projected/source counts and original-source checksums. The full website TSV is preserved. This is an implementation scope decision; no explicit user approval of the 100-row policy is inferred. Target search excludes disease-ID filtering to avoid misleading membership results from this selection.

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
