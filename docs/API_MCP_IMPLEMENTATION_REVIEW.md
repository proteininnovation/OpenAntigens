# API/MCP implementation review

The implementation follows `API_MCP_PLAN.md` in the `feat/public-api-mcp` worktree branch. The publication target is a PR in `proteininnovation/OpenAntigens`; production deployment and account changes are separate authorized steps.

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

Aggregate inspection of the archived October 17 public release found 9,949 targets and 78,125 reported constructs. Its declared frozen source reports total 3,745,998,719 bytes. The initial whole-report map and 8 MiB flat-manifest limit were unsuitable for this scale. Source compilation now releases each raw report before reading the next, validator checks stream per target, and portable fingerprint/export checks retain their safeguards. The release-manifest and combined retained-index budgets are separate 64 MiB engineering trial limits; individual index, target and small-shard limits remain 32, 8 and 1 MiB.

The output provenance is also compact: the response identifies the complete input-hash map in the immutable manifest instead of repeating thousands of hashes in every result. A 10,001-input regression protects the ordinary response budget. Source/data hashes remain complete in the release artifact.

A local full-catalogue build on an isolated copy of the archived site is in progress. Exact package size, Node startup, memory and retrieval results must be recorded before claiming full-catalogue fit. These local measurements will not establish the managed host's actual limits. The current GoDaddy account was inspected read-only; detailed account observations remain in ignored local evidence, while the tracked host record is public-safe.

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

The accepted R4 corrections are in progress. The same four seats will reread the final source and outputs. The full-catalogue disease-feed scope and package fit remain separate unresolved gates; no full-catalogue deployment claim is made from small fixtures.
