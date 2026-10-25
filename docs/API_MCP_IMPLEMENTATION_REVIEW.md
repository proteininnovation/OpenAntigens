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

Fixes and targeted regressions are in progress. The same four seats will perform a full rereview before the PR.
