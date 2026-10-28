# API/MCP planning review

## Scope and completion rule

Review the design in `API_MCP_PLAN.md` against canonical source baseline `8008192`, the October 17 hosting screenshot, and inspected official GoDaddy/MCP documentation. The requested deliverable is a reviewed plan. Implementation, penetration testing, hosting feasibility tests, public launch, and production approval are separate steps.

Use independent security and systems reviewers. Each reports an issue ID, severity, failure scenario, source location, and minimum corrective action. Revise the plan, then have both reviewers reread the full revised design and the resolution record. Continue until both report no unresolved Critical, High, or Medium plan defects and no new material defects in a complete review round. Resolve Low findings or record a specific reason for accepting them.

An empirical hosting unknown is acceptable in a plan only when it has a named measurement/check, a safe stop rule, and an explicit owner decision if it fails. It remains a launch prerequisite. Planning convergence establishes agreement about the design, not evidence that the host or future code satisfies it.

## Evidence preparation

- Read the shared static build, compact downloads, construct export/mapping helpers, mouse overrides, development/SQLite APIs, public release packaging, release validators, static deployer, and security/deployment tests.
- Independent research covered reusable functions/data boundaries, the current official GoDaddy/Node deployment model and MCP SDK/protocol, and the threat model.
- Before adversarial review, corrected a release-ID/construct-ID hash dependency cycle, separated canonical payload hashes from final download-byte checksums, specified retention of public downloads outside app runtime storage, added shared-account isolation and incident response, and recorded current Node/MCP evidence.
- Existing local untracked target-curation material in the base checkout was preserved. Planning work is isolated on `docs/api-mcp-plan` in the attached worktree.

## Review rounds

### Round 1: independent adversarial review

Security and systems reviewers each read the complete plan and checked relevant repository/official protocol material. A separate source specialist checked projection/export semantics against the actual builders. No Critical or High issue was reported. The author triaged the three material source findings as Medium.

| ID | Severity | Failure scenario | Plan correction |
| --- | --- | --- | --- |
| SEC-R1-1; SYS-R1-1 | Medium | A single 300 KiB evidence record fits storage limits but exceeds the 256 KiB response cap; index/construct downloads cannot supply it. | Require complete per-species `records.jsonl.gz`, with all projected target/evidence/reference/construct records. Overflow errors identify the exact record and final immutable artifact/hash. Validate decompressed completeness and shard equivalence. |
| SYS-R1-2; DATA-R1-2 | Medium | An implementer cannot derive the promised section/class/filter tokens; portal tab category and stored classification differ, which also changes construct identity. | Define exact section/filter enums and arguments. Preserve raw classification separately; use the shared seven-tab category mapping for filtering and construct identity. Freeze an identity-tuple fixture. |
| DATA-R1-1 | Medium | A mouse record can retain a UniProt target identifier while exporting a RefSeq fallback sequence; one accession field would mislabel that sequence. | Separate identifier accession/source from sequence accession/source/hash, preserve mapping provenance, and leave an unrecorded sequence accession null. |
| DATA-R1-3 | Medium | Reused homolog loaders search cwd caches; identical frozen reports can yield different exports on another machine. | Require snapshot-scoped caches/reference/mapping rows, no network or ambient fallback, declared input hashes, and poisoned-cwd/network-disabled reproducibility tests. |
| SYS-R1-3 | Low | MCP clients cannot enumerate retained releases even though HTTP clients can. | `get_release` includes current ID and retained release IDs/metadata, as well as the selected release. |

The author also clarified that catalogue track is an inclusion rationale distinct from analyzed topology. An unavailable analyzed field cannot inherit the catalogue bucket as a scientific result.

### Round 2: full rereview

Both independent adversarial reviewers reread the complete revised plan. The source specialist performed a second data-contract/source-semantics review.

| Review | Result | Limits |
| --- | --- | --- |
| Independent security, round 2 | SEC-R1-1 closed. No unresolved or new material findings after checking every plan section and the interactions among bounds, public fields, hashes, retention, privacy, maintenance and rollback. | No runtime, host-control, dependency-installation or real-payload testing. |
| Independent systems, round 2 | All SYS-R1 findings closed. No new material finding; release identity, complete overflow fallback, exact filters/categories, discovery and publication gates are coherent. | One Low editorial issue: the section-route enum reference pointed to section 3. Corrected to section 5. No deployment or client tests. |
| Source specialist, round 2 | All DATA-R1 findings closed. Identifier/sequence provenance, seven-tab category mapping, snapshot cache isolation, exact filters and hash construction agree with the inspected source. | No scientific reanalysis, real-release build or sequence validation. |

**Convergence reached after two adversarial review rounds.** No known Critical, High or Medium plan defect remains. The Low discovery finding and Low cross-reference issue are resolved. The only plan change after the second full rereads was that cross-reference correction; it does not change the design.

## Required work after planning

The implementation and launch gates in the plan remain open: signed-in hosting availability, selected-release size/memory/latency, actual proxy/MCP-client behavior, process isolation, log/rate controls, named operational ownership, and publication/rollback rehearsal. They have checks and safe stop rules; review convergence does not close them.

No service code, account setting, app, credential, DNS entry, scientific data, or production artifact was changed. No push, merge, app creation, paid action, or deployment was performed. The deliverable is the local worktree plan, this review record, and the README discovery link.
