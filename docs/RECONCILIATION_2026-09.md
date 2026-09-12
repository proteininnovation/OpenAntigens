# September 2026 repository reconciliation

The official repository is `proteininnovation/OpenAntigens`. This revision starts from its public `main` commit `49dcc77` (July 27, 2026) and applies later development changes as patches. It preserves the official history.

| Original change date | Changes carried forward |
|---|---|
| August 25 | Serve compressed assets under their original filenames; verify deployment and cache purge. |
| September 5-6 | Validate PDB/PAE pairs, use global alignment consistently, prevent stale exports, preserve invalid calculator inputs, and validate complete builds and frontend assets. |
| September 11-12 | Add construct tabs, section navigation, concise selection and export guidance, canonical-sequence labels, assembly/exposure and PDB reminders, and the copy-control layout fix. |

Some changes were recommitted on September 11. Their original author dates and file diffs were used to distinguish new work from changes already in the July public release.

The reconciliation retains the official Complex Portal bulk lookups, mouse-specific guidance and provenance, privacy page and analytics restrictions, pinned viewer assets, and public release exclusions. AF3 remains disabled by default; public builds reject AF3-derived reports. Release validation also checks compressed text assets. Local model catalogs, manuscripts, production records and private development history are excluded.

The compression update removes the obsolete `--include-html` option. HTML remains plain; JavaScript, CSS and SVG retain their public filenames with gzip encoding configured through `.htaccess`.

Validation covers the complete Python suite, static JSON and SQLite builds, live serving, browser navigation, mutations, exports, the 3D viewer and mobile layouts. The review portal uses saved July analyses; rebuilding its pages does not rerun the scientific pipeline. Production deployment is separate from this PR.

Local checks passed: 328 pytest tests and 41 subtests, 327 unittest tests, and the Chrome regression suite. The saved-data preview contains 60 human and 53 matching mouse reports, with 1,580 construct cards and 104 structures. Report anchors, tab counts and local frontend assets were checked across all reports; human and mouse GPCR reports were exercised at 1920, 390 and 320 pixels.

The offline preview has no populated PubTator links and therefore does not pass the production release gate. Restore cached literature metadata or fetch it during the production build, then run the complete release validator before deployment. The preview passed the separate frontend-asset and public-content exclusion checks.

## Retained correctness and performance fixes

A second comparison against the latest development checkout confirmed these changes are retained:

| Area | Preserved behavior |
|---|---|
| Alignment | Shared global alignment for portal mappings and family identities, optional Parasail/Biopython acceleration, and a 1,024-entry cache keyed by sequences and backend. |
| Build efficiency | Reuse one analyzer per worker, write batch summaries at increasing intervals, and checkpoint ortholog rows individually. |
| Viewer efficiency | Reuse the PAE heatmap, schedule selection rendering through animation frames, and keep the existing builder state. |
| Compressed assets | Keep JS/CSS/SVG filenames; serve gzip bytes with the correct encoding and MIME headers in both portals, without URL rewriting or double compression. Writes remain atomic and repeatable. |
| Freshness and integrity | Preserve asset-content hashes, stale-file pruning, snapshot-specific sequence reads, complete-build checks, validated downloads and PDB/PAE pairs, and corrected UniProt pagination. |

The alignment, ortholog-checkpoint and family-identity functions match the development versions, as do the compared batch-worker, compression and portal-mapping functions. Focused verification passed 95 tests and 21 subtests; pytest emitted one cache-write warning under the local sandbox. No missing implementation was found in these areas.
