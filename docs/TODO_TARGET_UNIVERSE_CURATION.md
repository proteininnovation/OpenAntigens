# TODO — target-universe curation sprint (deferred, opened 2026-10-15)

**Status: deferred by decision.** The 2026-10-14 page review surfaced targets whose page describes
them as surface or secreted antigens although their dominant annotation is intracellular. Membership
in the target universe is **out of scope for the current remediation work** — some of these proteins
may legitimately reach the cell surface, and deciding that is curation, not a bug fix. This file
records what was measured so a later sprint can start from evidence instead of repeating the review.

Nothing here blocks the remediation PRs. The one related item that *is* in scope is the rendering bug
where a catalogue bucket makes the page print an empty "Extracellular Region" heading for a record
with no ectodomain (133 human / 229 mouse pages) — that is a display defect, not a membership
question, and it is fixed in the remediation plan.

## What was measured

Evidence: `outputs/page-review-2026-10-14/remediation/` — plan §1 row 8, `work/trace_topology.md`,
`redteam/round2_counts.md`, `work/trace_f13.py`, `work/rt2_f13.py`.

| Observation | Size | Note |
|---|---|---|
| Reports that ship a design region and constructs although the record has no SIGNAL, no TM, no GPI and no extracellular TOPO_DOM | **263 pages** (131 human / 132 mouse) | the honest data definition; a reviewer keyword count gave 533, which is not reproducible from data |
| Of those, leaderless members of the curated "Secreted" catalogue, where the pipeline invents the region (`pipeline.py:894-926`) | 146 human / 132 mouse (130 human ship a construct) | includes real leaderless cytokines (IL-1β, FGF2 style) that belong in the universe |
| Universe rows whose UniProt entry no longer carries a "Secreted" location (release 2026_03) | 9 of the 146 leaderless; 50 of 1,670 Secreted-bucket rows | e.g. ABHD18 (Q0P651) is now annotated Mitochondrion only |
| Universe rows with "pseudogene" in the protein name | 20 | e.g. ANXA2P2 |
| Membrane-bucket rows resting on SURFY alone, with no non-isoform UniProt surface annotation | ~140 (400 rows total lack one; 260 are isoform-only) | for 72, SURFY's own subcellular field names only an organelle (ER/Golgi membrane) |
| Reports whose own `family_name` says "Intracellular" while the topology bucket is Secreted or membrane | 9 (6 ship a region) | cheap internal inconsistency signal; too small to be the remedy |
| Mouse-vs-human topology disagreements (mouse UniProt calls a region extracellular where the human ortholog calls it cytoplasmic) | 6, positional proxy | e.g. Atp11c |

The classifier that builds the universe lives in the separate **AKB2** repository
(`AKB2/docs/ACCESSIBLE_PROTEOME_METHODS.md`), not here. The rule is: a non-isoform UniProt location of
exactly "Secreted" → `secreted_soluble`. It does not require a signal peptide, which is deliberate.

## Suggested plan for the sprint

1. **Re-run the AKB2 classifier against the current UniProt release.** This alone removes the 9 stale
   leaderless rows and the 50 stale Secreted rows. Expect the 5,328-protein denominator to move, which
   changes every page count, the index and the manuscript's numbers — so do it as its own release.
2. **Decide the pseudogene rule** (20 rows). Either exclude, or keep with an explicit "putative
   pseudogene product" label on the page.
3. **Decide the SURFY-only rule** (~140 rows, 72 with organelle-only subcellular text). Options: require
   a second line of evidence; keep with a visible caveat; or keep only where an HPA/Uhlen secretome or
   surfaceome source agrees.
4. **Carry the inclusion evidence onto the page.** `target_sets.py:41-48` writes six columns and drops
   `classification_evidence`, `supporting_evidence_sources`, `nonisoform_*_locations` and
   `surfy_uniprot_subcellular`, so the portal cannot show *why* a protein is in the universe. Adding an
   "Inclusion evidence" row would let a reader judge the borderline cases instead of trusting the
   bucket — this is the highest-value item for a scientist and is independent of any membership change.
5. **Hand-curate the leaderless set** (49 human rows rest on UniProt alone with no HPA/secretome
   support). Non-classical secretion is real; the question is per protein.
6. **Add the internal consistency check** (`family_name` vs topology bucket) as a release assertion so
   new contradictions are caught at build time, whatever the curation decides.

## Related, already in the remediation plan (do not duplicate here)

- `pipeline.py:894-926` invents a design region for leaderless secreted targets and labels it
  "Complete extracellular domain based on topology annotations" (`pipeline.py:2159`) regardless of
  source. Honest labelling was **deferred with this file** — it describes a membership judgement, so it
  waits for the curation decision.
- `topology.py:23` counts every UniProt TRANSMEM without reading its description, so β-hairpin and
  β-barrel membrane proteins (MACPF: C6, C8A, C8B, C9, PRF1, MPEG1; gasdermins; VDAC1 — 19 reports) are
  treated as α-helical multipass and lose their ectodomain. **That one is a bug and is in the plan**
  (PR-7); it is not a membership question.
