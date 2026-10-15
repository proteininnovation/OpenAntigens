"""Empty states must say what was checked, not print a bare negative.

Every assertion runs against **both** portal titles, per AGENTS.md: the fixes
live in the shared renderer and a species-specific regression would be a
divergence, not a feature.

Families covered (``outputs/page-review-2026-10-14/remediation/`` rows 17, 16,
20, 33, 40 and the furin follow-through):

* F03 — "The completed BLAST search returned no sequence similarity hits."
  printed beside a full-length hit table on 1,658 human / 1,598 mouse pages,
  because with no design region the region search never ran.
* F04a — "No InterPro annotations." beside a cited family accession on 2,205 /
  1,916 pages: the stored list is already filtered to domain-like entries
  inside the design region.
* F05 — "No assembly requirements reported." on all 4,702 mouse pages, 588 of
  them twins of a human page that lists curated complexes.
* F27 — bare "None" list placeholders on 2,729 / 2,746 pages.
* F31 — "No cysteine analysis available." with no reason on 465 / 711 pages.
"""

from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from agdesign2.portal import _build_index_entry, render_detail_page  # noqa: E402

FIXTURES = ROOT / "tests" / "golden" / "fixtures"
PORTAL_TITLES = ("OpenAntigens", "OpenAntigens Mouse")


def _render(report: dict, portal_title: str) -> str:
    item = {
        "json_report": f"{report['target']['entry_name'].lower()}_report.json",
        "status": "ok",
        "query": report["target"].get("gene_symbol"),
    }
    entry = _build_index_entry(item, report, None)
    with TemporaryDirectory() as root:
        return render_detail_page(
            entry,
            report,
            batch_dir=Path(root),
            portal_title=portal_title,
            sibling_link=None,
        )


def _base_report() -> dict:
    """A report with **no** design region — 10 of the 11 human fixtures are like
    this, which is why F03 and F04a reach thousands of pages."""
    return json.loads((FIXTURES / "sfta3_human_report.json").read_text(encoding="utf-8"))


def _region_report() -> dict:
    """A report that does have a design region, for the counter-cases."""
    report = json.loads((FIXTURES / "itb5_human_report.json").read_text(encoding="utf-8"))
    assert report.get("ectodomain"), "fixture must have a design region"
    return report


class BlastEmptyStateTests(unittest.TestCase):
    """F03."""

    def test_no_design_region_does_not_claim_a_completed_region_search(self) -> None:
        for title in PORTAL_TITLES:
            with self.subTest(portal_title=title):
                report = _base_report()
                report["ectodomain"] = None
                report["cross_reactivity_hits"] = []
                report["full_length_cross_reactivity_hits"] = [
                    {"accession": f"P0000{i}", "species": "mouse", "bit_score": 100 - i}
                    for i in range(3)
                ]
                html = _render(report, title)
                self.assertNotIn("The completed BLAST search returned no sequence similarity hits.", html)
                self.assertIn("no region-restricted BLAST search was run", html)
                # The sidebar card must count — and name — the search that ran.
                self.assertIn("<dt>Full-length BLAST sequence similarity hits</dt><dd>3</dd>", html)

    def test_design_region_with_a_true_zero_keeps_the_completed_negative(self) -> None:
        # Counter-case: where the search did run, the honest negative stays.
        for title in PORTAL_TITLES:
            with self.subTest(portal_title=title):
                report = _region_report()
                report["cross_reactivity_hits"] = []
                report["full_length_cross_reactivity_hits"] = []
                html = _render(report, title)
                self.assertIn("The completed BLAST search returned no sequence similarity hits.", html)

    def test_a_design_region_target_does_not_claim_a_completed_full_length_search(self) -> None:
        """The full-length search runs only for multipass targets and targets
        with no design region (pipeline.py:667, module_refresh.py:132).

        3,237 human / 2,712 mouse pages have a design region, are not
        multipass, and never ran it — and none of them carries a BLAST note, so
        the unavailable check cannot see it either.
        """
        for title in PORTAL_TITLES:
            with self.subTest(portal_title=title):
                report = _region_report()
                report["topology"] = {"topology_class": "single_pass"}
                report["cross_reactivity_hits"] = [{"accession": "P1"}]
                report["full_length_cross_reactivity_hits"] = []
                report["notes"] = []
                html = _render(report, title)
                self.assertNotIn("completed full-length BLAST search returned no", html)
                self.assertIn("No full-length BLAST search was run for this target", html)

    def test_a_multipass_target_with_no_hits_does_claim_a_completed_search(self) -> None:
        for title in PORTAL_TITLES:
            with self.subTest(portal_title=title):
                report = _region_report()
                report["topology"] = {"topology_class": "multipass_mixed"}
                report["cross_reactivity_hits"] = []
                report["full_length_cross_reactivity_hits"] = []
                report["notes"] = []
                html = _render(report, title)
                self.assertIn("completed full-length BLAST search returned no", html)

    def test_full_length_block_has_an_empty_state(self) -> None:
        for title in PORTAL_TITLES:
            with self.subTest(portal_title=title):
                report = _base_report()
                report["ectodomain"] = None
                report["cross_reactivity_hits"] = []
                report["full_length_cross_reactivity_hits"] = []
                html = _render(report, title)
                self.assertIn("completed full-length BLAST search returned no sequence similarity hits", html)


class InterProEmptyStateTests(unittest.TestCase):
    """F04a."""

    def test_table_scope_is_stated(self) -> None:
        for title in PORTAL_TITLES:
            with self.subTest(portal_title=title):
                html = _render(_base_report(), title)
                self.assertIn("lie entirely inside the design region", html)
                self.assertIn("Family-level entries are shown under Family Context", html)

    def test_no_design_region_says_so(self) -> None:
        for title in PORTAL_TITLES:
            with self.subTest(portal_title=title):
                report = _base_report()
                report["ectodomain"] = None
                report["interpro_annotations"] = []
                # The fixture itself carries an InterPro "lookup failed" note
                # (it is one of the 157 human F04b reports); clear it so this
                # case exercises the no-design-region branch.
                report["notes"] = []
                report["canonical_family"] = {"accession": "IPR005461", "name": "A family"}
                html = _render(report, title)
                self.assertNotIn("No InterPro annotations.", html)
                self.assertIn("no region-restricted domain entries could be listed", html)

    def test_failed_lookup_says_so(self) -> None:
        for title in PORTAL_TITLES:
            with self.subTest(portal_title=title):
                report = _base_report()
                report["interpro_annotations"] = []
                report["notes"] = [
                    {
                        "severity": "warning",
                        "source": "InterPro",
                        "message": "InterPro/Pfam annotation lookup failed: Expecting value",
                    }
                ]
                html = _render(report, title)
                self.assertNotIn("No InterPro annotations.", html)
                self.assertIn("InterPro lookup failed for this report", html)

    def test_region_present_but_nothing_inside_it(self) -> None:
        for title in PORTAL_TITLES:
            with self.subTest(portal_title=title):
                report = _region_report()
                report["interpro_annotations"] = []
                report["notes"] = []
                html = _render(report, title)
                self.assertIn("lies entirely inside the design region", html)


class AssemblyEmptyStateTests(unittest.TestCase):
    """F05."""

    def test_scope_sentence_replaces_the_bare_negative(self) -> None:
        for title in PORTAL_TITLES:
            with self.subTest(portal_title=title):
                report = _base_report()
                report["assembly_requirements"] = []
                html = _render(report, title)
                self.assertNotIn("No assembly requirements reported.", html)
                self.assertIn("only evaluated for canonical integrin chains", html)

    def test_curated_complexes_are_not_followed_by_an_empty_result_claim(self) -> None:
        # 630 human pages list Complex Portal records; only 26 of them have an
        # obligatory assembly row, so keying on assembly_requirements alone
        # left 604 showing the sentence under a populated table.
        report = _base_report()
        report["assembly_requirements"] = []
        report["complex_portal_complexes"] = [
            {"complex_id": "CPX-1", "name": "A complex", "members": ["X", "Y"]}
        ]
        html = _render(report, "OpenAntigens")
        self.assertNotIn("an empty result here does not rule out", html)

    def test_a_populated_table_is_not_followed_by_an_empty_result_claim(self) -> None:
        # 26 human pages do list an obligatory requirement; the scope sentence
        # contradicted the table immediately above it.
        for title in PORTAL_TITLES:
            with self.subTest(portal_title=title):
                report = _base_report()
                report["assembly_requirements"] = [
                    {"classification": "obligatory_partner_requirement", "obligatory": True,
                     "confidence": "high", "summary": "Heterodimer", "partners": ["ITGB3"]}
                ]
                html = _render(report, title)
                self.assertIn("ITGB3", html)
                self.assertNotIn("an empty result here does not rule out", html)

    def test_non_obligatory_annotations_are_declared_not_denied(self) -> None:
        # 1,771 mouse reports carry non-obligatory SUBUNIT rows the table does
        # not show; claiming nothing was evaluated would be false for them.
        for title in PORTAL_TITLES:
            with self.subTest(portal_title=title):
                report = _base_report()
                report["assembly_requirements"] = [
                    {"classification": "interaction_context", "obligatory": False,
                     "confidence": "medium", "summary": "Interacts with X", "partners": ["X"]},
                    {"classification": "heteromeric_assembly", "obligatory": False,
                     "confidence": "medium", "summary": "Heterodimer with Y", "partners": ["Y"]},
                ]
                html = _render(report, title)
                self.assertIn("2 non-obligatory interaction annotation(s) are recorded", html)
                self.assertNotIn("only evaluated for canonical integrin chains", html)

    def test_mouse_says_complex_portal_was_not_queried(self) -> None:
        report = _base_report()
        report["assembly_requirements"] = []
        mouse = _render(report, "OpenAntigens Mouse")
        human = _render(report, "OpenAntigens")
        self.assertIn("Complex Portal is not queried for mouse targets", mouse)
        self.assertNotIn("Complex Portal is not queried for mouse targets", human)


class CysteineEmptyStateTests(unittest.TestCase):
    """F31."""

    def test_no_region_and_not_multipass_says_no_structure_was_loaded(self) -> None:
        for title in PORTAL_TITLES:
            with self.subTest(portal_title=title):
                report = _base_report()
                report["ectodomain"] = None
                report["topology"] = {"topology_class": "secreted_soluble"}
                report["cysteine_analysis"] = []
                html = _render(report, title)
                self.assertNotIn("No cysteine analysis available.", html)
                self.assertIn("no structure was loaded and no cysteine analysis was run", html)

    def test_region_present_says_the_region_has_no_cysteines(self) -> None:
        for title in PORTAL_TITLES:
            with self.subTest(portal_title=title):
                report = _region_report()
                report["cysteine_analysis"] = []
                html = _render(report, title)
                self.assertIn("The analyzed region contains no cysteines.", html)

    def test_structure_skipped_is_not_reported_as_no_cysteines(self) -> None:
        # The note the pipeline actually writes (pipeline.py:419) says
        # "skipped", not "failed". Reports in this state DO have cysteines in
        # the region — celr1_human has 81 — so the third branch would be false.
        for title in PORTAL_TITLES:
            with self.subTest(portal_title=title):
                report = _region_report()
                report["cysteine_analysis"] = []
                report["notes"] = [
                    {
                        "severity": "warning",
                        "source": "AlphaFold",
                        "message": (
                            "Structure-driven analyses were skipped: No AlphaFold model "
                            "matches the canonical UniProt sequence for Q9NYQ6."
                        ),
                    }
                ]
                html = _render(report, title)
                self.assertNotIn("The analyzed region contains no cysteines.", html)
                self.assertIn("Structure-driven analyses were skipped", html)

    def test_structure_failure_says_so(self) -> None:
        for title in PORTAL_TITLES:
            with self.subTest(portal_title=title):
                report = _region_report()
                report["cysteine_analysis"] = []
                report["notes"] = [
                    {
                        "severity": "warning",
                        "source": "AlphaFold",
                        "message": "Structure download failed for this target",
                    }
                ]
                html = _render(report, title)
                self.assertIn("Structure-driven analyses were skipped", html)


class FurinEmptyStateTests(unittest.TestCase):
    """The furin scan covers the design region only.

    ``_find_furin_sites`` returns [] immediately when there is no ectodomain
    (``pipeline.py:3865-3867``), so the predicate is the design region — the
    same one F03 uses — not the presence of a sequence. Every report in both
    snapshots has a sequence, so a sequence-based branch would be dead code.
    """

    def test_no_design_region_says_no_scan_was_run(self) -> None:
        for title in PORTAL_TITLES:
            with self.subTest(portal_title=title):
                report = _base_report()
                report["ectodomain"] = None
                report["furin_sites"] = []
                html = _render(report, title)
                self.assertNotIn("No candidate furin-like motifs found", html)
                self.assertIn("no furin-like motif search was run", html)

    def test_design_region_names_the_region_not_the_sequence(self) -> None:
        for title in PORTAL_TITLES:
            with self.subTest(portal_title=title):
                report = _region_report()
                report["furin_sites"] = []
                html = _render(report, title)
                self.assertIn("No candidate furin-like motif was found in the design region.", html)
                # The markdown report says "in the ectodomain"; the HTML must not
                # claim a wider scope than the scan actually has.
                self.assertNotIn("in the analyzed sequence", html)


class CrossReactivityCardLabelTests(unittest.TestCase):
    """The card's number must name which search it counts.

    ``_index_cross_reactivity_count`` returns the full-length count for
    multipass targets even when they have a design region and region hits, so
    an unlabelled number contradicted the region table below it on 263 human /
    244 mouse pages.
    """

    def test_multipass_card_says_full_length(self) -> None:
        for title in PORTAL_TITLES:
            with self.subTest(portal_title=title):
                report = _region_report()
                report["topology"] = {"topology_class": "multipass_mixed"}
                report["cross_reactivity_hits"] = [{"accession": "P1"}]
                report["full_length_cross_reactivity_hits"] = [{"accession": f"P{i}"} for i in range(6)]
                html = _render(report, title)
                self.assertIn("<dt>Full-length BLAST sequence similarity hits</dt><dd>6</dd>", html)

    def test_single_pass_card_names_the_design_region(self) -> None:
        for title in PORTAL_TITLES:
            with self.subTest(portal_title=title):
                report = _region_report()
                report["cross_reactivity_hits"] = [{"accession": f"P{i}"} for i in range(4)]
                report["full_length_cross_reactivity_hits"] = [{"accession": f"Q{i}"} for i in range(9)]
                html = _render(report, title)
                self.assertIn("BLAST sequence similarity hits</dt><dd>4</dd>", html)
                self.assertNotIn("<dt>Full-length BLAST sequence similarity hits</dt>", html)


class NonePlaceholderTests(unittest.TestCase):
    """F27 — a bare "None" reads as a finding; say whose checks found nothing."""

    def test_every_placeholder_site_is_covered(self) -> None:
        """Render a construct card whose evidence, warnings, PTM, metrics and
        GPCRdb-segment lists are all empty at once.

        The fixture-driven test below does not reach three of these sites, so
        reverting any one of them to the bare "None" left it green.
        """
        from tempfile import TemporaryDirectory

        from agdesign2.portal import _render_construct_card

        construct = {
            "name": "full_ectodomain",
            "start": 25,
            "end": 180,
            "sequence": "A" * 156,
            "evidence": [],
            "warnings": [],
            "ptms": [],
            "structural_metrics": {},
            "homologs": [],
            "gpcr_segments": [],
            "gpcr_generic_range": {"start": "5.50", "end": "6.50"},
        }
        target = {"entry_name": "IL7_HUMAN", "gene_symbol": "IL7",
                  "accession": "P10168", "taxon_id": "9606"}
        with TemporaryDirectory() as tmp:
            html = _render_construct_card(
                construct, target=target, batch_dir=Path(tmp), page_dir=Path(tmp),
                furin_sites=[], ptms=[],
            )
        self.assertNotIn("<li>None</li>", html)
        self.assertIn("None recorded by this report", html)          # evidence + warnings
        self.assertIn("No annotated PTM lies inside this construct.", html)
        self.assertIn("No structural metric was computed for this construct.", html)
        self.assertIn("No GPCRdb segment overlaps this construct.", html)

    def test_no_bare_none_list_item_survives(self) -> None:
        for title in PORTAL_TITLES:
            with self.subTest(portal_title=title):
                for name in ("itb5_human", "ghrhr_human", "vdac1_human"):
                    report = json.loads(
                        (FIXTURES / f"{name}_report.json").read_text(encoding="utf-8")
                    )
                    html = _render(report, title)
                    self.assertNotIn("<li>None</li>", html, name)


if __name__ == "__main__":
    unittest.main()
