"""PR-3: labels and exports that did not match the data behind them.

Every assertion runs against **both** portal titles (AGENTS.md parity); the
fixes all live in the shared renderer.

* **F13B/OA-017** (133 human / 229 mouse) — the target-universe bucket is an
  inclusion rationale, not a finding about the analysed record. All four
  consumers of ``_design_region_label`` printed a bucket-derived heading over
  an empty Boundary/Label/Source. Only the multipass branch consulted the
  record.
* **F15** (724 / 407 pages, 1,556 / 983 rows) — a conserved motif was reported
  "detected" whenever GPCRdb had *annotated* every position, with no residue
  test, so a receptor whose 3.50 is not an arginine still read "DRY motif
  detected".
* **SF4a** (92 / 1,681) — blank identity and coverage cells with no reason.
* **SF5** — four identity numbers on one page, computed four ways, under a
  column headed just "Identity".
* **F24** (86 / 67) — an accessible-identity value on a page that states it has
  no extracellular region.
* **OA-002** (3,333 / 2,874) — design-region BLAST ranges printed in
  query-local coordinates, unlabelled.
* **OA-019** — the workbench picker offered GPCR fusion variants, then exported
  the unmodified native loop under the fusion's name.
* **PDB id collisions** (113 / 7) — one PDB chain backing two constructs got
  one id.
* **SF10** — the payload's ``entryName`` held an accession.
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

from agdesign2.gpcrdb import motif_status  # noqa: E402
from agdesign2.portal import (  # noqa: E402
    NO_DESIGN_REGION_LABEL,
    _build_index_entry,
    _design_region_label,
    _motif_status_from_stored,
    _construct_workbench_js,
    _render_cross_reactivity_rows,
    _soluble_construct_label,
    render_detail_page,
)

FIXTURES = ROOT / "tests" / "golden" / "fixtures"
PORTAL_TITLES = ("OpenAntigens", "OpenAntigens Mouse")


def _identity_cell(html: str) -> str:
    """The Homology table's Identity cell, not "somewhere on the page".

    The reason also reaches the Entry and accessible-identity cells, so a test
    that greps the whole page passes even with the identity fix reverted.
    """
    import re

    section = re.search(r'id="homology".*?</section>', html, re.S)
    assert section, "no Homology section"
    row = re.search(r"<tbody>(.*?)</tbody>", section.group(0), re.S)
    assert row, "no Homology rows"
    cells = re.findall(r"<td>(.*?)</td>", row.group(1), re.S)
    return cells[2].strip()


def _report(name: str = "sfta3_human") -> dict:
    return json.loads((FIXTURES / f"{name}_report.json").read_text(encoding="utf-8"))


def _render(report: dict, portal_title: str, entry_extra: dict | None = None) -> str:
    item = {
        "json_report": f"{report['target']['entry_name'].lower()}_report.json",
        "status": "ok",
        "query": report["target"].get("gene_symbol"),
    }
    entry = _build_index_entry(item, report, None)
    entry.update(entry_extra or {})
    with TemporaryDirectory() as root:
        return render_detail_page(
            entry, report, batch_dir=Path(root), portal_title=portal_title, sibling_link=None
        )


class DesignRegionLabelTests(unittest.TestCase):
    """F13B/OA-017."""

    def test_bucket_does_not_invent_a_region_the_record_does_not_have(self) -> None:
        report = _report()
        report["ectodomain"] = None
        for bucket in ("Single-pass", "Secreted", "GPI"):
            with self.subTest(bucket=bucket):
                entry = {"topology_bucket": bucket}
                self.assertEqual(_design_region_label(entry, report), NO_DESIGN_REGION_LABEL)
                self.assertEqual(_soluble_construct_label(entry, report), "Soluble Constructs")

    def test_the_multipass_branch_is_left_alone(self) -> None:
        # It already consulted the record, and "Full-Length Membrane Protein"
        # is a real design scope, not an empty one. 1,532 human / 1,373 mouse
        # pages depend on it; the plan says explicitly to leave it.
        report = _report()
        report["ectodomain"] = None
        entry = {"topology_bucket": "Multipass"}
        self.assertEqual(_design_region_label(entry, report), "Full-Length Membrane Protein")
        self.assertNotIn(NO_DESIGN_REGION_LABEL, _soluble_construct_label(entry, report))

    def test_bucket_still_names_the_region_when_the_record_has_one(self) -> None:
        report = _report("itb5_human")
        self.assertEqual(
            _design_region_label({"topology_bucket": "GPI"}, report),
            "GPI-Anchored Extracellular Region",
        )
        self.assertEqual(
            _design_region_label({"topology_bucket": "Secreted"}, report), "Secreted Region"
        )
        self.assertEqual(
            _design_region_label({"topology_bucket": "Single-pass"}, report),
            "Extracellular Region",
        )

    def test_all_four_consumers_agree_on_the_page(self) -> None:
        for title in PORTAL_TITLES:
            with self.subTest(portal_title=title):
                report = _report()
                report["ectodomain"] = None
                html = _render(report, title, {"topology_bucket": "Single-pass"})
                self.assertIn(f"<h2>{NO_DESIGN_REGION_LABEL}</h2>", html)
                self.assertNotIn("<h2>Extracellular Region</h2>", html)
                self.assertNotIn("Soluble Extracellular-region Constructs", html)
                # The Track row must say the bucket is the catalogue's.
                self.assertIn("target-universe inclusion rationale", html)
                self.assertIn("Single-pass", html)


class ConservedMotifStatusTests(unittest.TestCase):
    """F15 — position presence is not motif presence."""

    def test_matching_residues_are_matched(self) -> None:
        self.assertEqual(
            motif_status("DRY / E/DRY activation motif",
                         {"3.49": "D", "3.50": "R", "3.51": "Y"}, True),
            "matched",
        )

    def test_a_substituted_residue_is_named(self) -> None:
        self.assertEqual(
            motif_status("DRY / E/DRY activation motif",
                         {"3.49": "D", "3.50": "Q", "3.51": "Y"}, True),
            "variant at 3.50",
        )

    def test_the_x_position_is_unconstrained(self) -> None:
        # NPxxY: 7.51 and 7.52 are the "xx" and may be anything.
        self.assertEqual(
            motif_status("NPxxY / helix-7 motif",
                         {"7.49": "N", "7.50": "P", "7.51": "V", "7.52": "I", "7.53": "Y"}, True),
            "matched",
        )

    def test_e_is_accepted_at_3_49(self) -> None:
        # The motif is "E/DRY": glutamate is the documented alternative.
        self.assertEqual(
            motif_status("DRY / E/DRY activation motif",
                         {"3.49": "E", "3.50": "R", "3.51": "Y"}, True),
            "matched",
        )

    def test_stored_detected_status_is_recomputed_from_the_residues(self) -> None:
        stored = {
            "name": "DRY / E/DRY activation motif",
            "generic_numbers": ["3.49x49", "3.50x50", "3.51x51"],
            "sequence": "DQY",
            "status": "detected",
        }
        self.assertEqual(_motif_status_from_stored(stored), "variant at 3.50")

    def test_unrecognised_shape_keeps_the_stored_status(self) -> None:
        self.assertEqual(
            _motif_status_from_stored({"name": "x", "generic_numbers": [], "status": "partial"}),
            "partial",
        )


class HomologyCellTests(unittest.TestCase):
    """SF4a, SF5, F24."""

    def _with_rows(self, rows: list[dict], *, ectodomain: dict | None) -> dict:
        report = _report("itb5_human")
        report["ectodomain_homology"] = rows
        if ectodomain is None:
            report["ectodomain"] = None
        return report

    def test_blank_identity_gets_a_reason(self) -> None:
        for title in PORTAL_TITLES:
            with self.subTest(portal_title=title):
                # The mouse live path leaves identity None when the target has
                # no ectodomain (pipeline.py:1600); the cell was simply blank.
                report = self._with_rows(
                    [{"species": "mouse", "accession": "P00001", "identity": None,
                      "coverage": None}],
                    ectodomain=None,
                )
                html = _render(report, title)
                self.assertEqual(_identity_cell(html), "no extracellular region")

    def test_a_projected_row_with_no_number_reads_n_a(self) -> None:
        for title in PORTAL_TITLES:
            with self.subTest(portal_title=title):
                report = self._with_rows(
                    [{"species": "mouse", "accession": "P00001", "identity": None,
                      "coverage": None, "ectodomain_start": 24, "ectodomain_end": 719}],
                    ectodomain={"start": 1},
                )
                html = _render(report, title)
                self.assertEqual(_identity_cell(html), "n/a")

    def test_a_failed_boundary_projection_says_so(self) -> None:
        """The cause behind all 92 human pages the plan counted: the ortholog
        resolved and the target has a region, but the projection failed the
        0.70 coverage gate. They all read "n/a" before."""
        for title in PORTAL_TITLES:
            with self.subTest(portal_title=title):
                report = self._with_rows(
                    [{"species": "mouse", "accession": "P00001", "identity": None,
                      "coverage": None, "ectodomain_start": None, "ectodomain_end": None}],
                    ectodomain={"start": 1},
                )
                html = _render(report, title)
                self.assertEqual(
                    _identity_cell(html), "boundaries could not be mapped confidently"
                )

    def test_entryless_row_says_no_ortholog_resolved(self) -> None:
        for title in PORTAL_TITLES:
            with self.subTest(portal_title=title):
                report = self._with_rows(
                    [{"species": "macaca_fascicularis", "identity": None, "coverage": None}],
                    ectodomain={"start": 1},
                )
                html = _render(report, title)
                self.assertEqual(_identity_cell(html), "no ortholog resolved")

    def test_withheld_row_says_so(self) -> None:
        # The status exactly as ortholog_table.py writes it. It once rendered
        # "withheld: withheld: ..." because the renderer prefixed it again.
        for title in PORTAL_TITLES:
            with self.subTest(portal_title=title):
                report = self._with_rows(
                    [{"species": "macaca_fascicularis",
                      "selection_status": "withheld: only 13 aligned residues",
                      "identity": None, "coverage": None}],
                    ectodomain={"start": 1},
                )
                html = _render(report, title)
                self.assertEqual(_identity_cell(html), "withheld: only 13 aligned residues")

    def test_identity_header_names_the_reference(self) -> None:
        for title in PORTAL_TITLES:
            with self.subTest(portal_title=title):
                html = _render(_report("itb5_human"), title)
                self.assertIn("Identity to this human protein", html)
                self.assertIn("their numbers are not directly comparable", html)

    def test_no_region_means_no_accessible_identity(self) -> None:
        # F24: the value is meaningless where the page says there is no region.
        for title in PORTAL_TITLES:
            with self.subTest(portal_title=title):
                report = self._with_rows(
                    [{"species": "mouse", "accession": "P00001", "identity": 90.0,
                      "coverage": 95.0, "extracellular_surface_identity": 88.0,
                      "extracellular_surface_aligned_positions": 22}],
                    ectodomain=None,
                )
                html = _render(report, title)
                self.assertNotIn("88.00% (22 aa)", html)


class CrossReactivityCoordinateTests(unittest.TestCase):
    """OA-002."""

    HIT = {
        "subject_id": "sp|P00001|TEST_MOUSE",
        "query_alignment": "ACDEF",
        "subject_alignment": "ACDEF",
        "query_start": 1,
        "query_end": 5,
        "subject_start": 10,
        "subject_end": 14,
    }

    def test_region_local_coordinates_are_converted(self) -> None:
        with TemporaryDirectory() as tmp:
            html = _render_cross_reactivity_rows(
                [dict(self.HIT)], batch_dir=Path(tmp), query_origin=201
            )
        self.assertIn("201-205", html)
        self.assertNotIn("human 1-5 |", html)

    def test_full_length_search_is_unchanged(self) -> None:
        with TemporaryDirectory() as tmp:
            html = _render_cross_reactivity_rows([dict(self.HIT)], batch_dir=Path(tmp))
        self.assertIn("1-5", html)


class WorkbenchPickerTests(unittest.TestCase):
    """OA-019 — a fusion variant must never be exported as a native slice."""

    def test_the_picker_reads_only_the_native_construct_section(self) -> None:
        # The picker ships inside the structure widget, which needs a local
        # AlphaFold model; no committed fixture has one, so assert on the
        # emitted script itself.
        js = _construct_workbench_js()
        self.assertIn('const titles = new Set(["Precomputed Construct Summary"]);', js)
        self.assertNotIn("GPCR Expression and Stabilization Strategy", js)


class PdbConstructIdTests(unittest.TestCase):
    """113 human / 7 mouse pages gave two different sequences one name."""

    def test_two_spans_of_one_chain_get_different_ids(self) -> None:
        from agdesign2.config import AnalysisConfig
        from agdesign2.models import ExperimentalConstruct, Region
        from agdesign2.pipeline import AntigenAnalyzer

        analyzer = AntigenAnalyzer(config=AnalysisConfig(generate_assets=False))
        construct = ExperimentalConstruct(pdb_id="6OL9", method="X-ray", resolution="2.0 A", chains=[])
        first = Region(start=21, end=224, label="Chain A", source="PDB", metadata={"chain": "A"})
        second = Region(start=431, end=532, label="Chain A", source="PDB", metadata={"chain": "A"})
        self.assertNotEqual(
            analyzer._pdb_construct_name(construct, first),
            analyzer._pdb_construct_name(construct, second),
        )
        self.assertEqual(analyzer._pdb_construct_name(construct, first), "pdb_6ol9_a_21_224")


class EntryNameTests(unittest.TestCase):
    """SF10 — the ortholog table has an accession and a gene symbol, never a
    UniProt entry name, so entry_name=accession asserted something untrue.

    Nothing covered this: reverting all three sites left the whole suite green,
    golden manifest included, because the fixtures never reach the precomputed
    homology path.
    """

    def test_table_sourced_homology_rows_carry_no_entry_name(self) -> None:
        import re

        source = (ROOT / "src" / "agdesign2" / "precomputed.py").read_text(encoding="utf-8")
        constructions = re.findall(r"HomologyRecord\((.*?)\n\s*\)", source, re.S)
        self.assertTrue(constructions, "no HomologyRecord constructions found")
        for block in constructions:
            with self.subTest(block=block.strip().splitlines()[0]):
                self.assertIn("entry_name=None", block)

    def test_the_display_falls_back_to_the_accession(self) -> None:
        for title in PORTAL_TITLES:
            with self.subTest(portal_title=title):
                report = _report("itb5_human")
                report["ectodomain_homology"] = [
                    {"species": "mouse", "accession": "NP_001.1", "entry_name": None,
                     "identity": 90.0, "coverage": 95.0,
                     "ectodomain_start": 24, "ectodomain_end": 719}
                ]
                html = _render(report, title)
                self.assertIn("NP_001.1", html)

    def test_the_markdown_report_falls_back_too(self) -> None:
        """SF10 emptied the Markdown Entry column: _render_homology_table has
        no accession column to fall back to."""
        from agdesign2.models import HomologyRecord
        from agdesign2.reporting import _render_homology_table

        rows = _render_homology_table(
            [
                HomologyRecord(
                    species="mouse", accession="NP_001.1", entry_name=None,
                    ectodomain_start=1, ectodomain_end=10, available=True,
                )
            ]
        )
        self.assertIn("NP_001.1", "\n".join(rows))


class CrossReactivityProducerTests(unittest.TestCase):
    """OA-002's producer half.

    BLAST's query_alignment starts at hit.query_start WITHIN the queried
    sequence, not at the region origin, so the accessible-identity position set
    was shifted by query_start - 1 on 117,043 of the 148,231 stored values.
    """

    def _annotate(self, query_start: int, positions: set[int]):
        from agdesign2.config import AnalysisConfig
        from agdesign2.models import BlastHit
        from agdesign2.pipeline import AntigenAnalyzer

        analyzer = AntigenAnalyzer(config=AnalysisConfig(generate_assets=False))
        hit = BlastHit(
            subject_id="sp|P1|X_MOUSE", description="", species="mouse",
            identity=100.0, coverage=100.0, alignment_length=9,
            evalue=0.0, bitscore=10.0,
            query_start=query_start, query_end=query_start + 8,
            subject_start=1, subject_end=9,
            query_alignment="ACDEFGHIK", subject_alignment="ACDEFGHIW",
        )
        return analyzer._annotate_blast_hits_with_surface_identity(
            [hit], query_region_start=101, extracellular_surface_positions=positions
        )[0]

    def test_the_offset_inside_the_region_is_applied(self) -> None:
        # Region starts at 101; the alignment starts 4 residues into it, so the
        # first aligned residue is protein position 105.
        self.assertEqual(self._annotate(5, {105, 106, 107}).extracellular_surface_matches, 3)

    def test_an_alignment_at_the_region_start_is_unchanged(self) -> None:
        self.assertEqual(self._annotate(1, {101, 102, 103}).extracellular_surface_matches, 3)

    def test_the_old_unshifted_positions_no_longer_match(self) -> None:
        self.assertEqual(self._annotate(5, {101, 102, 103}).extracellular_surface_matches, 0)


if __name__ == "__main__":
    unittest.main()
