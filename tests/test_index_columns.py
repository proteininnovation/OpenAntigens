"""PR-10: index values that contradicted their own pages (plan row 41).

* ``has_assembly_requirements`` was true on 1,771 mouse rows against 26 human,
  because non-obligatory UniProt SUBUNIT rows counted while the page shows
  obligatory requirements only.
* ``top_disease_name`` was an Open Targets therapeutic-area node on 5,189 of
  5,328 human rows: indirect scores propagate up the ontology, so a root area
  always outscores the diseases under it.

Deliberately NOT changed: the eight mouse/source columns that are blank on
every human row. ``downloads.html`` documents them as blank there, and plan §7
treats them as expected. ``has_assets = 0`` is also correct (construct-image
rendering is disabled at ``cli.py:656``) and must not be "fixed".
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from agdesign2.open_targets import _top_specific_association  # noqa: E402


class TopDiseaseTests(unittest.TestCase):
    ROWS = [
        {"disease_id": "MONDO_0045024", "disease_name": "cancer or benign tumor", "score": 0.9},
        {"disease_id": "EFO_0000616", "disease_name": "neoplasm", "score": 0.85},
        {"disease_id": "EFO_0000305", "disease_name": "breast carcinoma", "score": 0.8},
    ]

    def test_a_therapeutic_area_is_not_the_top_disease(self) -> None:
        top = _top_specific_association(self.ROWS, {"MONDO_0045024", "EFO_0000616"})
        self.assertEqual(top["disease_name"], "breast carcinoma")

    def test_without_area_metadata_the_plain_top_row_is_kept(self) -> None:
        self.assertEqual(_top_specific_association(self.ROWS, set())["disease_name"],
                         "cancer or benign tumor")

    def test_all_areas_falls_back_rather_than_blanking(self) -> None:
        rows = self.ROWS[:2]
        self.assertEqual(
            _top_specific_association(rows, {"MONDO_0045024", "EFO_0000616"})["disease_name"],
            "cancer or benign tumor",
        )

    def test_no_associations(self) -> None:
        self.assertEqual(_top_specific_association([], {"X"}), {})


class AssemblyFlagTests(unittest.TestCase):
    def _flag(self, requirements: list[dict]) -> bool:
        import json

        from agdesign2.portal import _build_index_entry

        report = json.loads(
            (ROOT / "tests" / "golden" / "fixtures" / "sfta3_human_report.json").read_text(encoding="utf-8")
        )
        report["assembly_requirements"] = requirements
        entry = _build_index_entry({"json_report": "x.json", "status": "ok"}, report, None)
        return bool(entry["has_assembly_requirements"])

    def test_non_obligatory_rows_do_not_set_the_flag(self) -> None:
        self.assertFalse(self._flag([{"classification": "interaction_context", "obligatory": False}]))

    def test_an_obligatory_row_does(self) -> None:
        self.assertTrue(self._flag([{"classification": "obligatory_partner_requirement", "obligatory": True}]))

    def test_none(self) -> None:
        self.assertFalse(self._flag([]))


if __name__ == "__main__":
    unittest.main()
