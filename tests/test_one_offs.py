"""PR-11: one-off corrections and the support-page findings.

* **mouse/llms.txt** was byte-identical to the human file, so the file an agent
  fetches first sent it to the human index, manifest, metadata, Methods and
  help — and the agent prompt on the mouse portal named the human llms.txt.
* **F28** (33 human / 9 mouse pages) — UniProt SITE "Not glycosylated" records
  that a sequon is NOT modified; a keyword test rendered it as a PTM.
* **F30** (30 / 14) — a glycosylation or disulfide inside an annotated
  transmembrane span, shown with no comment.
* **Row 46** — constructs up to 14,451 aa with no size warning.
* Stale text: help.html and constructs.html described the identity metric and
  the missing-homolog case as they were before PR-4 and PR-8.
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from agdesign2.portal import (  # noqa: E402
    CONSTRUCT_SIZE_WARNING_AA,
    _agent_guidance_panel,
    _construct_size_warning,
    _is_ptm_feature_dict,
    _ptm_membrane_conflict,
    _report_ptms,
    mouse_llms_txt,
)
from agdesign2.site_docs import agent_guide_markdown  # noqa: E402


class MouseLlmsTxtTests(unittest.TestCase):
    PER_PORTAL = (
        "downloads/agdesign2_portal_index.json",
        "downloads/agdesign2_portal_index.tsv",
        "downloads/download_manifest.json",
        "portal_metadata.json",
        "methods.html",
    )

    def test_the_mouse_file_is_not_the_human_file(self) -> None:
        self.assertNotEqual(mouse_llms_txt(agent_guide_markdown()), agent_guide_markdown())

    def test_per_portal_artifacts_point_at_mouse(self) -> None:
        mouse = mouse_llms_txt(agent_guide_markdown())
        for artifact in self.PER_PORTAL:
            with self.subTest(artifact=artifact):
                if f"https://openantigens.org/{artifact}" in agent_guide_markdown():
                    self.assertIn(f"https://openantigens.org/mouse/{artifact}", mouse)
                    self.assertNotIn(f"(https://openantigens.org/{artifact})", mouse)

    def test_the_shared_citation_files_stay_at_the_root(self) -> None:
        mouse = mouse_llms_txt(agent_guide_markdown())
        self.assertIn("https://openantigens.org/downloads/openantigens.bib", mouse)

    def test_the_release_gate_accepts_the_per_portal_file(self) -> None:
        # The gate used to require the two files to be byte-identical, which is the
        # defect this class describes; a rebuild failed validate-release-data on it.
        mouse = mouse_llms_txt(agent_guide_markdown())
        self.assertNotEqual(mouse, agent_guide_markdown())
        self.assertIn("https://openantigens.org/mouse/downloads/agdesign2_portal_index.json", mouse)

    def test_the_mouse_agent_prompt_names_the_mouse_llms_txt(self) -> None:
        mouse = _agent_guidance_panel(portal_title="OpenAntigens Mouse")
        human = _agent_guidance_panel()
        self.assertIn("https://openantigens.org/mouse/llms.txt", mouse)
        self.assertIn("read https://openantigens.org/llms.txt", human)


class NegativeSiteTests(unittest.TestCase):
    """F28."""

    def test_a_negative_annotation_is_not_a_ptm(self) -> None:
        self.assertFalse(_is_ptm_feature_dict({"type": "SITE", "description": "Not glycosylated"}))
        self.assertFalse(_is_ptm_feature_dict({"type": "SITE", "description": "Not N-glycosylated"}))

    def test_a_real_glycosylation_site_still_is(self) -> None:
        self.assertTrue(_is_ptm_feature_dict({"type": "SITE", "description": "Glycosylated in variant X"}))

    def test_already_built_reports_drop_it_on_render(self) -> None:
        report = {"ptms": [
            {"type": "SITE", "description": "Not glycosylated", "start": 5, "end": 5},
            {"type": "CARBOHYD", "description": "N-linked", "start": 9, "end": 9},
        ]}
        self.assertEqual([p["type"] for p in _report_ptms(report)], ["CARBOHYD"])


class MembranePtmTests(unittest.TestCase):
    """F30."""

    SPANS = [(39, 59)]

    def test_a_glycosylation_in_a_helix_is_flagged(self) -> None:
        self.assertIn("transmembrane span (39-59)",
                      _ptm_membrane_conflict({"type": "CARBOHYD", "start": 45}, self.SPANS))

    def test_one_outside_is_not(self) -> None:
        self.assertEqual(_ptm_membrane_conflict({"type": "CARBOHYD", "start": 70}, self.SPANS), "")

    def test_a_phosphorylation_is_not_judged(self) -> None:
        # Only modifications that cannot occur in a bilayer are flagged.
        self.assertEqual(_ptm_membrane_conflict({"type": "MOD_RES", "start": 45}, self.SPANS), "")


class ConstructSizeTests(unittest.TestCase):
    """Row 46 — MUC16's full ectodomain is 14,451 aa."""

    def test_a_very_large_construct_is_warned(self) -> None:
        self.assertIn("14,451 aa", _construct_size_warning({"start": 1, "end": 14451}))

    def test_an_ordinary_construct_is_not(self) -> None:
        self.assertEqual(_construct_size_warning({"start": 1, "end": CONSTRUCT_SIZE_WARNING_AA}), "")


class SupportPageTextTests(unittest.TestCase):
    def test_constructs_page_describes_the_current_metric(self) -> None:
        source = (SRC / "agdesign2" / "portal.py").read_text(encoding="utf-8")
        self.assertIn("BLOSUM62 alignment (gap open 11, extend 1) with free end gaps", source)

    def test_help_page_names_the_withheld_case(self) -> None:
        source = (SRC / "agdesign2" / "portal.py").read_text(encoding="utf-8")
        self.assertIn("the ortholog was withheld because it failed the acceptance gate", source)


if __name__ == "__main__":
    unittest.main()
