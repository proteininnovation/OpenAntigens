"""PR-2: three mouse-render defects, each with a human counter-case.

* **F06a** — a mouse page listed its own species as "Mapping unavailable" on
  4,359 of 4,702 pages (human: 0). ``_collect_species_matches`` already skips
  the target's own species (``pipeline.py:1485``), but
  ``_map_construct_to_homologs`` looped a hard-coded ``("mouse",
  "macaca_fascicularis")`` pair regardless. Fixed at the cause, and filtered in
  the renderer so the pages already built that way lose the row on re-render.
* **F19-4 / OA-008** — mouse Methods claimed every target is analyzed against
  reviewed mouse UniProt data; 642 of 4,702 reports (13.7%) are RefSeq-sourced.
* **F21** — the mouse agent guide's data-access links resolved to the *human*
  index, manifest and metadata, on the file an agent fetches first.
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from agdesign2.config import AnalysisConfig  # noqa: E402
from agdesign2.models import ConstructSuggestion, Region  # noqa: E402
from agdesign2.pipeline import AntigenAnalyzer  # noqa: E402
from agdesign2.portal import (  # noqa: E402
    _PER_PORTAL_GUIDE_ARTIFACTS,
    render_agent_guide_page,
    render_methods_page,
)


def _map(target_species: str):
    analyzer = AntigenAnalyzer(
        config=AnalysisConfig(
            generate_assets=False,
            prefer_precomputed_references=False,
            target_species=target_species,
        )
    )
    return analyzer._map_construct_to_homologs(
        construct=ConstructSuggestion(
            name="full_ectodomain", start=1, end=281, score=88.0,
            rationale="test", classification="full_ectodomain",
        ),
        target_sequence="A" * 330,
        homolog_context={},
        ectodomain=Region(start=25, end=180, label="Major extracellular region", source="topology"),
    )


class SelfSpeciesHomologTests(unittest.TestCase):
    """F06a, at the cause."""

    def test_mouse_target_gets_no_mouse_homolog_row(self) -> None:
        self.assertEqual([h.species for h in _map("mouse")], ["macaca_fascicularis"])

    def test_human_target_still_gets_both_rows(self) -> None:
        self.assertEqual(
            [h.species for h in _map("human")], ["mouse", "macaca_fascicularis"]
        )


class SelfSpeciesRenderTests(unittest.TestCase):
    """F06a, on the 4,359 pages already built with the row baked in."""

    def _card(self, entry_name: str, species: str) -> str:
        from tempfile import TemporaryDirectory

        from agdesign2.portal import _render_construct_card

        construct = {
            "name": "full_ectodomain",
            "start": 25,
            "end": 180,
            "sequence": "A" * 156,
            "homologs": [
                {"species": species, "available": False,
                 "notes": ["No same-name UniProt entry found for species."]},
                {"species": "macaca_fascicularis", "available": False,
                 "notes": ["No same-name UniProt entry found for species."]},
            ],
        }
        target = {"entry_name": entry_name, "gene_symbol": "Il7", "accession": "P10168",
                  "taxon_id": "10090" if entry_name.endswith("_MOUSE") else "9606"}
        with TemporaryDirectory() as tmp:
            return _render_construct_card(
                construct, target=target, batch_dir=Path(tmp), page_dir=Path(tmp),
                furin_sites=[], ptms=[],
            )

    def test_mouse_card_drops_the_self_species_row(self) -> None:
        html = self._card("IL7_MOUSE", "mouse")
        self.assertNotIn("_MOUSE_mapping", html)
        self.assertIn("_MACFA_mapping", html)

    def test_human_card_keeps_the_mouse_row(self) -> None:
        html = self._card("IL7_HUMAN", "mouse")
        self.assertIn("_MOUSE_mapping", html)
        self.assertIn("_MACFA_mapping", html)


class MouseMethodsProvenanceTests(unittest.TestCase):
    """F19-4 / Codex OA-008."""

    def test_mouse_methods_does_not_claim_every_target_is_reviewed_uniprot(self) -> None:
        mouse = render_methods_page(portal_title="OpenAntigens Mouse")
        self.assertNotIn("Each included target is analyzed against reviewed mouse UniProt data", mouse)
        self.assertIn("Most included targets are analyzed against reviewed mouse UniProt data", mouse)
        self.assertIn("RefSeq record named in the report", mouse)
        # The provenance table two paragraphs below said the opposite.
        self.assertNotIn("<td>Reviewed UniProt mouse records</td>", mouse)
        self.assertNotIn("<td>Reviewed UniProt homolog lookup</td>", mouse)
        self.assertIn("where no reviewed mouse entry exists", mouse)

    def test_human_methods_is_untouched(self) -> None:
        human = render_methods_page()
        self.assertIn("The human target universe is sourced from the refined cell-surface", human)
        self.assertNotIn("Most included targets are analyzed against reviewed mouse UniProt", human)


class AgentGuideLinkTests(unittest.TestCase):
    """F21."""

    # Listed explicitly, not read from _PER_PORTAL_GUIDE_ARTIFACTS: a test that
    # iterates the constant it is testing passes when the constant shrinks.
    PER_PORTAL = (
        "downloads.html",
        "help.html",
        "builder.html",
        "constructs.html",
        "methods.html",
        "downloads/agdesign2_portal_index.json",
        "downloads/agdesign2_portal_index.tsv",
        "downloads/download_manifest.json",
        "portal_metadata.json",
    )

    def test_the_constant_still_covers_every_per_portal_artifact(self) -> None:
        self.assertEqual(set(_PER_PORTAL_GUIDE_ARTIFACTS), set(self.PER_PORTAL))

    def test_mouse_guide_points_at_the_mouse_artifacts(self) -> None:
        mouse = render_agent_guide_page(portal_title="OpenAntigens Mouse")
        for artifact in self.PER_PORTAL:
            with self.subTest(artifact=artifact):
                self.assertIn(f'href="{artifact}"', mouse)
                self.assertNotIn(f'href="../{artifact}"', mouse)

    def test_mouse_guide_names_the_mouse_index_root(self) -> None:
        mouse = render_agent_guide_page(portal_title="OpenAntigens Mouse")
        self.assertIn("Mouse index paths are relative to https://openantigens.org/mouse/", mouse)
        self.assertNotIn("Human index paths are relative to", mouse)

    def test_mouse_guide_keeps_the_shared_site_links_relative_to_the_human_root(self) -> None:
        # Deliberate: the citation and the human portal itself are shared, and
        # tests/test_site_docs.py pins them.
        mouse = render_agent_guide_page(portal_title="OpenAntigens Mouse")
        self.assertIn('href="../index.html"', mouse)
        self.assertIn('href="../help.html#cite-openantigens"', mouse)

    def test_human_guide_is_unchanged(self) -> None:
        human = render_agent_guide_page()
        for artifact in self.PER_PORTAL:
            self.assertIn(f'href="{artifact}"', human)
        self.assertIn("Human index paths are relative to https://openantigens.org/", human)


if __name__ == "__main__":
    unittest.main()
