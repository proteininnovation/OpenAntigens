"""PR-6: cysteine pairing and exposure.

* **F12-A** (428 human / 324 mouse pages, 545 / 218 where the model picked a
  *different* partner) — pairing was geometry-only, so a curated UniProt
  DISULFID whose sulfurs sit further than 2.4 A apart in the AlphaFold model
  was reported "unpaired" and both ends were offered for Cys→Ser.
* **F12-B** (74 / 53) — "unpaired" where the curated partner is simply outside
  the construct. The offer stays (owner decision, plan §5.4); only the wording
  changes.
* **F12-C** (207 / 159) — an interchain bond (UniProt writes it ``n-n``) is
  unpaired in a monomer by construction, and Cys→Ser was offered on it.
* **F12-D** (124 / 118) — cysteines inside a SIGNAL span listed as liabilities.
* **F11** (1,473 / 1,263) — residues inside a membrane span called surface
  exposed. 71% of those calls come from computing SASA on an isolated helix.
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
from agdesign2.models import Feature, Region, TopologyAnnotation  # noqa: E402
from agdesign2.pipeline import AntigenAnalyzer  # noqa: E402
from agdesign2.portal import _NO_CYS_TO_SER_WARNING_PREFIXES  # noqa: E402
from agdesign2.structure_utils import INTERCHAIN_DISULFIDE_WARNING, analyze_cysteines  # noqa: E402


def _atoms(positions: dict[int, tuple[float, float, float]]):
    return {residue: {"SG": coords, "CA": coords} for residue, coords in positions.items()}


# Two cysteines 3.36 A apart — the ACVRL1 C46/C69 distance, outside the 2.4 A
# geometric cutoff but a curated bond.
FAR_PAIR = _atoms({46: (0.0, 0.0, 0.0), 69: (3.36, 0.0, 0.0)})
NAMES = {46: "CYS", 69: "CYS"}


class CuratedDisulfideTests(unittest.TestCase):
    """F12-A."""

    def test_a_curated_pair_beats_model_geometry(self) -> None:
        findings = analyze_cysteines(
            FAR_PAIR, NAMES, 1, 100,
            curated_disulfides=[Feature(type="DISULFID", start=46, end=69, description="")],
        )
        by_position = {f.position: f for f in findings}
        self.assertEqual(by_position[46].paired_with, 69)
        self.assertEqual(by_position[69].paired_with, 46)
        self.assertIsNone(by_position[46].warning)

    def test_without_the_annotation_geometry_is_not_relaxed(self) -> None:
        findings = analyze_cysteines(FAR_PAIR, NAMES, 1, 100)
        self.assertTrue(all(f.paired_with is None for f in findings))

    def test_a_curated_pair_the_model_pairs_differently_is_overridden(self) -> None:
        # 545 human / 218 mouse pages: the model pairs a different partner.
        atoms = _atoms({46: (0.0, 0.0, 0.0), 69: (3.36, 0.0, 0.0), 80: (2.0, 0.0, 0.0)})
        names = {46: "CYS", 69: "CYS", 80: "CYS"}
        findings = analyze_cysteines(
            atoms, names, 1, 100,
            curated_disulfides=[Feature(type="DISULFID", start=46, end=69, description="")],
        )
        by_position = {f.position: f for f in findings}
        self.assertEqual(by_position[46].paired_with, 69)
        self.assertIsNone(by_position[80].paired_with)

    def test_a_partner_outside_the_range_is_named_not_called_absent(self) -> None:
        """F12-B. The residue stays unpaired — it is, in this construct."""
        findings = analyze_cysteines(
            FAR_PAIR, NAMES, 1, 50,
            curated_disulfides=[Feature(type="DISULFID", start=46, end=69, description="")],
        )
        finding = next(f for f in findings if f.position == 46)
        self.assertIsNone(finding.paired_with)
        self.assertEqual(finding.partner_outside_range, 69)
        self.assertIn("outside this range", finding.warning)

    def test_an_interchain_bond_is_labelled(self) -> None:
        """F12-C. UniProt writes an interchain bond as n-n."""
        findings = analyze_cysteines(
            _atoms({46: (0.0, 0.0, 0.0)}), {46: "CYS"}, 1, 100,
            curated_disulfides=[Feature(type="DISULFID", start=46, end=46, description="Interchain")],
        )
        finding = findings[0]
        self.assertIsNone(finding.paired_with)
        self.assertEqual(finding.warning, INTERCHAIN_DISULFIDE_WARNING)

    def test_the_interchain_warning_blocks_cys_to_ser(self) -> None:
        self.assertTrue(
            INTERCHAIN_DISULFIDE_WARNING.startswith(_NO_CYS_TO_SER_WARNING_PREFIXES)
        )

    def test_an_orphan_with_no_annotation_still_offers_cys_to_ser(self) -> None:
        # Owner decision, plan §5.4: the offer is only withdrawn when BOTH ends
        # are in scope.
        findings = analyze_cysteines(_atoms({46: (0.0, 0.0, 0.0)}), {46: "CYS"}, 1, 100)
        finding = findings[0]
        self.assertIsNone(finding.paired_with)
        self.assertFalse(
            str(finding.warning or "").startswith(_NO_CYS_TO_SER_WARNING_PREFIXES)
        )


class MembraneExposureTests(unittest.TestCase):
    """F11, under the owner's decision to scope it to membrane spans (§5.2)."""

    def _withhold(self, exposure: dict[int, bool | None]) -> dict[int, bool | None]:
        analyzer = AntigenAnalyzer(config=AnalysisConfig(generate_assets=False))
        topology = TopologyAnnotation(
            topology_class="multipass_compact",
            extracellular_regions=[Region(start=1, end=38, label="e", source="t")],
            transmembrane_regions=[Region(start=39, end=59, label="TM1", source="UniProt")],
            cytoplasmic_regions=[Region(start=60, end=100, label="c", source="t")],
        )
        return analyzer._withhold_membrane_exposure(
            exposure, sequence_length=100, topology=topology
        )

    def test_a_membrane_residue_gets_no_call(self) -> None:
        # ACVRL1-style: 5HT1A C49 and C56 sit inside TM1 39-59.
        result = self._withhold({49: True, 56: True})
        self.assertIsNone(result[49])
        self.assertIsNone(result[56])

    def test_the_key_is_kept_so_the_second_heuristic_does_not_fire(self) -> None:
        # analyze_cysteines falls back to _is_surface_exposed for an ABSENT
        # residue, so the value has to be an explicit None.
        self.assertIn(49, self._withhold({49: True}))

    def test_extracellular_and_intracellular_residues_keep_their_call(self) -> None:
        result = self._withhold({10: True, 70: True, 80: False})
        self.assertIs(result[10], True)
        self.assertIs(result[70], True)
        self.assertIs(result[80], False)

    def test_no_exposure_means_no_exposure_warning(self) -> None:
        findings = analyze_cysteines(
            _atoms({49: (0.0, 0.0, 0.0)}), {49: "CYS"}, 1, 100,
            surface_exposure={49: None},
        )
        self.assertIsNone(findings[0].surface_exposed)
        self.assertIsNone(findings[0].warning)


if __name__ == "__main__":
    unittest.main()
