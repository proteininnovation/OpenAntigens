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


class ConstructCardTests(unittest.TestCase):
    """A cysteine that is not a free thiol must not be called unpaired.

    The Cys->Ser offer has filtered on the warning prefix since F12-C; the
    construct card did not, so interchain and ambiguously annotated cysteines
    were listed in "Contains unpaired cysteine(s)".
    """

    def _warnings(self, disulfides: list[Feature]) -> list[str]:
        from agdesign2.models import ConstructSuggestion, Region

        analyzer = AntigenAnalyzer(config=AnalysisConfig(generate_assets=False, prefer_precomputed_references=False))
        # C30-C40 pair by geometry, so C10 is the only cysteine in question.
        names = {10: "CYS", 30: "CYS", 40: "CYS"}
        atoms = _atoms({10: (0.0, 0.0, 0.0), 30: (50.0, 0.0, 0.0), 40: (52.0, 0.0, 0.0)})
        details = analyzer._build_construct_details(
            target_entry_name="TEST_HUMAN",
            target_sequence="A" * 100,
            constructs=[ConstructSuggestion(name="c", start=1, end=100, score=1.0, rationale="t")],
            homolog_context={},
            ectodomain=Region(start=1, end=100, label="t", source="t"),
            ligand_interactions=[],
            interpro_annotations=[],
            ptms=disulfides,
            structure_atoms=atoms,
            residue_names=names,
        )
        return [w for w in details[0].warnings if "unpaired cysteine" in w]

    def test_an_interchain_cysteine_is_not_listed_as_unpaired(self) -> None:
        self.assertEqual(
            self._warnings([Feature(type="DISULFID", start=10, end=10, description="Interchain")]),
            [],
        )

    def test_an_ambiguous_cysteine_is_not_listed_as_unpaired(self) -> None:
        self.assertEqual(
            self._warnings([Feature(type="DISULFID", start=10, end=30, description=""),
                            Feature(type="DISULFID", start=10, end=40, description="")]),
            [],
        )

    def test_a_genuinely_free_cysteine_still_is(self) -> None:
        self.assertTrue(self._warnings([]))

class PairingConsistencyTests(unittest.TestCase):
    """PR-6b: the curated override left the table contradicting itself.

    Writing paired[a] = b without revisiting whoever the model had paired to
    either end gave C10->C30, C20->C10, C30->C10. C20 claimed a partner that
    claimed someone else, so a genuinely free cysteine was never listed as
    unpaired and never offered for Cys->Ser. The configuration is the commit's
    own "model pairs a different partner" case: 509 human / 198 mouse pages.
    """

    def test_a_displaced_cysteine_becomes_unpaired(self) -> None:
        atoms = _atoms({10: (0.0, 0.0, 0.0), 20: (2.0, 0.0, 0.0), 30: (9.0, 0.0, 0.0)})
        names = {10: "CYS", 20: "CYS", 30: "CYS"}
        findings = {
            f.position: f
            for f in analyze_cysteines(
                atoms, names, 1, 100,
                curated_disulfides=[Feature(type="DISULFID", start=10, end=30, description="")],
            )
        }
        self.assertEqual(findings[10].paired_with, 30)
        self.assertEqual(findings[30].paired_with, 10)
        self.assertIsNone(findings[20].paired_with)

    def test_the_pairing_is_symmetric(self) -> None:
        atoms = _atoms({10: (0.0, 0.0, 0.0), 20: (2.0, 0.0, 0.0), 30: (9.0, 0.0, 0.0)})
        names = {10: "CYS", 20: "CYS", 30: "CYS"}
        findings = {
            f.position: f
            for f in analyze_cysteines(
                atoms, names, 1, 100,
                curated_disulfides=[Feature(type="DISULFID", start=10, end=30, description="")],
            )
        }
        for finding in findings.values():
            if finding.paired_with is not None:
                self.assertEqual(findings[finding.paired_with].paired_with, finding.position)

    def test_alternative_curated_connectivity_is_not_offered_for_cys_to_ser(self) -> None:
        # I18RA_MOUSE: 154-179, 174-220, 179-220. Last-in-file used to win, so
        # 179-220 stole C179 and C154 came out unpaired with no warning.
        names = {154: "CYS", 174: "CYS", 179: "CYS", 220: "CYS"}
        atoms = _atoms({r: (float(r) * 10, 0.0, 0.0) for r in names})  # no geometric pairs
        bonds = [Feature(type="DISULFID", start=a, end=b, description="")
                 for a, b in ((154, 179), (174, 220), (179, 220))]
        results = []
        for order in (bonds, list(reversed(bonds))):
            findings = {f.position: f for f in analyze_cysteines(
                atoms, names, 1, 300, curated_disulfides=order)}
            for position in names:
                self.assertTrue(
                    str(findings[position].warning or "").startswith(_NO_CYS_TO_SER_WARNING_PREFIXES),
                    f"C{position} would be offered for Cys->Ser",
                )
            results.append({p: (f.paired_with, f.warning) for p, f in findings.items()})
        self.assertEqual(results[0], results[1], "file order decided the pairing")

    def test_geometry_may_not_pick_a_partner_no_annotation_names(self) -> None:
        """18 cysteines asserted a pairing no DISULFID names, with no warning."""
        names = {10: "CYS", 20: "CYS", 30: "CYS", 40: "CYS"}
        # C10 is in two curated bonds (30 and 40); geometry puts it next to C20.
        atoms = _atoms({10: (0.0, 0.0, 0.0), 20: (2.0, 0.0, 0.0),
                        30: (50.0, 0.0, 0.0), 40: (90.0, 0.0, 0.0)})
        findings = {f.position: f for f in analyze_cysteines(
            atoms, names, 1, 100,
            curated_disulfides=[Feature(type="DISULFID", start=10, end=30, description=""),
                                Feature(type="DISULFID", start=10, end=40, description="")],
        )}
        self.assertIsNone(findings[10].paired_with)
        self.assertIsNone(findings[20].paired_with)
        self.assertTrue(
            str(findings[10].warning or "").startswith(_NO_CYS_TO_SER_WARNING_PREFIXES),
            findings[10].warning,
        )

    def test_a_paired_ambiguous_cysteine_still_says_the_bond_is_ambiguous(self) -> None:
        # Geometry picks C30, which IS one of the curated alternatives: the
        # pair is kept, but the page must not present it as settled.
        names = {10: "CYS", 30: "CYS", 40: "CYS"}
        atoms = _atoms({10: (0.0, 0.0, 0.0), 30: (2.0, 0.0, 0.0), 40: (90.0, 0.0, 0.0)})
        findings = {f.position: f for f in analyze_cysteines(
            atoms, names, 1, 100,
            curated_disulfides=[Feature(type="DISULFID", start=10, end=30, description=""),
                                Feature(type="DISULFID", start=10, end=40, description="")],
        )}
        self.assertEqual(findings[10].paired_with, 30)
        self.assertTrue(
            str(findings[10].warning or "").startswith(_NO_CYS_TO_SER_WARNING_PREFIXES),
            findings[10].warning,
        )

    def test_the_same_bond_listed_twice_is_one_bond(self) -> None:
        # IGHG1_HUMAN lists 144-204 twice (blank and "Intrachain").
        names = {10: "CYS", 30: "CYS"}
        atoms = _atoms({10: (0.0, 0.0, 0.0), 30: (90.0, 0.0, 0.0)})
        findings = {f.position: f for f in analyze_cysteines(
            atoms, names, 1, 100,
            curated_disulfides=[Feature(type="DISULFID", start=10, end=30, description=""),
                                Feature(type="DISULFID", start=10, end=30, description="Intrachain")],
        )}
        self.assertEqual(findings[10].paired_with, 30)
        self.assertIsNone(findings[10].warning)

    def test_a_curated_partner_that_is_not_a_cysteine_does_not_destroy_a_pairing(self) -> None:
        # Position 30 is in range but is not CYS in the model. The old code
        # called it "outside this range" (false) and threw away C10's correct
        # geometric pairing to C20.
        atoms = _atoms({10: (0.0, 0.0, 0.0), 20: (2.0, 0.0, 0.0)})
        names = {10: "CYS", 20: "CYS", 30: "SER"}
        findings = {
            f.position: f
            for f in analyze_cysteines(
                atoms, names, 1, 100,
                curated_disulfides=[Feature(type="DISULFID", start=10, end=30, description="")],
            )
        }
        self.assertEqual(findings[10].paired_with, 20)
        self.assertIn("not a cysteine in this model", findings[10].warning)

    def test_a_model_paired_interchain_residue_is_not_told_otherwise(self) -> None:
        atoms = _atoms({10: (0.0, 0.0, 0.0), 20: (2.0, 0.0, 0.0)})
        findings = {
            f.position: f
            for f in analyze_cysteines(
                atoms, {10: "CYS", 20: "CYS"}, 1, 100,
                curated_disulfides=[Feature(type="DISULFID", start=10, end=10, description="Interchain")],
            )
        }
        self.assertEqual(findings[10].paired_with, 20)
        self.assertNotIn("partner is on another chain", findings[10].warning)
        self.assertIn("the model pairs it with C20", findings[10].warning)


class MembraneKeyCoverageTests(unittest.TestCase):
    """PR-6b: every membrane residue needs an explicit key, not just the ones
    the SASA pass covered. The target set is bounded by the design region, so a
    membrane cysteine outside it had no key and fell through to the second
    heuristic — 1,476 construct-level findings on 165 human / 149 mouse pages.
    """

    def test_a_membrane_residue_outside_the_sasa_scope_still_gets_a_key(self) -> None:
        analyzer = AntigenAnalyzer(config=AnalysisConfig(generate_assets=False))
        topology = TopologyAnnotation(
            topology_class="multipass_compact",
            extracellular_regions=[Region(start=1, end=38, label="e", source="t")],
            transmembrane_regions=[Region(start=39, end=59, label="TM1", source="UniProt")],
            cytoplasmic_regions=[Region(start=60, end=100, label="c", source="t")],
        )
        # The SASA pass only covered residue 10; 49 is in TM1 and absent.
        result = analyzer._withhold_membrane_exposure(
            {10: True}, sequence_length=100, topology=topology
        )
        self.assertIn(49, result)
        self.assertIsNone(result[49])

    def test_the_stored_classification_is_withheld_too(self) -> None:
        source = (SRC / "agdesign2" / "pipeline.py").read_text(encoding="utf-8")
        self.assertIn("No solvent-exposure call is made inside a membrane span.", source)


class WarningLabelTests(unittest.TestCase):
    """PR-6b: F12-B's user-visible half and F11's second producer had no test."""

    def test_the_label_names_a_partner_outside_the_construct(self) -> None:
        from agdesign2.models import CysteineFinding

        analyzer = AntigenAnalyzer(config=AnalysisConfig(generate_assets=False))
        label = analyzer._unpaired_cysteine_warning_label(
            CysteineFinding(position=139, paired_with=None, surface_exposed=None,
                            partner_outside_range=162)
        )
        self.assertIn("curated partner C162 is outside this construct", label)

    def test_the_viewer_payload_withholds_membrane_exposure(self) -> None:
        """F11's second producer: _build_residue_annotations feeds the 3D
        viewer and has its own membrane check, separate from the dict."""
        from types import SimpleNamespace

        analyzer = AntigenAnalyzer(config=AnalysisConfig(generate_assets=False))
        topology = TopologyAnnotation(
            topology_class="multipass_compact",
            extracellular_regions=[Region(start=1, end=38, label="e", source="t")],
            transmembrane_regions=[Region(start=39, end=59, label="TM1", source="UniProt")],
        )
        accessibility = SimpleNamespace(
            surface_exposed=True, classification="exposed", reason="exposed in the model",
            full_model_sasa=1.0, full_model_relative_sasa=0.5,
            domain_context_sasa=1.0, domain_context_relative_sasa=0.5, pae_block_id=None,
        )
        annotations = analyzer._build_residue_annotations(
            sequence="M" * 60,
            topology=topology,
            residue_plddt={},
            surface_accessibility={10: accessibility, 49: accessibility},
        )
        by_position = {a.position: a for a in annotations}
        self.assertIs(by_position[10].surface_exposed, True)
        self.assertIsNone(by_position[49].surface_exposed, "residue 49 is inside TM1")
        self.assertIsNone(by_position[49].surface_accessibility)
