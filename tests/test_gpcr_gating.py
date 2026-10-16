"""PR-5: GPCR classification, fusion placement and the GPCRdb boundary ruling.

* **F01** (536 human / 508 mouse pages) — GPCR advice was gated on
  ``tm_count >= 7``, a proxy, not a classification. ABC and SLC transporters,
  CNG and TRP channels, Na+/H+ exchangers and PAQR receptors were handed a
  "Conservative 7TM expression series", a loop screen and — on 329 / 295 pages
  — a generated BRIL/T4L fusion sequence.
* **F02-A** (683 / 261) — the only guard on a fusion window was overlap with
  UniProt TRANSMEM, which GPCRdb's structure-based helices extend well past.
* **F02-B** (84 / 227 reports with a SURFY topology) — SURFY segments arrive
  after the tail and the loops are computed, so the tail stays None and a
  C-terminal segment is offered as an internal loop.
* **F02-C** (707 / 393) — the advisory loop screen ran with no blocker at all.
* **F17-3** (80 / 71) — the design region runs into GPCRdb TM1.
* **OA-009** — a cassette replacing part of an annotated DOMAIN.

The band constants are the owner's decision of 2026-10-15: the published A2A
BRIL window, generic 5.70 to 6.24, with "no cut inside an annotated TM segment"
as the fallback where no generic numbering exists.
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
from agdesign2.models import (  # noqa: E402
    DomainAnnotation,
    Feature,
    GPCRdbAnnotation,
    GPCRdbSegment,
    Region,
    TargetRecord,
)
from agdesign2.pipeline import (  # noqa: E402
    GPCR_ICL3_FUSION_DEPTH,
    AntigenAnalyzer,
    is_gpcr_by_name,
)
from agdesign2.topology import apply_surfy_topology_fallback, derive_ectodomain  # noqa: E402


def _target(protein_name: str, entry_name: str = "TEST_HUMAN") -> TargetRecord:
    return TargetRecord(
        accession="Q00000",
        entry_name=entry_name,
        gene_symbol="TEST",
        protein_name=protein_name,
        organism="Homo sapiens",
        taxon_id="9606",
        sequence="M" * 510,
    )


class GpcrClassificationTests(unittest.TestCase):
    """F01, under the owner's decision to keep known GPCRs (plan §5.1)."""

    GPCRS = (
        "Olfactory receptor 51A7",
        "Vomeronasal type-1 receptor 4",
        "Taste receptor type 2 member 38",
        "Rhodopsin",
        "Short-wave-sensitive opsin 1",
        "Beta-2 adrenergic receptor",
        "Frizzled-1",
        "Adhesion G protein-coupled receptor B1",
        "Trace amine-associated receptor 1",
        "Mas-related G-protein coupled receptor member X1",
        "5-hydroxytryptamine receptor 1A",
    )
    NOT_GPCRS = (
        "ATP-binding cassette sub-family A member 9",
        "Cyclic nucleotide-gated cation channel alpha-1",
        "Adiponectin receptor protein 1",
        "Sodium/hydrogen exchanger 4",
        "Adenylate cyclase type 3",
        "Excitatory amino acid transporter 4",
        "Transmembrane protein 245",
        # The three the plan names: a token match that is not a receptor.
        "G protein-coupled receptor-associated protein LMBRD2",
        "Integral membrane protein GPR137C",
        "Integral membrane protein GPR180",
        # Ligand-gated ion channels whose names look like GPCR families.
        "P2X purinoceptor 6",
        "5-hydroxytryptamine receptor 3A",
        "Modulator of smoothened protein",
    )

    def test_known_gpcrs_are_classified_without_a_gpcrdb_entry(self) -> None:
        for name in self.GPCRS:
            with self.subTest(protein_name=name):
                self.assertTrue(is_gpcr_by_name(_target(name), []))

    def test_non_gpcrs_are_not(self) -> None:
        for name in self.NOT_GPCRS:
            with self.subTest(protein_name=name):
                self.assertFalse(is_gpcr_by_name(_target(name), []))

    def test_a_named_olfactory_receptor_with_no_gpcrdb_entry_keeps_its_advice(self) -> None:
        """The owner's decision, asserted directly (plan §5.1)."""
        analyzer = AntigenAnalyzer(config=AnalysisConfig(generate_assets=False, prefer_precomputed_references=False))
        topology = _seven_tm_topology(analyzer.config)
        suggestions = analyzer._suggest_advanced_membrane_engineering(
            target=_target("Olfactory receptor 51A7", "OR51A7_HUMAN"),
            features=[],
            topology=topology.topology,
            interpro_annotations=[],
            gpcr_annotation=None,
        )
        self.assertTrue(any(s.category.startswith("gpcr_") for s in suggestions))

    def test_a_seven_tm_transporter_loses_it(self) -> None:
        analyzer = AntigenAnalyzer(config=AnalysisConfig(generate_assets=False, prefer_precomputed_references=False))
        topology = _seven_tm_topology(analyzer.config)
        target = _target("ATP-binding cassette sub-family A member 9", "ABCA9_HUMAN")
        self.assertEqual(
            analyzer._suggest_advanced_membrane_engineering(
                target=target, features=[], topology=topology.topology,
                interpro_annotations=[], gpcr_annotation=None,
            ),
            [],
        )
        self.assertEqual(
            analyzer._suggest_gpcr_engineering_variants(
                target=target, topology=topology.topology, features=[],
                ligand_interactions=[], ptms=[], gpcr_annotation=None,
            ),
            [],
        )


    def test_a_target_classified_only_by_interpro_still_gets_fusion_variants(self) -> None:
        """The variant gate once passed an empty InterPro list to the classifier."""
        analyzer = AntigenAnalyzer(config=AnalysisConfig(generate_assets=False, prefer_precomputed_references=False))
        topology = _seven_tm_topology(analyzer.config)
        target = _target("Uncharacterized protein", "UNCH_HUMAN")
        interpro = [DomainAnnotation(
            accession="IPR000276", name="G protein-coupled receptor, rhodopsin-like",
            type="family", source_database="interpro", start=25, end=500,
        )]
        self.assertEqual(
            analyzer._suggest_gpcr_engineering_variants(
                target=target, topology=topology.topology, features=[],
                ligand_interactions=[], ptms=[], gpcr_annotation=None,
                interpro_annotations=[],
            ),
            [],
        )
        self.assertTrue(
            analyzer._suggest_gpcr_engineering_variants(
                target=target, topology=topology.topology, features=[],
                ligand_interactions=[], ptms=[], gpcr_annotation=None,
                interpro_annotations=interpro,
            )
        )

def _seven_tm_topology(config: AnalysisConfig):
    features = [Feature(type="SIGNAL", start=1, end=24, description="Signal peptide")]
    position = 25
    for index in range(7):
        features.append(Feature(type="TOPO_DOM", start=position, end=position + 24,
                                description="Extracellular" if index % 2 == 0 else "Cytoplasmic"))
        position += 25
        features.append(Feature(type="TRANSMEM", start=position, end=position + 20, description="Helical"))
        position += 21
    features.append(Feature(type="TOPO_DOM", start=position, end=510, description="Cytoplasmic"))
    return derive_ectodomain(features, 510, config)


class FusionBandTests(unittest.TestCase):
    """F02-A: the band is in GPCRdb generic numbers, not raw overlap."""

    # A2A's own GPCRdb segments, retrieved 2026-10-15 (plan §5.6).
    A2A = GPCRdbAnnotation(
        entry_name="aa2ar_human",
        accession="P29274",
        name="Adenosine receptor A2a",
        segments=[
            GPCRdbSegment(name="TM5", start=173, end=213, generic_start="5.34x36", generic_end="5.74x74"),
            GPCRdbSegment(name="ICL3", start=214, end=218),
            GPCRdbSegment(name="TM6", start=219, end=259, generic_start="6.21x21", generic_end="6.61x61"),
        ],
    )

    def _blockers(self, cut_start: int, cut_end: int, annotation=None):
        """Blockers for a construct that CUTS at ``cut_start``/``cut_end``.

        The cassette keeps ``gpcr_loop_fusion_retain_residues`` at each end, so
        the loop is wider than the cut. The published A2A numbers below are cut
        points, so the loop is widened here to match.
        """
        config = AnalysisConfig(generate_assets=False, prefer_precomputed_references=False)
        retain = config.gpcr_loop_fusion_retain_residues
        analyzer = AntigenAnalyzer(config=config)
        return analyzer._gpcr_band_blockers(
            loop=Region(start=cut_start - retain, end=cut_end + retain, label="ICL3", source="test"),
            gpcr_annotation=self.A2A if annotation is None else annotation,
        )

    def test_the_published_window_is_allowed(self) -> None:
        # Lys209-Gly218: 209 is generic 5.70 in TM5, 218 is inside ICL3.
        self.assertEqual(self._blockers(209, 218), [])

    def test_both_published_readings_of_the_a2a_construct_are_allowed(self) -> None:
        """The literature reports the canonical construct two ways.

        Lys209-Gly218 (5.70) and L208/R222 (5.69 and 6.24). The owner's
        follow-up decision set the TM5 floor to 5.69 so a rule meant to say
        "do not cut deeper than a construct known to fold and crystallise"
        does not block that construct by one residue.
        """
        self.assertEqual(self._blockers(209, 218), [])
        self.assertEqual(self._blockers(208, 222), [])

    def test_the_portals_own_a2a_variant_is_blocked(self) -> None:
        # It replaces 204-229, ~14 residues of TM5 and ~16 of TM6 — taken from
        # the F02-A fallback path itself.
        blockers = self._blockers(204, 229)
        self.assertTrue(blockers)
        self.assertIn("residues into GPCRdb", blockers[0])

    def test_a_cut_one_residue_too_deep_in_tm5_is_blocked(self) -> None:
        # A2A's TM5 ends at 213, so 208 keeps 5 residues and is allowed;
        # 207 keeps 6 and is not.
        self.assertEqual(self._blockers(208, 218), [])
        self.assertTrue(self._blockers(207, 218))

    def test_a_cut_one_residue_too_deep_in_tm6_is_blocked(self) -> None:
        # A2A's TM6 starts at 219, so 222 keeps 3 residues and is allowed;
        # 223 keeps 4 and is not.
        self.assertEqual(self._blockers(209, 222), [])
        self.assertTrue(self._blockers(209, 223))

    def test_a_loop_clear_of_both_helices_is_allowed(self) -> None:
        self.assertEqual(self._blockers(214, 218), [])

    def test_a_receptor_without_generic_numbering_uses_the_same_depths(self) -> None:
        """Class B, C and F receptors carry no class-A generic numbering.

        The depth is measured from the annotated segment boundary, which every
        class has, so they no longer need the "block any helix cut" fallback
        the absolute window required (owner decision 2026-10-16).
        """
        class_b = GPCRdbAnnotation(
            entry_name="x", accession="Q1", name="A class-B receptor",
            segments=[GPCRdbSegment(name="TM5", start=173, end=213)],
        )
        self.assertEqual(self._blockers(209, 218, annotation=class_b), [])
        self.assertTrue(self._blockers(200, 218, annotation=class_b))

    def test_a_helix_with_no_published_construct_is_blocked(self) -> None:
        other = GPCRdbAnnotation(
            entry_name="x", accession="Q1", name="x",
            segments=[GPCRdbSegment(name="TM3", start=173, end=213)],
        )
        blockers = self._blockers(209, 218, annotation=other)
        self.assertTrue(blockers)
        self.assertIn("only for TM5 and TM6", blockers[0])

    def test_no_annotation_means_no_band_blocker(self) -> None:
        analyzer = AntigenAnalyzer(config=AnalysisConfig(generate_assets=False))
        self.assertEqual(
            analyzer._gpcr_band_blockers(
                loop=Region(start=1, end=10, label="l", source="t"), gpcr_annotation=None
            ),
            [],
        )

    def test_the_depth_constants_are_the_owners_decision(self) -> None:
        # A2A's L208 is 5 residues before the end of its TM5 (173-213) and
        # R222 is 3 after the start of its TM6 (219-259). Owner decision
        # 2026-10-16: the same construct, expressed as a depth so it travels.
        self.assertEqual(GPCR_ICL3_FUSION_DEPTH["5"], 5)
        self.assertEqual(GPCR_ICL3_FUSION_DEPTH["6"], 3)


class DomainBlockerTests(unittest.TestCase):
    """OA-009."""

    def test_a_loop_overlapping_a_domain_is_blocked(self) -> None:
        analyzer = AntigenAnalyzer(config=AnalysisConfig(generate_assets=False, prefer_precomputed_references=False))
        topology = _seven_tm_topology(analyzer.config).topology
        blockers = analyzer._gpcr_loop_engineering_blockers(
            loop=Region(start=300, end=330, label="ICL", source="t"),
            topology=topology, features=[
                Feature(type="DOMAIN", start=310, end=340, description="PDZ"),
            ],
            ligand_interactions=[], ptms=[],
            domains=[Feature(type="DOMAIN", start=310, end=340, description="PDZ")],
        )
        self.assertTrue(any("DOMAIN annotation 310-340" in b for b in blockers))


class SurfyTailTests(unittest.TestCase):
    """F02-B, with the or10ac1_mouse numbers."""

    def test_the_tail_is_recomputed_and_segments_are_clipped(self) -> None:
        config = AnalysisConfig(generate_assets=False, prefer_precomputed_references=False)
        topology = _seven_tm_topology(config)
        annotation = topology.topology
        annotation.cytoplasmic_regions.clear()
        annotation.extracellular_regions.clear()
        annotation.transmembrane_regions.clear()
        annotation.transmembrane_regions.append(
            Region(start=276, end=296, label="TM7", source="UniProt")
        )
        annotation.cytoplasmic_tail_region = None
        # The SURFY topology string is "KIND:start-end;..." (topology.py:447).
        changed = apply_surfy_topology_fallback(
            annotation, "NC:1-264;CY:265-317", sequence_length=317
        )
        self.assertTrue(changed)
        # The SURFY CY segment 265-317 straddles TM7 276-296; the longer half
        # (297-317) survives and becomes the tail.
        self.assertIsNotNone(annotation.cytoplasmic_tail_region)
        self.assertEqual(annotation.cytoplasmic_tail_region.start, 297)
        self.assertTrue(
            all(
                region.start > 296 or region.end < 276
                for region in annotation.cytoplasmic_regions
            ),
            "a cytoplasmic region still reaches inside TM7",
        )


if __name__ == "__main__":
    unittest.main()


class GpcrWiringTests(unittest.TestCase):
    """PR-5b: items a red-team round showed no test exercised.

    Reverting each of these left the whole suite green: the band blocker's
    wiring into the variant path, the F02-C advisory blockers, the F17-3 clip,
    and the GPCRdb non-404 raise. FusionBandTests called the band helper
    directly, so the helper was pinned but its use was not.
    """

    A2A = FusionBandTests.A2A

    def _analyzer(self) -> AntigenAnalyzer:
        return AntigenAnalyzer(config=AnalysisConfig(generate_assets=False, prefer_precomputed_references=False))

    def test_the_band_blocker_is_wired_into_the_loop_blockers(self) -> None:
        topology = _seven_tm_topology(AnalysisConfig()).topology
        blockers = self._analyzer()._gpcr_loop_engineering_blockers(
            loop=Region(start=199, end=234, label="ICL3", source="t"),
            topology=topology, features=[], ligand_interactions=[], ptms=[],
            gpcr_annotation=self.A2A,
        )
        self.assertTrue([b for b in blockers if "ICL3-fusion construct" in b], blockers)

    def test_a_helix_without_a_window_is_not_said_to_lack_numbering(self) -> None:
        # Red team: 106 messages claimed "no generic numbering" for TM3/TM4
        # segments that carry it. Only TM5/TM6 have a published window.
        annotation = GPCRdbAnnotation(
            entry_name="x", accession="Q1", name="x",
            segments=[GPCRdbSegment(name="TM3", start=100, end=130,
                                    generic_start="3.22x22", generic_end="3.52x52")],
        )
        blockers = self._analyzer()._gpcr_band_blockers(
            loop=Region(start=110, end=125, label="ICL", source="t"), gpcr_annotation=annotation,
        )
        self.assertTrue(blockers)
        self.assertNotIn("no generic numbering", blockers[0])
        self.assertIn("only for TM5 and TM6", blockers[0])

    def test_the_advisory_screen_says_when_nothing_is_eligible(self) -> None:
        topology = _seven_tm_topology(AnalysisConfig()).topology
        suggestions = self._analyzer()._suggest_advanced_membrane_engineering(
            target=_target("Beta-2 adrenergic receptor"),
            features=[Feature(type="DOMAIN", start=1, end=510, description="covers everything")],
            topology=topology, interpro_annotations=[], gpcr_annotation=None,
        )
        categories = [s.category for s in suggestions]
        self.assertIn("gpcr_loop_engineering_blocked", categories)
        self.assertNotIn("gpcr_loop_engineering", categories)

    def test_the_tm1_clip_trims_and_notes(self) -> None:
        from agdesign2.topology import TopologyResult

        annotation = GPCRdbAnnotation(
            entry_name="agrg3", accession="Q1", name="x",
            segments=[GPCRdbSegment(name="TM1", start=261, end=297)],
        )
        result = TopologyResult(
            ectodomain=Region(start=21, end=270, label="Extracellular", source="topology"),
            notes=[], topology=None,
        )
        notes: list = []
        self._analyzer()._clip_ectodomain_to_gpcrdb_tm1(result, annotation, notes)
        self.assertEqual(result.ectodomain.end, 260)
        self.assertTrue(notes and "GPCRdb places TM1 at 261-297" in notes[0].message)

    def test_the_tm1_clip_runs_before_the_interpro_filter(self) -> None:
        """Removing the call site left the suite green, and the ORDER matters.

        InterPro annotations surviving the design-region filter become domain_*
        and repeat_* constructs, so filtering against the un-clipped region put
        a domain reaching into GPCRdb TM1 back on offer. Deliberately a
        source-order assertion: the invariant is the order of two statements
        inside ``analyze_target``, and both sit past target resolution, so
        there is no cheaper behavioural probe.
        """
        source = (SRC / "agdesign2" / "pipeline.py").read_text(encoding="utf-8")
        clip = source.index("self._clip_ectodomain_to_gpcrdb_tm1(topology, gpcr_annotation, notes)")
        design_filter = source.index(
            "interpro_annotations = self._filter_design_annotations(interpro_annotations, topology.ectodomain)"
        )
        self.assertLess(clip, design_filter)

    def test_a_non_404_gpcrdb_error_raises(self) -> None:
        from agdesign2.exceptions import ExternalServiceError
        from agdesign2.gpcrdb import _is_not_found

        self.assertTrue(_is_not_found(ExternalServiceError("HTTP error for https://x: 404")))
        self.assertFalse(_is_not_found(ExternalServiceError("HTTP error for https://x: 503")))
        self.assertFalse(_is_not_found(ExternalServiceError("Connection error for https://x: offline")))

    def test_smoothened_is_a_gpcr_and_its_modulator_is_not(self) -> None:
        self.assertTrue(is_gpcr_by_name(_target("Protein smoothened"), []))
        self.assertFalse(is_gpcr_by_name(_target("Modulator of smoothened protein"), []))
        self.assertFalse(is_gpcr_by_name(_target("Small cell adhesion glycoprotein"), []))


class SurfyClipKeepsBothPiecesTests(unittest.TestCase):
    """PR-5b: a region split by a helix keeps every piece outside it."""

    def test_both_sides_of_a_helix_survive(self) -> None:
        from agdesign2.topology import _clip_regions_to_transmembrane

        regions = [Region(start=265, end=317, label="CY", source="SURFY")]
        _clip_regions_to_transmembrane(regions, [Region(start=276, end=296, label="TM7", source="UniProt")])
        self.assertEqual([(r.start, r.end) for r in regions], [(265, 275), (297, 317)])
