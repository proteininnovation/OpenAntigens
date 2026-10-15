"""PR-7: design regions, construct names and the annotations behind them.

* **F10a** (10 human / 9 mouse) — every UniProt TRANSMEM was counted as a
  helix and the description was never read, so MACPF and gasdermin
  beta-hairpins and the VDAC1 beta-barrel made C6/C8A/C8B/C9/PRF1/MPEG1,
  GSDMA/B/D and VDAC1 Track Multipass with a blank ectodomain — while C7, which
  has no TRANSMEM at all, is correctly Secreted 23-843.
* **F16b** — the GPI region ran past the omega site because UniProt gives some
  entries a CHAIN covering the whole sequence (CD14 20-367, site at 345).
* **F16c** (306 / 263) — "Secreted mature chain" containing the propeptide.
* **OA-018** — ``Propeptide``/``Lipidation`` normalised to PROPEPTIDE/LIPIDATION
  while the pipeline tests PROPEP/LIPID, so the propeptide fallback never fired.
* **F14a** (408 / 365) — UniProt REGION "Disordered" offered as a construct.
* **Row 29** (1,068 / 779) — two spans sharing a description got one name.
* **F32** (45 / 27 on a cleavage SITE, 189 / 158 on any processing boundary) —
  a furin motif on an annotated cleavage site offered for mutation with no sign
  that the site is curated and functional.
* **OA-012** (422 / 365 reports) — a full-length construct projected through an
  ectodomain-only alignment, with the coverage denominator hiding the gap.
* **RSM-002** — alternate-isoform feature coordinates applied to the canonical
  record (PPAP, a secreted enzyme, classified single-pass in both species).
* **F10b** (1,531 / 1,347) — a WARNING saying the protein is not what the same
  page has just classified it as.
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
from agdesign2.models import Feature  # noqa: E402
from agdesign2.sequence_utils import global_align, map_query_region_to_subject  # noqa: E402
from agdesign2.topology import (  # noqa: E402
    derive_ectodomain,
    feature_type,
    is_beta_stranded_transmem,
)


class BetaStrandedTransmemTests(unittest.TestCase):
    """F10a."""

    def test_a_beta_strand_is_not_a_helix(self) -> None:
        self.assertTrue(
            is_beta_stranded_transmem(Feature(type="TRANSMEM", start=1, end=10, description="Beta stranded"))
        )
        self.assertFalse(
            is_beta_stranded_transmem(Feature(type="TRANSMEM", start=1, end=10, description="Helical"))
        )
        self.assertFalse(
            is_beta_stranded_transmem(Feature(type="TRANSMEM", start=1, end=10, description="Helical; Name=3"))
        )

    def test_a_macpf_hairpin_protein_is_secreted_not_multipass(self) -> None:
        # C9-shaped: a signal peptide and two beta-hairpins, no helix.
        features = [
            Feature(type="SIGNAL", start=1, end=21, description="Signal peptide"),
            Feature(type="TRANSMEM", start=280, end=300, description="Beta stranded"),
            Feature(type="TRANSMEM", start=320, end=340, description="Beta stranded"),
        ]
        result = derive_ectodomain(features, 559, AnalysisConfig())
        self.assertIsNotNone(result.ectodomain)
        self.assertEqual(result.ectodomain.start, 22)
        self.assertNotIn("multipass", result.topology.topology_class)

    def test_a_real_helix_still_counts(self) -> None:
        features = [
            Feature(type="SIGNAL", start=1, end=21, description="Signal peptide"),
            Feature(type="TRANSMEM", start=280, end=300, description="Helical"),
            Feature(type="TRANSMEM", start=320, end=340, description="Helical"),
        ]
        result = derive_ectodomain(features, 559, AnalysisConfig())
        self.assertIn("multipass", result.topology.topology_class)


class FeatureTypeAliasTests(unittest.TestCase):
    """OA-018, read side."""

    def test_legacy_spellings_are_accepted(self) -> None:
        self.assertEqual(feature_type(Feature(type="PROPEPTIDE", start=1, end=2, description="")), "PROPEP")
        self.assertEqual(feature_type(Feature(type="LIPIDATION", start=1, end=2, description="")), "LIPID")
        self.assertEqual(feature_type(Feature(type="PROPEP", start=1, end=2, description="")), "PROPEP")

    def test_the_ingestion_alias_is_present(self) -> None:
        from agdesign2.uniprot import UniProtClient

        normalize = UniProtClient.__dict__["_normalize_feature_type"]
        self.assertEqual(normalize(None, "Propeptide"), "PROPEP")
        self.assertEqual(normalize(None, "Lipidation"), "LIPID")
        self.assertEqual(normalize(None, "Chain"), "CHAIN")


class GpiRegionCapTests(unittest.TestCase):
    """F16b, with the CD14 numbers."""

    def test_the_region_stops_at_the_omega_site(self) -> None:
        features = [
            Feature(type="SIGNAL", start=1, end=19, description="Signal peptide"),
            Feature(type="CHAIN", start=20, end=375, description="Monocyte differentiation antigen CD14"),
            Feature(type="LIPID", start=345, end=345, description="GPI-anchor amidated asparagine"),
            Feature(type="PROPEP", start=346, end=375, description="Removed in mature form"),
        ]
        result = derive_ectodomain(features, 375, AnalysisConfig())
        self.assertIsNotNone(result.ectodomain)
        self.assertEqual(result.ectodomain.end, 345)

    def test_a_region_already_inside_the_site_is_untouched(self) -> None:
        features = [
            Feature(type="SIGNAL", start=1, end=19, description="Signal peptide"),
            Feature(type="CHAIN", start=20, end=100, description="Mature chain"),
            Feature(type="LIPID", start=120, end=120, description="GPI-anchor amidated serine"),
        ]
        result = derive_ectodomain(features, 150, AnalysisConfig())
        # The anchor is outside every CHAIN, so the chain branch does not apply.
        self.assertIsNotNone(result.ectodomain)


class SecretedPrecursorLabelTests(unittest.TestCase):
    """F16c, with the BMP8B numbers."""

    def test_a_region_holding_a_propeptide_is_not_called_the_mature_chain(self) -> None:
        features = [
            Feature(type="SIGNAL", start=1, end=19, description="Signal peptide"),
            Feature(type="PROPEP", start=20, end=263, description="Removed in mature form"),
            Feature(type="CHAIN", start=264, end=402, description="Bone morphogenetic protein 8B"),
        ]
        result = derive_ectodomain(features, 402, AnalysisConfig())
        self.assertIsNotNone(result.ectodomain)
        self.assertEqual(result.ectodomain.label, "Post-signal precursor region")
        # The interval is deliberately kept: a pro-domain is often wanted.
        self.assertEqual((result.ectodomain.start, result.ectodomain.end), (20, 402))
        self.assertEqual(result.ectodomain.metadata["propeptide_spans"], [[20, 263]])
        self.assertEqual(result.ectodomain.metadata["processed_chain_spans"], [[264, 402]])

    def test_a_region_with_no_propeptide_keeps_the_old_label(self) -> None:
        features = [
            Feature(type="SIGNAL", start=1, end=19, description="Signal peptide"),
            Feature(type="CHAIN", start=20, end=402, description="A mature chain"),
        ]
        result = derive_ectodomain(features, 402, AnalysisConfig())
        self.assertEqual(result.ectodomain.label, "Secreted mature chain")


class DisorderRegionTests(unittest.TestCase):
    """F14a and row 29."""

    def _suggestions(self, features):
        from agdesign2.models import Region
        from agdesign2.pipeline import AntigenAnalyzer

        analyzer = AntigenAnalyzer(config=AnalysisConfig(generate_assets=False, min_construct_length=5))
        return analyzer._suggest_constructs(
            features=features,
            topology=None,
            ectodomain=Region(start=20, end=300, label="Extracellular", source="topology"),
            sequence_length=400,
            structural_regions=[],
            experimental_constructs=[],
            cysteine_analysis=[],
            interpro_annotations=[],
        )

    def test_a_disordered_region_is_not_offered_as_a_construct(self) -> None:
        names = [c.name for c in self._suggestions([
            Feature(type="REGION", start=30, end=120, description="Disordered"),
            Feature(type="DOMAIN", start=150, end=250, description="Ig-like C2-type"),
        ])]
        self.assertFalse([n for n in names if n.startswith("disordered")], names)
        self.assertTrue([n for n in names if n.startswith("ig_like_c2_type")], names)

    def test_two_spans_sharing_a_description_get_different_names(self) -> None:
        names = [c.name for c in self._suggestions([
            Feature(type="DOMAIN", start=30, end=100, description="Ig-like C2-type"),
            Feature(type="DOMAIN", start=150, end=250, description="Ig-like C2-type"),
        ])]
        domain_names = [n for n in names if n.startswith("ig_like_c2_type")]
        self.assertEqual(len(domain_names), len(set(domain_names)), domain_names)


class FurinCuratedSiteTests(unittest.TestCase):
    """F32."""

    def _sites(self, features):
        from agdesign2.models import Region
        from agdesign2.pipeline import AntigenAnalyzer

        analyzer = AntigenAnalyzer(config=AnalysisConfig(generate_assets=False))
        sequence = "M" * 209 + "RRRR" + "M" * 100
        return analyzer._find_furin_sites(
            sequence,
            Region(start=1, end=len(sequence), label="r", source="t"),
            features,
        )

    def test_a_motif_on_a_curated_cleavage_site_is_labelled(self) -> None:
        sites = self._sites([
            Feature(type="SITE", start=213, end=214, description="Cleavage; by furin and PCSK7"),
        ])
        self.assertTrue(sites)
        self.assertIn("Cleavage; by furin and PCSK7", sites[0].curated_annotation)
        # The edit is still offered: an uncleavable variant is a real design goal.
        self.assertTrue(sites[0].suggestions)

    def test_a_motif_with_no_annotation_is_not_labelled(self) -> None:
        self.assertIsNone(self._sites([])[0].curated_annotation)

    def test_a_propeptide_junction_counts(self) -> None:
        sites = self._sites([
            Feature(type="PROPEP", start=20, end=212, description="Removed in mature form"),
        ])
        self.assertIn("PROPEPTIDE", sites[0].curated_annotation)


class ProjectionDenominatorTests(unittest.TestCase):
    """OA-012, with the abca4_human shape."""

    def test_a_request_larger_than_the_alignment_is_rejected(self) -> None:
        # The alignment covers 100 residues; the construct asks for 500.
        alignment = global_align("A" * 100, "A" * 100)
        start, end, sequence, notes = map_query_region_to_subject(
            alignment, query_start=1, query_end=100,
            subject_sequence="A" * 100, requested_length=500,
        )
        self.assertEqual((start, end, sequence), (None, None, None))
        self.assertIn("100/500", notes[0])

    def test_a_contained_request_still_passes(self) -> None:
        alignment = global_align("A" * 100, "A" * 100)
        start, end, _, notes = map_query_region_to_subject(
            alignment, query_start=1, query_end=100,
            subject_sequence="A" * 100, requested_length=100,
        )
        self.assertEqual((start, end), (1, 100))

    def test_the_default_is_the_old_behaviour(self) -> None:
        alignment = global_align("A" * 100, "A" * 100)
        start, end, _, _ = map_query_region_to_subject(
            alignment, query_start=1, query_end=100, subject_sequence="A" * 100
        )
        self.assertEqual((start, end), (1, 100))


class IsoformFeatureTests(unittest.TestCase):
    """RSM-002, with the PPAP numbers."""

    def test_a_feature_past_the_sequence_end_is_dropped(self) -> None:
        from agdesign2.uniprot import UniProtClient

        entry = {
            "sequence": {"value": "M" * 386},
            "features": [
                {"type": "Chain", "location": {"start": {"value": 33}, "end": {"value": 386}}},
                {"type": "Transmembrane", "location": {"start": {"value": 382}, "end": {"value": 404}}},
                {"type": "Topological domain", "location": {"start": {"value": 405}, "end": {"value": 418}}},
            ],
        }
        features = UniProtClient.get_features(UniProtClient.__new__(UniProtClient), entry)
        self.assertEqual([f.type for f in features], ["CHAIN"])


class SingleTransmembraneNoteTests(unittest.TestCase):
    """F10b — the note must not contradict the classification on the page."""

    def test_the_note_describes_what_happened(self) -> None:
        features = [Feature(type="TRANSMEM", start=s, end=s + 20, description="Helical")
                    for s in (50, 100, 150, 200, 250, 300, 350)]
        result = derive_ectodomain(features, 500, AnalysisConfig())
        messages = [note.message for note in result.notes]
        self.assertFalse([m for m in messages if "Expected exactly one" in m], messages)
        self.assertTrue([m for m in messages if "no single extracellular domain" in m], messages)
        self.assertFalse(
            [n for n in result.notes if n.severity == "warning" and n.source == "topology"],
            "the note is not a warning about the report",
        )


if __name__ == "__main__":
    unittest.main()
