from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from agdesign2.config import AnalysisConfig
from agdesign2.models import Feature
from agdesign2.topology import apply_surfy_topology_fallback, derive_ectodomain


class TopologyTests(unittest.TestCase):
    def test_derives_n_terminal_ectodomain(self) -> None:
        features = [
            Feature(type="SIGNAL", start=1, end=24, description="Signal peptide"),
            Feature(type="TOPO_DOM", start=25, end=620, description="Extracellular"),
            Feature(type="TRANSMEM", start=621, end=643, description="Helical"),
            Feature(type="TOPO_DOM", start=644, end=700, description="Cytoplasmic"),
        ]
        result = derive_ectodomain(features, 700, AnalysisConfig())
        self.assertIsNotNone(result.ectodomain)
        self.assertEqual((result.ectodomain.start, result.ectodomain.end), (25, 620))

    def test_uses_signal_peptide_as_fallback(self) -> None:
        features = [
            Feature(type="SIGNAL", start=1, end=18, description="Signal peptide"),
            Feature(type="TRANSMEM", start=220, end=240, description="Helical"),
        ]
        result = derive_ectodomain(features, 300, AnalysisConfig())
        self.assertIsNotNone(result.ectodomain)
        self.assertEqual(result.ectodomain.start, 19)
        self.assertTrue(result.notes)

    def test_returns_none_for_non_single_pass(self) -> None:
        features = [
            Feature(type="TRANSMEM", start=10, end=30),
            Feature(type="TRANSMEM", start=50, end=70),
        ]
        result = derive_ectodomain(features, 100, AnalysisConfig())
        self.assertIsNone(result.ectodomain)
        self.assertIsNotNone(result.topology)
        self.assertEqual(result.topology.topology_class, "multipass_compact")

    def test_preserves_intramembrane_regions_without_changing_tm_count(self) -> None:
        features = [
            Feature(type="TOPO_DOM", start=1, end=20, description="Extracellular"),
            Feature(type="INTRAMEMBRANE", start=21, end=35, description="Pore-forming"),
            Feature(type="TOPO_DOM", start=36, end=45, description="Extracellular"),
            Feature(type="TRANSMEM", start=46, end=66, description="Helical"),
            Feature(type="TOPO_DOM", start=67, end=100, description="Cytoplasmic"),
        ]

        result = derive_ectodomain(features, 100, AnalysisConfig())

        self.assertIsNotNone(result.topology)
        self.assertEqual(result.topology.topology_class, "single_pass")
        self.assertEqual(len(result.topology.transmembrane_regions), 1)
        self.assertEqual(len(result.topology.intramembrane_regions), 1)
        self.assertEqual(
            (result.topology.intramembrane_regions[0].start, result.topology.intramembrane_regions[0].end),
            (21, 35),
        )

    def test_recovers_dominant_extracellular_region_from_conflicting_multi_tm_annotations(self) -> None:
        features = [
            Feature(type="SIGNAL", start=1, end=19, description="Signal peptide"),
            Feature(type="TOPO_DOM", start=20, end=847, description="Extracellular"),
            Feature(type="TRANSMEM", start=256, end=276, description="Helical"),
            Feature(type="TOPO_DOM", start=277, end=364, description="Lumenal"),
            Feature(type="TRANSMEM", start=848, end=868, description="Helical"),
            Feature(type="TOPO_DOM", start=869, end=956, description="Cytoplasmic"),
        ]
        result = derive_ectodomain(features, 956, AnalysisConfig())
        self.assertIsNotNone(result.ectodomain)
        self.assertEqual((result.ectodomain.start, result.ectodomain.end), (20, 847))
        self.assertTrue(any("dominant extracellular topological domain" in note.message.lower() for note in result.notes))
        self.assertIsNotNone(result.topology)
        self.assertEqual(result.topology.topology_class, "multipass_mixed")
        self.assertIsNotNone(result.topology.membrane_core_region)
        self.assertEqual((result.topology.membrane_core_region.start, result.topology.membrane_core_region.end), (256, 868))
        self.assertIsNotNone(result.topology.cytoplasmic_tail_region)
        self.assertEqual((result.topology.cytoplasmic_tail_region.start, result.topology.cytoplasmic_tail_region.end), (869, 956))

    def test_derives_secreted_mature_chain_from_signal_peptide(self) -> None:
        features = [
            Feature(type="SIGNAL", start=1, end=24, description="Signal peptide"),
        ]
        result = derive_ectodomain(features, 266, AnalysisConfig())
        self.assertIsNotNone(result.ectodomain)
        self.assertEqual((result.ectodomain.start, result.ectodomain.end), (25, 266))
        self.assertEqual(result.ectodomain.label, "Secreted mature chain")

    def test_gpi_anchored_target_uses_processed_chain_boundary(self) -> None:
        features = [
            Feature(type="SIGNAL", start=1, end=18, description="Signal peptide"),
            Feature(type="CHAIN", start=19, end=284, description="Carbonic anhydrase 4"),
            Feature(type="LIPIDATION", start=284, end=284, description="GPI-anchor amidated serine"),
        ]
        result = derive_ectodomain(features, 312, AnalysisConfig())
        self.assertEqual((result.ectodomain.start, result.ectodomain.end), (19, 284))
        self.assertEqual(result.ectodomain.label, "GPI-anchored mature chain")
        self.assertEqual(result.topology.topology_class, "gpi_anchored")

    def test_gpi_anchored_target_uses_full_processed_chain_envelope(self) -> None:
        features = [
            Feature(type="SIGNAL", start=1, end=24, description="Signal peptide"),
            Feature(type="CHAIN", start=25, end=358, description="Alpha chain"),
            Feature(type="CHAIN", start=359, end=554, description="Beta chain"),
            Feature(type="LIPIDATION", start=554, end=554, description="GPI-anchor amidated serine"),
        ]
        result = derive_ectodomain(features, 580, AnalysisConfig())
        self.assertEqual((result.ectodomain.start, result.ectodomain.end), (25, 554))

    def test_gpi_anchor_without_processed_chain_retains_untrimmed_construct(self) -> None:
        features = [
            Feature(type="SIGNAL", start=1, end=18, description="Signal peptide"),
            Feature(type="LIPIDATION", start=59, end=59, description="GPI-anchor amidated glycine"),
        ]
        result = derive_ectodomain(features, 80, AnalysisConfig())
        self.assertEqual((result.ectodomain.start, result.ectodomain.end), (19, 80))
        self.assertEqual(result.ectodomain.label, "GPI-anchored sequence (untrimmed)")
        self.assertEqual(result.topology.topology_class, "gpi_anchored_untrimmed")
        self.assertTrue(any("without trimming" in note.message for note in result.notes))

    def test_transmembrane_topology_takes_precedence_over_gpi_processing(self) -> None:
        features = [
            Feature(type="TRANSMEM", start=21, end=43, description="Helical"),
            Feature(type="TOPO_DOM", start=44, end=160, description="Extracellular"),
            Feature(type="CHAIN", start=1, end=161, description="Tetherin"),
            Feature(type="LIPIDATION", start=161, end=161, description="GPI-anchor amidated serine"),
        ]
        result = derive_ectodomain(features, 180, AnalysisConfig())
        self.assertEqual((result.ectodomain.start, result.ectodomain.end), (44, 180))
        self.assertEqual(result.topology.topology_class, "single_pass")

    def test_derives_secreted_chain_from_extracellular_topology_without_tm(self) -> None:
        features = [
            Feature(type="SIGNAL", start=1, end=18, description="Signal peptide"),
            Feature(type="TOPO_DOM", start=19, end=180, description="Extracellular"),
        ]
        result = derive_ectodomain(features, 180, AnalysisConfig())
        self.assertIsNotNone(result.ectodomain)
        self.assertEqual((result.ectodomain.start, result.ectodomain.end), (19, 180))
        self.assertEqual(result.ectodomain.label, "Secreted extracellular chain")

    def test_derives_major_extracellular_region_for_multipass_mixed_target(self) -> None:
        features = [
            Feature(type="SIGNAL", start=1, end=24, description="Signal peptide"),
            Feature(type="TOPO_DOM", start=25, end=180, description="Extracellular"),
            Feature(type="TRANSMEM", start=181, end=201, description="Helical"),
            Feature(type="TOPO_DOM", start=202, end=220, description="Cytoplasmic"),
            Feature(type="TRANSMEM", start=221, end=241, description="Helical"),
            Feature(type="TOPO_DOM", start=242, end=260, description="Cytoplasmic"),
            Feature(type="TRANSMEM", start=261, end=281, description="Helical"),
            Feature(type="TOPO_DOM", start=282, end=330, description="Cytoplasmic"),
        ]
        result = derive_ectodomain(features, 330, AnalysisConfig())
        self.assertIsNotNone(result.ectodomain)
        self.assertEqual((result.ectodomain.start, result.ectodomain.end), (25, 180))
        self.assertIsNotNone(result.topology)
        self.assertEqual(result.topology.topology_class, "multipass_mixed")
        self.assertEqual((result.topology.major_extracellular_region.start, result.topology.major_extracellular_region.end), (25, 180))
        self.assertEqual((result.topology.membrane_core_region.start, result.topology.membrane_core_region.end), (181, 281))
        self.assertEqual((result.topology.cytoplasmic_tail_region.start, result.topology.cytoplasmic_tail_region.end), (282, 330))

    def test_multipass_requires_large_extracellular_region_for_mixed_case(self) -> None:
        features = [
            Feature(type="SIGNAL", start=1, end=20, description="Signal peptide"),
            Feature(type="TOPO_DOM", start=21, end=99, description="Extracellular"),
            Feature(type="TRANSMEM", start=100, end=120, description="Helical"),
            Feature(type="TOPO_DOM", start=121, end=140, description="Cytoplasmic"),
            Feature(type="TRANSMEM", start=141, end=161, description="Helical"),
            Feature(type="TOPO_DOM", start=162, end=220, description="Cytoplasmic"),
        ]
        result = derive_ectodomain(features, 220, AnalysisConfig())
        self.assertIsNone(result.ectodomain)
        self.assertIsNotNone(result.topology)
        self.assertEqual(result.topology.topology_class, "multipass_compact")

    def test_surfy_topology_fallback_fills_missing_multipass_loops(self) -> None:
        features = [
            Feature(type="TRANSMEM", start=36, end=53, description="Helical"),
            Feature(type="TRANSMEM", start=221, end=242, description="Helical"),
        ]
        result = derive_ectodomain(features, 300, AnalysisConfig())
        self.assertIsNotNone(result.topology)
        self.assertEqual(result.topology.extracellular_regions, [])
        changed = apply_surfy_topology_fallback(
            result.topology,
            "CY:1-35;TM:36-53;NC:54-220;TM:221-242;CY:243-261",
            sequence_length=300,
        )
        self.assertTrue(changed)
        self.assertEqual((result.topology.extracellular_regions[0].start, result.topology.extracellular_regions[0].end), (54, 220))
        self.assertEqual((result.topology.cytoplasmic_regions[0].start, result.topology.cytoplasmic_regions[0].end), (1, 35))
        self.assertEqual(result.topology.extracellular_regions[0].source, "SURFY")


if __name__ == "__main__":
    unittest.main()
