from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from agdesign2.config import AnalysisConfig
from agdesign2.models import Region
from agdesign2.structure_utils import (
    analyze_cysteines,
    build_pae_stats,
    calculate_residue_sasa,
    classify_surface_accessibility,
    collect_region_split_diagnostics,
    find_structured_regions,
    split_domains_from_pae,
    summarize_construct_quality,
)


class StructureTests(unittest.TestCase):
    def test_cysteine_analysis_reports_nearest_unpaired_cysteine_distance(self) -> None:
        atoms = {
            3: {"CA": (0.0, 0.0, 0.0), "SG": (0.0, 0.0, 0.0)},
            8: {"CA": (8.0, 0.0, 0.0), "SG": (5.1, 0.0, 0.0)},
            20: {"CA": (20.0, 0.0, 0.0), "SG": (20.0, 0.0, 0.0)},
        }
        residue_names = {3: "CYS", 8: "CYS", 20: "CYS"}

        findings = analyze_cysteines(atoms, residue_names, 1, 25)
        by_position = {finding.position: finding for finding in findings}

        self.assertIsNone(by_position[3].paired_with)
        self.assertEqual(by_position[3].closest_unpaired_cysteine, 8)
        self.assertEqual(by_position[3].closest_unpaired_cysteine_distance, 5.1)
        self.assertEqual(by_position[3].closest_unpaired_cysteine_distance_atom, "SG-SG")

    def test_cysteine_analysis_uses_ca_distance_when_sulfur_atoms_are_missing(self) -> None:
        atoms = {
            3: {"CA": (0.0, 0.0, 0.0)},
            8: {"CA": (7.5, 0.0, 0.0)},
        }
        residue_names = {3: "CYS", 8: "CYS"}

        findings = analyze_cysteines(atoms, residue_names, 1, 10)
        by_position = {finding.position: finding for finding in findings}

        self.assertEqual(by_position[3].closest_unpaired_cysteine, 8)
        self.assertEqual(by_position[3].closest_unpaired_cysteine_distance, 7.5)
        self.assertEqual(by_position[3].closest_unpaired_cysteine_distance_atom, "CA-CA")

    def test_finds_structured_regions(self) -> None:
        config = AnalysisConfig(min_structured_segment=3, max_disordered_gap=1)
        plddt = {1: 85.0, 2: 82.0, 3: 81.0, 4: 40.0, 5: 90.0, 6: 91.0, 7: 92.0}
        regions = find_structured_regions(1, 7, plddt, config)
        self.assertEqual([(region.start, region.end) for region in regions], [(1, 7)])

    def test_splits_domains_from_pae(self) -> None:
        config = AnalysisConfig(min_domain_size=2, pae_domain_threshold=10.0, pae_separation_threshold=5.0)
        matrix = [
            [1, 1, 18, 18],
            [1, 1, 18, 18],
            [18, 18, 1, 1],
            [18, 18, 1, 1],
        ]
        regions = split_domains_from_pae(Region(start=1, end=4, label="segment", source="AlphaFold"), matrix, config)
        self.assertEqual([(region.start, region.end) for region in regions], [(1, 2), (3, 4)])

    def test_does_not_split_cohesive_region(self) -> None:
        config = AnalysisConfig(min_domain_size=2, pae_domain_threshold=10.0, pae_separation_threshold=5.0)
        matrix = [
            [1, 3, 5, 5],
            [3, 1, 5, 5],
            [5, 5, 1, 3],
            [5, 5, 3, 1],
        ]
        regions = split_domains_from_pae(Region(start=1, end=4, label="segment", source="AlphaFold"), matrix, config)
        self.assertEqual([(region.start, region.end) for region in regions], [(1, 4)])

    def test_prefers_cut_with_low_confidence_linker(self) -> None:
        config = AnalysisConfig(
            min_domain_size=2,
            pae_domain_threshold=10.0,
            pae_separation_threshold=3.0,
            domain_linker_window=1,
        )
        matrix = [
            [1, 1, 8, 18, 18, 18],
            [1, 1, 8, 18, 18, 18],
            [8, 8, 1, 18, 18, 18],
            [18, 18, 18, 1, 8, 8],
            [18, 18, 18, 8, 1, 1],
            [18, 18, 18, 8, 1, 1],
        ]
        plddt = {1: 90.0, 2: 89.0, 3: 55.0, 4: 58.0, 5: 91.0, 6: 90.0}
        regions = split_domains_from_pae(
            Region(start=1, end=6, label="segment", source="AlphaFold"),
            matrix,
            config,
            residue_plddt=plddt,
        )
        self.assertEqual([(region.start, region.end) for region in regions], [(1, 3), (4, 6)])

    def test_accepts_precomputed_pae_stats(self) -> None:
        config = AnalysisConfig(min_domain_size=2, pae_domain_threshold=10.0, pae_separation_threshold=5.0)
        matrix = [
            [1, 1, 18, 18],
            [1, 1, 18, 18],
            [18, 18, 1, 1],
            [18, 18, 1, 1],
        ]
        stats = build_pae_stats(matrix)
        regions = split_domains_from_pae(Region(start=1, end=4, label="segment", source="AlphaFold"), stats, config)
        self.assertEqual([(region.start, region.end) for region in regions], [(1, 2), (3, 4)])

    def test_construct_quality_and_split_diagnostics_accept_precomputed_pae_stats(self) -> None:
        config = AnalysisConfig(min_domain_size=2, pae_domain_threshold=10.0, pae_separation_threshold=5.0)
        matrix = [
            [1, 1, 18, 18],
            [1, 1, 18, 18],
            [18, 18, 1, 1],
            [18, 18, 1, 1],
        ]
        stats = build_pae_stats(matrix)
        plddt = {1: 90.0, 2: 88.0, 3: 91.0, 4: 89.0}
        metrics = summarize_construct_quality(1, 4, plddt, stats, config)
        diagnostics = collect_region_split_diagnostics(
            Region(start=1, end=4, label="construct", source="Construct"),
            stats,
            config,
            residue_plddt=plddt,
        )
        self.assertEqual(metrics.mean_intra_pae, 12.33)
        self.assertEqual(metrics.max_intra_pae, 18.0)
        self.assertEqual([(item.parent_start, item.boundary_after, item.parent_end) for item in diagnostics], [(1, 2, 4)])

    def test_construct_quality_gracefully_handles_regions_outside_pae_bounds(self) -> None:
        config = AnalysisConfig(min_domain_size=2, pae_domain_threshold=10.0, pae_separation_threshold=5.0)
        matrix = [
            [1, 1, 18, 18],
            [1, 1, 18, 18],
            [18, 18, 1, 1],
            [18, 18, 1, 1],
        ]
        stats = build_pae_stats(matrix)
        plddt = {4: 88.0, 5: 90.0, 6: 87.0}
        metrics = summarize_construct_quality(4, 6, plddt, stats, config)
        diagnostics = collect_region_split_diagnostics(
            Region(start=4, end=6, label="construct", source="Construct"),
            stats,
            config,
            residue_plddt=plddt,
        )
        self.assertEqual(metrics.mean_plddt, 88.33)
        self.assertEqual(metrics.mean_intra_pae, 0.0)
        self.assertEqual(metrics.max_intra_pae, 1.0)
        self.assertEqual(diagnostics, [])

    def test_sasa_increases_when_uncertain_neighboring_block_is_removed(self) -> None:
        config = AnalysisConfig(
            sasa_slices=20,
            sasa_relative_exposed_threshold=0.20,
            sasa_low_confidence_plddt_threshold=50.0,
        )
        atoms = {
            1: {"CA": (0.0, 0.0, 0.0)},
            2: {"CA": (3.0, 0.0, 0.0)},
            3: {"CA": (-3.0, 0.0, 0.0)},
            4: {"CA": (0.0, 3.0, 0.0)},
            5: {"CA": (0.0, -3.0, 0.0)},
            6: {"CA": (0.0, 0.0, 3.0)},
            7: {"CA": (0.0, 0.0, -3.0)},
        }
        residue_names = {index: "ALA" for index in atoms}
        full_sasa = calculate_residue_sasa(
            atoms,
            residue_names,
            residues={1},
            context_residues=set(atoms),
            probe_radius=config.sasa_probe_radius,
            slices=config.sasa_slices,
        )
        block_sasa = calculate_residue_sasa(
            atoms,
            residue_names,
            residues={1},
            context_residues={1},
            probe_radius=config.sasa_probe_radius,
            slices=config.sasa_slices,
        )
        self.assertLess(full_sasa[1], block_sasa[1])
        accessibility = classify_surface_accessibility(
            atoms,
            residue_names,
            {index: 90.0 for index in atoms},
            [Region(start=1, end=1, label="Block A", source="AlphaFold")],
            config,
        )
        self.assertEqual(accessibility[1].classification, "conditional")
        self.assertTrue(accessibility[1].surface_exposed)

    def test_low_plddt_residue_is_accessible_but_low_confidence(self) -> None:
        config = AnalysisConfig(sasa_slices=20, sasa_low_confidence_plddt_threshold=50.0)
        atoms = {
            1: {"CA": (0.0, 0.0, 0.0)},
            2: {"CA": (3.0, 0.0, 0.0)},
            3: {"CA": (-3.0, 0.0, 0.0)},
            4: {"CA": (0.0, 3.0, 0.0)},
            5: {"CA": (0.0, -3.0, 0.0)},
            6: {"CA": (0.0, 0.0, 3.0)},
            7: {"CA": (0.0, 0.0, -3.0)},
        }
        residue_names = {index: "ALA" for index in atoms}
        accessibility = classify_surface_accessibility(
            atoms,
            residue_names,
            {index: 40.0 for index in atoms},
            [Region(start=1, end=1, label="Block A", source="AlphaFold")],
            config,
        )
        self.assertEqual(accessibility[1].classification, "conditional")
        self.assertTrue(accessibility[1].surface_exposed)
        self.assertIn("low-confidence", accessibility[1].reason)

    def test_surface_accessibility_can_target_residue_subset(self) -> None:
        atoms = {
            1: {"CA": (0.0, 0.0, 0.0)},
            2: {"CA": (3.8, 0.0, 0.0)},
            3: {"CA": (7.6, 0.0, 0.0)},
        }
        residue_names = {index: "ALA" for index in atoms}

        accessibility = classify_surface_accessibility(
            atoms,
            residue_names,
            {index: 90.0 for index in atoms},
            [Region(start=1, end=3, label="Block", source="AlphaFold")],
            AnalysisConfig(),
            target_residues={2},
        )

        self.assertEqual(set(accessibility), {2})


if __name__ == "__main__":
    unittest.main()
