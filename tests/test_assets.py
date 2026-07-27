from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from agdesign2.assets import (
    _remove_disabled_asset_outputs,
    _report_assets_complete,
    _reuse_duplicate_boundary_quality_plots,
    generate_assets_for_report_json,
)
from agdesign2.config import AnalysisConfig
from agdesign2.structure_utils import pae_matches_sequence_length


class AssetResumeTests(unittest.TestCase):
    def test_resume_collapses_existing_duplicate_boundary_plots(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            report_path = root / "test_human_report.json"
            assets_dir = root / "test_human_report_assets"
            assets_dir.mkdir()
            primary = assets_dir / "pdb_1abc_quality.png"
            duplicate = assets_dir / "pdb_2def_quality.png"
            primary.write_bytes(b"primary")
            duplicate.write_bytes(b"duplicate")
            payload = {
                "target": {"entry_name": "TEST_HUMAN"},
                "construct_details": [
                    {
                        "name": "pdb_1abc",
                        "start": 2,
                        "end": 5,
                        "quality_plot": "test_human_report_assets/pdb_1abc_quality.png",
                    },
                    {
                        "name": "pdb_2def",
                        "start": 2,
                        "end": 5,
                        "quality_plot": "test_human_report_assets/pdb_2def_quality.png",
                    },
                ],
            }

            changed = _reuse_duplicate_boundary_quality_plots(
                payload=payload,
                report_path=report_path,
            )

            refs = {
                construct["quality_plot"]
                for construct in payload["construct_details"]
            }
            self.assertTrue(changed)
            self.assertEqual(
                refs,
                {"test_human_report_assets/pdb_1abc_quality.png"},
            )
            self.assertTrue(primary.exists())
            self.assertFalse(duplicate.exists())

    def test_incompatible_cached_model_removes_stale_assets(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            data_dir = root / "data"
            alphafold_dir = data_dir / "alphafold"
            alphafold_dir.mkdir(parents=True)
            report_path = root / "test_human_report.json"
            assets_dir = root / "test_human_report_assets"
            assets_dir.mkdir()
            quality_path = assets_dir / "full_ectodomain_quality.png"
            quality_path.write_bytes(b"png")
            payload = {
                "target": {
                    "accession": "P00001",
                    "entry_name": "TEST_HUMAN",
                    "gene_symbol": "TEST",
                    "protein_name": "Test",
                    "organism": "Homo sapiens",
                    "taxon_id": 9606,
                    "sequence": "CC",
                },
                "resolution": {},
                "features": [],
                "ectodomain": {"start": 1, "end": 2, "label": "test", "source": "test"},
                "structural_regions": [],
                "experimental_constructs": [],
                "construct_recommendations": [],
                "furin_sites": [],
                "cysteine_analysis": [],
                "assembly_requirements": [],
                "interpro_annotations": [],
                "species_name_matches": [],
                "ectodomain_homology": [],
                "family_context": None,
                "cross_reactivity_hits": [],
                "construct_details": [
                    {
                        "name": "full_ectodomain",
                        "start": 1,
                        "end": 2,
                        "length": 2,
                        "sequence": "CC",
                        "score": 1.0,
                        "rationale": "test",
                        "quality_plot": "test_human_report_assets/full_ectodomain_quality.png",
                    }
                ],
            }
            report_path.write_text(json.dumps(payload), encoding="utf-8")
            (alphafold_dir / "P00001.pdb").write_text(
                "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 90.00           C  \n",
                encoding="utf-8",
            )
            (alphafold_dir / "P00001.pae.json").write_text(
                '{"predicted_aligned_error":[[0.0,0.0],[0.0,0.0]]}',
                encoding="utf-8",
            )

            result = generate_assets_for_report_json(
                report_path,
                config=AnalysisConfig(
                    data_dir=data_dir,
                    render_structure_images=False,
                    render_quality_plots=True,
                ),
                resume=True,
            )

            updated = json.loads(report_path.read_text(encoding="utf-8"))
            self.assertEqual(result["status"], "skipped_incompatible_alphafold")
            self.assertIsNone(updated["construct_details"][0]["quality_plot"])
            self.assertFalse(quality_path.exists())

    def test_pae_dimensions_must_match_canonical_sequence(self) -> None:
        self.assertTrue(pae_matches_sequence_length([[0.0, 1.0], [1.0, 0.0]], 2))
        self.assertFalse(pae_matches_sequence_length([[0.0]], 2))
        self.assertFalse(pae_matches_sequence_length([[0.0, 1.0]], 1))

    def test_no_structure_mode_removes_stale_pymol_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            report_path = root / "test_human_report.json"
            assets_dir = root / "test_human_report_assets"
            assets_dir.mkdir()
            structure_path = assets_dir / "full_ectodomain_structure.png"
            structure_path.write_bytes(b"png")
            payload = {
                "target": {"entry_name": "TEST_HUMAN"},
                "construct_details": [
                    {
                        "name": "full ectodomain",
                        "start": 1,
                        "end": 100,
                        "structure_image": "test_human_report_assets/full_ectodomain_structure.png",
                    }
                ],
            }

            changed = _remove_disabled_asset_outputs(
                payload=payload,
                report_path=report_path,
                config=AnalysisConfig(
                    render_structure_images=False,
                    render_quality_plots=True,
                ),
            )

            self.assertTrue(changed)
            self.assertIsNone(payload["construct_details"][0]["structure_image"])
            self.assertFalse(structure_path.exists())

    def test_missing_requested_structure_images_are_not_treated_as_complete(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            report_path = root / "test_human_report.json"
            assets_dir = root / "test_human_report_assets"
            assets_dir.mkdir()
            (assets_dir / "full_ectodomain_quality.png").write_bytes(b"png")
            payload = {
                "target": {"entry_name": "TEST_HUMAN"},
                "ectodomain": {"start": 1, "end": 100},
                "construct_details": [
                    {
                        "name": "full ectodomain",
                        "start": 1,
                        "end": 100,
                        "quality_plot": "test_human_report_assets/full_ectodomain_quality.png",
                        "structure_image": None,
                    }
                ],
            }

            quality_only = AnalysisConfig(render_structure_images=False, render_quality_plots=True)
            with_structures = AnalysisConfig(render_structure_images=True, render_quality_plots=True)

            self.assertTrue(_report_assets_complete(payload=payload, report_path=report_path, config=quality_only))
            self.assertFalse(_report_assets_complete(payload=payload, report_path=report_path, config=with_structures))


if __name__ == "__main__":
    unittest.main()
