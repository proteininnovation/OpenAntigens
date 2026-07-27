from __future__ import annotations

import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from matplotlib.axes import Axes

from agdesign2.models import ConstructDetail, Region
from agdesign2.visualization import ConstructAssetGenerator, _topology_intervals


class VisualizationTests(unittest.TestCase):
    @patch.dict(os.environ, {}, clear=True)
    def test_resolves_common_pymol_location(self) -> None:
        def fake_which(candidate: str) -> str | None:
            if candidate == "/opt/homebrew/bin/pymol":
                return candidate
            return None

        with patch("agdesign2.visualization.shutil.which", side_effect=fake_which):
            generator = ConstructAssetGenerator()

        self.assertEqual(generator.pymol_executable, "/opt/homebrew/bin/pymol")
        with patch("agdesign2.visualization.shutil.which", side_effect=fake_which):
            self.assertTrue(generator.is_pymol_available())

    def test_topology_intervals_group_residue_annotations(self) -> None:
        annotations = [
            {"position": 1, "topology_location": "extracellular"},
            {"position": 2, "topology_location": "extracellular"},
            {"position": 3, "topology_location": "membrane"},
            {"position": 4, "topology_location": "intracellular"},
        ]

        self.assertEqual(
            _topology_intervals(annotations, 1, 4),
            [(1, 2, "extracellular"), (3, 3, "membrane"), (4, 4, "intracellular")],
        )

    def test_quality_plot_accepts_topology_annotations(self) -> None:
        generator = ConstructAssetGenerator(pymol_executable="missing-pymol")
        residue_plddt = {index: 80.0 for index in range(1, 7)}
        pae_matrix = [[1.0 if row == col else 5.0 for col in range(6)] for row in range(6)]
        annotations = [
            {"position": 1, "topology_location": "extracellular"},
            {"position": 2, "topology_location": "extracellular"},
            {"position": 3, "topology_location": "membrane"},
            {"position": 4, "topology_location": "membrane"},
            {"position": 5, "topology_location": "intracellular"},
            {"position": 6, "topology_location": "intracellular"},
        ]

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "quality.png"
            generator.render_quality_plot(
                residue_plddt=residue_plddt,
                pae_matrix=pae_matrix,
                ectodomain=Region(start=1, end=6, label="test", source="test"),
                construct=ConstructDetail(
                    name="test_construct",
                    start=2,
                    end=5,
                    length=4,
                    sequence="AAAA",
                    score=1.0,
                    rationale="test",
                ),
                residue_annotations=annotations,
                output_path=output_path,
            )

            self.assertTrue(os.path.exists(output_path))
            self.assertGreater(os.path.getsize(output_path), 0)

    def test_duplicate_boundaries_share_one_quality_plot(self) -> None:
        generator = ConstructAssetGenerator(pymol_executable="missing-pymol")
        constructs = [
            ConstructDetail(
                name="pdb_1abc",
                start=2,
                end=5,
                length=4,
                sequence="AAAA",
                score=1.0,
                rationale="test",
            ),
            ConstructDetail(
                name="pdb_2def",
                start=2,
                end=5,
                length=4,
                sequence="AAAA",
                score=1.0,
                rationale="test",
            ),
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)

            def render_once(*, output_path: Path, **_kwargs) -> None:
                output_path.write_bytes(b"png")

            with (
                patch.object(generator, "is_pymol_available", return_value=False),
                patch.object(
                    generator,
                    "_render_quality_plot_from_context",
                    side_effect=render_once,
                ) as render_quality,
            ):
                results = generator.generate_assets_for_constructs(
                    constructs=constructs,
                    pdb_path=output_dir / "model.pdb",
                    residue_plddt={index: 80.0 for index in range(1, 7)},
                    pae_matrix=[[1.0] * 6 for _ in range(6)],
                    ectodomain=Region(start=1, end=6, label="test", source="test"),
                    output_dir=output_dir,
                    render_structure_image=False,
                    render_quality_plot=True,
                )

            self.assertEqual(render_quality.call_count, 1)
            self.assertEqual(results["pdb_1abc"], results["pdb_2def"])
            self.assertEqual(len(list(output_dir.glob("*_quality.png"))), 1)

    def test_quality_plot_title_uses_shared_boundary_not_construct_name(self) -> None:
        generator = ConstructAssetGenerator(pymol_executable="missing-pymol")
        construct = ConstructDetail(
            name="pdb_1abc",
            start=2,
            end=5,
            length=4,
            sequence="AAAA",
            score=1.0,
            rationale="test",
        )
        with tempfile.TemporaryDirectory() as tmpdir, patch.object(
            Axes,
            "set_title",
            autospec=True,
        ) as set_title:
            generator.render_quality_plot(
                residue_plddt={index: 80.0 for index in range(1, 7)},
                pae_matrix=[[1.0] * 6 for _ in range(6)],
                ectodomain=Region(start=1, end=6, label="test", source="test"),
                construct=construct,
                output_path=Path(tmpdir) / "quality.png",
            )

        titles = [str(call.args[1]) for call in set_title.call_args_list]
        self.assertTrue(titles)
        self.assertTrue(all("Residues 2-5" in title for title in titles))
        self.assertTrue(all("pdb_1abc" not in title for title in titles))


if __name__ == "__main__":
    unittest.main()
