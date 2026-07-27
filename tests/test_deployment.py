from __future__ import annotations

import sys
import tempfile
import unittest
from contextlib import contextmanager
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from agdesign2.deployment import (
    FRESH_SNAPSHOT_STEPS,
    _ortholog_table_complete,
    _select_stratified_smoke_rows,
    build_fresh_snapshot,
    prepare_deployment_data,
)
from agdesign2.cli import build_parser


class OrthologResumeCompletenessTests(unittest.TestCase):
    def test_partial_ortholog_table_is_not_treated_as_complete(self) -> None:
        import json

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            target_tsv = root / "targets.tsv"
            target_tsv.write_text("uniprot_name\nAAA_HUMAN\nBBB_HUMAN\nCCC_HUMAN\n", encoding="utf-8")
            out_tsv = root / "ortholog_reference_table.tsv"
            json_path = out_tsv.with_suffix(".json")

            self.assertFalse(_ortholog_table_complete(out_tsv, target_tsv))  # no file
            json_path.write_text(json.dumps([{"input_index": 1}, {"input_index": 2}]), encoding="utf-8")
            self.assertFalse(_ortholog_table_complete(out_tsv, target_tsv))  # 2 of 3 -> resume
            json_path.write_text(json.dumps([{"input_index": i} for i in range(1, 4)]), encoding="utf-8")
            self.assertTrue(_ortholog_table_complete(out_tsv, target_tsv))  # complete


class DeploymentPrecomputeTests(unittest.TestCase):
    def test_fresh_snapshot_disables_structure_images_by_default(self) -> None:
        parser = build_parser()
        default_args = parser.parse_args(["build-fresh-snapshot", "--dry-run"])
        opt_in_args = parser.parse_args(
            ["build-fresh-snapshot", "--structure-images", "--dry-run"]
        )

        self.assertFalse(default_args.render_structure_images)
        self.assertTrue(opt_in_args.render_structure_images)

    def test_build_fresh_snapshot_dry_run_uses_snapshot_local_paths(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            with _chdir(root):
                result = build_fresh_snapshot(
                    source_csv=root / "surface.csv",
                    snapshot_root=root / "snapshots",
                    snapshot_name="test-snapshot",
                    limit=20,
                    uniprot_prefetch_mode="targets",
                    include_mouse=False,
                    dry_run=True,
                )

        self.assertEqual(result.snapshot_dir.name, "test-snapshot")
        self.assertIn("test-snapshot", str(result.target_tsv))
        self.assertIn("test-snapshot", str(result.output_dir))
        self.assertEqual([row["step"] for row in result.timings], list(FRESH_SNAPSHOT_STEPS))
        self.assertTrue(all(row["status"] == "planned" for row in result.timings))
        self.assertIsNone(result.selected_smoke_tsv)
        self.assertIsNone(result.mouse_public_site_dir)

    def test_full_smoke_dry_run_reports_mouse_and_sample_paths(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            with _chdir(root):
                result = build_fresh_snapshot(
                    source_csv=root / "surface.csv",
                    snapshot_root=root / "snapshots",
                    snapshot_name="openantigen-full-smoke-100",
                    sample_size=100,
                    sample_seed=20260530,
                    sample_mode="stratified-random",
                    include_mouse=True,
                    fetch_mouse_alphafold=True,
                    dry_run=True,
                )

        self.assertEqual(result.selected_smoke_tsv.name, "selected_smoke_targets.tsv")
        self.assertEqual(result.full_smoke_manifest_path.name, "full_smoke_manifest.json")
        self.assertEqual(result.mouse_output_dir.name, "mouse_data")
        self.assertEqual(result.mouse_public_site_dir.name, "mouse")
        self.assertEqual([row["step"] for row in result.timings], list(FRESH_SNAPSHOT_STEPS) + ["build-mouse-portal"])

    def test_fresh_snapshot_dry_run_includes_af3_catalog_step_when_requested(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            catalog = root / "af3_catalog.sqlite"
            catalog.write_text("", encoding="utf-8")
            result = build_fresh_snapshot(
                source_csv=root / "surface.csv",
                snapshot_root=root / "snapshots",
                snapshot_name="test-snapshot",
                af3_catalog=catalog,
                include_mouse=False,
                dry_run=True,
            )

        steps = [row["step"] for row in result.timings]
        self.assertEqual(steps[0], "build-target-tsv")
        self.assertEqual(steps[1], "prepare-af3-catalog-artifacts")
        self.assertEqual(steps[2:], list(FRESH_SNAPSHOT_STEPS)[1:])

    def test_fresh_snapshot_dry_run_ignores_unrequested_af3_catalog(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            catalog = root / "data" / "af3_catalog" / "af3_catalog.sqlite"
            catalog.parent.mkdir(parents=True)
            catalog.write_text("", encoding="utf-8")
            with _chdir(root):
                result = build_fresh_snapshot(
                    source_csv=root / "surface.csv",
                    snapshot_root=root / "snapshots",
                    snapshot_name="test-snapshot",
                    dry_run=True,
                )

        steps = [row["step"] for row in result.timings]
        self.assertNotIn("prepare-af3-catalog-artifacts", steps)

    def test_fresh_snapshot_includes_mouse_by_default(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            result = build_fresh_snapshot(
                source_csv=root / "surface.csv",
                snapshot_root=root / "snapshots",
                snapshot_name="test-snapshot",
                dry_run=True,
            )
        steps = [row["step"] for row in result.timings]
        self.assertIn("build-mouse-portal", steps)
        self.assertIsNotNone(result.mouse_public_site_dir)

    def test_fresh_snapshot_rejects_ambiguous_selection_modes(self) -> None:
        with self.assertRaisesRegex(ValueError, "--limit cannot be combined"):
            build_fresh_snapshot(limit=10, sample_size=100, dry_run=True)
        with self.assertRaisesRegex(ValueError, "--sample-mode requires --sample-size"):
            build_fresh_snapshot(sample_mode="stratified-random", dry_run=True)

    def test_stratified_smoke_sampler_is_deterministic_and_unique(self) -> None:
        rows = _sampler_rows()
        first, first_manifest = _select_stratified_smoke_rows(rows, sample_size=12, sample_seed=20260530)
        second, second_manifest = _select_stratified_smoke_rows(rows, sample_size=12, sample_seed=20260530)

        self.assertEqual([row["uniprot_name"] for row in first], [row["uniprot_name"] for row in second])
        self.assertEqual(first_manifest["selected_count"], 12)
        self.assertEqual(second_manifest["sample_seed"], 20260530)
        self.assertEqual(len({row["uniprot_name"] for row in first}), 12)

    def test_stratified_smoke_sampler_fills_when_strata_are_undersized(self) -> None:
        rows = [
            _target_row(index, "Secreted", "secreted protein")
            for index in range(18)
        ] + [
            _target_row(50, "GPI", "gpi anchored protein"),
            _target_row(51, "Multipass", "G protein-coupled receptor"),
        ]

        selected, manifest = _select_stratified_smoke_rows(rows, sample_size=10, sample_seed=7)

        self.assertEqual(len(selected), 10)
        self.assertEqual(len({row["uniprot_name"] for row in selected}), 10)
        self.assertGreater(manifest["random_fill_count"], 0)

    def test_prepare_deployment_data_dry_run_lists_bulk_steps(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            result = prepare_deployment_data(
                target_tsv=root / "targets.tsv",
                summary_path=root / "batch_summary.json",
                output_dir=root / "outputs",
                cache_dir=root / "cache",
                data_dir=root / "data",
                dry_run=True,
            )

        self.assertEqual(
            result.steps_planned,
            [
                "prefetch-uniprot",
                "setup-blast",
                "build-ortholog-table",
                "build-family-alignments",
                "build-paralog-reference",
                "build-open-targets",
                "refresh-cross-reactivity",
                "build-portal",
            ],
        )
        self.assertEqual(result.steps_completed, [])

    def test_prepare_deployment_data_dry_run_honors_skip_flags(self) -> None:
        result = prepare_deployment_data(
            dry_run=True,
            skip_uniprot=True,
            skip_blast_setup=True,
            skip_open_targets=True,
            skip_portal=True,
        )

        self.assertEqual(
            result.steps_planned,
            [
                "build-ortholog-table",
                "build-family-alignments",
                "build-paralog-reference",
                "refresh-cross-reactivity",
            ],
        )

def _sampler_rows() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    buckets = [
        ("Secreted", "secreted cytokine"),
        ("GPI", "gpi anchored receptor"),
        ("Single-pass", "single-pass receptor complex subunit"),
        ("Multipass", "multipass transporter"),
        ("Multipass", "G protein-coupled receptor"),
    ]
    index = 1
    for bucket, protein_name in buckets:
        for _ in range(5):
            rows.append(_target_row(index, bucket, protein_name))
            index += 1
    return rows


def _target_row(index: int, bucket: str, protein_name: str) -> dict[str, str]:
    return {
        "uniprot_id": f"P{index:05d}",
        "uniprot_name": f"TARGET{index}_HUMAN",
        "gene_symbol": f"T{index}",
        "topology_bucket": bucket,
        "primary_topology_class": bucket,
        "protein_name": protein_name,
    }


@contextmanager
def _chdir(path: Path):
    previous = Path.cwd()
    try:
        import os

        os.chdir(path)
        yield
    finally:
        os.chdir(previous)


if __name__ == "__main__":
    unittest.main()
