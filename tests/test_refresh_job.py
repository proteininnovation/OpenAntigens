from __future__ import annotations

import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from agdesign2.refresh_job import run_accessible_portal_refresh_job
from agdesign2.target_sets import DEFAULT_ACCESSIBLE_BUCKETS, build_accessible_target_tsv

from tests.test_batch import FakeAnalyzer


class TargetSetTests(unittest.TestCase):
    def test_build_accessible_target_tsv_filters_and_deduplicates(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            source_csv = tmpdir_path / "accessible.csv"
            output_tsv = tmpdir_path / "targets.tsv"
            with source_csv.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=[
                        "accession",
                        "uniprot_entry_name",
                        "gene_symbol",
                        "protein_name",
                        "primary_topology_class",
                        "simplified_topology_bucket",
                    ],
                )
                writer.writeheader()
                writer.writerow(
                    {
                        "accession": "P1",
                        "uniprot_entry_name": "AAA_HUMAN",
                        "gene_symbol": "AAA",
                        "protein_name": "A",
                        "primary_topology_class": "secreted_soluble",
                        "simplified_topology_bucket": "Secreted",
                    }
                )
                writer.writerow(
                    {
                        "accession": "P1",
                        "uniprot_entry_name": "AAA_HUMAN",
                        "gene_symbol": "AAA",
                        "protein_name": "A duplicate",
                        "primary_topology_class": "secreted_soluble",
                        "simplified_topology_bucket": "Secreted",
                    }
                )
                writer.writerow(
                    {
                        "accession": "P2",
                        "uniprot_entry_name": "BBB_HUMAN",
                        "gene_symbol": "BBB",
                        "protein_name": "B",
                        "primary_topology_class": "multipass_surface",
                        "simplified_topology_bucket": "Multipass",
                    }
                )
            rows = build_accessible_target_tsv(
                source_csv=source_csv,
                output_tsv=output_tsv,
                allowed_buckets=("Secreted", "GPI", "Single-pass"),
            )
            with output_tsv.open("r", encoding="utf-8", newline="") as handle:
                written_rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["uniprot_name"], "AAA_HUMAN")
        self.assertEqual(rows[0]["gene_symbol"], "AAA")
        self.assertEqual(len(written_rows), 1)

    def test_default_accessible_target_tsv_includes_multipass(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            source_csv = tmpdir_path / "accessible.csv"
            output_tsv = tmpdir_path / "targets.tsv"
            with source_csv.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=[
                        "accession",
                        "uniprot_entry_name",
                        "gene_symbol",
                        "protein_name",
                        "primary_topology_class",
                        "simplified_topology_bucket",
                    ],
                )
                writer.writeheader()
                writer.writerow(
                    {
                        "accession": "P2",
                        "uniprot_entry_name": "BBB_HUMAN",
                        "gene_symbol": "BBB",
                        "protein_name": "B",
                        "primary_topology_class": "multipass_surface",
                        "simplified_topology_bucket": "Multipass",
                    }
                )
            rows = build_accessible_target_tsv(source_csv=source_csv, output_tsv=output_tsv)
        self.assertEqual(DEFAULT_ACCESSIBLE_BUCKETS, ("Secreted", "GPI", "Single-pass", "Multipass"))
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["topology_bucket"], "Multipass")


class RefreshJobTests(unittest.TestCase):
    def test_refresh_job_resumes_from_state(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            source_csv = tmpdir_path / "accessible.csv"
            with source_csv.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=[
                        "accession",
                        "uniprot_entry_name",
                        "gene_symbol",
                        "protein_name",
                        "primary_topology_class",
                        "simplified_topology_bucket",
                    ],
                )
                writer.writeheader()
                writer.writerow(
                    {
                        "accession": "P1",
                        "uniprot_entry_name": "AAA_HUMAN",
                        "gene_symbol": "AAA",
                        "protein_name": "A",
                        "primary_topology_class": "secreted_soluble",
                        "simplified_topology_bucket": "Secreted",
                    }
                )
                writer.writerow(
                    {
                        "accession": "P2",
                        "uniprot_entry_name": "BBB_HUMAN",
                        "gene_symbol": "BBB",
                        "protein_name": "B",
                        "primary_topology_class": "single_pass_surface",
                        "simplified_topology_bucket": "Single-pass",
                    }
                )

            output_tsv = tmpdir_path / "targets.tsv"
            output_dir = tmpdir_path / "out"
            state_path = tmpdir_path / "state.json"

            with mock.patch("agdesign2.refresh_job.AntigenAnalyzer", return_value=FakeAnalyzer()):
                with mock.patch("agdesign2.refresh_job.build_portal") as build_portal_mock:
                    first = run_accessible_portal_refresh_job(
                        source_csv=source_csv,
                        output_tsv=output_tsv,
                        output_dir=output_dir,
                        state_path=state_path,
                        max_targets=1,
                    )
                    second = run_accessible_portal_refresh_job(
                        source_csv=source_csv,
                        output_tsv=output_tsv,
                        output_dir=output_dir,
                        state_path=state_path,
                    )

            state = json.loads(state_path.read_text(encoding="utf-8"))
            self.assertEqual(first.processed_this_run, 1)
            self.assertEqual(second.completed_targets, 2)
            self.assertEqual([item["status"] for item in state["items"]], ["ok", "ok"])
            self.assertTrue((output_dir / "aaa_human_report.json").exists())
            self.assertTrue((output_dir / "bbb_human_report.json").exists())
            self.assertTrue(build_portal_mock.called)

    def test_refresh_job_migrates_state_when_new_targets_are_added(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            source_csv = tmpdir_path / "accessible.csv"
            fieldnames = [
                "accession",
                "uniprot_entry_name",
                "gene_symbol",
                "protein_name",
                "primary_topology_class",
                "simplified_topology_bucket",
            ]
            with source_csv.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerow(
                    {
                        "accession": "P1",
                        "uniprot_entry_name": "AAA_HUMAN",
                        "gene_symbol": "AAA",
                        "protein_name": "A",
                        "primary_topology_class": "secreted_soluble",
                        "simplified_topology_bucket": "Secreted",
                    }
                )
            output_tsv = tmpdir_path / "targets.tsv"
            output_dir = tmpdir_path / "out"
            state_path = tmpdir_path / "state.json"
            with mock.patch("agdesign2.refresh_job.AntigenAnalyzer", return_value=FakeAnalyzer()):
                run_accessible_portal_refresh_job(
                    source_csv=source_csv,
                    output_tsv=output_tsv,
                    output_dir=output_dir,
                    state_path=state_path,
                    max_targets=1,
                    allowed_buckets=("Secreted",),
                )

            with source_csv.open("a", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fieldnames)
                writer.writerow(
                    {
                        "accession": "P2",
                        "uniprot_entry_name": "BBB_HUMAN",
                        "gene_symbol": "BBB",
                        "protein_name": "B",
                        "primary_topology_class": "multipass_surface",
                        "simplified_topology_bucket": "Multipass",
                    }
                )
            with mock.patch("agdesign2.refresh_job.AntigenAnalyzer", return_value=FakeAnalyzer()):
                run_accessible_portal_refresh_job(
                    source_csv=source_csv,
                    output_tsv=output_tsv,
                    output_dir=output_dir,
                    state_path=state_path,
                    max_targets=0,
                    allowed_buckets=("Secreted", "Multipass"),
                )

            state = json.loads(state_path.read_text(encoding="utf-8"))
        self.assertEqual([item["query"] for item in state["items"]], ["AAA_HUMAN", "BBB_HUMAN"])
        self.assertEqual([item["status"] for item in state["items"]], ["ok", "pending"])
        self.assertEqual([item["batch_index"] for item in state["items"]], [1, 2])


class _AlwaysErrorAnalyzer:
    def __init__(self) -> None:
        self.calls: list[str] = []

    def analyze_target(self, query, output_dir=None):
        self.calls.append(query)
        raise RuntimeError("boom")


class RefreshJobStateAndRetryTests(unittest.TestCase):
    def test_reconcile_resets_leftover_running_to_pending(self) -> None:
        from agdesign2.refresh_job import _reconcile_state_with_existing_reports

        class Row:
            def __init__(self, query, column):
                self.query = query
                self.source_column = column
                self.record = {"gene_symbol": query}

        state = {
            "items": [
                {
                    "query": "AAA",
                    "source_column": "gene_symbol",
                    "status": "running",
                    "json_report": "stale.json",
                    "markdown_report": "stale.md",
                }
            ]
        }
        with tempfile.TemporaryDirectory() as tmp:
            _reconcile_state_with_existing_reports(
                state, rows=[Row("AAA", "gene_symbol")], output_root=Path(tmp)
            )
        self.assertEqual(state["items"][0]["status"], "pending")
        self.assertIsNone(state["items"][0]["json_report"])

    def _write_single_row_csv(self, path: Path) -> None:
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "accession", "uniprot_entry_name", "gene_symbol", "protein_name",
                    "primary_topology_class", "simplified_topology_bucket",
                ],
            )
            writer.writeheader()
            writer.writerow({
                "accession": "P1", "uniprot_entry_name": "AAA_HUMAN", "gene_symbol": "AAA",
                "protein_name": "A", "primary_topology_class": "secreted_soluble",
                "simplified_topology_bucket": "Secreted",
            })

    def test_max_attempts_stops_retrying_failed_target(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            source_csv = tmpdir_path / "accessible.csv"
            self._write_single_row_csv(source_csv)
            analyzer = _AlwaysErrorAnalyzer()
            with mock.patch("agdesign2.refresh_job.AntigenAnalyzer", return_value=analyzer):
                with mock.patch("agdesign2.refresh_job.build_portal"):
                    common = dict(
                        source_csv=source_csv,
                        output_tsv=tmpdir_path / "targets.tsv",
                        output_dir=tmpdir_path / "out",
                        state_path=tmpdir_path / "state.json",
                        max_attempts=1,
                    )
                    run_accessible_portal_refresh_job(**common)
                    run_accessible_portal_refresh_job(**common)
            # First run errors once (attempts=1); second run must skip the
            # exhausted target rather than re-invoking analyze_target.
            self.assertEqual(analyzer.calls, ["AAA_HUMAN"])
