from __future__ import annotations

import json
import importlib
import sys
import tempfile
import threading
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from agdesign2.open_targets import (
    OpenTargetsClient,
    build_open_targets_associations_from_bulk_downloads,
    build_open_targets_associations_from_summary,
    download_open_targets_bulk_datasets,
    resolve_open_targets_release,
)
from agdesign2.exceptions import ExternalServiceError


class FakeOpenTargetsClient:
    def resolve_target(self, gene_symbol: str):
        if gene_symbol == "EGFR":
            return {"ensembl_id": "ENSG00000146648", "approved_symbol": "EGFR", "approved_name": "EGFR"}
        return None

    def fetch_associated_diseases(self, ensembl_id: str, *, page_size: int, max_pages: int):
        return 2, [
            {
                "score": 0.91,
                "disease": {"id": "EFO_0000616", "name": "neoplasm"},
                "datasourceScores": [{"id": "chembl", "score": 0.8}],
            },
            {
                "score": 0.72,
                "disease": {"id": "MONDO_0007254", "name": "breast cancer"},
                "datasourceScores": [{"id": "europepmc", "score": 0.6}],
            },
            {
                "score": 0.95,
                "disease": {"id": "EFO_0001444", "name": "measurement"},
                "datasourceScores": [{"id": "europepmc", "score": 0.9}],
            },
            {
                "score": 0.96,
                "disease": {"id": "OBA_2050465", "name": "level of alpha-1B-glycoprotein in blood"},
                "datasourceScores": [{"id": "europepmc", "score": 0.9}],
            },
        ]


class OpenTargetsBuildTests(unittest.TestCase):
    def test_graphql_rejects_non_http_endpoint(self) -> None:
        client = OpenTargetsClient(endpoint="file:///etc/passwd")
        with self.assertRaisesRegex(ExternalServiceError, "Unsupported Open Targets URL"):
            client.graphql("query { test }")

    def test_bulk_downloads_use_bounded_parallel_workers(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            started: list[str] = []
            two_started = threading.Event()

            def fake_download(url: str, destination: Path) -> None:
                started.append(url)
                if len(started) >= 2:
                    two_started.set()
                self.assertTrue(two_started.wait(timeout=1.0))
                destination.write_bytes(b"parquet")

            with (
                mock.patch(
                    "agdesign2.open_targets.OPEN_TARGETS_BULK_DATASETS",
                    ("target",),
                ),
                mock.patch(
                    "agdesign2.open_targets._list_open_targets_parquet_urls",
                    return_value=[
                        "https://example.test/part-00000.parquet",
                        "https://example.test/part-00001.parquet",
                    ],
                ),
                mock.patch(
                    "agdesign2.open_targets._download_file",
                    side_effect=fake_download,
                ),
            ):
                download_open_targets_bulk_datasets(
                    release="26.06",
                    ftp_base_url="https://example.test",
                    data_dir=root,
                    jobs=2,
                )

            self.assertEqual(len(started), 2)
            self.assertTrue((root / "target" / "part-00000.parquet").is_file())
            self.assertTrue((root / "target" / "part-00001.parquet").is_file())

    def test_builds_association_artifacts_from_batch_summary(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            report_json = root / "egfr_human_report.json"
            report_json.write_text(
                json.dumps(
                    {
                        "target": {
                            "accession": "P00533",
                            "entry_name": "EGFR_HUMAN",
                            "gene_symbol": "EGFR",
                            "protein_name": "Epidermal growth factor receptor",
                        }
                    }
                ),
                encoding="utf-8",
            )
            summary_path = root / "batch_summary.json"
            summary_path.write_text(
                json.dumps(
                    [
                        {
                            "query": "EGFR_HUMAN",
                            "status": "ok",
                            "batch_index": 1,
                            "json_report": str(report_json),
                        }
                    ]
                ),
                encoding="utf-8",
            )
            result = build_open_targets_associations_from_summary(
                summary_path,
                client=FakeOpenTargetsClient(),
                page_size=50,
                max_pages=1,
                verbose=False,
            )
            payload = json.loads(result.output_json.read_text(encoding="utf-8"))
            tsv = result.output_tsv.read_text(encoding="utf-8")
        self.assertEqual(payload["target_count"], 1)
        self.assertEqual(payload["targets"][0]["ensembl_id"], "ENSG00000146648")
        self.assertEqual(payload["targets"][0]["top_disease_name"], "neoplasm")
        self.assertIn("breast cancer", tsv)
        self.assertNotIn("measurement", tsv)
        self.assertNotIn("level of alpha-1B-glycoprotein", tsv)
        self.assertIn("chembl:0.8", tsv)

    def test_builds_association_artifacts_from_bulk_tables(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            report_json = root / "egfr_human_report.json"
            report_json.write_text(
                json.dumps(
                    {
                        "target": {
                            "accession": "P00533",
                            "entry_name": "EGFR_HUMAN",
                            "gene_symbol": "EGFR",
                        }
                    }
                ),
                encoding="utf-8",
            )
            summary_path = root / "batch_summary.json"
            summary_path.write_text(
                json.dumps(
                    [
                        {
                            "query": "EGFR_HUMAN",
                            "status": "ok",
                            "batch_index": 1,
                            "json_report": str(report_json),
                        }
                    ]
                ),
                encoding="utf-8",
            )

            def fake_read(path, *, columns, filter_target_ids=None, column_aliases=None):
                path_text = str(path)
                if path_text.endswith("target"):
                    return [{"id": "ENSG00000146648", "approvedSymbol": "EGFR", "approvedName": "EGFR"}]
                if path_text.endswith("association_overall_direct"):
                    self.assertEqual(filter_target_ids, {"ENSG00000146648"})
                    self.assertEqual(column_aliases, {"score": ["score", "associationScore"]})
                    return [
                        {"targetId": "ENSG00000146648", "diseaseId": "EFO_0000616", "score": 0.91},
                        {"targetId": "ENSG00000146648", "diseaseId": "MONDO_0007254", "score": 0.72},
                        {"targetId": "ENSG00000146648", "diseaseId": "EFO_0001444", "score": 0.99},
                        {"targetId": "ENSG00000146648", "diseaseId": "OBA_2050465", "score": 0.98},
                    ]
                if path_text.endswith("association_by_datasource_indirect"):
                    self.assertEqual(filter_target_ids, {"ENSG00000146648"})
                    self.assertEqual(
                        column_aliases,
                        {
                            "score": ["score", "associationScore"],
                            "datasourceId": ["datasourceId", "aggregationValue"],
                        },
                    )
                    return [
                        {
                            "targetId": "ENSG00000146648",
                            "diseaseId": "EFO_0000616",
                            "datasourceId": "europepmc",
                            "score": 0.88,
                        },
                        {
                            "targetId": "ENSG00000146648",
                            "diseaseId": "MONDO_0007254",
                            "datasourceId": "eva",
                            "score": 0.93,
                        },
                        {
                            "targetId": "ENSG00000146648",
                            "diseaseId": "EFO_0001444",
                            "datasourceId": "gwas_credible_sets",
                            "score": 0.99,
                        },
                        {
                            "targetId": "ENSG00000146648",
                            "diseaseId": "OBA_2050465",
                            "datasourceId": "gwas_credible_sets",
                            "score": 0.98,
                        },
                    ]
                raise AssertionError(path_text)

            def fake_disease_read(path):
                return [
                    {"id": "EFO_0000616", "name": "neoplasm", "ancestors": ["EFO_0000408"]},
                    {"id": "MONDO_0007254", "name": "breast cancer", "ancestors": []},
                    {"id": "EFO_0001444", "name": "measurement", "ancestors": []},
                    {"id": "OBA_2050465", "name": "level of alpha-1B-glycoprotein in blood", "ancestors": []},
                ]

            with mock.patch("agdesign2.open_targets._read_parquet_dataset", side_effect=fake_read), mock.patch(
                "agdesign2.open_targets._read_open_targets_disease_dataset", side_effect=fake_disease_read
            ):
                result = build_open_targets_associations_from_bulk_downloads(
                    summary_path,
                    release="26.03",
                    data_dir=root / "open_targets",
                    download=False,
                    verbose=False,
                )
            payload = json.loads(result.output_json.read_text(encoding="utf-8"))
            tsv = result.output_tsv.read_text(encoding="utf-8")

        self.assertEqual(payload["target_count"], 1)
        self.assertEqual(payload["open_targets_release"], "26.03")
        self.assertEqual(payload["targets"][0]["ensembl_id"], "ENSG00000146648")
        self.assertEqual(payload["targets"][0]["top_disease_name"], "breast cancer")
        self.assertEqual(payload["targets"][0]["top_disease_score"], 0.93)
        self.assertEqual(payload["associations"][0]["direct_score"], 0.72)
        self.assertEqual(payload["associations"][0]["indirect_score"], 0.93)
        self.assertIn("indirect_datasource_scores", payload["associations"][0])
        self.assertIn("breast cancer", tsv)
        self.assertNotIn("measurement", tsv)
        self.assertNotIn("level of alpha-1B-glycoprotein", tsv)
        self.assertIn("eva:0.93", tsv)

    def test_read_parquet_dataset_accepts_open_targets_2603_alias_columns(self) -> None:
        pyarrow = importlib.import_module("pyarrow")
        parquet = importlib.import_module("pyarrow.parquet")

        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            dataset_dir = root / "association_by_datasource_indirect"
            dataset_dir.mkdir()
            table = pyarrow.table(
                {
                    "targetId": ["ENSG00000146648", "ENSG00000123456"],
                    "diseaseId": ["EFO_0000616", "EFO_0000616"],
                    "aggregationValue": ["europepmc", "eva"],
                    "associationScore": [0.88, 0.44],
                }
            )
            parquet.write_table(table, dataset_dir / "part-00000.parquet")

            from agdesign2.open_targets import _read_parquet_dataset

            rows = _read_parquet_dataset(
                dataset_dir,
                columns=["targetId", "diseaseId", "datasourceId", "score"],
                filter_target_ids={"ENSG00000146648"},
                column_aliases={
                    "score": ["score", "associationScore"],
                    "datasourceId": ["datasourceId", "aggregationValue"],
                },
            )

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["targetId"], "ENSG00000146648")
        self.assertEqual(rows[0]["datasourceId"], "europepmc")
        self.assertEqual(rows[0]["score"], 0.88)

    def test_resolves_latest_open_targets_release_from_local_cache_without_download(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            (root / "25.09").mkdir()
            (root / "26.03").mkdir()
            (root / "not-a-release").mkdir()

            release = resolve_open_targets_release(
                release="latest",
                data_dir=root,
                allow_network=False,
            )

        self.assertEqual(release, "26.03")

    def test_resolves_latest_open_targets_release_from_remote_listing(self) -> None:
        class FakeResponse:
            def __enter__(self):
                return self

            def __exit__(self, exc_type, exc, tb):
                return False

            def read(self):
                return b'<a href="25.09/">25.09/</a><a href="26.03/">26.03/</a><a href="README">README</a>'

        with mock.patch("agdesign2.open_targets.urlopen", return_value=FakeResponse()):
            release = resolve_open_targets_release(
                release="latest",
                ftp_base_url="https://example.org/opentargets",
                data_dir=Path("/does/not/matter"),
                allow_network=True,
            )

        self.assertEqual(release, "26.03")


class FetchAssociatedDiseasesPaginationTests(unittest.TestCase):
    def _client_returning(self, pages):
        from agdesign2.open_targets import OpenTargetsClient

        client = OpenTargetsClient.__new__(OpenTargetsClient)
        calls = {"n": 0}

        def fake_graphql(query, variables=None):
            idx = variables["page"]["index"]
            calls["n"] += 1
            associated = pages[idx] if idx < len(pages) else {"count": None, "rows": []}
            return {"target": {"associatedDiseases": associated}}

        client.graphql = fake_graphql  # type: ignore[assignment]
        return client, calls

    def test_null_count_does_not_truncate_to_one_page(self) -> None:
        # count is null on every page; pagination must continue until a short
        # page rather than stopping after page 0 (the (total_count or 0) bug).
        pages = [
            {"count": None, "rows": [{"disease": {"id": f"D{i}"}} for i in range(50)]},
            {"count": None, "rows": [{"disease": {"id": f"D{i}"}} for i in range(50, 90)]},
        ]
        client, calls = self._client_returning(pages)
        total, rows = client.fetch_associated_diseases("ENSG", page_size=50, max_pages=10)
        self.assertIsNone(total)
        self.assertEqual(len(rows), 90)
        self.assertEqual(calls["n"], 2)  # full page, then short page ends it

    def test_known_count_stops_at_total(self) -> None:
        pages = [
            {"count": 3, "rows": [{"disease": {"id": "D1"}}, {"disease": {"id": "D2"}}]},
            {"count": 3, "rows": [{"disease": {"id": "D3"}}]},
        ]
        client, calls = self._client_returning(pages)
        total, rows = client.fetch_associated_diseases("ENSG", page_size=2, max_pages=10)
        self.assertEqual(total, 3)
        self.assertEqual(len(rows), 3)


if __name__ == "__main__":
    unittest.main()
