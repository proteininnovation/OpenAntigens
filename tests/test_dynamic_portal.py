from __future__ import annotations

import hashlib
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

from agdesign2.dynamic_portal import create_app, summarize_dynamic_entries


class DynamicPortalTests(unittest.TestCase):
    def test_summarize_dynamic_entries_reads_batch_and_report_files(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            alphafold_dir = root / "alphafold"
            alphafold_dir.mkdir()
            (alphafold_dir / "P00001.pdb").write_text(
                "\n".join(
                    f"ATOM  {index:5d}  CA  MET A{index:4d}    "
                    f"{float(index):8.3f}{13.207:8.3f}{2.100:8.3f}  1.00 91.00           C"
                    for index in range(1, 121)
                )
                + "\nEND\n",
                encoding="utf-8",
            )
            report_assets_dir = root / "test_human_report_assets"
            report_assets_dir.mkdir()
            (report_assets_dir / "quality.png").write_bytes(b"png")
            (root / "open_targets_disease_associations.tsv").write_text("targetId\n", encoding="utf-8")
            (root / "open_targets_disease_associations.json").write_text(
                '{"associations": []}\n',
                encoding="utf-8",
            )
            report_json = root / "test_human_report.json"
            report_json.write_text(
                json.dumps(
                    {
                        "target": {
                            "accession": "P00001",
                            "entry_name": "TEST_HUMAN",
                            "gene_symbol": "TEST",
                            "protein_name": "Test protein",
                            "organism": "Homo sapiens",
                            "taxon_id": 9606,
                            "sequence": "M" * 120,
                        },
                        "ectodomain": {"start": 25, "end": 90, "label": "ecto", "source": "topology"},
                        "construct_details": [],
                        "construct_recommendations": [],
                        "cross_reactivity_hits": [],
                        "family_context": {"source": "HGNC", "family_names": ["Test family"], "members": []},
                        "assembly_requirements": [],
                        "experimental_constructs": [],
                        "interpro_annotations": [],
                    }
                ),
                encoding="utf-8",
            )
            summary_path = root / "batch_summary.json"
            summary_path.write_text(
                json.dumps(
                    [
                        {
                            "query": "TEST_HUMAN",
                            "status": "ok",
                            "batch_index": 1,
                            "json_report": str(report_json),
                            "markdown_report": str(root / "test_human_report.md"),
                            "error": "/private/workstation/report failed",
                        }
                    ]
                ),
                encoding="utf-8",
            )
            entries = summarize_dynamic_entries(summary_path)
            from fastapi.testclient import TestClient

            client = TestClient(create_app(summary_path))
            public_entry = client.get("/api/genes").json()[0]
            self.assertNotIn("json_report", public_entry)
            self.assertNotIn("markdown_report", public_entry)
            self.assertNotIn("error", public_entry)
            self.assertEqual(client.post("/api/refresh/TEST").status_code, 404)
            self.assertEqual(client.get("/batch_summary.json").status_code, 404)
            self.assertEqual(client.get("/test_human_report.json").status_code, 404)
            self.assertEqual(client.get("/constructs.html").status_code, 200)
            self.assertIn("TEST_HUMAN", client.get("/portal-index-data.js").text)
            self.assertEqual(client.get("/portal-disease-index.js").status_code, 200)
            for asset_name in (
                "apple-touch-icon.png",
                "favicon-32.png",
                "favicon.ico",
                "icon-192.png",
                "icon-512.png",
                "site.webmanifest",
            ):
                with self.subTest(asset_name=asset_name):
                    self.assertEqual(client.get(f"/{asset_name}").status_code, 200)
            self.assertEqual(
                client.get("/report_assets/test_human_report_assets/quality.png").content,
                b"png",
            )
            self.assertEqual(client.get("/open_targets_disease_associations.tsv").status_code, 200)
            self.assertEqual(client.get("/open_targets_disease_associations.json").status_code, 200)
            privacy_response = client.get("/privacy.html")
            self.assertEqual(privacy_response.status_code, 200)
            self.assertIn("Privacy and Analytics", privacy_response.text)
            self.assertIn("Plausible Analytics", privacy_response.text)
            detail_html = client.get("/reports/test_human.html").text
            self.assertIn("cdn.jsdelivr.net/npm/3dmol@2.5.5", detail_html)
            self.assertIn("sha256-98x4khrnLnYj6JzdERQ09Ywu/d0v/aHNISZEtAb7gBY=", detail_html)

            local_vendor = root / "portal" / "3Dmol-min.js"
            local_vendor.write_bytes(b"test 3Dmol")
            local_hash = hashlib.sha256(local_vendor.read_bytes()).hexdigest()
            with mock.patch("agdesign2.portal._THREEDMOL_SHA256", local_hash):
                local_client = TestClient(create_app(summary_path))
                self.assertEqual(local_client.get("/3Dmol-min.js").content, local_vendor.read_bytes())
                local_detail_html = local_client.get("/reports/test_human.html").text
            self.assertIn('<script src="../3Dmol-min.js"></script>', local_detail_html)
            self.assertNotIn("cdn.jsdelivr.net/npm/3dmol@2.5.5", local_detail_html)
            self.assertEqual(len(entries), 1)
            self.assertEqual(entries[0]["entry_name"], "TEST_HUMAN")
            self.assertEqual(entries[0]["gene_symbol"], "TEST")
            self.assertTrue(entries[0]["has_family_context"])

    def test_invalid_report_json_is_not_silently_ignored(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            report_path = root / "broken_report.json"
            report_path.write_text("{", encoding="utf-8")
            summary_path = root / "batch_summary.json"
            summary_path.write_text(
                json.dumps([{"query": "BROKEN_HUMAN", "json_report": str(report_path)}]),
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "Could not read report JSON"):
                summarize_dynamic_entries(summary_path)


if __name__ == "__main__":
    unittest.main()
