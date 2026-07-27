from __future__ import annotations

import hashlib
import json
import sqlite3
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from agdesign2.dynamic_portal import summarize_dynamic_entries
from agdesign2.portal import portal_index_json, portal_index_tsv
from agdesign2.portal_db import (
    PortalDatabaseError,
    build_portal_database,
    build_portal_from_database,
    create_db_portal_app,
    load_index_entries_from_db,
    load_report_from_db,
)
from scripts.build_public_release import validate_portal


class PortalDatabaseTests(unittest.TestCase):
    def test_public_indexes_do_not_expose_raw_error_details(self) -> None:
        entry = {
            "query": "BROKEN_HUMAN",
            "status": "error",
            "error": "/private/workstation/report failed",
        }

        for output in (portal_index_json([entry]), portal_index_tsv([entry])):
            self.assertIn("analysis failed", output)
            self.assertNotIn("/private/workstation", output)

    def test_imports_batch_reports_and_matches_dynamic_entries(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            summary_path, report_payload = _write_fixture(root)
            db_path = root / "data" / "openantigen.sqlite"

            build_portal_database(summary_path, db_path)

            db_entries = load_index_entries_from_db(db_path)
            dynamic_entries = summarize_dynamic_entries(summary_path)
            self.assertEqual(db_entries, _without_json_paths(dynamic_entries))
            self.assertEqual(load_report_from_db(db_path, "TEST_HUMAN"), report_payload)
            self.assertEqual(load_report_from_db(db_path, "test_human.html"), report_payload)
            self.assertEqual(portal_index_tsv(db_entries), portal_index_tsv(_without_json_paths(dynamic_entries)))
            self.assertEqual(portal_index_json(db_entries), portal_index_json(_without_json_paths(dynamic_entries)))
            with sqlite3.connect(db_path) as conn:
                metadata = dict(conn.execute("SELECT key, value FROM metadata"))
            self.assertNotIn("source_summary_path", metadata)
            self.assertEqual(metadata["assets_dir"], "assets")
            self.assertNotIn(tmpdir, json.dumps(metadata))

    def test_import_fails_loudly_for_missing_report_json(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            summary_path = root / "batch_summary.json"
            summary_path.write_text(
                json.dumps([{"query": "TEST_HUMAN", "status": "ok", "json_report": str(root / "missing.json")}]),
                encoding="utf-8",
            )

            with self.assertRaisesRegex(PortalDatabaseError, "missing report JSON"):
                build_portal_database(summary_path, root / "openantigen.sqlite")

    def test_import_fails_loudly_for_duplicate_entry_names(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            summary_path, _ = _write_fixture(root)
            report_2 = root / "test_human_report_copy.json"
            report_2.write_text((root / "test_human_report.json").read_text(encoding="utf-8"), encoding="utf-8")
            rows = json.loads(summary_path.read_text(encoding="utf-8"))
            rows.append({**rows[0], "json_report": str(report_2)})
            summary_path.write_text(json.dumps(rows), encoding="utf-8")

            with self.assertRaisesRegex(PortalDatabaseError, "Duplicate entry_name"):
                build_portal_database(summary_path, root / "openantigen.sqlite")

    def test_import_fails_loudly_for_missing_required_target_fields(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            summary_path, report_payload = _write_fixture(root)
            report_payload["target"]["accession"] = ""
            (root / "test_human_report.json").write_text(json.dumps(report_payload), encoding="utf-8")

            with self.assertRaisesRegex(PortalDatabaseError, "target.accession"):
                build_portal_database(summary_path, root / "openantigen.sqlite")

    def test_records_asset_metadata_for_copied_report_assets(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            summary_path, _ = _write_fixture(root)
            asset_dir = root / "test_human_report_assets"
            asset_dir.mkdir()
            (asset_dir / "quality.png").write_bytes(b"not really a png")
            db_path = root / "openantigen.sqlite"

            build_portal_database(summary_path, db_path)

            with sqlite3.connect(db_path) as conn:
                rows = conn.execute("SELECT path, kind, size_bytes, sha256 FROM assets WHERE kind = 'report_asset'").fetchall()
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0][0], "report_assets/test_human_report_assets/quality.png")
            self.assertEqual(rows[0][1], "report_asset")
            self.assertGreater(rows[0][2], 0)
            self.assertEqual(len(rows[0][3]), 64)

    def test_db_backed_fastapi_serves_pages_api_and_no_refresh_endpoint(self) -> None:
        try:
            from fastapi.testclient import TestClient
        except Exception as exc:  # pragma: no cover - dependency guard
            self.skipTest(f"FastAPI TestClient unavailable: {exc}")
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            summary_path, _ = _write_fixture(root)
            db_path = root / "openantigen.sqlite"
            build_portal_database(summary_path, db_path)
            client = TestClient(create_db_portal_app(db_path))

            self.assertEqual(client.get("/index.html").status_code, 200)
            self.assertIn("portal-index-data.js", client.get("/index.html").text)
            self.assertIn("TEST_HUMAN", client.get("/portal-index-data.js").text)
            self.assertEqual(client.get("/reports/test_human.html").status_code, 200)
            self.assertEqual(client.get("/api/genes").json()[0]["entry_name"], "TEST_HUMAN")
            self.assertEqual(client.get("/api/gene/TEST_HUMAN").json()["target"]["entry_name"], "TEST_HUMAN")
            self.assertIn("entry_name", client.get("/downloads/agdesign2_portal_index.tsv").text)
            self.assertEqual(client.get("/portal_metadata.json").json()["target_count"], 1)
            self.assertEqual(client.post("/api/refresh/TEST_HUMAN").status_code, 405)
            privacy_response = client.get("/privacy.html")
            self.assertEqual(privacy_response.status_code, 200)
            self.assertIn("Privacy and Analytics", privacy_response.text)
            self.assertIn("Plausible Analytics", privacy_response.text)

    def test_static_export_from_db_writes_portal_files(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            summary_path, _ = _write_fixture(root)
            db_path = root / "openantigen.sqlite"
            build_portal_database(summary_path, db_path)

            index_path = build_portal_from_database(db_path, root / "public_site")

            self.assertTrue(index_path.exists())
            self.assertTrue((root / "public_site" / "reports" / "test_human.html").exists())
            self.assertTrue((root / "public_site" / "downloads" / "agdesign2_portal_index.tsv").exists())
            self.assertTrue((root / "public_site" / "portal_metadata.json").exists())
            self.assertTrue((root / "public_site" / "portal.css").exists())
            self.assertTrue((root / "public_site" / ".htaccess").exists())
            privacy_html = (root / "public_site" / "privacy.html").read_text(encoding="utf-8")
            self.assertIn("Privacy and Analytics", privacy_html)
            self.assertIn("Plausible Analytics", privacy_html)
            validate_portal(root / "public_site")

    def test_db_dynamic_and_static_outputs_share_payloads_assets_and_viewer_scripts(self) -> None:
        try:
            from fastapi.testclient import TestClient
        except Exception as exc:  # pragma: no cover - dependency guard
            self.skipTest(f"FastAPI TestClient unavailable: {exc}")
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            summary_path, _ = _write_fixture(root, with_structure=True, with_vendor=True)
            db_path = root / "openantigen.sqlite"
            vendor_path = root / "vendor" / "3Dmol-min.js"
            vendor_hash = hashlib.sha256(vendor_path.read_bytes()).hexdigest()
            with mock.patch("agdesign2.portal._THREEDMOL_SHA256", vendor_hash):
                build_portal_database(summary_path, db_path, bundle_vendor_assets=True)
            static_root = root / "public_site"
            with mock.patch("agdesign2.portal._THREEDMOL_SHA256", vendor_hash):
                build_portal_from_database(db_path, static_root)
            client = TestClient(create_db_portal_app(db_path))

            checks = {
                "/portal-index-data.js": static_root / "portal-index-data.js",
                "/downloads/agdesign2_portal_index.tsv": static_root / "downloads" / "agdesign2_portal_index.tsv",
                "/downloads/agdesign2_portal_index.json": static_root / "downloads" / "agdesign2_portal_index.json",
                "/downloads/download_manifest.json": static_root / "downloads" / "download_manifest.json",
                "/3Dmol-min.js": static_root / "3Dmol-min.js",
                "/structures/test_human.pdb": static_root / "structures" / "test_human.pdb",
            }
            for route, static_path in checks.items():
                response = client.get(route)
                self.assertEqual(response.status_code, 200, route)
                self.assertTrue(static_path.exists(), str(static_path))
                self.assertEqual(response.content, static_path.read_bytes(), route)

            with mock.patch("agdesign2.portal._THREEDMOL_SHA256", vendor_hash):
                dynamic_report = client.get("/reports/test_human.html")
            self.assertEqual(dynamic_report.status_code, 200)
            static_report = (static_root / "reports" / "test_human.html").read_text(encoding="utf-8")
            for expected in ("TEST_HUMAN", "structureViewer", "plddtPlot", "paePlot", "../3Dmol-min.js"):
                self.assertIn(expected, dynamic_report.text)
                self.assertIn(expected, static_report)
            self.assertIn("../report-viewer.js", static_report)
            self.assertIn("../report_scripts/test_human.js", static_report)

            report_viewer_js = (static_root / "report-viewer.js").read_text(encoding="utf-8")
            report_data_js = (static_root / "report_scripts" / "test_human.js").read_text(encoding="utf-8")
            self.assertIn("Interactive AlphaFold viewer ready", report_viewer_js)
            self.assertIn("const viewerData", report_data_js)
            self.assertEqual(report_viewer_js.count("const selectionFasta"), 1)
            self.assertEqual(report_data_js.count("const selectionFasta"), 0)


def _write_fixture(root: Path, *, with_structure: bool = False, with_vendor: bool = False) -> tuple[Path, dict[str, object]]:
    sequence = "M" * 10 if with_structure else "M" * 120
    report_payload: dict[str, object] = {
        "target": {
            "accession": "P00001",
            "entry_name": "TEST_HUMAN",
            "gene_symbol": "TEST",
            "protein_name": "Test protein",
            "organism": "Homo sapiens",
            "taxon_id": 9606,
            "sequence": sequence,
        },
        "ectodomain": {
            "start": 1 if with_structure else 25,
            "end": len(sequence) if with_structure else 90,
            "label": "ecto",
            "source": "topology",
        },
        "construct_details": [
            {
                "name": "TEST_HUMAN_1-10",
                "start": 1,
                "end": 10,
                "sequence": sequence,
                "kind": "full_design_region",
            }
        ]
        if with_structure
        else [],
        "construct_recommendations": [],
        "cross_reactivity_hits": [],
        "family_context": {"source": "HGNC", "family_names": ["Test family"], "members": []},
        "assembly_requirements": [],
        "experimental_constructs": [],
        "interpro_annotations": [],
    }
    report_json = root / "test_human_report.json"
    report_json.write_text(json.dumps(report_payload), encoding="utf-8")
    summary_path = root / "batch_summary.json"
    summary_path.write_text(
        json.dumps(
            [
                {
                    "query": "TEST_HUMAN",
                    "status": "ok",
                    "batch_index": 1,
                    "pubtator_count": 1,
                    "pubtator_query_url": "https://www.ncbi.nlm.nih.gov/research/pubtator3/docsum?text=%40GENE_TEST",
                    "json_report": str(report_json),
                    "markdown_report": str(root / "test_human_report.md"),
                }
            ]
        ),
        encoding="utf-8",
    )
    if with_structure:
        alphafold_dir = root / "alphafold"
        alphafold_dir.mkdir()
        (alphafold_dir / "P00001.pdb").write_text(_minimal_methionine_pdb(len(sequence)), encoding="utf-8")
    if with_vendor:
        vendor_dir = root / "vendor"
        vendor_dir.mkdir()
        (vendor_dir / "3Dmol-min.js").write_text("window.$3Dmol = window.$3Dmol || {};\n", encoding="utf-8")
    return summary_path, report_payload


def _minimal_methionine_pdb(length: int) -> str:
    lines = [
        f"ATOM  {index:5d}  CA  MET A{index:4d}    {float(index):8.3f}{13.207:8.3f}{2.100:8.3f}  1.00 91.00           C"
        for index in range(1, length + 1)
    ]
    return "\n".join(lines) + "\nEND\n"


def _without_json_paths(entries: list[dict[str, object]]) -> list[dict[str, object]]:
    normalized: list[dict[str, object]] = []
    for entry in entries:
        item = dict(entry)
        if item.get("json_report"):
            item["json_report"] = ""
        if item.get("markdown_report"):
            item["markdown_report"] = ""
        normalized.append(item)
    return normalized


if __name__ == "__main__":
    unittest.main()
