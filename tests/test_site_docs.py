"""Citation provenance and agent-document parity across delivery paths."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from fastapi.testclient import TestClient

from agdesign2.dynamic_portal import create_app
from agdesign2.portal import build_portal, render_agent_guide_page
from agdesign2.portal_db import build_portal_database, build_portal_from_database, create_db_portal_app
from agdesign2.site_docs import PAPER, agent_guide_markdown, citation_cff, citation_downloads, citation_text, site_document_files
from scripts.build_public_release import copy_public_portal, render_release_readme
from tests.test_portal_db import _write_fixture

ROOT = Path(__file__).resolve().parents[1]
DOI = "10.64898/2026.07.30.741735"


class SiteDocumentationTests(unittest.TestCase):
    def test_verified_reference_and_repository_copies(self):
        self.assertEqual(PAPER["doi"], DOI)
        self.assertEqual(PAPER["date"], "2026-08-04")
        self.assertEqual([a["family"] for a in PAPER["authors"]], ["Teixeira", "Zhu", "Kothiwal", "Cao", "Mills"])
        self.assertEqual((ROOT / "llms.txt").read_text(), agent_guide_markdown())
        self.assertNotIn("{{citation}}", agent_guide_markdown())
        self.assertEqual((ROOT / "CITATION.cff").read_text(), citation_cff())
        self.assertIn("preferred-citation:", citation_cff())
        self.assertIn("  doi: " + DOI, citation_cff())
        for name in ("README.md", "DATA_LICENSE.md"):
            self.assertIn(citation_text(), (ROOT / name).read_text())
        for content in citation_downloads().values():
            self.assertIn(DOI, content)
            self.assertIn(PAPER["title"], content)
            self.assertIn("Preprint", content)

    def test_static_and_live_routes_include_the_same_documents(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            summary, _ = _write_fixture(root, with_structure=False)
            site = root / "site"
            with patch("agdesign2.portal._get_pubtator_literature_info", return_value={"count": 1, "query_url": "https://www.ncbi.nlm.nih.gov/research/pubtator3/docsum?text=TEST"}):
                build_portal(summary, output_dir=site, refresh_stale_reports=False)
            db = root / "portal.sqlite"
            build_portal_database(summary, db)
            db_site = root / "db-site"
            build_portal_from_database(db, db_site)
            expected = {**site_document_files(), "agent-guide.html": render_agent_guide_page()}
            for app in (create_app(summary), create_db_portal_app(db)):
                with TestClient(app) as client:
                    for name, content in expected.items():
                        response = client.get("/" + name)
                        self.assertEqual(response.status_code, 200, name)
                        self.assertEqual(response.text, content, name)
                    self.assertEqual(client.get("/downloads/private.txt").status_code, 404)
                    self.assertIn("application/x-research-info-systems", client.get("/downloads/openantigens.ris").headers["content-type"])
            for folder in (site, db_site):
                for name, content in site_document_files().items():
                    self.assertEqual((folder / name).read_text(), content, name)
                self.assertIn(DOI, (folder / "agent-guide.html").read_text())
                for name in ("index.html", "help.html", "terms.html", "methods.html", "downloads.html", "reports/test_human.html"):
                    html = (folder / name).read_text()
                    self.assertIn(DOI, html, name)
                    prefix = "../" if name.startswith("reports/") else ""
                    self.assertIn(f'href="{prefix}help.html#cite-openantigens"', html)
                metadata = json.loads((folder / "portal_metadata.json").read_text())
                manifest = json.loads((folder / "downloads/download_manifest.json").read_text())
                self.assertEqual(metadata["citation"]["doi"], DOI)
                self.assertEqual(manifest["citation"]["doi"], DOI)
                self.assertTrue({"downloads/openantigens.bib", "downloads/openantigens.ris"} <= {f["path"] for f in manifest["files"]})
            release = root / "release"
            copy_public_portal(site, release)
            self.assertTrue((release / "downloads/openantigens.bib").is_file())
            readme = render_release_readme({"portal_metadata": metadata, "release_created_utc": "2026-09-17", "software_version": "test", "file_count": 1, "total_size_human": "1 KB", "html_report_count": 1, "download_file_count": 4, "structure_file_count": 0})
            self.assertIn(citation_text(), readme)

    def test_mouse_guide_links_keep_human_and_mouse_portals_distinct(self):
        html = render_agent_guide_page(portal_title="OpenAntigens Mouse")
        self.assertIn('href="../index.html"', html)
        self.assertIn('href="../mouse/index.html"', html)
        self.assertIn('href="../help.html#cite-openantigens"', html)
        self.assertIn('href="help.html#cite-openantigens"', html)  # Local footer.
        self.assertIn("Start at https://openantigens.org/mouse/index.html", html)


if __name__ == "__main__":
    unittest.main()
