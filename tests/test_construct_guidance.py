"""Construct selection help stays concise and reachable in every serving mode."""
import json
import re
import tempfile
import unittest
from html import unescape
from pathlib import Path

from fastapi.testclient import TestClient

from agdesign2.dynamic_portal import create_app
from agdesign2.portal import build_portal, render_constructs_page, render_methods_page
from agdesign2.portal_db import build_portal_database, build_portal_from_database, create_db_portal_app
from tests.test_portal_db import _write_fixture


class ConstructGuidanceTests(unittest.TestCase):
    def test_guide_is_first_concise_and_links_to_existing_details(self):
        html = render_constructs_page()
        guide = re.search(r'<section[^>]*id="choose-constructs".*?</section>', html, re.S).group()
        self.assertEqual(html.count('id="choose-constructs"'), 1)
        self.assertLess(html.index(guide), html.index("What the construct list contains"))
        self.assertLessEqual(len(unescape(re.sub(r"<[^>]+>", " ", guide)).split()), 285)
        self.assertEqual(guide.count("<li>"), 7)  # Four steps and three final checks.
        for anchor in re.findall(r'href="#([^" ]+)"', guide):
            self.assertEqual(html.count(f'id="{anchor}"'), 1)
        classes = re.search(r'<section[^>]*id="construct-classes".*?</section>', html, re.S).group()
        self.assertEqual(classes.count("<tr>"), 8)  # Header and seven categories.
        self.assertIn("Full design region when neither label applies", classes)
        self.assertEqual(html.count('id="after-export"'), 1)
        self.assertLess(html.index('id="after-export"'), html.index("What the construct list contains"))
        for portal_title in ("OpenAntigens", "OpenAntigens Mouse"):
            page = render_constructs_page(portal_title=portal_title)
            after_export = re.search(r'<section[^>]*id="after-export".*?</section>', page, re.S).group()
            self.assertEqual(
                re.findall(r'href="https://doi.org/([^"]+)', after_export),
                [
                    "10.1038/nmeth.f.202", "10.1016/bs.mie.2021.08.019",
                    "10.3389/fbioe.2025.1661193", "10.1007/978-1-0716-3878-1_6",
                    "10.1016/B978-0-12-420070-8.00013-1",
                    "10.1002/0471140864.ps0909s73", "10.1002/0471140864.ps0601s80",
                ],
            )
            self.assertIn("https://blog.addgene.org/plasmids-101-protein-tags", after_export)
            self.assertNotIn("Tegel", after_export)
            self.assertNotIn("https://info.addgene.org/plasmids-101-topic-page", after_export)
            self.assertIn("For an exposed unpaired cysteine", page)
            self.assertIn("Furin-like motifs warrant particular attention", page)
        self.assertIn("without partner chains", render_methods_page())
        self.assertIn("Other free-text interaction rows are not shown", render_methods_page())

    def test_help_link_and_guide_across_serving_modes_and_mouse_theme(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            summary, payload = _write_fixture(root, with_structure=True)
            payload["target"]["canonical_isoform_id"] = "P00001-2"
            payload["construct_details"][0]["name"] = "pdb_1abc_a"
            (root / "test_human_report.json").write_text(json.dumps(payload))
            pages = []
            db = root / "portal.sqlite"
            build_portal_database(summary, db)
            for name, app in (("live-json", create_app(summary)), ("sqlite", create_db_portal_app(db))):
                with self.subTest(mode=name), TestClient(app) as client:
                    report = client.get("/reports/test_human.html")
                    self.assertEqual(report.status_code, 200)
                    self.assertIn('href="../constructs.html#choose-constructs" target="_blank" rel="noopener"', report.text)
                    pages.append(report.text)
                    guide = client.get("/constructs.html")
                    self.assertEqual(guide.status_code, 200)
                    self.assertEqual(guide.text, render_constructs_page())
            guide_section = re.search(r'<section[^>]*id="choose-constructs".*?</section>', render_constructs_page(), re.S).group()
            for theme in ("human", "mouse"):
                site = root / theme
                build_portal(summary, output_dir=site, refresh_stale_reports=False,
                             fetch_literature=False, bundle_vendor_assets=False, theme=theme)
                self.assertIn(guide_section, (site / "constructs.html").read_text())
                self.assertIn('href="../constructs.html#choose-constructs"', (site / "reports/test_human.html").read_text())
                pages.append((site / "reports/test_human.html").read_text())
            site = root / "db-static"
            build_portal_from_database(db, site)
            self.assertIn(guide_section, (site / "constructs.html").read_text())
            self.assertIn('href="../constructs.html#choose-constructs"', (site / "reports/test_human.html").read_text())
            pages.append((site / "reports/test_human.html").read_text())
            for page in pages:
                self.assertIn('href="../constructs.html#after-export" target="_blank" rel="noopener"', page)
                self.assertIn("UniProt canonical (P00001-2)", page)
                self.assertIn("Alternative isoforms are not analyzed separately", page)
                self.assertEqual(page.count("Exposure is calculated from the target-chain model"), 1)
                self.assertIn('href="#interactions-assembly"', page)
                self.assertIn('id="interactions-assembly"', page)
                self.assertEqual(page.count("These sequences follow mapped PDB boundaries"), 1)
