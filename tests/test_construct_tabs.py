"""Construct category coverage and parity across report serving paths."""
import json
import re
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from fastapi.testclient import TestClient

from agdesign2.dynamic_portal import create_app
from agdesign2.portal import _render_construct_tabs, build_portal
from agdesign2.portal_db import build_portal_database, build_portal_from_database, create_db_portal_app
from tests.test_portal_db import _write_fixture
from tests.test_report_navigation import section_targets, with_optional_sections


def tab_constructs(base):
    # Identical boundaries must retain every card, including both PDB records.
    names = [
        "full_ectodomain", "pdb_1abc_a", "pdb_2abc_b", "domain_example", "repeat_example",
        "curated_uniprot_region", "strict_domain_1", "structural_region_2", "lenient_region_1",
        "membrane_expression_full_length", "membrane_expression_n_tail_trimmed",
        "membrane_expression_c_tail_trimmed", "membrane_expression_terminal_trimmed",
        "membrane_expression_gpcr_c_tail_keep_15",
    ]
    return [
        {**base, "name": name, "classification": "membrane_expression" if name == "pdb_2abc_b" or name.startswith("membrane_expression_") else "complete_domain"}
        for name in names
    ]


class ConstructTabTests(unittest.TestCase):
    def render_tabs(self, constructs, label="Extracellular Region"):
        with patch("agdesign2.portal._render_construct_card", side_effect=lambda c, **kwargs: f'<article>{c["name"]}</article>'):
            return _render_construct_tabs(
                constructs, design_region_label=label, target={}, batch_dir=Path("."),
                page_dir=Path("."), furin_sites=[], ptms=[], allow_structure_assets=False,
                identity_label="Identity to human",
            )

    def test_categories_preserve_every_card_once_and_use_generation_kind(self):
        constructs = tab_constructs({"start": 1, "end": 10})
        html = self.render_tabs(constructs)
        expected = [
            "Full ectodomain (1)", "PDB (2)", "Annotated domains / repeats (3)",
            "Strict (2)", "Lenient (1)", "Full-length multipass (1)", "Trimmed multipass (4)",
        ]
        self.assertEqual(re.findall(r'role="tab"[^>]*>(.*?)</button>', html), expected)
        self.assertEqual(sum(int(n) for n in re.findall(r'\((\d+)\)</button>', html)), len(constructs))
        for construct in constructs:
            self.assertEqual(html.count(f'<article>{construct["name"]}</article>'), 1)
        panels = re.findall(r'<div role="tabpanel"[^>]*>(.*?)</div></div>', html, re.S)
        self.assertEqual(len(panels), 7)
        self.assertIn("Soluble PDB constructs", panels[1])
        self.assertIn("Membrane PDB constructs", panels[1])
        self.assertEqual(html.count("These sequences follow mapped PDB boundaries"), 1)
        self.assertLess(panels[1].index("These sequences"), panels[1].index("Soluble PDB constructs"))
        self.assertIn("pdb_2abc_b", panels[1])
        self.assertNotIn("pdb_2abc_b", panels[6])
        self.assertIn("structural_region_2", panels[3])
        self.assertEqual(html.count('aria-selected="true"'), 1)
        self.assertEqual(html.count('tabindex="0" hidden'), 6)

    def test_empty_single_category_and_design_region_labels(self):
        self.assertEqual(self.render_tabs([]), "<p>No constructs available.</p>")
        for label, expected in [
            ("Secreted Region", "Full secreted region"),
            ("GPI-Anchored Extracellular Region", "Full ectodomain"),
            ("Design Region", "Full design region"),
        ]:
            with self.subTest(label=label):
                html = self.render_tabs([{"name": "full_ectodomain"}], label)
                self.assertIn(f">{expected} (1)</button>", html)
                self.assertEqual(html.count('<button type="button" role="tab"'), 1)
                self.assertNotIn('tabindex="0" hidden', html)
        html = self.render_tabs([{"name": "lenient_region_1"}])
        self.assertNotIn("These sequences follow mapped PDB boundaries", html)
        self.assertIn('id="construct-tab-4" aria-controls="construct-panel-4" aria-selected="true"', html)

    def test_static_json_live_json_and_sqlite_render_the_same_tabs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            summary, report = _write_fixture(root)
            report = with_optional_sections(report)
            report["construct_details"] = tab_constructs({"start": 1, "end": 10, "sequence": "M" * 10})
            (root / "test_human_report.json").write_text(json.dumps(report))
            build_portal(summary, output_dir=root / "static", refresh_stale_reports=False, fetch_literature=False)
            db = root / "portal.sqlite"
            build_portal_database(summary, db)
            build_portal_from_database(db, root / "db-static")
            pages = [(root / directory / "reports/test_human.html").read_text() for directory in ("static", "db-static")]
            for app in (create_app(summary), create_db_portal_app(db)):
                with TestClient(app) as client:
                    response = client.get("/reports/test_human.html")
                    self.assertEqual(response.status_code, 200)
                    pages.append(response.text)
            expected_tabs = re.findall(r'<button[^>]*role="tab".*?</button>', pages[0])
            self.assertEqual(len(expected_tabs), 7)
            for page in pages:
                summary_html = page.split('id="construct-summary"', 1)[1].split('</tbody>', 1)[0]
                self.assertEqual(re.findall(r'<tr><td><code>(.*?)</code>', summary_html), [
                    "full_ectodomain", "pdb_1abc_a", "strict_domain_1", "structural_region_2",
                    "lenient_region_1", "curated_uniprot_region", "domain_example", "repeat_example",
                ])
                self.assertEqual(section_targets(page), section_targets(pages[0]))
                self.assertEqual(len(section_targets(page)), 18)
                self.assertEqual(re.findall(r'<button[^>]*role="tab".*?</button>', page), expected_tabs)
                details = page.split('id="construct-details"', 1)[1].split('</section>', 1)[0]
                self.assertEqual(details.count('class="construct-card"'), len(report["construct_details"]))


if __name__ == "__main__":
    unittest.main()
