"""The release checker must catch each defect it names, and pass a clean page.

Run against the frozen 2026-09-30 build it fails with the plan's own counts
(blast-negative-beside-hits 1,658 human / 1,598 mouse, mouse-self-mapping-row
4,359, duplicated-note 4,060 ...), which is what shows it detects them.
"""

from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("release_checks", ROOT / "scripts" / "release_checks.py")
release_checks = importlib.util.module_from_spec(_spec)
sys.modules["release_checks"] = release_checks
_spec.loader.exec_module(release_checks)


def _site(pages: dict[str, str]) -> TemporaryDirectory:
    tmp = TemporaryDirectory()
    root = Path(tmp.name)
    for rel, html in pages.items():
        path = root / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(html, encoding="utf-8")
    return tmp


def _results(pages: dict[str, str]) -> dict[str, "release_checks.Check"]:
    with _site(pages) as tmp:
        return {c.key: c for c in release_checks.run(Path(tmp), None, None)}


CLEAN = "<section id=\"notes\"><ul><li>a</li><li>b</li></ul></section>"
INDEX = "gene\tidentity\nX\t90\n"
CLEAN_SITE = {
    "reports/x_human.html": CLEAN,
    "mouse/reports/x_mouse.html": CLEAN,
    "downloads/agdesign2_portal_index.tsv": INDEX,
    "mouse/downloads/agdesign2_portal_index.tsv": INDEX,
}


class ReleaseCheckTests(unittest.TestCase):
    def test_a_clean_page_passes(self) -> None:
        results = _results(CLEAN_SITE)
        failed = [k for k, c in results.items() if not (c.passed or c.skipped)]
        self.assertEqual(failed, [])

    def test_an_empty_site_fails(self) -> None:
        results = _results({})
        self.assertEqual(dict(results["portal-present"].failures), {"human": 1, "mouse": 1})
        self.assertFalse(results["index-column-blank"].passed)

    def test_a_site_missing_its_mouse_portal_fails(self) -> None:
        site = {k: v for k, v in CLEAN_SITE.items() if not k.startswith("mouse/")}
        self.assertEqual(dict(_results(site)["portal-present"].failures), {"mouse": 1})

    def test_an_empty_report_directory_fails(self) -> None:
        # Both given, but one holds no *_report.json: the JSON checks would
        # otherwise pass over nothing at all.
        with _site(CLEAN_SITE) as tmp, TemporaryDirectory() as empty:
            checks = {c.key: c for c in release_checks.run(Path(tmp), Path(empty), Path(empty))}
        self.assertEqual(dict(checks["portal-present"].failures), {"human": 1, "mouse": 1})

    def test_one_json_directory_alone_is_refused(self) -> None:
        # Running one portal's JSON checks reported the other's as PASS.
        with _site(CLEAN_SITE) as tmp:
            with self.assertRaises(ValueError):
                release_checks.run(Path(tmp), Path(tmp), None)

    def test_each_html_defect_is_caught(self) -> None:
        cases = {
            "blast-negative-beside-hits": (
                "reports/a_human.html",
                "The completed BLAST search returned no sequence similarity hits. Showing 50 full-length hit(s)",
            ),
            "interpro-negative-beside-family": (
                "reports/b_human.html", "No InterPro annotations. Canonical family: X (IPR005461)",
            ),
            "mouse-self-mapping-row": ("mouse/reports/c_mouse.html", "<code>Il7_MOUSE_mapping</code>"),
            "duplicated-note": (
                "mouse/reports/d_mouse.html",
                "<section id=\"notes\"><ul><li>same</li><li>same</li></ul></section>",
            ),
            "motif-detected-label": ("reports/e_human.html", "<td>detected</td>"),
            "homolog-failure-as-no-records": (
                "mouse/reports/f_mouse.html",
                "Species homolog lookup failed <section id=\"homology\">No homolog records.</section>",
            ),
            "duplicate-construct-name": (
                "reports/g_human.html",
                "<section id=\"construct-summary\"><tr><td>A</td><td>x</td><td>1-10</td></tr>"
                "<tr><td>A</td><td>x</td><td>1-20</td></tr></section>",
            ),
            "oversized-construct-unwarned": (
                "reports/h_human.html",
                "<section id=\"construct-summary\"><tr><td>A</td><td>x</td><td>1-1600</td></tr></section>",
            ),
            "negative-site-as-ptm": (
                "mouse/reports/i_mouse.html", "<td>SITE</td><td>57</td><td>Not glycosylated</td>",
            ),
        }
        self.assertEqual(set(cases), {c.key for c in release_checks._html_checks()})
        for key, (path, html) in cases.items():
            with self.subTest(check=key):
                self.assertFalse(_results({path: html})[key].passed)

    def test_json_checks_are_skipped_not_passed_without_json(self) -> None:
        results = _results(CLEAN_SITE)
        self.assertTrue(results["gpi-region-past-omega"].skipped)
        self.assertFalse(results["gpi-region-past-omega"].passed)


if __name__ == "__main__":
    unittest.main()
