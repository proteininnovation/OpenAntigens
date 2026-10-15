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


class ReleaseCheckTests(unittest.TestCase):
    def test_a_clean_page_passes(self) -> None:
        results = _results({"reports/x_human.html": CLEAN})
        failed = [k for k, c in results.items() if not (c.passed or c.skipped)]
        self.assertEqual(failed, [])

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
        }
        for key, (path, html) in cases.items():
            with self.subTest(check=key):
                self.assertFalse(_results({path: html})[key].passed)

    def test_json_checks_are_skipped_not_passed_without_json(self) -> None:
        results = _results({"reports/x_human.html": CLEAN})
        self.assertTrue(results["gpi-region-past-omega"].skipped)
        self.assertFalse(results["gpi-region-past-omega"].passed)


if __name__ == "__main__":
    unittest.main()
