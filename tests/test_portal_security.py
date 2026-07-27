"""XSS-regression tests for portal rendering.

The detail page's structure-viewer embeds report-derived strings (gene symbol,
sequence, accession, viewer payload) inside an executable inline ``<script>``.
``json.dumps`` does not escape ``<``/``/``, so a ``</script>`` substring in any
of those fields closes the script element and the rest is parsed as HTML. The
static/db builds externalize this runtime to a ``.js`` file (neutralizing it),
but the dynamic portal (``serve-portal``) serves ``render_detail_page`` output
inline, where the breakout is live.

The fix routes those interpolations through ``_js_json``, which escapes the
HTML-significant characters as JSON ``\\u`` escapes — losslessly, so the decoded
JavaScript value (and therefore the science/data shown) is unchanged. The
client-side construct-card ``redraw()`` similarly must HTML-escape values it
writes via ``innerHTML``.
"""

from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from agdesign2.portal import _build_index_entry, _construct_name_html, _js_json, render_detail_page  # noqa: E402

FIXTURES = ROOT / "tests" / "golden" / "fixtures"
BREAKOUT = "</script><img src=x onerror=alert(1)>"


class JsJsonTests(unittest.TestCase):
    def test_escapes_script_breakout(self) -> None:
        out = _js_json(BREAKOUT)
        self.assertNotIn("</script>", out)
        self.assertNotIn("<img", out)
        self.assertIn("\\u003c", out)  # '<' escaped

    def test_is_lossless(self) -> None:
        # Escaping must not change the decoded value — the construct sequences,
        # identities, and names round-trip exactly.
        for value in [
            BREAKOUT,
            "MKWVTFISLLFLFSSAYS",
            {"geneSymbol": "CD19<>&", "rows": [{"seq": "AC</script>DE"}]},
            123,
            True,
            None,
            ["a<b", "c>d", "e&f"],
        ]:
            self.assertEqual(json.loads(_js_json(value)), value)

    def test_escapes_comment_and_ampersand(self) -> None:
        out = _js_json("<!-- & -->")
        self.assertNotIn("<!--", out)
        self.assertNotIn("&", out)


class PortalRenderXssTests(unittest.TestCase):
    def _render(self, report: dict, batch_dir: Path) -> str:
        item = {
            "json_report": f"{report['target']['entry_name'].lower()}_report.json",
            "status": "ok",
            "query": report["target"].get("gene_symbol"),
        }
        entry = _build_index_entry(item, report, None)
        return render_detail_page(entry, report, batch_dir=batch_dir, sibling_link=None)

    def test_no_raw_breakout_in_rendered_page(self) -> None:
        report = json.loads((FIXTURES / "sfta3_human_report.json").read_text(encoding="utf-8"))
        report["target"]["gene_symbol"] = BREAKOUT
        report["target"]["protein_name"] = BREAKOUT
        with TemporaryDirectory() as root:
            html = self._render(report, Path(root))
        # No executable-context breakout anywhere on the page.
        self.assertNotIn("</script><img", html)

    def test_construct_card_runtime_escapes_html(self) -> None:
        # redraw() builds table rows via innerHTML; it must route interpolated
        # cell values through an HTML escaper rather than inject them raw.
        report = json.loads((FIXTURES / "spr2a_human_report.json").read_text(encoding="utf-8"))
        with TemporaryDirectory() as root:
            html = self._render(report, Path(root))
        self.assertIn("renderConstructMutationCard", html)
        self.assertIn("escapeHtml", html)

    def test_interactions_assembly_shows_complex_portal_records_and_lookup_states(self) -> None:
        report = json.loads((FIXTURES / "sfta3_human_report.json").read_text(encoding="utf-8"))
        report["complex_portal_lookup"] = {"human": {"status": "no_hits"}}
        with TemporaryDirectory() as root:
            no_hits_html = self._render(report, Path(root))
        self.assertIn("Interactions / Assembly", no_hits_html)
        self.assertIn("Curated Complex Portal records", no_hits_html)
        self.assertNotIn("Complex Portal context", no_hits_html)
        self.assertIn("No curated human Complex Portal records matched this target.", no_hits_html)

        report["complex_portal_lookup"] = {"human": {"status": "ok"}}
        report["complex_portal_complexes"] = [
            {
                "complex_ac": "CPX-TEST",
                "name": "Test alpha/beta complex",
                "species": "Homo sapiens; 9606",
                "complex_assemblies": ["Heterodimer"],
                "participants": [{"identifier": "P08648", "name": "ITGA5"}, {"identifier": "P05556", "name": "ITGB1"}],
            }
        ]
        report["assembly_requirements"] = [
            {
                "classification": "obligatory_partner_requirement",
                "obligatory": True,
                "confidence": "high",
                "summary": "Requires an integrin beta partner.",
                "partners": ["ITGB1"],
            },
            {
                "classification": "interaction_context",
                "obligatory": False,
                "confidence": "medium",
                "summary": "Noisy non-obligatory interaction context.",
                "partners": ["HIV"],
            },
        ]
        with TemporaryDirectory() as root:
            records_html = self._render(report, Path(root))
        self.assertIn("CPX-TEST", records_html)
        self.assertIn("ITGA5, ITGB1", records_html)
        self.assertIn("Assembly requirements", records_html)
        self.assertIn("Requires an integrin beta partner.", records_html)
        self.assertNotIn("Noisy non-obligatory interaction context.", records_html)


class ConstructNameLinkTests(unittest.TestCase):
    def test_pdb_construct_name_links_to_rcsb(self) -> None:
        html = _construct_name_html({"name": "pdb_1ira_x", "pdb_id": "1IRA"})
        self.assertEqual(
            html,
            '<a class="inline-link" href="https://www.rcsb.org/structure/1IRA" '
            'target="_blank" rel="noreferrer">pdb_1ira_x</a>',
        )

    def test_non_pdb_construct_name_is_plain_escaped_text(self) -> None:
        self.assertEqual(_construct_name_html({"name": "full_ectodomain", "pdb_id": None}), "full_ectodomain")

    def test_malicious_pdb_id_is_neutralised(self) -> None:
        html = _construct_name_html({"name": "pdb_x", "pdb_id": '"><script>alert(1)</script>'})
        self.assertNotIn("<script>", html)
        self.assertNotIn('"><', html)


if __name__ == "__main__":
    unittest.main()
