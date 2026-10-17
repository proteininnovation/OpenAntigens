"""The static portal ships a gzip-enabling .htaccess that survives the
public_site packaging step (which otherwise strips dotfiles)."""

from __future__ import annotations

import io
import shutil
import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from agdesign2.portal import _write_portal_htaccess


class PortalHtaccessTests(unittest.TestCase):
    def test_write_portal_htaccess_enables_gzip_for_text(self) -> None:
        with TemporaryDirectory() as tmp:
            portal = Path(tmp)
            _write_portal_htaccess(portal)
            text = (portal / ".htaccess").read_text(encoding="utf-8")
            self.assertIn("mod_deflate", text)
            self.assertIn("application/javascript", text)
            # Large JSON/TSV downloads are intentionally NOT compressed.
            self.assertNotIn("application/json", text)

    def test_only_the_mouse_portal_carries_the_removed_page_redirects(self) -> None:
        with TemporaryDirectory() as tmp:
            human, mouse = Path(tmp) / "human", Path(tmp) / "mouse"
            human.mkdir()
            mouse.mkdir()
            _write_portal_htaccess(human)
            _write_portal_htaccess(mouse, mouse_redirects=True)
            self.assertNotIn("RedirectMatch", (human / ".htaccess").read_text(encoding="utf-8"))
            text = (mouse / ".htaccess").read_text(encoding="utf-8")
            # ccl12_mouse is the case the remediation plan names: PR-9's winner
            # rule moved the page to ccl2_mouse.
            self.assertIn(
                "RedirectMatch 301 ^/mouse/reports/ccl12_mouse\\.html$ /mouse/reports/ccl2_mouse.html",
                text,
            )

    def test_the_precompressed_rewrite_keeps_the_redirects(self) -> None:
        # compress_public_site rewrites every .htaccess in place; before the
        # rules moved into the writer, that silently dropped them.
        with TemporaryDirectory() as tmp:
            mouse = Path(tmp) / "mouse"
            mouse.mkdir()
            _write_portal_htaccess(mouse, mouse_redirects=True)
            plain = (mouse / ".htaccess").read_text(encoding="utf-8").count("RedirectMatch 301")
            _write_portal_htaccess(mouse, precompressed=True, mouse_redirects=True)
            text = (mouse / ".htaccess").read_text(encoding="utf-8")
            self.assertEqual(text.count("RedirectMatch 301"), plain)
            self.assertIn("Content-Encoding", text)

    def test_every_redirect_target_is_a_page_the_map_does_not_remove(self) -> None:
        # A target that is itself redirected would chain; a duplicate source
        # would make the second rule dead.
        import csv
        from importlib import resources

        from agdesign2.portal import _MOUSE_REDIRECT_MAP

        text = resources.files("agdesign2.data").joinpath(_MOUSE_REDIRECT_MAP).read_text(encoding="utf-8")
        rows = [row for row in csv.DictReader(io.StringIO(text), delimiter="\t") if row.get("removed_page")]
        sources = [row["removed_page"] for row in rows]
        self.assertEqual(len(sources), len(set(sources)))
        removed = {f"/mouse/reports/{source}" for source in sources}
        self.assertEqual([row["redirect_to"] for row in rows if row["redirect_to"] in removed], [])

    def test_package_public_site_preserves_htaccess(self) -> None:
        from agdesign2.deployment import _package_public_site

        with TemporaryDirectory() as tmp:
            portal = Path(tmp) / "portal"
            portal.mkdir()
            (portal / "index.html").write_text("<html></html>", encoding="utf-8")
            _write_portal_htaccess(portal)
            # A different dotfile must still be stripped.
            (portal / ".secret").write_text("nope", encoding="utf-8")

            public = Path(tmp) / "public_site"
            _package_public_site(portal, public)

            self.assertTrue((public / ".htaccess").exists())
            self.assertFalse((public / ".secret").exists())


if __name__ == "__main__":
    unittest.main()
