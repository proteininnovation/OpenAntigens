"""The static portal ships a gzip-enabling .htaccess that survives the
public_site packaging step (which otherwise strips dotfiles)."""

from __future__ import annotations

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
