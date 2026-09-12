"""Tests for the public_site disk-size reduction pass (GoDaddy quota fit)."""
from __future__ import annotations

import gzip
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from agdesign2.deployment import compress_public_site
from agdesign2.portal import _PORTAL_HTACCESS, _PRECOMPRESSED_PORTAL_HTACCESS


def _make_site(root: Path) -> None:
    (root / "index.html").write_text("<html>root index</html>", encoding="utf-8")
    (root / "portal-index-data.js").write_text("var INDEX = " + "x" * 5000 + ";", encoding="utf-8")
    (root / "styles.css").write_text("body{color:red}" * 200, encoding="utf-8")
    (root / "logo.svg").write_text("<svg>" + "p" * 1000 + "</svg>", encoding="utf-8")
    (root / "portal_metadata.json").write_text('{"k": 1}', encoding="utf-8")
    rs = root / "report_scripts"
    rs.mkdir()
    (rs / "egfr_human.js").write_text("const PDB = '" + "A" * 8000 + "';", encoding="utf-8")
    reports = root / "reports"
    reports.mkdir()
    (reports / "egfr_human.html").write_text("<html>" + "r" * 4000 + "</html>", encoding="utf-8")
    assets = root / "report_assets"
    assets.mkdir()
    (assets / "egfr.png").write_bytes(b"\x89PNG\r\n" + b"\x00" * 2000)


class CompressPublicSiteTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp_path = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.tmp_path, ignore_errors=True)

    def test_compresses_js_css_svg_under_public_filenames(self) -> None:
        site = self.tmp_path / "public_site"
        site.mkdir()
        _make_site(site)
        original_js = (site / "report_scripts" / "egfr_human.js").read_bytes()

        summary = compress_public_site(site)

        # Public filenames remain stable while their contents become gzip data.
        report_js = site / "report_scripts" / "egfr_human.js"
        assert report_js.exists()
        assert not (site / "report_scripts" / "egfr_human.js.gz").exists()
        assert (site / "portal-index-data.js").exists()
        assert (site / "styles.css").exists()
        assert (site / "logo.svg").exists()
        assert gzip.decompress(report_js.read_bytes()) == original_js
        # HTML, PNG, the directory index, and JSON metadata are left untouched.
        assert (site / "reports" / "egfr_human.html").exists()
        assert (site / "report_assets" / "egfr.png").exists()
        assert (site / "index.html").exists()
        assert (site / "portal_metadata.json").exists()
        assert summary["files_compressed"] == 4
        assert summary["bytes_before"] > summary["bytes_after"] > 0

    def test_idempotent(self) -> None:
        site = self.tmp_path / "public_site"
        site.mkdir()
        _make_site(site)

        first = compress_public_site(site)
        assert first["files_compressed"] == 4
        second = compress_public_site(site)
        assert second["files_compressed"] == 0
        assert gzip.decompress((site / "report_scripts" / "egfr_human.js").read_bytes()).startswith(b"const PDB")

    def test_interrupted_compression_preserves_original(self) -> None:
        site = self.tmp_path / "public_site"
        site.mkdir()
        _make_site(site)
        path = site / "styles.css"
        original = path.read_bytes()

        with patch("agdesign2.deployment.shutil.copyfileobj", side_effect=OSError("stopped")):
            with self.assertRaises(OSError):
                compress_public_site(site)

        assert path.read_bytes() == original
        assert not path.with_name(path.name + ".gzip.tmp").exists()

    def test_drop_unreferenced_structures(self) -> None:
        site = self.tmp_path / "public_site"
        site.mkdir()
        _make_site(site)
        for rel in ("structures", "mouse/structures"):
            directory = site / rel
            directory.mkdir(parents=True)
            (directory / "egfr_human.pdb").write_text("ATOM  ...\n", encoding="utf-8")
            (directory / "egfr_human.meta.json").write_text("{}", encoding="utf-8")

        summary = compress_public_site(site, drop_unreferenced_structures=True)

        assert not (site / "structures").exists()
        assert not (site / "mouse" / "structures").exists()
        assert summary["structures_bytes_removed"] > 0

    def test_drop_redundant_downloads(self) -> None:
        site = self.tmp_path / "public_site"
        site.mkdir()
        _make_site(site)
        (site / "open_targets_disease_associations.json").write_text("[]", encoding="utf-8")
        (site / "open_targets_disease_associations.tsv").write_text("a\tb\n", encoding="utf-8")
        (site / "downloads.html").write_text(
            "<div>\n"
            '          <a class="button" href="open_targets_disease_associations.tsv">Open Targets TSV</a>\n'
            '          <a class="button" href="open_targets_disease_associations.json">Open Targets JSON</a>\n'
            "</div>\n",
            encoding="utf-8",
        )

        summary = compress_public_site(site, drop_redundant_downloads=True)

        assert not (site / "open_targets_disease_associations.json").exists()
        # The TSV holds the same data and stays a plain wget-friendly download.
        assert (site / "open_targets_disease_associations.tsv").exists()
        downloads = (site / "downloads.html").read_text(encoding="utf-8")
        assert "open_targets_disease_associations.json" not in downloads  # dead button stripped
        assert "open_targets_disease_associations.tsv" in downloads  # TSV button retained
        assert summary["downloads_bytes_removed"] > 0

    def test_refreshes_root_and_mouse_htaccess(self) -> None:
        site = self.tmp_path / "public_site"
        site.mkdir()
        _make_site(site)
        # A snapshot packaged before the precompressed rules existed carries a stale
        # .htaccess at the root *and* in mouse/ (mod_rewrite does not inherit into a
        # subdir that has its own .htaccess), so both must be refreshed.
        stale = "# old hosting hints, no rewrite rules\n"
        (site / ".htaccess").write_text(stale, encoding="utf-8")
        mouse = site / "mouse"
        mouse.mkdir()
        (mouse / "index.html").write_text("<html>mouse</html>", encoding="utf-8")
        (mouse / ".htaccess").write_text(stale, encoding="utf-8")

        compress_public_site(site)

        for htaccess in (site / ".htaccess", mouse / ".htaccess"):
            text = htaccess.read_text(encoding="utf-8")
            assert "RewriteRule" not in text
            assert '<FilesMatch "\\.js$">' in text
            assert "Header set Content-Encoding gzip" in text

    def test_requires_packaged_public_site(self) -> None:
        empty = self.tmp_path / "not_a_site"
        empty.mkdir()
        with self.assertRaises(FileNotFoundError):
            compress_public_site(empty)

    def test_plain_and_precompressed_htaccess_are_separate(self) -> None:
        assert "Content-Encoding gzip" not in _PORTAL_HTACCESS
        assert "RewriteRule" not in _PORTAL_HTACCESS

        h = _PRECOMPRESSED_PORTAL_HTACCESS
        assert "RewriteRule" not in h
        assert '<FilesMatch "\\.js$">' in h
        assert '<FilesMatch "\\.css$">' in h
        assert '<FilesMatch "\\.svg$">' in h
        assert "Header set Content-Encoding gzip" in h
        assert "Header append Vary Accept-Encoding" in h
        assert "SetEnv no-gzip 1" in h


if __name__ == "__main__":
    unittest.main()
