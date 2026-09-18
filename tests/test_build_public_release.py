from __future__ import annotations

import json
import gzip
import tempfile
import unittest
from pathlib import Path

from scripts.build_public_release import (
    REQUIRED_PORTAL_DIRECTORIES,
    REQUIRED_PORTAL_FILES,
    build_manifest,
    copy_public_license_files,
    copy_public_portal,
    read_portal_metadata,
    resolve_release_dir,
    validate_portal,
)


def _write_minimal_portal(portal_dir: Path) -> None:
    for name in REQUIRED_PORTAL_DIRECTORIES:
        (portal_dir / name).mkdir(parents=True, exist_ok=True)
    for name in REQUIRED_PORTAL_FILES:
        path = portal_dir / name
        path.parent.mkdir(parents=True, exist_ok=True)
        if name == "portal-index-data.js":
            content = "https://www.ncbi.nlm.nih.gov/research/pubtator3/docsum?text=%40GENE_TEST"
        elif name == "downloads/agdesign2_portal_index.json":
            content = '[{"has_alphafold_structure": 1}]'
        else:
            content = "{}" if path.suffix == ".json" else "content"
        path.write_text(content, encoding="utf-8")


class PublicReleaseTests(unittest.TestCase):
    def test_portal_validation_requires_complete_nonempty_output(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            portal_dir = Path(tmpdir) / "portal"
            portal_dir.mkdir()
            for name in ("index.html", "portal.css", "help.html", "methods.html", "downloads.html", "terms.html"):
                (portal_dir / name).write_text("", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "missing or empty required files"):
                validate_portal(portal_dir)

    def test_portal_validation_requires_privacy_page(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            portal_dir = Path(tmpdir) / "portal"
            _write_minimal_portal(portal_dir)
            (portal_dir / "privacy.html").unlink()

            with self.assertRaisesRegex(ValueError, "privacy.html"):
                validate_portal(portal_dir)

    def test_portal_validation_rejects_missing_pubtator_values(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            portal_dir = Path(tmpdir) / "portal"
            _write_minimal_portal(portal_dir)
            (portal_dir / "portal-index-data.js").write_text(
                'window.OpenAntigenIndexRows = [{"dataset":{"pubtator":"0"},'
                '"cells":["<span class=\\"muted\\">n/a</span>"]}];',
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "no populated PubTator links"):
                validate_portal(portal_dir)

    def test_portal_validation_rejects_low_local_structure_coverage(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            portal_dir = Path(tmpdir) / "portal"
            _write_minimal_portal(portal_dir)
            rows = [{"has_alphafold_structure": int(index < 790)} for index in range(1000)]
            (portal_dir / "downloads/agdesign2_portal_index.json").write_text(
                json.dumps(rows), encoding="utf-8"
            )

            with self.assertRaisesRegex(ValueError, "local AlphaFold coverage is 790/1000"):
                validate_portal(portal_dir)

    def test_portal_symlinks_cannot_copy_outside_content(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            portal_dir = root / "portal"
            _write_minimal_portal(portal_dir)
            outside = root / "outside-secret.txt"
            outside.write_text("PRIVATE", encoding="utf-8")
            (portal_dir / "reports" / "leak.txt").symlink_to(outside)

            with self.assertRaisesRegex(ValueError, "symbolic links"):
                copy_public_portal(portal_dir, root / "release")
            self.assertFalse((root / "release").exists())

    def test_public_copy_keeps_apache_configuration(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            portal_dir = root / "portal"
            release_dir = root / "release"
            _write_minimal_portal(portal_dir)

            copy_public_portal(portal_dir, release_dir)

            self.assertEqual((release_dir / ".htaccess").read_text(encoding="utf-8"), "content")

    def test_public_copy_includes_only_active_license_files(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            release_dir = root / "release"
            release_dir.mkdir()
            for name in ("LICENSE", "DATA_LICENSE.md", "LICENSE.AlphaFold-Server-Output-Terms.md"):
                (root / name).write_text(name, encoding="utf-8")

            copy_public_license_files(root, release_dir)

            self.assertTrue((release_dir / "LICENSE").is_file())
            self.assertTrue((release_dir / "DATA_LICENSE.md").is_file())
            self.assertFalse((release_dir / "LICENSE.AlphaFold-Server-Output-Terms.md").exists())

    def test_portal_validation_rejects_af3_content(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            portal_dir = Path(tmpdir) / "portal"
            _write_minimal_portal(portal_dir)
            (portal_dir / "reports" / "target.html").write_text(
                "Using an AlphaFold 3 model",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "excluded AF3 content"):
                validate_portal(portal_dir)

    def test_compressed_assets_preserve_public_release_validation(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            portal_dir = Path(tmpdir) / "portal"
            _write_minimal_portal(portal_dir)
            index = portal_dir / "portal-index-data.js"
            index.write_bytes(gzip.compress(index.read_bytes()))
            validate_portal(portal_dir)
            script = portal_dir / "report_scripts" / "target.js"
            script.parent.mkdir(exist_ok=True)
            script.write_bytes(gzip.compress(b'const source = "AlphaFold 3";'))
            with self.assertRaisesRegex(ValueError, "excluded AF3 content"):
                validate_portal(portal_dir)

    def test_portal_validation_rejects_static_pymol_images(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            portal_dir = Path(tmpdir) / "portal"
            _write_minimal_portal(portal_dir)
            (portal_dir / "report_assets").mkdir()
            (portal_dir / "report_assets" / "target_structure.png").write_bytes(b"png")

            with self.assertRaisesRegex(ValueError, "excluded static PyMOL"):
                validate_portal(portal_dir)

    def test_release_directory_cannot_escape_output_root(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            output_root = Path(tmpdir) / "releases"
            self.assertEqual(
                resolve_release_dir(output_root, "openantigen-test"),
                (output_root / "openantigen-test").resolve(),
            )
            for name in ("", ".", "..", "../outside", "nested/release"):
                with self.subTest(name=name), self.assertRaises(ValueError):
                    resolve_release_dir(output_root, name)

    def test_manifest_does_not_publish_local_paths(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            release_dir = root / "releases" / "openantigen-test"
            release_dir.mkdir(parents=True)
            (release_dir / "index.html").write_text("", encoding="utf-8")

            manifest = build_manifest(repo_root=root, release_dir=release_dir, checksums=False)

        self.assertNotIn("source_portal_dir", manifest)
        self.assertNotIn("release_root", manifest)
        self.assertNotIn(tmpdir, json.dumps(manifest))

    def test_invalid_portal_metadata_is_an_error(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            release_dir = Path(tmpdir)
            (release_dir / "portal_metadata.json").write_text("[]", encoding="utf-8")
            with self.assertRaises(ValueError):
                read_portal_metadata(release_dir)


if __name__ == "__main__":
    unittest.main()
