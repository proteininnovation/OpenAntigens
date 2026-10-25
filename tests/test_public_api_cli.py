"""CLI contracts for frozen public API compilation and validation."""
from __future__ import annotations

from contextlib import redirect_stdout
import io
from pathlib import Path
import unittest
from unittest.mock import patch

from agdesign2.cli import build_parser, main


class PublicApiCliTests(unittest.TestCase):
    def test_compile_requires_both_species_and_explicit_provenance(self):
        args = build_parser().parse_args([
            "build-public-api", "--human-summary", "human.json", "--mouse-summary", "mouse.json",
            "--public-site", "site", "--output-dir", "bundle", "--software-commit", "a" * 40,
            "--snapshot-id", "frozen-release", "--previous-bundle", "previous",
        ])
        self.assertEqual(args.human_summary, "human.json")
        self.assertEqual(args.mouse_summary, "mouse.json")
        self.assertEqual(args.previous_bundle, "previous")

    def test_compile_dispatches_to_shared_builder(self):
        argv = ["agdesign2", "build-public-api", "--human-summary", "human.json",
                "--mouse-summary", "mouse.json", "--public-site", "site", "--output-dir", "bundle",
                "--software-commit", "a" * 40, "--snapshot-id", "frozen-release"]
        with patch("sys.argv", argv), patch("agdesign2.public_api_data.build_public_api_bundle", return_value=Path("bundle/manifest.json")) as build:
            with redirect_stdout(io.StringIO()) as output:
                main()
        build.assert_called_once_with(
            "human.json", "mouse.json", "bundle", public_site_dir="site", software_commit="a" * 40,
            snapshot_id="frozen-release", previous_bundle_dir=None, static_origin="https://openantigens.org",
        )
        self.assertEqual(output.getvalue().strip(), "bundle/manifest.json")

    def test_validate_dispatches_to_bundle_validator(self):
        with patch("sys.argv", ["agdesign2", "validate-public-api", "bundle"]), patch(
            "agdesign2.public_api_data.validate_public_api_bundle", return_value={"current_release_id": "a" * 64}
        ) as validate:
            with redirect_stdout(io.StringIO()) as output:
                main()
        validate.assert_called_once_with("bundle")
        self.assertEqual(output.getvalue().strip(), "valid release=" + "a" * 64)
