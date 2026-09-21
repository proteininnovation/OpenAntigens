from __future__ import annotations

import json
import os
import shutil
import sys
from unittest import mock
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEPLOY_SCRIPT = ROOT / "scripts" / "deploy_snapshot_to_godaddy.sh"


class GoDaddyDeployScriptTests(unittest.TestCase):
    def _fixture(self, *, ssh_exit: int = 0) -> tuple[tempfile.TemporaryDirectory[str], Path, Path, dict[str, str]]:
        tmp = tempfile.TemporaryDirectory()
        root = Path(tmp.name)
        snapshot = root / "snapshot"
        public_site = snapshot / "public_site"
        public_site.mkdir(parents=True)
        (public_site / "index.html").write_text(
            '<link rel="stylesheet" href="portal.css?v=test">\n'
            '<script src="portal-index-data.js?v=test"></script>\n'
            '<script src="portal-disease-index.js?v=test"></script>\n'
            '<script src="portal-index.js?v=test"></script>\n',
            encoding="utf-8",
        )
        human_row = {
            "status": "ok",
            "detail_page": "reports/example.html",
            "cross_reactivity_count": 1,
            "top_disease_name": "Example disease",
            "disease_count": 1,
            "has_alphafold_structure": 1,
        }
        mouse_row = {
            "status": "ok",
            "detail_page": "reports/example.html",
            "cross_reactivity_count": 1,
            "has_alphafold_structure": 1,
        }
        (public_site / "downloads").mkdir()
        (public_site / "mouse/downloads").mkdir(parents=True)
        (public_site / "reports").mkdir()
        (public_site / "mouse/reports").mkdir()
        (public_site / "mouse/report_scripts").mkdir()
        (public_site / "reports/example.html").write_text("human", encoding="utf-8")
        (public_site / "mouse/reports/example.html").write_text("mouse", encoding="utf-8")
        (public_site / "mouse/report_scripts/example.js").write_text("const report = {};", encoding="utf-8")
        (public_site / "portal_metadata.json").write_text(
            json.dumps({"target_count": 5_000}) + "\n", encoding="utf-8"
        )
        (public_site / "downloads/agdesign2_portal_index.json").write_text(
            json.dumps([human_row] * 5_000), encoding="utf-8"
        )
        (public_site / "mouse/portal_metadata.json").write_text(
            json.dumps({"target_count": 4_500}) + "\n", encoding="utf-8"
        )
        (public_site / "mouse/downloads/agdesign2_portal_index.json").write_text(
            json.dumps([mouse_row] * 4_500), encoding="utf-8"
        )
        (public_site / "open_targets_disease_associations.tsv").write_text(
            "entry_name\tdisease_name\topen_targets_url\n"
            + "EXAMPLE\tExample disease\thttps://example.test\n" * 5_000,
            encoding="utf-8",
        )
        (public_site / "portal.css").write_text("body{}", encoding="utf-8")
        (public_site / "portal-index-data.js").write_text("const rows = [];", encoding="utf-8")
        (public_site / "portal-disease-index.js").write_text("const diseases = [];", encoding="utf-8")
        (public_site / "portal-index.js").write_text("const index = [];", encoding="utf-8")
        report_scripts = public_site / "report_scripts"
        report_scripts.mkdir()
        (report_scripts / "example.js").write_text("const report = {};", encoding="utf-8")
        (report_scripts / "egfr_human.js").write_text("const report = {};", encoding="utf-8")

        marker_dir = root / "markers"
        marker_dir.mkdir()
        bin_dir = root / "bin"
        bin_dir.mkdir()
        identity_file = root / "deploy-key"
        identity_file.touch(mode=0o600)
        self._write_command(
            bin_dir / "ssh",
            f'printf "%s\\n" "$@" > "$TEST_MARKER_DIR/ssh_args"\nexit {ssh_exit}\n',
        )
        self._write_command(
            bin_dir / "curl",
            'args="$*"\n'
            'printf "%s\\n" "$args" >> "$TEST_MARKER_DIR/curl_args"\n'
            'if [[ -n "${TEST_CURL_FAIL_ASSET:-}" && "$args" == *"$TEST_CURL_FAIL_ASSET"* ]]; then exit 22; fi\n'
            'if [[ "$args" == *"--config -"* ]]; then input="$(cat)"; printf "%s\\n" "$input" >> "$TEST_MARKER_DIR/curl_args"; printf \'{"status":1}\'; exit 0; fi\n'
            'output=""; url=""\n'
            'while [[ "$#" -gt 0 ]]; do case "$1" in --output) output="$2"; shift 2;; http*) url="$1"; shift;; *) shift;; esac; done\n'
            'if [[ -n "$output" && "$output" != "/dev/null" && -n "$url" ]]; then rel="${url#*://}"; rel="${rel#*/}"; cp "$TEST_PUBLIC_SITE/$rel" "$output"; fi\n',
        )
        self._write_command(
            bin_dir / "agdesign2",
            'touch "$TEST_MARKER_DIR/agdesign2"\n',
        )
        self._write_command(
            bin_dir / "rsync",
            'printf "%s\\n" "$@" > "$TEST_MARKER_DIR/rsync_args"\n',
        )

        env = os.environ.copy()
        env.update(
            {
                "PATH": f"{bin_dir}:{env['PATH']}",
                "GODADDY_TARGET": "deploy@example:public_html/",
                "GODADDY_IDENTITY_FILE": str(identity_file),
                "AGDESIGN2_BIN": str(bin_dir / "agdesign2"),
                "SUCURI_ENV_FILE": str(root / "missing-sucuri.env"),
                "TEST_MARKER_DIR": str(marker_dir),
                "TEST_PUBLIC_SITE": str(public_site),
            }
        )
        return tmp, snapshot, marker_dir, env

    @staticmethod
    def _write_command(path: Path, body: str) -> None:
        path.write_text("#!/usr/bin/env bash\nset -euo pipefail\n" + body, encoding="utf-8")
        path.chmod(0o755)

    @staticmethod
    def _run(snapshot: Path, env: dict[str, str], *args: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [str(DEPLOY_SCRIPT), "--existing-snapshot", str(snapshot), *args],
            cwd=ROOT,
            env=env,
            capture_output=True,
            text=True,
            check=False,
        )

    def test_missing_sucuri_credentials_fail_before_ssh_or_compression(self) -> None:
        tmp, snapshot, marker_dir, env = self._fixture()
        self.addCleanup(tmp.cleanup)

        result = self._run(snapshot, env)

        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Sucuri credentials are required", result.stderr)
        self.assertFalse((marker_dir / "ssh_args").exists())
        self.assertFalse((marker_dir / "agdesign2").exists())
        self.assertFalse((marker_dir / "rsync_args").exists())

    def test_ssh_failure_stops_before_compression_or_rsync(self) -> None:
        tmp, snapshot, marker_dir, env = self._fixture(ssh_exit=23)
        self.addCleanup(tmp.cleanup)
        env.update({"SUCURI_API_KEY": "test-key", "SUCURI_API_SECRET": "test-secret"})

        result = self._run(snapshot, env)

        self.assertNotEqual(result.returncode, 0)
        self.assertIn("non-interactive SSH preflight failed or the publish root is not writable", result.stderr)
        ssh_args = (marker_dir / "ssh_args").read_text(encoding="utf-8")
        self.assertIn(str(Path(env["GODADDY_IDENTITY_FILE"])), ssh_args)
        self.assertIn("IdentitiesOnly=yes", ssh_args)
        self.assertIn("test -d public_html/", ssh_args)
        self.assertIn("test -w public_html/", ssh_args)
        self.assertFalse((marker_dir / "agdesign2").exists())
        self.assertFalse((marker_dir / "rsync_args").exists())

    def test_dry_run_uses_safe_mirror_arguments(self) -> None:
        tmp, snapshot, marker_dir, env = self._fixture()
        self.addCleanup(tmp.cleanup)
        env.update(
            {
                "COMPRESS": "0",
                "PURGE_SUCURI": "0",
            }
        )

        result = self._run(snapshot, env, "--dry-run")

        self.assertEqual(result.returncode, 0, result.stderr)
        rsync_args = (marker_dir / "rsync_args").read_text(encoding="utf-8")
        self.assertIn("--dry-run", rsync_args)
        self.assertIn("--itemize-changes", rsync_args)
        self.assertIn("--stats", rsync_args)
        self.assertIn("--delete-before", rsync_args)
        self.assertIn("--exclude=.ftpquota", rsync_args)
        self.assertIn("--exclude=.well-known/", rsync_args)
        self.assertIn(str(Path(env["GODADDY_IDENTITY_FILE"])), rsync_args)
        self.assertIn("IdentitiesOnly=yes", rsync_args)
        ssh_args = (marker_dir / "ssh_args").read_text(encoding="utf-8")
        self.assertIn(str(Path(env["GODADDY_IDENTITY_FILE"])), ssh_args)
        self.assertIn("IdentitiesOnly=yes", ssh_args)
        self.assertIn("test -d public_html/", ssh_args)
        self.assertIn("test -w public_html/", ssh_args)
        self.assertIn("dry-run complete; no remote files were changed", result.stdout)
        self.assertFalse((marker_dir / "agdesign2").exists())

    def test_legacy_gzip_suffix_stops_before_rsync(self) -> None:
        tmp, snapshot, marker_dir, env = self._fixture()
        self.addCleanup(tmp.cleanup)
        env.update({"COMPRESS": "0", "PURGE_SUCURI": "0"})
        (snapshot / "public_site" / "legacy.js.gz").write_bytes(b"gzip")

        result = self._run(snapshot, env)

        self.assertNotEqual(result.returncode, 0)
        self.assertIn("obsolete .gz suffix", result.stderr)
        self.assertFalse((marker_dir / "rsync_args").exists())

    def test_missing_blast_coverage_stops_before_rsync(self) -> None:
        tmp, snapshot, marker_dir, env = self._fixture()
        self.addCleanup(tmp.cleanup)
        env.update({"COMPRESS": "0", "PURGE_SUCURI": "0"})
        index = snapshot / "public_site/downloads/agdesign2_portal_index.json"
        rows = json.loads(index.read_text(encoding="utf-8"))
        for row in rows:
            row["cross_reactivity_count"] = 0
        index.write_text(json.dumps(rows), encoding="utf-8")

        result = self._run(snapshot, env)

        self.assertNotEqual(result.returncode, 0)
        self.assertIn("human BLAST coverage is 0/5000", result.stderr)
        self.assertFalse((marker_dir / "rsync_args").exists())

    def test_missing_open_targets_stops_before_rsync(self) -> None:
        tmp, snapshot, marker_dir, env = self._fixture()
        self.addCleanup(tmp.cleanup)
        env.update({"COMPRESS": "0", "PURGE_SUCURI": "0"})
        (snapshot / "public_site/open_targets_disease_associations.tsv").unlink()

        result = self._run(snapshot, env)

        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Open Targets association TSV is missing or empty", result.stderr)
        self.assertFalse((marker_dir / "rsync_args").exists())

    def test_publish_purges_before_live_asset_checks(self) -> None:
        tmp, snapshot, marker_dir, env = self._fixture()
        self.addCleanup(tmp.cleanup)
        env.update(
            {
                "COMPRESS": "0",
                "SUCURI_API_KEY": "test-key",
                "SUCURI_API_SECRET": "test-secret",
            }
        )

        result = self._run(snapshot, env)

        self.assertEqual(result.returncode, 0, result.stderr)
        calls = (marker_dir / "curl_args").read_text(encoding="utf-8")
        self.assertLess(calls.index("waf.sucuri.net"), calls.index("portal.css?v=test"))
        self.assertIn("report_scripts/egfr_human.js", calls)
        self.assertIn("live asset verification passed", result.stdout)

    def test_live_asset_failure_fails_publish(self) -> None:
        tmp, snapshot, _marker_dir, env = self._fixture()
        self.addCleanup(tmp.cleanup)
        env.update(
            {
                "COMPRESS": "0",
                "SUCURI_API_KEY": "test-key",
                "SUCURI_API_SECRET": "test-secret",
                "TEST_CURL_FAIL_ASSET": "portal.css",
            }
        )

        result = self._run(snapshot, env)

        self.assertNotEqual(result.returncode, 0)
        self.assertNotIn("live asset verification passed", result.stdout)

    def test_default_module_entrypoint_compresses_without_copied_cli_shebang(self) -> None:
        tmp, snapshot, marker_dir, env = self._fixture()
        self.addCleanup(tmp.cleanup)
        env.pop("AGDESIGN2_BIN")
        env.update({"PURGE_SUCURI": "0"})

        checkout = Path(tmp.name) / "checkout"
        (checkout / "scripts").mkdir(parents=True)
        (checkout / ".venv/bin").mkdir(parents=True)
        (checkout / ".venv/bin/python").symlink_to(sys.executable)
        (checkout / "src").symlink_to(ROOT / "src", target_is_directory=True)
        script = checkout / "scripts" / DEPLOY_SCRIPT.name
        shutil.copy2(DEPLOY_SCRIPT, script)
        with mock.patch(__name__ + ".DEPLOY_SCRIPT", script):
            result = self._run(snapshot, env)

        self.assertEqual(result.returncode, 0, result.stderr)
        compressed = snapshot / "public_site" / "portal.css"
        self.assertEqual(compressed.read_bytes()[:2], b"\x1f\x8b")
        self.assertFalse((snapshot / "public_site" / "portal.css.gz").exists())
        self.assertFalse((marker_dir / "agdesign2").exists())


if __name__ == "__main__":
    unittest.main()
