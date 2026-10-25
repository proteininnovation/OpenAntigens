"""Exercise the actual public bundle validator and the final deployment ZIP."""
from __future__ import annotations

import hashlib
import io
import importlib.util
import json
from pathlib import Path
import shutil
import stat
import zipfile

import tempfile
import unittest
from unittest.mock import patch
from contextlib import redirect_stdout

from agdesign2.public_api_data import build_public_api_bundle
from agdesign2.public_api_release import SOURCE_FILES, build_api_service_release

REPO = Path(__file__).resolve().parents[1]


def load_script(name: str):
    spec = importlib.util.spec_from_file_location(name, REPO / f"scripts/{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def build_frozen(root):
    spec = importlib.util.spec_from_file_location("public_api_fixture", REPO / "tests/build_public_api_fixture.py")
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    helper.build(root)
    service = root / "service"
    for name in SOURCE_FILES:
        destination = service / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(REPO / "services/public-api" / name, destination)
    (service / ".env").write_text("PRIVATE_ENV_CANARY_DO_NOT_PUBLISH")
    (service / "node_modules").mkdir()
    (service / "node_modules/private.js").write_text("PRIVATE_MODULE_CANARY_DO_NOT_PUBLISH")


def package(root, output, **kwargs):
    return build_api_service_release(root / "service", root / "bundle", root / "site", output, repo_root=REPO, **kwargs)


def snapshot(root):
    return {path.relative_to(root).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest() for path in root.rglob("*") if path.is_file() and not path.is_symlink()}


def write_root_manifest(bundle, value):
    content = (json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")) + "\n").encode()
    (bundle / "manifest.json").write_bytes(content)
    (bundle / "manifest.sha256").write_text(hashlib.sha256(content).hexdigest() + "  manifest.json\n")


def replace_release_manifest(bundle, site, root, content, pointer_index=0):
    pointer = root["releases"][pointer_index]
    pointer["sha256"] = hashlib.sha256(content).hexdigest()
    (bundle / pointer["manifest_path"]).write_bytes(content)
    static = site / f"data/releases/{pointer['release_id']}/manifest.json"
    static.write_bytes(content)
    static.with_suffix(".sha256").write_text(pointer["sha256"] + "  manifest.json\n")
    write_root_manifest(bundle, root)


class PublicApiReleaseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.frozen_temp = tempfile.TemporaryDirectory(prefix="oa-release-source-")
        cls.frozen = Path(cls.frozen_temp.name).resolve()
        build_frozen(cls.frozen)

    @classmethod
    def tearDownClass(cls):
        cls.frozen_temp.cleanup()

    def setUp(self):
        self.scratch = tempfile.TemporaryDirectory(prefix="oa-release-test-")
        self.addCleanup(self.scratch.cleanup)
        self.tmp_path = Path(self.scratch.name).resolve()
        self.release = self.tmp_path / "inputs"
        shutil.copytree(self.frozen, self.release)
        for species in ("human", "mouse"):
            summary = self.release / species / "batch_summary.json"
            rows = json.loads(summary.read_text())
            for row in rows:
                row["json_report"] = str(self.release / species / Path(row["json_report"]).name)
            summary.write_text(json.dumps(rows))

    def test_real_zip_is_deterministic_allowlisted_and_measured(self):
        release, tmp_path = self.release, self.tmp_path
        before = snapshot(release)
        output = tmp_path / "release.zip"
        dry = package(release, output, dry_run=True)
        assert not output.exists() and not output.with_suffix(".artifact.json").exists()
        assert snapshot(release) == before

        artifact = package(release, output)
        duplicate = package(release, tmp_path / "duplicate.zip")
        assert artifact["sha256"] == dry["sha256"] == duplicate["sha256"]
        assert artifact["sha256"] == hashlib.sha256(output.read_bytes()).hexdigest()
        assert artifact["bytes"] == output.stat().st_size < 100_000_000
        assert artifact["local_budgets"]["host_capacity_verified"] is False
        assert artifact["gates"]["node22_runtime_and_port"] == "not_verified"
        assert artifact["dependency_audit"]["status"] == "not_run"
        assert output.with_suffix(".artifact.json").is_file()
        with zipfile.ZipFile(output) as archive:
            names = archive.namelist()
            assert names == sorted(names) == sorted(artifact["files"])
            assert set(name for name in names if not name.startswith("data/")) == {*SOURCE_FILES, "LICENSE", "DATA_LICENSE.md", "service-build.json"}
            assert artifact["expanded_bytes"] == sum(entry.file_size for entry in archive.infolist())
            for entry in archive.infolist():
                assert entry.date_time == (1980, 1, 1, 0, 0, 0)
                assert stat.S_ISREG(entry.external_attr >> 16)
                payload = archive.read(entry)
                assert hashlib.sha256(payload).hexdigest() == artifact["files"][entry.filename]["sha256"]
                assert b"PRIVATE_ENV_CANARY" not in payload and b"PRIVATE_MODULE_CANARY" not in payload
        assert snapshot(release) == before


    def test_cli_stdout_is_compact_and_points_to_complete_local_record(self):
        builder = load_script("build_api_service_release")
        output = self.tmp_path / "cli-service.zip"
        stdout = io.StringIO()
        with redirect_stdout(stdout):
            status = builder.main([
                "--service-dir", str(self.release / "service"),
                "--bundle-dir", str(self.release / "bundle"),
                "--public-site-dir", str(self.release / "site"),
                "--output", str(output),
            ])
        assert status == 0
        summary = json.loads(stdout.getvalue())
        assert "files" not in summary
        record_path = output.with_suffix(".artifact.json")
        assert summary["artifact_record"] == str(record_path.resolve())
        record = json.loads(record_path.read_text())
        assert len(record["files"]) == record["file_count"]
        assert "artifact_record" not in record
        for key, value in record.items():
            if key != "files":
                assert summary[key] == value
        with zipfile.ZipFile(output) as archive:
            assert set(archive.namelist()) == set(record["files"])
            for name, descriptor in record["files"].items():
                assert hashlib.sha256(archive.read(name)).hexdigest() == descriptor["sha256"]

    def test_symlinks_rejected_and_inputs_preserved(self):
        release, tmp_path = self.release, self.tmp_path
        for location in ("bundle", "service", "site"):
            with self.subTest(location=location):
                bundle_root = json.loads((release / "bundle/manifest.json").read_text())
                release_id = bundle_root["current_release_id"]
                if location == "bundle":
                    target = release / "bundle" / bundle_root["releases"][0]["manifest_path"]
                elif location == "service":
                    target = release / "service/src/server.js"
                else:
                    target = release / f"site/data/releases/{release_id}/human/records.jsonl.gz"
                private = tmp_path / "private-canary"
                private.write_bytes(target.read_bytes())
                target.unlink()
                target.symlink_to(private)
                before = snapshot(release)
                with self.assertRaises((ValueError, OSError)):
                    package(release, tmp_path / "blocked.zip")
                assert snapshot(release) == before
                assert not (tmp_path / "blocked.zip").exists()
                target.unlink()
                target.write_bytes(private.read_bytes())

    def test_large_file_count_uses_zip64_and_reports_requirement(self):
        output = self.tmp_path / "many-files.zip"
        with patch.object(zipfile, "ZIP_FILECOUNT_LIMIT", 40):
            artifact = package(self.release, output)
        assert artifact["zip64_required"] is True
        with zipfile.ZipFile(output) as archive:
            assert len(archive.namelist()) == artifact["file_count"]
            assert archive.testzip() is None


    def test_unknown_data_canary_and_changed_data_fail_closed(self):
        release, tmp_path = self.release, self.tmp_path
        private = release / "bundle/private-report.json"
        private.write_text('{"private":"PRIVATE_REPORT_CANARY"}')
        before = snapshot(release)
        with self.assertRaises(ValueError):
            package(release, tmp_path / "private.zip")
        assert snapshot(release) == before
        private.unlink()
        root = json.loads((release / "bundle/manifest.json").read_text())
        path = release / "bundle" / root["releases"][0]["manifest_path"]
        path.write_bytes(path.read_bytes() + b" ")
        before = snapshot(release)
        with self.assertRaises(ValueError):
            package(release, tmp_path / "changed.zip")
        assert snapshot(release) == before


    def test_size_failure_and_dry_run_create_no_outputs(self):
        release, tmp_path = self.release, self.tmp_path
        for budget in ({"max_zip_bytes": 100}, {"max_expanded_bytes": 100}, {"max_zip_bytes": 100_000_001}):
            with self.subTest(budget=budget):
                output = tmp_path / "absent/output.zip"
                before = snapshot(release)
                with self.assertRaises(ValueError):
                    package(release, output, dry_run=True, **budget)
                assert not output.parent.exists()
                assert snapshot(release) == before


    def test_existing_output_and_unsafe_inputs_preserve_user_files(self):
        release, tmp_path = self.release, self.tmp_path
        output = tmp_path / "existing.zip"
        output.write_bytes(b"USER_ARTIFACT")
        with self.assertRaisesRegex(ValueError, "already exists"):
            package(release, output)
        assert output.read_bytes() == b"USER_ARTIFACT"
        with self.assertRaisesRegex(ValueError, "outside"):
            package(release, release / "bundle/output.zip")
        module = __import__("agdesign2.public_api_release", fromlist=["_relative"])
        for value in ("../private", "/absolute", "releases/../private", "releases//private", "releases/%2e%2e/private", "releases\\private", "C:/private", ".env"):
            with self.assertRaises(ValueError):
                module._relative(value)

    def test_oversized_root_manifest_is_rejected_before_zip_writes(self):
        path = self.release / "bundle/manifest.json"
        content = path.read_bytes().ljust(1024 * 1024 + 1, b" ")
        path.write_bytes(content)
        (path.parent / "manifest.sha256").write_text(hashlib.sha256(content).hexdigest() + "  manifest.json\n")
        before = snapshot(self.release)
        with self.assertRaises(ValueError):
            package(self.release, self.tmp_path / "root-too-large.zip")
        assert snapshot(self.release) == before and not (self.tmp_path / "root-too-large.zip").exists()

    def test_oversized_release_manifest_is_rejected_before_zip_writes(self):
        bundle = self.release / "bundle"
        root = json.loads((bundle / "manifest.json").read_text())
        path = bundle / root["releases"][0]["manifest_path"]
        replace_release_manifest(bundle, self.release / "site", root, path.read_bytes().ljust(64 * 1024 * 1024 + 1, b" "))
        before = snapshot(self.release)
        with self.assertRaises(ValueError):
            package(self.release, self.tmp_path / "release-too-large.zip")
        assert snapshot(self.release) == before and not (self.tmp_path / "release-too-large.zip").exists()

    def test_aggregate_retained_indexes_obey_runtime_budget(self):
        bundle = self.release / "aggregate-bundle"
        build_public_api_bundle(
            self.release / "human/batch_summary.json", self.release / "mouse/batch_summary.json", bundle,
            public_site_dir=self.release / "site", software_commit="4" * 40,
            snapshot_id="aggregate-budget-fixture", previous_bundle_dir=self.release / "bundle",
        )
        root = json.loads((bundle / "manifest.json").read_text())
        for pointer_index, pointer in enumerate(root["releases"]):
            manifest = json.loads((bundle / pointer["manifest_path"]).read_text())
            for species in ("human", "mouse"):
                name = manifest["species"][species]["index_path"]
                path = bundle / name
                content = path.read_bytes().ljust(17 * 1024 * 1024, b" ")
                path.write_bytes(content)
                manifest["files"][name].update(bytes=len(content), sha256=hashlib.sha256(content).hexdigest())
            content = (json.dumps(manifest, ensure_ascii=False, sort_keys=True, separators=(",", ":")) + "\n").encode()
            replace_release_manifest(bundle, self.release / "site", root, content, pointer_index)
        before = snapshot(self.release)
        with self.assertRaises(ValueError):
            build_api_service_release(self.release / "service", bundle, self.release / "site", self.tmp_path / "indexes-too-large.zip", repo_root=REPO)
        assert snapshot(self.release) == before and not (self.tmp_path / "indexes-too-large.zip").exists()


    def test_retention_validates_final_static_copies_and_separates_code_id(self):
        release, tmp_path = self.release, self.tmp_path
        current_root = json.loads((release / "bundle/manifest.json").read_text())
        old_id = current_root["current_release_id"]
        old_manifest = (release / "bundle" / current_root["releases"][0]["manifest_path"]).read_bytes()
        original = package(release, tmp_path / "one.zip", dry_run=True)
        server = release / "service/src/server.js"
        server_bytes = server.read_bytes()
        server.write_bytes(b"// a code-only change\n" + server_bytes)
        code_only = package(release, tmp_path / "code-only.zip", dry_run=True)
        assert code_only["current_release_id"] == old_id
        assert code_only["service_build_id"] != original["service_build_id"]
        server.write_bytes(server_bytes)
        next_bundle = release / "next-bundle"
        build_public_api_bundle(
            release / "human/batch_summary.json", release / "mouse/batch_summary.json", next_bundle,
            public_site_dir=release / "site", software_commit="1" * 40, snapshot_id="synthetic-next-release",
            previous_bundle_dir=release / "bundle",
        )
        artifact = build_api_service_release(release / "service", next_bundle, release / "site", tmp_path / "two.zip", repo_root=REPO)
        assert len(artifact["releases"]) == 2
        assert artifact["releases"][1]["release_id"] == old_id
        assert artifact["service_build_id"] == original["service_build_id"]
        assert artifact["current_release_id"] != old_id
        assert (next_bundle / f"releases/{old_id}/manifest.json").read_bytes() == old_manifest
        retained_static = release / "site" / f"data/releases/{old_id}/mouse/records.jsonl.gz"
        retained_bytes = retained_static.read_bytes()
        retained_static.unlink()
        with self.assertRaises(ValueError):
            build_api_service_release(release / "service", next_bundle, release / "site", tmp_path / "missing.zip", repo_root=REPO)
        retained_static.write_bytes(retained_bytes)


    def test_both_final_portal_metadata_releases_must_match(self):
        for name in ("portal_metadata.json", "mouse/portal_metadata.json"):
            with self.subTest(name=name):
                path = self.release / "site" / name
                original = path.read_bytes()
                metadata = json.loads(original)
                metadata["public_api"]["release_id"] = "f" * 64
                path.write_text(json.dumps(metadata))
                with self.assertRaisesRegex(ValueError, "portal metadata"):
                    package(self.release, self.tmp_path / "mismatch.zip")
                assert not (self.tmp_path / "mismatch.zip").exists()
                path.write_bytes(original)


    def test_audit_report_is_recorded_without_claiming_a_pass(self):
        release, tmp_path = self.release, self.tmp_path
        audit = tmp_path / "audit.json"
        audit.write_text(json.dumps({"metadata": {"vulnerabilities": {"high": 0, "critical": 0}}}))
        artifact = package(release, tmp_path / "audit.zip", dry_run=True, dependency_audit_report=audit)
        assert artifact["dependency_audit"]["status"] == "provided_unverified"
        assert artifact["dependency_audit"]["report"]["sha256"] == hashlib.sha256(audit.read_bytes()).hexdigest()

    def test_retired_tombstones_package_without_retired_payload_files(self):
        first = json.loads((self.release / "bundle/manifest.json").read_text())
        previous = self.release / "bundle"
        for generation in (2, 3):
            destination = self.release / f"generation-{generation}"
            build_public_api_bundle(
                self.release / "human/batch_summary.json", self.release / "mouse/batch_summary.json", destination,
                public_site_dir=self.release / "site", software_commit=str(generation) * 40,
                snapshot_id=f"synthetic-generation-{generation}", previous_bundle_dir=previous,
            )
            previous = destination
        root = json.loads((previous / "manifest.json").read_text())
        tombstone = {"release_id": first["current_release_id"], "retired_by_release_id": root["current_release_id"], "archive_url": None}
        assert root["retired_releases"] == [tombstone]
        artifact = build_api_service_release(self.release / "service", previous, self.release / "site", self.tmp_path / "retired.zip", repo_root=REPO)
        assert artifact["retired_releases"] == [tombstone]
        with zipfile.ZipFile(self.tmp_path / "retired.zip") as archive:
            assert not any(name.startswith(f"data/releases/{first['current_release_id']}/") for name in archive.namelist())
        root["retired_releases"][0]["archive_url"] = "https://PRIVATE_CANARY@archive.example.org/file.zip"
        content = (json.dumps(root, ensure_ascii=False, sort_keys=True, separators=(",", ":")) + "\n").encode()
        (previous / "manifest.json").write_bytes(content)
        (previous / "manifest.sha256").write_text(hashlib.sha256(content).hexdigest() + "  manifest.json\n")
        before = snapshot(self.release)
        with self.assertRaises(ValueError):
            build_api_service_release(self.release / "service", previous, self.release / "site", self.tmp_path / "unsafe-retired.zip", repo_root=REPO)
        assert snapshot(self.release) == before and not (self.tmp_path / "unsafe-retired.zip").exists()


    def test_credential_bearing_lockfile_is_rejected(self):
        release, tmp_path = self.release, self.tmp_path
        path = release / "service/package-lock.json"
        lock = json.loads(path.read_text())
        entry = next(value for name, value in lock["packages"].items() if name)
        entry["resolved"] = "https://PRIVATE_TOKEN_CANARY@registry.npmjs.org/package.tgz"
        path.write_text(json.dumps(lock))
        before = snapshot(release)
        with self.assertRaisesRegex(ValueError, "public npm registry"):
            package(release, tmp_path / "secret.zip")
        assert snapshot(release) == before and not (tmp_path / "secret.zip").exists()


    def test_discovery_activation_updates_both_and_preserves_unrelated_metadata(self):
        release, tmp_path = self.release, self.tmp_path
        smoke = load_script("check_api_service")
        release_id = "a" * 64
        site = tmp_path / "site"
        for name in ("portal_metadata.json", "mouse/portal_metadata.json"):
            path = site / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps({"other": "keep", "public_api": {"schema_version": 1, "release_id": release_id, "state": "packaged", "manifest_url": "keep-url"}}))
        checked = {"base_url": "https://api.example.org", "release_id": release_id, "species": ["human", "mouse"], "http": "passed", "mcp_wire": "passed"}
        smoke.activate_local_discovery(site, checked)
        for name in ("portal_metadata.json", "mouse/portal_metadata.json"):
            value = json.loads((site / name).read_text())
            assert value["other"] == "keep" and value["public_api"]["manifest_url"] == "keep-url"
            assert value["public_api"]["state"] == "active"
            assert value["public_api"]["service_url"] == "https://api.example.org"
            assert value["public_api"]["mcp_url"] == "https://api.example.org/mcp"
        original = snapshot(site)
        with self.assertRaises(ValueError):
            smoke.activate_local_discovery(site, {**checked, "release_id": "b" * 64})
        assert snapshot(site) == original


    def test_failed_smoke_preserves_local_discovery(self):
        release, tmp_path = self.release, self.tmp_path
        smoke = load_script("check_api_service")
        site = tmp_path / "site"
        site.mkdir()
        (site / "portal_metadata.json").write_text("USER_METADATA_CANARY")
        before = snapshot(site)
        def fail(*args, **kwargs):
            raise ValueError("synthetic release mismatch")
        with patch.object(smoke, "check_api_service", fail):
            assert smoke.main(["--base-url", "https://api.example.org", "--release-id", "a" * 64, "--activate-discovery-local-site", str(site)]) == 2
        assert snapshot(site) == before

    def test_interrupted_discovery_write_restores_both_originals(self):
        smoke = load_script("check_api_service")
        site = self.tmp_path / "site"
        for name in ("portal_metadata.json", "mouse/portal_metadata.json"):
            path = site / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps({"public_api": {"schema_version": 1, "release_id": "a" * 64, "state": "packaged"}}))
        before = snapshot(site)
        checked = {"base_url": "https://api.example.org", "release_id": "a" * 64, "species": ["human", "mouse"], "http": "passed", "mcp_wire": "passed"}
        real_replace = smoke.os.replace
        calls = 0
        def interrupt_second(source, destination):
            nonlocal calls
            calls += 1
            if calls == 2:
                raise OSError("synthetic interrupted write")
            return real_replace(source, destination)
        with patch.object(smoke.os, "replace", interrupt_second):
            with self.assertRaises(OSError):
                smoke.activate_local_discovery(site, checked)
        assert snapshot(site) == before
        assert not list(site.rglob(".api-discovery-*"))


    def test_smoke_rejects_urls_with_credentials_and_remote_cleartext(self):
        release, tmp_path = self.release, self.tmp_path
        smoke = load_script("check_api_service")
        for url in ("https://token@example.org", "https://example.org/?key=secret", "http://example.org", "https://example.org/mcp"):
            with self.assertRaises(ValueError):
                smoke.check_api_service(url, "a" * 64)

    def test_maintenance_metadata_hides_current_and_retained_discovery(self):
        smoke = load_script("check_api_service")
        pinned = {"release_id": "a" * 64, "discovery_state": "publication_in_progress"}
        assert smoke._discovery(pinned, "b" * 64, True) == []
        for leaked in ({"current_release_id": "b" * 64}, {"releases": []}, {"current_release_id": None}, {"releases": None}):
            with self.subTest(leaked=leaked), self.assertRaises(ValueError):
                smoke._discovery({**pinned, **leaked}, "b" * 64, True)
        active = {"current_release_id": "a" * 64, "releases": [{"release_id": "a" * 64}, {"release_id": "b" * 64}]}
        assert smoke._discovery(active, "a" * 64, False) == ["a" * 64, "b" * 64]

    def test_leaked_maintenance_discovery_cannot_activate_local_files(self):
        smoke = load_script("check_api_service")
        site = self.tmp_path / "site"
        for name in ("portal_metadata.json", "mouse/portal_metadata.json"):
            path = site / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("USER_METADATA_CANARY")
        before = snapshot(site)
        def response(opener, base, path, **kwargs):
            headers = {"Cache-Control": "no-store", "Retry-After": "60"}
            if path == "/healthz":
                return 200, {"status": "ok"}, headers
            if path == "/readyz":
                return 200, {"status": "ready"}, headers
            if path == "/v1/releases/current":
                return 503, {"error": {"code": "publication_in_progress"}}, headers
            return 200, {"release_id": "a" * 64, "schema_version": 1, "discovery_state": "publication_in_progress", "current_release_id": "a" * 64}, headers
        with patch.object(smoke, "_request", response):
            assert smoke.main(["--base-url", "https://api.example.org", "--release-id", "a" * 64, "--expect-maintenance", "--activate-discovery-local-site", str(site)]) == 2
        assert snapshot(site) == before


if __name__ == "__main__":
    unittest.main()
