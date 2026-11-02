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
from email.message import Message
from types import SimpleNamespace

from agdesign2.public_api_data import build_public_api_bundle
from agdesign2.public_api_release import SOURCE_FILES, build_api_service_release, repack_api_service_release

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
        after = snapshot(release)
        assert all(after[name] == value for name, value in before.items())
        assert all(name.startswith("site/api-runtime/") for name in after.keys() - before.keys())


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
            if key not in {"files", "query_static_assets"}:
                assert summary[key] == value
        with zipfile.ZipFile(output) as archive:
            assert set(archive.namelist()) == set(record["files"])
            for name, descriptor in record["files"].items():
                assert hashlib.sha256(archive.read(name)).hexdigest() == descriptor["sha256"]

    def test_code_only_repack_reuses_exact_data_and_replaces_only_current_code(self):
        from agdesign2 import public_api_data, public_api_query
        original = self.tmp_path / "original.zip"
        receipt = package(self.release, original)
        (self.release / "service/src/server.js").write_text("// reviewed runtime-only patch\n" + (self.release / "service/src/server.js").read_text())
        before = snapshot(self.release)
        output = self.tmp_path / "code-update.zip"
        with patch.object(public_api_data, "validate_public_api_bundle", side_effect=AssertionError("scientific validation must be inherited")), patch.object(public_api_query, "_build_query_runtime", side_effect=AssertionError("data must not be recompiled")):
            dry = repack_api_service_release(original.with_suffix(".artifact.json"), receipt["sha256"], self.release / "service", self.release / "site", output, repo_root=REPO, dry_run=True)
            assert not output.exists() and snapshot(self.release) == before
            updated = repack_api_service_release(original.with_suffix(".artifact.json"), receipt["sha256"], self.release / "service", self.release / "site", output, repo_root=REPO)
        assert updated["sha256"] == dry["sha256"] and updated["service_build_id"] != receipt["service_build_id"]
        assert updated["data_transport"] == receipt["data_transport"] and updated["current_release_id"] == receipt["current_release_id"]
        assert updated["validation_lineage"]["parent_zip_sha256"] == receipt["sha256"]
        assert "not_rerun" in updated["validation_lineage"]["scientific_validation"]
        with zipfile.ZipFile(original) as old, zipfile.ZipFile(output) as new:
            old_names = {name for name in old.namelist() if name.startswith("data/")}
            assert old_names == {name for name in new.namelist() if name.startswith("data/")}
            for name in old_names:
                assert old.read(name) == new.read(name)
            assert new.read("src/server.js").startswith(b"// reviewed runtime-only patch")
        assert snapshot(self.release) == before

    def test_code_only_repack_requires_independent_hash_and_all_static_bytes(self):
        original = self.tmp_path / "original.zip"
        receipt = package(self.release, original)
        output = self.tmp_path / "repacked.zip"
        before = snapshot(self.release)
        with self.assertRaises(ValueError):
            repack_api_service_release(original.with_suffix(".artifact.json"), "f" * 64, self.release / "service", self.release / "site", output, repo_root=REPO)
        assert not output.exists() and snapshot(self.release) == before
        missing = self.release / "site" / next(iter(receipt["query_static_assets"]))
        saved = missing.read_bytes()
        missing.unlink()
        try:
            with self.assertRaises(ValueError):
                repack_api_service_release(original.with_suffix(".artifact.json"), receipt["sha256"], self.release / "service", self.release / "site", output, repo_root=REPO)
            assert not output.exists()
        finally:
            missing.write_bytes(saved)

    def test_code_only_repack_uses_current_query_policy_without_rewriting_data(self):
        original = self.tmp_path / "original.zip"
        receipt = package(self.release, original)
        from agdesign2 import public_api_query
        before = snapshot(self.release)
        with patch.object(public_api_query, "validate_query_runtime", side_effect=ValueError("current retirement policy rejects this data")):
            with self.assertRaises(ValueError):
                repack_api_service_release(original.with_suffix(".artifact.json"), receipt["sha256"], self.release / "service", self.release / "site", self.tmp_path / "no-fallback.zip", repo_root=REPO)
        assert snapshot(self.release) == before and not (self.tmp_path / "no-fallback.zip").exists()

    def test_query_catalog_is_sealed_and_static_assets_preserve_logical_bytes(self):
        from agdesign2.public_api_transport import decode_part
        output = self.tmp_path / "query.zip"
        artifact = package(self.release, output)
        root = json.loads((self.release / "bundle/manifest.json").read_text())
        manifest = json.loads((self.release / "bundle" / root["releases"][0]["manifest_path"]).read_text())
        with zipfile.ZipFile(output) as archive:
            build = json.loads(archive.read("service-build.json"))
            seal = build["data_transport"]
            assert seal == artifact["data_transport"] and seal["format"] == "openantigens-own-static-query-v1"
            assert seal["source_root_manifest_sha256"] == hashlib.sha256((self.release / "bundle/manifest.json").read_bytes()).hexdigest()
            assert seal["query_root_sha256"] == hashlib.sha256(archive.read("data/query-root.json")).hexdigest()
            assert "data/manifest.json" not in archive.namelist()
            assert "data/" + root["releases"][0]["manifest_path"] not in archive.namelist()
            transport = json.loads(archive.read("data/transport-index.json"))
            for target, parts in transport["targets"].items():
                elements = Path(target).parts
                for number, descriptor in enumerate(parts):
                    asset = f"api-runtime/{elements[1]}/{elements[2]}/{Path(target).stem}-{number}-{descriptor['sha256']}.oap"
                    assert asset in artifact["query_static_assets"]
                    members = decode_part((self.release / "site" / asset).read_bytes(), descriptor, manifest["files"])
                    for name, content in members.items():
                        assert content == (self.release / "bundle" / name).read_bytes()
            excluded = {"index.tsv", "constructs.jsonl", "constructs.tsv", "constructs.fasta", "reference.fasta", "records.jsonl.gz"}
            assert not any("/downloads/" in name and Path(name).name in excluded for name in archive.namelist())
            assert not any(name.endswith(".oap") or name.startswith("data/packs/") for name in archive.namelist())
            data_entries = [info for info in archive.infolist() if info.filename.startswith("data/")]
            assert artifact["physical_data_file_count"] == len(data_entries)
            assert artifact["physical_data_bytes"] == sum(info.file_size for info in data_entries)
        assert artifact["query_static_asset_file_count"] == len(artifact["query_static_assets"])
        assert artifact["external_bulk_download_verification"]["runtime"] == "metadata_only_external_static_files"

    def test_changed_query_outputs_are_rejected_before_zip(self):
        from agdesign2 import public_api_query
        real_build = public_api_query._build_query_runtime
        before = snapshot(self.release)
        for defect in ("runtime", "asset", "extra", "seal", "private seal field"):
            with self.subTest(defect=defect):
                def tamper(source, destination, assets, **kwargs):
                    result = real_build(source, destination, assets, **kwargs)
                    if defect == "runtime":
                        path = destination / "query-root.json"
                        path.write_bytes(path.read_bytes() + b"PRIVATE_CANARY")
                    elif defect == "asset":
                        path = assets / next(iter(result["static_assets"]))
                        path.write_bytes(path.read_bytes() + b"PRIVATE_CANARY")
                    elif defect == "extra":
                        (destination / "private.sqlite").write_bytes(b"PRIVATE_CANARY")
                    elif defect == "seal":
                        result["seal"]["source_root_manifest_sha256"] = "f" * 64
                    else:
                        result["seal"]["source_directory"] = "/Users/PRIVATE_CANARY"
                    return result
                with patch.object(public_api_query, "_build_query_runtime", tamper):
                    with self.assertRaises(ValueError):
                        package(self.release, self.tmp_path / "corrupt-query.zip")
                assert not (self.tmp_path / "corrupt-query.zip").exists()
                assert snapshot(self.release) == before

    def test_dry_run_stages_nothing_and_immutable_collisions_preserve_existing_bytes(self):
        output = self.tmp_path / "query.zip"
        before = snapshot(self.release)
        preview = package(self.release, output, dry_run=True)
        assert snapshot(self.release) == before and not (self.release / "site/api-runtime").exists()
        name = next(iter(preview["query_static_assets"]))
        target = self.release / "site" / name
        target.parent.mkdir(parents=True)
        target.write_bytes(b"USER_EXISTING_IMMUTABLE_BYTES")
        collision = snapshot(self.release)
        with self.assertRaises(ValueError):
            package(self.release, output)
        assert snapshot(self.release) == collision and not output.exists()

    def test_publisher_derives_assets_from_sealed_zip_not_hash_shaped_directory_files(self):
        from agdesign2.public_api_release import write_query_asset_list
        output = self.tmp_path / "query.zip"
        artifact = package(self.release, output)
        unknown = self.release / "site/api-runtime" / ("f" * 64) / "human" / (("e" * 64) + "-0-" + ("d" * 64) + ".oap")
        unknown.parent.mkdir(parents=True)
        unknown.write_bytes(b"PRIVATE_CANARY")
        file_list = self.tmp_path / "publish.list"
        count = write_query_asset_list(self.release / "site", output.with_suffix(".artifact.json"), file_list)
        names = file_list.read_bytes().rstrip(b"\0").decode().split("\0")
        assert count == len(artifact["query_static_assets"]) + 24
        assert {name for name in names if name.startswith("api-runtime/")} == set(artifact["query_static_assets"])
        assert len([name for name in names if name.startswith("data/releases/")]) == 24
        assert unknown.relative_to(self.release / "site").as_posix() not in names

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
        with patch.object(zipfile, "ZIP_FILECOUNT_LIMIT", 20):
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


    def test_packager_independently_rejects_manifest_listed_private_database(self):
        from agdesign2 import public_api_data
        bundle = self.release / "bundle"
        root = json.loads((bundle / "manifest.json").read_text())
        manifest = json.loads((bundle / root["releases"][0]["manifest_path"]).read_text())
        release_id = root["current_release_id"]
        name = f"releases/{release_id}/human/downloads/private.sqlite"
        content = b"PRIVATE_CANARY synthetic raw confidential input\n"
        (bundle / name).write_bytes(content)
        descriptor = {"sha256": hashlib.sha256(content).hexdigest(), "bytes": len(content), "media_type": "application/octet-stream", "encoding": None, "schema": None}
        manifest["files"][name] = descriptor
        manifest["species"]["human"]["downloads"].append({"name": "private.sqlite", "url": f"https://openantigens.org/data/releases/{release_id}/human/private.sqlite", "path": name, **descriptor, "row_count": None, "license": "synthetic test"})
        (self.release / f"site/data/releases/{release_id}/human/private.sqlite").write_bytes(content)
        replace_release_manifest(bundle, self.release / "site", root, json.dumps(manifest, sort_keys=True, separators=(",", ":")).encode() + b"\n")
        before = snapshot(self.release)
        with patch.object(public_api_data, "validate_public_api_bundle", return_value=root):
            with self.assertRaises(ValueError):
                package(self.release, self.tmp_path / "private-database.zip")
        assert snapshot(self.release) == before and not (self.tmp_path / "private-database.zip").exists()

    def test_packager_scans_untyped_canary_masquerading_as_approved_exports(self):
        from agdesign2 import public_api_data
        bundle = self.release / "bundle"
        for download_name in ("openantigens.bib", "index.tsv"):
            with self.subTest(download=download_name):
                root = json.loads((bundle / "manifest.json").read_text())
                manifest = json.loads((bundle / root["releases"][0]["manifest_path"]).read_text())
                item = next(row for row in manifest["species"]["human"]["downloads"] if row["name"] == download_name)
                path = bundle / item["path"]
                original = path.read_bytes()
                content = original + b"\nPRIVATE_CANARY synthetic untyped export\n"
                path.write_bytes(content)
                item.update(bytes=len(content), sha256=hashlib.sha256(content).hexdigest())
                manifest["files"][item["path"]].update(bytes=item["bytes"], sha256=item["sha256"])
                (self.release / f"site/data/releases/{root['current_release_id']}/human/{download_name}").write_bytes(content)
                replace_release_manifest(bundle, self.release / "site", root, json.dumps(manifest, sort_keys=True, separators=(",", ":")).encode() + b"\n")
                before = snapshot(self.release)
                with patch.object(public_api_data, "validate_public_api_bundle", return_value=root):
                    with self.assertRaises(ValueError):
                        package(self.release, self.tmp_path / "masked-canary.zip")
                assert snapshot(self.release) == before and not (self.tmp_path / "masked-canary.zip").exists()
                path.write_bytes(original)
                item.update(bytes=len(original), sha256=hashlib.sha256(original).hexdigest())
                manifest["files"][item["path"]].update(bytes=item["bytes"], sha256=item["sha256"])
                (self.release / f"site/data/releases/{root['current_release_id']}/human/{download_name}").write_bytes(original)
                replace_release_manifest(bundle, self.release / "site", root, json.dumps(manifest, sort_keys=True, separators=(",", ":")).encode() + b"\n")

    def test_packager_independently_rejects_unreachable_typed_storage(self):
        from agdesign2 import public_api_data
        bundle = self.release / "bundle"
        root = json.loads((bundle / "manifest.json").read_text())
        manifest = json.loads((bundle / root["releases"][0]["manifest_path"]).read_text())
        source = next(name for name, metadata in manifest["files"].items() if metadata["schema"] == "sequence_shard")
        name = f"releases/{root['current_release_id']}/human/sequences/{'f' * 64}.json"
        (bundle / name).write_bytes((bundle / source).read_bytes())
        manifest["files"][name] = dict(manifest["files"][source])
        replace_release_manifest(bundle, self.release / "site", root, json.dumps(manifest, sort_keys=True, separators=(",", ":")).encode() + b"\n")
        before = snapshot(self.release)
        with patch.object(public_api_data, "validate_public_api_bundle", return_value=root):
            with self.assertRaises(ValueError):
                package(self.release, self.tmp_path / "unreachable.zip")
        assert snapshot(self.release) == before and not (self.tmp_path / "unreachable.zip").exists()

    def test_packager_independently_rejects_hidden_duplicate_json_values(self):
        from agdesign2 import public_api_data
        for defect in ("sequence", "schema"):
            with self.subTest(defect=defect):
                bundle = self.release / "bundle"
                root = json.loads((bundle / "manifest.json").read_text())
                manifest = json.loads((bundle / root["releases"][0]["manifest_path"]).read_text())
                if defect == "sequence":
                    names = [next(name for name, metadata in manifest["files"].items() if metadata["schema"] == "sequence_shard")]
                    key = "schema_version"
                else:
                    names = [name for name, metadata in manifest["files"].items() if metadata["media_type"] == "application/schema+json"]
                    key = "$id"
                originals = {}
                for name in names:
                    original = (bundle / name).read_bytes()
                    originals[name] = original
                    content = ("{" + json.dumps(key) + ':"PRIVATE_CANARY",').encode() + original.lstrip()[1:]
                    (bundle / name).write_bytes(content)
                    manifest["files"][name].update(bytes=len(content), sha256=hashlib.sha256(content).hexdigest())
                    for species in ("human", "mouse"):
                        for item in manifest["species"][species]["downloads"]:
                            if item["path"] == name:
                                item.update(bytes=len(content), sha256=hashlib.sha256(content).hexdigest())
                                (self.release / f"site/data/releases/{root['current_release_id']}/{species}/{item['name']}").write_bytes(content)
                replace_release_manifest(bundle, self.release / "site", root, json.dumps(manifest, sort_keys=True, separators=(",", ":")).encode() + b"\n")
                before = snapshot(self.release)
                with patch.object(public_api_data, "validate_public_api_bundle", return_value=root):
                    with self.assertRaises(ValueError) as raised:
                        package(self.release, self.tmp_path / "duplicate-json.zip")
                assert "PRIVATE_CANARY" not in str(raised.exception)
                assert snapshot(self.release) == before and not (self.tmp_path / "duplicate-json.zip").exists()
                for name, original in originals.items():
                    (bundle / name).write_bytes(original)
                    metadata = {"bytes": len(original), "sha256": hashlib.sha256(original).hexdigest()}
                    manifest["files"][name].update(**metadata)
                    for species in ("human", "mouse"):
                        for item in manifest["species"][species]["downloads"]:
                            if item["path"] == name:
                                item.update(**metadata)
                                (self.release / f"site/data/releases/{root['current_release_id']}/{species}/{item['name']}").write_bytes(original)
                replace_release_manifest(bundle, self.release / "site", root, json.dumps(manifest, sort_keys=True, separators=(",", ":")).encode() + b"\n")

    def test_original_website_disease_tsv_stays_outside_app_bundle(self):
        original = self.release / "site/open_targets_disease_associations.tsv"
        original.write_text("ORIGINAL_WEBSITE_DISEASE_TSV_CANARY\n")
        before = snapshot(self.release)
        output = self.tmp_path / "projection-only.zip"
        artifact = package(self.release, output)
        assert not any("open_targets_disease_associations.tsv" in name for name in artifact["files"])
        with zipfile.ZipFile(output) as archive:
            assert not any(b"ORIGINAL_WEBSITE_DISEASE_TSV_CANARY" in archive.read(name) for name in archive.namelist())
        after = snapshot(self.release)
        assert all(after[name] == value for name, value in before.items())
        assert all(name.startswith("site/api-runtime/") for name in after.keys() - before.keys())

    def test_external_bulk_limits_are_distinct_from_runtime_and_tiny_references(self):
        from agdesign2.public_api_release import _public_file_limit
        prefix = "releases/" + "a" * 64 + "/human/downloads/"
        for name in ("index.jsonl", "index.tsv", "constructs.jsonl", "constructs.tsv", "constructs.fasta", "reference.fasta", "records.jsonl.gz"):
            with self.subTest(role=name):
                assert _public_file_limit(prefix + name, "application/x-ndjson") == 1024 ** 3
        for name in ("schema.json", "citation.json", "openantigens.bib", "openantigens.ris"):
            with self.subTest(role=name):
                assert _public_file_limit(prefix + name, "application/json") == 1024 ** 2
        with self.assertRaises(ValueError):
            _public_file_limit(prefix + "open_targets_disease_associations.tsv", "text/tab-separated-values")
        assert _public_file_limit("releases/a/human/targets/b.json", "application/json") == 8 * 1024 ** 2
        assert _public_file_limit("releases/a/human/sequences/b.json", "application/json") == 1024 ** 2
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
        replace_release_manifest(bundle, self.release / "site", root, path.read_bytes().ljust(96 * 1024 * 1024 + 1, b" "))
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
            assert smoke.main(["--base-url", "https://api.example.org", "--release-id", "a" * 64, "--human-query", "SYNTHETIC", "--mouse-query", "SYNTHETIC", "--activate-discovery-local-site", str(site)]) == 2
        assert snapshot(site) == before

    def test_duplicate_metadata_cannot_activate_or_expose_hidden_values(self):
        smoke = load_script("check_api_service")
        site = self.tmp_path / "site"
        public = {"schema_version": 1, "release_id": "a" * 64, "state": "packaged"}
        for name in ("portal_metadata.json", "mouse/portal_metadata.json"):
            path = site / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps({"public_api": public}))
        (site / "portal_metadata.json").write_text('{"public_api":{"hidden":"PRIVATE_CANARY"},"public_api":' + json.dumps(public) + '}')
        before = snapshot(site)
        checked = {"base_url": "https://api.example.org", "release_id": "a" * 64, "species": ["human", "mouse"], "http": "passed", "mcp_wire": "passed"}
        with self.assertRaises(ValueError) as raised:
            smoke.activate_local_discovery(site, checked)
        assert "PRIVATE_CANARY" not in str(raised.exception) and snapshot(site) == before
        with self.assertRaises(ValueError):
            smoke._json_loads('{"outer":{"score":"PRIVATE_CANARY","score":1}}')

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
                smoke.check_api_service(url, "a" * 64, human_query="SYNTHETIC", mouse_query="SYNTHETIC")

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
            manifest = {"url": f"https://openantigens.org/data/releases/{'a' * 64}/manifest.json", "sha256": "b" * 64}
            provenance = {"software_commit": "fixture", "snapshot_id": "fixture", "source_versions": {}, "source_dates": {}, "notes": [], "input_hashes": {"count": 10_001, "manifest_url": manifest["url"], "manifest_sha256": manifest["sha256"], "json_pointer": "/provenance/input_sha256"}}
            return 200, {"release_id": "a" * 64, "schema_version": 1, "manifest": manifest, "provenance": provenance, "discovery_state": "publication_in_progress", "current_release_id": "a" * 64}, headers
        with patch.object(smoke, "_request", response):
            assert smoke.main(["--base-url", "https://api.example.org", "--release-id", "a" * 64, "--human-query", "SYNTHETIC", "--mouse-query", "SYNTHETIC", "--expect-maintenance", "--activate-discovery-local-site", str(site)]) == 2
        assert snapshot(site) == before

    def test_smoke_provenance_uses_immutable_reference_without_raw_input_map(self):
        smoke = load_script("check_api_service")
        release_id = "a" * 64
        manifest = {"url": f"https://openantigens.org/data/releases/{release_id}/manifest.json", "sha256": "b" * 64}
        provenance = {"software_commit": "fixture", "snapshot_id": "fixture", "source_versions": {}, "source_dates": {}, "notes": [], "input_hashes": {"count": 10_001, "manifest_url": manifest["url"], "manifest_sha256": manifest["sha256"], "json_pointer": "/provenance/input_sha256"}}
        response = {"release_id": release_id, "schema_version": 1, "manifest": manifest, "provenance": provenance, "data": []}
        assert smoke._envelope(response, release_id) is response
        with self.assertRaises(ValueError):
            smoke._envelope({**response, "provenance": {**provenance, "input_sha256": {"private": "c" * 64}}}, release_id)
        with self.assertRaises(ValueError):
            smoke._envelope({**response, "provenance": {**provenance, "input_hashes": {**provenance["input_hashes"], "manifest_sha256": "c" * 64}}}, release_id)

    def test_smoke_paces_and_bounds_rate_limit_retries(self):
        smoke = load_script("check_api_service")
        class Response:
            def __init__(self, status, retry=None):
                self.status = status
                self.headers = Message()
                self.headers["Content-Type"] = "application/json"
                if retry:
                    self.headers["Retry-After"] = retry
            def __enter__(self):
                return self
            def __exit__(self, *args):
                return None
            def read(self, limit):
                return b"{}"
        cases = (
            ("pace", [Response(200)], 8, 0.0, False, 9, 3.0),
            ("retry", [Response(429, "2"), Response(200)], 0, 0.0, False, 2, 2.0),
            ("attempt cap", [Response(429, "1")] * 3, 0, 0.0, True, 3, 2.0),
            ("wait cap", [Response(429, "2")], 0, 29.0, True, 1, 29.0),
            ("read cap", [], 20, 0.0, True, 20, 0.0),
        )
        for name, replies, initial_reads, initial_wait, fails, reads, waits in cases:
            with self.subTest(case=name):
                clock = [0.0]
                responses = iter(replies)
                opener = SimpleNamespace(open=lambda *args, **kwargs: next(responses))
                budget = {"started_at": 0.0, "requests": initial_reads, "wait_seconds": initial_wait}
                def advance(seconds):
                    clock[0] += seconds
                with patch.object(smoke.time, "monotonic", lambda: clock[0]), patch.object(smoke.time, "sleep", advance):
                    if fails:
                        with self.assertRaises(ValueError):
                            smoke._request(opener, "https://api.example.org", "/readyz", budget=budget)
                    else:
                        assert smoke._request(opener, "https://api.example.org", "/readyz", budget=budget)[0] == 200
                assert budget["requests"] == reads and budget["wait_seconds"] == waits

    def test_smoke_references_exclude_bulk_download_metadata(self):
        smoke = load_script("check_api_service")
        root = json.loads((self.release / "bundle/manifest.json").read_text())
        manifest = json.loads((self.release / "bundle" / root["releases"][0]["manifest_path"]).read_text())
        metadata = {"species": {species: {"taxon_id": value["taxon_id"], "capabilities": value["capabilities"], "references": [item for item in value["downloads"] if item["name"] in smoke.REFERENCE_ROLES]} for species, value in manifest["species"].items()}}
        smoke._references(metadata, root["current_release_id"])
        metadata["species"]["human"]["references"].append(manifest["species"]["human"]["records_download"])
        with self.assertRaises(ValueError):
            smoke._references(metadata, root["current_release_id"])


if __name__ == "__main__":
    unittest.main()
