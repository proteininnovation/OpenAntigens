"""Offline query compilation preserves source bytes and sealed complete coverage."""
from __future__ import annotations

import copy
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

from agdesign2.public_api_data import PublicDataError, archive_url_policy, build_public_api_bundle, canonical_bytes
from agdesign2.public_api_query import _build_query_runtime, build_query_runtime, validate_query_runtime
from agdesign2.public_api_transport import build_runtime_data

REPO = Path(__file__).resolve().parents[1]
SCHEMA = REPO / "services/public-api/src/query-schema.json"


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def digest(content):
    return hashlib.sha256(content).hexdigest()


def write_json(path, value):
    content = canonical_bytes(value) + b"\n"
    path.write_bytes(content)
    return {"bytes": len(content), "sha256": digest(content)}


def snapshot(root):
    return {p.relative_to(root).as_posix(): digest(p.read_bytes()) for p in root.rglob("*") if p.is_file()}


class QueryCompilerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory(prefix="query-compiler-tests-")
        cls.root = Path(cls.temporary.name).resolve()
        cls.frozen = cls.root / "frozen"
        fixture = load_module("query_fixture", REPO / "tests/build_public_api_fixture.py")
        fixture.build(cls.frozen)
        cls.bundle = cls.frozen / "bundle"
        cls.before = snapshot(cls.bundle)
        cls.runtime, cls.static = cls.root / "runtime", cls.root / "static"
        cls.result = build_query_runtime(cls.bundle, cls.runtime, cls.static, schema_path=SCHEMA)

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def test_exact_transport_summaries_schema_and_source_preserved(self):
        self.assertEqual(self.before, snapshot(self.bundle))
        inventory = validate_query_runtime(self.runtime, schema_path=SCHEMA, seal=self.result["seal"])
        self.assertEqual(inventory, self.result["static_assets"])
        self.assertEqual(set(inventory), set(snapshot(self.static)))
        source_transport = self.frozen / "runtime"
        self.assertEqual((self.runtime / "transport-index.json").read_bytes(), (source_transport / "transport-index.json").read_bytes())
        transport = json.loads((self.runtime / "transport-index.json").read_bytes())
        for target, parts in transport["targets"].items():
            _, rid, species, _, filename = target.split("/")
            key = filename[:-5]
            for number, part in enumerate(parts):
                asset = f"api-runtime/{rid}/{species}/{key}-{number}-{part['sha256']}.oap"
                self.assertEqual((self.static / asset).read_bytes(), (source_transport / part["path"]).read_bytes())
        query = json.loads((self.runtime / "query-root.json").read_bytes())
        self.assertEqual(query["source_root_manifest_sha256"], self.before["manifest.json"])
        for pointer in query["releases"]:
            release = json.loads((self.runtime / pointer["path"]).read_bytes())
            self.assertNotIn("files", release)
            self.assertFalse((self.runtime / f"releases/{pointer['release_id']}/manifest.json").exists())
            for sm in release["species"].values():
                self.assertNotIn("target_paths", sm)
                for field in ("index", "summary_lines"):
                    path = sm[field]["path"]
                    self.assertEqual((self.runtime / path).read_bytes(), (self.bundle / path).read_bytes())
            schema = release["public_schema"]["path"]
            self.assertEqual((self.runtime / schema).read_bytes(), (self.bundle / schema).read_bytes())
        self.assertFalse((self.runtime / "packs").exists())
        again = build_query_runtime(self.bundle, self.root / "again-runtime", self.root / "again-static", schema_path=SCHEMA)
        self.assertEqual(again, self.result)
        self.assertEqual(snapshot(self.root / "again-runtime"), snapshot(self.runtime))

    def changed_runtime(self, name):
        root = self.root / name
        shutil.copytree(self.runtime, root)
        return root, copy.deepcopy(self.result["seal"])

    def rehash_release(self, root, seal, change):
        query_path = root / "query-root.json"
        query = json.loads(query_path.read_bytes())
        pointer = query["releases"][0]
        release_path = root / pointer["path"]
        release = json.loads(release_path.read_bytes())
        change(release)
        pointer.update(write_json(release_path, release))
        seal["query_root_sha256"] = write_json(query_path, query)["sha256"]

    def test_offset_gap_hash_and_record_order_fail_after_resealing(self):
        for kind in ("gap", "hash", "order"):
            root, seal = self.changed_runtime("offset-" + kind)
            def change(release):
                sm = release["species"]["human"]
                metadata = sm["summary_lines"]["offsets"]
                path = root / metadata["path"]
                offsets = json.loads(path.read_bytes())
                if kind == "gap":
                    offsets["rows"][0][0] = 1
                elif kind == "hash":
                    offsets["rows"][0][2] = "0" * 64
                else:
                    offsets["index_sha256"] = "0" * 64
                metadata.update(write_json(path, offsets))
            self.rehash_release(root, seal, change)
            with self.subTest(kind=kind), self.assertRaises(PublicDataError):
                validate_query_runtime(root, schema_path=SCHEMA, seal=seal)

    def test_missing_target_and_foreign_part_fail_after_resealing(self):
        for kind in ("missing", "foreign"):
            root, seal = self.changed_runtime("routing-" + kind)
            query_path = root / "query-root.json"
            query = json.loads(query_path.read_bytes())
            route_path = root / query["transport_index"]["path"]
            routing = json.loads(route_path.read_bytes())
            first = next(iter(routing["targets"]))
            if kind == "missing":
                del routing["targets"][first]
            else:
                routing["targets"][first][0]["path"] = "packs/../private.sql"
            query["transport_index"].update(write_json(route_path, routing))
            seal["query_root_sha256"] = write_json(query_path, query)["sha256"]
            with self.subTest(kind=kind), self.assertRaises(PublicDataError):
                validate_query_runtime(root, schema_path=SCHEMA, seal=seal)

    def test_shared_closed_schema_origin_duplicates_and_inventory(self):
        for kind in ("unknown", "origin", "duplicate", "extra", "directory", "symlink"):
            root, seal = self.changed_runtime("closed-" + kind)
            query_path = root / "query-root.json"
            query = json.loads(query_path.read_bytes())
            if kind == "unknown":
                query["extra"] = "not allowed"
                seal["query_root_sha256"] = write_json(query_path, query)["sha256"]
            elif kind == "origin":
                query["static_origin"] = "https://unapproved.invalid"
                seal["static_origin"] = query["static_origin"]
                seal["query_root_sha256"] = write_json(query_path, query)["sha256"]
            elif kind == "duplicate":
                raw = query_path.read_bytes().replace(b'"format":', b'"format":"PRIVATE_CANARY","format":', 1)
                query_path.write_bytes(raw); seal["query_root_sha256"] = digest(raw)
            elif kind == "extra":
                (root / "private.sql").write_text("PRIVATE_CANARY")
            elif kind == "directory":
                (root / "unexpected-directory").mkdir()
            else:
                (root / "alias").symlink_to(self.runtime / "query-root.json")
            with self.subTest(kind=kind), self.assertRaises(PublicDataError) as caught:
                validate_query_runtime(root, schema_path=SCHEMA, seal=seal)
            self.assertNotIn("PRIVATE_CANARY", str(caught.exception))

    def test_shared_query_schema_uses_exact_public_text_policy(self):
        schema = json.loads(SCHEMA.read_bytes())
        self.assertEqual(schema["x-public-text-policy"], archive_url_policy())
        altered = self.root / "altered-query-schema.json"
        schema["x-public-text-policy"] = {}
        write_json(altered, schema)
        with self.assertRaises(PublicDataError):
            validate_query_runtime(self.runtime, schema_path=altered, seal=self.result["seal"])

    def test_validated_packager_path_rechecks_retirement_text(self):
        source = self.root / "retirement-source"
        shutil.copytree(self.bundle, source)
        path = source / "manifest.json"
        root = json.loads(path.read_bytes())
        root["retired_releases"] = [{"release_id": "f" * 64,
                                     "retired_by_release_id": root["current_release_id"],
                                     "archive_url": "https://archive.org/PRIVATE_CANARY_RETIREMENT"}]
        checksum = write_json(path, root)["sha256"]
        (source / "manifest.sha256").write_text(checksum + "  manifest.json\n")
        runtime, static = self.root / "retirement-output", self.root / "retirement-assets"
        with self.assertRaises(PublicDataError) as caught:
            _build_query_runtime(source, runtime, static, schema_path=SCHEMA, validated_root=root)
        self.assertFalse(runtime.exists())
        self.assertFalse(static.exists())
        self.assertNotIn("PRIVATE_CANARY", str(caught.exception))

    def test_retirement_archive_policy_after_query_root_resealing(self):
        safe = [None, "https://archive.org/openantigens/release.tar.gz", "https://archive.org:443/openantigens/release.tar.gz", "https://doi.org/10.5281/zenodo.123"]
        unsafe = [
            "http://archive.org/release.tar.gz", "//archive.org/release.tar.gz",
            "https://user:password@archive.org/release.tar.gz",
            "https://user@archive.org/release.tar.gz",
            "https://@archive.org/release.tar.gz",
            "https://:@archive.org/release.tar.gz",
            "https://archive.org/PRIVATE_CANARY_RETIREMENT",
            "https://archive.org/private/archive.tar.gz",
            "https://archive.org/PRIVATE%5FCANARY_RETIREMENT",
            "https://archive.org/PRIVATE%255FCANARY_RETIREMENT",
            "https://archive.org/%2Fprivate%2Farchive.tar.gz",
            "https://archive.org/%0Aarchive.tar.gz",
            "https://archive.org/%C2%A0archive.tar.gz",
            "https://archive.org/%E2%80%A8archive.tar.gz",
            "https://archive.org/archive\n.tar.gz",
            "https://archive.org/releases/archive.tar.gz?token=secret",
            "https://archive.org/releases/archive.tar.gz#secret",
            "https://archive.org:8443/releases/archive.tar.gz",
            "https://archive.org:invalid/releases/archive.tar.gz",
            "https://127.0.0.1/releases/archive.tar.gz",
            "https://[::1]/releases/archive.tar.gz",
            "https://localhost/releases/archive.tar.gz",
            "https://archive.local/releases/archive.tar.gz",
            "https://archive.invalid/releases/archive.tar.gz",
            "https://archive.org/PRIVATE%252525255FCANARY_RETIREMENT",
        ]
        for number, url in enumerate(safe + unsafe):
            root, seal = self.changed_runtime(f"retirement-{number}")
            path = root / "query-root.json"
            query = json.loads(path.read_bytes())
            query["retired_releases"] = [{"release_id": "f" * 64,
                                          "retired_by_release_id": query["current_release_id"],
                                          "archive_url": url}]
            seal["query_root_sha256"] = write_json(path, query)["sha256"]
            if url in safe:
                with self.subTest(case=number):
                    self.assertEqual(validate_query_runtime(root, schema_path=SCHEMA, seal=seal), self.result["static_assets"])
            else:
                with self.subTest(case=number), self.assertRaises(PublicDataError) as caught:
                    validate_query_runtime(root, schema_path=SCHEMA, seal=seal)
                self.assertNotIn("PRIVATE_CANARY", str(caught.exception))

    def test_aggregate_budget_fails_whole_build_without_outputs(self):
        runtime, static = self.root / "budget-runtime", self.root / "budget-static"
        with patch("agdesign2.public_api_query.TARGET_COMPRESSED_LIMIT", 1), self.assertRaises(PublicDataError):
            build_query_runtime(self.bundle, runtime, static, schema_path=SCHEMA)
        self.assertFalse(runtime.exists()); self.assertFalse(static.exists())
        self.assertEqual(self.before, snapshot(self.bundle))

    def test_one_mebibyte_summary_bound_is_enforced(self):
        runtime, static = self.root / "record-runtime", self.root / "record-static"
        with patch("agdesign2.public_api_query.RECORD_LIMIT", 1), self.assertRaises(PublicDataError):
            build_query_runtime(self.bundle, runtime, static, schema_path=SCHEMA)
        self.assertFalse(runtime.exists()); self.assertFalse(static.exists())

    def test_two_part_repacking_preserves_ids_and_changes_hash_qualified_assets(self):
        routing = json.loads((self.frozen / "runtime/transport-index.json").read_bytes())
        maximum = max(p["decoded_bytes"] for parts in routing["targets"].values() for p in parts)
        limit = (maximum * 3 + 4) // 5
        def split(source, output, *, validate):
            return build_runtime_data(source, output, validate=validate, part_limit=limit)
        runtime, static = self.root / "split-runtime", self.root / "split-static"
        with patch("agdesign2.public_api_query.build_runtime_data", side_effect=split):
            result = build_query_runtime(self.bundle, runtime, static, schema_path=SCHEMA)
        split_index = json.loads((runtime / "transport-index.json").read_bytes())
        self.assertEqual(max(map(len, split_index["targets"].values())), 2)
        self.assertEqual(set(split_index["targets"]), set(routing["targets"]))
        self.assertEqual(result["seal"]["source_root_manifest_sha256"], self.result["seal"]["source_root_manifest_sha256"])
        self.assertNotEqual(result["seal"]["query_root_sha256"], self.result["seal"]["query_root_sha256"])
        self.assertNotEqual(set(result["static_assets"]), set(self.result["static_assets"]))
        for path in set(result["static_assets"]) & set(self.result["static_assets"]):
            self.assertEqual((static / path).read_bytes(), (self.static / path).read_bytes())

    def test_current_previous_and_retirement_history_are_preserved(self):
        previous = self.bundle
        for number in (2, 3):
            destination = self.root / f"source-release-{number}"
            build_public_api_bundle(self.frozen / "human/batch_summary.json", self.frozen / "mouse/batch_summary.json", destination,
                                    public_site_dir=self.frozen / "site", software_commit="0" * 40,
                                    snapshot_id=f"query-release-{number}", previous_bundle_dir=previous)
            previous = destination
        result = build_query_runtime(previous, self.root / "retained-runtime", self.root / "retained-static", schema_path=SCHEMA)
        query = json.loads((self.root / "retained-runtime/query-root.json").read_bytes())
        original = json.loads((previous / "manifest.json").read_bytes())
        self.assertEqual(query["current_release_id"], original["current_release_id"])
        self.assertEqual(query["retired_releases"], original["retired_releases"])
        self.assertEqual(len(query["releases"]), 2)
        self.assertEqual(result["metrics"]["target_count"], 4)

    def test_unresolved_target_remains_in_local_summary_and_routing(self):
        helper = load_module("query_source_fixture", REPO / "tests/test_public_api_data.py")
        root = self.root / "unresolved"
        human, mouse, site = helper.source_fixture(root)
        rows = json.loads(human.read_bytes()); rows.append({"query": "UNRESOLVED_HUMAN", "status": "error"}); write_json(human, rows)
        path = site / "downloads/agdesign2_portal_index.json"
        rows = json.loads(path.read_bytes()); rows.append({"entry_name": "UNRESOLVED_HUMAN", "status": "error"}); write_json(path, rows)
        source = root / "bundle"
        build_public_api_bundle(human, mouse, source, public_site_dir=site, software_commit="0" * 40, snapshot_id="query-unresolved")
        result = build_query_runtime(source, root / "remote-runtime", root / "remote-static", schema_path=SCHEMA)
        self.assertEqual(result["metrics"]["target_count"], 3)
        query = json.loads((root / "remote-runtime/query-root.json").read_bytes())
        release = json.loads((root / "remote-runtime" / query["releases"][0]["path"]).read_bytes())
        sm = release["species"]["human"]
        self.assertEqual(sm["counts"]["targets"], 2)
        self.assertEqual(sm["counts"]["resolved_targets"], 1)
        summaries = [json.loads(line) for line in (root / "remote-runtime" / sm["summary_lines"]["path"]).read_bytes().splitlines()]
        self.assertEqual(sum(row["target_id"] is None for row in summaries), 1)

    def test_public_entrypoint_validates_source_and_preserves_existing_outputs(self):
        source = self.root / "invalid-source"; shutil.copytree(self.bundle, source)
        (source / "private.sql").write_text("PRIVATE_CANARY")
        with self.assertRaises(PublicDataError):
            build_query_runtime(source, self.root / "invalid-runtime", self.root / "invalid-static", schema_path=SCHEMA)
        before = snapshot(self.runtime)
        with self.assertRaises(PublicDataError):
            build_query_runtime(self.bundle, self.runtime, self.root / "unused-static", schema_path=SCHEMA)
        self.assertEqual(before, snapshot(self.runtime))


if __name__ == "__main__":
    unittest.main()
