"""Physical packing preserves authoritative virtual bytes and rejects malformed input."""
from __future__ import annotations

import base64
import copy
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import tempfile
import unittest

from agdesign2.public_api_data import PublicDataError, canonical_bytes, build_public_api_bundle
from agdesign2.public_api_transport import (
    FORMAT, EXTERNAL_DOWNLOADS, STANDALONE_DOWNLOADS,
    build_runtime_data, decode_part,
)

REPO = Path(__file__).resolve().parents[1]


def snapshot(root):
    return {p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in root.rglob("*") if p.is_file()}


def encoded_part(members):
    raw = canonical_bytes({"format": FORMAT, "members": members}) + b"\n"
    compressed = gzip.compress(raw, compresslevel=9, mtime=0)
    descriptor = {"path": "packs/" + "a" * 64 + "/human/" + "b" * 64 + "-0.json.gz", "sha256": hashlib.sha256(compressed).hexdigest(),
                  "bytes": len(compressed), "decoded_bytes": len(raw), "members": len(members),
                  "min_path": members[0][0], "max_path": members[-1][0]}
    return compressed, descriptor


class PublicTransportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory(prefix="api-transport-tests-")
        cls.root = Path(cls.temporary.name).resolve()
        spec = importlib.util.spec_from_file_location("transport_fixture", REPO / "tests/build_public_api_fixture.py")
        helper = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(helper)
        helper.build(cls.root / "frozen")
        cls.bundle = cls.root / "frozen/bundle"
        cls.runtime = cls.root / "runtime"
        cls.before = snapshot(cls.bundle)
        cls.seal = build_runtime_data(cls.bundle, cls.runtime)
        cls.index = json.loads((cls.runtime / "transport-index.json").read_bytes())
        cls.source_files = {}
        root = json.loads((cls.bundle / "manifest.json").read_bytes())
        for pointer in root["releases"]:
            cls.source_files.update(json.loads((cls.bundle / pointer["manifest_path"]).read_bytes())["files"])

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def test_complete_roundtrip_source_unchanged_and_deterministic(self):
        decoded = {}
        physical = {"manifest.json", "manifest.sha256", "transport-index.json"}
        for target, parts in self.index["targets"].items():
            for number, part in enumerate(parts):
                self.assertTrue(part["path"].endswith(f"-{number}.json.gz"))
                physical.add(part["path"])
                members = decode_part((self.runtime / part["path"]).read_bytes(), part, self.source_files)
                self.assertFalse(decoded.keys() & members.keys())
                decoded.update(members)
            self.assertIn(target, decoded)
        for path, content in decoded.items():
            self.assertEqual(content, (self.bundle / path).read_bytes())
        for path in self.source_files:
            if "/downloads/" in path:
                name = path.rsplit("/", 1)[-1]
                self.assertEqual((self.runtime / path).exists(), name in STANDALONE_DOWNLOADS)
                self.assertIn(name, STANDALONE_DOWNLOADS | EXTERNAL_DOWNLOADS)
            if (self.runtime / path).exists():
                physical.add(path)
        physical.update(p for p in self.before if p.endswith("/manifest.json"))
        self.assertEqual(set(snapshot(self.runtime)), physical)
        self.assertEqual(len(decoded), self.seal["packed_logical_file_count"])
        self.assertEqual(len(physical), self.seal["physical_file_count"])
        self.assertEqual(self.before, snapshot(self.bundle))
        again = self.root / "again"
        seal = build_runtime_data(self.bundle, again)
        self.assertEqual(seal, self.seal)
        self.assertEqual(snapshot(again), snapshot(self.runtime))
        self.assertEqual(self.seal["root_manifest_sha256"], self.before["manifest.json"])
        self.assertEqual(self.seal["transport_index_sha256"], snapshot(self.runtime)["transport-index.json"])

    def test_oversized_target_splits_with_complete_nonoverlapping_ranges(self):
        max_single = 0
        for parts in self.index["targets"].values():
            for part in parts:
                payload = json.loads(gzip.decompress((self.runtime / part["path"]).read_bytes()))
                for member in payload["members"]:
                    max_single = max(max_single, len(canonical_bytes({"format": FORMAT, "members": [member]})) + 1)
        destination = self.root / "split"
        seal = build_runtime_data(self.bundle, destination, part_limit=max_single)
        index = json.loads((destination / "transport-index.json").read_bytes())
        self.assertGreater(sum(map(len, index["targets"].values())), len(index["targets"]))
        count = 0
        for parts in index["targets"].values():
            previous = None
            for part in parts:
                self.assertLessEqual(part["decoded_bytes"], max_single)
                if previous is not None:
                    self.assertLess(previous, part["min_path"])
                previous = part["max_path"]
                count += len(decode_part((destination / part["path"]).read_bytes(), part, self.source_files))
        self.assertEqual(count, seal["packed_logical_file_count"])
        with self.assertRaises(PublicDataError):
            build_runtime_data(self.bundle, self.root / "too-small", part_limit=32)

    def test_current_and_previous_releases_have_separate_containers(self):
        combined = self.root / "retained-source"
        frozen = self.root / "frozen"
        build_public_api_bundle(frozen / "human/batch_summary.json", frozen / "mouse/batch_summary.json",
                                combined, public_site_dir=frozen / "site", software_commit="0" * 40,
                                snapshot_id="transport-retention-check", previous_bundle_dir=self.bundle)
        destination = self.root / "retained-runtime"
        seal = build_runtime_data(combined, destination)
        index = json.loads((destination / "transport-index.json").read_bytes())
        self.assertEqual(len({p.split("/")[1] for p in index["targets"]}), 2)
        self.assertEqual(len(index["targets"]), 2 * len(self.index["targets"]))
        self.assertEqual(seal["packed_logical_file_count"], 2 * self.seal["packed_logical_file_count"])
        self.assertEqual(seal["external_download_count"], 2 * self.seal["external_download_count"])

    def test_output_must_be_fresh_and_separate(self):
        for destination in (self.bundle, self.bundle / "nested", self.runtime):
            with self.subTest(destination=destination), self.assertRaises(PublicDataError):
                build_runtime_data(self.bundle, destination, validate=False)

    def changed_source(self, name):
        destination = self.root / name
        shutil.copytree(self.bundle, destination)
        return destination

    def rehash_manifest(self, bundle, change):
        root = json.loads((bundle / "manifest.json").read_bytes())
        pointer = root["releases"][0]
        path = bundle / pointer["manifest_path"]
        manifest = json.loads(path.read_bytes())
        change(manifest)
        raw = canonical_bytes(manifest) + b"\n"
        path.write_bytes(raw)
        pointer["sha256"] = hashlib.sha256(raw).hexdigest()
        raw = canonical_bytes(root) + b"\n"
        (bundle / "manifest.json").write_bytes(raw)
        (bundle / "manifest.sha256").write_text(hashlib.sha256(raw).hexdigest() + "  manifest.json\n")

    def test_unlisted_private_file_rejected(self):
        source = self.changed_source("private-extra")
        (source / "private.sql").write_text("PRIVATE_CANARY")
        with self.assertRaises(PublicDataError):
            build_runtime_data(source, self.root / "private-output", validate=False)

    def test_rehashed_unknown_role_rejected(self):
        source = self.changed_source("private-role")
        def change(manifest):
            path = f"releases/{manifest['release_id']}/human/downloads/private.sql"
            raw = b"PRIVATE_CANARY"
            (source / path).write_bytes(raw)
            manifest["files"][path] = {"bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest(),
                                       "media_type": "application/octet-stream", "schema": None, "encoding": None}
        self.rehash_manifest(source, change)
        with self.assertRaises(PublicDataError):
            build_runtime_data(source, self.root / "private-role-output", validate=False)

    def test_rehashed_unreachable_query_file_rejected(self):
        source = self.changed_source("unreachable")
        def change(manifest):
            path = next(p for p, m in manifest["files"].items() if m["schema"] == "sequence_shard")
            orphan = path.rsplit("/", 1)[0] + "/" + "f" * 64 + ".json"
            shutil.copyfile(source / path, source / orphan)
            manifest["files"][orphan] = copy.deepcopy(manifest["files"][path])
        self.rehash_manifest(source, change)
        with self.assertRaises(PublicDataError):
            build_runtime_data(source, self.root / "orphan-output", validate=False)

    def test_shadowed_duplicate_source_member_rejected(self):
        source = self.changed_source("shadowed-member")
        def change(manifest):
            path = next(p for p, m in manifest["files"].items() if m["schema"] == "sequence_shard")
            raw = (source / path).read_bytes().replace(b'"sequence":', b'"sequence":"PRIVATE_CANARY","sequence":', 1)
            (source / path).write_bytes(raw)
            manifest["files"][path].update(bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest())
        self.rehash_manifest(source, change)
        with self.assertRaises(PublicDataError) as caught:
            build_runtime_data(source, self.root / "shadowed-output", validate=False)
        self.assertNotIn("PRIVATE_CANARY", str(caught.exception))

    def test_symlink_rejected(self):
        source = self.changed_source("symlink")
        (source / "link").symlink_to(self.bundle / "manifest.json")
        with self.assertRaises(PublicDataError):
            build_runtime_data(source, self.root / "symlink-output", validate=False)


class PartBoundaryTests(unittest.TestCase):
    def setUp(self):
        self.path = "releases/" + "a" * 64 + "/human/targets/" + "b" * 64 + ".json"
        self.raw = b'{"data":1.0}\n'
        self.expected = {self.path: {"bytes": len(self.raw), "sha256": hashlib.sha256(self.raw).hexdigest(), "schema": "target_shard"}}
        self.members = [[self.path, base64.b64encode(self.raw).decode("ascii")]]

    def test_preserves_raw_json_number_spelling(self):
        content, descriptor = encoded_part(self.members)
        self.assertEqual(decode_part(content, descriptor, self.expected)[self.path], self.raw)

    def test_untrusted_member_paths_and_duplicates_rejected(self):
        for members in ([self.members[0], self.members[0]], [["../private.sql", "eA=="]], [[self.path + "/extra", self.members[0][1]]]):
            with self.subTest(members=members), self.assertRaises(PublicDataError):
                decode_part(*encoded_part(members), self.expected)

    def test_bad_base64_member_hash_count_and_range_rejected(self):
        variants = []
        bad = copy.deepcopy(self.members); bad[0][1] = "!" * len(bad[0][1]); variants.append(encoded_part(bad))
        bad = copy.deepcopy(self.members); bad[0][1] = base64.b64encode(b'x' * len(self.raw)).decode(); variants.append(encoded_part(bad))
        for field, value in (("path", "../private.sql"), ("members", 2), ("min_path", "invalid"), ("max_path", "invalid"), ("decoded_bytes", 1)):
            raw, descriptor = encoded_part(self.members); descriptor[field] = value; variants.append((raw, descriptor))
        for raw, descriptor in variants:
            with self.subTest(descriptor=descriptor), self.assertRaises(PublicDataError):
                decode_part(raw, descriptor, self.expected)

    def test_duplicate_keys_inside_logical_member_rejected(self):
        raw = b'{"data":{"value":"PRIVATE_CANARY","value":1}}\n'
        members = [[self.path, base64.b64encode(raw).decode("ascii")]]
        expected = {self.path: {**self.expected[self.path], "bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()}}
        with self.assertRaises(PublicDataError) as caught:
            decode_part(*encoded_part(members), expected)
        self.assertNotIn("PRIVATE_CANARY", str(caught.exception))

    def test_gzip_bomb_and_concatenated_stream_rejected(self):
        for raw in (gzip.compress(b'x' * 1_000_000, mtime=0), encoded_part(self.members)[0] + gzip.compress(b'extra', mtime=0)):
            _, descriptor = encoded_part(self.members)
            descriptor.update(bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest())
            with self.assertRaises(PublicDataError):
                decode_part(raw, descriptor, self.expected)


if __name__ == "__main__":
    unittest.main()
