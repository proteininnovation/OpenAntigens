"""Bounded physical storage for validated immutable public API bundles."""
from __future__ import annotations

import base64
import gzip
import hashlib
import json
from pathlib import Path
import re
import zlib

from .public_api_data import (MANIFEST_LIMIT, PublicDataError, _safe_path, _sha_file, canonical_bytes,
                              public_file_role, validate_public_api_bundle)

FORMAT = "openantigens-target-gzip-v1"
INDEX_LIMIT = PART_LIMIT = 16 * 1024 * 1024
COMPRESSED_LIMIT = 17 * 1024 * 1024
TARGET_LIMIT = 8 * 1024 * 1024
PAGE_LIMIT = 1024 * 1024
STANDALONE_DOWNLOADS = frozenset({"schema.json", "citation.json", "openantigens.bib", "openantigens.ris"})
EXTERNAL_DOWNLOADS = frozenset({"index.jsonl", "index.tsv", "constructs.jsonl", "constructs.tsv", "constructs.fasta", "reference.fasta", "records.jsonl.gz"})
_TARGET = re.compile(r"releases/([a-f0-9]{64})/(human|mouse)/targets/([a-f0-9]{64})\.json")
_HASH = re.compile(r"[a-f0-9]{64}")
_PACK = re.compile(r"packs/[a-f0-9]{64}/(?:human|mouse)/[a-f0-9]{64}-(?:0|[1-9][0-9]*)\.json\.gz")


def _fail(message: str) -> None:
    raise PublicDataError(message)


def _json(content: bytes):
    def pairs(items):
        value = {}
        for key, item in items:
            if key in value:
                _fail("Duplicate transport JSON key")
            value[key] = item
        return value
    try:
        return json.loads(content.decode("utf-8"), object_pairs_hook=pairs)
    except (UnicodeError, ValueError) as exc:
        raise PublicDataError("Invalid transport JSON") from exc


def _read(root: Path, path: str, maximum: int, descriptor: dict | None = None) -> bytes:
    source = _safe_path(root, path)
    size = source.stat().st_size
    if size > maximum or (descriptor is not None and size != descriptor["bytes"]):
        _fail("Transport file size mismatch")
    with source.open("rb") as handle:
        content = handle.read(maximum + 1)
    if len(content) != size or len(content) > maximum:
        _fail("Transport file changed or exceeds its bound")
    if descriptor is not None and hashlib.sha256(content).hexdigest() != descriptor["sha256"]:
        _fail("Transport file checksum mismatch")
    return content


def _part_path(target_path: str, number: int) -> str:
    match = _TARGET.fullmatch(target_path)
    if not match:
        _fail("Invalid transport target path")
    return f"packs/{match[1]}/{match[2]}/{match[3]}-{number}.json.gz"


def _members_bytes(members: list[list[str]]) -> bytes:
    return canonical_bytes({"format": FORMAT, "members": members}) + b"\n"


def decode_part(content: bytes, descriptor: dict, expected: dict[str, dict]) -> dict[str, bytes]:
    """Validate one bounded container against source-authoritative member metadata."""
    keys = {"path", "sha256", "bytes", "decoded_bytes", "members", "min_path", "max_path"}
    if not isinstance(descriptor, dict) or set(descriptor) != keys:
        _fail("Invalid transport part descriptor")
    if not isinstance(descriptor["path"], str) or not _PACK.fullmatch(descriptor["path"]):
        _fail("Invalid transport physical path")
    if not _HASH.fullmatch(str(descriptor["sha256"])):
        _fail("Invalid transport checksum")
    if any(type(descriptor[k]) is not int or descriptor[k] <= 0 for k in ("bytes", "decoded_bytes", "members")):
        _fail("Invalid transport part bounds")
    if descriptor["bytes"] > COMPRESSED_LIMIT or descriptor["decoded_bytes"] > PART_LIMIT:
        _fail("Transport part exceeds its bound")
    if len(content) != descriptor["bytes"] or hashlib.sha256(content).hexdigest() != descriptor["sha256"]:
        _fail("Transport compressed checksum mismatch")
    try:
        decoder = zlib.decompressobj(31)
        raw = decoder.decompress(content, descriptor["decoded_bytes"] + 1)
        if len(raw) != descriptor["decoded_bytes"] or not decoder.eof or decoder.unused_data or decoder.unconsumed_tail:
            _fail("Transport gzip bounds or trailing data mismatch")
    except zlib.error as exc:
        raise PublicDataError("Invalid transport gzip") from exc
    payload = _json(raw)
    if not isinstance(payload, dict) or set(payload) != {"format", "members"} or payload["format"] != FORMAT or not isinstance(payload["members"], list):
        _fail("Invalid transport member container")
    entries = payload["members"]
    if not entries or len(entries) != descriptor["members"]:
        _fail("Transport member count mismatch")
    result, previous = {}, None
    for item in entries:
        if not isinstance(item, list) or len(item) != 2 or not all(isinstance(v, str) for v in item):
            _fail("Invalid transport member")
        path, encoded = item
        # Embedded paths cannot authorize members outside the validated source graph.
        if path not in expected or previous is not None and path <= previous:
            _fail("Unexpected, duplicate or unordered transport member")
        previous = path
        meta = expected[path]
        bound = TARGET_LIMIT if meta.get("schema") == "target_shard" else PAGE_LIMIT
        if type(meta.get("bytes")) is not int or not 0 <= meta["bytes"] <= bound:
            _fail("Transport logical member exceeds its bound")
        if len(encoded) != 4 * ((meta["bytes"] + 2) // 3):
            _fail("Transport base64 length mismatch")
        try:
            value = base64.b64decode(encoded, validate=True)
        except ValueError as exc:
            raise PublicDataError("Invalid transport base64") from exc
        if base64.b64encode(value).decode("ascii") != encoded or len(value) != meta["bytes"] or hashlib.sha256(value).hexdigest() != meta["sha256"]:
            _fail("Transport logical checksum mismatch")
        result[path] = value
    if entries[0][0] != descriptor["min_path"] or entries[-1][0] != descriptor["max_path"]:
        _fail("Transport member range mismatch")
    return result


def build_runtime_data(bundle_dir: str | Path, output_dir: str | Path, *, validate: bool = True, part_limit: int = PART_LIMIT) -> dict:
    """Write a fresh physical runtime copy; never modify the authoritative bundle."""
    source_arg, output_arg = Path(bundle_dir), Path(output_dir)
    if source_arg.is_symlink() or output_arg.is_symlink() or not source_arg.is_dir():
        _fail("Invalid transport input/output roots")
    source, output = source_arg.resolve(), output_arg.resolve()
    if output.is_relative_to(source) or source.is_relative_to(output):
        _fail("Transport output must be separate from the source")
    if output.exists() and any(output.iterdir()):
        _fail("Transport output must be absent or empty")
    if type(part_limit) is not int or not 1 <= part_limit <= PART_LIMIT:
        _fail("Invalid transport part budget")
    if validate:
        validate_public_api_bundle(source)
    root_bytes = _read(source, "manifest.json", PAGE_LIMIT)
    root_sha = hashlib.sha256(root_bytes).hexdigest()
    root = _json(root_bytes)
    if not isinstance(root, dict) or not isinstance(root.get("releases"), list) or not 1 <= len(root["releases"]) <= 2:
        _fail("Invalid transport source root manifest")
    root_sidecar = _read(source, "manifest.sha256", 256)
    if root_sidecar != f"{root_sha}  manifest.json\n".encode("ascii"):
        _fail("Transport source root checksum mismatch")
    standalone = {"manifest.json": root_bytes, "manifest.sha256": root_sidecar}
    allowed, packed = set(standalone), set()
    target_groups, manifests, manifest_checks = {}, [], {}
    external = 0
    for pointer in root["releases"]:
        rid, path = pointer.get("release_id"), pointer.get("manifest_path")
        if not isinstance(rid, str) or not _HASH.fullmatch(rid) or path != f"releases/{rid}/manifest.json":
            _fail("Invalid transport release pointer")
        raw = _read(source, path, MANIFEST_LIMIT)
        if hashlib.sha256(raw).hexdigest() != pointer.get("sha256"):
            _fail("Transport release checksum mismatch")
        manifest = _json(raw)
        if manifest.get("release_id") != rid or not isinstance(manifest.get("files"), dict):
            _fail("Invalid transport release manifest")
        standalone[path] = None
        manifest_checks[path] = {"bytes": len(raw), "sha256": pointer["sha256"]}
        del raw
        allowed.add(path)
        files = manifest["files"]
        roles = {name: public_file_role(name, meta, rid) for name, meta in files.items()}
        for name, meta in files.items():
            file = _safe_path(source, name)
            if file.stat().st_size != meta["bytes"] or _sha_file(file) != meta["sha256"]:
                _fail("Transport source file checksum mismatch")
            allowed.add(name)
            role = roles[name]
            if role == "index" or role == "download" and name.rsplit("/", 1)[-1] in STANDALONE_DOWNLOADS:
                standalone[name] = None
            elif role == "download":
                if name.rsplit("/", 1)[-1] not in EXTERNAL_DOWNLOADS:
                    _fail("Unapproved external transport download")
                external += 1
        for name in sorted(n for n, role in roles.items() if role == "target"):
            target = _json(_read(source, name, TARGET_LIMIT, files[name]))["data"]
            references = [(name, "target")]
            if target["reference_sequence_path"]:
                references.append((target["reference_sequence_path"], "sequence"))
            references.extend((p, "sequence") for p in target["construct_sequence_paths"].values() if p is not None)
            references.extend((p, "constructs") for p in target["constructs"]["pages"])
            references.extend((p, "section") for value in target["sections"].values() for p in value["pages"])
            names = [p for p, _ in references]
            species_prefix = name.split("/targets/", 1)[0] + "/"
            if len(names) != len(set(names)) or any(p in packed or not p.startswith(species_prefix) or roles.get(p) != role for p, role in references):
                _fail("Invalid or duplicate transport target reachability")
            packed.update(names)
            target_groups[name] = (sorted(names), files)
        if {p for p, role in roles.items() if role not in {"index", "download"}} != packed.intersection(files):
            _fail("Unreachable public query file in transport source")
        manifests.append(manifest)
    actual = set()
    for path in source.rglob("*"):
        if path.is_symlink():
            _fail("Transport source contains a symlink")
        if path.is_file():
            actual.add(path.relative_to(source).as_posix())
    if actual != allowed:
        _fail("Unexpected file in transport source")
    output.mkdir(parents=True, exist_ok=True)
    index = {"format": FORMAT, "root_manifest_sha256": root_sha, "targets": {}}
    physical = {}

    def write(name: str, content: bytes) -> None:
        destination = _safe_path(output, name, exists=False)
        destination.parent.mkdir(parents=True, exist_ok=True)
        with destination.open("xb") as handle:
            handle.write(content)
        physical[name] = len(content)

    for name, content in standalone.items():
        if content is None:
            meta = manifest_checks.get(name) or next(m["files"][name] for m in manifests if name in m["files"])
            bound = MANIFEST_LIMIT if name in manifest_checks else 32 * 1024 * 1024 if name.endswith("/index.json") else PAGE_LIMIT
            content = _read(source, name, bound, meta)
        write(name, content)
    empty_size = len(_members_bytes([]))
    for target_path, (names, files) in sorted(target_groups.items()):
        descriptors, members, size = [], [], empty_size

        def flush() -> None:
            if not members:
                return
            raw = _members_bytes(members)
            if len(raw) != size or len(raw) > part_limit:
                _fail("Transport part exceeds its decoded budget")
            compressed = gzip.compress(raw, compresslevel=9, mtime=0)
            if len(compressed) > COMPRESSED_LIMIT:
                _fail("Transport compressed part exceeds its budget")
            path = _part_path(target_path, len(descriptors))
            descriptor = {"path": path, "sha256": hashlib.sha256(compressed).hexdigest(), "bytes": len(compressed), "decoded_bytes": len(raw), "members": len(members), "min_path": members[0][0], "max_path": members[-1][0]}
            decode_part(compressed, descriptor, {p: files[p] for p, _ in members})
            write(path, compressed)
            descriptors.append(descriptor)

        for name in names:
            meta = files[name]
            content = _read(source, name, TARGET_LIMIT if meta["schema"] == "target_shard" else PAGE_LIMIT, meta)
            member = [name, base64.b64encode(content).decode("ascii")]
            added = len(canonical_bytes(member))
            if empty_size + added > part_limit:
                _fail("Single transport member exceeds the part budget")
            if members and size + 1 + added > part_limit:
                flush()
                members, size = [], empty_size
            size += added + bool(members)
            members.append(member)
        flush()
        index["targets"][target_path] = descriptors
    index_bytes = canonical_bytes(index) + b"\n"
    if len(index_bytes) > INDEX_LIMIT:
        _fail("Transport index exceeds its budget")
    write("transport-index.json", index_bytes)
    if _read(source, "manifest.json", PAGE_LIMIT) != root_bytes:
        _fail("Source root changed while preparing transport")
    return {"format": FORMAT, "root_manifest_sha256": root_sha, "transport_index_sha256": hashlib.sha256(index_bytes).hexdigest(), "transport_index_bytes": len(index_bytes), "physical_file_count": len(physical), "physical_bytes": sum(physical.values()), "packed_logical_file_count": len(packed), "external_download_count": external}
