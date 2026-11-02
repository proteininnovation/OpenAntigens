"""Compile sealed local query metadata and immutable own-host target assets."""
from __future__ import annotations

from collections import Counter
import hashlib
import os
from pathlib import Path
import re
import shutil
import tempfile

from jsonschema import Draft202012Validator

from .public_api_data import (
    EXTERNAL_DOWNLOAD_LIMIT, INDEX_LIMIT, INDEX_TOTAL_LIMIT, MANIFEST_LIMIT,
    PUBLIC_DOWNLOAD_ROLES, SECTIONS, PublicDataError, _approved_release_schema, _safe_path,
    _sha_file, canonical_bytes, public_file_role, semantic_hash, strict_json_loads,
    validate_public_api_bundle, validate_retirement_metadata, archive_url_policy,
)
from .public_api_transport import build_runtime_data, _read

FORMAT = "openantigens-own-static-query-v1"
STATIC_ORIGIN = "https://openantigens.org"
METADATA_LIMIT = RECORD_LIMIT = 1024 * 1024
ROUTING_LIMIT = OFFSETS_LIMIT = 16 * 1024 * 1024
TARGET_COMPRESSED_LIMIT = 4 * 1024 * 1024
TARGET_DECODED_LIMIT = 32 * 1024 * 1024
TARGET_PART_LIMIT = 2
_TARGET = re.compile(r"releases/([a-f0-9]{64})/(human|mouse)/targets/([a-f0-9]{64})\.json")
_MEMBER = re.compile(r"releases/([a-f0-9]{64})/(human|mouse)/(?:targets/[a-f0-9]{64}\.json|constructs/[a-f0-9]{64}-(?:0|[1-9][0-9]*)\.json|sequences/[a-f0-9]{64}(?:-reference)?\.json|sections/[a-f0-9]{64}/[a-z_]+-(?:0|[1-9][0-9]*)\.json)")


def _fail(message):
    raise PublicDataError(message)


def _bytes(value):
    return canonical_bytes(value) + b"\n"


def _sha(content):
    return hashlib.sha256(content).hexdigest()


def _validators(schema_path):
    path = Path(schema_path)
    schema = strict_json_loads(_read(path.parent.resolve(), path.name, METADATA_LIMIT))
    if schema.get("x-public-text-policy") != archive_url_policy():
        _fail("Query code schema public-text policy mismatch")
    Draft202012Validator.check_schema(schema)
    def local_refs(value):
        if isinstance(value, dict):
            if "$ref" in value and not value["$ref"].startswith("#/"):
                _fail("Query code schema must use local references")
            for child in value.values():
                local_refs(child)
        elif isinstance(value, list):
            for child in value:
                local_refs(child)
    local_refs(schema)
    return {name: Draft202012Validator({**schema, "$ref": f"#/$defs/{name}"})
            for name in ("file", "seal", "query_root", "query_release", "offsets", "transport_index")}


def _check(validator, value, label):
    if not validator.is_valid(value):
        _fail(f"Invalid {label}")


def _descriptor(path, relative):
    return {"path": relative, "bytes": path.stat().st_size, "sha256": _sha_file(path)}


def _write(root, relative, content):
    path = _safe_path(root, relative, exists=False)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("xb") as handle:
        handle.write(content)
    return _descriptor(path, relative)


def _target_path(row, release_id, species):
    key = _sha(canonical_bytes(row.get("target_id") or row["record_id"]))
    expected = f"releases/{release_id}/{species}/targets/{key}.json"
    if row.get("summary_path", expected) != expected:
        _fail("Query index target path mismatch")
    return expected


def _asset_inventory(transport):
    assets = {}
    for target, parts in transport["targets"].items():
        match = _TARGET.fullmatch(target)
        if not match or not 1 <= len(parts) <= TARGET_PART_LIMIT:
            _fail("Query target exceeds its part-count budget")
        if sum(p["bytes"] for p in parts) > TARGET_COMPRESSED_LIMIT or sum(p["decoded_bytes"] for p in parts) > TARGET_DECODED_LIMIT:
            _fail("Query target exceeds its aggregate byte budget")
        previous = ""
        for number, part in enumerate(parts):
            if part["path"] != f"packs/{match[1]}/{match[2]}/{match[3]}-{number}.json.gz":
                _fail("Query part path mismatch")
            for boundary in (part["min_path"], part["max_path"]):
                member = _MEMBER.fullmatch(boundary)
                if not member or member.group(1, 2) != match.group(1, 2):
                    _fail("Query part range has foreign ownership")
                tail = boundary.split("/", 3)[3]
                pieces = tail.split("/")
                owner = None
                if pieces[0] == "targets":
                    owner = pieces[1][:-5]
                elif pieces[0] == "constructs":
                    owner = pieces[1].split("-", 1)[0]
                elif pieces[0] == "sections":
                    owner = pieces[1]
                    if pieces[2].rsplit("-", 1)[0] not in SECTIONS:
                        _fail("Query part range has an unknown section role")
                elif pieces[1].endswith("-reference.json"):
                    owner = pieces[1][:-len("-reference.json")]
                if owner is not None and owner != match[3]:
                    _fail("Query part range belongs to another target")
            if part["min_path"] > part["max_path"] or previous and previous >= part["min_path"]:
                _fail("Query part ranges overlap or are unordered")
            previous = part["max_path"]
            relative = f"api-runtime/{match[1]}/{match[2]}/{match[3]}-{number}-{part['sha256']}.oap"
            if relative in assets:
                _fail("Duplicate query asset")
            assets[relative] = {"bytes": part["bytes"], "sha256": part["sha256"]}
    return assets


def _summary_table(path, index, validator, expected, *, destination=None, offsets=None):
    checksum, result, position = hashlib.sha256(), [], 0
    if path.stat().st_size != expected["bytes"]:
        _fail("Summary JSONL size mismatch")
    with path.open("rb") as source:
        for row_number, row in enumerate(index):
            line = source.readline(RECORD_LIMIT + 1)
            if not line or len(line) > RECORD_LIMIT or not line.endswith(b"\n"):
                _fail("Summary JSONL record exceeds its bound or is incomplete")
            summary = strict_json_loads(line)
            _check(validator, summary, "released full summary")
            if any(summary.get(key) != value for key, value in row.items() if key != "summary_path"):
                _fail("Summary JSONL and compact index disagree")
            entry = [position, len(line), _sha(line)]
            if offsets is not None and (row_number >= len(offsets) or entry != offsets[row_number]):
                _fail("Summary offset, length or checksum mismatch")
            result.append(entry)
            checksum.update(line)
            position += len(line)
            if destination is not None:
                destination.write(line)
        if source.read(1) or offsets is not None and len(offsets) != len(index):
            _fail("Summary JSONL coverage mismatch")
    if position != expected["bytes"] or checksum.hexdigest() != expected["sha256"]:
        _fail("Summary JSONL checksum mismatch")
    return result


def _inventory(root):
    files = {}
    for path in root.rglob("*"):
        if path.is_symlink():
            _fail("Query output contains a symlink")
        if path.is_file():
            files[path.relative_to(root).as_posix()] = {"bytes": path.stat().st_size, "sha256": _sha_file(path)}
    return files


def validate_query_runtime(runtime_dir, *, schema_path, seal):
    """Validate local pins and complete query coverage; return authorized asset hashes."""
    root_arg = Path(runtime_dir)
    if root_arg.is_symlink() or not root_arg.is_dir():
        _fail("Invalid query runtime root")
    root = root_arg.resolve()
    validators = _validators(schema_path)
    _check(validators["seal"], seal, "query seal")
    raw = _read(root, "query-root.json", METADATA_LIMIT)
    if _sha(raw) != seal["query_root_sha256"]:
        _fail("Query root pin mismatch")
    query = strict_json_loads(raw)
    _check(validators["query_root"], query, "query root")
    validate_retirement_metadata(query)
    if query["source_root_manifest_sha256"] != seal["source_root_manifest_sha256"] or query["static_origin"] != seal["static_origin"]:
        _fail("Query source or origin pin mismatch")
    allowed = {"query-root.json"}
    def load(descriptor, maximum):
        _check(validators["file"], {k: descriptor[k] for k in ("path", "bytes", "sha256")}, "query file descriptor")
        allowed.add(descriptor["path"])
        return strict_json_loads(_read(root, descriptor["path"], maximum, descriptor))
    if query["transport_index"]["path"] != "transport-index.json":
        _fail("Unexpected query routing path")
    transport = load(query["transport_index"], ROUTING_LIMIT)
    _check(validators["transport_index"], transport, "query routing inventory")
    if transport["root_manifest_sha256"] != query["source_root_manifest_sha256"]:
        _fail("Query transport source mismatch")
    assets = _asset_inventory(transport)
    release_ids = [p["release_id"] for p in query["releases"]]
    if len(set(release_ids)) != len(release_ids) or query["current_release_id"] not in release_ids:
        _fail("Query retained release identity mismatch")
    retired = set()
    for item in query["retired_releases"]:
        if item["release_id"] in retired or item["release_id"] in release_ids or item["release_id"] == item["retired_by_release_id"]:
            _fail("Query retirement identity mismatch")
        retired.add(item["release_id"])
    covered, total_indexes = set(), 0
    for pointer in query["releases"]:
        rid = pointer["release_id"]
        if pointer["path"] != f"releases/{rid}/query-release.json":
            _fail("Unexpected query release path")
        release = load(pointer, METADATA_LIMIT)
        _check(validators["query_release"], release, "query release metadata")
        if release["release_id"] != rid or release["source_manifest"]["url"] != f"{STATIC_ORIGIN}/data/releases/{rid}/manifest.json":
            _fail("Query source release mismatch")
        provenance = release["provenance"]["input_hashes"]
        if provenance["manifest_sha256"] != release["source_manifest"]["sha256"] or provenance["manifest_url"] != release["source_manifest"]["url"]:
            _fail("Query provenance manifest mismatch")
        if release["public_schema"]["path"] not in {f"releases/{rid}/{sp}/downloads/schema.json" for sp in ("human", "mouse")}:
            _fail("Unexpected public schema path")
        public_schema = load(release["public_schema"], METADATA_LIMIT)
        _approved_release_schema(public_schema)
        if semantic_hash(public_schema) != release["source_payload_sha256"]["schema"] or semantic_hash(release["citation"]) != release["source_payload_sha256"]["citation"]:
            _fail("Query schema or citation source identity mismatch")
        public = {name: Draft202012Validator({**public_schema, "$ref": f"#/$defs/{name}"}) for name in ("index", "summary", "download")}
        for species, taxon in (("human", 9606), ("mouse", 10090)):
            sm = release["species"][species]
            prefix = f"releases/{rid}/{species}"
            if sm["taxon_id"] != taxon or sm["index"]["path"] != prefix + "/index.json" or sm["summary_lines"]["path"] != prefix + "/downloads/index.jsonl" or sm["summary_lines"]["offsets"]["path"] != prefix + "/summary-offsets.json":
                _fail("Query species file identity mismatch")
            index = load(sm["index"], INDEX_LIMIT)
            _check(public["index"], index, "released compact index")
            if index["release_id"] != rid or index["species"] != species:
                _fail("Query index envelope mismatch")
            rows = index["data"]
            total_indexes += sm["index"]["bytes"]
            if total_indexes > INDEX_TOTAL_LIMIT:
                _fail("Query retained index budget exceeded")
            paths = {_target_path(row, rid, species) for row in rows}
            ids = [row["target_id"] for row in rows if row["target_id"] is not None]
            if len(paths) != len(rows) or len(set(ids)) != len(ids) or any(row["species"] != species or row["taxon_id"] != taxon for row in rows):
                _fail("Duplicate or foreign query target identity")
            if sm["counts"]["targets"] != len(rows) or sm["counts"]["resolved_targets"] != len(ids) or sm["counts"]["status"] != dict(Counter(row["status"] for row in rows)):
                _fail("Query index aggregate counts mismatch")
            covered.update(paths)
            names = set()
            for item in sm["downloads"]:
                _check(public["download"], item, "released download metadata")
                if item["name"] in names or item["path"] != prefix + "/downloads/" + item["name"] or item["url"] != f"{STATIC_ORIGIN}/data/releases/{rid}/{species}/{item['name']}":
                    _fail("Query download identity mismatch")
                names.add(item["name"])
                public_file_role(item["path"], item, rid)
                if item["name"] == "index.jsonl" and any(item[k] != sm["summary_lines"][k] for k in ("path", "bytes", "sha256")):
                    _fail("Query local summary and public download disagree")
            if names != set(PUBLIC_DOWNLOAD_ROLES):
                _fail("Query download roles incomplete")
            offsets = load(sm["summary_lines"]["offsets"], OFFSETS_LIMIT)
            _check(validators["offsets"], offsets, "query summary offsets")
            if offsets["index_sha256"] != sm["index"]["sha256"] or offsets["jsonl_sha256"] != sm["summary_lines"]["sha256"]:
                _fail("Query summary offset source mismatch")
            lines = sm["summary_lines"]
            if lines["bytes"] > EXTERNAL_DOWNLOAD_LIMIT:
                _fail("Query summary file budget exceeded")
            allowed.add(lines["path"])
            _summary_table(_safe_path(root, lines["path"]), rows, public["summary"], lines, offsets=offsets["rows"])
    if covered != set(transport["targets"]):
        _fail("Query routing and index target coverage differ")
    if set(_inventory(root)) != allowed:
        _fail("Unexpected query runtime file")
    directories = {p.as_posix() for name in allowed for p in Path(name).parents if p.as_posix() != "."}
    if {p.relative_to(root).as_posix() for p in root.rglob("*") if p.is_dir()} != directories:
        _fail("Unexpected query runtime directory")
    return assets


def build_query_runtime(bundle_dir, runtime_output, static_output, *, schema_path):
    """Fully validate source data and build fresh local app-data and static assets."""
    source = Path(bundle_dir)
    if source.is_symlink() or not source.is_dir():
        _fail("Invalid public source root")
    validated_root = validate_public_api_bundle(source)
    return _build_query_runtime(source, runtime_output, static_output, schema_path=schema_path, validated_root=validated_root)


def _build_query_runtime(bundle_dir, runtime_output, static_output, *, schema_path, validated_root):
    """Internal packager path after complete source/static validation; no CLI bypass."""
    source = Path(bundle_dir).resolve()
    outputs = [Path(runtime_output).absolute(), Path(static_output).absolute()]
    if any(p.exists() or p.is_symlink() or not p.parent.is_dir() for p in outputs):
        _fail("Query outputs must be absent with existing parent directories")
    outputs = [p.resolve() for p in outputs]
    if any(p.is_relative_to(source) or source.is_relative_to(p) for p in outputs) or outputs[0].is_relative_to(outputs[1]) or outputs[1].is_relative_to(outputs[0]):
        _fail("Query output roots must be separate from each other and source")
    validators = _validators(schema_path)
    source_root_bytes = _read(source, "manifest.json", METADATA_LIMIT)
    root = strict_json_loads(source_root_bytes)
    if root != validated_root:
        _fail("Validated source root changed before query compilation")
    validate_retirement_metadata(root)
    with tempfile.TemporaryDirectory(prefix="openantigens-query-", dir=outputs[0].parent) as temporary:
        work = Path(temporary).resolve()
        transport_root, runtime, static = work / "transport", work / "runtime", work / "static"
        build_runtime_data(source, transport_root, validate=False)
        transport_raw = _read(transport_root, "transport-index.json", ROUTING_LIMIT)
        transport = strict_json_loads(transport_raw)
        assets = _asset_inventory(transport)
        _check(validators["transport_index"], transport, "query routing inventory")
        runtime.mkdir(); static.mkdir()
        transport_descriptor = _write(runtime, "transport-index.json", transport_raw)
        query = {"format": FORMAT, "source_root_manifest_sha256": _sha(source_root_bytes), "current_release_id": root["current_release_id"], "retired_releases": root.get("retired_releases", []), "transport_index": transport_descriptor, "releases": [], "static_origin": STATIC_ORIGIN}
        for pointer in root["releases"]:
            raw = _read(source, pointer["manifest_path"], MANIFEST_LIMIT)
            if _sha(raw) != pointer["sha256"]:
                _fail("Source release manifest changed during query compilation")
            manifest = strict_json_loads(raw)
            del raw
            rid, provenance = pointer["release_id"], manifest["provenance"]
            if manifest["manifest_url"] != f"{STATIC_ORIGIN}/data/releases/{rid}/manifest.json":
                _fail("Query runtime requires the approved own static origin")
            compact_provenance = {k: provenance[k] for k in ("software_commit", "snapshot_id", "source_versions", "source_dates", "notes")}
            compact_provenance.update(projection_methods=provenance.get("projection_methods"), input_hashes={"count": len(provenance["input_sha256"]), "manifest_url": manifest["manifest_url"], "manifest_sha256": pointer["sha256"], "json_pointer": "/provenance/input_sha256"})
            compact = {"format": "openantigens-query-release-v1", "release_id": rid, "schema_version": manifest["schema_version"], "source_manifest": {"url": manifest["manifest_url"], "sha256": pointer["sha256"]}, "source_payload_sha256": manifest["core"]["payload_sha256"], "citation": manifest["citation"], "licenses": manifest["licenses"], "provenance": compact_provenance, "species": {}}
            schema_file = next(p for p, meta in manifest["files"].items() if meta["media_type"] == "application/schema+json")
            schema_raw = _read(source, schema_file, METADATA_LIMIT, manifest["files"][schema_file])
            compact["public_schema"] = _write(runtime, schema_file, schema_raw)
            public_schema = strict_json_loads(schema_raw)
            _approved_release_schema(public_schema)
            summary_validator = Draft202012Validator({**public_schema, "$ref": "#/$defs/summary"})
            for species in ("human", "mouse"):
                sm = manifest["species"][species]
                prefix = f"releases/{rid}/{species}"
                index_raw = _read(source, sm["index_path"], INDEX_LIMIT, manifest["files"][sm["index_path"]])
                index = strict_json_loads(index_raw)["data"]
                index_descriptor = _write(runtime, sm["index_path"], index_raw)
                summary_path = prefix + "/downloads/index.jsonl"
                destination = _safe_path(runtime, summary_path, exists=False)
                destination.parent.mkdir(parents=True, exist_ok=True)
                with destination.open("xb") as handle:
                    rows = _summary_table(_safe_path(source, summary_path), index, summary_validator, manifest["files"][summary_path], destination=handle)
                summaries = _descriptor(destination, summary_path)
                offset_data = {"format": "openantigens-summary-offsets-v1", "index_sha256": index_descriptor["sha256"], "jsonl_sha256": summaries["sha256"], "rows": rows}
                _check(validators["offsets"], offset_data, "compiled summary offsets")
                offset_raw = _bytes(offset_data)
                if len(offset_raw) > OFFSETS_LIMIT:
                    _fail("Query summary offsets exceed their byte budget")
                summaries["offsets"] = _write(runtime, prefix + "/summary-offsets.json", offset_raw)
                compact["species"][species] = {"taxon_id": sm["taxon_id"], "index": index_descriptor, "counts": sm["counts"], "capabilities": sm["capabilities"], "downloads": sm["downloads"], "records_download_name": "records.jsonl.gz", "summary_lines": summaries}
            _check(validators["query_release"], compact, "compiled query release")
            content = _bytes(compact)
            if len(content) > METADATA_LIMIT:
                _fail("Query release metadata exceeds its byte budget")
            descriptor = _write(runtime, f"releases/{rid}/query-release.json", content)
            query["releases"].append({"release_id": rid, **descriptor})
        _check(validators["query_root"], query, "compiled query root")
        content = _bytes(query)
        if len(content) > METADATA_LIMIT:
            _fail("Query root exceeds its byte budget")
        root_descriptor = _write(runtime, "query-root.json", content)
        seal = {"format": FORMAT, "query_root_sha256": root_descriptor["sha256"], "source_root_manifest_sha256": query["source_root_manifest_sha256"], "static_origin": STATIC_ORIGIN}
        for target, parts in transport["targets"].items():
            match = _TARGET.fullmatch(target)
            for number, part in enumerate(parts):
                relative = f"api-runtime/{match[1]}/{match[2]}/{match[3]}-{number}-{part['sha256']}.oap"
                content = _read(transport_root, part["path"], TARGET_COMPRESSED_LIMIT, part)
                _write(static, relative, content)
        checked_assets = validate_query_runtime(runtime, schema_path=schema_path, seal=seal)
        runtime_files, static_files = _inventory(runtime), _inventory(static)
        if checked_assets != assets or static_files != assets:
            _fail("Compiled query asset inventory mismatch")
        if _read(source, "manifest.json", METADATA_LIMIT) != source_root_bytes:
            _fail("Source root changed during query compilation")
        if any(p.exists() or p.is_symlink() for p in outputs):
            _fail("Query output appeared during compilation")
        completed = []
        try:
            for staged, destination in zip((runtime, static), outputs):
                os.rename(staged, destination)
                completed.append(destination)
        except OSError:
            for created in completed:
                shutil.rmtree(created)
            raise
    return {"seal": seal, "runtime_files": runtime_files, "static_assets": static_files, "metrics": {"runtime_file_count": len(runtime_files), "runtime_bytes": sum(v["bytes"] for v in runtime_files.values()), "static_asset_count": len(static_files), "static_asset_bytes": sum(v["bytes"] for v in static_files.values()), "target_count": len(transport["targets"]), "release_count": len(query["releases"]), "max_parts_per_target": max(map(len, transport["targets"].values()), default=0), "max_target_compressed_bytes": max((sum(p["bytes"] for p in parts) for parts in transport["targets"].values()), default=0), "max_target_decoded_bytes": max((sum(p["decoded_bytes"] for p in parts) for parts in transport["targets"].values()), default=0)}}
