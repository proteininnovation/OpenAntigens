"""Build a bounded, deterministic public-service ZIP from validated public data.

This module packages an existing data release. It never projects new scientific
records, contacts a host, or publishes an artifact.
"""

from __future__ import annotations

import hashlib
import codecs
import json
import os
import re
import shutil
import stat
import subprocess
import tempfile
import zipfile
from pathlib import Path, PurePosixPath
from typing import Any
from contextlib import contextmanager
from urllib.parse import unquote, urlsplit

HARD_ZIP_LIMIT = 100_000_000
SOURCE_FILES = (
    "package.json", "package-lock.json", "src/server.js", "src/mcp.js",
    "src/provider.js", "src/contracts.js", "src/openapi.js",
    "src/analytics.js",
    "src/client-identity.js",
    "src/storage.js",
    "src/json.js",
    "src/public-boundary.js",
    "src/query-storage.js",
    "src/query-catalog.js",
    "src/query-schema.json",
)
LICENSE_FILES = ("LICENSE", "DATA_LICENSE.md")
HOST_GATES = (
    "subscription_app_slot_and_cost", "node22_runtime_and_port", "domain_tls_proxy",
    "host_disk_inodes_build_memory_runtime_memory", "startup_and_replica_limits",
    "maintenance_control", "deployment_rollback_and_old_traffic",
    "preview_visibility_expiry_and_publish_terms",
    "proxy_rate_limits_and_logs", "static_capacity_and_rollback",
    "rights_review_and_named_owner", "final_domain_http_and_mcp",
)
_SHA = re.compile(r"[0-9a-f]{64}\Z")
DOWNLOAD_ROLES = {
    "index.jsonl": ("application/x-ndjson", "summary", None),
    "index.tsv": ("text/tab-separated-values", None, None),
    "constructs.jsonl": ("application/x-ndjson", "construct", None),
    "constructs.tsv": ("text/tab-separated-values", None, None),
    "constructs.fasta": ("text/x-fasta", None, None),
    "reference.fasta": ("text/x-fasta", None, None),
    "schema.json": ("application/schema+json", None, None),
    "citation.json": ("application/json", "citation", None),
    "openantigens.bib": ("application/x-bibtex", None, None),
    "openantigens.ris": ("application/x-research-info-systems", None, None),
    "records.jsonl.gz": ("application/x-ndjson", None, "gzip"),
}
SECTIONS = {"topology", "domains", "ptms", "orthologs", "paralogs", "family", "assembly", "structures", "similarity", "literature", "disease"}
EXTERNAL_BULK_ROLES = {"index.jsonl", "index.tsv", "constructs.jsonl", "constructs.tsv", "constructs.fasta", "reference.fasta", "records.jsonl.gz"}


def _json_bytes(value: Any) -> bytes:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")


def _hash_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _root(path: Path | str) -> Path:
    path = Path(path).expanduser().absolute()
    if path.is_symlink() or not path.is_dir():
        raise ValueError("input root must be a real directory")
    return path.resolve()


def _canonical_repo(repo: Path) -> None:
    origin = subprocess.run(["git", "-C", str(repo), "remote", "get-url", "origin"], capture_output=True, text=True, check=False)
    if origin.returncode or origin.stdout.strip() not in {
        "https://github.com/proteininnovation/OpenAntigens", "https://github.com/proteininnovation/OpenAntigens.git",
        "git@github.com:proteininnovation/OpenAntigens.git", "ssh://git@github.com/proteininnovation/OpenAntigens.git",
    }:
        raise ValueError("origin must resolve to proteininnovation/OpenAntigens")


def _relative(name: str) -> str:
    if not isinstance(name, str) or not name or "\\" in name or "%" in name or "\x00" in name:
        raise ValueError("unsafe manifest path")
    parsed = PurePosixPath(name)
    if parsed.is_absolute() or parsed.as_posix() != name or any(part in {".", ".."} or part.startswith(".") for part in parsed.parts):
        raise ValueError("unsafe manifest path")
    if ":" in name or any(ord(char) < 32 for char in name):
        raise ValueError("unsafe manifest path")
    return name


def _file(root: Path, name: str) -> Path:
    name = _relative(name)
    candidate = root / name
    for part in (candidate, *candidate.parents):
        if part == root:
            break
        if part.is_symlink():
            raise ValueError("symbolic links are forbidden in release inputs")
    if not candidate.is_file() or not stat.S_ISREG(candidate.stat().st_mode) or not candidate.resolve().is_relative_to(root):
        raise ValueError(f"missing or unsafe public file: {name}")
    return candidate


def _all_files(root: Path) -> set[str]:
    names = set()
    for path in root.rglob("*"):
        if path.is_symlink():
            raise ValueError("symbolic links are forbidden in the data bundle")
        if path.is_file():
            names.add(_relative(path.relative_to(root).as_posix()))
        elif not path.is_dir():
            raise ValueError("special files are forbidden in the data bundle")
    return names


def _read_json(path: Path) -> dict:
    from agdesign2.public_api_data import strict_json_loads

    value = strict_json_loads(path.read_bytes())
    if not isinstance(value, dict):
        raise ValueError("manifest must be a JSON object")
    return value


def _check_bytes(path: Path, entry: dict) -> None:
    if type(entry.get("bytes")) is not int or entry["bytes"] < 0 or not _SHA.fullmatch(str(entry.get("sha256", ""))):
        raise ValueError("invalid file checksum metadata")
    if path.stat().st_size != entry["bytes"] or _hash_file(path) != entry["sha256"]:
        raise ValueError("public file size or checksum mismatch")


def _scan_public_text(path: Path) -> None:
    from agdesign2.public_api_data import _scan_public

    decoder = codecs.getincrementaldecoder("utf-8")()
    previous = ""
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            decoded = decoder.decode(chunk)
            _scan_public(previous + decoded)
            previous = (previous + decoded)[-256:]
    _scan_public(previous + decoder.decode(b"", final=True))


def _public_file_limit(name: str, media_type: str) -> int:
    if "/downloads/" in name:
        role = PurePosixPath(name).name
        if role in EXTERNAL_BULK_ROLES:
            return 1024 * 1024 * 1024
        if role in {"schema.json", "citation.json", "openantigens.bib", "openantigens.ris"}:
            return 1024 * 1024
        raise ValueError("unapproved public download role")
    if name.endswith("/index.json"):
        return 32 * 1024 * 1024
    if "/targets/" in name:
        return 8 * 1024 * 1024
    if media_type == "application/schema+json" or media_type == "application/json" and "/downloads/" not in name:
        return 1024 * 1024
    return 1024 * 1024


def _check_public_file_roles(bundle: Path, manifest: dict) -> None:
    """Require approved file roles and actual reachability, independent of membership."""
    from agdesign2.public_api_data import _approved_release_schema, _scan_public, _validate_release_value, public_citation_downloads, strict_json_loads

    release_id = manifest["release_id"]
    _scan_public(manifest)
    expected = set()

    def role(path: str, media: str, schema_name: str | None, encoding: str | None = None) -> Path:
        entry = manifest["files"].get(path)
        if not isinstance(entry, dict) or (entry.get("media_type"), entry.get("schema"), entry.get("encoding")) != (media, schema_name, encoding):
            raise ValueError("public file does not have its approved artifact role")
        expected.add(path)
        file = _file(bundle, path)
        limit = _public_file_limit(path, media)
        size = file.stat().st_size
        if size > limit:
            label = PurePosixPath(path).name if "/downloads/" in path else schema_name
            raise ValueError(f"public {label} is {size} bytes; limit is {limit} bytes")
        return file

    schema_bytes = role(f"releases/{release_id}/human/downloads/schema.json", "application/schema+json", None).read_bytes()
    schema = strict_json_loads(schema_bytes)
    _approved_release_schema(schema)
    citations = public_citation_downloads(manifest["citation"])

    def typed(path: str, schema_name: str, species: str) -> object:
        value = _read_json(role(path, "application/json", schema_name))
        _validate_release_value(schema, schema_name, value)
        _scan_public(value)
        if value.get("release_id") != release_id or value.get("species") != species:
            raise ValueError("reachable public shard release/species mismatch")
        return value["data"]

    for species in ("human", "mouse"):
        prefix = f"releases/{release_id}/{species}"
        metadata = manifest["species"][species]
        if metadata["index_path"] != f"{prefix}/index.json":
            raise ValueError("public index has an unapproved path")
        rows = typed(metadata["index_path"], "index", species)
        for row in rows:
            key = hashlib.sha256(_json_bytes(row["target_id"] or row["record_id"])).hexdigest()
            path = f"{prefix}/targets/{key}.json"
            if row.get("summary_path", path) != path:
                raise ValueError("public summary reference has an unapproved path")
            if row["target_id"] is None and path not in manifest["files"] and schema.get("x-openantigens-index-version", 1) == 1:
                continue
            if row["target_id"] is not None and metadata["target_paths"].get(row["target_id"]) != path:
                raise ValueError("public target reference has an unapproved path")
            target = typed(path, "target_shard", species)
            if target["summary"]["record_id"] != row["record_id"] or target["summary"]["target_id"] != row["target_id"]:
                raise ValueError("public target does not match its indexed identity")
            if set(target["sections"]) != SECTIONS:
                raise ValueError("public target evidence roles differ from the approved set")
            for section, descriptor in target["sections"].items():
                for number, page in enumerate(descriptor["pages"]):
                    if page != f"{prefix}/sections/{key}/{section}-{number}.json":
                        raise ValueError("public evidence page has an unapproved path")
                    typed(page, f"section_{section}", species)
            for number, page in enumerate(target["constructs"]["pages"]):
                if page != f"{prefix}/constructs/{key}-{number}.json":
                    raise ValueError("public construct page has an unapproved path")
                typed(page, "construct_page", species)
            reference = target["reference_sequence_path"]
            if reference is not None:
                if reference != f"{prefix}/sequences/{key}-reference.json":
                    raise ValueError("public reference sequence has an unapproved path")
                typed(reference, "sequence_shard", species)
            for construct_id, sequence in target["construct_sequence_paths"].items():
                if not _SHA.fullmatch(construct_id):
                    raise ValueError("public construct sequence identity is invalid")
                if sequence is not None:
                    if sequence != f"{prefix}/sequences/{construct_id}.json":
                        raise ValueError("public construct sequence has an unapproved path")
                    typed(sequence, "sequence_shard", species)
        downloads = metadata["downloads"]
        if len(downloads) != len(DOWNLOAD_ROLES) or {item["name"] for item in downloads} != set(DOWNLOAD_ROLES):
            raise ValueError("public downloads differ from the approved artifact roles")
        for item in downloads:
            name = item["name"]
            media, schema_name, encoding = DOWNLOAD_ROLES[name]
            path = f"{prefix}/downloads/{name}"
            if item["path"] != path or (item["media_type"], item["schema"], item["encoding"]) != (media, schema_name, encoding):
                raise ValueError("public download has an unapproved role or path")
            file = role(path, media, schema_name, encoding)
            if name == "schema.json":
                if file.read_bytes() != schema_bytes:
                    raise ValueError("public schema copies disagree")
            elif name in citations:
                if file.read_bytes() != citations[name]:
                    raise ValueError("public citation export differs from its approved source")
            elif name == "citation.json":
                if file.read_bytes() != _json_bytes(manifest["citation"]) + b"\n":
                    raise ValueError("public citation JSON differs from its approved source")
            elif encoding is None:
                _scan_public_text(file)
                if name.endswith(".tsv"):
                    record_type = "summary" if name == "index.tsv" else "construct"
                    header = ("\t".join(schema["$defs"][record_type]["required"]) + "\n").encode()
                    with file.open("rb") as handle:
                        if handle.readline(65536) != header:
                            raise ValueError("public tabular export has an unapproved header")
    if set(manifest["files"]) != expected:
        raise ValueError("public manifest contains unreachable or unapproved artifact files")


def _bundle_files(bundle: Path, site: Path) -> tuple[dict, dict[str, Path], list[dict]]:
    # The producer validates schemas, canonical IDs, counts and complete-record equivalence.
    from agdesign2.public_api_data import validate_public_api_bundle

    if _file(bundle, "manifest.json").stat().st_size > 1024 * 1024:
        raise ValueError("root manifest exceeds the service's 1 MiB startup limit")
    root = validate_public_api_bundle(bundle)
    _check_root = _read_json(_file(bundle, "manifest.json"))
    if root != _check_root or root.get("schema_version") != 1:
        raise ValueError("unsupported or changed public manifest")
    sidecar = _file(bundle, "manifest.sha256")
    if sidecar.read_text(encoding="ascii") != f"{_hash_file(bundle / 'manifest.json')}  manifest.json\n":
        raise ValueError("root manifest checksum mismatch")
    releases = root.get("releases")
    if not isinstance(releases, list) or not 1 <= len(releases) <= 2:
        raise ValueError("bundle must retain current and at most one previous release")
    ids = [item.get("release_id") for item in releases]
    if len(set(ids)) != len(ids) or any(not _SHA.fullmatch(str(item)) for item in ids) or root.get("current_release_id") != ids[0]:
        raise ValueError("invalid retained release identifiers")
    retired = root.get("retired_releases", [])
    if not isinstance(retired, list) or len(retired) > 1000:
        raise ValueError("invalid retired release history")
    retired_ids = set()
    for item in retired:
        if not isinstance(item, dict) or set(item) != {"release_id", "retired_by_release_id", "archive_url"}:
            raise ValueError("invalid retired release tombstone")
        retired_id = item["release_id"]
        if not _SHA.fullmatch(str(retired_id)) or not _SHA.fullmatch(str(item["retired_by_release_id"])) or retired_id in ids or retired_id in retired_ids or retired_id == item["retired_by_release_id"]:
            raise ValueError("invalid retired release identifiers")
        retired_ids.add(retired_id)
        if item["archive_url"] is not None:
            url = item["archive_url"]
            if not isinstance(url, str) or any(ord(char) <= 32 for char in url):
                raise ValueError("unsafe retired archive URL")
            parsed = urlsplit(url)
            if parsed.scheme != "https" or not parsed.hostname or parsed.username or parsed.password:
                raise ValueError("unsafe retired archive URL")
    files = {"manifest.json": bundle / "manifest.json", "manifest.sha256": sidecar}
    summaries = []
    index_bytes = 0
    for pointer in releases:
        release_id = pointer["release_id"]
        manifest_name = f"releases/{release_id}/manifest.json"
        if pointer.get("manifest_path") != manifest_name:
            raise ValueError("release manifest is outside its immutable release")
        manifest_path = _file(bundle, manifest_name)
        if manifest_path.stat().st_size > 96 * 1024 * 1024:
            raise ValueError("release manifest exceeds the service's 96 MiB startup limit")
        if _hash_file(manifest_path) != pointer.get("sha256"):
            raise ValueError("release manifest checksum mismatch")
        manifest = _read_json(manifest_path)
        if manifest.get("release_id") != release_id or hashlib.sha256(_json_bytes(manifest.get("core"))).hexdigest() != release_id:
            raise ValueError("release identifier does not match its canonical core")
        if set(manifest.get("species", {})) != {"human", "mouse"}:
            raise ValueError("each release must contain human and mouse data")
        _check_public_file_roles(bundle, manifest)
        files[manifest_name] = manifest_path
        for name, entry in manifest["files"].items():
            if not name.startswith(f"releases/{release_id}/"):
                raise ValueError("release file is outside its immutable release")
            source = _file(bundle, name)
            _check_bytes(source, entry)
            maximum = _public_file_limit(name, entry["media_type"])
            if entry["bytes"] > maximum:
                label = PurePosixPath(name).name if "/downloads/" in name else entry["schema"]
                raise ValueError(f"public {label} is {entry['bytes']} bytes; limit is {maximum} bytes")
            files[name] = source
        species_summary = {}
        for species, taxon in (("human", 9606), ("mouse", 10090)):
            metadata = manifest["species"][species]
            index_bytes += manifest["files"][metadata["index_path"]]["bytes"]
            if index_bytes > 64 * 1024 * 1024:
                raise ValueError("retained compact indexes exceed the service's combined 64 MiB budget")
            if metadata.get("taxon_id") != taxon:
                raise ValueError("species taxon mismatch")
            if species == "mouse" and metadata.get("capabilities", {}).get("disease") is not False:
                raise ValueError("mouse disease capability must be explicitly false")
            records = metadata.get("records_download")
            if not isinstance(records, dict) or records.get("name") != "records.jsonl.gz":
                raise ValueError("complete-record fallback is required for each species")
            downloads = metadata.get("downloads", [])
            if records not in downloads:
                raise ValueError("complete-record fallback must be a listed download")
            for item in downloads:
                name = item.get("path", "")
                if name not in manifest["files"] or not name.startswith(f"releases/{release_id}/{species}/"):
                    raise ValueError("download is outside its species release")
                _check_bytes(files[name], item)
                url = urlsplit(item.get("url", ""))
                public_name = f"data/releases/{release_id}/{species}/{item['name']}"
                if url.scheme != "https" or not url.netloc or url.query or url.fragment or unquote(url.path) != "/" + public_name:
                    raise ValueError("download URL must identify its immutable static copy")
                _check_bytes(_file(site, public_name), item)
            species_summary[species] = {key: metadata[key] for key in ("taxon_id", "counts", "capabilities")}
            species_summary[species]["download_bytes"] = sum(item["bytes"] for item in downloads)
            species_summary[species]["download_files"] = len(downloads)
            if release_id == root["current_release_id"]:
                portal_path = "mouse/portal_metadata.json" if species == "mouse" else "portal_metadata.json"
                portal = _read_json(_file(site, portal_path))
                api = portal.get("public_api")
                if not isinstance(api, dict) or api.get("schema_version") != 1 or api.get("release_id") != release_id or api.get("state") not in {"packaged", "active"}:
                    raise ValueError("final portal metadata disagrees with the current API release")
                if api.get("manifest_url") != manifest.get("manifest_url") or portal.get("target_count") != metadata["counts"]["targets"]:
                    raise ValueError("final portal manifest or target count disagrees with the API release")
        static_manifest = _file(site, f"data/releases/{release_id}/manifest.json")
        if _hash_file(static_manifest) != pointer["sha256"]:
            raise ValueError("static release manifest differs from the service manifest")
        static_sidecar = _file(site, f"data/releases/{release_id}/manifest.sha256")
        if static_sidecar.read_text(encoding="ascii") != f"{pointer['sha256']}  manifest.json\n":
            raise ValueError("static release manifest checksum mismatch")
        summaries.append({"release_id": release_id, "manifest_sha256": pointer["sha256"], "species": species_summary, "static_metadata_bytes": static_manifest.stat().st_size + static_sidecar.stat().st_size})
    if _all_files(bundle) != set(files):
        raise ValueError("data bundle contains files outside its public manifest")
    return root, files, summaries


def _runtime_files(runtime: Path, logical_files: dict[str, Path], seal: dict) -> dict[str, Path]:
    """Read back only the sealed physical inventory generated from public inputs."""
    from agdesign2.public_api_data import _scan_public

    count_fields = {"transport_index_bytes", "physical_file_count", "physical_bytes", "packed_logical_file_count", "external_download_count"}
    if set(seal) != {"format", "root_manifest_sha256", "transport_index_sha256", *count_fields} or any(type(seal[name]) is not int or seal[name] < 0 for name in count_fields):
        raise ValueError("physical transport seal contains an unapproved field or count")
    if seal.get("format") != "openantigens-target-gzip-v1" or seal.get("root_manifest_sha256") != _hash_file(logical_files["manifest.json"]):
        raise ValueError("physical transport root seal disagrees with validated public data")
    index_file = _file(runtime, "transport-index.json")
    if index_file.stat().st_size > 16 * 1024 * 1024 or index_file.stat().st_size != seal["transport_index_bytes"] or _hash_file(index_file) != seal.get("transport_index_sha256"):
        raise ValueError("physical transport index seal mismatch")
    index = _read_json(index_file)
    if set(index) != {"format", "root_manifest_sha256", "targets"} or index["format"] != seal["format"] or index["root_manifest_sha256"] != seal["root_manifest_sha256"] or not isinstance(index["targets"], dict):
        raise ValueError("invalid physical transport index")
    _scan_public(index)
    targets = {name for name in logical_files if "/targets/" in name}
    if set(index["targets"]) != targets:
        raise ValueError("physical transport does not cover the validated target graph")
    expected = {"transport-index.json": index_file}
    tiny_references = {"schema.json", "citation.json", "openantigens.bib", "openantigens.ris"}
    for name, original in logical_files.items():
        standalone = name in {"manifest.json", "manifest.sha256"} or name.endswith("/manifest.json") or name.endswith("/index.json") or "/downloads/" in name and PurePosixPath(name).name in tiny_references
        if standalone:
            file = _file(runtime, name)
            if file.stat().st_size != original.stat().st_size or _hash_file(file) != _hash_file(original):
                raise ValueError("physical standalone file differs from its validated public source")
            expected[name] = file
    for target, parts in index["targets"].items():
        path = PurePosixPath(target)
        release_id, species, key = path.parts[1], path.parts[2], path.stem
        if not isinstance(parts, list) or not parts:
            raise ValueError("physical target has no transport parts")
        previous_max = None
        for number, descriptor in enumerate(parts):
            if not isinstance(descriptor, dict) or set(descriptor) != {"path", "sha256", "bytes", "decoded_bytes", "members", "min_path", "max_path"}:
                raise ValueError("invalid physical part descriptor")
            name = f"packs/{release_id}/{species}/{key}-{number}.json.gz"
            if descriptor["path"] != name or type(descriptor["decoded_bytes"]) is not int or not 0 < descriptor["decoded_bytes"] <= 16 * 1024 * 1024 or type(descriptor["members"]) is not int or descriptor["members"] < 1:
                raise ValueError("invalid physical target part range or size")
            low, high = descriptor["min_path"], descriptor["max_path"]
            if not isinstance(low, str) or not isinstance(high, str) or low not in logical_files or high not in logical_files or low > high or previous_max is not None and low <= previous_max:
                raise ValueError("physical target part ranges overlap or escape public inputs")
            previous_max = high
            file = _file(runtime, name)
            _check_bytes(file, descriptor)
            if file.stat().st_size > 17 * 1024 * 1024:
                raise ValueError("compressed physical part exceeds its service limit")
            expected[name] = file
    if _all_files(runtime) != set(expected) or len(expected) != seal.get("physical_file_count") or sum(path.stat().st_size for path in expected.values()) != seal.get("physical_bytes"):
        raise ValueError("physical runtime contains an unexpected file or changed inventory")
    return expected


def _query_files(runtime: Path, compiled: dict, schema_path: Path) -> tuple[dict[str, Path], dict]:
    from agdesign2.public_api_query import validate_query_runtime

    assets = validate_query_runtime(runtime, schema_path=schema_path, seal=compiled["seal"])
    if assets != compiled["static_assets"] or _all_files(runtime) != set(compiled["runtime_files"]):
        raise ValueError("query compiler inventory disagrees with its sealed runtime")
    files = {}
    for name, metadata in compiled["runtime_files"].items():
        file = _file(runtime, name)
        _check_bytes(file, metadata)
        files[name] = file
    return files, assets


def _asset_destination(site: Path, name: str) -> Path:
    _relative(name)
    if not re.fullmatch(r"api-runtime/[a-f0-9]{64}/(?:human|mouse)/[a-f0-9]{64}-[01]-[a-f0-9]{64}\.oap", name):
        raise ValueError("unapproved immutable query asset path")
    target = site / name
    for parent in (target, *target.parents):
        if parent == site:
            break
        if parent.is_symlink():
            raise ValueError("query asset destinations cannot contain symlinks")
    if not target.resolve().is_relative_to(site):
        raise ValueError("query asset destination escapes the static root")
    return target


def _stage_query_assets(site: Path, source: Path, inventory: dict, *, dry_run: bool) -> None:
    for name, metadata in inventory.items():
        _check_bytes(_file(source, name), metadata)
        target = _asset_destination(site, name)
        if target.exists():
            _check_bytes(_file(site, name), metadata)
    if dry_run:
        return
    for name, metadata in inventory.items():
        target = _asset_destination(site, name)
        if target.exists():
            _check_bytes(_file(site, name), metadata)
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        temporary = None
        try:
            with tempfile.NamedTemporaryFile(dir=target.parent, prefix=".query-asset-", delete=False) as destination:
                temporary = Path(destination.name)
                with _file(source, name).open("rb") as original:
                    shutil.copyfileobj(original, destination)
            _check_bytes(temporary, metadata)
            temporary.chmod(0o644)
            try:
                os.link(temporary, target)
            except FileExistsError:
                _check_bytes(_file(site, name), metadata)
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)


@contextmanager
def _reviewed_query_data(public_site_dir: Path | str, artifact_record: Path | str, *, expected_zip_sha256: str | None = None, schema_path: Path | None = None):
    """Verify a reviewed archive and current local/static boundaries without old code execution."""
    from agdesign2.public_api_query import validate_query_runtime

    site = _root(public_site_dir)
    record_path = Path(artifact_record).expanduser().resolve()
    record = _read_json(record_path)
    if expected_zip_sha256 is not None and (not _SHA.fullmatch(expected_zip_sha256) or record.get("sha256") != expected_zip_sha256):
        raise ValueError("reviewed archive does not match the independently supplied ZIP SHA")
    archive_name = record.get("artifact", "")
    if not isinstance(archive_name, str) or PurePosixPath(archive_name).name != archive_name or not archive_name.endswith(".zip"):
        raise ValueError("review record must name one local app ZIP")
    archive_path = _file(record_path.parent, archive_name)
    if _hash_file(archive_path) != record.get("sha256") or archive_path.stat().st_size != record.get("bytes"):
        raise ValueError("reviewed app ZIP checksum mismatch")
    schema_path = schema_path or Path(__file__).resolve().parents[2] / "services/public-api/src/query-schema.json"
    if archive_path.stat().st_size >= HARD_ZIP_LIMIT or archive_path.stat().st_size > 80_000_000:
        raise ValueError("reviewed ZIP exceeds the default publication budget")
    with tempfile.TemporaryDirectory(prefix="oa-query-publish-check-") as temporary:
        folder = Path(temporary).resolve()
        runtime = folder / "data"
        runtime.mkdir()
        with zipfile.ZipFile(archive_path) as archive:
            names = archive.namelist()
            if len(set(names)) != len(names) or set(names) != set(record.get("files", {})):
                raise ValueError("reviewed ZIP inventory mismatch")
            if sum(info.file_size for info in archive.infolist()) > 400_000_000:
                raise ValueError("reviewed ZIP exceeds the publication expanded budget")
            build = None
            for info in archive.infolist():
                name = _relative(info.filename)
                if not stat.S_ISREG(info.external_attr >> 16):
                    raise ValueError("reviewed ZIP contains a non-regular entry")
                expected = record["files"][name]
                if type(expected.get("bytes")) is not int or expected["bytes"] < 0 or info.file_size != expected["bytes"]:
                    raise ValueError("reviewed ZIP declared member size mismatch")
                digest = hashlib.sha256()
                size = 0
                selected = name.startswith("data/") or name == "service-build.json"
                destination = folder / name
                if selected:
                    destination.parent.mkdir(parents=True, exist_ok=True)
                with archive.open(info) as original:
                    target = destination.open("xb") if selected else None
                    try:
                        for chunk in iter(lambda: original.read(1024 * 1024), b""):
                            digest.update(chunk)
                            size += len(chunk)
                            if target is not None:
                                target.write(chunk)
                    finally:
                        if target is not None:
                            target.close()
                if {"sha256": digest.hexdigest(), "bytes": size} != expected:
                    raise ValueError("reviewed ZIP member checksum mismatch")
            build = _read_json(folder / "service-build.json")
            code_names = set(build.get("files", {}))
            if not code_names.issubset(set(SOURCE_FILES) | set(LICENSE_FILES)) or {name for name in names if not name.startswith("data/")} != code_names | {"service-build.json"}:
                raise ValueError("reviewed ZIP code inventory is not closed")
            if any(record["files"].get(name) != metadata for name, metadata in build["files"].items()):
                raise ValueError("reviewed source manifest disagrees with ZIP receipt")
        inventory = validate_query_runtime(runtime, schema_path=schema_path, seal=build["data_transport"])
        if inventory != record.get("query_static_assets"):
            raise ValueError("review record assets differ from the sealed app catalog")
        for name, metadata in inventory.items():
            _asset_destination(site, name)
            _check_bytes(_file(site, name), metadata)
        query_root = _read_json(runtime / "query-root.json")
        for pointer in query_root["releases"]:
            release = _read_json(_file(runtime, pointer["path"]))
            release_id = release["release_id"]
            prefix = f"data/releases/{release_id}"
            manifest = _file(site, prefix + "/manifest.json")
            manifest_sha = release["source_manifest"]["sha256"]
            if manifest.stat().st_size > 96 * 1024 * 1024 or _hash_file(manifest) != manifest_sha:
                raise ValueError("static source manifest differs from reviewed query catalog")
            inventory[prefix + "/manifest.json"] = {"bytes": manifest.stat().st_size, "sha256": manifest_sha}
            checksum = _file(site, prefix + "/manifest.sha256")
            expected_checksum = f"{manifest_sha}  manifest.json\n".encode("ascii")
            if checksum.read_bytes() != expected_checksum:
                raise ValueError("static source manifest sidecar differs from reviewed catalog")
            inventory[prefix + "/manifest.sha256"] = {"bytes": len(expected_checksum), "sha256": hashlib.sha256(expected_checksum).hexdigest()}
            for species in ("human", "mouse"):
                downloads = release["species"][species]["downloads"]
                if len(downloads) != len(DOWNLOAD_ROLES) or {item["name"] for item in downloads} != set(DOWNLOAD_ROLES):
                    raise ValueError("reviewed static exports have an unapproved role")
                for item in downloads:
                    name = item["name"]
                    if (item["media_type"], item["schema"], item["encoding"]) != DOWNLOAD_ROLES[name]:
                        raise ValueError("reviewed static export media/schema/encoding mismatch")
                    relative = f"{prefix}/{species}/{name}"
                    _check_bytes(_file(site, relative), item)
                    inventory[relative] = {"bytes": item["bytes"], "sha256": item["sha256"]}
        yield runtime, build, record, inventory, archive_path


def write_query_asset_list(public_site_dir: Path | str, artifact_record: Path | str, output_list: Path | str) -> int:
    """Authorize additive transfer only from reviewed, currently validated catalogs."""
    with _reviewed_query_data(public_site_dir, artifact_record) as (_, _, _, inventory, _):
        with Path(output_list).open("xb") as handle:
            handle.write(b"".join(name.encode("ascii") + b"\0" for name in sorted(inventory)))
        return len(inventory)


def build_api_service_release(
    service_dir: Path | str, bundle_dir: Path | str, public_site_dir: Path | str,
    output_zip: Path | str, *, repo_root: Path | str,
    dry_run: bool = False, max_zip_bytes: int = HARD_ZIP_LIMIT,
    max_expanded_bytes: int = 500_000_000, reserve_fraction: float = 0.20,
    dependency_audit_report: Path | str | None = None,
) -> dict:
    """Validate and package immutable data plus an exact source-file allowlist.

    Dry runs measure the real compressed archive in temporary storage and leave
    inputs and requested output paths unchanged. All output writes are exclusive.
    Limits are local packaging budgets; they do not assert a host's capacity.
    """
    if type(max_zip_bytes) is not int or not 0 < max_zip_bytes <= HARD_ZIP_LIMIT or type(max_expanded_bytes) is not int or max_expanded_bytes <= 0:
        raise ValueError("invalid local size budget")
    if not isinstance(reserve_fraction, (float, int)) or not 0 <= reserve_fraction < 1:
        raise ValueError("reserve fraction must be between zero and one")
    service, bundle, site, repo = map(_root, (service_dir, bundle_dir, public_site_dir, repo_root))
    _canonical_repo(repo)
    output = Path(output_zip).expanduser().absolute()
    if output.suffix != ".zip" or output.is_symlink():
        raise ValueError("output must be a new .zip file")
    artifact_path = output.with_suffix(".artifact.json")
    if output.exists() or artifact_path.exists() or artifact_path.is_symlink():
        raise ValueError("release output already exists; choose a new output name")
    for root_dir in (service, bundle, site):
        if output.resolve().is_relative_to(root_dir) or artifact_path.resolve().is_relative_to(root_dir):
            raise ValueError("output must be outside all packaging inputs")
    root, data_files, releases = _bundle_files(bundle, site)
    from agdesign2.public_api_query import _build_query_runtime

    logical_file_count = len(data_files)
    logical_bytes = sum(path.stat().st_size for path in data_files.values())
    with tempfile.TemporaryDirectory(prefix="oa-service-runtime-") as temporary:
        runtime = Path(temporary).resolve() / "data"
        assets_root = Path(temporary).resolve() / "assets"
        schema_path = _file(service, "src/query-schema.json")
        compiled = _build_query_runtime(bundle, runtime, assets_root, schema_path=schema_path, validated_root=root)
        data_files, asset_inventory = _query_files(runtime, compiled, schema_path)
        transport = compiled["seal"]
        _stage_query_assets(site, assets_root, asset_inventory, dry_run=True)
        return _write_service_archive(service, repo, output, site, data_files, transport, root, releases, logical_file_count, logical_bytes, asset_inventory, assets_root=assets_root, dry_run=dry_run, max_zip_bytes=max_zip_bytes, max_expanded_bytes=max_expanded_bytes, reserve_fraction=reserve_fraction, dependency_audit_report=dependency_audit_report)


def _write_service_archive(service, repo, output, site, data_files, transport, root, releases, logical_file_count, logical_bytes, asset_inventory, *, assets_root, dry_run, max_zip_bytes, max_expanded_bytes, reserve_fraction, dependency_audit_report, lineage=None):
    artifact_path = output.with_suffix(".artifact.json")
    sources = {name: _file(service, name) for name in SOURCE_FILES}
    sources.update({name: _file(repo, name) for name in LICENSE_FILES})
    package = _read_json(sources["package.json"])
    lock = _read_json(sources["package-lock.json"])
    if package.get("engines", {}).get("node") != ">=22 <23" or package.get("scripts", {}).get("start") != "node src/server.js":
        raise ValueError("service package must declare the supported Node 22 start command")
    if lock.get("lockfileVersion") != 3 or lock.get("packages", {}).get("", {}).get("dependencies") != package.get("dependencies"):
        raise ValueError("dependency lock does not match package.json")
    if not package.get("dependencies") or any(not re.fullmatch(r"[0-9]+\.[0-9]+\.[0-9]+(?:-[A-Za-z0-9.-]+)?", str(version)) for version in package["dependencies"].values()):
        raise ValueError("runtime dependencies must use exact pinned versions")
    for name, entry in lock["packages"].items():
        if not name:
            continue
        resolved = urlsplit(entry.get("resolved", ""))
        if not name.startswith("node_modules/") or entry.get("link") or resolved.scheme != "https" or resolved.netloc != "registry.npmjs.org" or resolved.query or resolved.fragment or not entry.get("integrity", "").startswith(("sha512-", "sha256-")):
            raise ValueError("lockfile dependencies must use integrity-pinned public npm registry URLs")
    code_files = {name: {"sha256": _hash_file(path), "bytes": path.stat().st_size} for name, path in sources.items()}
    build_id = hashlib.sha256(_json_bytes(code_files)).hexdigest()
    build = {"schema_version": 1, "service_build_id": build_id, "runtime": "node22", "files": code_files, "data_transport": transport}
    generated = {"service-build.json": _json_bytes(build) + b"\n"}
    inputs = {**sources, **{f"data/{name}": path for name, path in data_files.items()}}
    expected = {name: {"sha256": _hash_file(path), "bytes": path.stat().st_size} for name, path in inputs.items()}
    expected.update({name: {"sha256": hashlib.sha256(value).hexdigest(), "bytes": len(value)} for name, value in generated.items()})
    expanded = sum(entry["bytes"] for entry in expected.values())
    static_bytes = sum(item["static_metadata_bytes"] + sum(value["download_bytes"] for value in item["species"].values()) for item in releases)
    static_files = sum(2 + sum(value["download_files"] for value in item["species"].values()) for item in releases)
    if expanded > max_expanded_bytes * (1 - reserve_fraction):
        raise ValueError(f"expanded artifact is {expanded} bytes; local budget with reserved headroom is {int(max_expanded_bytes * (1 - reserve_fraction))} bytes")
    audit = {"status": "not_run", "command": "npm audit --omit=dev --audit-level=high", "report": None}
    if dependency_audit_report is not None:
        audit_path = Path(dependency_audit_report)
        if audit_path.is_symlink() or not audit_path.is_file():
            raise ValueError("dependency audit report must be a regular file")
        report = _read_json(audit_path)
        audit = {"status": "provided_unverified", "command": audit["command"], "report": {"name": audit_path.name, "sha256": _hash_file(audit_path), "bytes": audit_path.stat().st_size}, "reported_vulnerabilities": report.get("metadata", {}).get("vulnerabilities"), "note": "Report is recorded; command execution, freshness and the installed dependency tree require review."}
    with tempfile.SpooledTemporaryFile(max_size=8 * 1024 * 1024) as spool:
        with zipfile.ZipFile(spool, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
            for name in sorted(expected):
                info = zipfile.ZipInfo(_relative(name), date_time=(1980, 1, 1, 0, 0, 0))
                info.create_system = 3
                info.external_attr = (stat.S_IFREG | 0o644) << 16
                info.compress_type = zipfile.ZIP_DEFLATED
                info._compresslevel = 9
                info.file_size = expected[name]["bytes"]
                digest = hashlib.sha256()
                count = 0
                with archive.open(info, "w") as destination:
                    if name in generated:
                        chunks = (generated[name],)
                        for chunk in chunks:
                            destination.write(chunk)
                            digest.update(chunk)
                            count += len(chunk)
                    else:
                        with inputs[name].open("rb") as source:
                            for chunk in iter(lambda: source.read(1024 * 1024), b""):
                                destination.write(chunk)
                                digest.update(chunk)
                                count += len(chunk)
                if {"sha256": digest.hexdigest(), "bytes": count} != expected[name]:
                    raise ValueError("release input changed while packaging")
        size = spool.tell()
        if size >= HARD_ZIP_LIMIT or size > max_zip_bytes * (1 - reserve_fraction):
            raise ValueError(f"compressed artifact is {size} bytes; local budget with reserved headroom is {int(max_zip_bytes * (1 - reserve_fraction))} bytes and hard ceiling is below {HARD_ZIP_LIMIT} bytes")
        spool.seek(0)
        archive_hash = hashlib.file_digest(spool, "sha256").hexdigest()
        artifact = {
            "schema_version": 1, "artifact": output.name, "sha256": archive_hash,
            "bytes": size, "expanded_bytes": expanded, "file_count": len(expected),
            "physical_data_file_count": len(data_files), "physical_data_bytes": sum(path.stat().st_size for path in data_files.values()),
            "logical_data_file_count": logical_file_count, "logical_data_bytes": logical_bytes,
            "data_transport": transport,
            "query_static_assets": asset_inventory,
            "query_static_asset_file_count": len(asset_inventory),
            "query_static_asset_bytes": sum(item["bytes"] for item in asset_inventory.values()),
            "external_bulk_download_verification": {"staging": "complete_bytes_and_exports_verified", "runtime": "metadata_only_external_static_files"},
            "zip64_required": len(expected) > zipfile.ZIP_FILECOUNT_LIMIT or any(entry["bytes"] * 1.05 > zipfile.ZIP64_LIMIT for entry in expected.values()),
            "retained_static_download_bytes": static_bytes, "retained_static_download_file_count": static_files,
            "minimum_app_plus_retained_static_bytes": expanded + static_bytes + sum(item["bytes"] for item in asset_inventory.values()),
            "service_build_id": build_id, "current_release_id": root["current_release_id"],
            "releases": releases, "runtime": "node22", "files": expected,
            "retired_releases": root.get("retired_releases", []),
            "local_budgets": {"zip_bytes": max_zip_bytes, "expanded_bytes": max_expanded_bytes, "reserve_fraction": reserve_fraction, "zip_remaining_bytes": max_zip_bytes - size, "expanded_remaining_bytes": max_expanded_bytes - expanded, "host_capacity_verified": False},
            "dependency_audit": audit,
            "gates": {"public_bundle_and_static_downloads": "passed", "exact_source_allowlist": "passed", "local_package_size_with_headroom": "passed", "physical_transport_integrity": "passed", "dependency_audit": "not_verified", "service_contract_tests": "not_verified", **{name: "not_verified" for name in HOST_GATES}},
            "publication_state": "prepared_local_artifact", "dry_run": dry_run,
        }
        if lineage is not None:
            actual_data = {name[5:]: metadata for name, metadata in expected.items() if name.startswith("data/")}
            if actual_data != lineage["reused_data_files"]:
                raise ValueError("code-only repack changed the reviewed data inventory or bytes")
            artifact["validation_lineage"] = {key: value for key, value in lineage.items() if key != "reused_data_files"}
            artifact["gates"]["public_bundle_and_static_downloads"] = "inherited_science_current_query_and_static_checks_passed"
        if not dry_run:
            if assets_root is not None:
                _stage_query_assets(site, assets_root, asset_inventory, dry_run=False)
            output.parent.mkdir(parents=True, exist_ok=True)
            created = []
            try:
                with output.open("xb") as destination:
                    created.append(output)
                    spool.seek(0)
                    shutil.copyfileobj(spool, destination)
                with artifact_path.open("xb") as destination:
                    created.append(artifact_path)
                    destination.write(_json_bytes(artifact) + b"\n")
            except BaseException:
                for path in created:
                    path.unlink(missing_ok=True)
                raise
    return artifact


def repack_api_service_release(
    artifact_record: Path | str, expected_zip_sha256: str,
    service_dir: Path | str, public_site_dir: Path | str, output_zip: Path | str,
    *, repo_root: Path | str, dry_run: bool = False,
    max_zip_bytes: int = HARD_ZIP_LIMIT, max_expanded_bytes: int = 500_000_000,
    reserve_fraction: float = 0.20, dependency_audit_report: Path | str | None = None,
) -> dict:
    """Replace code only, inheriting a separately pinned reviewed scientific receipt.

    This route is ineligible for projection, schema, transport or data changes,
    or data-boundary fixes not covered by the current query validator.
    """
    if not isinstance(expected_zip_sha256, str) or not _SHA.fullmatch(expected_zip_sha256):
        raise ValueError("code-only repack requires an independently supplied reviewed ZIP SHA")
    if type(max_zip_bytes) is not int or not 0 < max_zip_bytes <= HARD_ZIP_LIMIT or type(max_expanded_bytes) is not int or max_expanded_bytes <= 0 or not isinstance(reserve_fraction, (float, int)) or not 0 <= reserve_fraction < 1:
        raise ValueError("invalid local repack budget")
    service, site, repo = map(_root, (service_dir, public_site_dir, repo_root))
    _canonical_repo(repo)
    output = Path(output_zip).expanduser().absolute()
    receipt = output.with_suffix(".artifact.json")
    record_path = Path(artifact_record).expanduser().resolve()
    if output.suffix != ".zip" or output.is_symlink() or output.exists() or receipt.exists() or receipt.is_symlink():
        raise ValueError("code-only output must be a new .zip and artifact record")
    if any(output.resolve().is_relative_to(root) or receipt.resolve().is_relative_to(root) for root in (service, site)) or output.resolve() == record_path or receipt.resolve() == record_path:
        raise ValueError("code-only output must be outside its inputs")
    schema = _file(service, "src/query-schema.json")
    with _reviewed_query_data(site, record_path, expected_zip_sha256=expected_zip_sha256, schema_path=schema) as (runtime, old_build, old_record, static_files, archive):
        transport = old_build["data_transport"]
        if transport.get("format") != "openantigens-own-static-query-v1" or transport != old_record.get("data_transport"):
            raise ValueError("code-only repack requires the reviewed own-static data format and identical seal")
        query_root = _read_json(runtime / "query-root.json")
        if query_root["current_release_id"] != old_record.get("current_release_id") or query_root["retired_releases"] != old_record.get("retired_releases", []):
            raise ValueError("reviewed data release identity differs from its receipt")
        data_files = {name: _file(runtime, name) for name in _all_files(runtime)}
        old_data = {name[5:]: metadata for name, metadata in old_record["files"].items() if name.startswith("data/")}
        if set(data_files) != set(old_data):
            raise ValueError("code-only data path set changed")
        for name, metadata in old_data.items():
            _check_bytes(data_files[name], metadata)
        releases = []
        for pointer in query_root["releases"]:
            value = _read_json(_file(runtime, pointer["path"]))
            species = {}
            for name in ("human", "mouse"):
                entry = value["species"][name]
                species[name] = {key: entry[key] for key in ("taxon_id", "counts", "capabilities")}
                species[name].update(download_bytes=sum(item["bytes"] for item in entry["downloads"]), download_files=len(entry["downloads"]))
            prefix = f"data/releases/{value['release_id']}"
            releases.append({"release_id": value["release_id"], "manifest_sha256": value["source_manifest"]["sha256"], "species": species, "static_metadata_bytes": static_files[prefix + "/manifest.json"]["bytes"] + static_files[prefix + "/manifest.sha256"]["bytes"]})
        for name in ("logical_data_file_count", "logical_data_bytes"):
            if type(old_record.get(name)) is not int or old_record[name] < 1:
                raise ValueError("reviewed scientific receipt lacks inherited dataset measurements")
        lineage = {
            "mode": "code_only_repack", "parent_zip_sha256": expected_zip_sha256,
            "parent_artifact_record_sha256": _hash_file(record_path),
            "scientific_validation": "inherited_from_independently_pinned_full_validation_receipt_not_rerun",
            "current_query_schema_retirement_and_static_checks": "rerun_passed",
            "data_inventory_and_bytes": "unchanged", "data_transport_and_scientific_ids": "unchanged",
            "eligibility": "runtime_code_only_no_projection_schema_transport_data_or_uncovered_data_boundary_change",
            "source_measurements": "inherited_parent_receipt", "reused_data_files": old_data,
        }
        return _write_service_archive(service, repo, output, site, data_files, transport, query_root, releases, old_record["logical_data_file_count"], old_record["logical_data_bytes"], old_record["query_static_assets"], assets_root=None, dry_run=dry_run, max_zip_bytes=max_zip_bytes, max_expanded_bytes=max_expanded_bytes, reserve_fraction=reserve_fraction, dependency_audit_report=dependency_audit_report, lineage=lineage)
