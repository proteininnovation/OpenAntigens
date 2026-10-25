"""Build a bounded, deterministic public-service ZIP from validated public data.

This module packages an existing data release. It never projects new scientific
records, contacts a host, or publishes an artifact.
"""

from __future__ import annotations

import hashlib
import json
import re
import shutil
import stat
import subprocess
import tempfile
import zipfile
from pathlib import Path, PurePosixPath
from typing import Any
from urllib.parse import unquote, urlsplit

HARD_ZIP_LIMIT = 100_000_000
SOURCE_FILES = (
    "package.json", "package-lock.json", "src/server.js", "src/mcp.js",
    "src/provider.js", "src/contracts.js", "src/openapi.js",
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
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError("manifest must be a JSON object")
    return value


def _check_bytes(path: Path, entry: dict) -> None:
    if type(entry.get("bytes")) is not int or entry["bytes"] < 0 or not _SHA.fullmatch(str(entry.get("sha256", ""))):
        raise ValueError("invalid file checksum metadata")
    if path.stat().st_size != entry["bytes"] or _hash_file(path) != entry["sha256"]:
        raise ValueError("public file size or checksum mismatch")


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
        if manifest_path.stat().st_size > 64 * 1024 * 1024:
            raise ValueError("release manifest exceeds the service's 64 MiB startup limit")
        if _hash_file(manifest_path) != pointer.get("sha256"):
            raise ValueError("release manifest checksum mismatch")
        manifest = _read_json(manifest_path)
        if manifest.get("release_id") != release_id or hashlib.sha256(_json_bytes(manifest.get("core"))).hexdigest() != release_id:
            raise ValueError("release identifier does not match its canonical core")
        if set(manifest.get("species", {})) != {"human", "mouse"}:
            raise ValueError("each release must contain human and mouse data")
        files[manifest_name] = manifest_path
        for name, entry in manifest["files"].items():
            if not name.startswith(f"releases/{release_id}/"):
                raise ValueError("release file is outside its immutable release")
            source = _file(bundle, name)
            _check_bytes(source, entry)
            if name.endswith("/index.json"):
                maximum = 32 * 1024 * 1024
            elif "/targets/" in name:
                maximum = 8 * 1024 * 1024
            elif entry["media_type"] == "application/schema+json" or entry["media_type"] == "application/json" and "/downloads/" not in name:
                maximum = 1024 * 1024
            else:
                maximum = 128 * 1024 * 1024
            if entry["bytes"] > maximum:
                raise ValueError("public file exceeds the service's startup read limit")
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
    origin = subprocess.run(
        ["git", "-C", str(repo), "remote", "get-url", "origin"],
        capture_output=True, text=True, check=False,
    )
    if origin.returncode or origin.stdout.strip() not in {
        "https://github.com/proteininnovation/OpenAntigens", "https://github.com/proteininnovation/OpenAntigens.git",
        "git@github.com:proteininnovation/OpenAntigens.git", "ssh://git@github.com/proteininnovation/OpenAntigens.git",
    }:
        raise ValueError("origin must resolve to proteininnovation/OpenAntigens")
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
    build = {"schema_version": 1, "service_build_id": build_id, "runtime": "node22", "files": code_files}
    generated = {"service-build.json": _json_bytes(build) + b"\n"}
    inputs = {**sources, **{f"data/{name}": path for name, path in data_files.items()}}
    expected = {name: {"sha256": _hash_file(path), "bytes": path.stat().st_size} for name, path in inputs.items()}
    expected.update({name: {"sha256": hashlib.sha256(value).hexdigest(), "bytes": len(value)} for name, value in generated.items()})
    expanded = sum(entry["bytes"] for entry in expected.values())
    static_bytes = sum(item["static_metadata_bytes"] + sum(value["download_bytes"] for value in item["species"].values()) for item in releases)
    static_files = sum(2 + sum(value["download_files"] for value in item["species"].values()) for item in releases)
    if expanded > max_expanded_bytes * (1 - reserve_fraction):
        raise ValueError("expanded artifact exceeds its local budget with reserved headroom")
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
            raise ValueError("compressed artifact exceeds its local budget with reserved headroom")
        spool.seek(0)
        archive_hash = hashlib.file_digest(spool, "sha256").hexdigest()
        artifact = {
            "schema_version": 1, "artifact": output.name, "sha256": archive_hash,
            "bytes": size, "expanded_bytes": expanded, "file_count": len(expected),
            "zip64_required": len(expected) > zipfile.ZIP_FILECOUNT_LIMIT or any(entry["bytes"] * 1.05 > zipfile.ZIP64_LIMIT for entry in expected.values()),
            "retained_static_download_bytes": static_bytes, "retained_static_download_file_count": static_files,
            "minimum_app_plus_retained_static_bytes": expanded + static_bytes,
            "service_build_id": build_id, "current_release_id": root["current_release_id"],
            "releases": releases, "runtime": "node22", "files": expected,
            "retired_releases": root.get("retired_releases", []),
            "local_budgets": {"zip_bytes": max_zip_bytes, "expanded_bytes": max_expanded_bytes, "reserve_fraction": reserve_fraction, "zip_remaining_bytes": max_zip_bytes - size, "expanded_remaining_bytes": max_expanded_bytes - expanded, "host_capacity_verified": False},
            "dependency_audit": audit,
            "gates": {"public_bundle_and_static_downloads": "passed", "exact_source_allowlist": "passed", "local_package_size_with_headroom": "passed", "dependency_audit": "not_verified", "service_contract_tests": "not_verified", **{name: "not_verified" for name in HOST_GATES}},
            "publication_state": "prepared_local_artifact", "dry_run": dry_run,
        }
        if not dry_run:
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
