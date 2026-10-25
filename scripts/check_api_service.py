#!/usr/bin/env python3
"""Bounded HTTP and independent MCP wire smoke checks for an explicit release."""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import tempfile
from pathlib import Path
from urllib.error import HTTPError
from urllib.parse import urlsplit
from urllib.request import HTTPRedirectHandler, Request, build_opener

SHA = re.compile(r"[0-9a-f]{64}\Z")
TOOLS = {"get_release", "search_targets", "get_target", "list_constructs", "get_sequence", "get_downloads"}
MAX_RESPONSE_BYTES = 1024 * 1024


class NoRedirect(HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        return None


def _base_url(value: str) -> str:
    if any(ord(char) <= 32 for char in value):
        raise ValueError("base URL contains whitespace or control characters")
    parsed = urlsplit(value)
    if parsed.scheme not in {"http", "https"} or not parsed.hostname or parsed.username or parsed.password or parsed.query or parsed.fragment or parsed.path not in {"", "/"}:
        raise ValueError("base URL must be an HTTP(S) origin without credentials, a path, query or fragment")
    if parsed.scheme == "http" and parsed.hostname not in {"localhost", "127.0.0.1", "::1"}:
        raise ValueError("remote smoke checks require HTTPS")
    return value.rstrip("/")


def _request(opener, base: str, path: str, *, payload: dict | None = None, protocol: str | None = None) -> tuple[int, object, object]:
    headers = {"Accept": "application/json, text/event-stream" if payload is not None else "application/json"}
    body = None
    if payload is not None:
        headers["Content-Type"] = "application/json"
        if protocol:
            headers["MCP-Protocol-Version"] = protocol
            headers["MCP-Method"] = payload["method"]
            if payload.get("params", {}).get("name"):
                headers["MCP-Name"] = payload["params"]["name"]
        body = json.dumps(payload, separators=(",", ":")).encode("utf-8")
    request = Request(base + path, data=body, headers=headers, method="POST" if body is not None else "GET")
    try:
        response = opener.open(request, timeout=5)
    except HTTPError as exc:
        response = exc
    with response:
        data = response.read(MAX_RESPONSE_BYTES + 1)
        if len(data) > MAX_RESPONSE_BYTES:
            raise ValueError("smoke response exceeds the byte limit")
        if response.headers.get_content_type() != "application/json":
            raise ValueError("smoke response must be JSON")
        return response.status, json.loads(data), response.headers


def _envelope(value: object, release_id: str, species: str | None = None) -> dict:
    if not isinstance(value, dict) or value.get("release_id") != release_id or value.get("schema_version") != 1:
        raise ValueError("response release/schema mismatch")
    if species is not None:
        taxon = 9606 if species == "human" else 10090
        if value.get("species") != species or value.get("taxon_id") != taxon:
            raise ValueError("response species/taxon mismatch")
        if not isinstance(value.get("data"), list):
            raise ValueError("species collection must contain an array")
        for row in value["data"]:
            if not isinstance(row, dict) or row.get("species", species) != species or row.get("taxon_id", taxon) != taxon:
                raise ValueError("collection contains a record from another species")
    return value


def _discovery(value: dict, expected_current: str, maintenance: bool) -> list[str]:
    if maintenance:
        if value.get("discovery_state") != "publication_in_progress" or "current_release_id" in value or "releases" in value:
            raise ValueError("pinned metadata leaks current discovery during maintenance")
        return []
    if value.get("current_release_id") != expected_current:
        raise ValueError("literal metadata identifies another current release")
    retained = value.get("releases")
    ids = [entry.get("release_id") for entry in retained] if isinstance(retained, list) and all(isinstance(entry, dict) for entry in retained) else []
    if not ids or len(ids) > 2 or len(set(ids)) != len(ids) or ids[0] != expected_current or any(not SHA.fullmatch(str(item)) for item in ids):
        raise ValueError("retained release discovery is invalid")
    return ids


def check_api_service(base_url: str, release_id: str, *, expect_maintenance: bool = False, retained_release_ids: list[str] | None = None) -> dict:
    """Perform at most twenty-two sequential reads; send no credentials or data."""
    base = _base_url(base_url)
    if not SHA.fullmatch(release_id) or any(not SHA.fullmatch(value) for value in retained_release_ids or []):
        raise ValueError("expected release IDs must be lowercase SHA-256 values")
    known_ids = list(dict.fromkeys([release_id, *(retained_release_ids or [])]))
    if len(known_ids) > 2:
        raise ValueError("smoke checks accept current plus at most one caller-known retained release")
    opener = build_opener(NoRedirect)
    for path, state in (("/healthz", "ok"), ("/readyz", "ready")):
        status, value, headers = _request(opener, base, path)
        if status != 200 or not isinstance(value, dict) or value.get("status") != state or headers.get("Cache-Control") != "no-store":
            raise ValueError("service health/readiness check failed")
    status, current, headers = _request(opener, base, "/v1/releases/current")
    if expect_maintenance:
        if status != 503 or not isinstance(current, dict) or current.get("error", {}).get("code") != "publication_in_progress" or headers.get("Cache-Control") != "no-store" or not headers.get("Retry-After"):
            raise ValueError("current discovery is not in the expected maintenance state")
    elif status != 200 or headers.get("Cache-Control") != "no-store":
        raise ValueError("current discovery check failed")
    else:
        _envelope(current, release_id)
    status, metadata, headers = _request(opener, base, f"/v1/releases/{release_id}")
    if status != 200 or headers.get("Cache-Control") != "no-store":
        raise ValueError("literal release metadata check failed")
    _envelope(metadata, release_id)
    discovered_ids = _discovery(metadata, release_id, expect_maintenance)
    ids = known_ids if expect_maintenance else discovered_ids
    if not expect_maintenance and retained_release_ids is not None and set(ids) != set(known_ids):
        raise ValueError("retained release IDs do not match the expected release set")
    released = metadata.get("data")
    if not isinstance(released, dict) or set(released.get("species", {})) != {"human", "mouse"}:
        raise ValueError("release metadata must describe both species")
    for species, taxon in (("human", 9606), ("mouse", 10090)):
        entry = released["species"][species]
        if not isinstance(entry, dict) or entry.get("taxon_id") != taxon or species == "mouse" and entry.get("capabilities", {}).get("disease") is not False:
            raise ValueError("release metadata species/capability mismatch")
    http = {}
    http_downloads = {}
    for species in ("human", "mouse"):
        status, collection, _ = _request(opener, base, f"/v1/releases/{release_id}/{species}/targets?limit=1")
        if status != 200:
            raise ValueError("HTTP species check failed")
        http[species] = _envelope(collection, release_id, species)
        status, downloads, _ = _request(opener, base, f"/v1/releases/{release_id}/{species}/downloads")
        if status != 200:
            raise ValueError("HTTP download discovery failed")
        http_downloads[species] = _envelope(downloads, release_id, species)
    request_id = 0
    protocol = "2026-07-28"

    def rpc(method: str, params: dict) -> dict:
        nonlocal request_id
        request_id += 1
        meta = {"io.modelcontextprotocol/protocolVersion": protocol, "io.modelcontextprotocol/clientCapabilities": {}, "io.modelcontextprotocol/clientInfo": {"name": "openantigens-python-smoke", "version": "1"}}
        status, value, headers = _request(opener, base, "/mcp", payload={"jsonrpc": "2.0", "id": request_id, "method": method, "params": {**params, "_meta": meta}}, protocol=protocol)
        if status != 200 or not isinstance(value, dict) or value.get("jsonrpc") != "2.0" or value.get("id") != request_id or "error" in value or not isinstance(value.get("result"), dict):
            raise ValueError("MCP wire check failed")
        if headers.get("Mcp-Session-Id"):
            raise ValueError("MCP endpoint unexpectedly requires a session")
        if headers.get("Cache-Control") != "no-store":
            raise ValueError("MCP responses require no-store caching")
        return value["result"]

    discovered = rpc("server/discover", {})
    if "tools" not in discovered.get("capabilities", {}):
        raise ValueError("MCP discovery did not advertise tools")
    listed = rpc("tools/list", {})
    inventory = listed.get("tools")
    if not isinstance(inventory, list) or not all(isinstance(tool, dict) for tool in inventory) or {tool.get("name") for tool in inventory} != TOOLS:
        raise ValueError("MCP tool inventory mismatch")

    def tool(name: str, arguments: dict) -> dict:
        result = rpc("tools/call", {"name": name, "arguments": arguments})
        if result.get("isError"):
            raise ValueError("MCP tool returned an error")
        return _envelope(result.get("structuredContent"), arguments.get("release_id", release_id), arguments.get("species"))

    selected = tool("get_release", {"release_id": release_id})
    _discovery(selected, release_id, expect_maintenance)
    if selected.get("data") != metadata.get("data") or selected.get("releases") != metadata.get("releases"):
        raise ValueError("HTTP and MCP release metadata differ")
    if expect_maintenance:
        current_result = rpc("tools/call", {"name": "get_release", "arguments": {}})
        content = current_result.get("content", [])
        error = json.loads(content[0].get("text", "null")) if isinstance(content, list) and content and isinstance(content[0], dict) else None
        if not current_result.get("isError") or not isinstance(error, dict) or error.get("code") != "publication_in_progress":
            raise ValueError("MCP current discovery is not in maintenance")
    elif tool("get_release", {}).get("data") != metadata.get("data"):
        raise ValueError("MCP current discovery differs from the expected release")
    for species in ("human", "mouse"):
        selected = tool("search_targets", {"release_id": release_id, "species": species, "limit": 1})
        if selected.get("data") != http[species].get("data") or selected.get("total") != http[species].get("total"):
            raise ValueError("HTTP and MCP species records differ")
        downloads = tool("get_downloads", {"release_id": release_id, "species": species})
        if downloads.get("data") != http_downloads[species].get("data"):
            raise ValueError("HTTP and MCP download metadata differ")
    extra_http_reads = 0
    for known_id in known_ids[1:]:
        status, pinned, headers = _request(opener, base, f"/v1/releases/{known_id}")
        extra_http_reads += 1
        if status != 200 or headers.get("Cache-Control") != "no-store":
            raise ValueError("caller-known retained release is unavailable")
        _envelope(pinned, known_id)
        _discovery(pinned, release_id, expect_maintenance)
        selected = tool("get_release", {"release_id": known_id})
        _discovery(selected, release_id, expect_maintenance)
        if selected.get("data") != pinned.get("data"):
            raise ValueError("retained release metadata differs across interfaces")
        for species in ("human", "mouse"):
            status, collection, _ = _request(opener, base, f"/v1/releases/{known_id}/{species}/targets?limit=1")
            extra_http_reads += 1
            if status != 200:
                raise ValueError("caller-known retained species is unavailable")
            _envelope(collection, known_id, species)
            selected = tool("search_targets", {"release_id": known_id, "species": species, "limit": 1})
            if selected.get("data") != collection.get("data") or selected.get("total") != collection.get("total"):
                raise ValueError("retained species records differ across interfaces")
    return {"base_url": base, "release_id": release_id, "retained_release_ids": ids, "retained_release_ids_source": "caller" if expect_maintenance else "service_discovery", "species": ["human", "mouse"], "http": "passed", "mcp_wire": "passed", "publication_state": "maintenance" if expect_maintenance else "active", "requests": 8 + extra_http_reads + request_id}


def activate_local_discovery(public_site_dir: Path, checked: dict) -> None:
    """Update both local metadata files after a successful, explicit smoke run."""
    base = _base_url(checked.get("base_url", ""))
    if urlsplit(base).scheme != "https" or checked.get("http") != "passed" or checked.get("mcp_wire") != "passed" or checked.get("species") != ["human", "mouse"] or not SHA.fullmatch(str(checked.get("release_id", ""))):
        raise ValueError("local discovery activation requires passed HTTPS checks for both species")
    if public_site_dir.is_symlink() or not public_site_dir.is_dir():
        raise ValueError("local static root must be a real directory")
    root = public_site_dir.resolve()
    originals = {}
    prepared = {}
    for name in ("portal_metadata.json", "mouse/portal_metadata.json"):
        path = root / name
        if path.is_symlink() or path.parent.is_symlink() or not path.is_file() or not path.resolve().is_relative_to(root):
            raise ValueError("both local portal metadata files must be regular files")
        original = path.read_bytes()
        metadata = json.loads(original)
        api = metadata.get("public_api") if isinstance(metadata, dict) else None
        if not isinstance(api, dict) or api.get("schema_version") != 1 or api.get("release_id") != checked["release_id"] or api.get("state") not in {"packaged", "active"}:
            raise ValueError("local portal API metadata release/state mismatch")
        originals[path] = original
        api.update(state="active", service_url=base, mcp_url=base + "/mcp")
        prepared[path] = (json.dumps(metadata, ensure_ascii=False, indent=2, allow_nan=False) + "\n").encode("utf-8")

    def replace(path: Path, value: bytes) -> None:
        with tempfile.NamedTemporaryFile(dir=path.parent, prefix=".api-discovery-", delete=False) as handle:
            temporary = Path(handle.name)
            handle.write(value)
        try:
            os.chmod(temporary, path.stat().st_mode & 0o777)
            os.replace(temporary, path)
        finally:
            temporary.unlink(missing_ok=True)

    updated = []
    try:
        for path, value in prepared.items():
            if path.read_bytes() != originals[path]:
                raise ValueError("local metadata changed during discovery activation")
            replace(path, value)
            updated.append(path)
    except BaseException:
        for path in updated:
            replace(path, originals[path])
        raise


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-url", required=True)
    parser.add_argument("--release-id", required=True)
    parser.add_argument("--retained-release-id", action="append", default=None, help="Expected previous retained ID; supply every retained ID to enforce the set.")
    parser.add_argument("--expect-maintenance", action="store_true")
    parser.add_argument("--activate-discovery-local-site", type=Path, help="After all checks pass, activate API links in both local portal metadata files; publishes nothing.")
    args = parser.parse_args(argv)
    try:
        result = check_api_service(args.base_url, args.release_id, expect_maintenance=args.expect_maintenance, retained_release_ids=args.retained_release_id)
        if args.activate_discovery_local_site is not None:
            activate_local_discovery(args.activate_discovery_local_site, result)
            result["local_discovery"] = "activated_for_both_species"
    except (ValueError, OSError, json.JSONDecodeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(result, sort_keys=True, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
