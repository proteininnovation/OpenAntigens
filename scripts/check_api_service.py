#!/usr/bin/env python3
"""Bounded HTTP and independent MCP wire smoke checks for an explicit release."""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import tempfile
import time
from pathlib import Path
from urllib.error import HTTPError
from urllib.parse import urlencode, urlsplit
from urllib.request import HTTPRedirectHandler, Request, build_opener

SHA = re.compile(r"[0-9a-f]{64}\Z")
TOOLS = {"get_release", "search_targets", "get_target", "list_constructs", "get_sequence"}
REFERENCE_ROLES = {
    "schema.json": ("application/schema+json", None),
    "citation.json": ("application/json", "citation"),
    "openantigens.bib": ("application/x-bibtex", None),
    "openantigens.ris": ("application/x-research-info-systems", None),
}
MAX_RESPONSE_BYTES = 1024 * 1024


def _json_loads(content):
    def unique(pairs):
        value = {}
        for key, item in pairs:
            if key in value:
                raise ValueError("JSON object keys must be unique")
            value[key] = item
        return value
    return json.loads(content, object_pairs_hook=unique)


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


def _wait(seconds: float, budget: dict) -> None:
    if seconds < 0 or budget["wait_seconds"] + seconds > 30:
        raise ValueError("smoke rate-limit wait budget exceeds thirty seconds")
    time.sleep(seconds)
    budget["wait_seconds"] += seconds


def _request(opener, base: str, path: str, *, budget: dict, payload: dict | None = None, protocol: str | None = None) -> tuple[int, object, object]:
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
    for attempt in range(3):
        if budget["requests"] >= 20:
            raise ValueError("smoke request budget exceeds twenty reads")
        scheduled = max(0, (budget["requests"] + 1 - 8) * 3)
        delay = max(0, scheduled - (time.monotonic() - budget["started_at"]))
        if delay:
            _wait(delay, budget)
        budget["requests"] += 1
        try:
            response = opener.open(request, timeout=5)
        except HTTPError as exc:
            response = exc
        with response:
            data = response.read(MAX_RESPONSE_BYTES + 1)
            if budget["requests"] == 1:
                budget["started_at"] = time.monotonic()
            if len(data) > MAX_RESPONSE_BYTES:
                raise ValueError("smoke response exceeds the byte limit")
            if response.status != 429:
                if response.headers.get_content_type() != "application/json":
                    raise ValueError("smoke response must be JSON")
                return response.status, _json_loads(data), response.headers
            retry = response.headers.get("Retry-After", "")
        if attempt == 2:
            raise ValueError("smoke rate limit persisted after three attempts")
        if not re.fullmatch(r"[1-9][0-9]*", retry):
            raise ValueError("rate-limited smoke response needs numeric Retry-After seconds")
        _wait(int(retry), budget)


def _provenance(value: object, release_id: str, manifest: object) -> None:
    if not isinstance(value, dict) or "input_sha256" in value or not isinstance(manifest, dict):
        raise ValueError("response provenance must use a compact immutable-manifest reference")
    if not isinstance(value.get("software_commit"), str) or not value["software_commit"] or not isinstance(value.get("snapshot_id"), str) or not value["snapshot_id"] or not isinstance(value.get("source_versions"), dict) or not isinstance(value.get("source_dates"), dict) or not isinstance(value.get("notes"), list):
        raise ValueError("response provenance identity is incomplete")
    hashes = value.get("input_hashes")
    if not isinstance(hashes, dict) or type(hashes.get("count")) is not int or hashes["count"] < 0 or hashes.get("json_pointer") != "/provenance/input_sha256" or hashes.get("manifest_sha256") != manifest.get("sha256") or not SHA.fullmatch(str(manifest.get("sha256", ""))):
        raise ValueError("response input hashes reference is invalid")
    url = urlsplit(str(hashes.get("manifest_url", "")))
    if hashes.get("manifest_url") != manifest.get("url") or url.scheme != "https" or not url.hostname or url.username or url.password or url.query or url.fragment or url.path != f"/data/releases/{release_id}/manifest.json":
        raise ValueError("response input hashes must reference their immutable release manifest")


def _envelope(value: object, release_id: str, species: str | None = None) -> dict:
    if not isinstance(value, dict) or value.get("release_id") != release_id or value.get("schema_version") != 1:
        raise ValueError("response release/schema mismatch")
    _provenance(value.get("provenance"), release_id, value.get("manifest"))
    if species is None and isinstance(value.get("data"), dict) and "species" in value["data"]:
        _provenance(value["data"].get("provenance"), release_id, value.get("manifest"))
    if species is not None:
        taxon = 9606 if species == "human" else 10090
        if value.get("species") != species or value.get("taxon_id") != taxon:
            raise ValueError("response species/taxon mismatch")
        if not isinstance(value.get("data"), list) or not 0 < len(value["data"]) <= 1 or type(value.get("total")) is not int or value["total"] < 1:
            raise ValueError("the explicit species smoke query must return a known target")
        for row in value["data"]:
            if not isinstance(row, dict) or row.get("species") != species or row.get("taxon_id") != taxon or not {"record_id", "target_id", "protein_name", "construct_count", "section_availability"}.issubset(row):
                raise ValueError("collection contains a record from another species")
            if not SHA.fullmatch(str(row["record_id"])) or row["target_id"] is not None and (not isinstance(row["target_id"], str) or not row["target_id"].startswith(f"{taxon}:")):
                raise ValueError("collection record identity does not match its species")
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
    if not ids or any(not isinstance(item, str) or not SHA.fullmatch(item) for item in ids) or len(ids) > 2 or len(set(ids)) != len(ids) or ids[0] != expected_current:
        raise ValueError("retained release discovery is invalid")
    return ids


def _references(metadata: dict, release_id: str) -> None:
    species_map = metadata.get("species")
    if not isinstance(species_map, dict) or set(species_map) != {"human", "mouse"}:
        raise ValueError("release metadata must describe both species")
    for species, taxon in (("human", 9606), ("mouse", 10090)):
        entry = species_map[species]
        capabilities = entry.get("capabilities") if isinstance(entry, dict) else None
        if not isinstance(entry, dict) or entry.get("taxon_id") != taxon or not isinstance(capabilities, dict) or species == "mouse" and capabilities.get("disease") is not False or "downloads" in entry:
            raise ValueError("release metadata species/capability/reference mismatch")
        references = entry.get("references")
        if not isinstance(references, list) or len(references) != len(REFERENCE_ROLES) or not all(isinstance(item, dict) and isinstance(item.get("name"), str) for item in references) or {item.get("name") for item in references} != set(REFERENCE_ROLES):
            raise ValueError("interactive release references must contain only schemas and citations")
        for item in references:
            media, schema = REFERENCE_ROLES[item["name"]]
            url = urlsplit(str(item.get("url", "")))
            if item.get("media_type") != media or item.get("schema") != schema or item.get("encoding") is not None or type(item.get("bytes")) is not int or item["bytes"] < 0 or not SHA.fullmatch(str(item.get("sha256", ""))) or item.get("path") != f"releases/{release_id}/{species}/downloads/{item['name']}" or url.scheme != "https" or not url.hostname or url.username or url.password or url.query or url.fragment or url.path != f"/data/releases/{release_id}/{species}/{item['name']}":
                raise ValueError("interactive reference metadata has an unapproved role or immutable URL")


def check_api_service(base_url: str, release_id: str, *, human_query: str, mouse_query: str, expect_maintenance: bool = False, retained_release_ids: list[str] | None = None) -> dict:
    """Bounded interactive checks: twenty attempts and thirty seconds of waiting."""
    base = _base_url(base_url)
    if not SHA.fullmatch(release_id) or any(not SHA.fullmatch(value) for value in retained_release_ids or []):
        raise ValueError("expected release IDs must be lowercase SHA-256 values")
    known_ids = list(dict.fromkeys([release_id, *(retained_release_ids or [])]))
    if len(known_ids) > 2:
        raise ValueError("smoke checks accept current plus at most one caller-known retained release")
    if not isinstance(human_query, str) or not isinstance(mouse_query, str):
        raise ValueError("both species smoke queries must be explicit strings")
    queries = {"human": human_query.strip(), "mouse": mouse_query.strip()}
    if any(not 2 <= len(query) <= 200 or "?" in query or "*" in query for query in queries.values()):
        raise ValueError("each species smoke query must contain two to two hundred literal characters")
    opener = build_opener(NoRedirect)
    budget = {"started_at": time.monotonic(), "requests": 0, "wait_seconds": 0.0}
    for path, state in (("/healthz", "ok"), ("/readyz", "ready")):
        status, value, headers = _request(opener, base, path, budget=budget)
        if status != 200 or not isinstance(value, dict) or value.get("status") != state or headers.get("Cache-Control") != "no-store":
            raise ValueError("service health/readiness check failed")
    status, current, headers = _request(opener, base, "/v1/releases/current", budget=budget)
    if expect_maintenance:
        if status != 503 or not isinstance(current, dict) or current.get("error", {}).get("code") != "publication_in_progress" or headers.get("Cache-Control") != "no-store" or not headers.get("Retry-After"):
            raise ValueError("current discovery is not in the expected maintenance state")
    elif status != 200 or headers.get("Cache-Control") != "no-store":
        raise ValueError("current discovery check failed")
    else:
        _envelope(current, release_id)
    status, metadata, headers = _request(opener, base, f"/v1/releases/{release_id}", budget=budget)
    if status != 200 or headers.get("Cache-Control") != "no-store":
        raise ValueError("literal release metadata check failed")
    _envelope(metadata, release_id)
    discovered_ids = _discovery(metadata, release_id, expect_maintenance)
    ids = known_ids if expect_maintenance else discovered_ids
    if not expect_maintenance and retained_release_ids is not None and set(ids) != set(known_ids):
        raise ValueError("retained release IDs do not match the expected release set")
    released = metadata.get("data")
    if not isinstance(released, dict):
        raise ValueError("release metadata is invalid")
    _references(released, release_id)
    http = {}
    for species in ("human", "mouse"):
        query = urlencode({"limit": 1, "q": queries[species]})
        status, collection, _ = _request(opener, base, f"/v1/releases/{release_id}/{species}/targets?{query}", budget=budget)
        if status != 200:
            raise ValueError("HTTP species check failed")
        http[species] = _envelope(collection, release_id, species)
    request_id = 0
    protocol = "2026-07-28"

    def rpc(method: str, params: dict) -> dict:
        nonlocal request_id
        request_id += 1
        meta = {"io.modelcontextprotocol/protocolVersion": protocol, "io.modelcontextprotocol/clientCapabilities": {}, "io.modelcontextprotocol/clientInfo": {"name": "openantigens-python-smoke", "version": "1"}}
        status, value, headers = _request(opener, base, "/mcp", budget=budget, payload={"jsonrpc": "2.0", "id": request_id, "method": method, "params": {**params, "_meta": meta}}, protocol=protocol)
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
    if not isinstance(inventory, list) or len(inventory) != len(TOOLS) or not all(isinstance(tool, dict) and isinstance(tool.get("name"), str) for tool in inventory) or {tool.get("name") for tool in inventory} != TOOLS:
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
        error = _json_loads(content[0].get("text", "null")) if isinstance(content, list) and content and isinstance(content[0], dict) else None
        if not current_result.get("isError") or not isinstance(error, dict) or error.get("code") != "publication_in_progress":
            raise ValueError("MCP current discovery is not in maintenance")
    elif tool("get_release", {}).get("data") != metadata.get("data"):
        raise ValueError("MCP current discovery differs from the expected release")
    for species in ("human", "mouse"):
        selected = tool("search_targets", {"release_id": release_id, "species": species, "q": queries[species], "limit": 1})
        if selected.get("data") != http[species].get("data") or selected.get("total") != http[species].get("total") or selected.get("refine_query") != http[species].get("refine_query"):
            raise ValueError("HTTP and MCP species records differ")
    for known_id in known_ids[1:]:
        status, pinned, headers = _request(opener, base, f"/v1/releases/{known_id}", budget=budget)
        if status != 200 or headers.get("Cache-Control") != "no-store":
            raise ValueError("caller-known retained release is unavailable")
        _envelope(pinned, known_id)
        _discovery(pinned, release_id, expect_maintenance)
        _references(pinned["data"], known_id)
        selected = tool("get_release", {"release_id": known_id})
        _discovery(selected, release_id, expect_maintenance)
        if selected.get("data") != pinned.get("data"):
            raise ValueError("retained release metadata differs across interfaces")
        for species in ("human", "mouse"):
            query = urlencode({"limit": 1, "q": queries[species]})
            status, collection, _ = _request(opener, base, f"/v1/releases/{known_id}/{species}/targets?{query}", budget=budget)
            if status != 200:
                raise ValueError("caller-known retained species is unavailable")
            _envelope(collection, known_id, species)
            selected = tool("search_targets", {"release_id": known_id, "species": species, "q": queries[species], "limit": 1})
            if selected.get("data") != collection.get("data") or selected.get("total") != collection.get("total") or selected.get("refine_query") != collection.get("refine_query"):
                raise ValueError("retained species records differ across interfaces")
    return {"base_url": base, "release_id": release_id, "retained_release_ids": ids, "retained_release_ids_source": "caller" if expect_maintenance else "service_discovery", "species": ["human", "mouse"], "http": "passed", "mcp_wire": "passed", "publication_state": "maintenance" if expect_maintenance else "active", "requests": budget["requests"], "wait_seconds": round(budget["wait_seconds"], 3)}


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
        metadata = _json_loads(original)
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
    parser.add_argument("--human-query", required=True, help="Explicit literal query for a known human target; at least two characters.")
    parser.add_argument("--mouse-query", required=True, help="Explicit literal query for a known mouse target; at least two characters.")
    parser.add_argument("--retained-release-id", action="append", default=None, help="Expected previous retained ID; supply every retained ID to enforce the set.")
    parser.add_argument("--expect-maintenance", action="store_true")
    parser.add_argument("--activate-discovery-local-site", type=Path, help="After all checks pass, activate API links in both local portal metadata files; publishes nothing.")
    args = parser.parse_args(argv)
    try:
        result = check_api_service(args.base_url, args.release_id, human_query=args.human_query, mouse_query=args.mouse_query, expect_maintenance=args.expect_maintenance, retained_release_ids=args.retained_release_id)
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
