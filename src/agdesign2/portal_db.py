from __future__ import annotations

import hashlib
import json
import os
import shutil
import sqlite3
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .portal import (
    _apply_open_targets_to_entry,
    _copy_brand_assets,
    _copy_portal_structure,
    _copy_report_assets,
    _externalize_inline_script,
    _externalize_plain_scripts,
    _ensure_vendor_assets,
    _load_open_targets_index,
    _portal_css,
    _render_index_data_js,
    _resolve_existing_path,
    _safe_filename,
    _split_report_viewer_runtime,
    _structure_source_label,
    _write_text_atomic,
    _write_portal_htaccess,
    portal_build_metadata,
    portal_download_manifest,
    portal_index_json,
    portal_index_tsv,
    render_builder_page,
    render_calculator_page,
    render_constructs_page,
    render_detail_page,
    render_downloads_page,
    render_help_page,
    render_index_page,
    render_methods_page,
    render_privacy_page,
    render_terms_page,
)
from .dynamic_portal import _build_dynamic_index_entry


SCHEMA_VERSION = 1


class PortalDatabaseError(RuntimeError):
    """Raised when a portal SQLite snapshot cannot be built or read safely."""


def build_portal_database(
    summary_path: str | Path,
    db_path: str | Path,
    *,
    assets_dir: str | Path | None = None,
    bundle_vendor_assets: bool = False,
) -> Path:
    summary_file = Path(summary_path).resolve()
    if not summary_file.exists():
        raise FileNotFoundError(f"Summary file not found: {summary_file}")
    database_path = Path(db_path).resolve()
    asset_root = Path(assets_dir).resolve() if assets_dir is not None else database_path.parent / "assets"
    batch_dir = summary_file.parent
    results = _load_summary_rows(summary_file)

    database_path.parent.mkdir(parents=True, exist_ok=True)
    asset_root.mkdir(parents=True, exist_ok=True)
    structures_dir = asset_root / "structures"
    report_assets_dir = asset_root / "report_assets"
    structures_dir.mkdir(parents=True, exist_ok=True)
    report_assets_dir.mkdir(parents=True, exist_ok=True)
    _copy_brand_assets(asset_root)
    _write_text_atomic(asset_root / "portal.css", _portal_css())
    if bundle_vendor_assets:
        _ensure_vendor_assets(portal_dir=asset_root, batch_dir=batch_dir)

    disease_index = _load_open_targets_index(batch_dir)
    entries: list[dict[str, Any]] = []
    reports: list[dict[str, Any]] = []
    seen_entries: set[str] = set()
    seen_details: set[str] = set()
    assets: dict[str, dict[str, Any]] = {}

    for row_number, item in enumerate(results, start=1):
        report_json_path = _resolve_existing_path(item.get("json_report"), batch_dir)
        if report_json_path is None or not report_json_path.exists():
            raise PortalDatabaseError(f"Row {row_number} is missing report JSON: {item.get('json_report')!r}")
        report_data = _load_report(report_json_path)
        _validate_report(report_data, report_json_path)
        entry = _build_dynamic_index_entry(item, report_data, report_json_path)
        _apply_open_targets_to_entry(entry, disease_index)

        entry_name = str(entry.get("entry_name") or "").strip()
        detail_page = str(entry.get("detail_page") or "").strip()
        if not entry_name:
            raise PortalDatabaseError(f"Row {row_number} produced an empty entry_name")
        if not detail_page:
            raise PortalDatabaseError(f"Row {row_number} for {entry_name} produced an empty detail_page")
        entry_key = entry_name.upper()
        detail_key = detail_page.lower()
        if entry_key in seen_entries:
            raise PortalDatabaseError(f"Duplicate entry_name in portal database: {entry_name}")
        if detail_key in seen_details:
            raise PortalDatabaseError(f"Duplicate detail_page in portal database: {detail_page}")
        seen_entries.add(entry_key)
        seen_details.add(detail_key)

        structure_path = _copy_portal_structure(report_data=report_data, structures_dir=structures_dir, batch_dir=batch_dir)
        if structure_path:
            entry["portal_structure_path"] = _asset_reference(Path(structure_path), asset_root)
            entry["has_alphafold_structure"] = True
            entry["has_any_structure"] = True
            entry["structure_source"] = _structure_source_label(Path(structure_path))
            _record_asset(assets, Path(structure_path), asset_root, "structure")
        copied_assets = _copy_report_assets(
            report_data=report_data,
            report_json_path=report_json_path,
            report_assets_dir=report_assets_dir,
            batch_dir=batch_dir,
        )
        if copied_assets is not None:
            for file_path in copied_assets.rglob("*"):
                if file_path.is_file():
                    _record_asset(assets, file_path, asset_root, "report_asset")

        entry["json_report"] = ""
        entry["markdown_report"] = ""
        entries.append(entry)
        reports.append({"entry_name": entry_name, "detail_page": detail_page, "payload": report_data})

    for path in asset_root.iterdir():
        if path.is_file():
            _record_asset(assets, path, asset_root, "portal_asset")
    disease_payloads = _load_disease_payloads(batch_dir)

    with sqlite3.connect(database_path) as conn:
        _initialize_schema(conn)
        conn.execute("DELETE FROM metadata")
        conn.execute("DELETE FROM entries")
        conn.execute("DELETE FROM reports")
        conn.execute("DELETE FROM assets")
        conn.execute("DELETE FROM derived_files")
        _insert_metadata(conn, database_path, asset_root, entries)
        for entry in entries:
            conn.execute(
                """
                INSERT INTO entries(entry_name, detail_page, gene_symbol, status, entry_json)
                VALUES (?, ?, ?, ?, ?)
                """,
                (
                    str(entry["entry_name"]),
                    str(entry["detail_page"]),
                    str(entry.get("gene_symbol") or ""),
                    str(entry.get("status") or ""),
                    json.dumps(entry, sort_keys=True),
                ),
            )
        for report in reports:
            conn.execute(
                """
                INSERT INTO reports(entry_name, detail_page, report_json)
                VALUES (?, ?, ?)
                """,
                (
                    report["entry_name"],
                    report["detail_page"],
                    json.dumps(report["payload"], sort_keys=True),
                ),
            )
        for asset in sorted(assets.values(), key=lambda item: item["path"]):
            conn.execute(
                """
                INSERT INTO assets(path, kind, size_bytes, sha256)
                VALUES (?, ?, ?, ?)
                """,
                (asset["path"], asset["kind"], asset["size_bytes"], asset["sha256"]),
            )
        for path, payload in disease_payloads.items():
            conn.execute(
                """
                INSERT INTO derived_files(path, media_type, content)
                VALUES (?, ?, ?)
                """,
                (path, _media_type_for_path(path), payload),
            )
    return database_path


def load_index_entries_from_db(db_path: str | Path) -> list[dict[str, Any]]:
    with _connect_readonly(db_path) as conn:
        rows = conn.execute("SELECT entry_json FROM entries ORDER BY rowid").fetchall()
    return [json.loads(row[0]) for row in rows]


def load_report_from_db(db_path: str | Path, entry_or_detail: str) -> dict[str, Any]:
    key = str(entry_or_detail or "").strip()
    if not key:
        raise PortalDatabaseError("Missing entry name or detail page")
    with _connect_readonly(db_path) as conn:
        row = conn.execute(
            """
            SELECT report_json FROM reports
            WHERE UPPER(entry_name) = UPPER(?) OR LOWER(detail_page) = LOWER(?)
            """,
            (key, key),
        ).fetchone()
    if row is None:
        raise KeyError(key)
    return json.loads(row[0])


def load_entry_from_db(db_path: str | Path, entry_or_detail: str) -> dict[str, Any]:
    key = str(entry_or_detail or "").strip()
    if not key:
        raise PortalDatabaseError("Missing entry name or detail page")
    with _connect_readonly(db_path) as conn:
        row = conn.execute(
            """
            SELECT entry_json FROM entries
            WHERE UPPER(entry_name) = UPPER(?) OR LOWER(detail_page) = LOWER(?)
            """,
            (key, key),
        ).fetchone()
    if row is None:
        raise KeyError(key)
    return json.loads(row[0])


def load_derived_file_from_db(db_path: str | Path, path: str) -> tuple[str, str] | None:
    normalized = str(path or "").strip().lstrip("/")
    with _connect_readonly(db_path) as conn:
        row = conn.execute("SELECT media_type, content FROM derived_files WHERE path = ?", (normalized,)).fetchone()
    if row is None:
        return None
    return str(row[0]), str(row[1])


def build_portal_from_database(
    db_path: str | Path,
    output_dir: str | Path,
    *,
    assets_dir: str | Path | None = None,
    bundle_vendor_assets: bool = False,
) -> Path:
    database_path = Path(db_path).resolve()
    asset_root = _resolve_assets_dir(database_path, assets_dir)
    portal_dir = Path(output_dir).resolve()
    reports_dir = portal_dir / "reports"
    structures_dir = portal_dir / "structures"
    report_assets_dir = portal_dir / "report_assets"
    report_scripts_dir = portal_dir / "report_scripts"
    for directory in (portal_dir, reports_dir, structures_dir, report_assets_dir, report_scripts_dir):
        directory.mkdir(parents=True, exist_ok=True)
    _copy_asset_tree(asset_root, portal_dir)
    _write_portal_htaccess(portal_dir)
    if bundle_vendor_assets:
        from .portal import _ensure_vendor_assets

        _ensure_vendor_assets(portal_dir=portal_dir, batch_dir=asset_root)

    entries = load_index_entries_from_db(database_path)
    _write_portal_download_files_from_db(portal_dir, database_path, entries)
    index_html, index_js = _externalize_inline_script(render_index_page(entries), script_src="portal-index.js")
    _write_text_atomic(portal_dir / "index.html", index_html)
    _write_text_atomic(portal_dir / "portal-index-data.js", _render_index_data_js(entries))
    _write_text_atomic(portal_dir / "portal-index.js", index_js)
    shared_report_viewer_js = ""
    for entry in entries:
        detail_page = str(entry.get("detail_page") or "")
        if not detail_page:
            continue
        report_data = load_report_from_db(database_path, detail_page)
        detail_html = render_detail_page(entry, report_data, batch_dir=asset_root, portal_dir=portal_dir)
        detail_script_name = _safe_filename(str(Path(detail_page).with_suffix("").name)) + ".js"
        detail_html, detail_js = _externalize_plain_scripts(detail_html, script_src=f"../report_scripts/{detail_script_name}")
        detail_js, shared_js = _split_report_viewer_runtime(detail_js)
        if shared_js:
            shared_report_viewer_js = shared_js
            detail_html = detail_html.replace(
                f'<script src="../report_scripts/{detail_script_name}"></script>',
                f'<script src="../report-viewer.js"></script>\n    <script src="../report_scripts/{detail_script_name}"></script>',
                1,
            )
        if detail_js:
            _write_text_atomic(report_scripts_dir / detail_script_name, detail_js)
        _write_text_atomic(reports_dir / detail_page, detail_html)
    if shared_report_viewer_js:
        _write_text_atomic(portal_dir / "report-viewer.js", shared_report_viewer_js)
    _write_text_atomic(portal_dir / "help.html", render_help_page())
    _write_text_atomic(portal_dir / "builder.html", render_builder_page())
    _write_text_atomic(portal_dir / "constructs.html", render_constructs_page())
    _write_text_atomic(portal_dir / "methods.html", render_methods_page())
    _write_text_atomic(portal_dir / "downloads.html", render_downloads_page(entries))
    _write_text_atomic(portal_dir / "calculator.html", render_calculator_page())
    _write_text_atomic(portal_dir / "terms.html", render_terms_page())
    _write_text_atomic(portal_dir / "privacy.html", render_privacy_page())
    _write_text_atomic(portal_dir / "portal_metadata.json", portal_build_metadata(entries))
    return portal_dir / "index.html"


def create_db_portal_app(db_path: str | Path, *, assets_dir: str | Path | None = None):
    try:
        from fastapi import FastAPI, HTTPException
        from fastapi.middleware.cors import CORSMiddleware
        from fastapi.responses import FileResponse, HTMLResponse, JSONResponse, PlainTextResponse, RedirectResponse
    except Exception as exc:  # pragma: no cover - dependency/runtime guard
        raise RuntimeError("DB-backed portal requires FastAPI and uvicorn. Reinstall the project dependencies first.") from exc

    database_path = Path(db_path).resolve()
    if not database_path.exists():
        raise FileNotFoundError(f"Portal database not found: {database_path}")
    asset_root = _resolve_assets_dir(database_path, assets_dir)
    app = FastAPI(title="OpenAntigens SQLite Portal")
    app.add_middleware(
        CORSMiddleware,
        allow_origins=["*"],
        allow_credentials=False,
        allow_methods=["*"],
        allow_headers=["*"],
    )

    def entries() -> list[dict[str, Any]]:
        return load_index_entries_from_db(database_path)

    def detail_payload(key: str) -> tuple[dict[str, Any], dict[str, Any]]:
        try:
            return load_entry_from_db(database_path, key), load_report_from_db(database_path, key)
        except KeyError:
            raise HTTPException(status_code=404, detail="Gene report not found")

    @app.get("/", include_in_schema=False)
    def root():
        return RedirectResponse(url="/index.html", status_code=307)

    @app.get("/index.html", response_class=HTMLResponse, include_in_schema=False)
    def index():
        return HTMLResponse(render_index_page(entries()))

    @app.get("/help.html", response_class=HTMLResponse, include_in_schema=False)
    def help_page():
        return HTMLResponse(render_help_page())

    @app.get("/builder.html", response_class=HTMLResponse, include_in_schema=False)
    def builder_page():
        return HTMLResponse(render_builder_page())

    @app.get("/constructs.html", response_class=HTMLResponse, include_in_schema=False)
    def constructs_page():
        return HTMLResponse(render_constructs_page())

    @app.get("/methods.html", response_class=HTMLResponse, include_in_schema=False)
    def methods_page():
        return HTMLResponse(render_methods_page())

    @app.get("/downloads.html", response_class=HTMLResponse, include_in_schema=False)
    def downloads_page():
        return HTMLResponse(render_downloads_page(entries()))

    @app.get("/calculator.html", response_class=HTMLResponse, include_in_schema=False)
    def calculator_page():
        return HTMLResponse(render_calculator_page())

    @app.get("/terms.html", response_class=HTMLResponse, include_in_schema=False)
    def terms_page():
        return HTMLResponse(render_terms_page())

    @app.get("/privacy.html", response_class=HTMLResponse, include_in_schema=False)
    def privacy_page():
        return HTMLResponse(render_privacy_page())

    @app.get("/reports/{detail_page}", response_class=HTMLResponse, include_in_schema=False)
    def report_page(detail_page: str):
        entry, report = detail_payload(detail_page)
        return HTMLResponse(render_detail_page(entry, report, batch_dir=asset_root, portal_dir=asset_root))

    @app.get("/portal.css", response_class=PlainTextResponse, include_in_schema=False)
    def portal_css():
        return PlainTextResponse(_portal_css(), media_type="text/css")

    @app.get("/portal-index-data.js", response_class=PlainTextResponse, include_in_schema=False)
    def portal_index_data():
        return PlainTextResponse(_render_index_data_js(entries()), media_type="application/javascript")

    @app.get("/api/genes")
    def api_genes():
        return JSONResponse(entries())

    @app.get("/api/gene/{entry_name}")
    def api_gene(entry_name: str):
        _, report = detail_payload(entry_name)
        return JSONResponse(report)

    @app.get("/downloads/agdesign2_portal_index.tsv", response_class=PlainTextResponse, include_in_schema=False)
    def download_index_tsv():
        return PlainTextResponse(portal_index_tsv(entries()), media_type="text/tab-separated-values")

    @app.get("/downloads/agdesign2_portal_index.json", response_class=PlainTextResponse, include_in_schema=False)
    def download_index_json():
        return PlainTextResponse(portal_index_json(entries()), media_type="application/json")

    @app.get("/downloads/download_manifest.json", response_class=PlainTextResponse, include_in_schema=False)
    def download_manifest():
        return PlainTextResponse(portal_download_manifest(entries()), media_type="application/json")

    @app.get("/portal_metadata.json", response_class=PlainTextResponse, include_in_schema=False)
    def portal_metadata():
        return PlainTextResponse(portal_build_metadata(entries()), media_type="application/json")

    @app.get("/{file_path:path}", include_in_schema=False)
    def asset_or_derived_file(file_path: str):
        normalized = str(file_path or "").strip().lstrip("/")
        derived = load_derived_file_from_db(database_path, normalized)
        if derived is not None:
            media_type, content = derived
            return PlainTextResponse(content, media_type=media_type)
        path = (asset_root / normalized).resolve()
        if not _is_relative_to(path, asset_root) or not path.exists() or not path.is_file():
            raise HTTPException(status_code=404, detail="File not found")
        return FileResponse(path)

    return app


def run_db_portal(
    db_path: str | Path,
    *,
    assets_dir: str | Path | None = None,
    host: str = "127.0.0.1",
    port: int = 8000,
) -> None:
    try:
        import uvicorn
    except Exception as exc:  # pragma: no cover - dependency/runtime guard
        raise RuntimeError("DB-backed portal requires uvicorn. Reinstall the project dependencies first.") from exc
    uvicorn.run(create_db_portal_app(db_path, assets_dir=assets_dir), host=host, port=port)


def _initialize_schema(conn: sqlite3.Connection) -> None:
    conn.executescript(
        """
        PRAGMA foreign_keys = ON;
        CREATE TABLE IF NOT EXISTS metadata (
            key TEXT PRIMARY KEY,
            value TEXT NOT NULL
        );
        CREATE TABLE IF NOT EXISTS entries (
            entry_name TEXT PRIMARY KEY,
            detail_page TEXT NOT NULL UNIQUE,
            gene_symbol TEXT NOT NULL,
            status TEXT NOT NULL,
            entry_json TEXT NOT NULL
        );
        CREATE TABLE IF NOT EXISTS reports (
            entry_name TEXT PRIMARY KEY REFERENCES entries(entry_name) ON DELETE CASCADE,
            detail_page TEXT NOT NULL UNIQUE,
            report_json TEXT NOT NULL
        );
        CREATE TABLE IF NOT EXISTS assets (
            path TEXT PRIMARY KEY,
            kind TEXT NOT NULL,
            size_bytes INTEGER NOT NULL,
            sha256 TEXT NOT NULL
        );
        CREATE TABLE IF NOT EXISTS derived_files (
            path TEXT PRIMARY KEY,
            media_type TEXT NOT NULL,
            content TEXT NOT NULL
        );
        """
    )


def _insert_metadata(conn: sqlite3.Connection, database_path: Path, asset_root: Path, entries: list[dict[str, Any]]) -> None:
    values = {
        "schema_version": str(SCHEMA_VERSION),
        "assets_dir": os.path.relpath(asset_root, database_path.parent),
        "generated_at_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "target_count": str(len(entries)),
    }
    conn.executemany("INSERT INTO metadata(key, value) VALUES (?, ?)", sorted(values.items()))


def _load_summary_rows(summary_file: Path) -> list[dict[str, Any]]:
    try:
        payload = json.loads(summary_file.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise PortalDatabaseError(f"Invalid batch summary JSON: {summary_file}") from exc
    if not isinstance(payload, list):
        raise PortalDatabaseError(f"Batch summary must be a list: {summary_file}")
    return payload


def _load_report(report_json_path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(report_json_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise PortalDatabaseError(f"Invalid report JSON: {report_json_path}") from exc
    if not isinstance(payload, dict):
        raise PortalDatabaseError(f"Report JSON must be an object: {report_json_path}")
    return payload


def _validate_report(report_data: dict[str, Any], report_json_path: Path) -> None:
    target = report_data.get("target")
    if not isinstance(target, dict):
        raise PortalDatabaseError(f"Report is missing target object: {report_json_path}")
    for field in ("entry_name", "protein_name", "accession"):
        if not str(target.get(field) or "").strip():
            raise PortalDatabaseError(f"Report {report_json_path} is missing target.{field}")


def _record_asset(assets: dict[str, dict[str, Any]], path: Path, asset_root: Path, kind: str) -> None:
    if not path.exists() or not path.is_file():
        raise PortalDatabaseError(f"Referenced asset is missing: {path}")
    relative = _asset_reference(path, asset_root)
    assets[relative] = {
        "path": relative,
        "kind": kind,
        "size_bytes": path.stat().st_size,
        "sha256": _sha256(path),
    }


def _asset_reference(path: Path, asset_root: Path) -> str:
    resolved = path.resolve()
    root = asset_root.resolve()
    if not _is_relative_to(resolved, root):
        raise PortalDatabaseError(f"Asset path is outside assets directory: {resolved}")
    return resolved.relative_to(root).as_posix()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_disease_payloads(batch_dir: Path) -> dict[str, str]:
    payloads: dict[str, str] = {}
    for suffix in ("tsv", "json"):
        path = batch_dir / f"open_targets_disease_associations.{suffix}"
        if path.exists():
            payloads[path.name] = path.read_text(encoding="utf-8")
    return payloads


def _media_type_for_path(path: str) -> str:
    if path.endswith(".json"):
        return "application/json"
    if path.endswith(".tsv"):
        return "text/tab-separated-values"
    if path.endswith(".css"):
        return "text/css"
    return "text/plain"


def _connect_readonly(db_path: str | Path) -> sqlite3.Connection:
    database_path = Path(db_path).resolve()
    if not database_path.exists():
        raise FileNotFoundError(f"Portal database not found: {database_path}")
    conn = sqlite3.connect(f"file:{database_path}?mode=ro", uri=True)
    conn.execute("PRAGMA foreign_keys = ON")
    return conn


def _resolve_assets_dir(db_path: Path, assets_dir: str | Path | None) -> Path:
    if assets_dir is not None:
        asset_root = Path(assets_dir).resolve()
    else:
        with _connect_readonly(db_path) as conn:
            row = conn.execute("SELECT value FROM metadata WHERE key = 'assets_dir'").fetchone()
        asset_root = (db_path.parent / str(row[0])).resolve() if row is not None else db_path.parent / "assets"
    if not asset_root.exists():
        raise FileNotFoundError(f"Portal assets directory not found: {asset_root}")
    return asset_root


def _copy_asset_tree(source_root: Path, destination_root: Path) -> None:
    for source in source_root.rglob("*"):
        relative = source.relative_to(source_root)
        destination = destination_root / relative
        if source.is_dir():
            destination.mkdir(parents=True, exist_ok=True)
        elif source.is_file():
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)


def _write_portal_download_files_from_db(portal_dir: Path, db_path: Path, entries: list[dict[str, Any]]) -> None:
    downloads_dir = portal_dir / "downloads"
    downloads_dir.mkdir(parents=True, exist_ok=True)
    _write_text_atomic(downloads_dir / "agdesign2_portal_index.tsv", portal_index_tsv(entries))
    _write_text_atomic(downloads_dir / "agdesign2_portal_index.json", portal_index_json(entries))
    _write_text_atomic(downloads_dir / "download_manifest.json", portal_download_manifest(entries))
    for name in ("open_targets_disease_associations.tsv", "open_targets_disease_associations.json"):
        derived = load_derived_file_from_db(db_path, name)
        if derived is not None:
            _write_text_atomic(portal_dir / name, derived[1])


def _is_relative_to(path: Path, parent: Path) -> bool:
    try:
        path.resolve().relative_to(parent.resolve())
    except ValueError:
        return False
    return True
