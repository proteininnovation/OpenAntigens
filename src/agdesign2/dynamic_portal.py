from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .site_docs import agent_guide_markdown, citation_downloads, site_document_files

from .portal import (
    _copy_brand_assets,
    _copy_portal_structure,
    _copy_report_assets,
    _format_region,
    _index_cross_reactivity_count,
    _index_family_name,
    _is_valid_threedmol_asset,
    _load_open_targets_index,
    _apply_open_targets_to_entry,
    _portal_css,
    _refresh_reports_if_needed,
    _render_disease_index_js,
    _render_index_data_js,
    _resolve_alphafold_artifact,
    _available_disease_downloads,
    _resolve_existing_path,
    _safe_filename,
    _structure_source_label,
    _target_aliases,
    portal_download_manifest,
    portal_build_metadata,
    portal_index_json,
    portal_index_tsv,
    render_builder_page,
    render_calculator_page,
    render_constructs_page,
    render_detail_page,
    render_downloads_page,
    render_help_page,
    render_agent_guide_page,
    render_index_page,
    render_methods_page,
    render_privacy_page,
    render_terms_page,
)

_PUBLIC_ROOT_ASSETS = {
    "apple-touch-icon.png",
    "favicon-32.png",
    "favicon.ico",
    "icon-192.png",
    "icon-512.png",
    "site.webmanifest",
}


def _is_within_directory(path: Path, directory: Path) -> bool:
    """True only if ``path`` resolves to a location inside ``directory``.

    Replaces a string ``startswith`` check, which let a sibling whose path is a
    string prefix of the served directory (e.g. ``/data/run`` vs
    ``/data/run-secrets``) slip through.
    """
    try:
        path.resolve().relative_to(directory.resolve())
        return True
    except ValueError:
        return False


def create_app(
    summary_path: str | Path,
    *,
    refresh_reports: bool = False,
    refresh_stale_reports: bool = False,
    verbose: bool = False,
    enable_complex_portal: bool = False,
):
    try:
        from fastapi import FastAPI, HTTPException
        from fastapi.responses import FileResponse, HTMLResponse, JSONResponse, PlainTextResponse, RedirectResponse
    except Exception as exc:  # pragma: no cover - dependency/runtime guard
        raise RuntimeError(
            "Dynamic portal requires FastAPI and uvicorn. Reinstall the project dependencies first."
        ) from exc

    summary_file = Path(summary_path).resolve()
    if not summary_file.exists():
        raise FileNotFoundError(f"Summary file not found: {summary_file}")
    batch_dir = summary_file.parent
    portal_dir = batch_dir / "portal"
    structures_dir = portal_dir / "structures"
    structures_dir.mkdir(parents=True, exist_ok=True)
    report_assets_dir = portal_dir / "report_assets"
    report_assets_dir.mkdir(parents=True, exist_ok=True)
    _copy_brand_assets(portal_dir)
    (portal_dir / "portal.css").write_text(_portal_css(), encoding="utf-8")

    app = FastAPI(title="OpenAntigens Dynamic Portal")
    cache: dict[str, Any] = {
        "summary_mtime_ns": None,
        "results": None,
        "entries": None,
        "index_by_detail": None,
        "disease_mtime_ns": None,
        "disease_index": None,
    }

    def maybe_refresh_results() -> list[dict[str, Any]]:
        needs_refresh = refresh_reports or refresh_stale_reports
        if not needs_refresh:
            return json.loads(summary_file.read_text(encoding="utf-8"))
        return _refresh_reports_if_needed(
            results=json.loads(summary_file.read_text(encoding="utf-8")),
            summary_file=summary_file,
            batch_dir=batch_dir,
            refresh_reports=refresh_reports,
            refresh_stale_reports=refresh_stale_reports,
            verbose=verbose,
            enable_complex_portal=enable_complex_portal,
        )

    def load_results() -> list[dict[str, Any]]:
        summary_mtime_ns = summary_file.stat().st_mtime_ns
        cached_mtime_ns = cache["summary_mtime_ns"]
        if cache["results"] is not None and cached_mtime_ns == summary_mtime_ns:
            return cache["results"]
        results = maybe_refresh_results()
        cache["summary_mtime_ns"] = summary_file.stat().st_mtime_ns
        cache["results"] = results
        cache["entries"] = None
        cache["index_by_detail"] = None
        return results

    def load_entries() -> list[dict[str, Any]]:
        results = load_results()
        disease_path = batch_dir / "open_targets_disease_associations.json"
        current_disease_mtime = disease_path.stat().st_mtime_ns if disease_path.exists() else None
        if cache["disease_mtime_ns"] != current_disease_mtime:
            cache["entries"] = None
            cache["index_by_detail"] = None
            cache["disease_index"] = None
        cached_entries = cache["entries"]
        if cached_entries is not None:
            return cached_entries
        entries: list[dict[str, Any]] = []
        index_by_detail: dict[str, dict[str, Any]] = {}
        for item in results:
            report_json_path = _resolve_existing_path(item.get("json_report"), batch_dir)
            report_data: dict[str, Any] | None = None
            if report_json_path is not None:
                try:
                    report_data = json.loads(report_json_path.read_text(encoding="utf-8"))
                except (OSError, json.JSONDecodeError) as exc:
                    raise ValueError(f"Could not read report JSON: {report_json_path}") from exc
            entry = _build_dynamic_index_entry(item, report_data, report_json_path)
            if report_data is not None:
                entry["portal_structure_path"] = _copy_portal_structure(
                    report_data=report_data,
                    structures_dir=structures_dir,
                    batch_dir=batch_dir,
                )
                _copy_report_assets(
                    report_data=report_data,
                    report_json_path=report_json_path,
                    report_assets_dir=report_assets_dir,
                    batch_dir=batch_dir,
                )
                if entry["portal_structure_path"]:
                    entry["has_alphafold_structure"] = True
                    entry["has_any_structure"] = True
            _apply_open_targets_to_entry(entry, load_disease_index())
            entries.append(entry)
            detail_page = entry.get("detail_page")
            if detail_page:
                index_by_detail[str(detail_page)] = entry
        cache["entries"] = entries
        cache["index_by_detail"] = index_by_detail
        return entries

    def load_disease_index() -> dict[str, dict[str, Any]]:
        disease_path = batch_dir / "open_targets_disease_associations.json"
        mtime = disease_path.stat().st_mtime_ns if disease_path.exists() else None
        if cache["disease_index"] is not None and cache["disease_mtime_ns"] == mtime:
            return cache["disease_index"]
        cache["disease_mtime_ns"] = mtime
        cache["disease_index"] = _load_open_targets_index(batch_dir)
        return cache["disease_index"]

    def load_detail_payload(detail_page: str) -> dict[str, Any] | None:
        load_entries()
        entry = (cache["index_by_detail"] or {}).get(detail_page)
        if entry is None:
            return None
        report_json_path = _resolve_existing_path(entry.get("json_report"), batch_dir)
        if report_json_path is None:
            return None
        try:
            report_data = json.loads(report_json_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise ValueError(f"Could not read report JSON: {report_json_path}") from exc
        entry = dict(entry)
        entry["portal_structure_path"] = _copy_portal_structure(
            report_data=report_data,
            structures_dir=structures_dir,
            batch_dir=batch_dir,
        )
        _copy_report_assets(
            report_data=report_data,
            report_json_path=report_json_path,
            report_assets_dir=report_assets_dir,
            batch_dir=batch_dir,
        )
        if entry["portal_structure_path"]:
            entry["has_alphafold_structure"] = True
            entry["has_any_structure"] = True
            entry["structure_source"] = _structure_source_label(Path(str(entry["portal_structure_path"])))
        return {"entry": entry, "report": report_data}

    @app.get("/", include_in_schema=False)
    def root():
        return RedirectResponse(url="/index.html", status_code=307)

    @app.get("/index.html", response_class=HTMLResponse, include_in_schema=False)
    def index():
        return HTMLResponse(render_index_page(load_entries()))

    @app.get("/agent-guide.html", response_class=HTMLResponse, include_in_schema=False)
    def agent_guide_page():
        return HTMLResponse(render_agent_guide_page())

    @app.get("/llms.txt", response_class=PlainTextResponse, include_in_schema=False)
    def agent_guide_text():
        return PlainTextResponse(agent_guide_markdown())

    @app.get("/agent-guide.js", response_class=PlainTextResponse, include_in_schema=False)
    def agent_guide_script():
        return PlainTextResponse(site_document_files()["agent-guide.js"], media_type="text/javascript")

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
        return HTMLResponse(render_downloads_page(load_entries(), disease_downloads=_available_disease_downloads(batch_dir)))

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
        payload = load_detail_payload(detail_page)
        if payload is None:
            raise HTTPException(status_code=404, detail="Report page not found")
        return HTMLResponse(render_detail_page(payload["entry"], payload["report"], batch_dir=batch_dir, portal_dir=portal_dir))

    @app.get("/portal.css", response_class=PlainTextResponse, include_in_schema=False)
    def portal_css():
        return PlainTextResponse(_portal_css(), media_type="text/css")

    @app.get("/portal-index-data.js", response_class=PlainTextResponse, include_in_schema=False)
    def portal_index_data():
        return PlainTextResponse(
            _render_index_data_js(load_entries()),
            media_type="application/javascript",
        )

    @app.get("/portal-disease-index.js", response_class=PlainTextResponse, include_in_schema=False)
    def portal_disease_index():
        return PlainTextResponse(
            _render_disease_index_js(load_entries()),
            media_type="application/javascript",
        )

    @app.get("/3Dmol-min.js", include_in_schema=False)
    def threedmol_script():
        path = portal_dir / "3Dmol-min.js"
        if not _is_valid_threedmol_asset(path):
            raise HTTPException(status_code=404, detail="Vendor asset not found")
        return FileResponse(path, media_type="application/javascript")

    @app.get("/report_assets/{asset_path:path}", include_in_schema=False)
    def report_asset(asset_path: str):
        path = (report_assets_dir / asset_path).resolve()
        if not _is_within_directory(path, report_assets_dir) or not path.is_file():
            raise HTTPException(status_code=404, detail="Report asset not found")
        return FileResponse(path)

    @app.get("/structures/{filename}", include_in_schema=False)
    def structure_file(filename: str):
        path = (structures_dir / filename).resolve()
        if not _is_within_directory(path, structures_dir) or not path.exists() or not path.is_file():
            raise HTTPException(status_code=404, detail="Structure not found")
        return FileResponse(path)

    @app.get("/ipi-logo-dark-800.png", include_in_schema=False)
    def logo_dark():
        path = _resolve_brand_asset("ipi-logo-dark-800.png", portal_dir=portal_dir)
        if path is None:
            raise HTTPException(status_code=404, detail="Logo not found")
        return FileResponse(path)

    @app.get("/ipi-logo-light-800.png", include_in_schema=False)
    def logo_light():
        path = _resolve_brand_asset("ipi-logo-light-800.png", portal_dir=portal_dir)
        if path is None:
            raise HTTPException(status_code=404, detail="Logo not found")
        return FileResponse(path)

    @app.get("/api/genes")
    def api_genes():
        entries = [
            {
                key: value
                for key, value in entry.items()
                if key not in {"json_report", "markdown_report", "error"}
            }
            for entry in load_entries()
        ]
        return JSONResponse(entries)

    @app.get("/downloads/agdesign2_portal_index.tsv", response_class=PlainTextResponse, include_in_schema=False)
    def download_index_tsv():
        return PlainTextResponse(portal_index_tsv(load_entries()), media_type="text/tab-separated-values")

    @app.get("/downloads/agdesign2_portal_index.json", response_class=PlainTextResponse, include_in_schema=False)
    def download_index_json():
        return PlainTextResponse(portal_index_json(load_entries()), media_type="application/json")

    @app.get("/downloads/download_manifest.json", response_class=PlainTextResponse, include_in_schema=False)
    def download_manifest():
        return PlainTextResponse(portal_download_manifest(load_entries(), disease_downloads=_available_disease_downloads(batch_dir)), media_type="application/json")

    @app.get("/open_targets_disease_associations.tsv", include_in_schema=False)
    def open_targets_tsv():
        path = batch_dir / "open_targets_disease_associations.tsv"
        if not path.is_file():
            raise HTTPException(status_code=404, detail="Open Targets download not found")
        return FileResponse(path, media_type="text/tab-separated-values")

    @app.get("/open_targets_disease_associations.json", include_in_schema=False)
    def open_targets_json():
        path = batch_dir / "open_targets_disease_associations.json"
        if not path.is_file():
            raise HTTPException(status_code=404, detail="Open Targets download not found")
        return FileResponse(path, media_type="application/json")

    @app.get("/downloads/{citation_file}", response_class=PlainTextResponse, include_in_schema=False)
    def download_citation(citation_file: str):
        content = citation_downloads().get(citation_file)
        if content is None:
            raise HTTPException(status_code=404, detail="Citation file not found")
        media_type = "application/x-research-info-systems" if citation_file.endswith(".ris") else "application/x-bibtex"
        return PlainTextResponse(content, media_type=media_type, headers={"Content-Disposition": f'attachment; filename="{citation_file}"'})

    @app.get("/portal_metadata.json", response_class=PlainTextResponse, include_in_schema=False)
    def portal_metadata():
        return PlainTextResponse(portal_build_metadata(load_entries()), media_type="application/json")

    @app.get("/api/gene/{entry_name}")
    def api_gene(entry_name: str):
        detail_page = f"{entry_name.lower()}.html"
        payload = load_detail_payload(detail_page)
        if payload is None:
            raise HTTPException(status_code=404, detail="Gene report not found")
        return JSONResponse(payload["report"])

    @app.get("/{filename}", include_in_schema=False)
    def root_asset(filename: str):
        if filename not in _PUBLIC_ROOT_ASSETS:
            raise HTTPException(status_code=404, detail="File not found")
        path = portal_dir / filename
        if not path.is_file():
            raise HTTPException(status_code=404, detail="File not found")
        return FileResponse(path)

    return app


def run_dynamic_portal(
    summary_path: str | Path,
    *,
    host: str = "127.0.0.1",
    port: int = 8000,
    refresh_reports: bool = False,
    refresh_stale_reports: bool = False,
    verbose: bool = False,
    enable_complex_portal: bool = False,
) -> None:
    try:
        import uvicorn
    except Exception as exc:  # pragma: no cover - dependency/runtime guard
        raise RuntimeError(
            "Dynamic portal requires uvicorn. Reinstall the project dependencies first."
        ) from exc

    app = create_app(
        summary_path,
        refresh_reports=refresh_reports,
        refresh_stale_reports=refresh_stale_reports,
        verbose=verbose,
        enable_complex_portal=enable_complex_portal,
    )
    uvicorn.run(app, host=host, port=port)


def summarize_dynamic_entries(summary_path: str | Path) -> list[dict[str, Any]]:
    summary_file = Path(summary_path).resolve()
    batch_dir = summary_file.parent
    results = json.loads(summary_file.read_text(encoding="utf-8"))
    disease_index = _load_open_targets_index(batch_dir)
    entries: list[dict[str, Any]] = []
    for item in results:
        report_json_path = _resolve_existing_path(item.get("json_report"), batch_dir)
        report_data: dict[str, Any] | None = None
        if report_json_path is not None:
            try:
                report_data = json.loads(report_json_path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError) as exc:
                raise ValueError(f"Could not read report JSON: {report_json_path}") from exc
        entry = _build_dynamic_index_entry(item, report_data, report_json_path)
        _apply_open_targets_to_entry(entry, disease_index)
        entries.append(entry)
    return entries


def _resolve_brand_asset(name: str, *, portal_dir: Path) -> Path | None:
    for candidate in (portal_dir / name, Path.cwd() / "assets" / "branding" / name, Path.cwd() / name):
        if candidate.exists():
            return candidate.resolve()
    return None


def _build_dynamic_index_entry(
    item: dict[str, Any],
    report_data: dict[str, Any] | None,
    report_json_path: Path | None,
) -> dict[str, Any]:
    entry = dict(item)
    if report_data is None:
        query = item.get("query")
        entry["detail_page"] = None
        entry["gene_symbol"] = item.get("gene_symbol") or query
        entry["protein_name"] = item.get("protein_name") or ""
        entry["family_name"] = item.get("family_name") or ""
        entry["construct_count"] = ""
        entry["cross_reactivity_count"] = ""
        entry["ectodomain"] = ""
        entry["entry_name"] = item.get("resolved_entry_name") or query
        entry["has_alphafold_structure"] = False
        entry["structure_source"] = ""
        entry["has_pdb"] = False
        entry["has_any_structure"] = False
        entry["pubtator_count"] = item.get("pubtator_count")
        entry["pubtator_query_url"] = item.get("pubtator_query_url")
        return entry

    target = report_data.get("target", {})
    ectodomain = report_data.get("ectodomain") or {}
    entry_name = target.get("entry_name") or item.get("resolved_entry_name") or item.get("query")
    assets_dir = None
    if report_json_path is not None:
        assets_dir = report_json_path.parent / f"{str(entry_name).lower()}_report_assets"
    has_assets = any(
        (detail.get("structure_image") or detail.get("quality_plot"))
        for detail in (report_data.get("construct_details") or [])
    )
    if not has_assets and assets_dir is not None and assets_dir.exists():
        has_assets = any(assets_dir.glob("*.png"))
    accession = str(target.get("accession") or "").strip()
    alphafold_path = (
        _resolve_alphafold_artifact(
            accession,
            ".pdb",
            batch_dir=report_json_path.parent if report_json_path is not None else Path.cwd(),
        )
        if accession
        else None
    )
    has_alphafold = alphafold_path is not None
    has_pdb = bool(report_data.get("experimental_constructs"))
    literature = report_data.get("literature_context") or {}
    entry.update(
        {
            "entry_name": entry_name,
            "gene_symbol": target.get("gene_symbol") or item.get("query"),
            "protein_name": target.get("protein_name") or "",
            "target_aliases": " ".join(_target_aliases(target, entry=entry)),
            "family_name": _index_family_name(report_data),
            "construct_count": len(report_data.get("construct_details") or report_data.get("construct_recommendations") or []),
            "cross_reactivity_count": _index_cross_reactivity_count(report_data),
            "ectodomain": _format_region(ectodomain),
            "ectodomain_start": ectodomain.get("start"),
            "ectodomain_end": ectodomain.get("end"),
            "detail_page": f"{_safe_filename(str(entry_name).lower())}.html",
            "has_family_context": bool(report_data.get("family_context")),
            "has_assembly_requirements": bool(report_data.get("assembly_requirements")),
            "has_alphafold_structure": has_alphafold,
            "structure_source": _structure_source_label(alphafold_path) if has_alphafold and alphafold_path is not None else "",
            "has_pdb": has_pdb,
            "has_any_structure": has_alphafold or has_pdb,
            "has_interpro": bool(report_data.get("interpro_annotations")),
            "has_assets": has_assets,
            "pubtator_count": item.get("pubtator_count", literature.get("count")),
            "pubtator_query_url": item.get("pubtator_query_url", literature.get("query_url")),
        }
    )
    if report_json_path is not None:
        entry["json_report"] = str(report_json_path)
    return entry
