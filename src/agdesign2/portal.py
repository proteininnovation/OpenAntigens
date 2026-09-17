from __future__ import annotations

import base64
import csv
from collections import Counter
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import io
import os
import re
import shutil
from time import perf_counter
from html import escape
from pathlib import Path
from typing import Any
from urllib.parse import quote
from urllib.request import urlopen

from .cache import FileCache
from .blast import extract_accession, extract_entry_name
from .af3 import is_af3_artifact
from .config import AnalysisConfig
from .http import HttpClient, _atomic_write
from .pipeline import AntigenAnalyzer
from .sequence_utils import global_align, slice_sequence
from .structure_utils import load_pae_matrix, pae_matches_sequence_length, parse_alphafold_pdb
from .site_docs import PAPER, agent_guide_html, citation_html, citation_short_html, citation_metadata, site_document_files


# accession -> refseq target dict (or None) per ortholog table, keyed by
# (path, mtime_ns) so the multi-MB sequence-bearing TSV is parsed once per build
# instead of re-scanned on every homolog cache miss.
_ORTHOLOG_REFSEQ_INDEX_CACHE: dict[tuple[str, int], dict[str, dict[str, Any] | None]] = {}
_PUBTATOR_COUNT_CACHE: dict[tuple[str, str], dict[str, Any] | None] = {}
_ALPHAFOLD_SEQUENCE_MATCH_CACHE: dict[tuple[str, int, int, str], bool] = {}
_PORTAL_PROCESS_STARTED_UTC = datetime.now(timezone.utc).replace(microsecond=0)
_PORTAL_VERSION_CACHE: str | None = None
_PORTAL_CSS_VERSION = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()[:16]
_PORTAL_JS_VERSION = _PORTAL_CSS_VERSION

_PLAUSIBLE_SCRIPT_URL = "https://plausible.io/js/pa-Mda6DF7h4_8b17ZUlD_9E.js"
_ANALYTICS_EVENT_NAMES = (
    "Download PDB Full",
    "Download PDB Selection",
    "Download PNG",
    "Copy TSV",
    "Copy FASTA",
)
_THREEDMOL_URL = "https://cdn.jsdelivr.net/npm/3dmol@2.5.5/build/3Dmol-min.js"
_THREEDMOL_SHA256 = "f7cc78921ae72e7623e89cdd111434f58c2efddd2ffda1cd212644b406fb8016"
_THREEDMOL_SRI = "sha256-" + base64.b64encode(bytes.fromhex(_THREEDMOL_SHA256)).decode("ascii")
_THREE_TO_ONE = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "ASX": "B",
    "GLX": "Z",
    "SEC": "U",
    "PYL": "O",
    "UNK": "X",
}


def build_portal(
    summary_path: str | Path,
    output_dir: str | Path | None = None,
    *,
    refresh_reports: bool = False,
    refresh_stale_reports: bool = True,
    verbose: bool = False,
    enable_complex_portal: bool = False,
    bundle_vendor_assets: bool = False,
    portal_title: str = "OpenAntigens",
    portal_subtitle: str | None = None,
    portal_intro: str | None = None,
    theme: str = "human",
    include_disease_context: bool = True,
    # When False, skip the live PubTator literature fetch and read pubtator_count
    # from the committed report data instead. Keeps rendering hermetic and
    # deterministic (e.g. the golden test); deploys leave this True.
    fetch_literature: bool = True,
    sibling_link: tuple[str, str] | None = ("OpenAntigens Mouse", "mouse/index.html"),
    detail_sibling_link: tuple[str, str] | None = ("OpenAntigens Mouse", "../mouse/index.html"),
    allow_af3_structures: bool = False,
) -> Path:
    summary_file = Path(summary_path)
    batch_dir = summary_file.parent
    results = json.loads(summary_file.read_text(encoding="utf-8"))
    results = _refresh_reports_if_needed(
        results=results,
        summary_file=summary_file,
        batch_dir=batch_dir,
        refresh_reports=refresh_reports,
        refresh_stale_reports=refresh_stale_reports,
        verbose=verbose,
        enable_complex_portal=enable_complex_portal,
    )
    portal_dir = Path(output_dir) if output_dir is not None else batch_dir / "portal"
    portal_dir.mkdir(parents=True, exist_ok=True)
    reports_dir = portal_dir / "reports"
    reports_dir.mkdir(parents=True, exist_ok=True)
    structures_dir = portal_dir / "structures"
    structures_dir.mkdir(parents=True, exist_ok=True)
    report_assets_dir = portal_dir / "report_assets"
    if report_assets_dir.exists():
        shutil.rmtree(report_assets_dir)
    report_assets_dir.mkdir(parents=True, exist_ok=True)
    report_scripts_dir = portal_dir / "report_scripts"
    report_scripts_dir.mkdir(parents=True, exist_ok=True)
    _copy_brand_assets(portal_dir)
    _write_portal_htaccess(portal_dir)
    if bundle_vendor_assets:
        _ensure_vendor_assets(portal_dir=portal_dir, batch_dir=batch_dir)

    report_entries: list[dict[str, Any]] = []
    disease_index = _load_open_targets_index(batch_dir) if include_disease_context else {}
    (portal_dir / "portal.css").write_text(_portal_css(theme=theme), encoding="utf-8")
    shared_report_viewer_js = ""
    for item in results:
        report_json_path = _resolve_existing_path(item.get("json_report"), batch_dir)
        report_data = None
        if report_json_path is not None:
            try:
                report_data = json.loads(report_json_path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError) as exc:
                raise ValueError(f"Could not read report JSON: {report_json_path}") from exc
        if (
            report_data is not None
            and not allow_af3_structures
            and _report_uses_af3_structure(report_data)
        ):
            raise ValueError(
                f"Public portal build rejected an AF3-derived report: {report_json_path}"
            )
        entry = _build_index_entry(
            item,
            report_data,
            report_json_path,
            fetch_literature=fetch_literature,
            allow_af3_structures=allow_af3_structures,
        )
        if report_data is not None:
            entry["portal_structure_path"] = _copy_portal_structure(
                report_data=report_data,
                structures_dir=structures_dir,
                batch_dir=batch_dir,
                allow_af3_structures=allow_af3_structures,
            )
            if entry["portal_structure_path"]:
                entry["has_alphafold_structure"] = True
                entry["has_any_structure"] = True
                entry["structure_source"] = _structure_source_label(Path(str(entry["portal_structure_path"])))
            _copy_report_assets(
                report_data=report_data,
                report_json_path=report_json_path,
                report_assets_dir=report_assets_dir,
                batch_dir=batch_dir,
            )
        if include_disease_context:
            _apply_open_targets_to_entry(entry, disease_index)
        report_entries.append(entry)
        if report_data is not None and entry.get("detail_page"):
            detail_html = render_detail_page(
                entry,
                report_data,
                batch_dir=batch_dir,
                portal_dir=portal_dir,
                portal_title=portal_title,
                include_disease_context=include_disease_context,
                sibling_link=detail_sibling_link,
                allow_af3_structures=allow_af3_structures,
            )
            detail_script_name = _safe_filename(str(Path(str(entry["detail_page"])).with_suffix("").name)) + ".js"
            detail_html, detail_js = _externalize_plain_scripts(
                detail_html,
                script_src=f"../report_scripts/{detail_script_name}?v={_PORTAL_JS_VERSION}",
            )
            detail_js, shared_js = _split_report_viewer_runtime(detail_js)
            if shared_js:
                shared_report_viewer_js = shared_js
                detail_html = detail_html.replace(
                    f'<script src="../report_scripts/{detail_script_name}?v={_PORTAL_JS_VERSION}"></script>',
                    f'<script src="../report-viewer.js?v={_PORTAL_JS_VERSION}"></script>\n    <script src="../report_scripts/{detail_script_name}?v={_PORTAL_JS_VERSION}"></script>',
                    1,
                )
            if detail_js:
                _write_text_atomic(report_scripts_dir / detail_script_name, detail_js)
            (reports_dir / str(entry["detail_page"])).write_text(detail_html, encoding="utf-8")
    _write_portal_download_files(
        portal_dir,
        report_entries,
        batch_dir=batch_dir,
        include_disease_context=include_disease_context,
    )
    index_html, index_js = _externalize_inline_script(
        render_index_page(
            report_entries,
            portal_title=portal_title,
            portal_subtitle=portal_subtitle,
            portal_intro=portal_intro,
            include_disease_context=include_disease_context,
            sibling_link=sibling_link,
        ),
        script_src=f"portal-index.js?v={_PORTAL_JS_VERSION}",
    )
    _write_text_atomic(portal_dir / "index.html", index_html)
    _write_text_atomic(
        portal_dir / "portal-index-data.js",
        _render_index_data_js(report_entries, include_disease_context=include_disease_context),
    )
    # The per-target disease-name haystack is large (the bulk of the old index
    # payload) but only needed for the disease filter, not first paint. Ship it
    # as a separate deferred script so the protein list renders from the much
    # smaller row payload first.
    if include_disease_context:
        _write_text_atomic(
            portal_dir / "portal-disease-index.js",
            _render_disease_index_js(report_entries),
        )
    _write_text_atomic(portal_dir / "portal-index.js", index_js)
    if shared_report_viewer_js:
        _write_text_atomic(portal_dir / "report-viewer.js", shared_report_viewer_js)
    _write_text_atomic(
        portal_dir / "help.html",
        render_help_page(portal_title=portal_title, include_disease_context=include_disease_context),
    )
    _write_text_atomic(portal_dir / "builder.html", render_builder_page(portal_title=portal_title))
    _write_text_atomic(portal_dir / "constructs.html", render_constructs_page(portal_title=portal_title))
    _write_text_atomic(portal_dir / "methods.html", render_methods_page(portal_title=portal_title))
    _write_text_atomic(
        portal_dir / "downloads.html",
        render_downloads_page(
            report_entries,
            include_disease_context=include_disease_context,
            portal_title=portal_title,
            disease_downloads=_available_disease_downloads(portal_dir),
        ),
    )
    _write_text_atomic(portal_dir / "calculator.html", render_calculator_page(portal_title=portal_title))
    _write_text_atomic(portal_dir / "terms.html", render_terms_page(portal_title=portal_title))
    _write_text_atomic(portal_dir / "privacy.html", render_privacy_page(portal_title=portal_title))
    _write_site_documents(portal_dir, portal_title=portal_title)
    _write_text_atomic(portal_dir / "portal_metadata.json", portal_build_metadata(report_entries, portal_title=portal_title))
    _prune_generated_target_files(portal_dir, report_entries)
    _version_portal_assets(portal_dir)
    return portal_dir / "index.html"


def _version_portal_assets(portal_dir: Path) -> None:
    """Version local CSS/JS references by the bytes actually served."""
    from urllib.parse import urlsplit

    versions: dict[Path, str] = {}
    for page in [*portal_dir.glob("*.html"), *(portal_dir / "reports").glob("*.html")]:
        def replace(match: re.Match) -> str:
            url = urlsplit(match.group("url"))
            if url.scheme or url.netloc or not url.path.endswith((".js", ".css")):
                return match.group(0)
            asset = (page.parent / url.path).resolve()
            if not asset.is_relative_to(portal_dir.resolve()) or not asset.is_file():
                return match.group(0)
            if asset not in versions:
                versions[asset] = hashlib.sha256(asset.read_bytes()).hexdigest()[:16]
            return f'{match.group("attribute")}="{url.path}?v={versions[asset]}"'

        html = re.sub(r'(?P<attribute>src|href)="(?P<url>[^"]+)"', replace, page.read_text())
        _write_text_atomic(page, html)


def _prune_generated_target_files(portal_dir: Path, entries: list[dict[str, Any]]) -> None:
    reports = {Path(str(entry["detail_page"])).name for entry in entries}
    scripts = {Path(name).with_suffix(".js").name for name in reports}
    structures = {Path(str(entry["portal_structure_path"])).stem for entry in entries if entry.get("portal_structure_path")}
    for directory, suffix, keep in (("reports", ".html", reports), ("report_scripts", ".js", scripts)):
        for path in (portal_dir / directory).glob("*" + suffix):
            if path.name not in keep:
                path.unlink()
    for path in (portal_dir / "structures").glob("*"):
        if path.name.endswith((".pdb", ".meta.json")):
            stem = path.name.removesuffix(".meta.json").removesuffix(".pdb")
            if stem not in structures:
                path.unlink()


def _externalize_inline_script(html: str, *, script_src: str) -> tuple[str, str]:
    """Move the final inline script block into a cacheable sibling JS file."""

    start = html.rfind("  <script>")
    if start < 0:
        return html, ""
    content_start = start + len("  <script>")
    end = html.find("  </script>", content_start)
    if end < 0:
        return html, ""
    script = html[content_start:end].strip() + "\n"
    replacement = f'  <script src="{escape(script_src)}" defer></script>'
    return html[:start] + replacement + html[end + len("  </script>") :], script


def _externalize_plain_scripts(html: str, *, script_src: str) -> tuple[str, str]:
    """Move plain inline JavaScript blocks into one external script.

    Script tags with attributes are intentionally left untouched so JSON payloads
    and third-party library includes stay in their original positions.
    """

    pattern = re.compile(r"(?P<indent>[ \t]*)<script>(?P<body>.*?)</script>", re.DOTALL)
    scripts: list[str] = []
    pieces: list[str] = []
    cursor = 0
    inserted = False
    for match in pattern.finditer(html):
        body = match.group("body").strip()
        if not body:
            continue
        scripts.append(body)
        pieces.append(html[cursor : match.start()])
        if not inserted:
            indent = match.group("indent")
            pieces.append(f'{indent}<script src="{escape(script_src)}"></script>')
            inserted = True
        cursor = match.end()
    if not scripts:
        return html, ""
    pieces.append(html[cursor:])
    return "".join(pieces), "\n\n".join(scripts) + "\n"


def _split_report_viewer_runtime(script: str) -> tuple[str, str]:
    """Split a generated report viewer script into data and shared runtime JS."""

    marker = "        const viewerEl = document.getElementById('structureViewer');"
    if marker not in script:
        return script, ""
    stripped = script.strip()
    opener = "(() => {"
    closer = "      })();"
    if not stripped.startswith(opener) or not stripped.endswith(closer):
        return script, ""
    body = stripped[len(opener) : -len(closer)]
    marker_index = body.find(marker)
    if marker_index < 0:
        return script, ""
    data_block = body[:marker_index].rstrip()
    runtime_block = body[marker_index:].rstrip()
    data_script = f"""(() => {{{data_block}
        window.OpenAntigenReportViewer?.init({{
          sequence,
          pdbText,
          ectodomain,
          designRegionLabel,
          designRegionLower,
          plotScopeFullLength,
          targetInfo,
          viewerData
        }});
      }})();
"""
    shared_script = f"""window.OpenAntigenReportViewer = window.OpenAntigenReportViewer || {{}};
window.OpenAntigenReportViewer.init = function(payload) {{
        const sequence = payload.sequence || '';
        const pdbText = payload.pdbText || '';
        const ectodomain = payload.ectodomain || {{}};
        const designRegionLabel = payload.designRegionLabel || 'Design Region';
        const designRegionLower = payload.designRegionLower || 'design region';
        const plotScopeFullLength = Boolean(payload.plotScopeFullLength);
        const targetInfo = payload.targetInfo || {{}};
        const viewerData = payload.viewerData || {{}};
{runtime_block}
}};
"""
    return data_script, shared_script


def _js_json(value: Any) -> str:
    """``json.dumps`` safe for embedding inside an executable inline ``<script>``.

    ``json.dumps`` leaves ``<``, ``>`` and ``&`` unescaped, so a ``</script>``
    or ``<!--`` substring in a string value can terminate the script element and
    let following bytes be parsed as HTML. These characters only occur inside
    JSON string tokens, so replacing them with their ``\\u`` escapes keeps the
    decoded value byte-for-byte identical while making the output safe in an
    HTML context. ``ensure_ascii`` already escapes U+2028/U+2029.
    """
    return (
        json.dumps(value)
        .replace("<", "\\u003c")
        .replace(">", "\\u003e")
        .replace("&", "\\u0026")
    )


def _render_index_data_js(entries: list[dict[str, Any]], *, include_disease_context: bool = True) -> str:
    payload = [_index_row_payload(item, include_disease_context=include_disease_context) for item in entries]
    return "window.OpenAntigenIndexRows = " + json.dumps(payload, separators=(",", ":")) + ";\n"


def _render_disease_index_js(entries: list[dict[str, Any]]) -> str:
    """Disease-name filter haystack keyed by lower-cased entry name.

    Loaded as a deferred script so the (large) disease data does not block the
    protein list from rendering. Only targets with disease names are included.
    """
    index = {}
    for item in entries:
        haystack = _entry_disease_haystack(item)
        if haystack:
            index[str(item.get("entry_name") or "").lower()] = haystack
    return "window.OpenAntigenDiseaseIndex = " + json.dumps(index, separators=(",", ":")) + ";\n"


def _index_row_payload(item: dict[str, Any], *, include_disease_context: bool = True) -> dict[str, Any]:
    detail_label = str(item.get("entry_name") or item.get("query") or "view")
    if item.get("detail_page"):
        detail_path = str(item["detail_page"])
        detail_href = detail_path if detail_path.startswith("reports/") else f"reports/{detail_path}"
        detail_link = f"<a href=\"{escape(detail_href)}\">{escape(detail_label)}</a>"
    else:
        detail_link = f'<span class="muted">{escape(detail_label)}</span>'
    dataset = {
        "query": str(item.get("query") or "").lower(),
        "gene": str(item.get("gene_symbol") or "").lower(),
        "status": str(item.get("status") or ""),
        "index": str(item.get("batch_index") or 0),
        "constructs": str(item.get("construct_count") or 0),
        "cross": str(item.get("cross_reactivity_count") or 0),
        "entry": str(item.get("entry_name") or "").lower(),
        "protein": str(item.get("protein_name") or "").lower(),
        "family": str(item.get("family_name") or "").lower(),
        "aliases": str(item.get("target_aliases") or "").lower(),
        "track": str(item.get("topology_bucket") or "").lower(),
        # `diseases` (the disease-name filter haystack) is shipped separately in
        # portal-disease-index.js and merged in client-side after first paint.
        "diseaseScore": str(item.get("top_disease_score") or 0),
        "pubtator": str(item.get("pubtator_count") or 0),
        "ectoStart": str(item.get("ectodomain_start") or 0),
        "ectoEnd": str(item.get("ectodomain_end") or 0),
        "hasFamily": "1" if item.get("has_family_context") else "0",
        "hasAssembly": "1" if item.get("has_assembly_requirements") else "0",
        "hasStructure": "1" if item.get("has_any_structure") else "0",
        "hasAlphafold": "1" if item.get("has_alphafold_structure") else "0",
        "hasPdb": "1" if item.get("has_pdb") else "0",
        "hasInterpro": "1" if item.get("has_interpro") else "0",
        "hasAssets": "1" if item.get("has_assets") else "0",
    }
    cells = [
        escape(str(item.get("batch_index", ""))),
        detail_link,
        escape(str(item.get("gene_symbol") or "")),
        escape(str(item.get("protein_name") or "")),
        escape(str(item.get("family_name") or "")),
        escape(str(item.get("topology_bucket") or "")),
        _render_pubtator_link(item.get("pubtator_count"), item.get("pubtator_query_url")),
        escape(str(item.get("ectodomain") or "")),
        _structure_badge(bool(item.get("has_alphafold_structure")), label="AF"),
        _structure_badge(bool(item.get("has_pdb")), label="PDB"),
        escape(str(item.get("construct_count") or "")),
        escape(str(item.get("cross_reactivity_count") or "")),
    ]
    if include_disease_context:
        cells.insert(5, _render_top_disease(item))
    return {"dataset": dataset, "cells": cells}


def _refresh_reports_if_needed(
    *,
    results: list[dict[str, Any]],
    summary_file: Path,
    batch_dir: Path,
    refresh_reports: bool,
    refresh_stale_reports: bool,
    verbose: bool,
    enable_complex_portal: bool,
) -> list[dict[str, Any]]:
    if not results:
        return results
    if not refresh_reports and not refresh_stale_reports:
        return results

    from .module_refresh import _module_refresh_config

    config = _module_refresh_config(batch_dir=batch_dir, verbose=verbose, enable_complex_portal=enable_complex_portal)
    dependencies = _portal_refresh_dependencies(config)
    analyzer: AntigenAnalyzer | None = None
    refreshed = False

    for item in results:
        if not _should_refresh_report(
            item=item,
            batch_dir=batch_dir,
            refresh_reports=refresh_reports,
            refresh_stale_reports=refresh_stale_reports,
            dependencies=dependencies,
        ):
            continue
        if analyzer is None:
            analyzer = AntigenAnalyzer(config=config)
        query = str(item.get("query") or "").strip()
        if not query:
            continue
        if verbose:
            print(f"[agdesign2-portal] refreshing report for {query}", flush=True)
        started = perf_counter()
        try:
            report = analyzer.analyze_target(query, output_dir=batch_dir)
            stem = f"{report.target.entry_name.lower()}_report"
            item["status"] = "ok"
            item["resolved_entry_name"] = report.target.entry_name
            item["json_report"] = str(batch_dir / f"{stem}.json")
            item["markdown_report"] = str(batch_dir / f"{stem}.md")
            item["error"] = None
        except Exception as exc:
            item["status"] = "error"
            item["error"] = str(exc)
        item["duration_seconds"] = round(perf_counter() - started, 2)
        refreshed = True

    if refreshed:
        _write_text_atomic(summary_file, json.dumps(results, indent=2) + "\n")
    return results


def _portal_refresh_dependencies(config: AnalysisConfig) -> list[Path]:
    dependencies: list[Path] = []
    for path in (
        config.precomputed_ortholog_table_path,
        config.precomputed_family_alignment_index_path,
    ):
        if path.exists():
            dependencies.append(path)
    return dependencies


def _should_refresh_report(
    *,
    item: dict[str, Any],
    batch_dir: Path,
    refresh_reports: bool,
    refresh_stale_reports: bool,
    dependencies: list[Path],
) -> bool:
    if refresh_reports:
        return True
    if not refresh_stale_reports:
        return False
    if item.get("status") == "error":
        return True
    report_json_path = _resolve_existing_path(item.get("json_report"), batch_dir)
    report_md_path = _resolve_existing_path(item.get("markdown_report"), batch_dir)
    if report_json_path is None or report_md_path is None:
        return True
    if not report_json_path.exists() or not report_md_path.exists():
        return True
    if not dependencies:
        return False
    report_mtime = min(report_json_path.stat().st_mtime, report_md_path.stat().st_mtime)
    return any(path.stat().st_mtime > report_mtime for path in dependencies if path.exists())


def _write_text_atomic(path: Path, text: str) -> None:
    _atomic_write(path, text.encode("utf-8"))


def _get_pubtator_literature_info(gene_symbol: str | None, *, cache_dir: Path) -> dict[str, Any] | None:
    symbol = str(gene_symbol or "").strip().upper()
    if not symbol:
        return None
    cache_key = (str(cache_dir.resolve()), symbol)
    if cache_key in _PUBTATOR_COUNT_CACHE:
        return _PUBTATOR_COUNT_CACHE[cache_key]
    query_text = f"@GENE_{symbol}"
    query_url = f"https://www.ncbi.nlm.nih.gov/research/pubtator3/docsum?text={quote(query_text)}"
    api_url = f"https://www.ncbi.nlm.nih.gov/research/pubtator3-api/search/?text={quote(query_text)}"
    info: dict[str, Any] = {
        "gene_symbol": symbol,
        "query_text": query_text,
        "query_url": query_url,
        "count": None,
    }
    cache = FileCache(cache_dir.resolve())
    cache_path = cache.path_for("pubtator_search", api_url, ".json")
    http = HttpClient(cache=cache, timeout=20)
    for attempt in range(2):
        try:
            payload = http.fetch_json(api_url, cache_namespace="pubtator_search", cache_key=api_url)
            count = _pubtator_payload_count(payload)
            if count is not None:
                info["count"] = count
                break
            if attempt == 0 and cache_path.exists():
                cache_path.unlink()
                continue
        except Exception:
            if attempt == 0 and cache_path.exists():
                try:
                    cache_path.unlink()
                except OSError:
                    pass
                continue
            break
    _PUBTATOR_COUNT_CACHE[cache_key] = info
    return info


def _pubtator_payload_count(payload: Any) -> int | None:
    if not isinstance(payload, dict):
        return None
    count = payload.get("count")
    if count is None:
        return None
    try:
        return int(count)
    except (TypeError, ValueError):
        return None


def _coerce_pubtator_count(count: Any) -> int | None:
    if isinstance(count, bool):
        return None
    if isinstance(count, int):
        return count
    if isinstance(count, str):
        normalized = count.strip().replace(",", "")
        if normalized.isdigit():
            return int(normalized)
    return None


def _pubtator_count_or_existing(pubtator_info: dict[str, Any] | None, item: dict[str, Any]) -> int | None:
    count = _coerce_pubtator_count(pubtator_info.get("count") if pubtator_info else None)
    if count is not None:
        return count
    return _coerce_pubtator_count(item.get("pubtator_count"))


_PORTAL_HTACCESS_COMMON = """\
<IfModule mod_deflate.c>
  AddOutputFilterByType DEFLATE text/html text/css text/plain text/javascript application/javascript image/svg+xml
</IfModule>
<IfModule mod_expires.c>
  ExpiresActive On
  ExpiresByType application/javascript "access plus 7 days"
  ExpiresByType text/javascript "access plus 7 days"
  ExpiresByType text/css "access plus 7 days"
  ExpiresByType image/png "access plus 30 days"
  ExpiresByType text/html "access plus 0 seconds"
</IfModule>
"""


_PORTAL_HTACCESS = """\
# OpenAntigens static portal — Apache hosting hints for plain assets.
""" + _PORTAL_HTACCESS_COMMON


_PRECOMPRESSED_PORTAL_HTACCESS = """\
# OpenAntigens static portal — Apache hosting hints for precompressed assets.
# JS/CSS/SVG retain their public filenames but contain gzip bytes, avoiding any
# dependency on mod_rewrite while keeping the deployment within its disk quota.
<IfModule mod_headers.c>
  <FilesMatch "\\.js$">
    ForceType text/javascript
    Header set Content-Encoding gzip
    Header append Vary Accept-Encoding
  </FilesMatch>
  <FilesMatch "\\.css$">
    ForceType text/css
    Header set Content-Encoding gzip
    Header append Vary Accept-Encoding
  </FilesMatch>
  <FilesMatch "\\.svg$">
    ForceType image/svg+xml
    Header set Content-Encoding gzip
    Header append Vary Accept-Encoding
  </FilesMatch>
</IfModule>
# The files already contain gzip bytes; never let mod_deflate encode them again.
<FilesMatch "\\.(?:js|css|svg)$">
  SetEnv no-gzip 1
</FilesMatch>
""" + _PORTAL_HTACCESS_COMMON


def _write_portal_htaccess(portal_dir: Path, *, precompressed: bool = False) -> None:
    """Write a conservative Apache .htaccess enabling gzip for text assets.

    The published portal is large (the index data blob alone is ~19 MB; gzip
    takes it to ~2.5 MB). Directives are <IfModule>-guarded so an unsupported
    module is a no-op rather than a server error. Verify on the target server
    because AllowOverride may restrict .htaccess.
    """
    content = _PRECOMPRESSED_PORTAL_HTACCESS if precompressed else _PORTAL_HTACCESS
    (portal_dir / ".htaccess").write_text(content, encoding="utf-8")


def _copy_brand_assets(portal_dir: Path) -> None:
    repo_brand_dir = Path(__file__).resolve().parents[2] / "assets" / "branding"
    for asset_name in (
        "ipi-logo-dark-800.png",
        "ipi-logo-light-800.png",
        "favicon.ico",
        "favicon-32.png",
        "apple-touch-icon.png",
        "icon-192.png",
        "icon-512.png",
        "site.webmanifest",
    ):
        sources = (
            Path.cwd() / "assets" / "branding" / asset_name,
            repo_brand_dir / asset_name,
            Path.cwd() / asset_name,
        )
        source = next((candidate for candidate in sources if candidate.exists()), None)
        if source is not None:
            destination = portal_dir / asset_name
            if not destination.exists() or source.stat().st_mtime > destination.stat().st_mtime:
                shutil.copyfile(source, destination)


def _portal_favicon_links(*, prefix: str = "") -> str:
    return f"""  <link rel="icon" href="{escape(prefix)}favicon.ico" sizes="any">
  <link rel="icon" type="image/png" sizes="32x32" href="{escape(prefix)}favicon-32.png?v={_PORTAL_CSS_VERSION}">
  <link rel="apple-touch-icon" sizes="180x180" href="{escape(prefix)}apple-touch-icon.png?v={_PORTAL_CSS_VERSION}">
  <link rel="manifest" href="{escape(prefix)}site.webmanifest?v={_PORTAL_CSS_VERSION}">"""


def _portal_analytics_head() -> str:
    allowed_events = json.dumps(list(_ANALYTICS_EVENT_NAMES), separators=(",", ":"))
    return f"""  <!-- Privacy-friendly analytics by Plausible -->
  <script async src="{_PLAUSIBLE_SCRIPT_URL}"></script>
  <script data-openantigens-analytics>
    window.plausible = window.plausible || function () {{
      (plausible.q = plausible.q || []).push(arguments);
    }};
    plausible.init = plausible.init || function (options) {{
      plausible.o = options || {{}};
    }};
    plausible.init({{
      fileDownloads: {{ fileExtensions: ["tsv", "json"] }},
      formSubmissions: false,
      outboundLinks: false
    }});
    const openAntigensAnalyticsEvents = new Set({allowed_events});
    window.trackOpenAntigensEvent = function (name) {{
      if (!openAntigensAnalyticsEvents.has(name)) {{
        throw new Error(`Unknown OpenAntigens analytics event: ${{name}}`);
      }}
      window.plausible(name);
    }};
  </script>"""


def _ensure_vendor_assets(*, portal_dir: Path, batch_dir: Path) -> None:
    """Copy or download browser libraries required for a standalone snapshot."""

    destination = portal_dir / "3Dmol-min.js"
    if destination.exists():
        if _is_valid_threedmol_asset(destination):
            return
        raise RuntimeError(f"Existing 3Dmol-min.js did not match the pinned SHA-256: {destination}")
    candidates = [
        Path.cwd() / "assets" / "vendor" / "3Dmol-min.js",
        batch_dir / "vendor" / "3Dmol-min.js",
        batch_dir.parent / "vendor" / "3Dmol-min.js",
        batch_dir.parent.parent / "data" / "vendor" / "3Dmol-min.js",
    ]
    for source in candidates:
        if not source.exists():
            continue
        if not _is_valid_threedmol_asset(source):
            raise RuntimeError(f"Local 3Dmol-min.js did not match the pinned SHA-256: {source}")
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)
        return
    try:
        with urlopen(_THREEDMOL_URL, timeout=30) as response:  # nosec B310
            payload = response.read()
    except Exception as exc:
        raise RuntimeError("Could not bundle 3Dmol-min.js for a self-contained portal snapshot") from exc
    if hashlib.sha256(payload).hexdigest() != _THREEDMOL_SHA256:
        raise RuntimeError("Downloaded 3Dmol-min.js did not match the pinned SHA-256")
    destination.write_bytes(payload)


def _is_valid_threedmol_asset(path: Path) -> bool:
    return path.is_file() and hashlib.sha256(path.read_bytes()).hexdigest() == _THREEDMOL_SHA256


def _portal_nav(*, prefix: str = "", active: str = "", sibling_link: tuple[str, str] | None = None) -> str:
    links = (
        ("index", "Browse", "index.html"),
        ("builder", "Builder", "builder.html"),
        ("constructs", "Constructs", "constructs.html"),
        ("calculator", "Calculator", "calculator.html"),
        ("help", "Help", "help.html"),
        ("methods", "Methods", "methods.html"),
        ("downloads", "Downloads", "downloads.html"),
        ("terms", "Terms", "terms.html"),
        ("agents", "For AI agents", "agent-guide.html"),
    )
    rendered = "".join(
        '<a class="{active}" href="{href}">{label}</a>'.format(
            active="active" if key == active else "",
            href=escape(f"{prefix}{href}"),
            label=escape(label),
        )
        for key, label, href in links
    )
    return "<nav class=\"site-nav\" aria-label=\"Portal navigation\">" + rendered + "</nav>"


def _portal_species_switch(portal_title: str, sibling_link: tuple[str, str] | None) -> str:
    if sibling_link is None:
        return ""
    label, href = sibling_link
    current_species = "Mouse" if "mouse" in portal_title.lower() else "Human"
    target_species = "Mouse" if "mouse" in label.lower() else "Human"
    return f"""
            <div class="species-switcher" aria-label="Species portal">
              <span class="species-current">{escape(current_species)}</span>
              <a class="species-link" href="{escape(href)}" aria-label="{escape(label)}">{escape(target_species)}</a>
            </div>
"""


def _portal_brand_header(
    *,
    prefix: str = "",
    active: str = "",
    chip_html: str = "",
    portal_title: str = "OpenAntigens",
    sibling_link: tuple[str, str] | None = None,
) -> str:
    chip = chip_html or '<div class="hero-chip">Antigen constructs and supporting evidence</div>'
    species_switch = _portal_species_switch(portal_title, sibling_link)
    return f"""
      <div class="brand-bar">
        <div class="brand-lockup">
          <a class="brand-home-link" href="{escape(prefix)}index.html" aria-label="OpenAntigens home" style="color:#fff;text-decoration:none;display:inline-flex;align-items:center;">
            <img class="brand-logo" src="{escape(prefix)}ipi-logo-light-800.png" alt="Institute for Protein Innovation">
          </a>
          <div class="brand-copy">
            <a class="brand-kicker" href="https://proteininnovation.org/" target="_blank" rel="noreferrer" style="color:#fff;text-decoration:none;">Institute for Protein Innovation</a>
            <a class="brand-product" href="{escape(prefix)}index.html" style="color:#fff;text-decoration:none;">{escape(portal_title)}</a>
          </div>
        </div>
        <div class="brand-actions">
          {_portal_nav(prefix=prefix, active=active)}
          <div class="brand-secondary-actions">
          {species_switch}
          {chip}
          </div>
        </div>
      </div>
"""


def _portal_footer(*, prefix: str = "") -> str:
    version = _openantigen_version()
    build_date = _portal_build_date()
    return f"""
    <footer class="site-footer">
      <div class="footer-grid">
        <div class="footer-brand">
          <a class="footer-logo-link" href="https://proteininnovation.org/" target="_blank" rel="noreferrer" aria-label="Institute for Protein Innovation website">
            <img class="footer-logo" src="{escape(prefix)}ipi-logo-dark-800.png" alt="Institute for Protein Innovation">
          </a>
          <div>
            <strong>OpenAntigens</strong>
            <p>Antigen construct design for cell-surface and secreted proteins.</p>
          </div>
        </div>
        <nav class="footer-links" aria-label="Footer navigation">
          <a href="{escape(prefix)}builder.html">Builder guide</a>
          <a href="{escape(prefix)}constructs.html">Construct methodology</a>
          <a href="{escape(prefix)}methods.html">Methods</a>
          <a href="{escape(prefix)}downloads.html">Downloads</a>
          <a href="{escape(prefix)}terms.html">Terms</a>
          <a href="{escape(prefix)}privacy.html">Privacy</a>
          <a href="{escape(prefix)}agent-guide.html">For AI agents</a>
          <a href="{escape(prefix)}llms.txt">llms.txt</a>
        </nav>
      </div>
      <p class="footer-citation">Using OpenAntigens in your research? Cite
        {citation_short_html()} (preprint).
        <a href="{escape(prefix)}help.html#cite-openantigens">How to cite</a>
      </p>
      <div class="footer-meta">
        <a href="https://proteininnovation.org/" target="_blank" rel="noreferrer">Institute for Protein Innovation</a>
        <span>OpenAntigens v{escape(version)}</span>
        <span>Release snapshot {escape(build_date)} UTC</span>
        <span>Created by Andre A. R. Teixeira</span>
        <a href="mailto:andre.teixeira@proteininnovation.org">andre.teixeira@proteininnovation.org</a>
      </div>
</footer>
"""


def _openantigen_version() -> str:
    global _PORTAL_VERSION_CACHE
    if _PORTAL_VERSION_CACHE is not None:
        return _PORTAL_VERSION_CACHE
    try:
        _PORTAL_VERSION_CACHE = importlib.metadata.version("agdesign2")
    except importlib.metadata.PackageNotFoundError:
        _PORTAL_VERSION_CACHE = _version_from_pyproject()
    return _PORTAL_VERSION_CACHE


def _version_from_pyproject() -> str:
    pyproject = Path.cwd() / "pyproject.toml"
    if not pyproject.is_file():
        return "unknown"
    try:
        for line in pyproject.read_text(encoding="utf-8").splitlines():
            stripped = line.strip()
            if stripped.startswith("version"):
                return stripped.split("=", 1)[1].strip().strip('"').strip("'")
    except Exception:
        pass
    return "unknown"


def _entry_name_base(entry_name: str) -> str:
    text = str(entry_name or "").strip()
    if "_" not in text:
        return text
    base, species = text.rsplit("_", 1)
    if species.upper() in {"HUMAN", "MOUSE", "MACFA"}:
        return base
    return text


def _portal_build_date() -> str:
    return os.environ.get("OPENANTIGEN_BUILD_DATE") or _PORTAL_PROCESS_STARTED_UTC.date().isoformat()


def _topology_counts(entries: list[dict[str, Any]]) -> dict[str, int]:
    counts = {"secreted": 0, "single_pass": 0, "multipass": 0, "gpi": 0, "other": 0}
    for item in entries:
        raw = str(item.get("topology_bucket") or item.get("track") or "").strip().lower()
        normalized = raw.replace("_", "-").replace(" ", "-")
        if "secreted" in normalized:
            counts["secreted"] += 1
        elif "single-pass" in normalized or normalized == "singlepass":
            counts["single_pass"] += 1
        elif "multi-pass" in normalized or "multipass" in normalized:
            counts["multipass"] += 1
        elif "gpi" in normalized:
            counts["gpi"] += 1
        else:
            counts["other"] += 1
    return counts


def _citation_download_links() -> str:
    return '<a class="inline-link" href="downloads/openantigens.bib" download>BibTeX</a> · <a class="inline-link" href="downloads/openantigens.ris" download>RIS</a>'


def _citation_section() -> str:
    return f"""
    <section class="card doc-card" id="cite-openantigens">
      <h2>How to cite OpenAntigens</h2>
      <p>If you use OpenAntigens in research, please cite:</p>
      <p class="paper-citation">{citation_html()}</p>
      <p class="section-note">Preprint, version {escape(PAPER['version'])}; posted {escape(PAPER['date'])}. Not peer reviewed.</p>
      <p>Import the reference: {_citation_download_links()}</p>
      <p>Also record the portal build date, software version, report or downloaded artifact URL, and your access date. These identify the release you used; the paper describes the database.</p>
      <p>Cite underlying databases or primary studies where you use their evidence. See <a class="inline-link" href="terms.html">Terms</a> for data attribution and reuse requirements.</p>
    </section>
"""


def _agent_guidance_panel(*, portal_title: str = "OpenAntigens", guide: bool = False) -> str:
    site = "https://openantigens.org/mouse/" if portal_title == "OpenAntigens Mouse" else "https://openantigens.org/"
    prompt = (
        "Before using OpenAntigens, read https://openantigens.org/llms.txt. "
        f"Start at {site}index.html. Confirm each target's species and accession, preserve its residue numbering, "
        "distinguish predicted from experimental evidence, and cite the inspected report and release. "
        "For research use, include the OpenAntigens preprint citation given in the guide.\n\nMy question: "
    )
    label, href = ("Open plain text /llms.txt", "llms.txt") if guide else ("Read the agent guide", "agent-guide.html")
    return f"""
    <section class="agent-banner" aria-labelledby="agent-banner-title">
      <div><p class="agent-eyebrow">For AI assistants</p>
        <h2 id="agent-banner-title">Using Claude, ChatGPT, or another agent?</h2>
        <p>Start with the guide for search, downloads, residue numbering, and evidence interpretation.</p></div>
      <div class="agent-actions"><a class="agent-guide-link" href="{href}">{label}</a>
        <button id="agent-copy" class="agent-copy" type="button">Copy AI prompt</button></div>
      <details id="agent-prompt" class="agent-prompt"><summary>View prompt</summary>
        <label for="agent-prompt-text">Paste this into your assistant and add your question.</label>
        <textarea id="agent-prompt-text" rows="5" readonly>{escape(prompt)}</textarea>
      </details>
      <span id="agent-copy-status" class="agent-copy-status" role="status" aria-live="polite"></span>
    </section>
    <script src="agent-guide.js" defer></script>
"""


def render_agent_guide_page(*, portal_title: str = "OpenAntigens") -> str:
    human_prefix = "../" if portal_title == "OpenAntigens Mouse" else ""
    content = agent_guide_html().replace('href="https://openantigens.org/"', f'href="{human_prefix}index.html"')
    content = content.replace('href="https://openantigens.org/', f'href="{human_prefix}')
    headings = re.findall(r"<h2>([^<]+)</h2>", content)
    links = []
    for heading in headings:
        anchor = re.sub(r"[^a-z0-9]+", "-", heading.lower()).strip("-")
        content = content.replace(f"<h2>{heading}</h2>", f'<h2 id="{anchor}">{heading}</h2>')
        links.append(f'<a href="#{anchor}">{heading}</a>')
    return _render_info_page(
        title="Agent guide", active="agents", eyebrow="For AI assistants",
        heading="Use OpenAntigens with an AI assistant.",
        intro="Find the right record, interpret its evidence, and keep the source with your answer.",
        portal_title=portal_title,
        body_html=_agent_guidance_panel(portal_title=portal_title, guide=True) + f"""
        <div class="guide-layout"><nav class="guide-toc" aria-label="Guide contents">
          <strong>In this guide</strong>{''.join(links)}</nav>
          <article class="card guide-content">{content}</article></div>""",
    )


def _write_site_documents(portal_dir: Path, *, portal_title: str = "OpenAntigens") -> None:
    _write_text_atomic(portal_dir / "agent-guide.html", render_agent_guide_page(portal_title=portal_title))
    for name, content in site_document_files().items():
        destination = portal_dir / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        _write_text_atomic(destination, content)


def render_index_page(
    entries: list[dict[str, Any]],
    *,
    portal_title: str = "OpenAntigens",
    portal_subtitle: str | None = None,
    portal_intro: str | None = None,
    include_disease_context: bool = True,
    sibling_link: tuple[str, str] | None = None,
) -> str:
    total = len(entries)
    ready_entries = [item for item in entries if item.get("status") == "ok"]
    ready_total = len(ready_entries)
    alphafold_count = sum(1 for item in entries if item.get("has_alphafold_structure"))
    alphafold_missing = total - alphafold_count
    topology_counts = _topology_counts(ready_entries)

    subtitle = portal_subtitle or "Antigen constructs for human cell-surface and secreted proteins."
    intro = portal_intro or "Open a protein report to compare proposed construct boundaries and inspect sequence and structural evidence. Export selected sequences as FASTA or TSV."
    disease_controls = (
        """
      <input id="diseaseSearchBox" type="search" list="diseaseSuggestions" placeholder="Filter by disease" aria-label="Filter by disease">
      <datalist id="diseaseSuggestions"></datalist>"""
        if include_disease_context
        else ""
    )
    disease_header = (
        '\n            <th><button type="button" class="sort-header" data-sort-key="disease">Top disease <span class="sort-arrow">↕</span></button></th>'
        if include_disease_context
        else ""
    )
    # Deferred so the disease-name filter haystack (large) loads after the
    # protein list has rendered from the smaller row payload.
    disease_index_script = (
        f'\n  <script defer src="portal-disease-index.js?v={_PORTAL_JS_VERSION}"></script>'
        if include_disease_context
        else ""
    )
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{escape(portal_title)}</title>
{_portal_favicon_links()}
  <link rel="stylesheet" href="portal.css?v={_PORTAL_CSS_VERSION}">
{_portal_analytics_head()}
</head>
<body>
  <main class="page">
    <header class="site-hero site-hero-index">
      {_portal_brand_header(active="index", portal_title=portal_title, sibling_link=sibling_link)}
      <div class="hero-grid">
        <div class="hero-copy">
          <p class="eyebrow">Recombinant antigen design</p>
          <h1>{escape(subtitle)}</h1>
          <p class="hero-text">{escape(intro)}</p>
        </div>
        <div class="hero-panel">
          <div class="stats">
            <div class="stat stat-primary"><span>{ready_total:,}</span><label>Completed reports</label></div>
            <div class="stat"><span>{topology_counts['secreted']:,}</span><label>Secreted</label></div>
            <div class="stat"><span>{topology_counts['single_pass']:,}</span><label>Single-pass</label></div>
            <div class="stat"><span>{topology_counts['multipass']:,}</span><label>Multipass</label></div>
            <div class="stat"><span>{topology_counts['gpi']:,}</span><label>GPI-anchored</label></div>
          </div>
        </div>
      </div>
    </header>

    {_agent_guidance_panel(portal_title=portal_title)}
    <section class="controls">
      <input id="searchBox" type="search" placeholder="Search gene, protein, or UniProt entry" aria-label="Search gene, entry, protein, alias, family, or topology">
      {disease_controls}
      <span id="sortSummary" class="sort-summary">Sorted by PubTator hits ↓</span>
      <button id="resetFilters" type="button">Reset</button>
    </section>

    <p class="agent-essentials">Confirm species and accession. Preserve residue numbering. Proposed constructs require experimental validation.
      <a href="agent-guide.html#interpret-supporting-evidence">How to read the evidence</a></p>
    <section class="pagination-bar" aria-label="Index pagination">
      <label class="page-size-control" for="pageSizeSelect">
        Rows per page
        <select id="pageSizeSelect">
          <option value="10">10</option>
          <option value="25" selected>25</option>
          <option value="50">50</option>
          <option value="100">100</option>
        </select>
      </label>
      <button id="pageFirst" type="button" class="mini-button">First</button>
      <button id="pagePrev" type="button" class="mini-button">Previous</button>
      <span id="pageSummary" class="pagination-summary"></span>
      <button id="pageNext" type="button" class="mini-button">Next</button>
      <button id="pageLast" type="button" class="mini-button">Last</button>
    </section>

    <section class="card table-card">
      <table id="resultsTable">
        <thead>
          <tr>
            <th><button type="button" class="sort-header" data-sort-key="index"># <span class="sort-arrow">↕</span></button></th>
            <th><button type="button" class="sort-header" data-sort-key="entry">Entry <span class="sort-arrow">↕</span></button></th>
            <th><button type="button" class="sort-header" data-sort-key="gene">Gene <span class="sort-arrow">↕</span></button></th>
            <th><button type="button" class="sort-header" data-sort-key="protein">Protein <span class="sort-arrow">↕</span></button></th>
            <th><button type="button" class="sort-header" data-sort-key="family">Family <span class="sort-arrow">↕</span></button></th>
            {disease_header}
            <th><button type="button" class="sort-header" data-sort-key="track">Track <span class="sort-arrow">↕</span></button></th>
            <th><button type="button" class="sort-header active" data-sort-key="pubtator">PubTator hits <span class="sort-arrow">↓</span></button></th>
            <th><button type="button" class="sort-header" data-sort-key="ecto">Design Region <span class="sort-arrow">↕</span></button></th>
            <th><button type="button" class="sort-header" data-sort-key="alphafold">AlphaFold <span class="sort-arrow">↕</span></button></th>
            <th><button type="button" class="sort-header" data-sort-key="pdb">PDB <span class="sort-arrow">↕</span></button></th>
            <th><button type="button" class="sort-header" data-sort-key="constructs">Constructs <span class="sort-arrow">↕</span></button></th>
            <th><button type="button" class="sort-header" data-sort-key="cross">BLAST hits <span class="sort-arrow">↕</span></button></th>
          </tr>
        </thead>
        <tbody></tbody>
      </table>
    </section>
    {_portal_footer()}
  </main>
  <script src="portal-index-data.js?v={_PORTAL_JS_VERSION}"></script>{disease_index_script}
  <script>
    const searchBox = document.getElementById('searchBox');
    const diseaseSearchBox = document.getElementById('diseaseSearchBox');
    const diseaseSuggestions = document.getElementById('diseaseSuggestions');
    const sortSummary = document.getElementById('sortSummary');
    const resetFilters = document.getElementById('resetFilters');
    const pageFirst = document.getElementById('pageFirst');
    const pagePrev = document.getElementById('pagePrev');
    const pageSummary = document.getElementById('pageSummary');
    const pageNext = document.getElementById('pageNext');
    const pageLast = document.getElementById('pageLast');
    const pageSizeSelect = document.getElementById('pageSizeSelect');
    const tbody = document.querySelector('#resultsTable tbody');
    // Keep rows as lightweight data objects and build a <tr> only for the rows
    // on the visible page. Materializing all rows up front (createElement +
    // innerHTML parse per row) was the dominant first-paint cost on large
    // snapshots — sorting, filtering and search only read `dataset`/`_haystack`,
    // never the DOM node, so the nodes are pure render output.
    function enrichPayload(payload) {{
      const d = payload.dataset || (payload.dataset = {{}});
      // Disease names ship separately (portal-disease-index.js) and are merged
      // in after first paint; default to empty so filters/sorts stay safe.
      if (d.diseases === undefined) d.diseases = '';
      payload._haystack = `${{d.query || ''}} ${{d.gene || ''}} ${{d.entry || ''}} ${{d.protein || ''}} ${{d.aliases || ''}} ${{d.family || ''}} ${{d.track || ''}}`;
      return payload;
    }}
    function buildRowElement(payload) {{
      const row = document.createElement('tr');
      const d = payload.dataset || {{}};
      for (const key in d) {{ row.dataset[key] = String(d[key] ?? ''); }}
      row.innerHTML = (payload.cells || []).map((cell) => `<td>${{cell || ''}}</td>`).join('');
      return row;
    }}
    const rows = (window.OpenAntigenIndexRows || []).map(enrichPayload);
    const sortHeaders = Array.from(document.querySelectorAll('.sort-header'));
    let pageSize = Number(pageSizeSelect.value) || 25;
    let sortState = {{ key: 'pubtator', dir: 'desc' }};
    let currentPage = 1;
    let filteredRows = rows.slice();

    function toNumber(value) {{
      const parsed = Number(value);
      return Number.isFinite(parsed) ? parsed : 0;
    }}

    function compareRows(a, b, key, dir) {{
      let comparison = 0;
      switch (key) {{
        case 'entry':
          comparison = a.dataset.entry.localeCompare(b.dataset.entry);
          break;
        case 'gene':
          comparison = a.dataset.gene.localeCompare(b.dataset.gene);
          break;
        case 'protein':
          comparison = a.dataset.protein.localeCompare(b.dataset.protein);
          break;
        case 'family':
          comparison = a.dataset.family.localeCompare(b.dataset.family);
          break;
        case 'disease':
          comparison = toNumber(a.dataset.diseaseScore) - toNumber(b.dataset.diseaseScore) || a.dataset.diseases.localeCompare(b.dataset.diseases);
          break;
        case 'track':
          comparison = a.dataset.track.localeCompare(b.dataset.track);
          break;
        case 'pubtator':
          comparison = toNumber(a.dataset.pubtator) - toNumber(b.dataset.pubtator);
          break;
        case 'ecto':
          comparison = toNumber(a.dataset.ectoStart) - toNumber(b.dataset.ectoStart) || toNumber(a.dataset.ectoEnd) - toNumber(b.dataset.ectoEnd);
          break;
        case 'alphafold':
          comparison = toNumber(a.dataset.hasAlphafold) - toNumber(b.dataset.hasAlphafold);
          break;
        case 'pdb':
          comparison = toNumber(a.dataset.hasPdb) - toNumber(b.dataset.hasPdb);
          break;
        case 'constructs':
          comparison = toNumber(a.dataset.constructs) - toNumber(b.dataset.constructs);
          break;
        case 'cross':
          comparison = toNumber(a.dataset.cross) - toNumber(b.dataset.cross);
          break;
        case 'index':
        default:
          comparison = toNumber(a.dataset.index) - toNumber(b.dataset.index);
      }}
      if (comparison === 0) {{
        comparison = a.dataset.gene.localeCompare(b.dataset.gene) || a.dataset.entry.localeCompare(b.dataset.entry);
      }}
      return dir === 'desc' ? -comparison : comparison;
    }}

    function updateSortHeaders() {{
      const labels = {{
        index: 'batch order',
        entry: 'entry',
        gene: 'gene',
        protein: 'protein',
        family: 'family',
        disease: 'top disease score',
        track: 'track',
        pubtator: 'PubTator hits',
        ecto: 'ectodomain',
        alphafold: 'AlphaFold availability',
        pdb: 'PDB availability',
        constructs: 'construct count',
        cross: 'BLAST hit count',
      }};
      sortHeaders.forEach((button) => {{
        const arrow = button.querySelector('.sort-arrow');
        if (button.dataset.sortKey === sortState.key) {{
          button.classList.add('active');
          arrow.textContent = sortState.dir === 'asc' ? '↑' : '↓';
        }} else {{
          button.classList.remove('active');
          arrow.textContent = '↕';
        }}
      }});
      const directionLabel = sortState.dir === 'asc' ? '↑' : '↓';
      sortSummary.textContent = `Sorted by ${{labels[sortState.key] || sortState.key}} ${{directionLabel}}`;
    }}

    function paginationPageCount() {{
      return Math.max(1, Math.ceil(filteredRows.length / pageSize));
    }}

    function renderCurrentPage() {{
      const pageCount = paginationPageCount();
      currentPage = Math.min(Math.max(currentPage, 1), pageCount);
      const start = (currentPage - 1) * pageSize;
      const end = Math.min(start + pageSize, filteredRows.length);
      tbody.replaceChildren(...filteredRows.slice(start, end).map(buildRowElement));
      if (filteredRows.length) {{
        pageSummary.textContent = `Page ${{currentPage}} of ${{pageCount}} · showing ${{start + 1}}-${{end}} of ${{filteredRows.length}} matching / ${{rows.length}} total`;
      }} else {{
        pageSummary.textContent = 'No rows match the current filters';
      }}
      pageFirst.disabled = currentPage <= 1;
      pagePrev.disabled = currentPage <= 1;
      pageNext.disabled = currentPage >= pageCount || filteredRows.length === 0;
      pageLast.disabled = currentPage >= pageCount || filteredRows.length === 0;
      updateSortHeaders();
    }}

    function applyFilters(resetPage = false) {{
      if (resetPage) currentPage = 1;
      const term = searchBox.value.toLowerCase().trim();
      const diseaseTerm = diseaseSearchBox ? diseaseSearchBox.value.toLowerCase().trim() : '';
      const matchingRows = [];
      rows.forEach((row) => {{
        const haystack = row._haystack;
        const termMatch = !term || haystack.includes(term);
        const diseaseMatch = !diseaseTerm || row.dataset.diseases.includes(diseaseTerm);
        const show = termMatch && diseaseMatch;
        if (show) matchingRows.push(row);
      }});
      filteredRows = matchingRows.sort((a, b) => compareRows(a, b, sortState.key, sortState.dir));
      renderCurrentPage();
    }}
    function debounce(fn, wait) {{
      let timer = null;
      return (...args) => {{
        if (timer) clearTimeout(timer);
        timer = setTimeout(() => fn(...args), wait);
      }};
    }}
    const debouncedApplyFilters = debounce(applyFilters, 150);
    searchBox.addEventListener('input', () => debouncedApplyFilters(true));
    if (diseaseSearchBox) {{
      diseaseSearchBox.addEventListener('input', () => {{
        updateDiseaseSuggestions();
        debouncedApplyFilters(true);
      }});
    }}
    pageFirst.addEventListener('click', () => {{
      currentPage = 1;
      renderCurrentPage();
    }});
    pagePrev.addEventListener('click', () => {{
      currentPage -= 1;
      renderCurrentPage();
    }});
    pageNext.addEventListener('click', () => {{
      currentPage += 1;
      renderCurrentPage();
    }});
    pageLast.addEventListener('click', () => {{
      currentPage = paginationPageCount();
      renderCurrentPage();
    }});
    pageSizeSelect.addEventListener('change', () => {{
      pageSize = Number(pageSizeSelect.value) || 25;
      currentPage = 1;
      renderCurrentPage();
    }});
    sortHeaders.forEach((button) => {{
      button.addEventListener('click', () => {{
        const key = button.dataset.sortKey;
        if (sortState.key === key) {{
          sortState.dir = sortState.dir === 'asc' ? 'desc' : 'asc';
        }} else {{
          sortState = {{ key, dir: key === 'gene' || key === 'entry' || key === 'protein' || key === 'family' || key === 'track' ? 'asc' : 'desc' }};
          if (key === 'index' || key === 'ecto') {{
            sortState.dir = 'asc';
          }}
        }}
        applyFilters(true);
      }});
    }});
    resetFilters.addEventListener('click', () => {{
      searchBox.value = '';
      if (diseaseSearchBox) diseaseSearchBox.value = '';
      sortState = {{ key: 'pubtator', dir: 'desc' }};
      applyFilters(true);
    }});
    let allDiseaseNames = Array.from(new Set(rows.flatMap((row) => (row.dataset.diseases || '').split('|').map((name) => name.trim()).filter(Boolean)))).sort((a, b) => a.localeCompare(b));
    function rebuildDiseaseSuggestions() {{
      allDiseaseNames = Array.from(new Set(rows.flatMap((row) => (row.dataset.diseases || '').split('|').map((name) => name.trim()).filter(Boolean)))).sort((a, b) => a.localeCompare(b));
    }}
    function updateDiseaseSuggestions() {{
      if (!diseaseSearchBox || !diseaseSuggestions) return;
      const term = diseaseSearchBox.value.toLowerCase().trim();
      const matches = allDiseaseNames.filter((name) => !term || name.toLowerCase().includes(term)).slice(0, 25);
      diseaseSuggestions.replaceChildren(...matches.map((name) => {{
        const option = document.createElement('option');
        option.value = name;
        return option;
      }}));
    }}
    // Merge the deferred disease-name index (portal-disease-index.js) once it
    // has loaded, then refresh suggestions and the active filter. Works over
    // http(s) and file:// alike since it is a plain script, not a fetch.
    let diseaseIndexApplied = false;
    function applyDiseaseIndex() {{
      if (diseaseIndexApplied) return;
      const index = window.OpenAntigenDiseaseIndex;
      if (!index) return;
      diseaseIndexApplied = true;
      rows.forEach((row) => {{
        const names = index[row.dataset.entry];
        if (names) row.dataset.diseases = names;
      }});
      rebuildDiseaseSuggestions();
      updateDiseaseSuggestions();
      applyFilters(false);
    }}
    updateDiseaseSuggestions();
    applyFilters();
    if (window.OpenAntigenDiseaseIndex) {{
      applyDiseaseIndex();
    }} else {{
      document.addEventListener('DOMContentLoaded', applyDiseaseIndex);
    }}
  </script>
</body>
</html>
"""


def render_help_page(*, portal_title: str = "OpenAntigens", include_disease_context: bool = True) -> str:
    is_mouse = portal_title == "OpenAntigens Mouse"
    target_scope = (
        "mouse orthologs resolved from OpenAntigens human targets"
        if is_mouse
        else "human cell-surface and secreted protein targets"
    )
    quick_search_scope = ""
    evidence_columns = (
        "topology, AlphaFold, PDB, construct-count, and BLAST-hit"
        if not include_disease_context
        else "topology, disease, AlphaFold, PDB, construct-count, and BLAST-hit"
    )
    search_description = (
        "The main search box scans gene symbol, UniProt entry, protein name, aliases, family name, and topology track."
        if not include_disease_context
        else "The main search box scans gene symbol, UniProt entry, protein name, aliases, family name, and topology track. Use the disease-focused filter to search indexed Open Targets disease associations."
    )
    sort_fields = (
        "entry, gene, protein name, family, topology track, PubTator hits, design region, AlphaFold availability, PDB availability, construct count, or BLAST-hit count"
        if not include_disease_context
        else "entry, gene, protein name, family, disease score, topology track, PubTator hits, design region, AlphaFold availability, PDB availability, construct count, or BLAST-hit count"
    )
    export_scope = (
        "mouse construct name, boundaries, sequence, and available cynomolgus monkey equivalent regions; the originating human target is retained as provenance"
        if is_mouse
        else "human construct name, boundaries, sequence, and equivalent mouse and cynomolgus monkey regions"
    )
    disease_row = (
        "<tr><td>Open Targets Disease Associations</td><td>Indirect and direct disease associations and scores from Open Targets.</td><td>Review the underlying evidence before drawing conclusions about a target's therapeutic value.</td></tr>"
        if include_disease_context
        else ""
    )
    homology_description = (
        "Mouse-to-cynomolgus monkey equivalent regions when reference mappings are available. The originating human target is listed separately as provenance."
        if is_mouse
        else "Human-to-mouse and human-to-cynomolgus monkey equivalent regions when reference mappings are available."
    )
    example_searches = (
        "<code>EGFR</code>, <code>ITGAV</code>, <code>DKK1</code>, <code>OPRM</code>, <code>cadherin</code>"
        if not include_disease_context
        else "<code>EGFR</code>, <code>ITGAV</code>, <code>DKK1</code>, <code>OPRM</code>, <code>cadherin</code>"
    )
    search_goal = "Search the gene" if not include_disease_context else "Search the gene, then use the disease-focused filter when needed"
    target_suffix = "MOUSE" if is_mouse else "HUMAN"
    primary_only_design = "mouse-only" if is_mouse else "human-only"
    cross_species_workflow = (
        "Check the mouse and available cynomolgus monkey sequences in the Live Selected Region, then review the originating human target separately as provenance."
        if is_mouse
        else "Check human, mouse, and cynomolgus monkey sequences in the Live Selected Region."
    )
    return _render_info_page(
        title="Help and Examples",
        active="help",
        eyebrow="First-time user guide",
        heading="How to use OpenAntigens reports",
        intro=(
            f"OpenAntigens assembles sequence, annotation, and structure evidence for {target_scope}. "
            "Search for a target, review its proposed construct boundaries, and copy the sequences you choose to test."
        ),
        portal_title=portal_title,
        body_html=f"""
    {_citation_section()}
    <section class="card doc-card">
      <h2>Quick start</h2>
      <ol class="doc-list">
        <li>Search for a gene, UniProt entry, protein name, alias, family, topology track{quick_search_scope} on the <a class="inline-link" href="index.html">Browse</a> page.</li>
        <li>Scan the {evidence_columns} columns to identify reports with the evidence needed for review.</li>
        <li>Open a target report by clicking the UniProt entry name. Public snapshots show the results of a completed build. Pending or failed rows retain their recorded status and error details.</li>
        <li>Use the <strong>Interactive Construct Builder</strong> to adjust boundaries and inspect linked sequences, warnings, and available structure panels. You can also export sequences from the precomputed construct cards.</li>
        <li>Review the ordered construct sections, copy TSV or FASTA outputs, and download the release index or manifest for downstream analysis.</li>
      </ol>
    </section>

    <section class="doc-grid">
      <article class="card doc-card">
        <h2>Search effectively</h2>
        <p>{search_description}</p>
        <p><strong>Example searches:</strong> {example_searches}.</p>
      </article>
      <article class="card doc-card">
        <h2>Filter and sort</h2>
        <p>The Browse page defaults to all rows sorted by PubTator hits from highest to lowest. Search terms narrow the table; <strong>Reset</strong> clears both search boxes and restores PubTator-hit sorting.</p>
        <p>Click table headers to sort by {sort_fields}.</p>
        <p>Choose how many rows to show, then use the page controls to browse the remaining results.</p>
      </article>
      <article class="card doc-card">
        <h2>Understand topology tracks</h2>
        <p><strong>Secreted</strong>, <strong>GPI-anchored</strong>, and <strong>single-pass</strong> proteins usually emphasize soluble antigen constructs from the mature secreted protein or extracellular region. <strong>Multipass</strong> proteins include a membrane-expression track, and mixed cases with large extracellular regions also include soluble extracellular-region constructs.</p>
        <p>For multipass targets, do not interpret the soluble construct list as a complete expression strategy. Review the membrane-expression suggestions and full-length context.</p>
      </article>
      <article class="card doc-card">
        <h2>Static release snapshots</h2>
        <p>Public releases contain the results recorded during their build. The pages can be browsed locally or on the website.</p>
        <p>Opening a page does not update its annotations.</p>
      </article>
    </section>

    <section class="card doc-card">
      <h2>Report section map</h2>
      <div class="table-scroll">
        <table>
          <thead><tr><th>Section</th><th>What it tells you</th><th>How to use it</th></tr></thead>
          <tbody>
            <tr><td>Target / Design Region</td><td>Resolved UniProt target, organism, sequence length, topology-derived design scope, links, PubTator count, and source identifiers.</td><td>Confirm target resolution and compare the extracellular, secreted, or membrane scope with the planned design.</td></tr>
            <tr><td>Interactive Construct Builder</td><td>Linked sequence, cysteine/furin warnings, copy-ready sequence exports, and structure/pLDDT/PAE panels when compatible local structure assets exist.</td><td>Adjust boundaries and export {export_scope}.</td></tr>
            {disease_row}
            <tr><td>Homology</td><td>{homology_description}</td><td>Review conservation, boundary transfer, and cross-species reagent risk before animal screening or validation.</td></tr>
            <tr><td>Sequence similarity and cross-reactivity context</td><td>Local BLAST hits ranked by bit score, with identity, coverage, E-value, coordinates, and alignment text.</td><td>Inspect high-scoring non-self hits as possible proteins to counter-screen.</td></tr>
            <tr><td>Cysteines / PTMs / Furin sites</td><td>Potential expression liabilities, curated UniProt PTM or processing annotations, and optional mutation guidance.</td><td>Review warnings before ordering constructs. Optional Cys-to-Ser and furin-site edits in the builder update names and sequences in real time.</td></tr>
            <tr><td>Family Context</td><td>Canonical family/paralog context and directional identity matrix when available.</td><td>Identify close paralogs that may create specificity risks or useful comparison antigens.</td></tr>
            <tr><td>Interactions / Assembly</td><td>Known interaction or complex context, including possible obligatory partners when evidence supports it.</td><td>Flag proteins that may require co-expression, partner chains, or extra caution when designing isolated soluble constructs.</td></tr>
            <tr><td>Construct Summary / Constructs</td><td>Precomputed construct classes, boundaries, homolog-equivalent sequences, images, and sequence exports.</td><td>Review construct classes in display order, then refine boundaries in the builder when needed.</td></tr>
          </tbody>
        </table>
      </div>
    </section>

    <section class="card doc-card">
      <h2>Recommended construct-review order</h2>
      <p>The summary and builder list soluble constructs in this order: <strong>full design region</strong>, <strong>PDB-backed</strong>, <strong>strict calculated</strong>, <strong>lenient calculated</strong>, then <strong>domain annotated</strong> constructs. For multipass proteins, review <strong>membrane-expression constructs</strong> separately from soluble extracellular-region constructs. The <a class="inline-link" href="constructs.html">Pre-generated Constructs page</a> explains how their boundaries are generated. Display order is not a validated ranking of expression performance.</p>
    </section>

    <section class="doc-grid">
      <article class="card doc-card">
        <h2>pLDDT in one sentence</h2>
        <p>pLDDT tells you whether AlphaFold is locally confident at each residue. High pLDDT supports local fold confidence, but it does not prove that two domains have a fixed relative orientation.</p>
      </article>
      <article class="card doc-card">
        <h2>PAE in one sentence</h2>
        <p>PAE reports AlphaFold's estimated uncertainty in relative placement between residues or regions. Low PAE within a block supports a cohesive structural unit; high PAE between blocks indicates uncertain relative placement.</p>
      </article>
    </section>

    <section class="card doc-card">
      <h2>Common workflows</h2>
      <div class="table-scroll">
        <table>
          <thead><tr><th>Goal</th><th>Recommended workflow</th></tr></thead>
          <tbody>
            <tr><td>Find antibody-discovery antigen candidates</td><td>{search_goal}, filter for structure availability, inspect full design-region and PDB-backed constructs, then review BLAST similarity and homolog-equivalent sequences.</td></tr>
            <tr><td>Design a compact domain antigen</td><td>Open the report, start with domain annotated and strict calculated constructs, inspect pLDDT/PAE boundaries, and export TSV/FASTA from the builder.</td></tr>
            <tr><td>Prioritize cross-species screening constructs</td><td>{cross_species_workflow} Prioritize mapped constructs with high identity and conserved boundaries.</td></tr>
            <tr><td>Review paralog cross-reactivity risk</td><td>Use Family Context and BLAST similarity together. The family matrix compares paralogs across their sequence regions; BLAST alignments locate similar stretches.</td></tr>
            <tr><td>Handle a multipass protein</td><td>Look for a large extracellular region if present, then separately review membrane-expression constructs, full-length AlphaFold context, full-length BLAST, and GPCR/membrane engineering suggestions.</td></tr>
            <tr><td>Review release-state gaps</td><td>Check the index status and error details. Public pages do not rerun an analysis when opened.</td></tr>
          </tbody>
        </table>
      </div>
    </section>

    <section class="doc-grid">
      <article class="card doc-card">
        <h2>Example: EGFR</h2>
        <p>Search <code>EGFR</code>, open <code>EGFR_{target_suffix}</code>, compare full extracellular-region and domain constructs, inspect PDB-backed boundaries, and check ERBB family context before choosing an antigen for antibody discovery.</p>
      </article>
      <article class="card doc-card">
        <h2>Example: ITGAV</h2>
        <p>Search <code>ITGAV</code>, review interaction and assembly context carefully, and treat isolated integrin-alpha constructs cautiously because integrins often require beta-chain context for native assembly.</p>
      </article>
      <article class="card doc-card">
        <h2>Example: DKK1</h2>
        <p>Search <code>DKK1</code> to see how secreted proteins are handled. Focus on mature extracellular sequence, cysteine-rich domains, disulfide context, and compact domain constructs.</p>
      </article>
      <article class="card doc-card">
        <h2>Example: OPRM</h2>
        <p>Search <code>OPRM</code> for a multipass example. Prioritize full-length membrane-expression context and treat short loops differently from soluble extracellular-region antigens.</p>
      </article>
    </section>

    <section class="card doc-card">
      <h2>Troubleshooting</h2>
      <div class="table-scroll">
        <table>
          <thead><tr><th>Problem</th><th>Likely reason</th><th>What to do</th></tr></thead>
          <tbody>
            <tr><td>Report is pending or failed</td><td>The target could not be completed for the current release because of source-data limitations, structure mismatch, or processing failure.</td><td>Inspect the release index error field. The target may be resolved in a future release.</td></tr>
            <tr><td>AlphaFold is missing</td><td>No canonical AlphaFold model was found or the available model did not match the canonical UniProt sequence.</td><td>Review available PDB and domain evidence. A new analysis or release is needed to incorporate a later structure.</td></tr>
            <tr><td>Homolog sequence is missing</td><td>A usable reference sequence or alignment was unavailable when the report was built.</td><td>Use {primary_only_design} until the ortholog table is updated, or inspect source identifiers manually.</td></tr>
            <tr><td>Family matrix is absent</td><td>No canonical family was identified, the family exceeded the configured size limit, or precomputed paralog data are unavailable.</td><td>Review BLAST local similarity and InterPro domain annotations for specificity context.</td></tr>
            <tr><td>BLAST hit list is long</td><td>The query region is highly conserved, contains common domains, or the full-length track captured broad family similarity.</td><td>Inspect the list, which is ranked by bit score, and check which alignments overlap the selected construct.</td></tr>
            <tr><td>Construct looks biologically wrong</td><td>Automated heuristics can miss ligand sites, partner requirements, topology edge cases, or literature-specific constraints.</td><td>Use the builder to adjust boundaries and check Methods for the exact evidence sources and limitations.</td></tr>
          </tbody>
        </table>
      </div>
    </section>

    <section class="doc-grid">
      <article class="card doc-card">
        <h2>Where to go next</h2>
        <p>Use the <a class="inline-link" href="builder.html">Builder guide</a> for detailed instructions on sequence/structure selection, pLDDT, PAE, Cys-to-Ser edits, and FASTA/TSV exports.</p>
      </article>
      <article class="card doc-card">
        <h2>Methods</h2>
        <p>Use the <a class="inline-link" href="methods.html">Methods page</a> for source databases, construct heuristics, BLAST setup, family/paralog logic, update cadence, and limitations.</p>
      </article>
      <article class="card doc-card">
        <h2>Construct methodology</h2>
        <p>Use the <a class="inline-link" href="constructs.html">Constructs page</a> for the dedicated explanation of pre-generated construct classes, structural diagnostics, homolog transfer, PTM and processing review, and limitations.</p>
      </article>
      <article class="card doc-card">
        <h2>Downloads</h2>
        <p>Use the <a class="inline-link" href="downloads.html">Downloads page</a> for the portal index, manifest, and flat files that support reproducible downstream analysis.</p>
      </article>
      <article class="card doc-card">
        <h2>Terms</h2>
        <p>Use the <a class="inline-link" href="terms.html">Terms page</a> for software license, generated annotation license, source-data terms, and contact information.</p>
      </article>
    </section>
""",
    )


def render_builder_page(*, portal_title: str = "OpenAntigens") -> str:
    is_mouse = portal_title == "OpenAntigens Mouse"
    selected_region_scope = (
        "mouse construct name, boundaries, sequence, and available cynomolgus monkey equivalent regions; originating human target identifiers are provenance metadata"
        if is_mouse
        else "human construct name, boundaries, sequence, and equivalent mouse and cynomolgus monkey regions"
    )
    equivalent_sequence_check = (
        "Does the available cynomolgus monkey equivalent sequence preserve the selected region if cross-species antibody screening is planned?"
        if is_mouse
        else "Do mouse and cynomolgus monkey equivalent sequences preserve the selected region if cross-species antibody screening is planned?"
    )
    return _render_info_page(
        title="Interactive Construct Builder Guide",
        active="builder",
        eyebrow="Construct design support",
        heading="How to use the Interactive Construct Builder",
        intro=(
            "Use the builder to inspect and adjust antigen boundaries from each protein report. "
            "Select a construct or enter start and end positions, inspect the available evidence, then copy its sequence."
        ),
        portal_title=portal_title,
        body_html=f"""
    <section class="doc-grid">
      <article class="card doc-card">
        <h2>What the builder shows</h2>
        <p>Reports with a compatible local AlphaFold model include the interactive builder and Live Selected Region table. The pLDDT trace and PAE matrix appear when their data are available. Reports without a compatible model retain precomputed construct cards and their sequence exports.</p>
        <p>Selections stay synchronized across the panels present in the report. Selecting residues in the sequence or plots highlights the same residue window in the structure and updates the TSV/FASTA outputs.</p>
      </article>
      <article class="card doc-card">
        <h2>Basic workflow</h2>
        <ol class="doc-list">
          <li>Start from the construct classes present in the report: full design region, PDB-backed constructs, strict calculated regions, lenient calculated regions, then annotated domains.</li>
          <li>Choose a construct, enter start and end positions and click Apply boundaries, or select residues in the sequence or available plots.</li>
          <li>Inspect the model and annotations at both boundaries.</li>
          <li>Review warnings for unpaired cysteines and furin-like cleavage sites.</li>
          <li>Copy the TSV or FASTA output for cloning, ordering, or downstream analysis.</li>
        </ol>
      </article>
      <article class="card doc-card">
        <h2>Live Selected Region</h2>
        <p>The table reports sequences for the current selection: {selected_region_scope} where ortholog mappings exist.</p>
        <p>The construct name follows the format <code>GENE_SPECIES_start-end</code>. Selecting optional Cys-to-Ser mutations adds mutation suffixes to the name and updates the sequence immediately.</p>
      </article>
      <article class="card doc-card">
        <h2>Structure panel</h2>
        <p>Reports with a local AlphaFold model show the structure as a cartoon colored by AlphaFold confidence. Selected residues are highlighted for boundary review against folded domains, linkers, cysteines, and visible surface regions.</p>
        <p>Cysteines are emphasized in yellow with side chains where present. Exposed or unpaired cysteines can complicate recombinant soluble antigen expression and purification.</p>
      </article>
    </section>

    <section class="card doc-card">
      <h2>How to interpret pLDDT</h2>
      <p><strong>pLDDT</strong> is AlphaFold's per-residue local confidence estimate. It describes confidence in the immediate geometry around a residue. Domain-domain orientation requires PAE review.</p>
      <div class="table-scroll">
        <table>
          <thead><tr><th>pLDDT range</th><th>Practical interpretation</th><th>Construct-design value</th></tr></thead>
          <tbody>
            <tr><td>90-100</td><td>Very high local confidence.</td><td>Check these residues against domain and topology annotations when choosing boundaries.</td></tr>
            <tr><td>70-90</td><td>Generally confident local structure.</td><td>Often acceptable for domain cores and boundary-adjacent residues.</td></tr>
            <tr><td>50-70</td><td>Low confidence or flexible geometry.</td><td>Review manually. This range often marks flexible linkers, uncertain loops, or boundary regions.</td></tr>
            <tr><td>Below 50</td><td>Very low confidence, often disordered.</td><td>Consider trimming unless the region has a required biological role.</td></tr>
          </tbody>
        </table>
      </div>
      <p>For soluble antigen design, high mean pLDDT and a high fraction of residues above the structured threshold support local fold confidence. Use PAE to evaluate relative orientation between domains.</p>
    </section>

    <section class="card doc-card">
      <h2>How to interpret PAE</h2>
      <p><strong>PAE</strong>, predicted aligned error, estimates how uncertain AlphaFold is about the position of one residue or region when another residue or region is aligned. Lower PAE means AlphaFold expects the relative placement of those residues to be reliable; higher PAE means the relative placement is uncertain.</p>
      <div class="table-scroll">
        <table>
          <thead><tr><th>PAE pattern</th><th>What it suggests</th><th>Construct-design value</th></tr></thead>
          <tbody>
            <tr><td>Low PAE block within a region</td><td>Low predicted error within the region.</td><td>Compare the block with annotated domains before choosing a construct.</td></tr>
            <tr><td>High PAE between two neighboring blocks</td><td>Uncertain relative placement of the blocks.</td><td>Consider splitting constructs at or near the linker.</td></tr>
            <tr><td>Low pLDDT and high boundary PAE</td><td>This is a candidate boundary between regions.</td><td>Review as a possible trimming or domain-separation point.</td></tr>
            <tr><td>High pLDDT but high inter-domain PAE</td><td>Each domain has local fold support, with uncertain domain orientation.</td><td>Expression and relative domain orientation require experimental testing.</td></tr>
          </tbody>
        </table>
      </div>
      <p>PAE can differ when the aligned and evaluated residues are reversed. It describes model uncertainty; it does not measure the stability of an expressed construct.</p>
    </section>

    <section class="doc-grid">
      <article class="card doc-card">
        <h2>Choosing boundaries</h2>
        <p>Inspect both termini and avoid cutting through a helix, strand, or required domain. A flexible linker may provide a suitable boundary. For soluble constructs, review whether signal peptides, transmembrane segments, cytoplasmic tails, or disordered extensions should be excluded, and account for annotated processing sites.</p>
        <p>Compare candidate boundaries with domain annotations, PDB coverage, and processing sites. Use pLDDT and PAE to inspect the model near each cut.</p>
      </article>
      <article class="card doc-card">
        <h2>Strict versus lenient constructs</h2>
        <p>Strict constructs start from pLDDT segments and may be split further using PAE. Lenient constructs use a lower pLDDT threshold and allow longer gaps. See Constructs for the thresholds and filters.</p>
      </article>
      <article class="card doc-card">
        <h2>Warnings</h2>
        <p>Cysteine warnings report missing or ambiguous partners in the model. PTM annotations identify recorded modifications and processing sites. Furin warnings mark sequence motifs for review; a motif match does not establish cleavage in your expression system.</p>
      </article>
      <article class="card doc-card">
        <h2>Manual overrides</h2>
        <p>Use manual boundary choices for ligand-binding sites, obligate partners, required disulfides, epitope constraints, glycosylation, and literature-supported constructs.</p>
      </article>
    </section>

    <section class="card doc-card">
      <h2>Practical checklist</h2>
      <ol class="doc-list">
        <li>Confirm the intended topology: soluble extracellular region or membrane protein.</li>
        <li>Inspect pLDDT across the selected sequence.</li>
        <li>Check PAE within the region and across its boundaries.</li>
        <li>Review terminal residues and any linker retained in the construct.</li>
        <li>Check cysteine partners, furin motifs, and annotated binding sites.</li>
        <li>Check the selected sequence against PTM and processing annotations.</li>
        <li>{equivalent_sequence_check}</li>
      </ol>
      <p>For the full methodology behind the pre-generated construct list, see the <a class="inline-link" href="constructs.html">Pre-generated Constructs page</a>.</p>
    </section>
""",
    )


def render_constructs_page(*, portal_title: str = "OpenAntigens") -> str:
    is_mouse = portal_title == "OpenAntigens Mouse"
    primary_species = "mouse" if is_mouse else "human"
    homolog_species = "Cynomolgus monkey" if is_mouse else "Mouse and cynomolgus monkey"
    identity_label = "Identity to mouse source" if is_mouse else "Identity to human"
    identity_description = (
        "Construct-specific sequence identity between homolog-equivalent and mouse regions."
        if is_mouse
        else "Construct-specific sequence identity between species-equivalent and human regions."
    )
    transfer_description = (
        "Mouse construct boundaries are transferred to available cynomolgus monkey homologs by aligning the full design scope and projecting the mouse boundary positions. The originating human target is retained as provenance rather than exported as a mapped construct."
        if is_mouse
        else "Human construct boundaries are transferred to mouse and cynomolgus monkey by aligning the full design scope to each ortholog sequence. The alignment projects human boundary positions onto ortholog sequences."
    )
    return _render_info_page(
        title="Pre-generated Constructs",
        active="constructs",
        eyebrow="Construct methodology",
        heading="Pre-generated antigen constructs and design heuristics",
        intro=(
            "OpenAntigens proposes boundaries using the target sequence and topology, domain annotations, AlphaFold confidence, PDB chain coverage, "
            "and ortholog mappings when those inputs are available."
        ),
        portal_title=portal_title,
        body_html=f"""
    <section class="card doc-card" id="choose-constructs">
      <h2>How to choose among pre-generated constructs</h2>
      <p>Start with the protein region your experiment requires, then compare the evidence supporting its boundaries. The categories describe how boundaries were chosen; they do not rank expected expression success.</p>
      <ol class="construct-choice-workflow" role="list" aria-label="Construct selection workflow">
        <li><strong>1. Define your region</strong><span>Decide whether you need the complete soluble region, a specific domain, or a membrane protein.</span></li>
        <li><strong>2. Check experimental precedent</strong><span>Look for a relevant PDB structure and review the construct described in the associated study.</span></li>
        <li><strong>3. Compare candidates</strong><span>If the best boundaries are unclear, compare a larger region with a smaller candidate.</span></li>
        <li><strong>4. Inspect and export</strong><span>Use the builder to inspect the boundaries and warnings before exporting your sequence.</span></li>
      </ol>
      <p>A relevant PDB structure can guide your first candidate. Check the original study for fusion proteins, affinity tags, mutations, and required partners that may be absent from the mapped target sequence.</p>
      <table class="construct-choice-table">
        <thead><tr><th scope="col">Your goal</th><th scope="col">Where to start</th></tr></thead>
        <tbody>
          <tr><td>Complete soluble region</td><td><strong>Full ectodomain / Full secreted region</strong></td></tr>
          <tr><td>Specific domain or repeat</td><td><strong>Annotated domains / repeats</strong>; compare <strong>Strict</strong> boundaries</td></tr>
          <tr><td>Neighboring domains together</td><td><strong>Lenient</strong> or the full soluble region</td></tr>
          <tr><td>Membrane protein</td><td><strong>Full-length multipass</strong>; compare <strong>Trimmed multipass</strong> and membrane <strong>PDB</strong> candidates</td></tr>
        </tbody>
      </table>
      <h3>Before export</h3>
      <ul class="doc-list">
        <li>Inspect both ends of the construct and retain any required domains or partners. Predicted boundaries may cut through helices or strands and do not establish that the isolated region will fold independently.</li>
        <li>Check whether truncation removed a cysteine's disulfide partner. Surface exposure alone does not justify mutating a cysteine, and burial does not establish that a substitution is safe.</li>
        <li>Review PTMs, processing sites, and relevant ortholog mappings. Use warnings to identify features that need closer inspection before deciding whether to keep or change the construct.</li>
      </ul>
      <p class="construct-choice-more">More detail: <a class="inline-link" href="#construct-classes">construct categories</a> · <a class="inline-link" href="builder.html">using the builder</a> · <a class="inline-link" href="#construct-feature-review">cysteines and processing</a>.</p>
    </section>

    <section class="card doc-card" id="after-export">
      <h2>After export</h2>
      <p>The exported sequence contains your selected target region and chosen edits. Choose a vector and host appropriate for your experiment, adding a secretion leader or tags where needed. Check whether native signal or processing sequences are already included.</p>
      <p>Practical resources: <a class="inline-link" href="https://blog.addgene.org/plasmids-101-protein-tags">Addgene&rsquo;s protein-tag guide</a> · <a class="inline-link" href="https://doi.org/10.1016/j.nbt.2020.05.002">Tegel et al.: mammalian protein production</a></p>
    </section>

    <section class="card doc-card">
      <h2>What the construct list contains</h2>
      <p>Each report lists proposed constructs with their sequences and the evidence used to choose the boundaries. The available classes depend on the annotations and structures found for that target.</p>
      <p>Before ordering a sequence, check the target literature, required partners, and the assay you plan to run.</p>
    </section>

    <section class="card doc-card">
      <h2>Strict and lenient construct algorithm</h2>
      <p>Calculated constructs use a compatible AlphaFold model and PAE matrix. Strict constructs use conservative pLDDT seeds followed by recursive PAE boundary scoring; lenient constructs use a lower pLDDT threshold to preserve larger structured regions.</p>
      <p>Both tracks begin with intervals selected by pLDDT. The lenient track keeps these intervals; the strict track can split them further using PAE.</p>
      <p>The steps below use the default configuration. A release built with different settings can produce different boundaries.</p>
    </section>

    <section class="card doc-card">
      <h2>Step-by-step strict and lenient generation</h2>
      <ol class="doc-list">
        <li><strong>Define the eligible design region.</strong> The pipeline first determines the region that can be used for soluble construct design. For secreted proteins this is usually the mature secreted chain; for single-pass proteins it is usually the extracellular region; for multipass proteins, large extracellular regions can be treated as mixed soluble/membrane cases, while membrane-expression constructs are handled separately. Strict and lenient soluble constructs stay inside this region.</li>
        <li><strong>Retrieve canonical AlphaFold confidence data when available.</strong> The canonical UniProt sequence is matched to the AlphaFold structure and PAE matrix. pLDDT is used as a per-residue local-confidence signal. PAE is used as a pairwise confidence signal that indicates whether two parts of the model have a confident relative orientation.</li>
        <li><strong>Generate lenient pLDDT regions.</strong> The lenient track marks residues with pLDDT at least 60, allows low-confidence gaps up to 12 residues, and discards resulting regions shorter than 25 residues. This track preserves larger structured units and can include more than one domain when the intervening linker is short or only moderately uncertain.</li>
        <li><strong>Generate strict pLDDT seed regions.</strong> The strict track marks residues with pLDDT at least 70, allows low-confidence gaps up to 8 residues, and discards seed regions shorter than 25 residues. These high-confidence intervals then feed PAE-based domain splitting.</li>
        <li><strong>Evaluate PAE split candidates inside each strict seed.</strong> For every possible cut, both resulting sides must remain at least 40 amino acids. The pipeline computes intrablock PAE within each side, interblock PAE between the two sides, and separation, defined as interblock PAE minus mean intrablock PAE.</li>
        <li><strong>Apply the PAE thresholds.</strong> A cut is eligible only when interblock PAE is at least 12 angstroms and separation is at least 4 angstroms. These thresholds select cuts where the model is less certain about placement between the blocks than within them.</li>
        <li><strong>Score eligible strict split boundaries.</strong> Candidate boundaries are ranked using <code>split score = separation + 0.35 x max(0, boundary PAE - intrablock PAE) + linker bonus</code>. The linker bonus increases when local pLDDT near the cut is low, because low-confidence linker residues support splitting adjacent structured blocks.</li>
        <li><strong>Recursively split strict regions.</strong> The highest-scoring eligible cut is applied, then the same PAE split search is repeated on the left and right child regions. Recursion stops when no valid cut remains, child regions would become too small, or the configured maximum recursion depth is reached.</li>
        <li><strong>Apply final construct filters.</strong> Calculated soluble constructs must stay inside the design region and must be at least 50 amino acids. Exact duplicate boundaries from different calculated or annotated sources are collapsed to keep reports concise; PDB-backed constructs are exempt from deduplication but still pass the scope and length filters.</li>
        <li><strong>Annotate each construct for review.</strong> Construct cards list sequence, boundaries, exports, included UniProt/InterPro annotations, PTMs, furin-like motifs, cysteine warnings, homolog-equivalent regions, and construct-specific conservation. Structural metrics, PAE split diagnostics, ligand-interaction annotations, and images appear when available.</li>
      </ol>
    </section>

    <section class="card doc-card" id="construct-classes">
      <h2>Construct classes</h2>
      <p>Constructs are grouped by how their boundaries were chosen. The table follows the Construct Details tab order.</p>
      <div class="table-scroll">
        <table>
          <thead><tr><th>Display order</th><th>Construct class</th><th>How boundaries are generated</th><th>Review focus</th></tr></thead>
          <tbody>
            <tr><td>1</td><td>Full ectodomain / Full secreted region (Full design region when neither label applies)</td><td>Uses the inferred mature secreted or extracellular design scope from topology, signal peptide, propeptide, chain, and transmembrane annotations.</td><td>Comparing the complete inferred soluble region with shorter constructs.</td></tr>
            <tr><td>2</td><td>PDB</td><td>Uses PDB chain coverage mapped to UniProt numbering. Eligible PDB chains are kept individually even when their boundaries duplicate another construct. Soluble and membrane candidates appear separately within the PDB tab. The construct name links to its RCSB entry.</td><td>Reviewing deposited structure boundaries.</td></tr>
            <tr><td>3</td><td>Annotated domains / repeats</td><td>Uses UniProt and InterPro domain-like annotations that lie fully inside the design scope.</td><td>Designing named biological domains or domain modules.</td></tr>
            <tr><td>4</td><td>Strict</td><td>Uses more conservative AlphaFold pLDDT/PAE segmentation to isolate compact, cohesive structural regions.</td><td>Smaller regions for experimental testing.</td></tr>
            <tr><td>5</td><td>Lenient</td><td>Uses more permissive structural thresholds to retain larger coherent or partially coherent regions.</td><td>Multi-domain surfaces and larger conformational epitopes.</td></tr>
            <tr><td>6</td><td>Full-length multipass</td><td>Retains the full target sequence, including transmembrane regions.</td><td>Reviewing the complete membrane-protein context.</td></tr>
            <tr><td>7</td><td>Trimmed multipass</td><td>Retains transmembrane regions with N-terminal, C-terminal, both-terminal, or GPCR C-terminal tail-series trims.</td><td>Comparing terminal variants with the full-length candidate; review each variant's warnings.</td></tr>
          </tbody>
        </table>
      </div>
      <p>The default minimum construct length is 50 amino acids, including PDB-derived suggestions. Experimental PDB records also appear separately in the report.</p>
    </section>

    <section class="doc-grid">
      <article class="card doc-card">
        <h2>Design scope comes first</h2>
        <p>OpenAntigens first determines the design scope: mature secreted sequence, GPI/single-pass extracellular region, a large multipass extracellular region, or a full-length membrane-protein expression region. Soluble constructs stay inside that inferred scope.</p>
        <p>For secreted proteins, signal peptide and propeptide annotations help define mature sequence. For single-pass proteins, transmembrane and topology annotations define extracellular boundaries. For multipass proteins, large extracellular regions are treated as mixed cases; short loops usually route to membrane-expression target review.</p>
      </article>
      <article class="card doc-card">
        <h2>Boundary filtering</h2>
        <p>Recommended soluble constructs must fit inside the design scope. InterPro or UniProt annotations that extend outside the mature secreted or extracellular region remain descriptive context.</p>
        <p>Soluble PDB-derived suggestions must lie fully inside the design region and pass the minimum-length filter. Multipass targets can also retain PDB chains as membrane-expression suggestions. The experimental-structure section provides the broader record of overlapping PDB evidence.</p>
      </article>
    </section>

    <section class="card doc-card">
      <h2>AlphaFold-derived structural heuristics</h2>
      <p>Construct cards report the following metrics when the required structure data are available.</p>
      <div class="table-scroll">
        <table>
          <thead><tr><th>Diagnostic</th><th>What it measures</th><th>How to interpret it</th></tr></thead>
          <tbody>
            <tr><td>Mean pLDDT</td><td>Average local AlphaFold confidence across the construct.</td><td>Higher values support local fold confidence. Low values suggest disorder, flexible tails, or uncertain loops.</td></tr>
            <tr><td>Structured coverage</td><td>Fraction of construct residues above the configured pLDDT structured threshold.</td><td>Higher coverage supports a more consistently structured construct.</td></tr>
            <tr><td>Mean intra-domain PAE</td><td>Average PAE among residues inside the construct or block.</td><td>Lower values suggest a coherent structural unit.</td></tr>
            <tr><td>Interblock PAE</td><td>PAE between two candidate neighboring blocks around a possible split.</td><td>Higher interblock PAE can support separating blocks into distinct constructs.</td></tr>
            <tr><td>Separation</td><td>Difference between interblock and intrablock PAE around a candidate boundary.</td><td>Larger separation supports a domain boundary or flexible hinge.</td></tr>
            <tr><td>Local boundary PAE</td><td>PAE near a candidate split boundary.</td><td>High boundary PAE can indicate uncertain relative orientation near the split.</td></tr>
            <tr><td>Linker pLDDT</td><td>pLDDT around residues connecting blocks.</td><td>Low linker pLDDT supports trimming or splitting near that linker.</td></tr>
          </tbody>
        </table>
      </div>
      <p>Interpret these diagnostics with domain annotations, PDB evidence, PTMs, ligand sites, cysteines, and the intended antibody-discovery strategy.</p>
    </section>

    <section class="card doc-card">
      <h2>Domain and PDB evidence</h2>
      <div class="table-scroll">
        <table>
          <thead><tr><th>Evidence</th><th>How it is used</th><th>Important limitation</th></tr></thead>
          <tbody>
            <tr><td>UniProt domains/features</td><td>Curated features are merged into the residue map and can define domain annotated constructs when boundaries fit inside the design scope.</td><td>Feature granularity varies by protein and curator coverage.</td></tr>
            <tr><td>InterPro domains</td><td>Specific domain signatures can define domain annotated constructs and help classify construct content.</td><td>Broad whole-protein signatures are filtered from construct recommendations when they extend outside the relevant design scope.</td></tr>
            <tr><td>Pfam domains</td><td>Useful as domain evidence and descriptive annotations.</td><td>Paralog matrices use canonical family identifiers; Pfam-only family IDs remain domain evidence.</td></tr>
            <tr><td>PDB mapped regions</td><td>Mapped experimental boundaries can become PDB-backed constructs and provide precedent for construct boundaries.</td><td>PDB constructs document deposited structure boundaries; antibody-discovery context still determines final suitability.</td></tr>
          </tbody>
        </table>
      </div>
    </section>

    <section class="card doc-card">
      <h2>What each construct card includes</h2>
      <p>Each card contains the primary {primary_species} sequence and boundaries, with TSV/FASTA exports. Mapped homologs, annotations, and structural metrics appear when available.</p>
      <div class="table-scroll">
        <table>
          <thead><tr><th>Field</th><th>Meaning</th><th>Why it matters</th></tr></thead>
          <tbody>
            <tr><td>Name</td><td><code>GENE_SPECIES_start-end</code>, with optional mutation suffixes when edits are applied in the builder.</td><td>Identifies the sequence and any selected edits.</td></tr>
            <tr><td>Boundary</td><td>Amino-acid start/end in the source protein sequence.</td><td>Defines the exact construct region.</td></tr>
            <tr><td>Sequence</td><td>Construct amino-acid sequence.</td><td>Can be copied directly into cloning or ordering workflows.</td></tr>
            <tr><td>Homologs</td><td>{homolog_species} equivalent regions when mapping is available.</td><td>Supports cross-species reagent and screening strategy.</td></tr>
            <tr><td>{identity_label}</td><td>{identity_description}</td><td>Helps triage conservation and cross-species reagent risk; antibody behavior still depends on epitope, structure, glycans, and assay context.</td></tr>
            <tr><td>PTMs included</td><td>Curated UniProt PTM or processing annotations overlapping the construct.</td><td>Flags glycosylation, cleavage, chain, propeptide, lipidation, or other features that may affect design.</td></tr>
            <tr><td>Cysteines</td><td>Construct-specific paired/unpaired cysteine analysis.</td><td>Flags potential disulfide, aggregation, or expression liabilities.</td></tr>
            <tr><td>Furin motifs</td><td>Basic furin-like motifs inside the construct.</td><td>Identifies motifs to check before choosing an expression system.</td></tr>
            <tr><td>Images/plots</td><td>Static structure PNGs when generated, plus pLDDT/PAE quality plots when available.</td><td>Provides a quick check of model confidence and boundary placement.</td></tr>
          </tbody>
        </table>
      </div>
    </section>

    <section class="doc-grid" id="construct-feature-review">
      <article class="card doc-card">
        <h2>Homolog transfer</h2>
        <p>{transfer_description} Residue numbers can differ across species.</p>
        <p>Missing mappings are marked unavailable. An available mapping still needs review, especially near gaps or repeated domains.</p>
      </article>
      <article class="card doc-card">
        <h2>Cysteine edits</h2>
        <p>Pre-generated construct cards report unpaired cysteine warnings. In the Interactive Construct Builder, checkboxes can apply Cys-to-Ser edits to selected unpaired cysteines; names, TSV, and FASTA update immediately.</p>
        <p>The builder excludes ambiguous cysteine contacts from its suggested edits. Check the model and curated disulfide annotations before applying a substitution.</p>
      </article>
      <article class="card doc-card">
        <h2>PTM and processing review</h2>
        <p>Construct cards list overlapping PTMs and processing features for boundary review.</p>
        <p>Cleavage, propeptide, signal peptide, and chain annotations deserve special attention because they can define mature protein boundaries.</p>
      </article>
      <article class="card doc-card">
        <h2>Multipass proteins</h2>
        <p>Large extracellular regions on multipass proteins are handled as mixed cases with soluble extracellular-region constructs and full-length membrane-protein context. Proteins with only short extracellular loops are treated primarily as membrane-expression targets.</p>
        <p>Review membrane-expression suggestions separately; they retain transmembrane sequence.</p>
      </article>
    </section>

    <section class="card doc-card">
      <h2>Limitations</h2>
      <p>The generator proposes sequence boundaries; it does not predict expression yield, folding, or antibody binding. Review the target literature and test the selected constructs.</p>
    </section>
""",
    )


def render_methods_page(*, portal_title: str = "OpenAntigens") -> str:
    is_mouse = portal_title == "OpenAntigens Mouse"
    target_identity_source = "Reviewed UniProt mouse records" if is_mouse else "Reviewed UniProt human records"
    target_universe_text = (
        "OpenAntigens Mouse contains unique mouse orthologs resolved from the OpenAntigens human target universe. "
        "Targets without a resolved mouse record are omitted, and duplicate mouse records are included once. "
        "Each included target is analyzed against reviewed mouse UniProt data; the originating human target remains in report metadata."
        if is_mouse
        else "The human target universe is sourced from the refined cell-surface and secreted protein topology table under <code>data/inputs/</code>. "
        "That table combines public evidence from UniProt, Thera-SAbDab, the SURFY surfaceome publication, and the Uhlén secretome publication. "
        "The active portal includes proteins categorized as secreted, GPI-anchored, single-pass membrane, or multipass membrane."
    )
    target_resolution_text = (
        "Mouse records are resolved through reviewed mouse UniProt entries. The build records skipped human targets, including missing mouse orthologs and duplicate mouse targets, in a separate manifest."
        if is_mouse
        else "Each input row is resolved primarily through reviewed human UniProt records. Resolution prefers stable UniProt entry names or accessions when available and preserves the original input identifiers for traceability. Ambiguous or unresolved rows remain represented in batch outputs with explicit error status."
    )
    ortholog_text = (
        "The mouse portal keeps the originating human target in report metadata and maps additional homolog regions when suitable reference sequences are available."
        if is_mouse
        else "Precomputed ortholog reference tables provide human, mouse, and cynomolgus monkey accessions, sequences, and pairwise identity to the human sequence. Mouse resolution ranks matches across the HCOP and human fallback symbols together, always preferring reviewed UniProt over unreviewed TrEMBL. Reports use these tables first. If no precomputed row covers the target, the analyzer runs live species lookup."
    )
    homolog_projection_text = (
        "Construct-level homolog sequences are inferred by aligning the mouse design scope to available homolog sequences, projecting the selected mouse boundaries through that alignment, and extracting the corresponding homolog region. A projection is rejected when fewer than 70% of the target-region residues align to the homolog."
        if is_mouse
        else "Construct-level homolog sequences are inferred by aligning the full human design scope to species ortholog sequences, projecting the human construct boundaries through that alignment, and extracting the corresponding species region. A projection is rejected when fewer than 70% of the human-region residues align to the homolog. Construct-specific identity is then calculated against the human construct sequence."
    )
    blast_database_text = (
        "Mouse target sequences are searched against the configured local human, mouse, and cynomolgus monkey protein databases."
        if is_mouse
        else "Human, mouse, and cynomolgus monkey primary databases use reviewed UniProt protein sets. For cynomolgus monkey, targets without hits in the reviewed set are searched against a broader local UniProt proteome database when that database is configured."
    )
    ortholog_source = (
        "Reviewed UniProt homolog lookup"
        if is_mouse
        else "Precomputed RefSeq/HCOP-derived ortholog table"
    )
    ortholog_use = (
        "Mouse and available cynomolgus monkey sequence mapping and construct-level identity; the originating human target is retained as provenance."
        if is_mouse
        else "Human, mouse, and cynomolgus monkey sequence mapping and construct-level identity."
    )
    disease_row = (
        '<tr><td>Disease context</td><td>Open Targets Platform</td><td>Target-disease association scores and disease-aware portal search.</td></tr>'
        if not is_mouse
        else ""
    )
    disease_methods = (
        "Open Targets direct target-disease associations are precomputed as a separate resumable artifact. The portal shows the top disease association in the index, loads a disease-name index for disease-aware filtering, and lists disease score plus datasource-score context on each target page."
        if not is_mouse
        else "The mouse portal does not display Open Targets disease associations or disease-search controls."
    )
    return _render_info_page(
        title="Methods",
        active="methods",
        eyebrow="Data generation",
        heading="How OpenAntigens generates portal annotations",
        intro=(
            "OpenAntigens combines UniProt annotations, AlphaFold models, PDB mappings, "
            "ortholog and paralog references, and local BLAST searches to annotate proposed constructs."
        ),
        portal_title=portal_title,
        body_html=f"""
    <section class="card doc-card">
      <h2>Overview</h2>
      <p>Database paper: {citation_short_html()} (preprint). See <a class="inline-link" href="help.html#cite-openantigens">How to cite OpenAntigens</a>.</p>
      <p>Each release records the results of a local analysis. Reports link proposed construct boundaries to the sequence, annotations, and structural evidence used in that analysis.</p>
      <p>The methods below describe the current pipeline. Previously generated releases retain the results of the software and reference data used for their build.</p>
      <div class="table-scroll">
        <table>
          <thead><tr><th>Evidence class</th><th>Primary source</th><th>Used for</th></tr></thead>
          <tbody>
            <tr><td>Target identity and sequence</td><td>{target_identity_source}</td><td>Canonical protein sequence, gene symbol, entry name, protein name, feature annotations.</td></tr>
            <tr><td>Target-universe support</td><td>UniProt, Thera-SAbDab, SURFY, and the Uhlén secretome publication</td><td>Public evidence retained in the canonical human input table used to define the source target universe.</td></tr>
            <tr><td>Topology and extracellular scope</td><td>Input topology context plus UniProt features</td><td>Secreted, GPI, single-pass, and multipass track assignment; soluble design-region or full-length membrane design scope.</td></tr>
            <tr><td>Predicted structure confidence</td><td>AlphaFold DB</td><td>pLDDT, PAE, structure visualization, structural segmentation, construct diagnostics.</td></tr>
            <tr><td>Experimental structure precedent</td><td>RCSB PDB mappings</td><td>Deposited construct boundaries and experimentally observed extracellular regions.</td></tr>
            <tr><td>GPCR-specific annotations</td><td>GPCRdb</td><td>GPCR class/family, segment boundaries, generic residue numbering, conserved motif context, and GPCR-focused construct interpretation.</td></tr>
            <tr><td>Domains and families</td><td>InterPro, UniProt, HGNC, Ensembl-derived references</td><td>Domain annotation, canonical family selection, paralog identity matrices.</td></tr>
            <tr><td>Orthologs</td><td>{ortholog_source}</td><td>{ortholog_use}</td></tr>
            <tr><td>Sequence similarity</td><td>Local blastp databases</td><td>Ranked non-self sequence similarity hits and original alignment text for cross-reactivity review.</td></tr>
            {disease_row}
            <tr><td>Literature context</td><td>PubTator3</td><td>Approximate target-linked PubMed hit counts and literature search links.</td></tr>
            <tr><td>PTMs and processing</td><td>UniProt feature annotations</td><td>Glycosylation, modified residues, lipidation, disulfides, cross-links, signal/propeptide/chain annotations, and cleavage-related sites.</td></tr>
          </tbody>
        </table>
      </div>
    </section>

    <section class="doc-grid">
      <article class="card doc-card">
        <h2>Target universe and resolution</h2>
        <p>{target_universe_text}</p>
        <p>{target_resolution_text}</p>
      </article>
      <article class="card doc-card">
        <h2>Topology and construct scope</h2>
        <p>Topology controls which design track is used. Secreted proteins use the mature extracellular protein as the soluble design scope. GPI-anchored proteins without a conventional transmembrane helix use the envelope of processed UniProt chain annotations. If no processed chain contains the annotated GPI anchor, OpenAntigens retains the post-signal sequence without trimming the GPI propeptide and marks that limitation. Single-pass proteins use the extracellular region defined by topology and transmembrane annotations.</p>
        <p>Multipass proteins are treated separately. When an eligible extracellular region has at least 80 amino acids by default, OpenAntigens includes a mixed track with both soluble extracellular-region constructs and full-length membrane-expression context. With shorter extracellular loops, the report emphasizes membrane-protein expression and full-length context.</p>
      </article>
      <article class="card doc-card">
        <h2>Sequence canonicality</h2>
        <p>OpenAntigens requires the report sequence, UniProt sequence, and AlphaFold model sequence to represent the same canonical isoform. Mismatching AlphaFold models are excluded from structure-dependent metrics because residue-number mismatches corrupt construct boundaries, pLDDT coloring, and sequence-to-structure selection.</p>
        <p>AlphaFold DB models are used only when they match the canonical report sequence.</p>
        <p>When a target lacks a matching AlphaFold DB model, non-structure-dependent annotations remain available; structure-derived segmentation and interactive structure features are limited.</p>
      </article>
      <article class="card doc-card">
        <h2>AlphaFold confidence</h2>
        <p>pLDDT is used as local residue-confidence evidence. High pLDDT supports local fold confidence, while low pLDDT is consistent with flexible linkers, disordered tails, or uncertain loops. PAE is used as regional confidence evidence: low internal PAE supports a coherent structural unit, while high PAE between neighboring blocks indicates uncertain relative orientation.</p>
        <p>Construct diagnostics summarize mean pLDDT, structured fraction, mean intrablock PAE, interblock PAE, local boundary PAE, separation, and linker pLDDT when available. These metrics are evidence features, not absolute pass/fail criteria.</p>
      </article>
      <article class="card doc-card">
        <h2>Residue topology and accessibility</h2>
        <p>Residues are classified by topology first, using UniProt topology features, transmembrane segments, signal peptides, propeptides, mature-chain annotations, and the input topology track. The reported location categories are extracellular, membrane, intracellular, or unknown. For multipass proteins, extracellular termini and loops are treated as extracellular regions inside a membrane-protein model.</p>
        <p>Surface accessibility is computed from the canonical AlphaFold structure with FreeSASA using the Lee-Richards algorithm. The default probe radius is 1.4 angstroms, the default slice count is 20, and residues are called solvent-accessible when relative SASA is at least 0.20.</p>
      </article>
    </section>

    <section class="card doc-card">
      <h2>Extracellular accessible residues</h2>
      <p>OpenAntigens combines topology and model accessibility to select extracellular residues for comparison. This set is used in ortholog and family identity calculations, BLAST summaries, and alignment highlighting.</p>
      <div class="table-scroll">
        <table>
          <thead><tr><th>Accessibility class</th><th>How it is assigned</th><th>Design interpretation</th></tr></thead>
          <tbody>
            <tr><td>Exposed</td><td>Relative SASA in the full AlphaFold model is at least the configured threshold.</td><td>Residue is directly solvent-accessible in the full model.</td></tr>
            <tr><td>Conditional</td><td>The residue is buried in the full model but exposed when SASA is recomputed inside a PAE-coherent block, or it has low pLDDT and is treated as accessible but low-confidence/disordered.</td><td>Flag for flexible or uncertain regions where full-model AlphaFold packing may hide residues exposed in solution or on cells.</td></tr>
            <tr><td>Buried</td><td>Relative SASA is below threshold in the full model and in available PAE-block contexts.</td><td>Lower antibody accessibility unless the biological conformation differs from the model.</td></tr>
            <tr><td>Uncertain</td><td>No structure or SASA context is available.</td><td>Manual accessibility review required.</td></tr>
          </tbody>
        </table>
      </div>
      <p>The pipeline also calculates SASA within PAE blocks to inspect residues that may be hidden by uncertain packing between blocks. Recomputing SASA within a block shows how its calculated accessibility changes when other blocks are removed; it does not establish accessibility in solution or on cells.</p>
      <p>SASA calculations use the target-chain model and selected structural blocks, without partner chains. The pipeline does not reconstruct biological homo-oligomeric assemblies or incorporate their interfaces into boundary selection. A residue exposed in the model may be buried against another subunit, including an identical copy, in the biological assembly.</p>
      <p>The extracellular accessible identity reported in homolog, family, and BLAST sections is calculated only over aligned primary-target residues that are both extracellular and classified as exposed or conditional.</p>
    </section>

    <section class="card doc-card">
      <h2>Construct generation</h2>
      <p>The analysis engine scores construct candidates, removes exact duplicate boundaries across calculated and annotated candidate sources (PDB-backed constructs are kept individually as evidence), filters soluble candidates outside the inferred design scope or below the minimum length, and the portal displays the surviving classes in a fixed review order.</p>
      <div class="table-scroll">
        <table>
          <thead><tr><th>Construct class</th><th>How it is generated</th><th>Primary use</th></tr></thead>
          <tbody>
            <tr><td>Full design region</td><td>Uses the inferred mature secreted or extracellular design region after signal peptide/propeptide/topology processing.</td><td>Provides the full inferred soluble region for comparison with shorter constructs.</td></tr>
            <tr><td>PDB-backed</td><td>Uses PDB chain coverage mapped to UniProt numbering. Soluble suggestions must fit inside the design region; membrane suggestions are handled separately. Eligible chains retain their PDB links even when boundaries repeat.</td><td>Captures experimentally observed boundaries.</td></tr>
            <tr><td>Domain annotated</td><td>Uses curated UniProt/InterPro domain boundaries that lie fully inside the design scope.</td><td>Captures biologically named domain units.</td></tr>
            <tr><td>Strict calculated</td><td>Uses more conservative pLDDT/PAE segmentation to favor compact single-domain-like units.</td><td>Emphasizes smaller constructs with stronger cohesive-domain support.</td></tr>
            <tr><td>Lenient calculated</td><td>Uses more permissive segmentation to retain larger structured units and multi-domain regions when support is acceptable.</td><td>Captures bigger surfaces and possible conformational epitopes.</td></tr>
            <tr><td>Membrane-expression</td><td>For multipass targets, trims obviously disordered cytoplasmic regions while retaining membrane-spanning context.</td><td>Supports whole membrane-protein expression planning.</td></tr>
          </tbody>
        </table>
      </div>
      <p>The minimum construct size defaults to 50 amino acids. Soluble constructs outside the inferred design scope are excluded. Exact duplicate boundaries from different calculated or annotated sources are collapsed during construct deduplication; PDB-backed constructs are exempt from deduplication, but still pass the scope and length filters.</p>
      <p>The <a class="inline-link" href="constructs.html">Constructs page</a> gives the default pLDDT thresholds, PAE split criteria, and boundary score. These classes describe the algorithm used, not measured expression or binding.</p>
    </section>

    <section class="doc-grid">
      <article class="card doc-card">
        <h2>Domain annotations</h2>
        <p>UniProt feature annotations and InterPro records are merged into the residue annotation map. Domain annotations must lie inside the relevant design scope to be used as construct recommendations. Large whole-protein signatures that extend beyond the mature secreted or extracellular region are retained as descriptive context when appropriate. They are excluded from soluble construct recommendations.</p>
        <p>InterPro provides integrated domain and family annotations. Pfam entries support domain annotation. Pfam-only family IDs are excluded as canonical family identifiers for precomputed family matrices.</p>
      </article>
      <article class="card doc-card">
        <h2>Family and paralog context</h2>
        <p>Canonical family assignment prefers InterPro family annotations when they are specific and informative. HGNC and Ensembl-derived paralog references are used to recover biologically relevant families when InterPro is incomplete or too broad. If multiple families are available, smaller and more specific families are preferred over large generic families.</p>
        <p>Large families can be skipped from detailed matrix rendering to keep pages readable. For each globally aligned pair, directional identity divides the number of identical residues by the length of the row protein's sequence region. Reversing the comparison changes the denominator when the lengths differ.</p>
      </article>
      <article class="card doc-card">
        <h2>Orthologs and homology</h2>
        <p>{ortholog_text}</p>
        <p>{homolog_projection_text}</p>
      </article>
      <article class="card doc-card">
        <h2>Experimental structures</h2>
        <p>PDB/RCSB evidence is mapped back to UniProt residue numbering where possible. Experimental construct boundaries that overlap the design scope are summarized and can become PDB-backed construct suggestions. PDB entries fully outside the mature secreted, extracellular, or relevant membrane-protein scope are ignored for construct recommendation.</p>
        <p>Use the linked PDB record to review the mapped chain and original experiment.</p>
      </article>
    </section>

    <section class="card doc-card">
      <h2>BLAST sequence similarity and cross-reactivity context</h2>
      <p>Local <code>blastp</code> searches provide sequence similarity context for reviewing possible cross-reactivity. {blast_database_text}</p>
      <div class="table-scroll">
        <table>
          <thead><tr><th>Reported field</th><th>Meaning</th><th>Interpretation</th></tr></thead>
          <tbody>
            <tr><td>Bit score</td><td>BLAST alignment score normalized for database-independent comparison.</td><td>Primary ranking field. Higher values indicate stronger sequence similarity.</td></tr>
            <tr><td>E-value</td><td>Expected number of alignments of similar quality by chance.</td><td>Lower values indicate more statistically significant alignments.</td></tr>
            <tr><td>Identity</td><td>Fraction of identical residues across the aligned region.</td><td>Useful for local similarity, but must be interpreted with coverage.</td></tr>
            <tr><td>Coverage</td><td>Aligned query span divided by query length, capped at 100%.</td><td>Distinguishes full-region similarity from short conserved motifs.</td></tr>
            <tr><td>Alignment text</td><td>Original BLAST alignment block.</td><td>Inspect this to determine whether similarity overlaps the selected construct or only a small motif.</td></tr>
          </tbody>
        </table>
      </div>
      <p>Self hits are filtered by accession and entry name. BLAST results are ranked sequence similarity signals, not measurements of antibody cross-reactivity. No universal threshold is applied because acceptable risk depends on antigen class, screening goal, and downstream assay design.</p>
    </section>

    <section class="doc-grid">
      <article class="card doc-card">
        <h2>PTMs and processing features</h2>
        <p>OpenAntigens extracts curated UniProt PTM-like features including glycosylation, modified residues, lipidation, disulfide bonds, cross-links, initiator methionine processing, signal peptides, propeptides, mature chains, peptides, and cleavage-related site annotations.</p>
        <p>Target pages list available annotations in a PTM section. Construct cards list which PTMs overlap each construct, and the Interactive Construct Builder updates the included PTMs in real time as the selected residue window changes. Processing annotations are included in boundary warnings.</p>
      </article>
      <article class="card doc-card">
        <h2>Cysteine analysis</h2>
        <p>Model pairing uses cysteine sulfur atoms within 2.4 angstroms. A pair is assigned only when each cysteine has exactly one candidate partner and the choices are reciprocal. Competing contacts are marked ambiguous. Curated UniProt disulfides are shown separately as annotations.</p>
        <p>The Interactive Construct Builder can optionally apply Cys-to-Ser edits to selected unpaired cysteines and updates copied TSV/FASTA sequences in real time.</p>
      </article>
      <article class="card doc-card">
        <h2>Furin-site scanning</h2>
        <p>The default scan finds R-X-[K/R]-R motifs, including overlapping matches, in the design region and constructs. The pipeline reports suggested edits as annotations and leaves application to manual review because cleavage risk is context-dependent and some basic motifs may be biologically important.</p>
      </article>
      <article class="card doc-card">
        <h2>Ligand and interaction regions</h2>
        <p>Curated UniProt features, PDB mappings, and annotation text are used to report ligand-binding or interaction-related regions when available. Check these annotations before trimming the sequence.</p>
      </article>
      <article class="card doc-card">
        <h2>Assembly and partner requirements</h2>
        <p>Interaction and complex information is summarized separately from obligatory-partner classification. Known interactions provide context; obligatory-partner warnings require stronger curated complex evidence.</p>
        <p>Obligatory-partner warnings are anchored to curated Complex Portal records when that module is enabled. UniProt <code>SUBUNIT</code> comments are retained as interaction and assembly context without independently creating obligatory-partner warnings.</p>
        <div class="table-scroll">
          <table>
            <thead><tr><th>Class</th><th>Evidence rule</th><th>How to use it</th></tr></thead>
            <tbody>
              <tr><td>Obligatory partner requirement</td><td>Human Complex Portal evidence supports a biologically relevant mandatory partner relationship, currently including curated integrin alpha/beta complexes.</td><td>Plan co-expression, partner-chain inclusion, or extra validation before treating an isolated construct as biologically complete.</td></tr>
              <tr><td>Mouse complex support</td><td>Mouse Complex Portal evidence supports a comparable ortholog complex. Mouse-only evidence stays contextual for the human report.</td><td>Cross-species context for reagent strategy and validation planning.</td></tr>
              <tr><td>Conditional assembly</td><td>Language such as ligand-induced, upon ligand, or ligand-triggered dimerization.</td><td>Context for conditional biology; partner inclusion depends on the assay and construct goal.</td></tr>
              <tr><td>Stable complex context</td><td>Curated human Complex Portal membership without mandatory-partner status for antigen design.</td><td>Biological context for deciding whether to co-express partners.</td></tr>
              <tr><td>Interaction context</td><td>Generic phrases such as interacts with or forms a complex with.</td><td>Interaction awareness without mandatory-partner status.</td></tr>
            </tbody>
          </table>
        </div>
        <p>A listed interaction alone does not produce an obligatory-partner warning. Check the evidence class and source record before deciding whether to include a partner.</p>
      </article>
    </section>

    <section class="card doc-card">
      <h2>Disease, literature, and antibody-resource context</h2>
      <p>{disease_methods}</p>
      <p>PubTator3 is queried for approximate gene-linked PubMed hit counts. The counts help locate literature; they do not grade its quality or validate a target. CiteAb links are provided as external antibody-resource shortcuts; CiteAb catalog content remains external.</p>
    </section>

    <section class="card doc-card">
      <h2>Static assets and visualization</h2>
      <p>When assets are generated, each construct can include structure images, pLDDT plots, and PAE plots with the construct region highlighted. Structure renderings use AlphaFold models colored by pLDDT and emphasize cysteines when possible. The portal also includes an Interactive Construct Builder that links sequence selection, structure highlighting, pLDDT, PAE, warnings, and copy-ready sequence outputs.</p>
      <p>If the interactive builder is unavailable, the precomputed cards still provide sequences and any saved plots.</p>
    </section>

    <section class="card doc-card">
      <h2>Update cadence and limitations</h2>
      <p>OpenAntigens is released as a static, versioned snapshot. Each release records the build date, target universe, and OpenAntigens version when available. Future releases may incorporate updated annotations.</p>
      <p>Because each release is static, source databases may have changed after the displayed release date. Cite the release date and downloaded artifact names when using OpenAntigens data in downstream analyses.</p>
      <div class="table-scroll">
        <table>
          <thead><tr><th>Limitation</th><th>Consequence</th></tr></thead>
          <tbody>
            <tr><td>AlphaFold confidence is computational.</td><td>High pLDDT/low PAE supports design reasoning. Expression, folding, and epitope presentation require experimental evidence; experimental structures remain separate evidence.</td></tr>
            <tr><td>Topology annotations can be incomplete or inconsistent.</td><td>Mature secreted and extracellular-region boundaries require manual review for important targets.</td></tr>
            <tr><td>BLAST is sequence-based.</td><td>It cannot fully predict antibody cross-reactivity, conformational epitope overlap, glycosylation effects, or cell-surface accessibility.</td></tr>
            <tr><td>Family assignment can be imperfect.</td><td>Some proteins have broad, overlapping, or missing family annotations; paralog matrices provide context and require biological review.</td></tr>
            <tr><td>Ortholog mapping depends on available canonical sequences.</td><td>Missing or predicted RefSeq records can limit mouse or cynomolgus monkey construct transfer.</td></tr>
            <tr><td>Construct heuristics are general-purpose.</td><td>Known biology, literature constructs, ligand-binding requirements, and assay-specific constraints can override automated suggestions.</td></tr>
          </tbody>
        </table>
      </div>
    </section>
""",
    )


def render_downloads_page(
    entries: list[dict[str, Any]],
    *,
    include_disease_context: bool = True,
    disease_downloads: tuple[str, ...] = (),
    portal_title: str = "OpenAntigens",
) -> str:
    total = len(entries)
    ok = sum(1 for item in entries if item.get("status") == "ok")
    errors = sum(1 for item in entries if item.get("status") == "error")
    is_mouse = portal_title == "OpenAntigens Mouse"
    open_targets_buttons = "".join(
        f'<a class="button" href="open_targets_disease_associations.{suffix}">Open Targets {suffix.upper()}</a>'
        for suffix in ("tsv", "json") if include_disease_context and suffix in disease_downloads
    )
    disease_schema_row = (
        '<tr><td><code>top_disease_name</code>, <code>top_disease_score</code>, <code>disease_count</code></td><td>Top Open Targets indirect disease association in the index and number of downloaded associations for the target. Per-target report pages also list direct overall scores.</td></tr>'
        if include_disease_context
        else ""
    )
    open_targets_note = (
        "<li>Open Targets TSV/JSON files contain downloaded indirect and direct target-disease association scores when included in the release.</li>"
        if include_disease_context
        else ""
    )
    mouse_schema_row = (
        "<tr><td><code>mouse_ortholog_*</code>, <code>mouse_refseq_accession</code>, <code>mouse_sequence_source</code></td><td>Mouse target identifiers and sequence provenance.</td></tr>"
        "<tr><td><code>human_source_*</code></td><td>Identifiers for the originating human OpenAntigens target.</td></tr>"
        if is_mouse
        else "<tr><td><code>mouse_ortholog_*</code>, <code>mouse_refseq_accession</code>, <code>mouse_sequence_source</code></td><td>Mouse ortholog identifiers and sequence provenance when present; blank in human rows without mouse-derived metadata.</td></tr>"
        "<tr><td><code>human_source_*</code></td><td>Human source identifiers for derived portals or provenance-aware exports when present.</td></tr>"
    )
    return _render_info_page(
        title="Downloads",
        active="downloads",
        eyebrow="Data access",
        heading="Download static release data",
        intro=(
            "OpenAntigens release snapshots include static flat files for the processed target index, release metadata, and Open Targets disease associations when available."
        ),
        portal_title=portal_title,
        body_html=f"""
    <section class="download-grid">
      <article class="card doc-card">
        <h2>Release flat files</h2>
      <p>For research use, please cite {citation_short_html()} (preprint) and record the release used. <a class="inline-link" href="help.html#cite-openantigens">How to cite</a> · {_citation_download_links()}</p>
        <p>Current index snapshot: {total:,} portal rows, {ok:,} successful reports, {errors:,} error rows.</p>
        <div class="link-group link-group-wrap">
          <a class="button" href="downloads/agdesign2_portal_index.tsv">Download TSV index</a>
          <a class="button" href="downloads/agdesign2_portal_index.json">Download JSON index</a>
          <a class="button" href="downloads/download_manifest.json">Download manifest</a>
          <a class="button" href="portal_metadata.json">Portal metadata</a>
          {open_targets_buttons}
        </div>
      </article>
    </section>

    <section class="card doc-card">
      <h2>Download contents</h2>
      <p>These files contain the current portal index and available release data. Open Targets links appear only when those files are included.</p>
        <ul class="doc-list">
          <li><code>agdesign2_portal_index.tsv</code> is the compact index for spreadsheet and scripting workflows.</li>
          <li><code>agdesign2_portal_index.json</code> contains the same compact index as JSON objects.</li>
          <li><code>download_manifest.json</code> describes included release files, formats, source notes, and index columns.</li>
          <li><code>portal_metadata.json</code> records release-level metadata such as build date, row counts, status counts, and topology counts.</li>
          {open_targets_note}
        </ul>
    </section>
    <section class="card doc-card">
      <h2>TSV schema</h2>
      <p>The TSV index is compact by design. Use the rendered target pages for construct-level evidence and the release index for programmatic analysis.</p>
      <div class="table-scroll">
        <table>
          <thead><tr><th>Column</th><th>Meaning</th></tr></thead>
          <tbody>
            <tr><td><code>entry_name</code>, <code>gene_symbol</code>, <code>protein_name</code></td><td>Resolved portal target identifiers.</td></tr>
            {mouse_schema_row}
            {disease_schema_row}
            <tr><td><code>status</code>, <code>error</code></td><td>Batch processing state and failure message when present.</td></tr>
            <tr><td><code>topology_bucket</code>, <code>ectodomain</code></td><td>Portal track and mature secreted or extracellular design boundary.</td></tr>
            <tr><td><code>construct_count</code>, <code>cross_reactivity_count</code></td><td>Number of reported constructs and BLAST hits used in the portal index.</td></tr>
            <tr><td><code>has_alphafold_structure</code>, <code>has_pdb</code>, <code>has_interpro</code></td><td>Availability flags for major evidence classes.</td></tr>
            <tr><td><code>detail_page</code></td><td>Rendered report page path within the static release.</td></tr>
          </tbody>
        </table>
      </div>
    </section>
""",
    )


def render_calculator_page(*, portal_title: str = "OpenAntigens") -> str:
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Protein Calculator | {escape(portal_title)}</title>
{_portal_favicon_links()}
  <link rel="stylesheet" href="portal.css?v={_PORTAL_CSS_VERSION}">
{_portal_analytics_head()}
</head>
<body>
  <main class="page">
    <header class="site-hero site-hero-doc">
      """ + _portal_brand_header(active="calculator", portal_title=portal_title) + """
      <div class="hero-grid">
        <div class="hero-copy">
          <p class="eyebrow">Bench calculator</p>
          <h1>Protein concentration calculator</h1>
          <p class="hero-text">Enter the protein's molecular weight, then a concentration or amount to convert between mass and moles. Add a value from a second category, such as volume, to calculate the remaining quantities.</p>
        </div>
        <div class="hero-panel calculator-hero-panel">
          <span>Molecular weight is required.</span>
          <span>Calculations run in your browser.</span>
        </div>
      </div>
    </header>

    <section class="calculator-layout" aria-label="Protein concentration calculator">
      <article class="card calculator-card">
        <div class="calculator-card-header">
          <h2>Inputs and calculated values</h2>
          <div class="calc-status waiting" id="calcStatus">Waiting for values</div>
        </div>

        <div class="calculator-form">
          <div class="calculator-row">
            <label for="molecularWeight">
              Molecular weight <span class="required-tag">required</span>
              <span>Enter a positive value and check the unit. The default is kDa.</span>
            </label>
            <input id="molecularWeight" data-kind="mw" type="number" min="0" step="any" inputmode="decimal" placeholder="e.g. 50">
            <select id="molecularWeightUnit" aria-label="Molecular weight unit">
              <option value="kDa" selected>kDa</option>
              <option value="Da">Da / g mol^-1</option>
            </select>
          </div>

          <div class="calculator-row">
            <label for="massConcentration">
              Mass concentration
              <span>Optional. With molecular weight, this fills molar concentration.</span>
            </label>
            <input id="massConcentration" data-kind="massConcentration" type="number" min="0" step="any" inputmode="decimal" placeholder="optional">
            <select id="massConcentrationUnit" data-kind="massConcentration" aria-label="Mass concentration unit">
              <option value="mg_per_ml">mg/mL</option>
              <option value="ug_per_ml">ug/mL</option>
            </select>
          </div>

          <div class="calculator-row">
            <label for="molarConcentration">
              Molar concentration
              <span>Optional. With molecular weight, this fills mass concentration.</span>
            </label>
            <input id="molarConcentration" data-kind="molarConcentration" type="number" min="0" step="any" inputmode="decimal" placeholder="optional">
            <select id="molarConcentrationUnit" data-kind="molarConcentration" aria-label="Molar concentration unit">
              <option value="M">M</option>
              <option value="mM">mM</option>
              <option value="uM" selected>uM</option>
              <option value="nM">nM</option>
              <option value="pM">pM</option>
            </select>
          </div>

          <div class="calculator-row">
            <label for="massQuantity">
              Total quantity by mass
              <span>Optional. With molecular weight, this fills total moles.</span>
            </label>
            <input id="massQuantity" data-kind="massQuantity" type="number" min="0" step="any" inputmode="decimal" placeholder="optional">
            <select id="massQuantityUnit" data-kind="massQuantity" aria-label="Mass quantity unit">
              <option value="g">g</option>
              <option value="mg">mg</option>
              <option value="ug" selected>ug</option>
              <option value="ng">ng</option>
            </select>
          </div>

          <div class="calculator-row">
            <label for="moleQuantity">
              Total quantity by moles
              <span>Optional. With molecular weight, this fills total mass.</span>
            </label>
            <input id="moleQuantity" data-kind="moleQuantity" type="number" min="0" step="any" inputmode="decimal" placeholder="optional">
            <select id="moleQuantityUnit" data-kind="moleQuantity" aria-label="Mole quantity unit">
              <option value="mol">mol</option>
              <option value="mmol">mmol</option>
              <option value="umol">umol</option>
              <option value="nmol" selected>nmol</option>
              <option value="pmol">pmol</option>
            </select>
          </div>

          <div class="calculator-row">
            <label for="volume">
              Total volume
              <span>Final solution volume.</span>
            </label>
            <input id="volume" data-kind="volume" type="number" min="0" step="any" inputmode="decimal" placeholder="optional">
            <select id="volumeUnit" data-kind="volume" aria-label="Volume unit">
              <option value="ml">mL</option>
              <option value="ul" selected>uL</option>
            </select>
          </div>
        </div>

        <div class="calculator-actions">
          <button type="button" class="mini-button" id="clearCalculated">Clear calculated</button>
          <button type="button" class="button" id="resetAll">Reset</button>
        </div>
      </article>

      <aside class="card calculator-card">
        <h2>Calculation state</h2>
        <div class="calc-message" id="calcMessage">
          <strong>Molecular weight is required.</strong> Then enter one concentration or quantity to infer its equivalent unit, or add volume to solve the full concentration/quantity/volume set.
        </div>
        <dl class="calculator-facts" aria-label="Normalized calculation values">
          <dt>Mass concentration</dt><dd id="factMassConc">-</dd>
          <dt>Molar concentration</dt><dd id="factMolarConc">-</dd>
          <dt>Total mass</dt><dd id="factMass">-</dd>
          <dt>Total moles</dt><dd id="factMoles">-</dd>
          <dt>Total volume</dt><dd id="factVolume">-</dd>
        </dl>
        <ul class="doc-list calculator-notes">
          <li>Enter molecular weight first. Use any two categories to solve the full calculation: concentration, amount, or volume.</li>
          <li>Fields you type into are source values. Automatically filled fields use a pale background.</li>
          <li>Changing a unit converts the displayed value so the underlying amount stays constant.</li>
          <li>If you enter both mass and molar values in one category, the latest edit is used. If you supply concentration, amount, and volume, they must agree.</li>
        </ul>
      </aside>
    </section>
    """ + _portal_footer() + """
  </main>

  <script>
    const fields = {
      mw: document.getElementById("molecularWeight"),
      massConcentration: document.getElementById("massConcentration"),
      molarConcentration: document.getElementById("molarConcentration"),
      massQuantity: document.getElementById("massQuantity"),
      moleQuantity: document.getElementById("moleQuantity"),
      volume: document.getElementById("volume")
    };

    const units = {
      mw: document.getElementById("molecularWeightUnit"),
      massConcentration: document.getElementById("massConcentrationUnit"),
      molarConcentration: document.getElementById("molarConcentrationUnit"),
      massQuantity: document.getElementById("massQuantityUnit"),
      moleQuantity: document.getElementById("moleQuantityUnit"),
      volume: document.getElementById("volumeUnit")
    };

    const statusEl = document.getElementById("calcStatus");
    const messageEl = document.getElementById("calcMessage");
    const sourceKinds = ["massConcentration", "molarConcentration", "massQuantity", "moleQuantity", "volume"];
    const sourceGroup = {
      massConcentration: "concentration",
      molarConcentration: "concentration",
      massQuantity: "quantity",
      moleQuantity: "quantity",
      volume: "volume"
    };
    const precision = 6;
    let userSources = new Set();
    let lastEdited = [];
    let isProgrammatic = false;
    let previousUnits = {};

    const massConcentrationFactors = {
      mg_per_ml: 1,
      ug_per_ml: 0.001
    };
    const molecularWeightFactors = {
      Da: 1,
      kDa: 1000
    };
    const molarConcentrationFactors = {
      M: 1,
      mM: 1e-3,
      uM: 1e-6,
      nM: 1e-9,
      pM: 1e-12
    };
    const massQuantityFactors = {
      g: 1,
      mg: 1e-3,
      ug: 1e-6,
      ng: 1e-9
    };
    const moleQuantityFactors = {
      mol: 1,
      mmol: 1e-3,
      umol: 1e-6,
      nmol: 1e-9,
      pmol: 1e-12
    };
    const volumeFactors = {
      ml: 1,
      ul: 0.001
    };

    Object.keys(units).forEach((kind) => {
      previousUnits[kind] = units[kind].value;
    });

    function parsePositiveNumber(input) {
      const value = Number(input.value);
      return Number.isFinite(value) && value > 0 ? value : null;
    }

    function parseMolecularWeightDa() {
      const value = parsePositiveNumber(fields.mw);
      if (value === null) return null;
      return value * molecularWeightFactors[units.mw.value];
    }

    function toBase(kind, value) {
      if (value === null) return null;
      if (kind === "mw") return value * molecularWeightFactors[units.mw.value];
      if (kind === "massConcentration") return value * massConcentrationFactors[units.massConcentration.value];
      if (kind === "molarConcentration") return value * molarConcentrationFactors[units.molarConcentration.value];
      if (kind === "massQuantity") return value * massQuantityFactors[units.massQuantity.value];
      if (kind === "moleQuantity") return value * moleQuantityFactors[units.moleQuantity.value];
      if (kind === "volume") return value * volumeFactors[units.volume.value];
      return value;
    }

    function fromBase(kind, value) {
      if (!Number.isFinite(value) || value <= 0) return "";
      if (kind === "massConcentration") return value / massConcentrationFactors[units.massConcentration.value];
      if (kind === "molarConcentration") return value / molarConcentrationFactors[units.molarConcentration.value];
      if (kind === "massQuantity") return value / massQuantityFactors[units.massQuantity.value];
      if (kind === "moleQuantity") return value / moleQuantityFactors[units.moleQuantity.value];
      if (kind === "volume") return value / volumeFactors[units.volume.value];
      return value;
    }

    function valueIsPresent(value) {
      return Number.isFinite(value) && value > 0;
    }

    function formatInput(value) {
      if (!Number.isFinite(value) || value <= 0) return "";
      const abs = Math.abs(value);
      if (abs >= 1e5 || abs < 1e-4) {
        return Number(value.toExponential(precision - 1)).toString();
      }
      return Number(value.toPrecision(precision)).toString();
    }

    function formatFact(value, unit) {
      if (!Number.isFinite(value) || value <= 0) return "-";
      const abs = Math.abs(value);
      const displayed = abs >= 1e5 || abs < 1e-3
        ? value.toExponential(4)
        : Number(value.toPrecision(5)).toString();
      return `${displayed} ${unit}`;
    }

    function setCalculated(kind, baseValue) {
      const input = fields[kind];
      input.value = formatInput(fromBase(kind, baseValue));
      input.classList.toggle("calculated", input.value !== "");
      input.classList.remove("user-set");
    }

    function clearCalculatedFields() {
      sourceKinds.forEach((kind) => {
        if (!userSources.has(kind)) {
          fields[kind].value = "";
          fields[kind].classList.remove("calculated", "user-set");
        }
      });
    }

    function markUserFields() {
      sourceKinds.forEach((kind) => {
        const isSource = userSources.has(kind);
        fields[kind].classList.toggle("user-set", isSource);
        if (isSource) fields[kind].classList.remove("calculated");
      });
    }

    function updateMessage(text, mode) {
      messageEl.innerHTML = text;
      messageEl.classList.remove("error", "ready");
      statusEl.classList.remove("ready", "waiting");
      if (mode === "error") {
        messageEl.classList.add("error");
        statusEl.textContent = "Check inputs";
        statusEl.classList.add("waiting");
      } else if (mode === "ready") {
        messageEl.classList.add("ready");
        statusEl.textContent = "Calculated";
        statusEl.classList.add("ready");
      } else {
        statusEl.textContent = "Waiting for values";
        statusEl.classList.add("waiting");
      }
    }

    function resetFacts() {
      document.getElementById("factMassConc").textContent = "-";
      document.getElementById("factMolarConc").textContent = "-";
      document.getElementById("factMass").textContent = "-";
      document.getElementById("factMoles").textContent = "-";
      document.getElementById("factVolume").textContent = "-";
    }

    function updateFacts(base) {
      document.getElementById("factMassConc").textContent = formatFact(base.massConcentration, "mg/mL");
      document.getElementById("factMolarConc").textContent = formatFact(base.molarConcentration, "M");
      document.getElementById("factMass").textContent = formatFact(base.massQuantity, "g");
      document.getElementById("factMoles").textContent = formatFact(base.moleQuantity, "mol");
      document.getElementById("factVolume").textContent = formatFact(base.volume, "mL");
    }

    function latestSourceByGroup(groups) {
      for (const kind of lastEdited) {
        if (userSources.has(kind) && groups.includes(sourceGroup[kind])) return kind;
      }
      return null;
    }

    function baseForSource(kind, mw) {
      const raw = parsePositiveNumber(fields[kind]);
      const value = toBase(kind, raw);
      if (value === null) return null;

      if (kind === "massConcentration") {
        return { massConcentration: value, molarConcentration: value / mw };
      }
      if (kind === "molarConcentration") {
        return { molarConcentration: value, massConcentration: value * mw };
      }
      if (kind === "massQuantity") {
        return { massQuantity: value, moleQuantity: value / mw };
      }
      if (kind === "moleQuantity") {
        return { moleQuantity: value, massQuantity: value * mw };
      }
      if (kind === "volume") {
        return { volume: value };
      }
      return {};
    }

    function syncEquivalentUserField(sourceKind, targetKind, baseValue) {
      if (userSources.has(targetKind) && !userSources.has(sourceKind)) return;
      if (sourceGroup[sourceKind] !== sourceGroup[targetKind]) return;
      if (!Number.isFinite(baseValue) || baseValue <= 0) return;
      const sourceIndex = lastEdited.indexOf(sourceKind);
      const targetIndex = lastEdited.indexOf(targetKind);
      if (targetIndex !== -1 && targetIndex < sourceIndex) return;
      setCalculated(targetKind, baseValue);
      userSources.delete(targetKind);
    }

    function solve() {
      if (isProgrammatic) return;
      isProgrammatic = true;

      const mw = parseMolecularWeightDa();
      clearCalculatedFields();
      resetFacts();
      markUserFields();

      if (Object.values(fields).some((input) => input.validity.badInput || (input.value !== "" && parsePositiveNumber(input) === null))) {
        updateMessage("<strong>Invalid input.</strong> Enter positive numbers or clear the invalid field.", "error");
        isProgrammatic = false;
        return;
      }

      if (mw === null) {
        updateMessage("<strong>Molecular weight is required.</strong> Enter a positive molecular weight before calculating.", "idle");
        isProgrammatic = false;
        return;
      }

      const validSources = [...userSources].filter((kind) => parsePositiveNumber(fields[kind]) !== null);
      userSources = new Set(validSources);
      lastEdited = lastEdited.filter((kind) => userSources.has(kind));
      markUserFields();

      const concentrationSource = latestSourceByGroup(["concentration"]);
      const quantitySource = latestSourceByGroup(["quantity"]);
      const volumeSource = latestSourceByGroup(["volume"]);
      const presentGroups = new Set(validSources.map((kind) => sourceGroup[kind]));

      let base = {
        massConcentration: null,
        molarConcentration: null,
        massQuantity: null,
        moleQuantity: null,
        volume: null
      };

      if (concentrationSource) Object.assign(base, baseForSource(concentrationSource, mw));
      if (quantitySource) Object.assign(base, baseForSource(quantitySource, mw));
      if (volumeSource) Object.assign(base, baseForSource(volumeSource, mw));

      if (valueIsPresent(base.massConcentration) && !valueIsPresent(base.molarConcentration)) base.molarConcentration = base.massConcentration / mw;
      if (valueIsPresent(base.molarConcentration) && !valueIsPresent(base.massConcentration)) base.massConcentration = base.molarConcentration * mw;
      if (valueIsPresent(base.massQuantity) && !valueIsPresent(base.moleQuantity)) base.moleQuantity = base.massQuantity / mw;
      if (valueIsPresent(base.moleQuantity) && !valueIsPresent(base.massQuantity)) base.massQuantity = base.moleQuantity * mw;

      const hasAnyEquivalentPair =
        (valueIsPresent(base.massConcentration) && valueIsPresent(base.molarConcentration)) ||
        (valueIsPresent(base.massQuantity) && valueIsPresent(base.moleQuantity));

      if (validSources.length < 2 || presentGroups.size < 2) {
        sourceKinds.forEach((kind) => {
          if (!userSources.has(kind)) setCalculated(kind, base[kind]);
        });
        markUserFields();
        updateFacts(base);
        if (hasAnyEquivalentPair) {
          updateMessage("<strong>Equivalent value calculated.</strong> Add volume or another value type to solve concentration, quantity, and volume together.", "ready");
        } else {
          updateMessage("<strong>Add a source value.</strong> Enter one concentration or quantity to infer its equivalent unit, or add two compatible value types for a full solution calculation.", "idle");
        }
        isProgrammatic = false;
        return;
      }

      if (presentGroups.size === 3) {
        const expectedMass = base.massConcentration * base.volume / 1000;
        if (Math.abs(expectedMass - base.massQuantity) > 1e-6 * Math.max(expectedMass, base.massQuantity)) {
          updateMessage("<strong>Inputs are inconsistent.</strong> Concentration multiplied by volume must equal mass. Correct or clear one of the three source categories.", "error");
          isProgrammatic = false;
          return;
        }
      }

      if (valueIsPresent(base.massConcentration) && valueIsPresent(base.volume) && !valueIsPresent(base.massQuantity)) {
        base.massQuantity = base.massConcentration * base.volume / 1000;
        base.moleQuantity = base.massQuantity / mw;
      }

      if (valueIsPresent(base.molarConcentration) && valueIsPresent(base.volume) && !valueIsPresent(base.moleQuantity)) {
        base.moleQuantity = base.molarConcentration * base.volume / 1000;
        base.massQuantity = base.moleQuantity * mw;
      }

      if (valueIsPresent(base.massQuantity) && valueIsPresent(base.volume) && !valueIsPresent(base.massConcentration)) {
        base.massConcentration = base.massQuantity * 1000 / base.volume;
        base.molarConcentration = base.massConcentration / mw;
      }

      if (valueIsPresent(base.moleQuantity) && valueIsPresent(base.volume) && !valueIsPresent(base.molarConcentration)) {
        base.molarConcentration = base.moleQuantity * 1000 / base.volume;
        base.massConcentration = base.molarConcentration * mw;
      }

      if (valueIsPresent(base.massQuantity) && valueIsPresent(base.massConcentration) && !valueIsPresent(base.volume)) {
        base.volume = base.massQuantity * 1000 / base.massConcentration;
      }

      if (valueIsPresent(base.moleQuantity) && valueIsPresent(base.molarConcentration) && !valueIsPresent(base.volume)) {
        base.volume = base.moleQuantity * 1000 / base.molarConcentration;
      }

      if (valueIsPresent(base.massConcentration) && !valueIsPresent(base.molarConcentration)) base.molarConcentration = base.massConcentration / mw;
      if (valueIsPresent(base.molarConcentration) && !valueIsPresent(base.massConcentration)) base.massConcentration = base.molarConcentration * mw;
      if (valueIsPresent(base.massQuantity) && !valueIsPresent(base.moleQuantity)) base.moleQuantity = base.massQuantity / mw;
      if (valueIsPresent(base.moleQuantity) && !valueIsPresent(base.massQuantity)) base.massQuantity = base.moleQuantity * mw;

      if (!valueIsPresent(base.massConcentration) || !valueIsPresent(base.molarConcentration) || !valueIsPresent(base.massQuantity) || !valueIsPresent(base.moleQuantity) || !valueIsPresent(base.volume)) {
        updateMessage("<strong>Calculation is incomplete.</strong> Check that all source values are positive numbers.", "error");
        isProgrammatic = false;
        return;
      }

      sourceKinds.forEach((kind) => {
        if (!userSources.has(kind)) setCalculated(kind, base[kind]);
      });

      if (concentrationSource === "massConcentration") syncEquivalentUserField("massConcentration", "molarConcentration", base.molarConcentration);
      if (concentrationSource === "molarConcentration") syncEquivalentUserField("molarConcentration", "massConcentration", base.massConcentration);
      if (quantitySource === "massQuantity") syncEquivalentUserField("massQuantity", "moleQuantity", base.moleQuantity);
      if (quantitySource === "moleQuantity") syncEquivalentUserField("moleQuantity", "massQuantity", base.massQuantity);

      markUserFields();
      updateFacts(base);
      updateMessage("<strong>Values calculated.</strong> The latest input in each category is used.", "ready");
      isProgrammatic = false;
    }

    function rememberSource(kind) {
      if (!sourceKinds.includes(kind)) return;
      if (fields[kind].value === "" && !fields[kind].validity.badInput) {
        userSources.delete(kind);
        lastEdited = lastEdited.filter((item) => item !== kind);
      } else {
        userSources.add(kind);
        lastEdited = [kind, ...lastEdited.filter((item) => item !== kind)];
      }
    }

    fields.mw.addEventListener("input", solve);

    sourceKinds.forEach((kind) => {
      fields[kind].addEventListener("input", () => {
        rememberSource(kind);
        solve();
      });
    });

    Object.entries(units).forEach(([kind, select]) => {
      select.addEventListener("change", () => {
        const input = fields[kind];
        const oldUnit = previousUnits[kind];
        const newUnit = select.value;
        const value = Number(input.value);
        const factorMap = kind === "mw"
          ? molecularWeightFactors
          : kind === "massConcentration"
          ? massConcentrationFactors
          : kind === "molarConcentration"
            ? molarConcentrationFactors
            : kind === "massQuantity"
              ? massQuantityFactors
              : kind === "moleQuantity"
                ? moleQuantityFactors
                : volumeFactors;

        if (Number.isFinite(value) && value > 0) {
          const baseValue = value * factorMap[oldUnit];
          input.value = formatInput(baseValue / factorMap[newUnit]);
        }
        previousUnits[kind] = newUnit;
        solve();
      });
    });

    document.getElementById("clearCalculated").addEventListener("click", () => {
      clearCalculatedFields();
      resetFacts();
      updateMessage("<strong>Calculated fields cleared.</strong> Your source values are still in place.", "idle");
    });

    document.getElementById("resetAll").addEventListener("click", () => {
      isProgrammatic = true;
      Object.values(fields).forEach((input) => {
        input.value = "";
        input.classList.remove("calculated", "user-set");
      });
      userSources = new Set();
      lastEdited = [];
      resetFacts();
      updateMessage("<strong>Molecular weight is required.</strong> Then enter one concentration or quantity to infer its equivalent unit, or add volume to solve the full concentration/quantity/volume set.", "idle");
      isProgrammatic = false;
    });

    solve();
  </script>
</body>
</html>
"""


def render_terms_page(*, portal_title: str = "OpenAntigens") -> str:
    return _render_info_page(
        title="License and Terms",
        active="terms",
        eyebrow="Use conditions",
        heading="Terms for OpenAntigens data and software",
        intro=(
            "OpenAntigens combines public biological data with "
            "OpenAntigens-generated computational annotations. Separate licenses apply to the software, generated "
            "annotations, and third-party source data."
        ),
        portal_title=portal_title,
        body_html=f"""
    <section class="card doc-card">
      <h2>Summary</h2>
      <div class="table-scroll">
        <table>
          <thead><tr><th>Material</th><th>Terms</th><th>Notes</th></tr></thead>
          <tbody>
            <tr><td>OpenAntigens software</td><td>Apache License 2.0</td><td>Applies to the project source code unless a file states otherwise.</td></tr>
            <tr><td>OpenAntigens-generated annotations</td><td>CC BY 4.0</td><td>Applies to OpenAntigens-created construct recommendations, classifications, computed summaries, derived report text, portal index metadata, and related derived annotations.</td></tr>
            <tr><td>Third-party source data</td><td>Original source terms apply</td><td>Sequences, structures, identifiers, disease associations, literature counts, and source annotations retain the terms, licenses, and citation expectations of the databases from which they were obtained.</td></tr>
          </tbody>
        </table>
      </div>
      <p>These terms describe research-data reuse for this release. Legal review and source-database compliance remain the responsibility of the person or organization reusing the data.</p>
    </section>

    <section class="card doc-card">
      <h2>Website analytics</h2>
      <p>OpenAntigens uses aggregate, cookieless analytics to understand how the portal and its downloadable resources are used. The <a class="inline-link" href="privacy.html">Privacy and Analytics page</a> describes what is measured, what is excluded, and how aggregate results may be reported to funders.</p>
    </section>

    <section class="doc-grid">
      <article class="card doc-card">
        <h2>OpenAntigens code</h2>
        <p>The OpenAntigens software is released under the <a class="inline-link" href="https://www.apache.org/licenses/LICENSE-2.0" target="_blank" rel="noreferrer">Apache License 2.0</a>. This includes the portal renderer, command-line tools, data-retrieval code, analysis logic, and supporting scripts unless a specific file states otherwise.</p>
      </article>
      <article class="card doc-card">
        <h2>OpenAntigens annotations</h2>
        <p>OpenAntigens-generated annotations are released under <a class="inline-link" href="https://creativecommons.org/licenses/by/4.0/" target="_blank" rel="noreferrer">Creative Commons Attribution 4.0 International (CC BY 4.0)</a>. Reusers may share and adapt these annotations, including for commercial use, provided appropriate attribution is given.</p>
      </article>
      <article class="card doc-card">
        <h2>Third-party source data</h2>
        <p>OpenAntigens reports may include or derive from third-party identifiers, protein sequences, AlphaFold structures and confidence metrics, PDB mappings, domain annotations, disease associations, literature counts, BLAST databases, and source links. Those materials retain the terms and citation requirements of their original providers.</p>
      </article>
      <article class="card doc-card">
        <h2>Research use</h2>
        <p>OpenAntigens is a research-use computational resource. Diagnosis, clinical decision-making, regulatory submissions, therapeutic selection, manufacturability claims, and expression-performance claims require independent expert review and experimental validation.</p>
      </article>
    </section>

    <section class="card doc-card">
      <h2>Source database terms and attribution</h2>
      <p>The table summarizes major sources used by OpenAntigens. Official source terms control reuse requirements.</p>
      <div class="table-scroll">
        <table>
          <thead><tr><th>Source</th><th>Used for</th><th>License / terms summary</th><th>Link</th></tr></thead>
          <tbody>
            <tr><td>UniProt</td><td>Reviewed protein records, sequences, features, topology, domains, PTMs, entry names</td><td>Copyrightable database content is under CC BY 4.0; attribution required.</td><td><a href="https://www.uniprot.org/help/license" target="_blank" rel="noreferrer">UniProt license</a></td></tr>
            <tr><td>AlphaFold DB</td><td>Predicted structures, pLDDT, PAE, structure visualizations</td><td>Data are available under CC BY 4.0; AlphaFold predictions are theoretical models and provided as-is.</td><td><a href="https://alphafold.com/" target="_blank" rel="noreferrer">AlphaFold DB</a></td></tr>
            <tr><td>Thera-SAbDab</td><td>Therapeutic-antibody target evidence retained in the canonical human target table</td><td>Publication and source terms apply; cite the Thera-SAbDab database publication.</td><td><a href="https://doi.org/10.1093/nar/gkz827" target="_blank" rel="noreferrer">Thera-SAbDab publication</a></td></tr>
            <tr><td>SURFY</td><td>Surfaceome support and topology context retained in the canonical human target table</td><td>Publication and source terms apply; cite the in silico human surfaceome publication.</td><td><a href="https://doi.org/10.1073/pnas.1808790115" target="_blank" rel="noreferrer">SURFY publication</a></td></tr>
            <tr><td>Uhlén human secretome</td><td>Secretome support retained in the canonical human target table</td><td>Publication and source terms apply; cite the human secretome publication.</td><td><a href="https://doi.org/10.1126/scisignal.aaz0274" target="_blank" rel="noreferrer">Human secretome publication</a></td></tr>
            <tr><td>RCSB PDB / wwPDB</td><td>Experimental structure evidence and construct-boundary precedent</td><td>PDB archive and RCSB API data are available under CC0 1.0; attribution to original structure authors is encouraged.</td><td><a href="https://www.rcsb.org/pages/usage-policy" target="_blank" rel="noreferrer">RCSB usage policy</a></td></tr>
            <tr><td>Complex Portal</td><td>Curated human complex membership and assembly context</td><td>EMBL-EBI terms and Complex Portal citation expectations apply.</td><td><a href="https://www.ebi.ac.uk/about/terms-of-use/" target="_blank" rel="noreferrer">EMBL-EBI terms</a></td></tr>
            <tr><td>GPCRdb</td><td>GPCR class/family annotations, segment boundaries, generic residue numbering, conserved motif context, and GPCR-focused interpretation</td><td>GPCRdb source-specific terms and citation expectations apply. Cite the main GPCRdb reference and any page-specific GPCRdb references relevant to reused annotations.</td><td><a href="https://docs.gpcrdb.org/citing.html" target="_blank" rel="noreferrer">GPCRdb citing</a> <a href="https://docs.gpcrdb.org/web_services.html" target="_blank" rel="noreferrer">GPCRdb web services</a></td></tr>
            <tr><td>Open Targets Platform</td><td>Indirect and direct target-disease association scores, disease links, and disease-aware search</td><td>Platform data are marked CC0 1.0; citation of Open Targets is expected in accordance with good scientific practice.</td><td><a href="https://platform-docs.opentargets.org/licence" target="_blank" rel="noreferrer">Open Targets licence</a></td></tr>
            <tr><td>NCBI RefSeq / NCBI Protein / NCBI Gene</td><td>Ortholog protein accessions and sequences in RefSeq-backed reference tables, including mouse and cynomolgus fallback/provenance records</td><td>NCBI places no restrictions on molecular data but cannot transfer rights that may belong to submitters or third parties; acknowledge NCBI/NLM where appropriate.</td><td><a href="https://www.ncbi.nlm.nih.gov/home/about/policies/" target="_blank" rel="noreferrer">NCBI policies</a></td></tr>
            <tr><td>InterPro, Pfam, EMBL-EBI resources</td><td>Domain/family annotations and family context</td><td>EMBL-EBI expects attribution and may apply source-specific terms; Pfam is CC0.</td><td><a href="https://www.ebi.ac.uk/about/terms-of-use/" target="_blank" rel="noreferrer">EMBL-EBI terms</a></td></tr>
            <tr><td>HGNC and Ensembl</td><td>Human gene symbols, family/paralog context, gene links</td><td>Source-specific terms and attribution expectations apply.</td><td><a href="https://www.genenames.org/" target="_blank" rel="noreferrer">HGNC</a> <a href="https://www.ensembl.org/" target="_blank" rel="noreferrer">Ensembl</a></td></tr>
            <tr><td>PubTator3 / PubMed</td><td>Literature-hit counts and literature links</td><td>NCBI/NLM policies apply; article text and abstracts may have separate copyright conditions.</td><td><a href="https://www.ncbi.nlm.nih.gov/research/pubtator3/" target="_blank" rel="noreferrer">PubTator3</a></td></tr>
            <tr><td>CiteAb</td><td>External antibody-search link for each target</td><td>CiteAb catalog content stays external to OpenAntigens downloads and reports.</td><td><a href="https://www.citeab.com/" target="_blank" rel="noreferrer">CiteAb</a></td></tr>
          </tbody>
        </table>
      </div>
    </section>

    <section class="doc-grid">
      <article class="card doc-card">
        <h2>Download scope</h2>
        <p>OpenAntigens downloads contain derived annotations alongside source-derived fields such as accessions, sequences, structure-confidence values, BLAST hit descriptions, disease associations, and links. Original source terms remain attached to source-derived fields included in an OpenAntigens file.</p>
      </article>
      <article class="card doc-card">
        <h2>Recommended scholarly citation</h2>
        <p class="paper-citation">{citation_html()}</p>
        <p>This is a recommended scholarly citation. Data reuse remains governed by the licenses and attribution requirements on this page. See <a class="inline-link" href="help.html#cite-openantigens">How to cite</a> for citation files and release details.</p>
        <h2>Attribution</h2>
        <p>Recommended attribution: <em>OpenAntigens, Institute for Protein Innovation, created by Andre A. R. Teixeira</em>. For publications or redistributed datasets, include the portal build date, downloaded artifact name, OpenAntigens version when available, and the relevant source databases.</p>
      </article>
      <article class="card doc-card">
        <h2>Endorsement</h2>
        <p>References or links to third-party databases, products, antibodies, structures, publications, or disease associations identify source material and external resources. They carry no endorsement by OpenAntigens, the Institute for Protein Innovation, the source database providers, or any listed institution.</p>
      </article>
      <article class="card doc-card">
        <h2>Warranty</h2>
        <p>OpenAntigens is provided as-is and without warranties of accuracy, completeness, non-infringement, merchantability, fitness for a particular purpose, expression performance, antigenicity, developability, or suitability for any specific use.</p>
      </article>
    </section>

    <section class="card doc-card">
      <h2>Limitations and user responsibility</h2>
      <p>Construct recommendations, BLAST sequence similarity rankings, paralog/family matrices, disease associations, pLDDT/PAE interpretation, and cysteine/furin warnings are computational annotations. Independent review, experimental validation, biosafety review, intellectual-property review, and compliance with institutional, funder, journal, clinical, regulatory, and commercial requirements remain the responsibility of the person or organization using the data.</p>
      <p>Results depend on the source versions, sequence mappings, structural models, and search databases used for the release. Check those records when interpreting an unexpected result.</p>
    </section>

    <section class="card doc-card">
      <h2>Corrections, licensing questions, and takedown requests</h2>
      <p>Contact: Andre A. R. Teixeira, Institute for Protein Innovation, <a class="inline-link" href="mailto:andre.teixeira@proteininnovation.org">andre.teixeira@proteininnovation.org</a>.</p>
      <p>Contact us to report incorrect attribution, problematic redistribution of source data, stale or incorrect annotations, broken links, or data requiring correction or removal.</p>
    </section>
""",
    )


def render_privacy_page(*, portal_title: str = "OpenAntigens") -> str:
    return _render_info_page(
        title="Privacy and Analytics",
        active="",
        eyebrow="Privacy",
        heading="How OpenAntigens measures site use",
        intro=(
            "OpenAntigens uses aggregate analytics to maintain the portal, decide what to improve, "
            "and report its reach to current and prospective funders."
        ),
        portal_title=portal_title,
        body_html="""
    <section class="card doc-card">
      <h2>What is measured</h2>
      <p>OpenAntigens uses Plausible Analytics to measure aggregate site use. Plausible processes page paths, approximate visitor counts, referring sites, approximate location, device, browser, and operating-system categories, and interactions with downloadable resources. OpenAntigens records named events when a visitor starts a browser-generated PDB or PNG download or successfully copies a TSV or FASTA export.</p>
    </section>

    <section class="card doc-card">
      <h2>What is not collected</h2>
      <p>The integration does not use analytics cookies, advertising identifiers, logins, session replay, or cross-site tracking. Custom events are limited to a fixed allowlist of names. They do not send site-search terms, protein sequences, selected residues, clipboard contents, names, or email addresses to analytics. OpenAntigens does not use analytics to identify or profile individual visitors.</p>
    </section>

    <section class="card doc-card">
      <h2>How measurements are used</h2>
      <p>We use these measurements to maintain OpenAntigens, decide what to improve, and report aggregate reach and use to current and prospective funders. Counts are approximate and may omit visitors who block analytics. For static TSV and JSON files, analytics record a link click. Custom PDB and PNG events record that the browser initiated an export, while copy events record that the browser reported a successful clipboard write. None of these events shows that the resource was used in an experiment.</p>
    </section>

    <section class="card doc-card">
      <h2>Processing, access, and retention</h2>
      <p>Analytics are processed by <a class="inline-link" href="https://plausible.io/data-policy" target="_blank" rel="noreferrer">Plausible Analytics under its data policy</a>. Dashboard access is limited to authorized OpenAntigens administrators at the Institute for Protein Innovation. Only aggregate summaries are shared with funders; visitor-level data and breakdowns with fewer than five unique visitors are not shared.</p>
      <p>Reviewed quarterly summaries are retained for historical reporting. Provider analytics are retained while OpenAntigens uses Plausible. OpenAntigens will delete the provider-held analytics if the site is removed from the service.</p>
    </section>

    <section class="card doc-card">
      <h2>Questions</h2>
      <p>Questions about OpenAntigens analytics can be sent to Andre A. R. Teixeira at <a class="inline-link" href="mailto:andre.teixeira@proteininnovation.org">andre.teixeira@proteininnovation.org</a>.</p>
    </section>
""",
    )


def _render_info_page(
    *,
    title: str,
    active: str,
    eyebrow: str,
    heading: str,
    intro: str,
    body_html: str,
    portal_title: str = "OpenAntigens",
) -> str:
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{escape(title)} | {escape(portal_title)}</title>
{_portal_favicon_links()}
  <link rel="stylesheet" href="portal.css?v={_PORTAL_CSS_VERSION}">
{_portal_analytics_head()}
</head>
<body>
  <main class="page">
    <header class="site-hero site-hero-doc">
      {_portal_brand_header(active=active, portal_title=portal_title)}
      <div class="hero-grid">
        <div class="hero-copy">
          <p class="eyebrow">{escape(eyebrow)}</p>
          <h1>{escape(heading)}</h1>
          <p class="hero-text">{escape(intro)}</p>
        </div>
      </div>
    </header>
    <div class="doc-layout">
      {body_html}
    </div>
    {_portal_footer()}
  </main>
</body>
</html>
"""


def _portal_download_rows(entries: list[dict[str, Any]], *, include_disease_context: bool = True) -> list[dict[str, Any]]:
    fields = _portal_download_fields(include_disease_context=include_disease_context)
    rows: list[dict[str, Any]] = []
    for item in entries:
        rows.append(
            {
                "batch_index": item.get("batch_index", ""),
                "query": item.get("query", ""),
                "source_column": item.get("source_column", ""),
                "entry_name": item.get("entry_name", ""),
                "gene_symbol": item.get("gene_symbol", ""),
                "protein_name": item.get("protein_name", ""),
                "mouse_ortholog_accession": item.get("mouse_ortholog_accession", ""),
                "mouse_ortholog_entry_name": item.get("mouse_ortholog_entry_name", ""),
                "mouse_ortholog_gene_symbol": item.get("mouse_ortholog_gene_symbol", ""),
                "mouse_refseq_accession": item.get("mouse_refseq_accession", ""),
                "mouse_sequence_source": item.get("mouse_sequence_source", ""),
                "human_source_accession": item.get("human_source_accession", ""),
                "human_source_entry_name": item.get("human_source_entry_name", ""),
                "human_source_gene_symbol": item.get("human_source_gene_symbol", ""),
                "family_name": item.get("family_name", ""),
                "top_disease_name": item.get("top_disease_name", ""),
                "top_disease_score": item.get("top_disease_score", ""),
                "disease_count": item.get("disease_count", ""),
                "status": item.get("status", ""),
                "error": "analysis failed" if item.get("error") else item.get("error", ""),
                "topology_bucket": item.get("topology_bucket", ""),
                "ectodomain": item.get("ectodomain", ""),
                "construct_count": item.get("construct_count", ""),
                "cross_reactivity_count": item.get("cross_reactivity_count", ""),
                "has_alphafold_structure": int(bool(item.get("has_alphafold_structure"))),
                "structure_source": item.get("structure_source", ""),
                "has_pdb": int(bool(item.get("has_pdb"))),
                "has_any_structure": int(bool(item.get("has_any_structure"))),
                "has_interpro": int(bool(item.get("has_interpro"))),
                "has_family_context": int(bool(item.get("has_family_context"))),
                "has_assembly_requirements": int(bool(item.get("has_assembly_requirements"))),
                "has_assets": int(bool(item.get("has_assets"))),
                "pubtator_count": item.get("pubtator_count", ""),
                "detail_page": f"reports/{item.get('detail_page')}" if item.get("detail_page") else "",
                "duration_seconds": item.get("duration_seconds", ""),
            }
        )
    return [{field: row.get(field, "") for field in fields} for row in rows]


def _portal_download_fields(*, include_disease_context: bool = True) -> list[str]:
    fields = [
        "batch_index",
        "query",
        "source_column",
        "entry_name",
        "gene_symbol",
        "protein_name",
        "mouse_ortholog_accession",
        "mouse_ortholog_entry_name",
        "mouse_ortholog_gene_symbol",
        "mouse_refseq_accession",
        "mouse_sequence_source",
        "human_source_accession",
        "human_source_entry_name",
        "human_source_gene_symbol",
        "family_name",
        "status",
        "error",
        "topology_bucket",
        "ectodomain",
        "construct_count",
        "cross_reactivity_count",
        "has_alphafold_structure",
        "structure_source",
        "has_pdb",
        "has_any_structure",
        "has_interpro",
        "has_family_context",
        "has_assembly_requirements",
        "has_assets",
        "pubtator_count",
        "detail_page",
        "duration_seconds",
    ]
    if include_disease_context:
        insertion = fields.index("status")
        fields[insertion:insertion] = ["top_disease_name", "top_disease_score", "disease_count"]
    return fields


def portal_index_tsv(entries: list[dict[str, Any]], *, include_disease_context: bool = True) -> str:
    output = io.StringIO()
    writer = csv.DictWriter(
        output,
        fieldnames=_portal_download_fields(include_disease_context=include_disease_context),
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(_portal_download_rows(entries, include_disease_context=include_disease_context))
    return output.getvalue()


def portal_index_json(entries: list[dict[str, Any]], *, include_disease_context: bool = True) -> str:
    return json.dumps(_portal_download_rows(entries, include_disease_context=include_disease_context), indent=2) + "\n"


def portal_download_manifest(entries: list[dict[str, Any]], *, include_disease_context: bool = True, disease_downloads: tuple[str, ...] = ()) -> str:
    statuses: dict[str, int] = {}
    for item in entries:
        status = str(item.get("status") or "unknown")
        statuses[status] = statuses.get(status, 0) + 1
    manifest = {
        "name": "OpenAntigens Downloads",
        "generated_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "row_count": len(entries),
        "status_counts": statuses,
        "files": [
            {
                "path": "downloads/agdesign2_portal_index.tsv",
                "format": "TSV",
                "description": "Compact portal index with one row per batch target.",
                "columns": _portal_download_fields(include_disease_context=include_disease_context),
            },
            {
                "path": "downloads/agdesign2_portal_index.json",
                "format": "JSON",
                "description": "Same compact portal index represented as JSON objects.",
            },
        ],
        "citation": citation_metadata(),
        "source_notes": [
            "Files in this manifest are static release artifacts and can be downloaded directly.",
            "Third-party source data retain their own terms and citation requirements.",
        ],
    }
    manifest["files"].extend([
        {"path": "downloads/openantigens.bib", "format": "BibTeX", "description": "OpenAntigens preprint citation."},
        {"path": "downloads/openantigens.ris", "format": "RIS", "description": "OpenAntigens preprint citation."},
    ])
    if include_disease_context:
        manifest["files"].extend(
            {
                "path": f"open_targets_disease_associations.{suffix}",
                "format": suffix.upper(),
                "description": "Open Targets target-disease association scores.",
            }
            for suffix in ("tsv", "json") if suffix in disease_downloads
        )
    return json.dumps(manifest, indent=2) + "\n"


def portal_build_metadata(entries: list[dict[str, Any]], *, portal_title: str = "OpenAntigens") -> str:
    rows = _portal_download_rows(entries)
    status_counts = Counter(str(row.get("status") or "unknown") for row in rows)
    topology_counts = Counter(str(row.get("topology_bucket") or "unknown") for row in rows)
    metadata = {
        "database": portal_title,
        "citation": citation_metadata(),
        "software_version": _openantigen_version(),
        "portal_build_date": _portal_build_date(),
        "portal_build_started_utc": _PORTAL_PROCESS_STARTED_UTC.isoformat(),
        "generated_at_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "target_count": len(rows),
        "status_counts": dict(sorted(status_counts.items())),
        "topology_counts": dict(sorted(topology_counts.items())),
        "license": {
            "software": "Apache-2.0",
            "openantigen_generated_annotations": "CC BY 4.0",
            "third_party_source_data": "Retains original source database terms; see DATA_LICENSE.md.",
        },
    }
    return json.dumps(metadata, indent=2) + "\n"


def _write_portal_download_files(
    portal_dir: Path,
    entries: list[dict[str, Any]],
    batch_dir: Path | None = None,
    *,
    include_disease_context: bool = True,
) -> None:
    downloads_dir = portal_dir / "downloads"
    downloads_dir.mkdir(parents=True, exist_ok=True)
    _write_text_atomic(
        downloads_dir / "agdesign2_portal_index.tsv",
        portal_index_tsv(entries, include_disease_context=include_disease_context),
    )
    _write_text_atomic(
        downloads_dir / "agdesign2_portal_index.json",
        portal_index_json(entries, include_disease_context=include_disease_context),
    )
    for suffix in ("tsv", "json"):
        destination = portal_dir / f"open_targets_disease_associations.{suffix}"
        source = batch_dir / destination.name if batch_dir is not None else None
        if include_disease_context and source is not None and source.is_file():
            shutil.copy2(source, destination)
        else:
            destination.unlink(missing_ok=True)
    _write_text_atomic(downloads_dir / "download_manifest.json", portal_download_manifest(
        entries, include_disease_context=include_disease_context,
        disease_downloads=_available_disease_downloads(portal_dir),
    ))


def _available_disease_downloads(root: Path) -> tuple[str, ...]:
    return tuple(suffix for suffix in ("tsv", "json") if (root / f"open_targets_disease_associations.{suffix}").is_file())


def _load_open_targets_index(batch_dir: Path) -> dict[str, dict[str, Any]]:
    path = batch_dir / "open_targets_disease_associations.json"
    if not path.exists():
        return {}
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}
    index: dict[str, dict[str, Any]] = {}
    associations_by_entry: dict[str, list[dict[str, Any]]] = {}
    for row in payload.get("associations", []):
        key = str(row.get("entry_name") or "").upper()
        if key:
            associations_by_entry.setdefault(key, []).append(row)
    for target in payload.get("targets", []):
        keys = {
            str(target.get("entry_name") or "").upper(),
            str(target.get("gene_symbol") or "").upper(),
            str(target.get("uniprot_accession") or "").upper(),
        }
        disease_payload = dict(target)
        entry_key = str(target.get("entry_name") or "").upper()
        disease_payload["associations"] = associations_by_entry.get(entry_key, [])
        for key in keys:
            if key:
                index[key] = disease_payload
    return index


def _apply_open_targets_to_entry(entry: dict[str, Any], disease_index: dict[str, dict[str, Any]]) -> None:
    keys = (
        str(entry.get("entry_name") or "").upper(),
        str(entry.get("resolved_entry_name") or "").upper(),
        str(entry.get("gene_symbol") or "").upper(),
        str(entry.get("query") or "").split("_", 1)[0].upper(),
    )
    payload = next((disease_index[key] for key in keys if key and key in disease_index), None)
    if not payload:
        entry.setdefault("top_disease_name", "")
        entry.setdefault("top_disease_score", "")
        entry.setdefault("disease_count", "")
        entry.setdefault("open_targets_diseases", [])
        return
    associations = payload.get("associations") or []
    entry["top_disease_name"] = payload.get("top_disease_name") or ""
    entry["top_disease_score"] = payload.get("top_disease_score") or ""
    entry["top_disease_id"] = payload.get("top_disease_id") or ""
    entry["disease_count"] = payload.get("association_count_downloaded") or len(associations)
    entry["open_targets_url"] = payload.get("open_targets_url") or ""
    entry["open_targets_diseases"] = associations


def _entry_disease_haystack(item: dict[str, Any]) -> str:
    names = [str(item.get("top_disease_name") or "")]
    # Keep the static index compact. Per-target pages still render the top
    # Open Targets rows, but the main-page search only needs representative
    # disease names rather than thousands of ontology descendants per target.
    names.extend(str(row.get("disease_name") or "") for row in (item.get("open_targets_diseases") or [])[:100])
    return "|".join(sorted({name.strip().lower() for name in names if name.strip()}))


def _render_top_disease(item: dict[str, Any]) -> str:
    name = str(item.get("top_disease_name") or "").strip()
    if not name:
        return ""
    score = _fmt_number(item.get("top_disease_score"))
    count = item.get("disease_count")
    url = str(item.get("open_targets_url") or "").strip()
    label = escape(name)
    if url:
        label = f'<a href="{escape(url)}" target="_blank" rel="noreferrer">{label}</a>'
    pieces = [label]
    if score:
        pieces.append(f'<span class="muted-small">indirect score {escape(score)}</span>')
    if count:
        pieces.append(f'<span class="muted-small">{escape(str(count))} disease(s)</span>')
    return "<div class=\"disease-cell\">" + "".join(pieces) + "</div>"


def _render_open_targets_disease_rows(rows: list[dict[str, Any]]) -> str:
    if not rows:
        return "<tr><td colspan=\"6\">No Open Targets disease associations loaded.</td></tr>"
    rendered: list[str] = []
    for index, row in enumerate(rows[:100], start=1):
        disease_id = str(row.get("disease_id") or "")
        disease_name = str(row.get("disease_name") or "")
        url = str(row.get("open_targets_url") or "")
        disease = escape(disease_name)
        if url:
            disease = f'<a href="{escape(url)}" target="_blank" rel="noreferrer">{disease}</a>'
        direct_score = row.get("direct_score")
        indirect_score = row.get("indirect_score")
        if direct_score in (None, ""):
            direct_score = row.get("score")
        indirect_sources = str(row.get("indirect_datasource_scores") or "")
        direct_sources = str(row.get("datasource_scores") or "")
        sources = indirect_sources or direct_sources
        rendered.append(
            "<tr><td>{rank}</td><td>{disease}<br><span class=\"muted-small\">{disease_id}</span></td><td>{indirect_score}</td><td>{direct_score}</td><td>{sources}</td><td>{links}</td></tr>".format(
                rank=escape(str(row.get("rank") or index)),
                disease=disease,
                disease_id=escape(disease_id),
                indirect_score=escape(_fmt_number(indirect_score)),
                direct_score=escape(_fmt_number(direct_score)),
                sources=escape(sources),
                links=f'<a href="{escape(url)}" target="_blank" rel="noreferrer">Open Targets</a>' if url else "",
            )
        )
    return "".join(rendered)


def _render_obligatory_partner_warning(requirements: list[dict[str, Any]]) -> str:
    obligatory = [item for item in requirements if item.get("obligatory")]
    if not obligatory:
        return ""
    partners: list[str] = []
    summaries: list[str] = []
    for item in obligatory:
        partners.extend(str(partner) for partner in item.get("partners") or [] if str(partner).strip())
        summary = str(item.get("summary") or "").strip()
        if summary:
            summaries.append(summary)
    partner_text = ""
    unique_partners = sorted({partner.strip() for partner in partners if partner.strip()})
    paragraph_style = "max-width:100ch;margin:0 0 10px;color:#153f48;line-height:1.55;"
    if unique_partners:
        shown = ", ".join(unique_partners[:8])
        if len(unique_partners) > 8:
            shown += f", +{len(unique_partners) - 8} more"
        partner_text = f'<p style="{paragraph_style}"><strong>Annotated partner(s):</strong> {escape(shown)}</p>'
    evidence_text = ""
    if summaries:
        evidence_text = f'<p style="max-width:100ch;margin:0;color:#153f48;line-height:1.55;"><strong>Evidence summary:</strong> {escape(summaries[0])}</p>'
    return f"""
    <section class="partner-warning-card" id="obligatory-complex-context" data-report-section="Obligatory complex context" aria-label="Obligatory protein-complex context" style="margin:0 0 18px;padding:20px 22px;border-radius:22px;border:1px solid rgba(20,120,130,0.28);background:linear-gradient(135deg,rgba(231,249,247,0.96),rgba(240,249,255,0.92)),#eefaf8;box-shadow:0 18px 42px rgba(12,84,94,0.1);">
      <div>
        <p class="warning-kicker" style="margin:0 0 10px;color:#0f7a7f;font-size:0.78rem;font-weight:800;letter-spacing:0.14em;text-transform:uppercase;">Complex-aware design note</p>
        <h2 style="margin:0 0 8px;color:#0b4f5a;">Obligatory protein-complex context</h2>
        <p style="{paragraph_style}">Curated annotations indicate that the biologically relevant form of this target may be an obligatory protein complex rather than an isolated chain. Treat construct designs as starting points and review literature, structures, and assay goals to decide whether co-expression, partner-chain inclusion, or other complex-preserving strategies are needed.</p>
        {partner_text}
        {evidence_text}
      </div>
    </section>
    """


def _design_region_label(entry: dict[str, Any], report: dict[str, Any]) -> str:
    bucket = str(entry.get("topology_bucket") or entry.get("track") or "").strip().lower()
    if bucket:
        if "gpi" in bucket:
            return "GPI-Anchored Extracellular Region"
        if "secret" in bucket:
            return "Secreted Region"
        if "single" in bucket:
            return "Extracellular Region"
        if "multi" in bucket:
            ectodomain = report.get("ectodomain") or {}
            if ectodomain.get("start") is not None and ectodomain.get("end") is not None:
                return "Extracellular Region"
            return "Full-Length Membrane Protein"

    topology = report.get("topology")
    track = " ".join(
        str(value or "").lower()
        for value in (
            entry.get("primary_topology_class"),
            entry.get("topology_class"),
            report.get("topology_class"),
            topology.get("topology_class") if isinstance(topology, dict) else "",
        )
    )
    if "gpi" in track or "single" in track or "multi" in track or "membrane" in track or "surface" in track:
        return "Extracellular Region"
    if "secret" in track:
        return "Secreted Region"
    ectodomain = report.get("ectodomain") or {}
    source = str(ectodomain.get("source") or ectodomain.get("label") or "").lower()
    if "secret" in source:
        return "Secreted Region"
    if "transmembrane" in source or "topology" in source or "extracellular" in source:
        return "Extracellular Region"
    return "Design Region"


def _soluble_construct_label(entry: dict[str, Any], report: dict[str, Any]) -> str:
    bucket = str(entry.get("topology_bucket") or entry.get("track") or "").strip().lower()
    if "gpi" in bucket:
        return "Soluble GPI-Anchored Extracellular-region Constructs"
    if "secret" in bucket:
        return "Soluble Secreted-region Constructs"
    if "single" in bucket:
        return "Soluble Extracellular-region Constructs"
    label = _design_region_label(entry, report)
    if label == "Secreted Region":
        return "Soluble Secreted-region Constructs"
    if "Extracellular Region" in label:
        return "Soluble Extracellular-region Constructs"
    return "Soluble Design-region Constructs"


def _render_entry_track(entry: dict[str, Any]) -> str:
    track = str(entry.get("topology_bucket") or entry.get("track") or "").strip()
    if not track:
        return ""
    return f"<dt>Track</dt><dd>{escape(track)}</dd>"


def _render_target_aliases(target: dict[str, Any], entry: dict[str, Any] | None = None) -> str:
    aliases = _target_aliases(target, entry=entry)
    if not aliases:
        return ""
    alias_text = "; ".join(aliases[:20])
    if len(aliases) > 20:
        alias_text += f"; +{len(aliases) - 20} more"
    return f"<dt>Also known as</dt><dd>{escape(alias_text)}</dd>"


def _target_aliases(target: dict[str, Any], *, entry: dict[str, Any] | None = None) -> list[str]:
    names: list[str] = []
    primary_names = {
        str(target.get("protein_name") or "").strip().lower(),
        str(target.get("gene_symbol") or "").strip().lower(),
        str(target.get("entry_name") or "").strip().lower(),
    }

    def add(value: Any) -> None:
        cleaned = str(value or "").strip()
        if not cleaned:
            return
        key = cleaned.lower()
        if key in primary_names or any(key == item.lower() for item in names):
            return
        names.append(cleaned)

    for key in ("alternative_names", "aliases", "synonyms", "protein_synonyms"):
        raw = target.get(key)
        if isinstance(raw, list):
            for item in raw:
                add(item)
        elif isinstance(raw, str):
            for item in re.split(r"\s*[;|]\s*", raw):
                add(item)
    for item in target.get("gene_synonyms") or []:
        add(item)

    # Compatibility for existing reports: older JSON did not store explicit UniProt
    # alternative names, but many useful short names are embedded parenthetically in
    # the recommended protein label.
    protein_labels = [str(target.get("protein_name") or "")]
    if entry:
        protein_labels.append(str(entry.get("protein_name") or ""))
    for protein_name in protein_labels:
        for match in re.findall(r"\(([^()]+)\)", protein_name):
            for part in re.split(r"\s*;\s*", match):
                add(part)
    combined_labels = " ".join(protein_labels + names).lower()
    for cd_name in re.findall(r"\bCD antigen (CD\d+[A-Za-z0-9-]*)\b", " ".join(protein_labels + names)):
        add(cd_name)
        if "heavy chain" in combined_labels and not cd_name.lower().endswith("hc"):
            add(f"{cd_name}hc")
    return names


def _blast_search_unavailable(report: dict[str, Any], *, full_length: bool) -> bool:
    for note in report.get("notes") or []:
        if not isinstance(note, dict):
            continue
        if str(note.get("source") or "").strip().upper() != "BLAST":
            continue
        message = str(note.get("message") or "").strip().lower()
        if full_length and "full-length" not in message:
            continue
        if not full_length and "full-length" in message:
            continue
        if "failed" in message or "not run" in message or "disabled" in message:
            return True
    return False


def render_detail_page(
    entry: dict[str, Any],
    report: dict[str, Any],
    *,
    batch_dir: Path,
    portal_dir: Path | None = None,
    portal_title: str = "OpenAntigens",
    include_disease_context: bool = True,
    sibling_link: tuple[str, str] | None = None,
    allow_af3_structures: bool = False,
) -> str:
    page_dir = (portal_dir or batch_dir / "portal") / "reports"
    target = report.get("target", {})
    ectodomain = report.get("ectodomain") or {}
    design_region_label = _design_region_label(entry, report)
    design_region_lower = design_region_label.lower()
    soluble_construct_label = _soluble_construct_label(entry, report)
    target_aliases_html = _render_target_aliases(target, entry)
    constructs = _sorted_constructs_for_display(report.get("construct_details") or [])
    homology = report.get("ectodomain_homology") or []
    cross_hits = report.get("cross_reactivity_hits") or []
    full_length_cross_hits = report.get("full_length_cross_reactivity_hits") or []
    family_context = report.get("family_context")
    full_length_family_context = report.get("full_length_family_context")
    canonical_family = report.get("canonical_family")
    assembly_requirements = report.get("assembly_requirements") or []
    complex_portal_complexes = report.get("complex_portal_complexes") or []
    complex_portal_lookup = report.get("complex_portal_lookup") or {}
    advanced_membrane_suggestions = report.get("advanced_membrane_suggestions") or []
    gpcr_engineering_variants = report.get("gpcr_engineering_variants") or []
    gpcr_annotation = report.get("gpcr_annotation")
    experimental_constructs = report.get("experimental_constructs") or []
    interpro_annotations = report.get("interpro_annotations") or []
    ptms = _report_ptms(report)
    notes = report.get("notes") or []
    has_matching_alphafold = bool(entry.get("has_alphafold_structure"))
    report_json_path = _resolve_existing_path(entry.get("json_report"), batch_dir)
    report_md_path = _resolve_existing_path(entry.get("markdown_report"), batch_dir)
    structure_widget = _render_structure_widget(
        entry,
        report,
        batch_dir=batch_dir,
        page_dir=page_dir,
        allow_af3_structures=allow_af3_structures,
    )
    has_interactive_structure = 'id="structureViewer"' in structure_widget

    soluble_constructs, membrane_constructs = _partition_construct_dicts(constructs)

    homology_rows_list: list[str] = []
    for item in homology:
        homolog_target = _fetch_homolog_target(item, batch_dir=batch_dir) or {}
        homology_rows_list.append(
            "<tr><td>{species}</td><td>{entry}</td><td>{identity}</td><td>{surface_identity}</td><td>{coverage}</td><td>{links}</td></tr>".format(
                species=escape(_species_display_label(item.get("species"))),
                entry=escape(str(item.get("entry_name") or item.get("accession") or "")),
                identity=escape(_fmt_pct(item.get("identity"))),
                surface_identity=escape(_fmt_surface_identity(item)),
                coverage=escape(_fmt_pct(item.get("coverage"))),
                links=_render_resource_links(
                    accession=str(homolog_target.get("accession") or item.get("accession") or ""),
                    entry_name=str(homolog_target.get("entry_name") or item.get("entry_name") or ""),
                    gene_symbol=str(homolog_target.get("gene_symbol") or target.get("gene_symbol") or ""),
                    species=str(item.get("species") or ""),
                ),
            )
        )
    homology_rows = "".join(homology_rows_list) or "<tr><td colspan=\"6\">No homolog records.</td></tr>"

    target_sequence = str(target.get("sequence") or "")
    ecto_start = ectodomain.get("start")
    ecto_end = ectodomain.get("end")
    human_ecto_sequence = (
        slice_sequence(target_sequence, int(ecto_start), int(ecto_end))
        if target_sequence and ecto_start is not None and ecto_end is not None
        else ""
    )

    query_species_label = _target_species_key(target)
    cross_rows = (
        _render_cross_reactivity_rows(
            cross_hits,
            batch_dir=batch_dir,
            query_label=query_species_label,
        )
        if cross_hits
        else ""
    )
    full_length_cross_rows = _render_cross_reactivity_rows(
        full_length_cross_hits,
        batch_dir=batch_dir,
        query_label=query_species_label,
    )
    blast_unavailable = _blast_search_unavailable(report, full_length=False)
    if cross_hits:
        design_region_cross_html = f"""
      <h3>{escape(design_region_label)} BLAST sequence similarity hits</h3>
      <p class="section-note">Showing {len(cross_hits)} hit(s), ranked by bit score. Scroll within the table to inspect the full list and expand rows for BLAST alignments.</p>
      <div class="table-scroll cross-reactivity-scroll">
        <table>
          <thead><tr><th>Rank</th><th>Hit / Alignment</th><th>Species</th><th>Bit score</th><th>Identity</th><th>Extracellular accessible identity</th><th>Coverage</th><th>E-value</th><th>Links</th></tr></thead>
          <tbody>{cross_rows}</tbody>
        </table>
      </div>
        """
    elif blast_unavailable:
        design_region_cross_html = f"""
      <h3>{escape(design_region_label)} BLAST sequence similarity hits</h3>
      <p class="section-note">This BLAST sequence similarity search was not completed for the current report.</p>
        """
    else:
        design_region_cross_html = f"""
      <h3>{escape(design_region_label)} BLAST sequence similarity hits</h3>
      <p class="section-note">The completed BLAST search returned no sequence similarity hits.</p>
        """
    full_length_cross_html = ""
    if full_length_cross_hits:
        full_length_cross_html = f"""
        <h3>Full-length BLAST sequence similarity hits</h3>
        <p class="section-note">Showing {len(full_length_cross_hits)} full-length hit(s), ranked by bit score.</p>
        <div class="table-scroll cross-reactivity-scroll">
          <table>
            <thead><tr><th>Rank</th><th>Hit / Alignment</th><th>Species</th><th>Bit score</th><th>Identity</th><th>Extracellular accessible identity</th><th>Coverage</th><th>E-value</th><th>Links</th></tr></thead>
            <tbody>{full_length_cross_rows}</tbody>
          </table>
        </div>
        """

    cysteine_rows = "".join(
        "<tr><td>{position}</td><td>{paired}</td><td>{surface}</td><td>{nearest}</td><td>{warning}</td></tr>".format(
            position=escape(f"C{item.get('position')}") if item.get("position") is not None else "",
            paired=escape(f"C{item.get('paired_with')}" if item.get("paired_with") is not None else "unpaired"),
            surface=escape("yes" if item.get("surface_exposed") else "no" if item.get("surface_exposed") is not None else "n/a"),
            nearest=escape(_closest_unpaired_cysteine_text(item)),
            warning=escape(str(item.get("warning") or "")),
        )
        for item in (report.get("cysteine_analysis") or [])
    ) or "<tr><td colspan=\"5\">No cysteine analysis available.</td></tr>"

    furin_rows = "".join(
        "<tr><td>{boundary}</td><td><code>{motif}</code></td><td>{suggestions}</td></tr>".format(
            boundary=escape(f"{item.get('start')}-{item.get('end')}" if item.get("start") is not None and item.get("end") is not None else ""),
            motif=escape(str(item.get("motif") or "")),
            suggestions=escape(", ".join(item.get("suggestions") or [])),
        )
        for item in (report.get("furin_sites") or [])
    ) or "<tr><td colspan=\"3\">No candidate furin-like motifs found.</td></tr>"

    ptm_rows = _render_ptm_rows(ptms)

    obligatory_assembly_requirements = [
        item for item in assembly_requirements if item.get("obligatory")
    ]
    assembly_rows = "".join(
        "<tr><td>{classification}</td><td>{obligatory}</td><td>{confidence}</td><td>{summary}</td><td>{partners}</td></tr>".format(
            classification=escape(str(item.get("classification") or "")),
            obligatory=escape("yes" if item.get("obligatory") else "no"),
            confidence=escape(str(item.get("confidence") or "")),
            summary=escape(str(item.get("summary") or "")),
            partners=escape(", ".join(item.get("partners") or [])),
        )
        for item in obligatory_assembly_requirements
    )
    obligatory_partner_warning_html = _render_obligatory_partner_warning(assembly_requirements)
    assembly_requirements_html = (
        f"""
        <h3>Assembly requirements</h3>
        <table>
          <thead><tr><th>Classification</th><th>Obligatory</th><th>Confidence</th><th>Summary</th><th>Partners</th></tr></thead>
          <tbody>{assembly_rows}</tbody>
        </table>"""
        if assembly_rows
        else ""
    )
    complex_portal_html = (
        f"""
        <h3>Curated Complex Portal records</h3>
        {_render_complex_portal_context(complex_portal_complexes, complex_portal_lookup)}"""
        if portal_title != "OpenAntigens Mouse"
        else ""
    )
    interaction_assembly_html = (
        complex_portal_html
        + assembly_requirements_html
        or '<p class="section-note">No assembly requirements reported.</p>'
    )

    soluble_construct_summary_rows = _render_precomputed_construct_summary_rows(
        soluble_constructs,
        empty_message=f"No soluble {design_region_lower} constructs.",
        allow_structure_actions=has_interactive_structure,
    )
    membrane_construct_summary_rows = _render_precomputed_construct_summary_rows(
        membrane_constructs,
        empty_message="No membrane-expression constructs.",
        allow_structure_actions=has_interactive_structure,
    )
    membrane_construct_summary_html = (
        f"""
      <h3>Native Membrane-Expression Constructs</h3>
      <table>
        <thead><tr><th>Construct</th><th>Group</th><th>Boundary</th><th>Length</th><th>Support</th><th>Action</th></tr></thead>
        <tbody>{membrane_construct_summary_rows}</tbody>
      </table>"""
        if membrane_constructs
        else ""
    )
    construct_details_html = _render_construct_tabs(
        constructs,
        design_region_label=design_region_label,
        target=target,
        batch_dir=batch_dir,
        page_dir=page_dir,
        furin_sites=report.get("furin_sites") or [],
        ptms=ptms,
        allow_structure_assets=has_matching_alphafold,
        identity_label=_identity_label(report),
    )

    interpro_rows = "".join(
        "<tr><td>{start}</td><td>{end}</td><td>{source}</td><td>{acc}</td><td>{type}</td><td>{name}</td></tr>".format(
            start=escape(str(item.get("start") or "")),
            end=escape(str(item.get("end") or "")),
            source=escape(str(item.get("source_database") or "")),
            acc=escape(str(item.get("accession") or "")),
            type=escape(str(item.get("type") or "")),
            name=escape(str(item.get("name") or "")),
        )
        for item in interpro_annotations[:50]
    ) or "<tr><td colspan=\"6\">No InterPro annotations.</td></tr>"

    family_html = _render_family_context_section(
        family_context,
        canonical_family=canonical_family,
        matrix_scope_label=design_region_lower,
    )
    # No-ectodomain targets (e.g. multipass GPCRs like OPRM) are absent from the
    # ectodomain-surface comparison by design. Point readers to the full-length
    # matrix, which does include the target, so its absence above isn't mistaken
    # for missing data.
    if family_context and report.get("ectodomain") is None and full_length_family_context:
        family_html += (
            "<p class=\"section-note\">This target has no extracellular ectodomain, so it is not a row in the "
            "comparison above. See the Full-Length Family Matrix below for its full-length comparison to family members.</p>"
        )
    full_length_family_html = ""
    if full_length_family_context:
        full_length_family_html = (
            "<details class=\"details-panel\"><summary>Full-Length Family Matrix</summary>"
            + _render_family_context_section(
                full_length_family_context,
                canonical_family=None,
                matrix_scope_label="full-length sequence",
            )
            + "</details>"
        )
    advanced_membrane_html = _render_advanced_membrane_suggestions_section(advanced_membrane_suggestions)
    advanced_membrane_section_html = (
        f"""
    <section class="card" id="membrane-engineering" data-report-section="Membrane Engineering">
      <h2>Membrane Engineering Test Variants</h2>
      {advanced_membrane_html}
    </section>
"""
        if advanced_membrane_suggestions
        else ""
    )
    gpcr_engineering_html = _render_gpcr_engineering_variants_section(gpcr_engineering_variants)
    gpcr_annotation_html = _render_gpcr_annotation_section(gpcr_annotation)

    note_items = "".join(
        f"<li><strong>{escape(str(note.get('severity') or '').upper())}</strong>: {escape(str(note.get('message') or ''))}</li>"
        for note in notes
    ) or "<li>No notes.</li>"

    json_link = _relative_link(report_json_path, page_dir) if report_json_path else None
    md_link = _relative_link(report_md_path, page_dir) if report_md_path else None
    construct_mutation_script = (
        ""
        if has_interactive_structure
        else f"<script>{_construct_mutation_card_script(target=target)}</script>"
    )
    construct_summary_note = (
        "Use the <strong>View in 3D structure</strong> button on any construct to highlight it in the Interactive Construct Builder above. "
        if has_interactive_structure
        else ""
    )
    target_species = _target_species_key(target)
    report_metadata = report.get("metadata") if isinstance(report.get("metadata"), dict) else {}
    human_source_target = report_metadata.get("human_source_target") if isinstance(report_metadata, dict) else {}
    human_source_html = ""
    if target_species == "mouse" and isinstance(human_source_target, dict):
        human_source_entry = str(human_source_target.get("entry_name") or "").strip()
        human_source_accession = str(human_source_target.get("accession") or "").strip()
        if human_source_entry or human_source_accession:
            human_source_html = (
                f"<dt>Originating human target</dt><dd>{escape(human_source_entry or human_source_accession)}"
                + (f" ({escape(human_source_accession)})" if human_source_entry and human_source_accession else "")
                + "</dd>"
            )
    target_links_html = _render_resource_links(
        accession=str(target.get("accession") or ""),
        entry_name=str(target.get("entry_name") or ""),
        gene_symbol=str(target.get("gene_symbol") or ""),
        species=target_species,
    )
    pubtator_html = _render_pubtator_link(entry.get("pubtator_count"), entry.get("pubtator_query_url"))
    disease_rows = _render_open_targets_disease_rows(entry.get("open_targets_diseases") or [])
    disease_section = (
        f"""
    <section class="card" id="disease-associations" data-report-section="Disease Associations">
      <h2>Open Targets Disease Associations</h2>
      <p class="section-note">Disease associations from Open Targets. The portal index is ranked and searched using the broader indirect score; direct overall scores are shown here for comparison.</p>
      <div class="table-scroll cross-reactivity-scroll">
        <table>
          <thead><tr><th>Rank</th><th>Disease</th><th>Indirect score</th><th>Direct score</th><th>Indirect datasource scores</th><th>Links</th></tr></thead>
          <tbody>{disease_rows}</tbody>
        </table>
      </div>
    </section>
"""
        if include_disease_context
        else ""
    )

    return f"""<!DOCTYPE html>
<html lang="en" class="report-page">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{escape(str(target.get('entry_name') or entry.get('query') or 'OpenAntigens report'))}</title>
{_portal_favicon_links(prefix="../")}
  <link rel="stylesheet" href="../portal.css?v={_PORTAL_CSS_VERSION}">
{_portal_analytics_head()}
</head>
<body>
  <main class="page detail">
    <header class="site-hero site-hero-detail">
      {_portal_brand_header(prefix="../", portal_title=portal_title, sibling_link=sibling_link, chip_html='<a class="hero-chip hero-chip-link" href="../index.html">&larr; Back to portal index</a>')}
      <div class="hero-grid">
        <div class="hero-copy">
          <p class="eyebrow">Target report</p>
          <h1>{escape(str(target.get('entry_name') or 'Unknown target'))}</h1>
          <p class="hero-text">{escape(str(target.get('protein_name') or ''))}</p>
        </div>
      </div>
    </header>

    <section class="grid two" id="overview" data-report-section="Overview">
      <article class="card">
        <h2>Target</h2>
        <dl class="meta">
          <dt>Gene</dt><dd>{escape(str(target.get('gene_symbol') or ''))}</dd>
          <dt>Accession</dt><dd>{escape(str(target.get('accession') or ''))}</dd>
          <dt>Organism</dt><dd>{escape(str(target.get('organism') or ''))}</dd>
          <dt>Sequence length</dt><dd>{len(target.get('sequence') or '')} aa</dd>
          {human_source_html}
          <dt>Sequence analyzed</dt><dd>UniProt canonical{(' (' + escape(str(target['canonical_isoform_id'])) + ')') if target.get('canonical_isoform_id') else ''}</dd>
          {_render_entry_track(entry)}
          {target_aliases_html}
          <dt>Status</dt><dd>{escape(str(entry.get('status') or ''))}</dd>
          <dt>Links</dt><dd>{target_links_html}</dd>
          <dt>PubTator hits</dt><dd>{pubtator_html}</dd>
        </dl>
        <p class="section-note">Alternative isoforms are not analyzed separately; the coordinates and annotations shown refer to this sequence.</p>
      </article>
      <article class="card">
        <h2>{escape(design_region_label)}</h2>
        <dl class="meta">
          <dt>Boundary</dt><dd>{escape(_format_region(ectodomain))}</dd>
          <dt>Label</dt><dd>{escape(str(ectodomain.get('label') or ''))}</dd>
          <dt>Source</dt><dd>{escape(str(ectodomain.get('source') or ''))}</dd>
          <dt>Constructs</dt><dd>{len(constructs)}</dd>
          <dt>BLAST sequence similarity hits</dt><dd>{len(cross_hits)}</dd>
        </dl>
      </article>
    </section>

    {obligatory_partner_warning_html}

    <section class="card" id="interactive-construct-builder" data-report-section="Construct Builder" tabindex="-1">
      <h2>Construct Builder</h2>
      <p class="section-note construct-choice-help"><a class="inline-link" href="../constructs.html#choose-constructs" target="_blank" rel="noopener">Don&rsquo;t know how to pick your construct?</a> <span class="muted">(opens in new tab)</span></p>
      {structure_widget}
    </section>

    <section class="card" id="construct-summary" data-report-section="Construct Summary">
      <h2>Precomputed Construct Summary</h2>
      <p class="section-note">These are algorithm-generated starting points. {construct_summary_note}Detailed construct evidence, warnings, ortholog sequences, copy-ready exports, and images are provided in the <a class="inline-link" href="#construct-details">Construct Details</a> section at the end of this page.</p>
      <h3>{escape(soluble_construct_label)}</h3>
      <table>
        <thead><tr><th>Construct</th><th>Group</th><th>Boundary</th><th>Length</th><th>Support</th><th>Action</th></tr></thead>
        <tbody>{soluble_construct_summary_rows}</tbody>
      </table>
      {membrane_construct_summary_html}
    </section>

    {gpcr_engineering_html}

    {disease_section}

    <section class="card" id="homology" data-report-section="Homology">
      <h2>Homology</h2>
      <table>
        <thead><tr><th>Species</th><th>Entry</th><th>Identity</th><th>Extracellular accessible identity</th><th>Coverage</th><th>Links</th></tr></thead>
        <tbody>{homology_rows}</tbody>
      </table>
    </section>

    <section class="card" id="cross-reactivity" data-report-section="Cross-reactivity">
      <h2>Sequence similarity and cross-reactivity context</h2>
      {design_region_cross_html}
      {full_length_cross_html}
    </section>

    <section class="stack">
      <article class="card" id="cysteines" data-report-section="Cysteines">
        <h2>Cysteines</h2>
        <p class="section-note">Exposure is calculated from the target-chain model and selected structural blocks, without partner chains. A residue exposed here may be buried against another subunit, including an identical copy. Review <a class="inline-link" href="#interactions-assembly">Interactions / Assembly</a> when interpreting exposure.</p>
        <table>
          <thead><tr><th>Residue</th><th>Paired with</th><th>Surface exposed</th><th>Closest unpaired cysteine</th><th>Warning</th></tr></thead>
          <tbody>{cysteine_rows}</tbody>
        </table>
      </article>
      <article class="card" id="furin-cleavage-sites" data-report-section="Furin Cleavage Sites">
        <h2>Candidate Furin-like Motifs</h2>
        <table>
          <thead><tr><th>Boundary</th><th>Motif</th><th>Suggested edits</th></tr></thead>
          <tbody>{furin_rows}</tbody>
        </table>
      </article>
    </section>

    <section class="card" id="post-translational-modifications" data-report-section="Post-translational Modifications">
      <h2>Post-translational Modifications</h2>
      <p class="section-note">Curated UniProt PTM and molecule-processing annotations. Processing and cleavage annotations should be reviewed before finalizing construct boundaries.</p>
      <div class="table-scroll">
        <table>
          <thead><tr><th>Start</th><th>End</th><th>Type</th><th>Category</th><th>Description</th><th>Warning</th></tr></thead>
          <tbody>{ptm_rows}</tbody>
        </table>
      </div>
    </section>

    <section class="stack">
      <article class="card" id="family-context" data-report-section="Family Context">
        <h2>Family Context</h2>
        {family_html}
        {full_length_family_html}
      </article>
      <article class="card" id="interactions-assembly" data-report-section="Interactions / Assembly">
        <h2>Interactions / Assembly</h2>
        {interaction_assembly_html}
      </article>
    </section>

    {advanced_membrane_section_html}

    {gpcr_annotation_html}

    <section class="card" id="interpro-pfam" data-report-section="InterPro / Pfam">
      <h2>InterPro / Pfam</h2>
      <table>
        <thead><tr><th>Start</th><th>End</th><th>Source</th><th>Accession</th><th>Type</th><th>Name</th></tr></thead>
        <tbody>{interpro_rows}</tbody>
      </table>
    </section>

    <section class="card" id="construct-details" data-report-section="Construct Details">
      <h2>Construct Details</h2>
      {construct_details_html}
    </section>

    <section class="card" id="notes" data-report-section="Notes">
      <h2>Notes</h2>
      <ul>{note_items}</ul>
    </section>
    {_portal_footer(prefix="../")}
  </main>
  {construct_mutation_script}
  {_report_section_navigation()}
</body>
</html>
"""


def _render_complex_portal_context(complexes: list[dict[str, Any]], lookup: Any) -> str:
    human_lookup = lookup.get("human") if isinstance(lookup, dict) else None
    human_status = str((human_lookup or {}).get("status") or "not_recorded")
    status_html = {
        "no_hits": "<p class=\"section-note\">No curated human Complex Portal records matched this target.</p>",
        "error": "<p class=\"section-note\">The human Complex Portal lookup failed for this release.</p>",
        "disabled": "<p class=\"section-note\">Complex Portal lookup was disabled for this release.</p>",
        "not_recorded": "<p class=\"section-note\">Complex Portal lookup was not recorded for this report.</p>",
    }.get(human_status, "")
    if not complexes:
        return status_html or "<p class=\"section-note\">No curated Complex Portal records were returned.</p>"
    rows = "".join(
        "<tr><td><code>{complex_ac}</code></td><td>{name}</td><td>{species}</td><td>{assembly}</td><td>{participants}</td></tr>".format(
            complex_ac=escape(str(item.get("complex_ac") or "")),
            name=escape(str(item.get("name") or "")),
            species=escape(str(item.get("species") or "")),
            assembly=escape(", ".join(item.get("complex_assemblies") or []) or "n/a"),
            participants=escape(
                ", ".join(
                    str(participant.get("name") or participant.get("identifier") or "")
                    for participant in (item.get("participants") or [])
                    if isinstance(participant, dict)
                )
            ),
        )
        for item in complexes
        if isinstance(item, dict)
    )
    if not rows:
        rows = '<tr><td colspan="5">No readable Complex Portal records.</td></tr>'
    return f"""{status_html}
      <div class="table-scroll">
        <table>
          <thead><tr><th>Complex Portal ID</th><th>Complex</th><th>Species</th><th>Assembly</th><th>Participants</th></tr></thead>
          <tbody>{rows}</tbody>
        </table>
      </div>"""


def _report_section_navigation() -> str:
    return """
    <button type="button" class="report-sections-button" popovertarget="report-sections-menu">Sections</button>
    <div id="report-sections-menu" popover>
      <nav aria-label="Report sections"></nav>
    </div>
    <script type="module">
    const menu = document.getElementById('report-sections-menu');
    const sections = Array.from(document.querySelectorAll('[data-report-section]'));
    function focusHeading(section) {
      requestAnimationFrame(() => {
        const heading = section.querySelector('h2');
        heading.tabIndex = -1;
        heading.focus({ preventScroll: true });
      });
    }
    function focusHashHeading() {
      const section = sections.find(section => '#' + section.id === window.location.hash);
      if (section) focusHeading(section);
    }
    const links = sections.map((section) => {
      const link = document.createElement('a');
      link.href = '#' + section.id;
      link.textContent = section.dataset.reportSection;
      link.addEventListener('click', (event) => {
        if (event.button !== 0 || event.metaKey || event.ctrlKey || event.shiftKey || event.altKey) return;
        menu.hidePopover();
        focusHeading(section);
      });
      menu.querySelector('nav').append(link);
      return link;
    });
    let scheduled = false;
    function markCurrentSection() {
      scheduled = false;
      let current = 0;
      sections.forEach((section, index) => {
        if (section.getBoundingClientRect().top <= 32) current = index;
      });
      if (window.scrollY + window.innerHeight >= document.documentElement.scrollHeight - 2) current = sections.length - 1;
      links.forEach((link, index) => {
        if (index === current) link.setAttribute('aria-current', 'location');
        else link.removeAttribute('aria-current');
      });
    }
    function scheduleCurrentSection() {
      if (!scheduled) {
        scheduled = true;
        requestAnimationFrame(markCurrentSection);
      }
    }
    window.addEventListener('scroll', scheduleCurrentSection, { passive: true });
    window.addEventListener('resize', scheduleCurrentSection);
    window.addEventListener('load', scheduleCurrentSection, { once: true });
    window.addEventListener('load', focusHashHeading, { once: true });
    window.addEventListener('hashchange', focusHashHeading);
    menu.addEventListener('toggle', scheduleCurrentSection);
    scheduleCurrentSection();
    </script>
    """


def _render_construct_tabs(
    constructs: list[dict[str, Any]],
    *,
    design_region_label: str,
    target: dict[str, Any],
    batch_dir: Path,
    page_dir: Path,
    furin_sites: list[dict[str, Any]],
    ptms: list[dict[str, Any]],
    allow_structure_assets: bool,
    identity_label: str,
) -> str:
    if not constructs:
        return "<p>No constructs available.</p>"
    full_label = "Full design region"
    if design_region_label == "Secreted Region":
        full_label = "Full secreted region"
    elif "Extracellular Region" in design_region_label:
        full_label = "Full ectodomain"
    labels = [full_label, "PDB", "Annotated domains / repeats", "Strict", "Lenient", "Full-length multipass", "Trimmed multipass"]
    groups: list[list[dict[str, Any]]] = [[] for _ in labels]
    for construct in constructs:
        group = _construct_display_group(construct)[0]
        if group == 5 and construct.get("name") != "membrane_expression_full_length":
            group = 6
        groups[group].append(construct)
    tabs: list[str] = []
    panels: list[str] = []
    for index, (label, group) in enumerate(zip(labels, groups)):
        if not group:
            continue
        selected = not tabs
        tab_id = f"construct-tab-{index}"
        panel_id = f"construct-panel-{index}"
        tabs.append(
            f'<button type="button" role="tab" id="{tab_id}" aria-controls="{panel_id}" '
            f'aria-selected="{str(selected).lower()}" tabindex="{0 if selected else -1}">{escape(label)} ({len(group)})</button>'
        )
        sections = [("", group)]
        if index == 1:
            soluble, membrane = _partition_construct_dicts(group)
            sections = [("Soluble PDB constructs", soluble), ("Membrane PDB constructs", membrane)]
        content: list[str] = []
        if index == 1:
            content.append('<p class="section-note">These sequences follow mapped PDB boundaries. Check the original study for tags, fusion partners, mutations, ligands and other components of the experimental construct.</p>')
        for heading, members in sections:
            if not members:
                continue
            if heading:
                content.append(f"<h3>{heading}</h3>")
            cards = "".join(
                _render_construct_card(
                    construct, target=target, batch_dir=batch_dir, page_dir=page_dir,
                    furin_sites=furin_sites, ptms=ptms, allow_structure_assets=allow_structure_assets,
                    identity_label=identity_label,
                )
                for construct in members
            )
            content.append(f'<div class="construct-list">{cards}</div>')
        panels.append(
            f'<div role="tabpanel" id="{panel_id}" aria-labelledby="{tab_id}" tabindex="0"'
            f'{"" if selected else " hidden"}>{"".join(content)}</div>'
        )
    return (
        '<div class="construct-tabs" role="tablist" aria-label="Construct categories">'
        + "".join(tabs) + "</div>" + "".join(panels)
        + """
        <script type="module">
        (() => {
          const section = document.getElementById('construct-details');
          const tabs = Array.from(section.querySelectorAll('[role="tab"]'));
          const activate = (active) => tabs.forEach((tab) => {
            const selected = tab === active;
            tab.setAttribute('aria-selected', String(selected));
            tab.tabIndex = selected ? 0 : -1;
            document.getElementById(tab.getAttribute('aria-controls')).hidden = !selected;
          });
          tabs.forEach((tab, index) => {
            tab.addEventListener('click', () => activate(tab));
            tab.addEventListener('keydown', (event) => {
              let next;
              if (event.key === 'ArrowRight') next = (index + 1) % tabs.length;
              else if (event.key === 'ArrowLeft') next = (index + tabs.length - 1) % tabs.length;
              else if (event.key === 'Home') next = 0;
              else if (event.key === 'End') next = tabs.length - 1;
              else return;
              event.preventDefault();
              activate(tabs[next]);
              tabs[next].focus();
            });
          });
        })();
        </script>
        """
    )


def _render_construct_card(
    construct: dict[str, Any],
    *,
    target: dict[str, Any],
    batch_dir: Path,
    page_dir: Path,
    furin_sites: list[dict[str, Any]],
    ptms: list[dict[str, Any]],
    allow_structure_assets: bool = True,
    identity_label: str = "Identity to human",
) -> str:
    human_entry = construct.get("human_entry_name") or target.get("entry_name") or ""
    evidence = construct.get("evidence") or []
    warnings = construct.get("warnings") or []
    homologs = construct.get("homologs") or []
    structural_metrics = construct.get("structural_metrics") or {}
    safe_name = _safe_filename(str(construct.get("name") or "construct"))
    assets_dir = batch_dir / f"{str(target.get('entry_name') or '').lower()}_report_assets"
    portal_assets_dir = page_dir.parent / "report_assets" / assets_dir.name
    use_portal_assets = portal_assets_dir.exists()
    if use_portal_assets:
        assets_dir = portal_assets_dir
    structure_reference = construct.get("structure_image")
    quality_reference = construct.get("quality_plot")
    if use_portal_assets:
        if structure_reference:
            structure_reference = str(assets_dir / Path(str(structure_reference)).name)
        if quality_reference:
            quality_reference = str(assets_dir / Path(str(quality_reference)).name)
    structure_path = None
    quality_path = None
    if allow_structure_assets:
        structure_path = _resolve_asset_for_construct(
            structure_reference,
            assets_dir / f"{safe_name}_structure.png",
            base_dir=batch_dir,
        )
        quality_path = _resolve_asset_for_construct(
            quality_reference,
            assets_dir / f"{safe_name}_quality.png",
            base_dir=batch_dir,
        )
    structure_img = (
        f'<img loading="lazy" src="{escape(_relative_link(structure_path, page_dir))}" alt="{escape(safe_name)} structure">'
        if structure_path
        else ""
    )
    quality_img = (
        f'<img loading="lazy" src="{escape(_relative_link(quality_path, page_dir))}" alt="{escape(safe_name)} quality plot">'
        if quality_path
        else ""
    )

    evidence_list = "".join(f"<li>{escape(_scientific_region_text(item))}</li>" for item in evidence) or "<li>None</li>"
    warning_list = "".join(f"<li>{escape(_scientific_region_text(item))}</li>" for item in warnings) or "<li>None</li>"
    construct_ptms = _features_in_region(ptms, construct.get("start"), construct.get("end")) if ptms else construct.get("ptms") or []
    ptm_list = "".join(
        "<li>{label} <strong>{category}</strong>: {description}</li>".format(
            label=escape(_ptm_display_label(item)),
            category=escape(_ptm_category(item)),
            description=escape(str(item.get("description") or item.get("type") or "")),
        )
        for item in construct_ptms
    ) or "<li>None</li>"
    gene_symbol = str(target.get("gene_symbol") or target.get("entry_name") or "TARGET")
    primary_species = _target_species_key(target)
    primary_entry = str(target.get("entry_name") or human_entry)
    primary_identity = str(construct.get("identity_to_source_human") or "100.00%")
    if isinstance(construct.get("identity_to_source_human"), (int, float)):
        primary_identity = _fmt_pct(construct.get("identity_to_source_human"))
    export_rows = [
        {
            "name": _format_export_name(gene_symbol, primary_species, construct.get("start"), construct.get("end")),
            "species": _species_code(primary_species),
            "species_raw": primary_species,
            "boundary": f"{construct.get('start')}-{construct.get('end')}",
            "start": construct.get("start"),
            "end": construct.get("end"),
            "sequence": str(construct.get("sequence") or ""),
            "target_region_start": construct.get("start"),
            "homolog_region_start": construct.get("start"),
            "aligned_target": str(construct.get("sequence") or ""),
            "aligned_homolog": str(construct.get("sequence") or ""),
            "identity": primary_identity,
            "surface_identity": "100.00%",
            "accession": str(target.get("accession") or ""),
            "entry_name": primary_entry,
            "gene_symbol": gene_symbol,
            "links": _render_resource_links(
                accession=str(target.get("accession") or ""),
                entry_name=primary_entry,
                gene_symbol=gene_symbol,
                species=primary_species,
            ),
        }
    ]
    for item in homologs:
        if not (item.get("available") and item.get("start") is not None and item.get("end") is not None):
            continue
        species = str(item.get("species") or "")
        homolog_sequence = str(item.get("sequence") or "")
        alignment = _portal_align_sequences(
            str(construct.get("sequence") or ""),
            homolog_sequence,
        )
        export_rows.append(
            {
                "name": _format_export_name(gene_symbol, species, item.get("start"), item.get("end")),
                "species": _species_code(species),
                "species_raw": species,
                "boundary": f"{item.get('start')}-{item.get('end')}",
                "start": item.get("start"),
                "end": item.get("end"),
                "sequence": homolog_sequence,
                "target_region_start": construct.get("start"),
                "homolog_region_start": item.get("start"),
                "aligned_target": alignment["aligned_query"],
                "aligned_homolog": alignment["aligned_subject"],
                "identity": _fmt_pct(item.get("identity_to_human")),
                "surface_identity": _fmt_surface_identity(item),
                "accession": str(item.get("accession") or ""),
                "entry_name": str(item.get("entry_name") or item.get("accession") or ""),
                "gene_symbol": gene_symbol,
                "links": _render_resource_links(
                    accession=str(item.get("accession") or ""),
                    entry_name=str(item.get("entry_name") or item.get("accession") or ""),
                    gene_symbol=gene_symbol,
                    species=species,
                ),
            }
        )
    unavailable_rows = [
        {
            "name": f"{gene_symbol}_{_species_code(str(item.get('species') or ''))}_mapping",
            "species": _species_code(str(item.get("species") or "")),
            "species_raw": str(item.get("species") or ""),
            "reason": "; ".join(str(note) for note in (item.get("notes") or []) if str(note).strip())
            or "No usable homolog boundary was available.",
            "accession": str(item.get("accession") or ""),
            "entry_name": str(item.get("entry_name") or item.get("accession") or ""),
            "gene_symbol": gene_symbol,
            "links": _render_resource_links(
                accession=str(item.get("accession") or ""),
                entry_name=str(item.get("entry_name") or item.get("accession") or ""),
                gene_symbol=gene_symbol,
                species=str(item.get("species") or ""),
            ),
        }
        for item in homologs
        if not (item.get("available") and item.get("start") is not None and item.get("end") is not None)
    ]
    sequence_rows_html = "".join(
        "<tr><td><code>{name}</code></td><td>{species}</td><td>{boundary}</td><td><code class=\"sequence\">{sequence}</code></td><td>{identity}</td><td>{surface_identity}</td><td>{links}</td></tr>".format(
            name=escape(row["name"]),
            species=escape(row["species"]),
            boundary=escape(row["boundary"]),
            sequence=escape(row["sequence"]),
            identity=escape(row["identity"]),
            surface_identity=escape(row["surface_identity"]),
            links=row["links"],
        )
        for row in export_rows
    )
    sequence_rows_html += "".join(
        "<tr><td><code>{name}</code></td><td>{species}</td><td>unavailable</td>"
        '<td><span class="muted">Mapping unavailable: {reason}</span></td><td>n/a</td><td>n/a</td><td>{links}</td></tr>'.format(
            name=escape(row["name"]),
            species=escape(row["species"]),
            reason=escape(row["reason"]),
            links=row["links"],
        )
        for row in unavailable_rows
    )
    construct_tsv = _export_tsv(
        export_rows,
        identity_header=_identity_header(identity_label),
        surface_identity_header=_surface_identity_header(identity_label),
    )
    construct_fasta = _export_fasta(
        export_rows,
        identity_header=_identity_header(identity_label),
        surface_identity_header=_surface_identity_header(identity_label),
    )
    construct_mutation_candidates = sorted(
        {
            int(finding.get("position"))
            for finding in (construct.get("cysteine_analysis") or [])
            if finding.get("position") is not None and finding.get("paired_with") is None
            and not str(finding.get("warning") or "").startswith("Ambiguous")
        }
    )
    construct_start = construct.get("start")
    construct_end = construct.get("end")
    try:
        construct_start_int = int(construct_start)
        construct_end_int = int(construct_end)
    except (TypeError, ValueError):
        construct_start_int = None
        construct_end_int = None
    construct_furin_sites = []
    if construct_start_int is not None and construct_end_int is not None:
        construct_furin_sites = [
            site
            for site in furin_sites
            if site.get("start") is not None
            and site.get("end") is not None
            and construct_start_int <= int(site.get("start"))
            and int(site.get("end")) <= construct_end_int
        ]
    construct_payload = {
        "geneSymbol": gene_symbol,
        "identityHeader": _identity_header(identity_label),
        "surfaceIdentityHeader": _surface_identity_header(identity_label),
        "rows": [
            {
                "species": row["species"],
                "speciesRaw": row["species_raw"],
                "start": row["start"],
                "end": row["end"],
                "sequence": row["sequence"],
                "targetRegionStart": row["target_region_start"],
                "homologRegionStart": row["homolog_region_start"],
                "alignedTarget": row["aligned_target"],
                "alignedHomolog": row["aligned_homolog"],
                "identity": row["identity"],
                "surfaceIdentity": row["surface_identity"],
                "accession": row["accession"],
                "entryName": row["entry_name"],
                "geneSymbol": row["gene_symbol"],
                "links": row["links"],
            }
            for row in export_rows
        ],
        "unavailableRows": [
            {
                "name": row["name"],
                "species": row["species"],
                "speciesRaw": row["species_raw"],
                "reason": row["reason"],
                "accession": row["accession"],
                "entryName": row["entry_name"],
                "geneSymbol": row["gene_symbol"],
                "links": row["links"],
            }
            for row in unavailable_rows
        ],
        "mutationCandidates": construct_mutation_candidates,
        "furinSites": construct_furin_sites,
        "ptms": construct_ptms,
    }

    metrics = [] if not allow_structure_assets else [
        ("Mean pLDDT", _fmt_number(structural_metrics.get("mean_plddt"))),
        ("Structured coverage", _fmt_fraction_pct(structural_metrics.get("structured_fraction"))),
        ("Mean intra-PAE", _fmt_number(structural_metrics.get("mean_intra_pae"))),
    ]
    metrics_html = "".join(f"<li><strong>{label}:</strong> {escape(value)}</li>" for label, value in metrics if value)
    gpcr_segments = construct.get("gpcr_segments") or []
    gpcr_segment_list = "".join(
        "<li>{name}: {start}-{end}{generic}</li>".format(
            name=escape(str(segment.get("name") or "")),
            start=escape(str(segment.get("start") or "")),
            end=escape(str(segment.get("end") or "")),
            generic=(
                " ("
                + escape(str(segment.get("generic_start") or ""))
                + "-"
                + escape(str(segment.get("generic_end") or ""))
                + ")"
                if segment.get("generic_start") or segment.get("generic_end")
                else ""
            ),
        )
        for segment in gpcr_segments
    ) or "<li>None</li>"
    gpcr_range = construct.get("gpcr_generic_range")
    gpcr_construct_html = ""
    if gpcr_segments or gpcr_range:
        gpcr_construct_html = f"""
          <div>
            <h4>GPCR Segments</h4>
            <p><strong>Generic range:</strong> {escape(str(gpcr_range or 'n/a'))}</p>
            <ul>{gpcr_segment_list}</ul>
          </div>
        """
    viewer_button = (
        _viewer_focus_button(construct.get("start"), construct.get("end"), construct.get("name"))
        if allow_structure_assets
        else ""
    )
    viewer_action_row = (
        f'<div class="construct-actions">{viewer_button}</div>' if viewer_button else ""
    )

    return f"""
    <article class="construct-card">
      <div class="construct-copy">
        <h3>{_construct_name_html(construct)}</h3>
        {viewer_action_row}
        <p><strong>Length:</strong> {escape(str(construct.get("length") or ""))} aa</p>
        <p><strong>Boundary:</strong> {escape(str(construct.get("start") or ""))}-{escape(str(construct.get("end") or ""))}</p>
        <p><strong>Classification:</strong> <code>{escape(str(construct.get("classification") or ""))}</code></p>
        <p><strong>Rationale:</strong> {escape(_scientific_region_text(construct.get("rationale") or ""))}</p>
        <p><strong>Domain summary:</strong> {escape(str(construct.get("domain_summary") or "None"))}</p>
        <div class="construct-export-card">
          <h4>Construct Sequences</h4>
          <div class="construct-mutation-controls"></div>
          <div class="selection-table-wrap">
            <table class="construct-sequence-table">
              <thead><tr><th>Name</th><th>Species</th><th>Boundary</th><th>Sequence</th><th>{escape(identity_label)}</th><th>Extracellular accessible identity</th><th>Links</th></tr></thead>
              <tbody class="construct-sequence-table-body">{sequence_rows_html}</tbody>
            </table>
          </div>
          <div class="export-box-grid">
            <div>
              <label class="export-label">TSV</label>
              <textarea class="selection-tsv construct-export-box construct-export-tsv" readonly>{escape(construct_tsv)}</textarea>
            </div>
            <div>
              <label class="export-label">FASTA</label>
              <textarea class="selection-tsv construct-export-box construct-export-fasta" readonly>{escape(construct_fasta)}</textarea>
            </div>
          </div>
          <script type="application/json" class="construct-payload">{_js_json(construct_payload)}</script>
        </div>
        <div class="construct-grid">
          <div>
            <h4>Evidence</h4>
            <ul>{evidence_list}</ul>
          </div>
          <div>
            <h4>Warnings</h4>
            <ul>{warning_list}</ul>
          </div>
          <div>
            <h4>Included PTMs</h4>
            <ul>{ptm_list}</ul>
          </div>
          <div>
            <h4>Structural Metrics</h4>
            <ul>{metrics_html or '<li>None</li>'}</ul>
          </div>
          {gpcr_construct_html}
        </div>
      </div>
      <div class="construct-images">
        {structure_img}
        {quality_img}
      </div>
    </article>
    """


def _construct_mutation_card_script(*, target: dict[str, Any]) -> str:
    return (
        f"const targetSpecies = {_js_json(_target_species_key(target))};\n"
        + """
      (() => {
        const homologPositionLookupCache = new Map();

        function speciesCode(species) {
          const normalized = String(species || '').toLowerCase();
          if (normalized === 'human') return 'HUMAN';
          if (normalized === 'mouse') return 'MOUSE';
          if (normalized === 'macaca_fascicularis') return 'MACFA';
          return String(species || 'SPECIES').toUpperCase().replace(/[^A-Z0-9]+/g, '_');
        }

        function regionName(geneSymbol, species, start, end) {
          const geneToken = String(geneSymbol || 'GENE').toUpperCase().replace(/[^A-Z0-9]+/g, '_');
          return `${geneToken}_${speciesCode(species)}_${start}-${end}`;
        }

        function parseMutationSuggestion(text) {
          const match = String(text || '').trim().match(/^([A-Za-z*])(\\d+)([A-Za-z*])$/);
          if (!match) return null;
          const position = Number(match[2]);
          if (!Number.isFinite(position)) return null;
          return { from: match[1].toUpperCase(), position, to: match[3].toUpperCase(), label: `${match[1].toUpperCase()}${position}${match[3].toUpperCase()}` };
        }

        function mutationEditKey(edit) {
          return `${String(edit.from || 'X').toUpperCase()}${Number(edit.position)}${String(edit.to || '').toUpperCase()}`;
        }

        function mutationEditLabel(edit) {
          return edit.label || mutationEditKey(edit);
        }

        function normalizeMutationEdits(edits) {
          const deduped = new Map();
          (edits || []).forEach((edit) => {
            if (!edit) return;
            const position = Number(edit.position);
            const to = String(edit.to || '').toUpperCase();
            if (!Number.isFinite(position) || !to) return;
            const from = String(edit.from || '').toUpperCase();
            const normalized = { ...edit, position, from, to, label: edit.label || `${from || 'X'}${position}${to}` };
            deduped.set(mutationEditKey(normalized), normalized);
          });
          return Array.from(deduped.values()).sort((a, b) => a.position - b.position || mutationEditLabel(a).localeCompare(mutationEditLabel(b)));
        }

        function cysteineEditsFromPositions(positions) {
          return normalizeMutationEdits((positions || []).map((position) => ({ position, from: 'C', to: 'S', label: `C${Number(position)}S`, kind: 'cysteine' })));
        }

        function furinSiteKey(site) {
          return `${Number(site.start)}-${Number(site.end)}:${String(site.motif || '')}`;
        }

        function furinEditsForSite(site) {
          const siteKey = furinSiteKey(site);
          const edits = (site.suggestions || [])
            .map(parseMutationSuggestion)
            .filter(Boolean)
            .map((edit) => ({ ...edit, kind: 'furin', siteKey, siteStart: Number(site.start), siteEnd: Number(site.end), motif: site.motif || '' }));
          return normalizeMutationEdits(edits);
        }

        function nameWithMutationEdits(geneSymbol, species, start, end, edits) {
          const suffix = normalizeMutationEdits(edits).map(mutationEditLabel).join('_');
          return `${regionName(geneSymbol, species, start, end)}${suffix ? `_${suffix}` : ''}`;
        }

        function mutateSequenceByGlobalEdits(sequenceText, boundaryStart, edits) {
          const chars = String(sequenceText || '').split('');
          const appliedEdits = [];
          normalizeMutationEdits(edits).forEach((edit) => {
            const localIndex = edit.position - Number(boundaryStart);
            if (localIndex < 0 || localIndex >= chars.length) return;
            if (edit.from && String(chars[localIndex] || '').toUpperCase() !== edit.from) return;
            chars[localIndex] = edit.to;
            appliedEdits.push(edit);
          });
          return { sequence: chars.join(''), appliedEdits };
        }

        function homologPositionLookup(row) {
          const key = `${String(row.speciesRaw || '').toLowerCase()}:${row.start}-${row.end}`;
          if (homologPositionLookupCache.has(key)) return homologPositionLookupCache.get(key);
          const targetRegionStart = Number(row.targetRegionStart || 1);
          const homologRegionStart = Number(row.homologRegionStart || row.start || 1);
          const lookup = new Map();
          let targetPosition = 0;
          let homologPosition = 0;
          const alignedTarget = String(row.alignedTarget || '');
          const alignedHomolog = String(row.alignedHomolog || '');
          for (let index = 0; index < alignedTarget.length; index += 1) {
            const targetResidue = alignedTarget[index];
            const homologResidue = alignedHomolog[index];
            if (targetResidue !== '-') targetPosition += 1;
            if (homologResidue !== '-') homologPosition += 1;
            if (targetResidue === '-' || homologResidue === '-') continue;
            lookup.set(
              targetRegionStart + targetPosition - 1,
              homologRegionStart + homologPosition - 1
            );
          }
          homologPositionLookupCache.set(key, lookup);
          return lookup;
        }

        function mapTargetMutationEditsToRow(row, edits) {
          const normalized = String(row.speciesRaw || '').toLowerCase();
          const mutationEdits = normalizeMutationEdits(edits);
          if (normalized === targetSpecies) return mutationEdits;
          const lookup = homologPositionLookup(row);
          if (!lookup) return [];
          return normalizeMutationEdits(
            mutationEdits
              .map((edit) => {
                const mappedPosition = lookup.get(edit.position);
                if (!Number.isFinite(mappedPosition)) return null;
                return {
                  ...edit,
                  position: mappedPosition,
                  label: `${edit.from || 'X'}${mappedPosition}${edit.to}`,
                  mappedFromTargetPosition: edit.position,
                };
              })
              .filter(Boolean)
          );
        }

        function renderConstructMutationCard(card) {
          const payloadEl = card.querySelector('.construct-payload');
          if (!payloadEl) return;
          let payload = null;
          try {
            payload = JSON.parse(payloadEl.textContent || '{}');
          } catch (error) {
            throw new Error(`Invalid construct payload: ${error.message}`);
          }
          const controlsEl = card.querySelector('.construct-mutation-controls');
          const tableBody = card.querySelector('.construct-sequence-table-body');
          const tsvEl = card.querySelector('.construct-export-tsv');
          const fastaEl = card.querySelector('.construct-export-fasta');
          if (!controlsEl || !tableBody || !tsvEl || !fastaEl) return;
          const activeMutations = new Set();
          const activeFurinSites = new Set();
          const furinMutationCandidates = (payload.furinSites || []).filter((site) => furinEditsForSite(site).length);

          function redraw() {
            const targetMutationEdits = normalizeMutationEdits(
              cysteineEditsFromPositions(Array.from(activeMutations)).concat(
                furinMutationCandidates.filter((site) => activeFurinSites.has(furinSiteKey(site))).flatMap(furinEditsForSite)
              )
            );
            const rows = (payload.rows || []).map((row) => {
              const rowMutationEdits = mapTargetMutationEditsToRow(
                row,
                targetMutationEdits
              );
              const mutated = mutateSequenceByGlobalEdits(row.sequence || '', Number(row.start), rowMutationEdits);
              return {
                name: nameWithMutationEdits(payload.geneSymbol, row.speciesRaw, row.start, row.end, mutated.appliedEdits),
                species: row.species || speciesCode(row.speciesRaw),
                boundary: `${row.start}-${row.end}`,
                sequence: mutated.sequence,
                identity: row.identity || '',
                surfaceIdentity: row.surfaceIdentity || '',
                links: row.links || '<span class="muted">n/a</span>',
              };
            });
            const escapeHtml = (value) => String(value == null ? '' : value).replace(/[&<>"']/g, (ch) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[ch]));
            tableBody.innerHTML = rows.map((row) => `
              <tr>
                <td><code>${escapeHtml(row.name)}</code></td>
                <td>${escapeHtml(row.species)}</td>
                <td>${escapeHtml(row.boundary)}</td>
                <td><code class="sequence">${escapeHtml(row.sequence)}</code></td>
                <td>${escapeHtml(row.identity)}</td>
                <td>${escapeHtml(row.surfaceIdentity)}</td>
                <td>${row.links}</td>
              </tr>
            `).join('') + (payload.unavailableRows || []).map((row) => `
              <tr>
                <td><code>${escapeHtml(row.name)}</code></td>
                <td>${escapeHtml(row.species)}</td>
                <td>unavailable</td>
                <td><span class="muted">Mapping unavailable: ${escapeHtml(row.reason)}</span></td>
                <td>n/a</td><td>n/a</td><td>${row.links || '<span class="muted">n/a</span>'}</td>
              </tr>
            `).join('');
            const identityHeader = payload.identityHeader || 'identity_to_human';
            const surfaceIdentityHeader = payload.surfaceIdentityHeader || 'extracellular_surface_identity_to_human';
            const header = ['name', 'species', 'boundary', 'sequence', identityHeader, surfaceIdentityHeader];
            tsvEl.value = [header.join('\\t')].concat(rows.map((row) => [row.name, row.species, row.boundary, row.sequence, row.identity, row.surfaceIdentity].join('\\t'))).join('\\n');
            fastaEl.value = rows.map((row) => `>${row.name} species=${row.species} boundary=${row.boundary} ${identityHeader}=${row.identity || 'n/a'} ${surfaceIdentityHeader}=${row.surfaceIdentity || 'n/a'}\\n${row.sequence.replace(/(.{1,80})/g, '$1\\n').trim()}`).join('\\n');
            if (!(payload.mutationCandidates || []).length && !furinMutationCandidates.length) {
              controlsEl.innerHTML = '<p class="viewer-info">No unpaired cysteines or furin sites in this construct.</p>';
              return;
            }
            const cysteineControls = (payload.mutationCandidates || []).length
              ? `
                <div>
                  <strong>Optional Cys→Ser edits</strong>
                  <div class="mutation-checkboxes">
                    ${payload.mutationCandidates.map((position) => `
                      <label class="mutation-checkbox">
                        <input type="checkbox" class="construct-mutation-checkbox" data-position="${position}" ${activeMutations.has(position) ? 'checked' : ''}>
                        C${position}S
                      </label>
                    `).join('')}
                  </div>
                </div>
              `
              : '';
            const furinControls = furinMutationCandidates.length
              ? `
                <div>
                  <strong>Optional furin-site edits</strong>
                  <div class="mutation-checkboxes">
                    ${furinMutationCandidates.map((site) => {
                      const siteKey = furinSiteKey(site);
                      const edits = furinEditsForSite(site).map(mutationEditLabel).join(', ');
                      return `
                        <label class="mutation-checkbox">
                          <input type="checkbox" class="construct-furin-mutation-checkbox" data-site-key="${siteKey}" ${activeFurinSites.has(siteKey) ? 'checked' : ''}>
                          ${escapeHtml(site.motif || 'furin site')} ${site.start}-${site.end} (${edits})
                        </label>
                      `;
                    }).join('')}
                  </div>
                </div>
              `
              : '';
            controlsEl.innerHTML = `${cysteineControls}${furinControls}`;
            controlsEl.querySelectorAll('.construct-mutation-checkbox').forEach((checkbox) => {
              checkbox.addEventListener('change', () => {
                const position = Number(checkbox.dataset.position);
                if (checkbox.checked) activeMutations.add(position);
                else activeMutations.delete(position);
                redraw();
              });
            });
            controlsEl.querySelectorAll('.construct-furin-mutation-checkbox').forEach((checkbox) => {
              checkbox.addEventListener('change', () => {
                const siteKey = checkbox.dataset.siteKey || '';
                if (checkbox.checked) activeFurinSites.add(siteKey);
                else activeFurinSites.delete(siteKey);
                redraw();
              });
            });
          }

          redraw();
        }

        const initializeCards = () => document.querySelectorAll('.construct-card').forEach(renderConstructMutationCard);
        if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', initializeCards, { once: true });
        else initializeCards();
      })();
    """
    )


def _render_cross_reactivity_rows(
    hits: list[dict[str, Any]],
    *,
    batch_dir: Path,
    query_label: str = "human",
) -> str:
    rows: list[str] = []
    for index, item in enumerate(hits, start=1):
        cross_target = _fetch_cross_reactivity_target(item, batch_dir=batch_dir) or {}
        alignment_html = _render_cross_reactivity_alignment(
            item,
            rank=index,
            query_label=query_label,
        )
        hit_links = _render_resource_links(
            accession=str(cross_target.get("accession") or extract_accession(str(item.get("subject_id") or "")) or ""),
            entry_name=str(cross_target.get("entry_name") or extract_entry_name(str(item.get("subject_id") or "")) or ""),
            gene_symbol=str(cross_target.get("gene_symbol") or ""),
            species=str(item.get("species") or ""),
        )
        hit_label = escape(str(item.get("subject_id") or ""))
        if item.get("description"):
            hit_label += f"<div class=\"muted hit-description\">{escape(str(item.get('description') or ''))}</div>"
        rows.append(
            "<tr><td>{rank}</td><td><details><summary>{hit}</summary>{alignment}</details></td><td>{species}</td><td>{bitscore}</td><td>{identity}</td><td>{surface_identity}</td><td>{coverage}</td><td>{evalue}</td><td>{links}</td></tr>".format(
                rank=escape(str(index)),
                hit=hit_label,
                alignment=alignment_html,
                species=escape(_species_display_label(item.get("species"))),
                bitscore=escape(_fmt_number(item.get("bitscore"))),
                identity=escape(_fmt_pct(item.get("identity"))),
                surface_identity=escape(_fmt_surface_identity(item)),
                coverage=escape(_fmt_pct(item.get("coverage"))),
                evalue=escape(_fmt_evalue(item.get("evalue"))),
                links=hit_links,
            )
        )
    return "".join(rows) or "<tr><td colspan=\"9\">No sequence similarity hits were found in this BLAST search.</td></tr>"


def _report_ptms(report: dict[str, Any]) -> list[dict[str, Any]]:
    explicit = report.get("ptms")
    if explicit:
        return list(explicit)
    return [feature for feature in (report.get("features") or []) if _is_ptm_feature_dict(feature)]


def _is_ptm_feature_dict(feature: dict[str, Any]) -> bool:
    feature_type = str(feature.get("type") or "").upper()
    description = str(feature.get("description") or "").lower()
    if feature_type in {
        "CARBOHYD",
        "MOD_RES",
        "LIPID",
        "DISULFID",
        "CROSSLNK",
        "INIT_MET",
        "SIGNAL",
        "PROPEP",
        "PEPTIDE",
        "CHAIN",
    }:
        return True
    return feature_type == "SITE" and any(
        keyword in description
        for keyword in ("cleavage", "glycosyl", "phospho", "acetyl", "methyl", "ubiquitin", "sumoyl", "lipid")
    )


def _ptm_category(feature: dict[str, Any]) -> str:
    metadata = feature.get("metadata") or {}
    if metadata.get("ptm_category"):
        return str(metadata["ptm_category"])
    feature_type = str(feature.get("type") or "").upper()
    description = str(feature.get("description") or "").lower()
    if feature_type == "CARBOHYD":
        return "glycosylation"
    if feature_type == "MOD_RES":
        return "modified residue"
    if feature_type == "LIPID":
        return "lipidation"
    if feature_type == "DISULFID":
        return "disulfide"
    if feature_type == "CROSSLNK":
        return "cross-link"
    if feature_type == "INIT_MET":
        return "initiator methionine processing"
    if feature_type in {"SIGNAL", "PROPEP", "PEPTIDE", "CHAIN"}:
        return "molecule processing"
    if "cleavage" in description:
        return "cleavage site"
    return "post-translational modification"


def _ptm_warning(feature: dict[str, Any]) -> str:
    feature_type = str(feature.get("type") or "").upper()
    category = _ptm_category(feature).lower()
    description = str(feature.get("description") or "").lower()
    if feature_type in {"SIGNAL", "PROPEP", "PEPTIDE"} or "cleavage" in category or "cleavage" in description:
        return "Review cleavage/processing boundary."
    return ""


def _closest_unpaired_cysteine_text(item: dict[str, Any]) -> str:
    if item.get("paired_with") is not None:
        return ""
    partner = item.get("closest_unpaired_cysteine")
    distance = item.get("closest_unpaired_cysteine_distance")
    if partner is None or distance is None:
        return "none detected"
    atom = str(item.get("closest_unpaired_cysteine_distance_atom") or "distance")
    return f"C{partner}, {_fmt_number(distance)} A ({atom})"


def _ptm_display_label(feature: dict[str, Any]) -> str:
    feature_type = str(feature.get("type") or "").upper()
    try:
        start = int(feature.get("start"))
        end = int(feature.get("end"))
    except (TypeError, ValueError):
        start = feature.get("start")
        end = feature.get("end")
    if feature_type == "DISULFID":
        return f"C{start}-C{end}"
    return f"{start}-{end}"


def _render_ptm_rows(ptms: list[dict[str, Any]]) -> str:
    rows = []
    for feature in sorted(ptms, key=lambda item: (int(item.get("start") or 0), int(item.get("end") or 0), str(item.get("type") or ""))):
        rows.append(
            "<tr><td>{start}</td><td>{end}</td><td>{type}</td><td>{category}</td><td>{description}</td><td>{warning}</td></tr>".format(
                start=escape(str(feature.get("start") or "")),
                end=escape(str(feature.get("end") or "")),
                type=escape(str(feature.get("type") or "")),
                category=escape(_ptm_category(feature)),
                description=escape(str(feature.get("description") or "")),
                warning=escape(_ptm_warning(feature)),
            )
        )
    return "".join(rows) or "<tr><td colspan=\"6\">No curated UniProt PTM annotations available.</td></tr>"


def _features_in_region(features: list[dict[str, Any]], start: int | None, end: int | None) -> list[dict[str, Any]]:
    if start is None or end is None:
        return []
    selected = []
    for feature in features:
        try:
            feature_start = int(feature.get("start"))
            feature_end = int(feature.get("end"))
        except (TypeError, ValueError):
            continue
        if _feature_overlaps_region(feature, feature_start=feature_start, feature_end=feature_end, start=start, end=end):
            selected.append(feature)
    return selected


def _feature_overlaps_region(
    feature: dict[str, Any],
    *,
    feature_start: int,
    feature_end: int,
    start: int,
    end: int,
) -> bool:
    if str(feature.get("type") or "").upper() == "DISULFID":
        # Disulfide annotations are cysteine-pair endpoints, not continuous spans.
        return start <= feature_start <= end or start <= feature_end <= end
    return feature_start <= end and feature_end >= start


def _render_family_context_section(
    family_context: dict[str, Any] | None,
    *,
    canonical_family: dict[str, Any] | None = None,
    matrix_scope_label: str = "ectodomain",
) -> str:
    if not family_context and not canonical_family:
        return "<p>No family context available.</p>"
    sections: list[str] = []
    if canonical_family:
        coverage = canonical_family.get("coverage_fraction")
        coverage_text = _fmt_pct(coverage) if coverage is not None else ""
        span = ""
        if canonical_family.get("start") is not None and canonical_family.get("end") is not None:
            span = f"{canonical_family.get('start')}-{canonical_family.get('end')}"
        sections.append(
            "<p><strong>Canonical family:</strong> {name} ({accession})</p>"
            "<p><strong>Source:</strong> {source} | <strong>Fragments:</strong> {fragments} | <strong>Span:</strong> {span} | <strong>Coverage:</strong> {coverage}</p>".format(
                name=escape(str(canonical_family.get("name") or "")),
                accession=escape(str(canonical_family.get("accession") or "")),
                source=escape(str(canonical_family.get("source_database") or "")),
                fragments=escape(str(canonical_family.get("fragment_count") or "")),
                span=escape(span or "n/a"),
                coverage=escape(coverage_text or "n/a"),
            )
        )
    if not family_context:
        return "".join(sections)
    family_names = family_context.get("family_names") or []
    members = family_context.get("members") or []
    labels = family_context.get("identity_matrix_labels") or []
    matrix = family_context.get("identity_matrix") or []
    coverage_labels = family_context.get("coverage_matrix_labels") or []
    coverage_matrix = family_context.get("coverage_matrix") or []
    metadata = family_context.get("metadata") or {}
    accessible_identity_rows = metadata.get("target_extracellular_accessible_identity") if isinstance(metadata, dict) else None
    mouse_matrix_caveat = metadata.get("mouse_ortholog_matrix_caveat") if isinstance(metadata, dict) else None
    member_rows = "".join(
        "<tr><td>{gene}</td><td>{entry}</td><td>{ecto}</td></tr>".format(
            gene=escape(str(item.get("gene_symbol") or "")),
            entry=escape(str(item.get("entry_name") or item.get("accession") or "")),
            ecto=escape(
                f"{item.get('ectodomain_start')}-{item.get('ectodomain_end')}"
                if item.get("ectodomain_start") is not None and item.get("ectodomain_end") is not None
                else ""
            ),
        )
        for item in members[:30]
    ) or "<tr><td colspan=\"3\">No family members.</td></tr>"
    sections.extend([
        f"<p><strong>Source:</strong> {escape(str(family_context.get('source') or ''))}</p>",
        f"<p><strong>Family:</strong> {escape(', '.join(family_names) if family_names else 'n/a')}</p>",
        "<table><thead><tr><th>Gene</th><th>Entry</th><th>Design region</th></tr></thead><tbody>{}</tbody></table>".format(member_rows),
    ])
    if mouse_matrix_caveat:
        sections.append(f'<p class="section-note">{escape(str(mouse_matrix_caveat))}</p>')
    if labels and matrix:
        header = "".join(f"<th>{escape(str(label))}</th>" for label in labels)
        rows = []
        for label, row in zip(labels, matrix):
            values = "".join(_render_identity_matrix_cell(value) for value in row)
            rows.append(f"<tr><th>{escape(str(label))}</th>{values}</tr>")
        sections.append("<h3>Directional Identity Matrix</h3>")
        sections.append(
            f"<p class=\"section-note\">Directional matrix: each cell reports row protein → column protein identity, normalized by the row {escape(matrix_scope_label)} length.</p>"
        )
        sections.append(
            "<div class=\"matrix-wrap\"><table><thead><tr><th>Protein</th>{}</tr></thead><tbody>{}</tbody></table></div>".format(
                header,
                "".join(rows),
            )
        )
    if isinstance(accessible_identity_rows, list) and accessible_identity_rows:
        rows = "".join(
            "<tr><td>{member}</td><td>{identity}</td><td>{aligned}</td><td>{matches}</td></tr>".format(
                member=escape(str(item.get("member_label") or "")),
                identity=escape(_fmt_pct(item.get("identity"))),
                aligned=escape(str(item.get("aligned_positions") or "")),
                matches=escape(str(item.get("matches") or "")),
            )
            for item in accessible_identity_rows
            if isinstance(item, dict)
        )
        if rows:
            sections.append("<h3>Target Extracellular Accessible Identity</h3>")
            sections.append(
                "<p class=\"section-note\">Target-perspective identity over residues classified as extracellular and accessible "
                "(SASA-exposed or conditionally exposed within a PAE-coherent block).</p>"
            )
            sections.append(
                "<table><thead><tr><th>Family member</th><th>Identity</th><th>Aligned accessible residues</th><th>Matches</th></tr></thead>"
                f"<tbody>{rows}</tbody></table>"
            )
    if coverage_labels and coverage_matrix:
        header = "".join(f"<th>{escape(str(label))}</th>" for label in coverage_labels)
        rows = []
        for label, row in zip(coverage_labels, coverage_matrix):
            values = "".join(_render_identity_matrix_cell(value) for value in row)
            rows.append(f"<tr><th>{escape(str(label))}</th>{values}</tr>")
        sections.append("<h3>Directional Coverage Matrix</h3>")
        sections.append(
            f"<p class=\"section-note\">Directional coverage matrix: each cell reports the percent of the row {escape(matrix_scope_label)} aligned to the column {escape(matrix_scope_label)}.</p>"
        )
        sections.append(
            "<div class=\"matrix-wrap\"><table><thead><tr><th>Protein</th>{}</tr></thead><tbody>{}</tbody></table></div>".format(
                header,
                "".join(rows),
            )
        )
    return "".join(sections)


def _render_gpcr_annotation_section(annotation: dict[str, Any] | None) -> str:
    if not annotation:
        return ""
    gpcr_url = str(annotation.get("url") or "")
    link = f'<a class="inline-link" href="{escape(gpcr_url)}">GPCRdb</a>' if gpcr_url else '<span class="muted">n/a</span>'
    family_path = " > ".join(str(item) for item in annotation.get("family_path") or []) or str(annotation.get("family_name") or "n/a")
    summary_rows = """
      <tr><th>GPCRdb entry</th><td><code>{entry}</code></td></tr>
      <tr><th>Name</th><td>{name}</td></tr>
      <tr><th>Class / family</th><td>{family}</td></tr>
      <tr><th>Numbering</th><td>{numbering}</td></tr>
      <tr><th>Link</th><td>{link}</td></tr>
    """.format(
        entry=escape(str(annotation.get("entry_name") or "")),
        name=escape(str(annotation.get("name") or "n/a")),
        family=escape(family_path),
        numbering=escape(str(annotation.get("residue_numbering_scheme") or "n/a")),
        link=link,
    )
    segment_rows = "".join(
        "<tr><td>{name}</td><td>{start}-{end}</td><td>{generic_start}</td><td>{generic_end}</td></tr>".format(
            name=escape(str(segment.get("name") or "")),
            start=escape(str(segment.get("start") or "")),
            end=escape(str(segment.get("end") or "")),
            generic_start=escape(str(segment.get("generic_start") or "")),
            generic_end=escape(str(segment.get("generic_end") or "")),
        )
        for segment in annotation.get("segments") or []
    ) or '<tr><td colspan="4">No GPCRdb segment mapping available.</td></tr>'
    motif_rows = "".join(
        "<tr><td>{name}</td><td>{positions}</td><td>{generic}</td><td>{sequence}</td><td>{status}</td></tr>".format(
            name=escape(str(motif.get("name") or "")),
            positions=escape(", ".join(str(item) for item in motif.get("positions") or [])),
            generic=escape(", ".join(str(item) for item in motif.get("generic_numbers") or [])),
            sequence=escape(str(motif.get("sequence") or "")),
            status=escape(str(motif.get("status") or "")),
        )
        for motif in annotation.get("conserved_motifs") or []
    ) or '<tr><td colspan="5">No conserved motif positions detected from GPCRdb numbering.</td></tr>'
    warnings = annotation.get("warnings") or []
    warning_html = ""
    if warnings:
        warning_html = "<p class=\"section-note\">" + " ".join(escape(str(item)) for item in warnings) + "</p>"
    return f"""
    <section class="card" id="gpcr-annotation" data-report-section="GPCR Annotation">
      <h2>GPCR Annotation</h2>
      <p class="section-note">GPCRdb annotations provide receptor-specific segment names and generic numbering. These annotations improve interpretation of membrane constructs but do not imply expression success.</p>
      {warning_html}
      <div class="grid two">
        <div>
          <table class="compact-table"><tbody>{summary_rows}</tbody></table>
        </div>
        <div>
          <h3>Conserved Motifs</h3>
          <table>
            <thead><tr><th>Motif</th><th>Residues</th><th>Generic numbers</th><th>Sequence</th><th>Status</th></tr></thead>
            <tbody>{motif_rows}</tbody>
          </table>
        </div>
      </div>
      <details class="details-panel">
        <summary>GPCRdb segment table</summary>
        <div class="table-scroll">
          <table>
            <thead><tr><th>Segment</th><th>Boundary</th><th>Generic start</th><th>Generic end</th></tr></thead>
            <tbody>{segment_rows}</tbody>
          </table>
        </div>
      </details>
    </section>
    """


def _render_advanced_membrane_suggestions_section(suggestions: list[dict[str, Any]]) -> str:
    if not suggestions:
        return "<p>No advanced membrane-engineering suggestions available.</p>"
    cards: list[str] = []
    for suggestion in suggestions:
        rationale = "".join(
            f"<li>{escape(str(item))}</li>"
            for item in (suggestion.get("rationale") or [])
        ) or "<li>None</li>"
        actions = "".join(
            f"<li>{escape(str(item))}</li>"
            for item in (suggestion.get("suggested_actions") or [])
        ) or "<li>None</li>"
        evidence = "".join(
            f"<li>{escape(str(item))}</li>"
            for item in (suggestion.get("evidence") or [])
        ) or "<li>None</li>"
        warnings = "".join(
            f"<li>{escape(str(item))}</li>"
            for item in (suggestion.get("warnings") or [])
        ) or "<li>None</li>"
        region = ""
        if suggestion.get("start") is not None and suggestion.get("end") is not None:
            region = f"<p><strong>Region:</strong> {escape(str(suggestion.get('start')))}-{escape(str(suggestion.get('end')))}</p>"
        cards.append(
            """
            <article class="construct-card">
              <div class="construct-copy">
                <h3>{title}</h3>
                <p><strong>Category:</strong> <code>{category}</code></p>
                {region}
                <p><strong>Summary:</strong> {summary}</p>
                <div class="construct-grid">
                  <div>
                    <h4>Rationale</h4>
                    <ul>{rationale}</ul>
                  </div>
                  <div>
                    <h4>Experiments to Consider</h4>
                    <ul>{actions}</ul>
                  </div>
                  <div>
                    <h4>Evidence</h4>
                    <ul>{evidence}</ul>
                  </div>
                  <div>
                    <h4>Warnings</h4>
                    <ul>{warnings}</ul>
                  </div>
                </div>
              </div>
            </article>
            """.format(
                title=escape(str(suggestion.get("title") or "")),
                category=escape(str(suggestion.get("category") or "")),
                region=region,
                summary=escape(str(suggestion.get("summary") or "")),
                rationale=rationale,
                actions=actions,
                evidence=evidence,
                warnings=warnings,
            )
        )
    return '<div class="construct-list">' + "".join(cards) + "</div>"


def _public_gpcr_variant(variant: dict[str, Any]) -> dict[str, Any]:
    public_variant = dict(variant)
    loop_source = str(variant.get("loop_source") or "").strip()
    evidence = [str(item) for item in (variant.get("evidence") or [])]
    topology_only = (
        loop_source == "topology-derived internal cytoplasmic loop"
        or any(
            item.startswith("topology-derived internal cytoplasmic loop:")
            for item in evidence
        )
    )
    if topology_only:
        public_variant["name"] = str(variant.get("name") or "").replace("_ICL3_", "_ICL_")
        if variant.get("strategy") == "ICL3 cassette replacement":
            public_variant["strategy"] = "Intracellular-loop cassette replacement"
    return public_variant


def _render_gpcr_engineering_variants_section(variants: list[dict[str, Any]]) -> str:
    if not variants:
        return ""
    rows: list[str] = []
    cards: list[str] = []
    for variant in variants:
        variant = _public_gpcr_variant(variant)
        name = str(variant.get("name") or "")
        start = variant.get("start")
        end = variant.get("end")
        length = variant.get("engineered_length")
        cassette = variant.get("cassette_name") or variant.get("cassette_id") or "n/a"
        sequence = str(variant.get("sequence") or "")
        focus_button = (
            _viewer_focus_button(start, end, name, compact=True)
            or '<span class="muted">n/a</span>'
        )
        rows.append(
            "<tr><td><code>{name}</code></td><td>{category}</td><td>{strategy}</td><td>{region}</td><td>{cassette}</td><td>{length}</td><td>{action}</td></tr>".format(
                name=escape(name),
                category=escape(str(variant.get("category") or "")),
                strategy=escape(str(variant.get("strategy") or "")),
                region=escape(f"{start}-{end}" if start is not None and end is not None else "n/a"),
                cassette=escape(str(cassette)),
                length=escape(str(length or "")),
                action=focus_button,
            )
        )
        warnings = "".join(f"<li>{escape(str(item))}</li>" for item in (variant.get("warnings") or [])) or "<li>None</li>"
        evidence = "".join(f"<li>{escape(str(item))}</li>" for item in (variant.get("evidence") or [])) or "<li>None</li>"
        actions = "".join(f"<li>{escape(str(item))}</li>" for item in (variant.get("suggested_actions") or [])) or "<li>None</li>"
        replacement = ""
        if variant.get("replaced_start") is not None and variant.get("replaced_end") is not None:
            replacement = (
                f"<p><strong>Replaced receptor residues:</strong> "
                f"{escape(str(variant.get('replaced_start')))}-{escape(str(variant.get('replaced_end')))}</p>"
            )
        linker = ""
        if variant.get("n_linker") is not None or variant.get("c_linker") is not None:
            linker = (
                f"<p><strong>Linkers:</strong> N <code>{escape(str(variant.get('n_linker') or ''))}</code>; "
                f"C <code>{escape(str(variant.get('c_linker') or ''))}</code></p>"
            )
        generic = ""
        if variant.get("gpcr_generic_range"):
            generic = f"<p><strong>GPCRdb generic range:</strong> {escape(str(variant.get('gpcr_generic_range')))}</p>"
        tag_suggestions = variant.get("tag_suggestions") or {}
        tag_items = []
        for label, values in tag_suggestions.items():
            joined = "; ".join(str(item) for item in values) if isinstance(values, list) else str(values)
            if joined:
                tag_items.append(f"<li><strong>{escape(str(label).replace('_', ' ').title())}:</strong> {escape(joined)}</li>")
        tags_html = "".join(tag_items) or "<li>No tag metadata.</li>"
        fasta_text = f">{name}\n{_wrap_plain_sequence(sequence)}" if sequence else ""
        tsv_text = _gpcr_variant_tsv(variant)
        export_html = ""
        if sequence:
            export_html = f"""
              <div class="construct-export">
                <h4>Copy-ready exports</h4>
                <label>TSV</label>
                <textarea class="selection-tsv construct-export-box construct-export-tsv" readonly>{escape(tsv_text)}</textarea>
                <label>FASTA</label>
                <textarea class="selection-tsv construct-export-box construct-export-fasta" readonly>{escape(fasta_text)}</textarea>
              </div>
            """
        else:
            export_html = '<p class="section-note">No sequence generated for this warning/advisory variant.</p>'
        cards.append(
            f"""
            <article class="construct-card">
              <div class="construct-copy">
                <h3><code>{escape(name)}</code></h3>
                <p><strong>Category:</strong> {escape(str(variant.get('category') or ''))}</p>
                <p><strong>Summary:</strong> {escape(str(variant.get('summary') or ''))}</p>
                {replacement}
                {linker}
                {generic}
                <div class="construct-grid">
                  <div><h4>Evidence</h4><ul>{evidence}</ul></div>
                  <div><h4>Experiments to Consider</h4><ul>{actions}</ul></div>
                  <div><h4>Warnings</h4><ul>{warnings}</ul></div>
                  <div><h4>Tag Metadata</h4><ul>{tags_html}</ul></div>
                </div>
                {export_html}
              </div>
            </article>
            """
        )
    return f"""
    <section class="card" id="gpcr-engineering" data-report-section="GPCR Engineering">
      <h2>GPCR Expression and Stabilization Strategy</h2>
      <p class="section-note">These are experimental GPCR expression/stability strategies for antibody-discovery workflows. They are intentionally separated from native OpenAntigens construct recommendations because loop fusions and partner cassettes can change receptor conformation, signaling state, trafficking, and epitope presentation.</p>
      <h3>Experimental GPCR Engineering Variants</h3>
      <div class="table-scroll">
        <table>
          <thead><tr><th>Variant</th><th>Category</th><th>Strategy</th><th>Focus region</th><th>Cassette</th><th>Length</th><th>Action</th></tr></thead>
          <tbody>{''.join(rows)}</tbody>
        </table>
      </div>
      <div class="construct-list">{''.join(cards)}</div>
    </section>
    """


def _wrap_plain_sequence(sequence: str, width: int = 80) -> str:
    return "\n".join(sequence[index : index + width] for index in range(0, len(sequence), width))


def _gpcr_variant_tsv(variant: dict[str, Any]) -> str:
    header = [
        "name",
        "category",
        "strategy",
        "sequence_role",
        "region",
        "cassette",
        "replaced_region",
        "n_linker",
        "c_linker",
        "engineered_length",
        "sequence",
    ]
    region = (
        f"{variant.get('start')}-{variant.get('end')}"
        if variant.get("start") is not None and variant.get("end") is not None
        else ""
    )
    replaced = (
        f"{variant.get('replaced_start')}-{variant.get('replaced_end')}"
        if variant.get("replaced_start") is not None and variant.get("replaced_end") is not None
        else ""
    )
    row = [
        str(variant.get("name") or ""),
        str(variant.get("category") or ""),
        str(variant.get("strategy") or ""),
        str(variant.get("sequence_role") or ""),
        region,
        str(variant.get("cassette_name") or variant.get("cassette_id") or ""),
        replaced,
        str(variant.get("n_linker") or ""),
        str(variant.get("c_linker") or ""),
        str(variant.get("engineered_length") or ""),
        str(variant.get("sequence") or ""),
    ]
    return "\t".join(header) + "\n" + "\t".join(row)


def _construct_workbench_js() -> str:
    """Return the shared DOM adapter for the integrated construct workbench."""

    return r"""
        let workbenchSelectionRows = [];
        let renderWorkbenchSpecies = () => {};
        function scheduleConstructWorkbench() {
          if (document.readyState === "loading") {
            document.addEventListener("DOMContentLoaded", installConstructWorkbench, { once: true });
            return;
          }
          installConstructWorkbench();
        }

        function installConstructWorkbench() {
          const builder = document.getElementById("interactive-construct-builder");
          const viewerLayout = builder?.querySelector(".viewer-layout");
          const viewerPanel = viewerLayout?.querySelector(".viewer-panel");
          const sequencePanel = viewerLayout?.querySelector(".sequence-panel");
          const plotCards = Array.from(builder?.querySelectorAll(".interactive-plots .plot-card") || []);
          const livePanel = builder?.querySelector(".selection-live-card");
          if (!builder || !viewerPanel || !sequencePanel || plotCards.length !== 2 || !livePanel) {
            throw new Error("The construct workbench could not find its report controls.");
          }
          if (builder.classList.contains("construct-workbench")) return;

          builder.classList.add("construct-workbench");
          document.body.classList.add("construct-workbench-active");

          const header = document.createElement("div");
          header.className = "construct-workbench-header";
          const heading = document.createElement("h2");
          heading.textContent = "Construct Builder";
          header.append(heading, builder.querySelector(".construct-choice-help"));

          const picker = makeConstructPicker();
          const sequenceMetadata = readSequenceMetadata();
          const boundaryEditor = makeBoundaryEditor(sequenceMetadata);
          const speciesView = makeCompactSpeciesView(boundaryEditor);
          compactLivePanel(boundaryEditor, speciesView);
          labelStructureTools();

          plotCards[0].querySelector("h3").textContent = "pLDDT confidence";
          plotCards[1].querySelector("h3").textContent = "PAE matrix";
          installWorkbenchPlddtRenderer(plotCards[0]);

          const center = document.createElement("main");
          center.className = "construct-workbench-center";
          const evidence = document.createElement("div");
          evidence.className = "construct-workbench-evidence";
          evidence.append(...plotCards);
          center.append(viewerPanel, evidence, makeSequenceDetails(sequenceMetadata));

          const shell = document.createElement("div");
          shell.className = "construct-workbench-shell";
          shell.append(picker, center, livePanel);
          builder.replaceChildren(header, shell);

          function findConstructSources() {
            const titles = new Set([
              "Precomputed Construct Summary",
              "GPCR Expression and Stabilization Strategy",
            ]);
            return Array.from(document.querySelectorAll("section.card")).filter(
              (section) => titles.has(section.querySelector(":scope > h2")?.textContent.trim()),
            );
          }

          function makeConstructPicker() {
            const pickerElement = document.createElement("aside");
            pickerElement.className = "construct-workbench-picker";
            const pickerHeading = document.createElement("div");
            pickerHeading.className = "construct-workbench-picker-heading";
            const pickerTitle = document.createElement("strong");
            pickerTitle.textContent = "Suggested constructs";
            const pickerNote = document.createElement("span");
            pickerNote.textContent = "Selecting one updates the workbench";
            pickerHeading.append(pickerTitle, pickerNote);

            const filter = document.createElement("input");
            filter.type = "search";
            filter.placeholder = "Filter constructs";
            filter.setAttribute("aria-label", "Filter suggested constructs");

            const list = document.createElement("div");
            list.className = "construct-workbench-list";
            const rows = findConstructSources()
              .flatMap((section) => Array.from(section.querySelectorAll("tbody tr")))
              .filter((row) => row.querySelector(".construct-highlight"));
            const groupOrder = ["Topology", "PDB-backed", "AlphaFold strict", "AlphaFold lenient"];
            const rowValue = (row, headings) => {
              const headers = Array.from(
                row.closest("table")?.querySelectorAll("thead th") || [],
                (cell) => cell.textContent.trim(),
              );
              const cells = row.querySelectorAll("td");
              const index = headings.map((candidate) => headers.indexOf(candidate)).find((candidate) => candidate >= 0);
              return index === undefined ? "" : cells[index]?.textContent.trim() || "";
            };
            const rank = (row) => {
              const index = groupOrder.indexOf(rowValue(row, ["Group", "Category"]));
              return index < 0 ? groupOrder.length : index;
            };
            rows.sort((a, b) => rank(a) - rank(b));

            rows.forEach((row, index) => {
              const sourceButton = row.querySelector(".construct-highlight");
              const constructName = rowValue(row, ["Construct", "Variant"]) || "Construct";
              const group = rowValue(row, ["Group", "Category"]);
              const boundary = rowValue(row, ["Boundary", "Focus region"]);
              const button = document.createElement("button");
              button.type = "button";
              button.className = "construct-workbench-option";
              button.dataset.search = row.textContent.toLowerCase();
              button.title = constructName;

              const name = document.createElement("strong");
              name.textContent = constructName;
              const category = document.createElement("small");
              category.textContent = group;
              const range = document.createElement("span");
              range.className = "construct-workbench-boundary";
              range.textContent = boundary;
              button.append(name, category, range);
              button.addEventListener("click", () => {
                list.querySelectorAll("button").forEach((item) => item.removeAttribute("aria-pressed"));
                button.setAttribute("aria-pressed", "true");
                sourceButton.click();
              });
              if (index === 0) button.setAttribute("aria-pressed", "true");
              list.append(button);
            });

            filter.addEventListener("input", () => {
              const query = filter.value.trim().toLowerCase();
              list.querySelectorAll("button").forEach((button) => {
                button.hidden = Boolean(query && !button.dataset.search.includes(query));
              });
            });

            const summary = findConstructSources()[0];
            const viewAll = document.createElement("button");
            viewAll.type = "button";
            viewAll.className = "construct-workbench-view-all";
            viewAll.textContent = `View all ${rows.length} constructs below`;
            viewAll.addEventListener("click", () => summary?.scrollIntoView({ behavior: "smooth", block: "start" }));
            pickerElement.append(pickerHeading, filter, list, viewAll);
            return pickerElement;
          }

          function readSequenceMetadata() {
            const targetName = sequencePanel.querySelector(".sequence-header strong")?.textContent.trim();
            const description = sequencePanel.querySelector(".sequence-header span")?.textContent.trim();
            const length = Number(description?.match(/^(\d+)\s+aa total/)?.[1]);
            const boundary = description?.match(/(\d+)-(\d+)/);
            if (!targetName || !Number.isInteger(length)) {
              throw new Error("The construct workbench could not read the report sequence.");
            }
            return {
              targetName,
              length,
              start: boundary ? Number(boundary[1]) : 1,
              end: boundary ? Number(boundary[2]) : length,
            };
          }

          function makeBoundaryEditor(sequenceMetadata) {
            const editor = document.createElement("div");
            editor.className = "construct-workbench-boundary-editor";
            const startLabel = document.createElement("label");
            const endLabel = document.createElement("label");
            const startInput = document.createElement("input");
            const endInput = document.createElement("input");
            startInput.type = endInput.type = "number";
            startInput.min = endInput.min = "1";
            startInput.value = String(sequenceMetadata.start);
            endInput.value = String(sequenceMetadata.end);
            startLabel.append("Start", startInput);
            endLabel.append("End", endInput);

            const apply = document.createElement("button");
            apply.type = "button";
            apply.textContent = "Apply boundaries";
            const message = document.createElement("span");
            message.className = "construct-workbench-boundary-message";
            apply.addEventListener("click", () => {
              const start = Number(startInput.value);
              const end = Number(endInput.value);
              if (!Number.isInteger(start) || !Number.isInteger(end) || start < 1 || end > sequenceMetadata.length || start > end) {
                message.textContent = `Enter a valid range from 1 to ${sequenceMetadata.length}.`;
                return;
              }
              focusRange(start, end, 'Selected range');
              message.textContent = `Applied ${start}-${end}.`;
            });
            editor.append(startLabel, endLabel, apply, message);
            return { element: editor, startInput, endInput };
          }

          function makeCompactSpeciesView(boundaryEditor) {
            const compact = document.createElement("div");
            compact.className = "construct-workbench-species";
            const render = () => {
              const populated = workbenchSelectionRows;
              compact.replaceChildren();
              if (!populated.length) {
                const empty = document.createElement("p");
                empty.className = "construct-workbench-species-empty";
                empty.textContent = "Select a residue range or a suggested construct to see live species sequences.";
                compact.append(empty);
                return;
              }
              populated.forEach((row, index) => {
                const sourceSpecies = row.species;
                const species = ["MACFA", "MACACA_FASCICULARIS"].includes(sourceSpecies.toUpperCase())
                  ? "CYNO. MONKEY"
                  : sourceSpecies;
                const item = document.createElement("div");
                item.className = "construct-workbench-species-row";
                const metadata = document.createElement("div");
                const label = document.createElement("strong");
                label.textContent = species;
                const boundary = document.createElement("span");
                boundary.textContent = row.boundary;
                const sequence = document.createElement("code");
                sequence.textContent = row.sequence;
                metadata.append(label, boundary);
                item.append(metadata, sequence);
                compact.append(item);

                if (index === 0) {
                  const selectedBoundary = boundary.textContent.match(/(\d+)-(\d+)/);
                  if (selectedBoundary) {
                    boundaryEditor.startInput.value = selectedBoundary[1];
                    boundaryEditor.endInput.value = selectedBoundary[2];
                  }
                }
              });
            };
            render();
            renderWorkbenchSpecies = render;
            return compact;
          }

          function compactLivePanel(boundaryEditor, speciesView) {
            livePanel.classList.add("construct-workbench-current");
            const title = livePanel.querySelector(".selection-live-header h3");
            if (title) title.textContent = "Current construct";
            livePanel.querySelector(".selection-live-header")?.after(boundaryEditor.element);
            livePanel.querySelector("#selectionSummary")?.after(speciesView);
            livePanel.querySelector(".selection-table-wrap")?.classList.add("construct-workbench-original-table");
            const tsv = livePanel.querySelector("#selectionTsv");
            const fasta = livePanel.querySelector("#selectionFasta");
            if (tsv && fasta) {
              const raw = document.createElement("details");
              raw.className = "construct-workbench-raw-exports";
              const summary = document.createElement("summary");
              summary.textContent = "Raw export text";
              raw.append(summary, tsv, fasta);
              livePanel.append(raw);
            }
          }

          function labelStructureTools() {
            const toolbar = viewerPanel.querySelector(".viewer-toolbar");
            if (!toolbar) return;
            const label = document.createElement("strong");
            label.className = "construct-workbench-toolbar-label";
            label.textContent = "Structure tools";
            toolbar.prepend(label);
            viewerPanel.classList.add("construct-workbench-viewer");

            const focusButton = toolbar.querySelector("#viewerEctodomain");
            if (focusButton) {
              focusButton.textContent = focusButton.textContent.includes("full-length")
                ? "Focus full length"
                : "Focus extracellular";
            }
            const selectedDownload = toolbar.querySelector("#viewerDownloadSelectedPdb");
            if (selectedDownload) selectedDownload.textContent = "Download selection";
            const screenshot = toolbar.querySelector("#viewerScreenshotPng");
            if (screenshot) screenshot.textContent = "Screenshot";

            const controls = document.createElement("div");
            controls.className = "construct-workbench-toolbar-controls";
            [
              ["View", ["#viewerReset", "#viewerEctodomain", "#viewerToggleSelectionOnly"]],
              ["Export", ["#viewerDownloadPdb", "#viewerDownloadSelectedPdb"]],
              ["Capture", ["#viewerScreenshotPng"]],
            ].forEach(([name, selectors]) => {
              const group = document.createElement("div");
              group.className = "construct-workbench-toolbar-group";
              group.setAttribute("aria-label", `${name} tools`);
              group.append(...selectors.map((selector) => toolbar.querySelector(selector)).filter(Boolean));
              controls.append(group);
            });
            toolbar.append(controls);

            const structureCanvas = viewerPanel.querySelector(".structure-viewer");
            const legend = toolbar.querySelector(".viewer-note");
            if (legend && structureCanvas) {
              legend.classList.add("construct-workbench-structure-legend");
              legend.textContent = "Color: AlphaFold pLDDT (B-factor)";
              structureCanvas.append(legend);
            }
            const status = viewerPanel.querySelector("#viewerInfo");
            if (status && structureCanvas) {
              status.classList.add("construct-workbench-structure-status");
              structureCanvas.append(status);
            }
          }

          function makeSequenceDetails(sequenceMetadata) {
            const details = document.createElement("details");
            details.className = "construct-workbench-sequence-details";
            const summary = document.createElement("summary");
            summary.textContent = `Open full canonical sequence · ${sequenceMetadata.targetName} · ${sequenceMetadata.length} aa`;
            details.append(summary, sequencePanel);
            return details;
          }

          function installWorkbenchPlddtRenderer(card) {
            const source = card.querySelector("#plddtPlot");
            if (!source) return;
            source.setAttribute("preserveAspectRatio", "none");
            source.setAttribute("shape-rendering", "geometricPrecision");
          }
        }
"""


def _render_structure_widget(
    entry: dict[str, Any],
    report: dict[str, Any],
    *,
    batch_dir: Path,
    page_dir: Path | None = None,
    allow_af3_structures: bool = False,
) -> str:
    local_vendor_asset = page_dir.parent / "3Dmol-min.js" if page_dir is not None else None
    if local_vendor_asset is not None and _is_valid_threedmol_asset(local_vendor_asset):
        threedmol_script = '<script src="../3Dmol-min.js"></script>'
    else:
        threedmol_script = (
            f'<script src="{_THREEDMOL_URL}" integrity="{_THREEDMOL_SRI}" '
            'crossorigin="anonymous"></script>'
        )
    structure_path = entry.get("portal_structure_path")
    if not structure_path and page_dir is not None:
        inferred_portal_dir = page_dir.parent
        inferred_structures_dir = inferred_portal_dir / "structures"
        inferred_structures_dir.mkdir(parents=True, exist_ok=True)
        structure_path = _copy_portal_structure(
            report_data=report,
            structures_dir=inferred_structures_dir,
            batch_dir=batch_dir,
            allow_af3_structures=allow_af3_structures,
        )
        if structure_path:
            entry["portal_structure_path"] = structure_path
    target = report.get("target", {})
    sequence = str(target.get("sequence") or "")
    ectodomain = report.get("ectodomain") or {}
    design_region_label = _design_region_label(entry, report)
    design_region_lower = design_region_label.lower()
    topology = report.get("topology") or {}
    topology_class = str(topology.get("topology_class") or "")
    plot_scope_full_length = topology_class.startswith("multipass")
    plot_scope_label = "full-length protein" if plot_scope_full_length else design_region_lower
    identity_label = _identity_label(report)
    if not structure_path or not sequence:
        return "<p>No local AlphaFold structure is available for interactive viewing.</p>"
    resolved_structure_path = _resolve_existing_path(str(structure_path), batch_dir)
    if resolved_structure_path is None:
        return "<p>No local AlphaFold structure is available for interactive viewing.</p>"
    if not _alphafold_pdb_matches_target_sequence(report, resolved_structure_path):
        return "<p>No matching canonical AlphaFold structure is available for interactive viewing.</p>"
    try:
        structure_text = resolved_structure_path.read_text(encoding="utf-8")
    except Exception:
        return "<p>No local AlphaFold structure is available for interactive viewing.</p>"
    viewer_payload = _load_viewer_payload(
        report=report,
        batch_dir=batch_dir,
        pdb_path=resolved_structure_path,
        allow_af3_structures=allow_af3_structures,
    )
    structure_source_label = _structure_source_label(resolved_structure_path)

    entry_name = str(target.get("entry_name") or entry.get("entry_name") or "target")
    sequence_rows = _render_sequence_rows(sequence, ectodomain=ectodomain)
    return f"""
    <div class="viewer-layout">
      <div class="viewer-panel">
        <div class="viewer-toolbar">
          <button type="button" class="mini-button" id="viewerReset">Reset view</button>
          <button type="button" class="mini-button" id="viewerEctodomain">Focus {escape(design_region_lower)}</button>
          <button type="button" class="mini-button" id="viewerToggleSelectionOnly" disabled>Hide unselected</button>
          <button type="button" class="mini-button" id="viewerDownloadPdb">Download PDB</button>
          <button type="button" class="mini-button" id="viewerDownloadSelectedPdb" disabled>Download selected PDB</button>
          <button type="button" class="mini-button" id="viewerScreenshotPng">Screenshot PNG</button>
          <span class="viewer-note">Colored by pLDDT from the {escape(structure_source_label)} B-factor field.</span>
        </div>
        <div id="structureViewer" class="structure-viewer"></div>
        <div id="viewerInfo" class="viewer-info">Click a residue in the sequence or structure to inspect it.</div>
      </div>
      <div class="sequence-panel">
        <div class="sequence-header">
          <strong>{escape(entry_name)}</strong>
          <span>{len(sequence)} aa total; {escape(design_region_lower)} {escape(_format_region(ectodomain) or 'n/a')}</span>
        </div>
        <div id="sequenceViewer" class="sequence-viewer">
          {sequence_rows}
        </div>
      </div>
    </div>
    <div class="interactive-plots">
      <div class="plot-card">
        <div class="plot-card-header">
          <h3>Interactive pLDDT</h3>
          <div class="plot-controls">
            <button type="button" class="mini-button" id="plddtZoomIn">Zoom in</button>
            <button type="button" class="mini-button" id="plddtZoomOut">Zoom out</button>
            <button type="button" class="mini-button" id="plddtZoomReset">Reset</button>
          </div>
        </div>
        <svg id="plddtPlot" class="plot-svg" viewBox="0 0 720 240" preserveAspectRatio="none"></svg>
        <div class="topology-legend" id="plotTopologyLegend"></div>
        <div id="plddtInfo" class="viewer-info">Viewing the {escape(plot_scope_label)} pLDDT trace. Use the mouse wheel or buttons to zoom.</div>
      </div>
      <div class="plot-card">
        <div class="plot-card-header">
          <h3>Interactive PAE</h3>
          <div class="plot-controls">
            <button type="button" class="mini-button" id="paeZoomIn">Zoom in</button>
            <button type="button" class="mini-button" id="paeZoomOut">Zoom out</button>
            <button type="button" class="mini-button" id="paeZoomReset">Reset</button>
          </div>
        </div>
        <canvas id="paePlot" class="pae-canvas" width="520" height="520"></canvas>
        <div class="topology-legend" id="paeTopologyLegend"></div>
        <div id="paeInfo" class="viewer-info">The {escape(plot_scope_label)} PAE block updates with the same selection. Use the mouse wheel or buttons to zoom.</div>
      </div>
    </div>
    <div class="selection-live-card">
      <div class="selection-live-header">
        <h3>Live Selected Region</h3>
        <div class="button-row">
          <button type="button" class="mini-button" id="copySelectionTsv">Copy TSV</button>
          <button type="button" class="mini-button" id="copySelectionFasta">Copy FASTA</button>
        </div>
      </div>
      <p class="selection-export-help"><a class="inline-link" href="../constructs.html#after-export" target="_blank" rel="noopener">After export <span class="muted">(opens in new tab)</span></a></p>
      <div id="selectionSummary" class="viewer-info">Select a residue or contiguous range to populate the table.</div>
      <div id="selectionWarnings" class="selection-warnings">
        <p class="viewer-info">Warnings for the selected region will appear here.</p>
      </div>
      <div class="selection-table-wrap">
        <table id="selectionTable">
          <thead>
            <tr><th>Name</th><th>Species</th><th>Boundary</th><th>Sequence</th><th>{escape(identity_label)}</th><th>Extracellular accessible identity</th><th>Links</th></tr>
          </thead>
          <tbody id="selectionTableBody">
            <tr><td colspan="7">No selection yet.</td></tr>
          </tbody>
        </table>
      </div>
      <textarea id="selectionTsv" class="selection-tsv" readonly></textarea>
      <textarea id="selectionFasta" class="selection-tsv" readonly></textarea>
    </div>
    {threedmol_script}
    <script>
      (() => {{
        const sequence = {_js_json(sequence)};
        const pdbText = {_js_json(structure_text)};
        const ectodomain = {{ start: {_js_json(ectodomain.get("start"))}, end: {_js_json(ectodomain.get("end"))} }};
        const designRegionLabel = {_js_json(design_region_label)};
        const designRegionLower = {_js_json(design_region_lower)};
        const plotScopeFullLength = {_js_json(plot_scope_full_length)};
        const targetInfo = {{
          accession: {_js_json(target.get("accession") or "")},
          entryName: {_js_json(target.get("entry_name") or "")},
          geneSymbol: {_js_json(target.get("gene_symbol") or target.get("entry_name") or "TARGET")},
          nameSymbol: {_js_json(_entry_name_base(str(target.get("entry_name") or "")) or target.get("gene_symbol") or "TARGET")},
          species: {_js_json(_target_species_key(target))},
          identityHeader: {_js_json(_identity_header(identity_label))},
          surfaceIdentityHeader: {_js_json(_surface_identity_header(identity_label))}
        }};
        const viewerData = {_js_json(viewer_payload)};
        const viewerEl = document.getElementById('structureViewer');
        const sequenceEl = document.getElementById('sequenceViewer');
        const infoEl = document.getElementById('viewerInfo');
        const plddtSvg = document.getElementById('plddtPlot');
        const plddtInfo = document.getElementById('plddtInfo');
        const paeCanvas = document.getElementById('paePlot');
        const paeInfo = document.getElementById('paeInfo');
        const plotTopologyLegend = document.getElementById('plotTopologyLegend');
        const paeTopologyLegend = document.getElementById('paeTopologyLegend');
        const selectionSummary = document.getElementById('selectionSummary');
        const selectionWarnings = document.getElementById('selectionWarnings');
        const selectionTableBody = document.getElementById('selectionTableBody');
        const selectionTsv = document.getElementById('selectionTsv');
        const selectionFasta = document.getElementById('selectionFasta');
        const copySelectionTsv = document.getElementById('copySelectionTsv');
        const copySelectionFasta = document.getElementById('copySelectionFasta');
        const plddtZoomIn = document.getElementById('plddtZoomIn');
        const plddtZoomOut = document.getElementById('plddtZoomOut');
        const plddtZoomReset = document.getElementById('plddtZoomReset');
        const paeZoomIn = document.getElementById('paeZoomIn');
        const paeZoomOut = document.getElementById('paeZoomOut');
        const paeZoomReset = document.getElementById('paeZoomReset');
        const resetButton = document.getElementById('viewerReset');
        const ectoButton = document.getElementById('viewerEctodomain');
        const selectionOnlyButton = document.getElementById('viewerToggleSelectionOnly');
        const downloadPdbButton = document.getElementById('viewerDownloadPdb');
        const downloadSelectedPdbButton = document.getElementById('viewerDownloadSelectedPdb');
        const screenshotPngButton = document.getElementById('viewerScreenshotPng');
{_construct_workbench_js()}
        const residueSpans = Array.from(sequenceEl.querySelectorAll('.residue'));
        const residueMap = new Map(residueSpans.map((el) => [Number(el.dataset.resi), el]));
        const sequenceToStructureResidues = Array.isArray(viewerData.sequenceToStructureResidues) ? viewerData.sequenceToStructureResidues : [];
        const structureToSequenceResidues = viewerData.structureToSequenceResidues || {{}};
        let viewer = null;
        let viewerInitStarted = false;
        let selectedResidue = null;
        let selectedRange = null;
        let selectedResidueSet = null;
        let anchorResidue = null;
        let dragStartResidue = null;
        let isDragging = false;
        let plddtDragStart = null;
        let plddtDragMoved = false;
        let paeDragStart = null;
        let paeDragMoved = false;
        let plotResize = null;
        let showSelectionOnly = false;
        const plddtLayout = {{ width: 720, height: 240, margin: {{ top: 16, right: 20, bottom: 30, left: 42 }} }};
        const plotEdgeThresholdPx = 10;
        const fullSequenceLength = Math.max(sequence.length, (viewerData.plddt || []).length || 0, 1);
        const plotDomain = plotScopeFullLength
          ? {{ start: 1, end: fullSequenceLength }}
          : (ectodomain.start && ectodomain.end)
          ? {{ start: Number(ectodomain.start), end: Number(ectodomain.end) }}
          : {{ start: 1, end: fullSequenceLength }};
        let plddtView = {{ start: plotDomain.start, end: plotDomain.end }};
        let paeView = {{ start: plotDomain.start, end: plotDomain.end }};
        const minPlddtWindow = 12;
        const minPaeWindow = 12;
        const selectionMutationPositions = new Set();
        const selectionFurinMutationSites = new Set();
        const homologPositionLookupCache = new Map();

        function alphaFoldColor(atom) {{
          const score = Number(atom.b) || 0;
          if (score >= 90) return '#0053d6';
          if (score >= 70) return '#65cbf3';
          if (score >= 50) return '#ffdb13';
          return '#ff7d45';
        }}

        function plddtColor(score) {{
          const value = Number(score) || 0;
          if (value >= 90) return '#0053d6';
          if (value >= 70) return '#65cbf3';
          if (value >= 50) return '#ffdb13';
          return '#ff7d45';
        }}

        function normalizeTopologyLocation(value) {{
          const text = String(value || '').toLowerCase();
          if (text.includes('membrane') || text.includes('transmembrane')) return 'membrane';
          if (text.includes('intra') || text.includes('cyto')) return 'intracellular';
          if (text.includes('extra') || text.includes('lumen')) return 'extracellular';
          return 'unknown';
        }}

        function topologyDisplayName(location) {{
          if (location === 'extracellular') return 'Extracellular';
          if (location === 'membrane') return 'Membrane';
          if (location === 'intracellular') return 'Intracellular';
          return 'Unknown';
        }}

        function topologyColor(location, alpha = 1) {{
          const colors = {{
            extracellular: [18, 145, 126],
            membrane: [240, 158, 38],
            intracellular: [105, 82, 168],
            unknown: [150, 150, 150],
          }};
          const rgb = colors[location] || colors.unknown;
          return `rgba(${{rgb[0]}}, ${{rgb[1]}}, ${{rgb[2]}}, ${{alpha}})`;
        }}

        const topologyByPosition = new Map((viewerData.residueAnnotations || []).map((item) => [
          Number(item.position),
          normalizeTopologyLocation(item.topology_location),
        ]));
        const gpcrSegments = ((viewerData.gpcrAnnotation || {{}}).segments || []).map((segment) => ({{
          name: String(segment.name || ''),
          start: Number(segment.start),
          end: Number(segment.end),
          genericStart: segment.generic_start || '',
          genericEnd: segment.generic_end || '',
        }})).filter((segment) => Number.isFinite(segment.start) && Number.isFinite(segment.end));

        function topologyAtPosition(position) {{
          return topologyByPosition.get(Number(position)) || 'unknown';
        }}

        function topologyIntervals(start, end) {{
          const intervals = [];
          let currentLocation = null;
          let currentStart = null;
          for (let residue = Number(start); residue <= Number(end); residue += 1) {{
            const location = topologyAtPosition(residue);
            if (currentLocation === null) {{
              currentLocation = location;
              currentStart = residue;
            }} else if (location !== currentLocation) {{
              intervals.push({{ start: currentStart, end: residue - 1, location: currentLocation }});
              currentLocation = location;
              currentStart = residue;
            }}
          }}
          if (currentLocation !== null) intervals.push({{ start: currentStart, end: Number(end), location: currentLocation }});
          return intervals;
        }}

        function renderTopologyLegend(targetEl) {{
          if (!targetEl) return;
          const topologyItems = ['extracellular', 'membrane', 'intracellular', 'unknown'].map((location) => `
            <span class="topology-legend-item">
              <span class="topology-swatch" style="background:${{topologyColor(location, 0.95)}}"></span>${{topologyDisplayName(location)}}
            </span>
          `).join('');
          const gpcrItem = gpcrSegments.length
            ? '<span class="topology-legend-item"><span class="topology-swatch gpcr-swatch"></span>GPCRdb segments</span>'
            : '';
          targetEl.innerHTML = topologyItems + gpcrItem;
        }}

        function applyCysteineStyle(targetViewer, selection = {{}}) {{
          targetViewer.addStyle({{ ...selection, resn: 'CYS' }}, {{
            stick: {{ color: 'yellow', radius: 0.3 }},
          }});
        }}

        function makeRange(start, end) {{
          const values = [];
          const safeStart = Number(start) || 0;
          const safeEnd = Number(end) || 0;
          for (let i = safeStart; i <= safeEnd; i += 1) values.push(i);
          return values;
        }}

        function mapSequenceResidueToStructure(residue) {{
          const residueIndex = Number(residue);
          const mapped = Number(sequenceToStructureResidues[residueIndex - 1]);
          if (Number.isFinite(mapped) && mapped > 0) return mapped;
          return residueIndex;
        }}

        function mapStructureResidueToSequence(residue) {{
          const residueIndex = Number(residue);
          const mapped = Number(structureToSequenceResidues[String(residueIndex)] ?? structureToSequenceResidues[residueIndex]);
          if (Number.isFinite(mapped) && mapped > 0) return mapped;
          return residueIndex;
        }}

        function structureResiduesForSequenceRange(start, end) {{
          const residues = [];
          const seen = new Set();
          for (let residue = Math.min(Number(start), Number(end)); residue <= Math.max(Number(start), Number(end)); residue += 1) {{
            const mapped = mapSequenceResidueToStructure(residue);
            if (!Number.isFinite(mapped) || mapped <= 0 || seen.has(mapped)) continue;
            seen.add(mapped);
            residues.push(mapped);
          }}
          return residues;
        }}

        function clampToRange(value, minValue, maxValue) {{
          return Math.max(minValue, Math.min(maxValue, Number(value) || minValue));
        }}

        function viewLength(view) {{
          return Math.max(1, Number(view.end) - Number(view.start) + 1);
        }}

        function normalizeView(start, end, domainStart, domainEnd, minWidth) {{
          const safeDomainStart = Number(domainStart);
          const safeDomainEnd = Number(domainEnd);
          const domainWidth = Math.max(1, safeDomainEnd - safeDomainStart + 1);
          const widthFloor = Math.min(domainWidth, Math.max(1, Number(minWidth) || 1));
          const orderedStart = Math.min(Number(start), Number(end));
          const orderedEnd = Math.max(Number(start), Number(end));
          let width = Math.max(widthFloor, orderedEnd - orderedStart + 1);
          width = Math.min(width, domainWidth);
          let nextStart = clampToRange(orderedStart, safeDomainStart, safeDomainEnd - width + 1);
          let nextEnd = nextStart + width - 1;
          if (nextEnd > safeDomainEnd) {{
            nextEnd = safeDomainEnd;
            nextStart = Math.max(safeDomainStart, nextEnd - width + 1);
          }}
          return {{ start: nextStart, end: nextEnd }};
        }}

        function selectionCenter() {{
          const bounds = selectionBounds();
          if (bounds) {{
            return Math.round((bounds.start + bounds.end) / 2);
          }}
          return Math.round((plotDomain.start + plotDomain.end) / 2);
        }}

        function setPlotView(source, start, end) {{
          if (source === 'plddt') {{
            plddtView = normalizeView(start, end, plotDomain.start, plotDomain.end, minPlddtWindow);
          }} else if (source === 'pae' && paeView) {{
            paeView = normalizeView(start, end, plotDomain.start, plotDomain.end, minPaeWindow);
          }}
          renderInteractivePlots();
        }}

        function resetPlotViews() {{
          plddtView = {{ start: plotDomain.start, end: plotDomain.end }};
          if (paeView) {{
            paeView = {{ start: plotDomain.start, end: plotDomain.end }};
          }}
        }}

        function zoomPlot(source, factor, centerResidue) {{
          const current = source === 'plddt' ? plddtView : paeView;
          if (!current) return;
          const minWidth = source === 'plddt' ? minPlddtWindow : minPaeWindow;
          const domainWidth = plotDomain.end - plotDomain.start + 1;
          const currentWidth = viewLength(current);
          const nextWidth = Math.max(minWidth, Math.min(domainWidth, Math.round(currentWidth * factor)));
          const center = clampToRange(centerResidue ?? selectionCenter(), plotDomain.start, plotDomain.end);
          const nextStart = Math.round(center - ((nextWidth - 1) / 2));
          setPlotView(source, nextStart, nextStart + nextWidth - 1);
        }}

        function ensureRangeVisibleInView(source, start, end) {{
          const current = source === 'plddt' ? plddtView : paeView;
          if (!current) return;
          const clippedStart = clampToRange(start, plotDomain.start, plotDomain.end);
          const clippedEnd = clampToRange(end, plotDomain.start, plotDomain.end);
          if (clippedStart >= current.start && clippedEnd <= current.end) return;
          const currentWidth = viewLength(current);
          if ((clippedEnd - clippedStart + 1) >= currentWidth) {{
            setPlotView(source, clippedStart, clippedEnd);
            return;
          }}
          let nextStart = current.start;
          let nextEnd = current.end;
          if (clippedStart < current.start) {{
            nextStart = clippedStart;
            nextEnd = clippedStart + currentWidth - 1;
          }} else if (clippedEnd > current.end) {{
            nextEnd = clippedEnd;
            nextStart = clippedEnd - currentWidth + 1;
          }}
          setPlotView(source, nextStart, nextEnd);
        }}

        function ensureSelectionVisible(start, end) {{
          if (end >= plotDomain.start && start <= plotDomain.end) {{
            ensureRangeVisibleInView('plddt', start, end);
            if (paeView) ensureRangeVisibleInView('pae', start, end);
          }}
        }}

        function clearSequenceSelection() {{
          residueSpans.forEach((el) => el.classList.remove('active-residue', 'active-range'));
        }}

        function updateSequenceSelection() {{
          clearSequenceSelection();
          if (selectedResidueSet && selectedResidueSet.positions.length) {{
            selectedResidueSet.positions.forEach((resi) => {{
              const span = residueMap.get(resi);
              if (span) span.classList.add('active-range');
            }});
          }}
          if (selectedRange) {{
            const range = makeRange(selectedRange.start, selectedRange.end);
            range.forEach((resi) => {{
              const span = residueMap.get(resi);
              if (span) span.classList.add('active-range');
            }});
          }}
          if (selectedResidue !== null) {{
            const span = residueMap.get(selectedResidue);
            if (span) {{
              span.classList.add('active-residue');
              span.scrollIntoView({{ block: 'nearest', inline: 'center', behavior: 'smooth' }});
            }}
          }}
        }}

        function applyBaseStyle(targetViewer) {{
          targetViewer.setStyle({{}}, {{ cartoon: {{ colorfunc: alphaFoldColor }} }});
          applyCysteineStyle(targetViewer);
        }}

        function selectedStructureResidues() {{
          if (selectedResidueSet && selectedResidueSet.positions.length) {{
            return selectedResidueSet.positions
              .map((position) => mapSequenceResidueToStructure(position))
              .filter((position) => Number.isFinite(Number(position)));
          }}
          if (selectedRange) {{
            return structureResiduesForSequenceRange(selectedRange.start, selectedRange.end);
          }}
          if (selectedResidue !== null) {{
            const structureResidue = mapSequenceResidueToStructure(selectedResidue);
            return Number.isFinite(Number(structureResidue)) ? [structureResidue] : [];
          }}
          return [];
        }}

        function selectedSequenceResidues() {{
          if (selectedResidueSet && selectedResidueSet.positions.length) {{
            return selectedResidueSet.positions
              .map((position) => Number(position))
              .filter((position) => Number.isFinite(position) && position > 0);
          }}
          if (selectedRange) return makeRange(selectedRange.start, selectedRange.end);
          if (selectedResidue !== null) return [selectedResidue];
          return [];
        }}

        function safeFilenamePart(value) {{
          return String(value || 'target')
            .trim()
            .replace(/[^A-Za-z0-9_.-]+/g, '_')
            .replace(/^_+|_+$/g, '')
            .toLowerCase() || 'target';
        }}

        function downloadTextFile(filename, text, mimeType = 'chemical/x-pdb') {{
          const blob = new Blob([text], {{ type: `${{mimeType}};charset=utf-8` }});
          const url = URL.createObjectURL(blob);
          const link = document.createElement('a');
          link.href = url;
          link.download = filename;
          document.body.appendChild(link);
          link.click();
          link.remove();
          URL.revokeObjectURL(url);
        }}

        function downloadDataUri(filename, dataUri) {{
          const link = document.createElement('a');
          link.href = dataUri;
          link.download = filename;
          document.body.appendChild(link);
          link.click();
          link.remove();
        }}

        function pdbResidueNumber(line) {{
          const fixed = Number.parseInt(line.slice(22, 26).trim(), 10);
          if (Number.isFinite(fixed)) return fixed;
          const fields = line.trim().split(/\\s+/);
          const fallback = Number.parseInt(fields[5], 10);
          return Number.isFinite(fallback) ? fallback : null;
        }}

        function pdbForSequenceResidues(sequenceResidues, label) {{
          const wanted = new Set(sequenceResidues.map((position) => Number(position)));
          const atomLines = [];
          pdbText.split(/\\r?\\n/).forEach((line) => {{
            if (!/^(ATOM  |HETATM)/.test(line)) return;
            const structureResidue = pdbResidueNumber(line);
            if (structureResidue === null) return;
            const sequenceResidue = mapStructureResidueToSequence(structureResidue);
            if (wanted.has(sequenceResidue)) atomLines.push(line);
          }});
          const remarks = [
            'REMARK OpenAntigens selected-region PDB export',
            `REMARK Target: ${{targetInfo.entryName || targetInfo.geneSymbol || 'target'}}`,
            `REMARK Selection: ${{label || 'selected residues'}}`,
            `REMARK Residues exported: ${{atomLines.length ? sequenceResidues.length : 0}} sequence position(s)`,
          ];
          if (!atomLines.length) {{
            remarks.push('REMARK No ATOM/HETATM records matched the current selection.');
          }}
          return remarks.concat(atomLines, ['TER', 'END']).join('\\n') + '\\n';
        }}

        function downloadFullPdb() {{
          downloadTextFile(`${{safeFilenamePart(targetInfo.entryName || targetInfo.geneSymbol)}}_alphafold.pdb`, pdbText.endsWith('\\n') ? pdbText : `${{pdbText}}\\n`);
          window.trackOpenAntigensEvent('Download PDB Full');
          setInfo('Full model PDB downloaded.');
        }}

        function downloadSelectedPdb() {{
          const bounds = selectionBounds();
          const residues = selectedSequenceResidues();
          if (!bounds || !residues.length) {{
            setInfo('Select a residue or region before downloading a selected-region PDB.');
            return;
          }}
          const suffix = bounds.nonContiguous ? 'selected_residues' : `${{bounds.start}}-${{bounds.end}}`;
          const filename = `${{safeFilenamePart(targetInfo.entryName || targetInfo.geneSymbol)}}_${{safeFilenamePart(suffix)}}.pdb`;
          downloadTextFile(filename, pdbForSequenceResidues(residues, bounds.label || suffix));
          window.trackOpenAntigensEvent('Download PDB Selection');
          setInfo(`Selected-region PDB downloaded for ${{bounds.label || suffix}}.`);
        }}

        function downloadViewerScreenshot() {{
          if (!viewer) {{
            setInfo('3D viewer is not ready; screenshot export is unavailable.');
            return;
          }}
          if (typeof viewer.pngURI !== 'function') {{
            setInfo('Screenshot export is unavailable because this 3Dmol build does not expose PNG export.');
            return;
          }}
          const originalWidth = Math.max(1, Math.round(viewerEl.clientWidth || (typeof viewer.getWidth === 'function' ? viewer.getWidth() : 0) || 900));
          const originalHeight = Math.max(1, Math.round(viewerEl.clientHeight || (typeof viewer.getHeight === 'function' ? viewer.getHeight() : 0) || 520));
          const exportScale = Math.max(1, Math.min(3, 3000 / originalWidth, 3000 / originalHeight));
          const exportWidth = Math.round(originalWidth * exportScale);
          const exportHeight = Math.round(originalHeight * exportScale);
          try {{
            if (typeof viewer.setWidth === 'function' && typeof viewer.setHeight === 'function') {{
              viewer.setWidth(exportWidth);
              viewer.setHeight(exportHeight);
            }}
            viewer.render();
            const pngUri = viewer.pngURI();
            if (!pngUri || !pngUri.startsWith('data:image/png')) {{
              throw new Error('3Dmol did not return PNG image data.');
            }}
            const filename = `${{safeFilenamePart(targetInfo.entryName || targetInfo.geneSymbol)}}_3d_view_${{exportWidth}}x${{exportHeight}}.png`;
            downloadDataUri(filename, pngUri);
            window.trackOpenAntigensEvent('Download PNG');
            setInfo(`3D viewer screenshot downloaded as ${{exportWidth}} x ${{exportHeight}} PNG.`);
          }} catch (error) {{
            setInfo(`Screenshot export failed: ${{error}}`);
          }} finally {{
            if (typeof viewer.setWidth === 'function' && typeof viewer.setHeight === 'function') {{
              viewer.setWidth(originalWidth);
              viewer.setHeight(originalHeight);
              viewer.render();
            }}
          }}
        }}

        function updateSelectionOnlyControl() {{
          if (!selectionOnlyButton) return;
          const hasSelection = Boolean(selectedRange || selectedResidue !== null || (selectedResidueSet && selectedResidueSet.positions.length));
          selectionOnlyButton.disabled = !hasSelection;
          selectionOnlyButton.textContent = showSelectionOnly ? 'Show all residues' : 'Hide unselected';
          selectionOnlyButton.classList.toggle('active', showSelectionOnly && hasSelection);
          if (downloadSelectedPdbButton) downloadSelectedPdbButton.disabled = !hasSelection;
        }}

        let selectionRenderPending = false;
        function renderSelection() {{
          updateSelectionOnlyControl();
          renderSelectionTable();
          if (selectionRenderPending) return;
          selectionRenderPending = true;
          requestAnimationFrame(() => {{
            selectionRenderPending = false;
            renderSelectionViews();
          }});
        }}

        function renderSelectionViews() {{
          if (viewer) {{
            viewer.setStyle({{}}, {{}});
            const isolatedResidues = selectedStructureResidues();
            if (showSelectionOnly && isolatedResidues.length) {{
              viewer.setStyle({{ resi: isolatedResidues }}, {{ cartoon: {{ colorfunc: alphaFoldColor }} }});
              applyCysteineStyle(viewer, {{ resi: isolatedResidues }});
            }} else {{
              applyBaseStyle(viewer);
            }}
            if (selectedRange) {{
              const structureResidues = structureResiduesForSequenceRange(selectedRange.start, selectedRange.end);
              viewer.addStyle({{ resi: structureResidues }}, {{
                cartoon: {{ color: '#4c956c', opacity: 0.95 }},
                stick: {{ color: '#4c956c', radius: 0.18 }}
              }});
            }}
            if (selectedResidueSet && selectedResidueSet.positions.length) {{
              const structureResidues = selectedStructureResidues();
              viewer.addStyle({{ resi: structureResidues }}, {{
                cartoon: {{ color: '#4c956c', opacity: 0.95 }},
                stick: {{ color: '#4c956c', radius: 0.18 }}
              }});
            }}
            if (selectedResidue !== null) {{
              viewer.addStyle({{ resi: [mapSequenceResidueToStructure(selectedResidue)] }}, {{
                stick: {{ color: '#c1121f', radius: 0.24 }},
                sphere: {{ color: '#c1121f', radius: 0.42 }}
              }});
            }}
            if (showSelectionOnly && isolatedResidues.length) {{
              applyCysteineStyle(viewer, {{ resi: isolatedResidues }});
            }} else {{
              applyCysteineStyle(viewer);
            }}
            viewer.render();
          }}
          updateSequenceSelection();
          renderInteractivePlots();
        }}

        function setInfo(message) {{
          infoEl.textContent = message;
        }}

        function focusRange(start, end, label) {{
          const rangeStart = Math.min(Number(start), Number(end));
          const rangeEnd = Math.max(Number(start), Number(end));
          const structureResidues = structureResiduesForSequenceRange(rangeStart, rangeEnd);
          selectedRange = {{ start: rangeStart, end: rangeEnd, label: label || 'Selected range' }};
          selectedResidue = null;
          selectedResidueSet = null;
          ensureSelectionVisible(rangeStart, rangeEnd);
          if (viewer) viewer.zoomTo({{ resi: structureResidues }});
          renderSelection();
          setInfo(`${{selectedRange.label}}: residues ${{selectedRange.start}}-${{selectedRange.end}}`);
        }}

        function selectResidue(resi) {{
          const residueIndex = Number(resi);
          const structureResidue = mapSequenceResidueToStructure(residueIndex);
          selectedResidue = residueIndex;
          selectedRange = null;
          selectedResidueSet = null;
          ensureSelectionVisible(residueIndex, residueIndex);
          if (viewer) viewer.zoomTo({{ resi: [structureResidue] }});
          renderSelection();
          setInfo(`Residue ${{residueIndex}}: ${{sequence[residueIndex - 1] || '?'}}`);
        }}

        function setRangeFromSelection(start, end, label) {{
          focusRange(start, end, label || 'Selected range');
          anchorResidue = Number(start);
        }}

        function extracellularTopologyPositions() {{
          return (viewerData.residueAnnotations || [])
            .filter((item) => normalizeTopologyLocation(item.topology_location) === 'extracellular')
            .map((item) => Number(item.position))
            .filter((position) => Number.isFinite(position))
            .filter((position, index, array) => array.indexOf(position) === index)
            .sort((a, b) => a - b);
        }}

        function focusResidueSet(positions, label) {{
          if (!viewer || !positions.length) return false;
          const sorted = Array.from(new Set(positions.map(Number).filter(Number.isFinite))).sort((a, b) => a - b);
          if (!sorted.length) return false;
          selectedResidueSet = {{ positions: sorted, label: label || 'Selected residues' }};
          selectedRange = null;
          selectedResidue = null;
          const minResidue = sorted[0];
          const maxResidue = sorted[sorted.length - 1];
          ensureSelectionVisible(minResidue, maxResidue);
          renderSelection();
          const structureResidues = selectedStructureResidues();
          if (structureResidues.length) viewer.zoomTo({{ resi: structureResidues }});
          viewer.render();
          setInfo(`${{selectedResidueSet.label}}: ${{sorted.length}} non-contiguous residue(s) focused.`);
          return true;
        }}

        function speciesCode(species) {{
          const normalized = String(species || '').toLowerCase();
          if (normalized === 'human') return 'HUMAN';
          if (normalized === 'mouse') return 'MOUSE';
          if (normalized === 'macaca_fascicularis') return 'MACFA';
          return String(species || 'SPECIES').toUpperCase().replace(/[^A-Z0-9]+/g, '_');
        }}

        function targetSpeciesKey() {{
          return String(targetInfo.species || 'human').toLowerCase();
        }}

        function proteinResourceUrl(accession, entryName) {{
          const accessionText = String(accession || '').trim();
          const entryText = String(entryName || '').trim();
          if (accessionText.startsWith('NP_') || accessionText.startsWith('XP_') || accessionText.startsWith('YP_') || accessionText.startsWith('WP_')) {{
            return `https://www.ncbi.nlm.nih.gov/protein/${{encodeURIComponent(accessionText)}}`;
          }}
          if (accessionText || entryText) {{
            return `https://www.uniprot.org/uniprotkb/${{encodeURIComponent(accessionText || entryText)}}`;
          }}
          return '';
        }}

        function proteinResourceLabel(accession) {{
          const accessionText = String(accession || '').trim();
          return accessionText.startsWith('NP_') || accessionText.startsWith('XP_') || accessionText.startsWith('YP_') || accessionText.startsWith('WP_')
            ? 'RefSeq'
            : 'UniProt';
        }}

        function uniprotResourceUrl(accession, entryName, geneSymbol, species) {{
          const accessionText = String(accession || '').trim();
          const entryText = String(entryName || '').trim();
          if (accessionText && !accessionText.startsWith('NP_') && !accessionText.startsWith('XP_') && !accessionText.startsWith('YP_') && !accessionText.startsWith('WP_')) {{
            return `https://www.uniprot.org/uniprotkb/${{encodeURIComponent(accessionText)}}`;
          }}
          if (entryText.includes('_') && !accessionText.startsWith('NP_') && !accessionText.startsWith('XP_') && !accessionText.startsWith('YP_') && !accessionText.startsWith('WP_')) {{
            return `https://www.uniprot.org/uniprotkb/${{encodeURIComponent(entryText)}}`;
          }}
          const symbol = String(geneSymbol || '').trim();
          const normalized = String(species || '').toLowerCase();
          const taxonMap = {{ human: '9606', mouse: '10090', macaca_fascicularis: '9541' }};
          const taxonId = taxonMap[normalized];
          if (!symbol || !taxonId) return '';
          const query = `(gene_exact:${{symbol}}) AND (organism_id:${{taxonId}})`;
          return `https://www.uniprot.org/uniprotkb?query=${{encodeURIComponent(query)}}`;
        }}

        function geneResourceUrl(geneSymbol, species) {{
          const symbol = String(geneSymbol || '').trim();
          if (!symbol) return '';
          const normalized = String(species || '').toLowerCase();
          if (normalized === 'human') {{
            return `https://www.genenames.org/tools/search/#!/?query=${{encodeURIComponent(symbol)}}`;
          }}
          const taxonMap = {{ mouse: '10090', macaca_fascicularis: '9541' }};
          const taxonId = taxonMap[normalized];
          if (!taxonId) return '';
          return `https://www.ncbi.nlm.nih.gov/gene/?term=${{encodeURIComponent(`${{symbol}}[Gene Name] AND ${{taxonId}}[Taxonomy ID]`)}}`;
        }}

        function geneResourceLabel(species) {{
          return String(species || '').toLowerCase() === 'human' ? 'HGNC' : 'Gene';
        }}

        function resourceLinksHtml(resource) {{
          const links = [];
          const proteinUrl = proteinResourceUrl(resource.accession, resource.entryName);
          if (proteinUrl) {{
            links.push(`<a class="inline-link" href="${{proteinUrl}}" target="_blank" rel="noreferrer">${{proteinResourceLabel(resource.accession)}}</a>`);
          }}
          const uniprotUrl = uniprotResourceUrl(resource.accession, resource.entryName, resource.geneSymbol, resource.species);
          if (uniprotUrl && (!proteinUrl || !proteinUrl.includes('uniprot.org/uniprotkb/'))) {{
            links.push(`<a class="inline-link" href="${{uniprotUrl}}" target="_blank" rel="noreferrer">UniProt</a>`);
          }}
          const geneUrl = geneResourceUrl(resource.geneSymbol, resource.species);
          if (geneUrl) {{
            links.push(`<a class="inline-link" href="${{geneUrl}}" target="_blank" rel="noreferrer">${{geneResourceLabel(resource.species)}}</a>`);
          }}
          return links.join(' ') || '<span class="muted">n/a</span>';
        }}

        function clampResidue(residue) {{
          return Math.max(1, Math.min(sequence.length, Number(residue) || 1));
        }}

        function residueFromPlddtPointer(event) {{
          const rect = plddtSvg.getBoundingClientRect();
          const margin = plddtLayout.margin;
          const plotWidth = plddtLayout.width - margin.left - margin.right;
          const relativeX = ((event.clientX - rect.left) / Math.max(rect.width, 1)) * plddtLayout.width;
          const clampedX = Math.max(margin.left, Math.min(margin.left + plotWidth, relativeX));
          const fraction = (clampedX - margin.left) / Math.max(plotWidth, 1);
          const residue = plddtView.start + Math.round(fraction * Math.max(viewLength(plddtView) - 1, 1));
          return clampToRange(residue, plddtView.start, plddtView.end);
        }}

        function residueFromPaePointer(event) {{
          if (!paeView) return null;
          const rect = paeCanvas.getBoundingClientRect();
          const x = Math.max(0, Math.min(1, (event.clientX - rect.left) / Math.max(rect.width, 1)));
          const y = Math.max(0, Math.min(1, (event.clientY - rect.top) / Math.max(rect.height, 1)));
          const fraction = (x + y) / 2;
          const residue = paeView.start + Math.round(fraction * Math.max(viewLength(paeView) - 1, 1));
          return clampToRange(residue, paeView.start, paeView.end);
        }}

        function plddtEdgeAtEvent(event) {{
          if (!selectedRange) return null;
          if (selectedRange.end < plddtView.start || selectedRange.start > plddtView.end) return null;
          const rect = plddtSvg.getBoundingClientRect();
          const margin = plddtLayout.margin;
          const plotWidth = plddtLayout.width - margin.left - margin.right;
          const relativeX = ((event.clientX - rect.left) / Math.max(rect.width, 1)) * plddtLayout.width;
          const leftX = margin.left + ((selectedRange.start - plddtView.start) / Math.max(viewLength(plddtView) - 1, 1)) * plotWidth;
          const rightX = margin.left + ((selectedRange.end - plddtView.start) / Math.max(viewLength(plddtView) - 1, 1)) * plotWidth;
          if (Math.abs(relativeX - leftX) <= plotEdgeThresholdPx) return 'left';
          if (Math.abs(relativeX - rightX) <= plotEdgeThresholdPx) return 'right';
          return null;
        }}

        function paeEdgeAtEvent(event) {{
          if (!selectedRange || !paeView) return null;
          if (selectedRange.end < paeView.start || selectedRange.start > paeView.end) return null;
          const rect = paeCanvas.getBoundingClientRect();
          const paeLen = viewLength(paeView);
          const left = ((selectedRange.start - paeView.start) / Math.max(paeLen, 1)) * rect.width;
          const right = ((selectedRange.end - paeView.start + 1) / Math.max(paeLen, 1)) * rect.width;
          const x = event.clientX - rect.left;
          const y = event.clientY - rect.top;
          const nearLeft = (Math.abs(x - left) <= plotEdgeThresholdPx || Math.abs(y - left) <= plotEdgeThresholdPx)
            && x <= right + plotEdgeThresholdPx
            && y <= right + plotEdgeThresholdPx;
          const nearRight = (Math.abs(x - right) <= plotEdgeThresholdPx || Math.abs(y - right) <= plotEdgeThresholdPx)
            && x >= left - plotEdgeThresholdPx
            && y >= left - plotEdgeThresholdPx;
          if (nearLeft) return 'left';
          if (nearRight) return 'right';
          return null;
        }}

        function updatePlotCursor(target, edge) {{
          target.style.cursor = edge ? 'ew-resize' : 'crosshair';
        }}

        function paeMatrixIndexForResidue(residue, size) {{
          if (size <= 1) return 0;
          const domainLength = plotDomain.end - plotDomain.start + 1;
          const clipped = clampToRange(residue, plotDomain.start, plotDomain.end);
          const fraction = domainLength > 1 ? (clipped - plotDomain.start) / (domainLength - 1) : 0;
          return Math.max(0, Math.min(size - 1, Math.round(fraction * (size - 1))));
        }}

        function plasmaColor(t) {{
          const stops = [
            [0.0, [13, 8, 135]],
            [0.13, [75, 3, 161]],
            [0.25, [125, 3, 168]],
            [0.38, [168, 34, 150]],
            [0.5, [203, 70, 121]],
            [0.63, [229, 107, 93]],
            [0.75, [248, 148, 65]],
            [0.88, [253, 195, 40]],
            [1.0, [240, 249, 33]],
          ];
          const clamped = Math.max(0, Math.min(1, t));
          for (let i = 0; i < stops.length - 1; i += 1) {{
            const [leftStop, leftColor] = stops[i];
            const [rightStop, rightColor] = stops[i + 1];
            if (clamped <= rightStop) {{
              const local = (clamped - leftStop) / Math.max(rightStop - leftStop, 1e-6);
              return leftColor.map((value, index) => Math.round(value + (rightColor[index] - value) * local));
            }}
          }}
          return stops[stops.length - 1][1];
        }}

        function renderInteractivePlots() {{
          renderPlddtPlot();
          renderPaePlot();
        }}

        function selectionBounds() {{
          if (selectedResidueSet && selectedResidueSet.positions.length) {{
            return {{
              start: selectedResidueSet.positions[0],
              end: selectedResidueSet.positions[selectedResidueSet.positions.length - 1],
              label: selectedResidueSet.label || 'Selected residues',
              nonContiguous: true,
              count: selectedResidueSet.positions.length,
            }};
          }}
          if (selectedRange) {{
            return {{ start: selectedRange.start, end: selectedRange.end, label: selectedRange.label || 'Selected range' }};
          }}
          if (selectedResidue !== null) {{
            return {{ start: selectedResidue, end: selectedResidue, label: `Residue ${{selectedResidue}}` }};
          }}
          return null;
        }}

        function plddtStats(start, end) {{
          const plddt = viewerData.plddt || [];
          const values = [];
          for (let residue = Math.min(Number(start), Number(end)); residue <= Math.max(Number(start), Number(end)); residue += 1) {{
            const value = Number(plddt[residue - 1]);
            if (Number.isFinite(value)) values.push(value);
          }}
          if (!values.length) return null;
          const min = Math.min(...values);
          const max = Math.max(...values);
          const mean = values.reduce((total, value) => total + value, 0) / values.length;
          return {{ min, max, mean }};
        }}

        function extracellularSurfacePositionsInRange(start, end) {{
          return (viewerData.residueAnnotations || [])
            .filter((item) => Number(item.position) >= start && Number(item.position) <= end)
            .filter((item) => String(item.topology_location || '').toLowerCase() === 'extracellular' && item.surface_exposed === true)
            .map((item) => Number(item.position))
            .filter((position) => Number.isFinite(position));
        }}

        function surfaceIdentityText(identity, alignedPositions) {{
          if (identity === null || identity === undefined || !Number.isFinite(Number(identity))) return '';
          const formatted = `${{Number(identity).toFixed(2)}}%`;
          if (!Number.isFinite(Number(alignedPositions))) return formatted;
          return `${{formatted}} (${{Number(alignedPositions)}} aa)`;
        }}

        function mapSelectionToHomolog(homolog, start, end) {{
          const humanRegionStart = Number(homolog.humanRegionStart || ectodomain.start || 1);
          const humanRegionEnd = Number(homolog.humanRegionEnd || ectodomain.end || sequence.length);
          const homologRegionStart = Number(homolog.homologRegionStart || homolog.ectodomainStart || 1);
          if (!(humanRegionStart && humanRegionEnd) || start < humanRegionStart || end > humanRegionEnd) {{
            return null;
          }}
          const queryStart = start - humanRegionStart + 1;
          const queryEnd = end - humanRegionStart + 1;
          let queryPos = 0;
          let subjectPos = 0;
          const mapped = [];
          let matches = 0;
          let alignedPositions = 0;
          let surfaceMatches = 0;
          let surfaceAlignedPositions = 0;
          const surfacePositions = new Set(extracellularSurfacePositionsInRange(start, end));
          const alignedHuman = homolog.alignedHumanEctodomain || '';
          const alignedHomolog = homolog.alignedHomologEctodomain || '';
          for (let i = 0; i < alignedHuman.length; i += 1) {{
            const q = alignedHuman[i];
            const s = alignedHomolog[i];
            if (q !== '-') queryPos += 1;
            if (s !== '-') subjectPos += 1;
            if (q === '-' || queryPos < queryStart || queryPos > queryEnd) continue;
            if (s === '-') continue;
            mapped.push(subjectPos);
            alignedPositions += 1;
            if (q === s) matches += 1;
            const humanGlobalPosition = humanRegionStart + queryPos - 1;
            if (surfacePositions.has(humanGlobalPosition)) {{
              surfaceAlignedPositions += 1;
              if (q === s) surfaceMatches += 1;
            }}
          }}
          if (!mapped.length) return null;
          const requestedPositions = Math.max(0, queryEnd - queryStart + 1);
          if (!requestedPositions || alignedPositions / requestedPositions < 0.70) return null;
          const mappedStart = Math.min(...mapped);
          const mappedEnd = Math.max(...mapped);
          const homologSequence = (homolog.ectodomainSequence || '').slice(mappedStart - 1, mappedEnd);
          const identity = alignedPositions ? (100 * matches / alignedPositions) : null;
          const surfaceIdentity = surfaceAlignedPositions ? (100 * surfaceMatches / surfaceAlignedPositions) : null;
          return {{
            start: homologRegionStart + mappedStart - 1,
            end: homologRegionStart + mappedEnd - 1,
            sequence: homologSequence,
            identity,
            surfaceIdentity,
            surfaceAlignedPositions,
          }};
        }}

        function regionName(geneSymbol, species, start, end) {{
          const geneToken = String(geneSymbol || 'GENE').toUpperCase().replace(/[^A-Z0-9]+/g, '_');
          return `${{geneToken}}_${{speciesCode(species)}}_${{start}}-${{end}}`;
        }}

        function mutationEditKey(edit) {{
          return `${{String(edit.from || 'X').toUpperCase()}}${{Number(edit.position)}}${{String(edit.to || '').toUpperCase()}}`;
        }}

        function mutationEditLabel(edit) {{
          return edit.label || mutationEditKey(edit);
        }}

        function parseMutationSuggestion(text) {{
          const match = String(text || '').trim().match(/^([A-Za-z*])(\\d+)([A-Za-z*])$/);
          if (!match) return null;
          const position = Number(match[2]);
          if (!Number.isFinite(position)) return null;
          return {{
            from: match[1].toUpperCase(),
            position,
            to: match[3].toUpperCase(),
            label: `${{match[1].toUpperCase()}}${{position}}${{match[3].toUpperCase()}}`,
          }};
        }}

        function normalizeMutationEdits(edits) {{
          const deduped = new Map();
          (edits || []).forEach((edit) => {{
            if (!edit) return;
            const position = Number(edit.position);
            const to = String(edit.to || '').toUpperCase();
            if (!Number.isFinite(position) || !to) return;
            const from = String(edit.from || '').toUpperCase();
            const normalized = {{
              ...edit,
              position,
              from,
              to,
              label: edit.label || `${{from || 'X'}}${{position}}${{to}}`,
            }};
            deduped.set(mutationEditKey(normalized), normalized);
          }});
          return Array.from(deduped.values()).sort((a, b) => a.position - b.position || mutationEditLabel(a).localeCompare(mutationEditLabel(b)));
        }}

        function cysteineEditsFromPositions(positions) {{
          return normalizeMutationEdits((positions || []).map((position) => ({{
            position,
            from: 'C',
            to: 'S',
            label: `C${{Number(position)}}S`,
            kind: 'cysteine',
          }})));
        }}

        function furinSiteKey(site) {{
          return `${{Number(site.start)}}-${{Number(site.end)}}:${{String(site.motif || '')}}`;
        }}

        function furinEditsForSite(site) {{
          const siteKey = furinSiteKey(site);
          const parsed = (site.suggestions || [])
            .map(parseMutationSuggestion)
            .filter(Boolean)
            .map((edit) => ({{
              ...edit,
              kind: 'furin',
              siteKey,
              siteStart: Number(site.start),
              siteEnd: Number(site.end),
              motif: site.motif || '',
            }}));
          if (parsed.length) return normalizeMutationEdits(parsed);
          const motif = String(site.motif || '');
          const start = Number(site.start);
          const end = Number(site.end);
          const fallback = [];
          if (motif && Number.isFinite(start)) {{
            fallback.push({{ position: start, from: motif[0], to: 'A', kind: 'furin', siteKey, siteStart: start, siteEnd: end, motif }});
          }}
          if (motif.length > 1 && Number.isFinite(end) && end !== start) {{
            fallback.push({{ position: end, from: motif[motif.length - 1], to: 'A', kind: 'furin', siteKey, siteStart: start, siteEnd: end, motif }});
          }}
          return normalizeMutationEdits(fallback);
        }}

        function mutationSuffixFromEdits(edits) {{
          const values = normalizeMutationEdits(edits);
          return values.length ? `_${{values.map(mutationEditLabel).join('_')}}` : '';
        }}

        function nameWithMutationEdits(geneSymbol, species, start, end, edits) {{
          return `${{regionName(geneSymbol, species, start, end)}}${{mutationSuffixFromEdits(edits)}}`;
        }}

        function mutateSequenceByGlobalEdits(sequenceText, boundaryStart, edits) {{
          const chars = String(sequenceText || '').split('');
          const appliedEdits = [];
          normalizeMutationEdits(edits).forEach((edit) => {{
            const localIndex = edit.position - Number(boundaryStart);
            if (localIndex < 0 || localIndex >= chars.length) return;
            if (edit.from && String(chars[localIndex] || '').toUpperCase() !== edit.from) return;
            chars[localIndex] = edit.to;
            appliedEdits.push(edit);
          }});
          return {{
            sequence: chars.join(''),
            appliedEdits,
          }};
        }}

        function homologPositionLookup(species) {{
          const key = String(species || '').toLowerCase();
          if (homologPositionLookupCache.has(key)) return homologPositionLookupCache.get(key);
          const homolog = (viewerData.homologMappings || []).find((item) => String(item.species || '').toLowerCase() === key);
          if (!homolog) {{
            homologPositionLookupCache.set(key, null);
            return null;
          }}
          const humanRegionStart = Number(homolog.humanRegionStart || ectodomain.start || 1);
          const homologRegionStart = Number(homolog.homologRegionStart || homolog.ectodomainStart || 1);
          const lookup = new Map();
          let humanPos = 0;
          let homologPos = 0;
          const alignedHuman = String(homolog.alignedHumanEctodomain || '');
          const alignedHomolog = String(homolog.alignedHomologEctodomain || '');
          for (let index = 0; index < alignedHuman.length; index += 1) {{
            const humanResidue = alignedHuman[index];
            const homologResidue = alignedHomolog[index];
            if (humanResidue !== '-') humanPos += 1;
            if (homologResidue !== '-') homologPos += 1;
            if (humanResidue === '-' || homologResidue === '-') continue;
            lookup.set(humanRegionStart + humanPos - 1, homologRegionStart + homologPos - 1);
          }}
          homologPositionLookupCache.set(key, lookup);
          return lookup;
        }}

        function mapTargetMutationEditsToSpecies(species, edits) {{
          const normalized = String(species || '').toLowerCase();
          const mutationEdits = normalizeMutationEdits(edits);
          if (normalized === targetSpeciesKey()) return mutationEdits;
          const lookup = homologPositionLookup(species);
          if (!lookup) return [];
          return normalizeMutationEdits(
            mutationEdits
              .map((edit) => {{
                const mappedPosition = lookup.get(edit.position);
                if (!Number.isFinite(mappedPosition)) return null;
                return {{
                  ...edit,
                  position: mappedPosition,
                  label: `${{edit.from || 'X'}}${{mappedPosition}}${{edit.to}}`,
                  mappedFromTargetPosition: edit.position,
                }};
              }})
              .filter(Boolean)
          );
        }}

        function selectionMutationCandidatePositions(start, end) {{
          const cysteines = viewerData.cysteineAnalysis || [];
          return cysteines
            .filter((finding) => Number(finding.position) >= start && Number(finding.position) <= end)
            .filter((finding) => !String(finding.warning || '').startsWith('Ambiguous'))
            .filter((finding) => finding.paired_with === null || finding.paired_with === undefined || Number(finding.paired_with) < start || Number(finding.paired_with) > end)
            .map((finding) => Number(finding.position))
            .filter((value) => Number.isFinite(value))
            .sort((a, b) => a - b);
        }}

        function selectionFurinMutationCandidateSites(start, end) {{
          const furinSites = viewerData.furinSites || [];
          return furinSites
            .filter((site) => Number(site.start) >= start && Number(site.end) <= end)
            .filter((site) => furinEditsForSite(site).length)
            .sort((a, b) => Number(a.start) - Number(b.start));
        }}

        function ptmCategory(feature) {{
          const metadata = feature.metadata || {{}};
          if (metadata.ptm_category) return metadata.ptm_category;
          const type = String(feature.type || '').toUpperCase();
          const description = String(feature.description || '').toLowerCase();
          if (type === 'CARBOHYD') return 'glycosylation';
          if (type === 'MOD_RES') return 'modified residue';
          if (type === 'LIPID') return 'lipidation';
          if (type === 'DISULFID') return 'disulfide';
          if (type === 'CROSSLNK') return 'cross-link';
          if (type === 'INIT_MET') return 'initiator methionine processing';
          if (['SIGNAL', 'PROPEP', 'PEPTIDE', 'CHAIN'].includes(type)) return 'molecule processing';
          if (description.includes('cleavage')) return 'cleavage site';
          return 'post-translational modification';
        }}

        function ptmWarning(feature) {{
          const type = String(feature.type || '').toUpperCase();
          const category = String(ptmCategory(feature) || '').toLowerCase();
          const description = String(feature.description || '').toLowerCase();
          if (['SIGNAL', 'PROPEP', 'PEPTIDE'].includes(type) || category.includes('cleavage') || description.includes('cleavage')) {{
            return 'Review cleavage/processing boundary.';
          }}
          return '';
        }}

        function closestUnpairedCysteineText(finding) {{
          if (finding.paired_with !== null && finding.paired_with !== undefined) return '';
          const partner = finding.closest_unpaired_cysteine;
          const distance = Number(finding.closest_unpaired_cysteine_distance);
          if (partner === null || partner === undefined || !Number.isFinite(distance)) return 'nearest unpaired cysteine: none detected';
          const atom = finding.closest_unpaired_cysteine_distance_atom || 'distance';
          return `nearest unpaired cysteine: C${{partner}}, ${{distance.toFixed(2)}} A (${{atom}})`;
        }}

        function ptmOverlapsSelection(feature, start, end) {{
          const featureStart = Number(feature.start);
          const featureEnd = Number(feature.end);
          if (!Number.isFinite(featureStart) || !Number.isFinite(featureEnd)) return false;
          const type = String(feature.type || '').toUpperCase();
          if (type === 'DISULFID') {{
            return (featureStart >= start && featureStart <= end) || (featureEnd >= start && featureEnd <= end);
          }}
          return featureStart <= end && featureEnd >= start;
        }}

        function ptmDisplayLabel(feature) {{
          const type = String(feature.type || '').toUpperCase();
          if (type === 'DISULFID') return `C${{feature.start}}-C${{feature.end}}`;
          return `${{feature.start}}-${{feature.end}}`;
        }}

        function selectionPtms(start, end) {{
          return (viewerData.ptms || [])
            .filter((feature) => ptmOverlapsSelection(feature, start, end))
            .sort((a, b) => Number(a.start) - Number(b.start) || Number(a.end) - Number(b.end));
        }}

        function escapeHtml(value) {{
          const element = document.createElement('span');
          element.textContent = String(value ?? '');
          return element.innerHTML.replaceAll('"', '&quot;').replaceAll("'", '&#39;');
        }}

        function renderSelectionTable() {{
          const bounds = selectionBounds();
          if (!bounds) {{
            selectionSummary.textContent = 'Select a residue or contiguous range to populate the table.';
            selectionWarnings.innerHTML = '<p class="viewer-info">Warnings for the selected region will appear here.</p>';
            selectionTableBody.innerHTML = '<tr><td colspan="7">No selection yet.</td></tr>';
            selectionTsv.value = '';
            selectionFasta.value = '';
            workbenchSelectionRows = [];
            renderWorkbenchSpecies();
            return;
          }}
          if (bounds.nonContiguous) {{
            selectionSummary.textContent = `${{bounds.label}} focused in the structure (${{bounds.count}} non-contiguous residue(s)). Select a contiguous range to export construct sequences.`;
            selectionWarnings.innerHTML = '<p class="viewer-info">Non-contiguous topology selections are for structural inspection only.</p>';
            selectionTableBody.innerHTML = '<tr><td colspan="7">Non-contiguous extracellular residues are focused in the 3D viewer; no construct sequence export is generated for this selection.</td></tr>';
            selectionTsv.value = '';
            selectionFasta.value = '';
            workbenchSelectionRows = [];
            renderWorkbenchSpecies();
            return;
          }}
          const mutationCandidates = selectionMutationCandidatePositions(bounds.start, bounds.end);
          const furinMutationCandidates = selectionFurinMutationCandidateSites(bounds.start, bounds.end);
          const allowedMutationSet = new Set(mutationCandidates);
          Array.from(selectionMutationPositions).forEach((position) => {{
            if (!allowedMutationSet.has(position)) selectionMutationPositions.delete(position);
          }});
          const allowedFurinSiteSet = new Set(furinMutationCandidates.map(furinSiteKey));
          Array.from(selectionFurinMutationSites).forEach((siteKey) => {{
            if (!allowedFurinSiteSet.has(siteKey)) selectionFurinMutationSites.delete(siteKey);
          }});
          const activeTargetMutationEdits = normalizeMutationEdits(
            cysteineEditsFromPositions(Array.from(selectionMutationPositions)).concat(
              furinMutationCandidates
                .filter((site) => selectionFurinMutationSites.has(furinSiteKey(site)))
                .flatMap(furinEditsForSite)
            )
          );
          const targetSequenceText = sequence.slice(bounds.start - 1, bounds.end);
          const mutatedTarget = mutateSequenceByGlobalEdits(targetSequenceText, bounds.start, activeTargetMutationEdits);
          const targetGeneSymbol = targetInfo.geneSymbol || targetInfo.entryName || 'TARGET';
          const targetNameSymbol = targetInfo.nameSymbol || String(targetInfo.entryName || '').replace(/_(HUMAN|MOUSE|MACFA)$/i, '') || targetGeneSymbol;
          const targetSpecies = targetSpeciesKey();
          const rows = [
            {{
              name: nameWithMutationEdits(targetNameSymbol, targetSpecies, bounds.start, bounds.end, mutatedTarget.appliedEdits),
              species: speciesCode(targetSpecies),
              boundary: `${{bounds.start}}-${{bounds.end}}`,
              sequence: mutatedTarget.sequence,
              identity: '100.00%',
              surfaceIdentity: '100.00%',
              accession: targetInfo.accession || '',
              entryName: targetInfo.entryName || '',
              geneSymbol: targetGeneSymbol,
              links: resourceLinksHtml({{
                accession: targetInfo.accession || '',
                entryName: targetInfo.entryName || '',
                geneSymbol: targetGeneSymbol,
                species: targetSpecies,
              }}),
            }}
          ];
          (viewerData.homologMappings || []).forEach((homolog) => {{
            const mapped = mapSelectionToHomolog(homolog, bounds.start, bounds.end);
            if (!mapped) return;
            const homologMutationEdits = mapTargetMutationEditsToSpecies(homolog.species, activeTargetMutationEdits);
            const mutatedHomolog = mutateSequenceByGlobalEdits(mapped.sequence || '', mapped.start, homologMutationEdits);
            rows.push({{
              name: nameWithMutationEdits(targetNameSymbol, homolog.species, mapped.start, mapped.end, mutatedHomolog.appliedEdits),
              species: speciesCode(homolog.species),
              boundary: `${{mapped.start}}-${{mapped.end}}`,
              sequence: mutatedHomolog.sequence,
              identity: mapped.identity === null ? '' : `${{mapped.identity.toFixed(2)}}%`,
              surfaceIdentity: surfaceIdentityText(mapped.surfaceIdentity, mapped.surfaceAlignedPositions),
              accession: homolog.accession || '',
              entryName: homolog.entryName || homolog.accession || '',
              geneSymbol: targetGeneSymbol,
              links: resourceLinksHtml({{
                accession: homolog.accession || '',
                entryName: homolog.entryName || homolog.accession || '',
                geneSymbol: targetGeneSymbol,
                species: homolog.species,
              }}),
            }});
          }});
          workbenchSelectionRows = rows;
          renderWorkbenchSpecies();
          selectionSummary.textContent = `${{bounds.label}} ready for copy / paste.`;
          selectionTableBody.innerHTML = rows.map((row) => `
            <tr>
              <td><code>${{escapeHtml(row.name)}}</code></td>
              <td>${{escapeHtml(row.species)}}</td>
              <td>${{escapeHtml(row.boundary)}}</td>
              <td><code class="sequence">${{escapeHtml(row.sequence)}}</code></td>
              <td>${{escapeHtml(row.identity)}}</td>
              <td>${{escapeHtml(row.surfaceIdentity || '')}}</td>
              <td>${{row.links}}</td>
            </tr>
          `).join('');
          const identityHeader = targetInfo.identityHeader || 'identity_to_human';
          const surfaceIdentityHeader = targetInfo.surfaceIdentityHeader || 'extracellular_surface_identity_to_human';
          const header = ['name', 'species', 'boundary', 'sequence', identityHeader, surfaceIdentityHeader];
          selectionTsv.value = [header.join('\\t')].concat(
            rows.map((row) => [row.name, row.species, row.boundary, row.sequence, row.identity, row.surfaceIdentity || ''].join('\\t'))
          ).join('\\n');
          selectionFasta.value = rows.map((row) => `>${{row.name}} species=${{row.species}} boundary=${{row.boundary}} ${{identityHeader}}=${{row.identity || 'n/a'}} ${{surfaceIdentityHeader}}=${{row.surfaceIdentity || 'n/a'}}\\n${{row.sequence.replace(/(.{{1,80}})/g, '$1\\n').trim()}}`).join('\\n');
          renderSelectionWarnings(bounds.start, bounds.end, mutationCandidates, furinMutationCandidates);
        }}

        function renderSelectionWarnings(start, end, mutationCandidates, furinMutationCandidates) {{
          const warnings = [];
          const cysteines = viewerData.cysteineAnalysis || [];
          const furinSites = viewerData.furinSites || [];
          const ptms = selectionPtms(start, end);
          const selectedCysteines = cysteines.filter((finding) => Number(finding.position) >= start && Number(finding.position) <= end);
          const unpairedForSelection = selectedCysteines.filter((finding) => {{
            if (finding.paired_with === null || finding.paired_with === undefined) return true;
            return Number(finding.paired_with) < start || Number(finding.paired_with) > end;
          }});
          if (unpairedForSelection.length) {{
            const parts = unpairedForSelection.map((finding) => {{
              const position = `C${{finding.position}}`;
              if (finding.paired_with === null || finding.paired_with === undefined) {{
                const details = [finding.warning || 'no paired cysteine detected', closestUnpairedCysteineText(finding)].filter(Boolean).join('; ');
                return `${{position}} (${{details}})`;
              }}
              return `${{position}} (partner C${{finding.paired_with}} lies outside the selected region)`;
            }});
            warnings.push({{
              kind: 'cysteine',
              title: 'Unpaired cysteine warning',
              items: parts,
            }});
          }}
          const containedFurinSites = furinSites.filter((site) => Number(site.start) >= start && Number(site.end) <= end);
          if (containedFurinSites.length) {{
            warnings.push({{
              kind: 'furin',
              title: 'Candidate furin-like motif warning',
              items: containedFurinSites.map((site) => `${{site.motif}} at ${{site.start}}-${{site.end}}; suggested edits: ${{(site.suggestions || []).join(', ')}}`),
            }});
          }}
          if (ptms.length) {{
            warnings.push({{
              kind: 'ptm',
              title: 'PTMs included in selected region',
              items: ptms.map((feature) => {{
                const base = `${{ptmDisplayLabel(feature)}} ${{ptmCategory(feature)}}: ${{feature.description || feature.type || ''}}`;
                const warning = ptmWarning(feature);
                return warning ? `${{base}} (${{warning}})` : base;
              }}),
            }});
          }}
          const warningHtml = warnings.length
            ? warnings.map((warning) => `
            <div class="selection-warning-block ${{warning.kind}}">
              <strong>${{warning.title}}</strong>
              <ul>${{warning.items.map((item) => `<li>${{escapeHtml(item)}}</li>`).join('')}}</ul>
            </div>
          `).join('')
            : '<p class="viewer-info">No unpaired-cysteine or furin-site warnings for the current selection.</p>';
          const mutationHtml = (mutationCandidates || []).length
            ? `
              <div class="mutation-controls">
                <strong>Optional Cys→Ser edits</strong>
                <div class="mutation-checkboxes">
                  ${{mutationCandidates.map((position) => `
                    <label class="mutation-checkbox">
                      <input type="checkbox" class="selection-mutation-checkbox" data-position="${{position}}" ${{selectionMutationPositions.has(position) ? 'checked' : ''}}>
                      C${{position}}S
                    </label>
                  `).join('')}}
                </div>
              </div>
            `
            : '';
          const furinMutationHtml = (furinMutationCandidates || []).length
            ? `
              <div class="mutation-controls">
                <strong>Optional furin-site edits</strong>
                <div class="mutation-checkboxes">
                  ${{furinMutationCandidates.map((site) => {{
                    const siteKey = furinSiteKey(site);
                    const edits = furinEditsForSite(site).map(mutationEditLabel).join(', ');
                    return `
                    <label class="mutation-checkbox">
                      <input type="checkbox" class="selection-furin-mutation-checkbox" data-site-key="${{siteKey}}" ${{selectionFurinMutationSites.has(siteKey) ? 'checked' : ''}}>
                      ${{escapeHtml(site.motif || 'furin site')}} ${{site.start}}-${{site.end}} (${{edits}})
                    </label>
                  `;
                  }}).join('')}}
                </div>
              </div>
            `
            : '';
          selectionWarnings.innerHTML = `${{warningHtml}}${{mutationHtml}}${{furinMutationHtml}}`;
          selectionWarnings.querySelectorAll('.selection-mutation-checkbox').forEach((checkbox) => {{
            checkbox.addEventListener('change', () => {{
              const position = Number(checkbox.dataset.position);
              if (checkbox.checked) {{
                selectionMutationPositions.add(position);
              }} else {{
                selectionMutationPositions.delete(position);
              }}
              renderSelectionTable();
            }});
          }});
          selectionWarnings.querySelectorAll('.selection-furin-mutation-checkbox').forEach((checkbox) => {{
            checkbox.addEventListener('change', () => {{
              const siteKey = checkbox.dataset.siteKey || '';
              if (checkbox.checked) {{
                selectionFurinMutationSites.add(siteKey);
              }} else {{
                selectionFurinMutationSites.delete(siteKey);
              }}
              renderSelectionTable();
            }});
          }});
        }}

        function renderConstructMutationCard(card) {{
          const payloadEl = card.querySelector('.construct-payload');
          if (!payloadEl) return;
          let payload = null;
          try {{
            payload = JSON.parse(payloadEl.textContent || '{{}}');
          }} catch (error) {{
            throw new Error(`Invalid construct payload: ${{error.message}}`);
          }}
          const controlsEl = card.querySelector('.construct-mutation-controls');
          const tableBody = card.querySelector('.construct-sequence-table-body');
          const tsvEl = card.querySelector('.construct-export-tsv');
          const fastaEl = card.querySelector('.construct-export-fasta');
          if (!controlsEl || !tableBody || !tsvEl || !fastaEl) return;
          const activeMutations = new Set();
          const activeFurinSites = new Set();
          const furinMutationCandidates = (payload.furinSites || []).filter((site) => furinEditsForSite(site).length);

          function redraw() {{
            const targetMutationEdits = normalizeMutationEdits(
              cysteineEditsFromPositions(Array.from(activeMutations)).concat(
                furinMutationCandidates
                  .filter((site) => activeFurinSites.has(furinSiteKey(site)))
                  .flatMap(furinEditsForSite)
              )
            );
            const rows = (payload.rows || []).map((row) => {{
              const rowMutationEdits = String(row.speciesRaw || '').toLowerCase() === targetSpeciesKey()
                ? targetMutationEdits
                : mapTargetMutationEditsToSpecies(row.speciesRaw, targetMutationEdits);
              const mutated = mutateSequenceByGlobalEdits(row.sequence || '', Number(row.start), rowMutationEdits);
              return {{
                name: nameWithMutationEdits(payload.geneSymbol, row.speciesRaw, row.start, row.end, mutated.appliedEdits),
                species: row.species || speciesCode(row.speciesRaw),
                boundary: `${{row.start}}-${{row.end}}`,
                sequence: mutated.sequence,
                identity: row.identity || '',
                surfaceIdentity: row.surfaceIdentity || '',
                accession: row.accession || '',
                entryName: row.entryName || '',
                geneSymbol: row.geneSymbol || payload.geneSymbol || '',
                speciesRaw: row.speciesRaw || '',
              }};
            }});
            tableBody.innerHTML = rows.map((row) => `
              <tr>
                <td><code>${{escapeHtml(row.name)}}</code></td>
                <td>${{escapeHtml(row.species)}}</td>
                <td>${{escapeHtml(row.boundary)}}</td>
                <td><code class="sequence">${{escapeHtml(row.sequence)}}</code></td>
                <td>${{escapeHtml(row.identity)}}</td>
                <td>${{escapeHtml(row.surfaceIdentity)}}</td>
                <td>${{resourceLinksHtml({{
                  accession: row.accession,
                  entryName: row.entryName,
                  geneSymbol: row.geneSymbol,
                  species: row.speciesRaw,
                }})}}</td>
              </tr>
            `).join('') + (payload.unavailableRows || []).map((row) => `
              <tr>
                <td><code>${{escapeHtml(row.name)}}</code></td>
                <td>${{escapeHtml(row.species)}}</td>
                <td>unavailable</td>
                <td><span class="muted">Mapping unavailable: ${{escapeHtml(row.reason)}}</span></td>
                <td>n/a</td><td>n/a</td>
                <td>${{resourceLinksHtml({{
                  accession: row.accession,
                  entryName: row.entryName,
                  geneSymbol: row.geneSymbol,
                  species: row.speciesRaw,
                }})}}</td>
              </tr>
            `).join('');
            const identityHeader = targetInfo.identityHeader || 'identity_to_human';
            const surfaceIdentityHeader = targetInfo.surfaceIdentityHeader || 'extracellular_surface_identity_to_human';
            const header = ['name', 'species', 'boundary', 'sequence', identityHeader, surfaceIdentityHeader];
            tsvEl.value = [header.join('\\t')].concat(
              rows.map((row) => [row.name, row.species, row.boundary, row.sequence, row.identity, row.surfaceIdentity].join('\\t'))
            ).join('\\n');
            fastaEl.value = rows.map((row) => `>${{row.name}} species=${{row.species}} boundary=${{row.boundary}} ${{identityHeader}}=${{row.identity || 'n/a'}} ${{surfaceIdentityHeader}}=${{row.surfaceIdentity || 'n/a'}}\\n${{row.sequence.replace(/(.{{1,80}})/g, '$1\\n').trim()}}`).join('\\n');
            if (!(payload.mutationCandidates || []).length && !furinMutationCandidates.length) {{
              controlsEl.innerHTML = '<p class="viewer-info">No unpaired cysteines or furin sites in this construct.</p>';
              return;
            }}
            const cysteineControls = (payload.mutationCandidates || []).length
              ? `
                <div>
                  <strong>Optional Cys→Ser edits</strong>
                  <div class="mutation-checkboxes">
                    ${{payload.mutationCandidates.map((position) => `
                      <label class="mutation-checkbox">
                        <input type="checkbox" class="construct-mutation-checkbox" data-position="${{position}}" ${{activeMutations.has(position) ? 'checked' : ''}}>
                        C${{position}}S
                      </label>
                    `).join('')}}
                  </div>
                </div>
              `
              : '';
            const furinControls = furinMutationCandidates.length
              ? `
                <div>
                  <strong>Optional furin-site edits</strong>
                  <div class="mutation-checkboxes">
                    ${{furinMutationCandidates.map((site) => {{
                      const siteKey = furinSiteKey(site);
                      const edits = furinEditsForSite(site).map(mutationEditLabel).join(', ');
                      return `
                        <label class="mutation-checkbox">
                          <input type="checkbox" class="construct-furin-mutation-checkbox" data-site-key="${{siteKey}}" ${{activeFurinSites.has(siteKey) ? 'checked' : ''}}>
                          ${{escapeHtml(site.motif || 'furin site')}} ${{site.start}}-${{site.end}} (${{edits}})
                        </label>
                      `;
                    }}).join('')}}
                  </div>
                </div>
              `
              : '';
            controlsEl.innerHTML = `${{cysteineControls}}${{furinControls}}`;
            controlsEl.querySelectorAll('.construct-mutation-checkbox').forEach((checkbox) => {{
              checkbox.addEventListener('change', () => {{
                const position = Number(checkbox.dataset.position);
                if (checkbox.checked) {{
                  activeMutations.add(position);
                }} else {{
                  activeMutations.delete(position);
                }}
                redraw();
              }});
            }});
            controlsEl.querySelectorAll('.construct-furin-mutation-checkbox').forEach((checkbox) => {{
              checkbox.addEventListener('change', () => {{
                const siteKey = checkbox.dataset.siteKey || '';
                if (checkbox.checked) {{
                  activeFurinSites.add(siteKey);
                }} else {{
                  activeFurinSites.delete(siteKey);
                }}
                redraw();
              }});
            }});
          }}

          redraw();
        }}

        function renderPlddtPlot() {{
          const width = plddtLayout.width;
          const height = plddtLayout.height;
          const margin = plddtLayout.margin;
          const plotWidth = width - margin.left - margin.right;
          const plotHeight = height - margin.top - margin.bottom;
          const plddt = viewerData.plddt || [];
          if (!plddt.length) {{
            plddtSvg.innerHTML = '';
            plddtInfo.textContent = 'No pLDDT values available.';
            return;
          }}
          const viewStart = plddtView.start;
          const viewEnd = plddtView.end;
          const xFor = (resi) => margin.left + ((resi - viewStart) / Math.max(viewLength(plddtView) - 1, 1)) * plotWidth;
          const yFor = (value) => margin.top + (1 - Math.max(0, Math.min(100, value)) / 100) * plotHeight;
          const residues = makeRange(viewStart, viewEnd);
          const pathSegments = [];
          for (let index = 0; index < residues.length - 1; index += 1) {{
            const residue = residues[index];
            const nextResidue = residues[index + 1];
            const x1 = xFor(residue).toFixed(2);
            const y1 = yFor(plddt[residue - 1] ?? 0).toFixed(2);
            const x2 = xFor(nextResidue).toFixed(2);
            const y2 = yFor(plddt[nextResidue - 1] ?? 0).toFixed(2);
            const segmentScore = ((Number(plddt[residue - 1]) || 0) + (Number(plddt[nextResidue - 1]) || 0)) / 2;
            pathSegments.push(`<line x1="${{x1}}" y1="${{y1}}" x2="${{x2}}" y2="${{y2}}" stroke="${{plddtColor(segmentScore)}}" stroke-width="3" stroke-linecap="round"/>`);
          }}
          let highlight = '';
          if (selectedRange) {{
            const visibleStart = Math.max(selectedRange.start, viewStart);
            const visibleEnd = Math.min(selectedRange.end, viewEnd);
            if (visibleStart <= visibleEnd) {{
              const x1 = xFor(visibleStart);
              const x2 = xFor(visibleEnd);
              const residueWidth = plotWidth / Math.max(viewLength(plddtView), 1);
              highlight = `<rect x="${{Math.min(x1, x2)}}" y="${{margin.top}}" width="${{Math.max(2, Math.abs(x2 - x1) + residueWidth)}}" height="${{plotHeight}}" fill="rgba(76,149,108,0.18)" stroke="#4c956c" stroke-width="2"/>`;
            }}
            const stats = plddtStats(selectedRange.start, selectedRange.end);
            const statsText = stats ? ` | pLDDT min/mean/max ${{stats.min.toFixed(1)}} / ${{stats.mean.toFixed(1)}} / ${{stats.max.toFixed(1)}}` : '';
            plddtInfo.textContent = `${{selectedRange.label || 'Selected range'}}: residues ${{selectedRange.start}}-${{selectedRange.end}} | view ${{viewStart}}-${{viewEnd}}${{statsText}}`;
          }} else if (selectedResidue !== null && selectedResidue >= viewStart && selectedResidue <= viewEnd) {{
            const x = xFor(selectedResidue);
            highlight = `<line x1="${{x}}" y1="${{margin.top}}" x2="${{x}}" y2="${{margin.top + plotHeight}}" stroke="#c1121f" stroke-width="2"/>`;
            plddtInfo.textContent = `Residue ${{selectedResidue}} pLDDT: ${{(plddt[selectedResidue - 1] ?? 0).toFixed(2)}} | view ${{viewStart}}-${{viewEnd}}`;
          }} else {{
            plddtInfo.textContent = `Viewing ${{designRegionLower}} residues ${{viewStart}}-${{viewEnd}}. Use the mouse wheel or buttons to zoom.`;
          }}
          const bands = [
            {{ start: 0, end: 50, color: 'rgba(255,125,69,0.14)' }},
            {{ start: 50, end: 70, color: 'rgba(255,219,19,0.16)' }},
            {{ start: 70, end: 90, color: 'rgba(101,203,243,0.16)' }},
            {{ start: 90, end: 100, color: 'rgba(0,83,214,0.12)' }},
          ].map((band) => {{
            const top = yFor(band.end);
            const bottom = yFor(band.start);
            return `<rect x="${{margin.left}}" y="${{top}}" width="${{plotWidth}}" height="${{Math.max(0, bottom - top)}}" fill="${{band.color}}"/>`;
          }}).join('');
          const topologyBands = topologyIntervals(viewStart, viewEnd).map((interval) => {{
            const x1 = xFor(interval.start);
            const x2 = xFor(interval.end);
            const residueWidth = plotWidth / Math.max(viewLength(plddtView), 1);
            const left = Math.min(x1, x2);
            const bandWidth = Math.max(2, Math.abs(x2 - x1) + residueWidth);
            return `
              <rect x="${{left}}" y="${{margin.top}}" width="${{bandWidth}}" height="${{plotHeight}}" fill="${{topologyColor(interval.location, 0.07)}}"/>
              <rect x="${{left}}" y="${{margin.top + plotHeight - 8}}" width="${{bandWidth}}" height="8" fill="${{topologyColor(interval.location, 0.78)}}"/>
            `;
          }}).join('');
          const gpcrSegmentBands = gpcrSegments
            .filter((segment) => segment.end >= viewStart && segment.start <= viewEnd)
            .map((segment) => {{
              const visibleStart = Math.max(segment.start, viewStart);
              const visibleEnd = Math.min(segment.end, viewEnd);
              const x1 = xFor(visibleStart);
              const x2 = xFor(visibleEnd);
              const residueWidth = plotWidth / Math.max(viewLength(plddtView), 1);
              const left = Math.min(x1, x2);
              const bandWidth = Math.max(2, Math.abs(x2 - x1) + residueWidth);
              const labelX = left + bandWidth / 2;
              const label = segment.name;
              const bandY = margin.top + plotHeight - 26;
              return `
                <rect x="${{left}}" y="${{bandY}}" width="${{bandWidth}}" height="16" rx="5" fill="rgba(31,31,31,0.78)"/>
                <text x="${{labelX}}" y="${{bandY + 12}}" text-anchor="middle" fill="white" font-size="9" font-weight="700">${{label}}</text>
              `;
            }}).join('');
          const guides = [50, 70, 90].map((value) => {{
            const y = yFor(value);
            return `<line x1="${{margin.left}}" y1="${{y}}" x2="${{margin.left + plotWidth}}" y2="${{y}}" stroke="#e7ddd1" stroke-dasharray="4 4"/><text x="6" y="${{y + 4}}" fill="#6f675b" font-size="11">${{value}}</text>`;
          }}).join('');
          const ticks = Array.from(new Set([viewStart, Math.round((viewStart + viewEnd) / 2), viewEnd])).map((resi) => {{
            const x = xFor(resi);
            return `<line x1="${{x}}" y1="${{margin.top + plotHeight}}" x2="${{x}}" y2="${{margin.top + plotHeight + 5}}" stroke="#bcae99"/><text x="${{x}}" y="${{height - 6}}" text-anchor="middle" fill="#6f675b" font-size="11">${{resi}}</text>`;
          }}).join('');
          plddtSvg.innerHTML = `
            <rect x="0" y="0" width="${{width}}" height="${{height}}" fill="white" rx="12"/>
            ${{bands}}
            ${{topologyBands}}
            ${{gpcrSegmentBands}}
            <line x1="${{margin.left}}" y1="${{margin.top + plotHeight}}" x2="${{margin.left + plotWidth}}" y2="${{margin.top + plotHeight}}" stroke="#bcae99"/>
            <line x1="${{margin.left}}" y1="${{margin.top}}" x2="${{margin.left}}" y2="${{margin.top + plotHeight}}" stroke="#bcae99"/>
            ${{guides}}
            ${{highlight}}
            ${{pathSegments.join('')}}
            ${{ticks}}
          `;
        }}

        let paeHeatmapCache = null;
        function renderPaePlot() {{
          const ctx = paeCanvas.getContext('2d');
          const matrix = viewerData.paeMatrix || [];
          if (!matrix.length || !paeView) {{
            ctx.clearRect(0, 0, paeCanvas.width, paeCanvas.height);
            paeInfo.textContent = 'No PAE matrix available.';
            return;
          }}
          const width = paeCanvas.width;
          const height = paeCanvas.height;
          const size = matrix.length;
          const viewStart = paeView.start;
          const viewEnd = paeView.end;
          const paeWindow = viewLength(paeView);
          const cacheKey = `${{width}}:${{height}}:${{viewStart}}:${{viewEnd}}:${{viewerData.paeMax}}`;
          let image;
          if (paeHeatmapCache && paeHeatmapCache.key === cacheKey && paeHeatmapCache.matrix === matrix) {{
            image = paeHeatmapCache.image;
          }} else {{
            image = ctx.createImageData(width, height);
            for (let y = 0; y < height; y += 1) {{
              const rowResidue = viewStart + (y / Math.max(height - 1, 1)) * Math.max(paeWindow - 1, 1);
              const row = paeMatrixIndexForResidue(rowResidue, size);
              for (let x = 0; x < width; x += 1) {{
                const colResidue = viewStart + (x / Math.max(width - 1, 1)) * Math.max(paeWindow - 1, 1);
                const col = paeMatrixIndexForResidue(colResidue, size);
                const value = Number(matrix[row][col]) || 0;
                const t = Math.max(0, Math.min(1, value / (viewerData.paeMax || 30)));
                const [r, g, b] = plasmaColor(t);
                const idx = (y * width + x) * 4;
                image.data[idx] = r;
                image.data[idx + 1] = g;
                image.data[idx + 2] = b;
                image.data[idx + 3] = 255;
              }}
            }}
            paeHeatmapCache = {{ key: cacheKey, matrix, image }};
          }}
          ctx.putImageData(image, 0, 0);
          const stripWidth = Math.max(10, Math.round(Math.min(width, height) * 0.026));
          for (const interval of topologyIntervals(viewStart, viewEnd)) {{
            const startPixel = ((interval.start - viewStart) / Math.max(paeWindow, 1)) * width;
            const endPixel = ((interval.end - viewStart + 1) / Math.max(paeWindow, 1)) * width;
            const left = Math.max(0, Math.min(width, startPixel));
            const right = Math.max(0, Math.min(width, endPixel));
            const span = Math.max(1, right - left);
            ctx.fillStyle = topologyColor(interval.location, 0.92);
            ctx.fillRect(left, 0, span, stripWidth);
            ctx.fillRect(0, left, stripWidth, span);
            ctx.strokeStyle = topologyColor(interval.location, 0.18);
            ctx.lineWidth = 1;
            ctx.strokeRect(left, left, span, span);
          }}
          ctx.strokeStyle = '#bcae99';
          ctx.lineWidth = 1;
          ctx.strokeRect(0.5, 0.5, width - 1, height - 1);
          if (selectedRange) {{
            const visibleStart = Math.max(selectedRange.start, viewStart);
            const visibleEnd = Math.min(selectedRange.end, viewEnd);
            if (visibleStart <= visibleEnd) {{
              const x1 = ((visibleStart - viewStart) / Math.max(paeWindow, 1)) * width;
              const x2 = ((visibleEnd - viewStart + 1) / Math.max(paeWindow, 1)) * width;
              const left = Math.max(0, Math.min(width, Math.min(x1, x2)));
              const right = Math.max(0, Math.min(width, Math.max(x1, x2)));
              ctx.fillStyle = 'rgba(76,149,108,0.12)';
              ctx.fillRect(left, left, Math.max(2, right - left), Math.max(2, right - left));
              ctx.strokeStyle = '#4c956c';
              ctx.lineWidth = 3;
              ctx.strokeRect(left, left, Math.max(2, right - left), Math.max(2, right - left));
              ctx.fillStyle = '#4c956c';
              ctx.beginPath(); ctx.arc(left, left, 5, 0, Math.PI * 2); ctx.fill();
              ctx.beginPath(); ctx.arc(right, right, 5, 0, Math.PI * 2); ctx.fill();
            }}
            paeInfo.textContent = `${{selectedRange.label || 'Selected range'}} on ${{designRegionLower}} PAE | view ${{viewStart}}-${{viewEnd}}`;
          }} else if (selectedResidue !== null && selectedResidue >= viewStart && selectedResidue <= viewEnd) {{
            const x = ((selectedResidue - viewStart + 0.5) / Math.max(paeWindow, 1)) * width;
            ctx.strokeStyle = '#c1121f';
            ctx.lineWidth = 2;
            ctx.beginPath(); ctx.moveTo(x, 0); ctx.lineTo(x, height); ctx.stroke();
            ctx.beginPath(); ctx.moveTo(0, x); ctx.lineTo(width, x); ctx.stroke();
            paeInfo.textContent = `Residue ${{selectedResidue}} crosshair on ${{designRegionLower}} PAE | view ${{viewStart}}-${{viewEnd}}`;
          }} else {{
            paeInfo.textContent = `Viewing ${{designRegionLower}} residues ${{viewStart}}-${{viewEnd}}. Use the mouse wheel or buttons to zoom.`;
          }}
        }}

        residueSpans.forEach((span) => {{
          span.addEventListener('mousedown', (event) => {{
            event.preventDefault();
            const residueIndex = Number(span.dataset.resi);
            if (event.shiftKey && anchorResidue !== null) {{
              setRangeFromSelection(anchorResidue, residueIndex, 'Sequence range');
              return;
            }}
            isDragging = true;
            dragStartResidue = residueIndex;
            setRangeFromSelection(residueIndex, residueIndex, 'Sequence range');
          }});
          span.addEventListener('mouseenter', () => {{
            if (!isDragging || dragStartResidue === null) return;
            const residueIndex = Number(span.dataset.resi);
            setRangeFromSelection(dragStartResidue, residueIndex, 'Sequence range');
          }});
        }});

        plddtSvg.addEventListener('mousedown', (event) => {{
          event.preventDefault();
          const resizeEdge = plddtEdgeAtEvent(event);
          if (resizeEdge) {{
            plotResize = {{ source: 'plddt', edge: resizeEdge }};
            return;
          }}
          const residueIndex = residueFromPlddtPointer(event);
          if (event.shiftKey && anchorResidue !== null) {{
            setRangeFromSelection(anchorResidue, residueIndex, 'pLDDT range');
            return;
          }}
          plddtDragStart = residueIndex;
          plddtDragMoved = false;
          setRangeFromSelection(residueIndex, residueIndex, 'pLDDT range');
        }});

        plddtSvg.addEventListener('mousemove', (event) => {{
          if (plotResize && plotResize.source === 'plddt' && selectedRange) {{
            const residueIndex = residueFromPlddtPointer(event);
            if (plotResize.edge === 'left') {{
              setRangeFromSelection(residueIndex, selectedRange.end, selectedRange.label || 'pLDDT range');
            }} else {{
              setRangeFromSelection(selectedRange.start, residueIndex, selectedRange.label || 'pLDDT range');
            }}
            updatePlotCursor(plddtSvg, plotResize.edge);
            return;
          }}
          if (plddtDragStart === null) return;
          const residueIndex = residueFromPlddtPointer(event);
          plddtDragMoved = plddtDragMoved || residueIndex !== plddtDragStart;
          setRangeFromSelection(plddtDragStart, residueIndex, 'pLDDT range');
        }});

        plddtSvg.addEventListener('mouseleave', () => {{
          if (!plotResize || plotResize.source !== 'plddt') {{
            updatePlotCursor(plddtSvg, null);
          }}
        }});

        plddtSvg.addEventListener('mousemove', (event) => {{
          if (plotResize && plotResize.source === 'plddt') return;
          updatePlotCursor(plddtSvg, plddtEdgeAtEvent(event));
        }});

        plddtSvg.addEventListener('wheel', (event) => {{
          event.preventDefault();
          zoomPlot('plddt', event.deltaY < 0 ? 0.75 : 1.35, residueFromPlddtPointer(event));
        }}, {{ passive: false }});

        paeCanvas.addEventListener('mousedown', (event) => {{
          event.preventDefault();
          const resizeEdge = paeEdgeAtEvent(event);
          if (resizeEdge) {{
            plotResize = {{ source: 'pae', edge: resizeEdge }};
            return;
          }}
          const residueIndex = residueFromPaePointer(event);
          if (residueIndex === null) return;
          if (event.shiftKey && anchorResidue !== null) {{
            setRangeFromSelection(anchorResidue, residueIndex, 'PAE range');
            return;
          }}
          paeDragStart = residueIndex;
          paeDragMoved = false;
          setRangeFromSelection(residueIndex, residueIndex, 'PAE range');
        }});

        paeCanvas.addEventListener('mousemove', (event) => {{
          if (plotResize && plotResize.source === 'pae' && selectedRange) {{
            const residueIndex = residueFromPaePointer(event);
            if (residueIndex === null) return;
            if (plotResize.edge === 'left') {{
              setRangeFromSelection(residueIndex, selectedRange.end, selectedRange.label || 'PAE range');
            }} else {{
              setRangeFromSelection(selectedRange.start, residueIndex, selectedRange.label || 'PAE range');
            }}
            updatePlotCursor(paeCanvas, plotResize.edge);
            return;
          }}
          if (paeDragStart === null) return;
          const residueIndex = residueFromPaePointer(event);
          if (residueIndex === null) return;
          paeDragMoved = paeDragMoved || residueIndex !== paeDragStart;
          setRangeFromSelection(paeDragStart, residueIndex, 'PAE range');
        }});

        paeCanvas.addEventListener('mouseleave', () => {{
          if (!plotResize || plotResize.source !== 'pae') {{
            updatePlotCursor(paeCanvas, null);
          }}
        }});

        paeCanvas.addEventListener('mousemove', (event) => {{
          if (plotResize && plotResize.source === 'pae') return;
          updatePlotCursor(paeCanvas, paeEdgeAtEvent(event));
        }});

        paeCanvas.addEventListener('wheel', (event) => {{
          event.preventDefault();
          const residueIndex = residueFromPaePointer(event);
          zoomPlot('pae', event.deltaY < 0 ? 0.75 : 1.35, residueIndex ?? selectionCenter());
        }}, {{ passive: false }});

        document.addEventListener('mouseup', () => {{
          if (isDragging) {{
            isDragging = false;
            dragStartResidue = null;
          }}
          if (plddtDragStart !== null) {{
            if (!plddtDragMoved) {{
              selectResidue(plddtDragStart);
              anchorResidue = plddtDragStart;
            }}
            plddtDragStart = null;
            plddtDragMoved = false;
          }}
          if (paeDragStart !== null) {{
            if (!paeDragMoved) {{
              selectResidue(paeDragStart);
              anchorResidue = paeDragStart;
            }}
            paeDragStart = null;
            paeDragMoved = false;
          }}
          if (plotResize) {{
            if (plotResize.source === 'plddt') updatePlotCursor(plddtSvg, null);
            if (plotResize.source === 'pae') updatePlotCursor(paeCanvas, null);
            plotResize = null;
          }}
        }});

        document.addEventListener('click', (event) => {{
          const button = event.target.closest('.construct-highlight');
          if (!button) return;
          event.preventDefault();
          focusRange(button.dataset.start, button.dataset.end, button.dataset.label);
          const builderSection = document.getElementById('interactive-construct-builder');
          if (builderSection) {{
            builderSection.scrollIntoView({{ behavior: 'smooth', block: 'start' }});
            builderSection.focus({{ preventScroll: true }});
          }}
        }});

        resetButton.addEventListener('click', () => {{
          selectedResidue = null;
          selectedRange = null;
          selectedResidueSet = null;
          showSelectionOnly = false;
          resetPlotViews();
          renderSelection();
          if (viewer) {{ viewer.zoomTo(); viewer.render(); }}
          setInfo('Reset to full-structure view.');
        }});

        ectoButton.addEventListener('click', () => {{
          if (ectodomain.start && ectodomain.end) {{
            focusRange(ectodomain.start, ectodomain.end, designRegionLabel);
            return;
          }}
          const extracellularPositions = extracellularTopologyPositions();
          if (extracellularPositions.length) {{
            focusResidueSet(extracellularPositions, 'Extracellular topology residues');
            return;
          }}
          setInfo('No extracellular region or extracellular topology residues are available for focusing.');
        }});

        selectionOnlyButton.addEventListener('click', () => {{
          if (!(selectedRange || selectedResidue !== null || (selectedResidueSet && selectedResidueSet.positions.length))) {{
            showSelectionOnly = false;
            updateSelectionOnlyControl();
            return;
          }}
          showSelectionOnly = !showSelectionOnly;
          renderSelection();
          if (viewer && showSelectionOnly) {{
            const residues = selectedStructureResidues();
            if (residues.length) viewer.zoomTo({{ resi: residues }});
          }}
          if (viewer) viewer.render();
          setInfo(showSelectionOnly ? 'Showing only the selected region.' : 'Showing the full structure.');
        }});

        if (downloadPdbButton) {{
          downloadPdbButton.addEventListener('click', downloadFullPdb);
        }}

        if (downloadSelectedPdbButton) {{
          downloadSelectedPdbButton.addEventListener('click', downloadSelectedPdb);
        }}

        if (screenshotPngButton) {{
          screenshotPngButton.addEventListener('click', downloadViewerScreenshot);
        }}

        copySelectionTsv.addEventListener('click', async () => {{
          if (!selectionTsv.value) return;
          try {{
            await navigator.clipboard.writeText(selectionTsv.value);
            window.trackOpenAntigensEvent('Copy TSV');
            selectionSummary.textContent = 'Selection table copied as TSV.';
          }} catch (error) {{
            selectionSummary.textContent = 'Copy failed; the TSV is still available in the text box below.';
          }}
        }});

        copySelectionFasta.addEventListener('click', async () => {{
          if (!selectionFasta.value) return;
          try {{
            await navigator.clipboard.writeText(selectionFasta.value);
            window.trackOpenAntigensEvent('Copy FASTA');
            selectionSummary.textContent = 'Selection table copied as FASTA.';
          }} catch (error) {{
            selectionSummary.textContent = 'Copy failed; the FASTA is still available in the text box below.';
          }}
        }});

        plddtZoomIn.addEventListener('click', () => zoomPlot('plddt', 0.75, selectionCenter()));
        plddtZoomOut.addEventListener('click', () => zoomPlot('plddt', 1.35, selectionCenter()));
        plddtZoomReset.addEventListener('click', () => setPlotView('plddt', plotDomain.start, plotDomain.end));
        paeZoomIn.addEventListener('click', () => zoomPlot('pae', 0.75, selectionCenter()));
        paeZoomOut.addEventListener('click', () => zoomPlot('pae', 1.35, selectionCenter()));
        paeZoomReset.addEventListener('click', () => {{
          if (paeView) setPlotView('pae', plotDomain.start, plotDomain.end);
        }});

        const initializeCards = () => document.querySelectorAll('.construct-card').forEach(renderConstructMutationCard);
        if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', initializeCards, {{ once: true }});
        else initializeCards();
        renderTopologyLegend(plotTopologyLegend);
        renderTopologyLegend(paeTopologyLegend);
        updateSelectionOnlyControl();
        scheduleConstructWorkbench();
        renderInteractivePlots();
        renderSelectionTable();

        function initStructureViewer() {{
          if (viewerInitStarted) return;
          viewerInitStarted = true;
          attemptStructureViewerInit(0);
        }}

        function attemptStructureViewerInit(attempt) {{
          const mol3d = window.$3Dmol || window['3Dmol'];
          if (!mol3d || typeof mol3d.createViewer !== 'function') {{
            // 3Dmol may still be loading (e.g. the async CDN fallback when the
            // bundled library is absent). Retry briefly rather than failing
            // permanently; give up after ~10s.
            if (attempt < 40) {{
              setInfo('Loading the structure viewer library…');
              setTimeout(() => attemptStructureViewerInit(attempt + 1), 250);
            }} else {{
              setInfo('Unable to load the AlphaFold structure: 3Dmol viewer library did not load.');
            }}
            return;
          }}
          try {{
            viewer = mol3d.createViewer(viewerEl, {{ backgroundColor: 'white', antialias: true, preserveDrawingBuffer: true }});
            viewer.addModel(pdbText, 'pdb');
            applyBaseStyle(viewer);
            viewer.setClickable({{}}, true, (atom) => {{
              if (atom && atom.resi !== undefined) {{
                selectResidue(mapStructureResidueToSequence(atom.resi));
              }}
            }});
            viewer.zoomTo();
            viewer.render();
            renderSelection();
            setInfo('Interactive AlphaFold viewer ready. Click a residue or a construct button.');
          }} catch (error) {{
            setInfo(`Unable to load the AlphaFold structure: ${{error}}`);
          }}
        }}

        // Defer the WebGL viewer (3Dmol context + model parse + render) off the
        // initial page-load path: initialize when the viewer scrolls near the
        // viewport, or on first interaction with its section (capture phase, so
        // init runs before the triggering control's handler). Every viewer.*
        // call site guards on a null viewer, so any earlier interaction degrades
        // to a no-op rather than an error. Immediate init where IntersectionObserver
        // is unavailable preserves the original eager behavior.
        const viewerSection = viewerEl ? (viewerEl.closest('section') || viewerEl) : null;
        if (viewerEl && 'IntersectionObserver' in window) {{
          const viewerObserver = new IntersectionObserver((entries, observer) => {{
            if (entries.some((entry) => entry.isIntersecting)) {{
              observer.disconnect();
              initStructureViewer();
            }}
          }}, {{ rootMargin: '400px' }});
          viewerObserver.observe(viewerEl);
          if (viewerSection) {{
            ['pointerdown', 'focusin'].forEach((eventName) => {{
              viewerSection.addEventListener(eventName, initStructureViewer, {{ once: true, capture: true }});
            }});
          }}
        }} else {{
          initStructureViewer();
        }}
      }})();
    </script>
    """


def _render_sequence_rows(sequence: str, *, ectodomain: dict[str, Any], width: int = 50) -> str:
    ecto_start = ectodomain.get("start")
    ecto_end = ectodomain.get("end")
    rows: list[str] = []
    for row_start in range(0, len(sequence), width):
        row_end = min(len(sequence), row_start + width)
        residues = []
        for index, residue in enumerate(sequence[row_start:row_end], start=row_start + 1):
            classes = ["residue"]
            if ecto_start is not None and ecto_end is not None and ecto_start <= index <= ecto_end:
                classes.append("ecto-residue")
            residues.append(
                f'<button type="button" class="{" ".join(classes)}" data-resi="{index}" title="{residue}{index}">{escape(residue)}</button>'
            )
        rows.append(
            '<div class="sequence-row"><span class="sequence-index">{start:>4}</span><div class="sequence-residues">{residues}</div></div>'.format(
                start=row_start + 1,
                residues="".join(residues),
            )
        )
    return "".join(rows)


def _build_index_entry(
    item: dict[str, Any],
    report_data: dict[str, Any] | None,
    report_json_path: Path | None,
    *,
    fetch_literature: bool = True,
    allow_af3_structures: bool = False,
) -> dict[str, Any]:
    entry = dict(item)
    inferred_gene_symbol = None
    if report_data is not None:
        inferred_gene_symbol = (report_data.get("target") or {}).get("gene_symbol")
    if not inferred_gene_symbol:
        query_text = str(item.get("query") or "").strip()
        inferred_gene_symbol = query_text.split("_", 1)[0] if "_" in query_text else query_text
    pubtator_info = None
    if fetch_literature:
        from .module_refresh import _module_refresh_config

        batch_dir = report_json_path.parent if report_json_path else Path("outputs/surfy_batch")
        cache_dir = _module_refresh_config(batch_dir=batch_dir, verbose=False, enable_complex_portal=False).cache_dir
        pubtator_info = _get_pubtator_literature_info(inferred_gene_symbol, cache_dir=cache_dir)
    pubtator_count = _pubtator_count_or_existing(pubtator_info, item)
    if report_data is None:
        entry["detail_page"] = None
        entry["gene_symbol"] = item.get("query")
        entry["protein_name"] = ""
        entry["family_name"] = ""
        entry["construct_count"] = ""
        entry["cross_reactivity_count"] = ""
        entry["ectodomain"] = ""
        entry["entry_name"] = item.get("resolved_entry_name") or item.get("query")
        entry["has_alphafold_structure"] = False
        entry["structure_source"] = ""
        entry["has_pdb"] = False
        entry["has_any_structure"] = False
        entry["pubtator_count"] = pubtator_count
        entry["pubtator_query_url"] = pubtator_info.get("query_url") if pubtator_info else None
        return entry

    target = report_data.get("target", {})
    ectodomain = report_data.get("ectodomain") or {}
    family_name = _index_family_name(report_data)
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
            allow_af3_structures=allow_af3_structures,
        )
        if accession
        else None
    )
    has_alphafold = alphafold_path is not None and _alphafold_pdb_matches_target_sequence(report_data, alphafold_path)
    has_pdb = bool(report_data.get("experimental_constructs"))
    entry.update(
        {
            "entry_name": entry_name,
            "gene_symbol": target.get("gene_symbol") or item.get("query"),
            "protein_name": target.get("protein_name") or "",
            "target_aliases": " ".join(_target_aliases(target, entry=entry)),
            "family_name": family_name,
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
            "pubtator_count": pubtator_count,
            "pubtator_query_url": pubtator_info.get("query_url") if pubtator_info else None,
        }
    )
    if report_json_path is not None:
        entry["json_report"] = str(report_json_path)
    return entry


def _index_family_name(report_data: dict[str, Any]) -> str:
    canonical_family = report_data.get("canonical_family") or {}
    canonical_name = str(canonical_family.get("name") or "").strip()
    if canonical_name:
        return canonical_name
    family_context = report_data.get("family_context") or {}
    family_names = family_context.get("family_names") or []
    for name in family_names:
        cleaned = str(name or "").strip()
        if cleaned:
            return cleaned
    return ""


def _index_cross_reactivity_count(report_data: dict[str, Any]) -> int:
    topology = report_data.get("topology") or {}
    topology_class = str(topology.get("topology_class") or "")
    if topology_class.startswith("multipass") or report_data.get("ectodomain") is None:
        return len(report_data.get("full_length_cross_reactivity_hits") or [])
    return len(report_data.get("cross_reactivity_hits") or [])


def _structure_badge(available: bool, *, label: str) -> str:
    badge_class = "yes" if available else "no"
    text = "Yes" if available else "No"
    return f'<span class="structure-badge {badge_class}" aria-label="{escape(label)} {escape(text)}">{escape(text)}</span>'


def _render_pubtator_link(count: Any, query_url: str | None) -> str:
    numeric_count = _coerce_pubtator_count(count)
    label = f"{numeric_count:,}" if numeric_count is not None else "n/a"
    if query_url:
        return f'<a class="inline-link pubtator-link" href="{escape(query_url)}" target="_blank" rel="noreferrer">{escape(label)}</a>'
    return f'<span class="muted">{escape(label)}</span>'


def _species_code(species: str | None) -> str:
    normalized = str(species or "").lower()
    if normalized == "human":
        return "HUMAN"
    if normalized == "mouse":
        return "MOUSE"
    if normalized == "macaca_fascicularis":
        return "MACFA"
    cleaned = "".join(ch if ch.isalnum() else "_" for ch in str(species or "species").upper()).strip("_")
    return cleaned or "SPECIES"


def _species_display_label(species: str | None) -> str:
    normalized = str(species or "").strip().lower().replace(" ", "_")
    if normalized in {"human", "homo_sapiens"}:
        return "human"
    if normalized in {"mouse", "mus_musculus"}:
        return "mouse"
    if normalized in {"macaca_fascicularis", "macfa", "macaque", "cynomolgus_monkey"}:
        return "cynomolgus monkey"
    return str(species or "")


def _target_species_key(target: dict[str, Any]) -> str:
    taxon = str(target.get("taxon_id") or "").strip()
    organism = str(target.get("organism") or "").lower()
    entry_name = str(target.get("entry_name") or "").upper()
    if taxon == "10090" or "mus musculus" in organism or entry_name.endswith("_MOUSE"):
        return "mouse"
    if taxon == "9541" or "macaca fascicularis" in organism or entry_name.endswith("_MACFA"):
        return "macaca_fascicularis"
    return "human"


def _identity_label(report: dict[str, Any]) -> str:
    metadata = report.get("metadata") if isinstance(report.get("metadata"), dict) else {}
    label = str(metadata.get("construct_identity_label") or "").strip()
    return label or "Identity to human"


def _format_export_name(gene_symbol: str | None, species: str | None, start: int | None, end: int | None) -> str:
    gene_token = "".join(ch if ch.isalnum() else "_" for ch in str(gene_symbol or "GENE").upper()).strip("_") or "GENE"
    start_text = str(start or "")
    end_text = str(end or "")
    return f"{gene_token}_{_species_code(species)}_{start_text}-{end_text}"


def _protein_resource_url(accession: str | None, entry_name: str | None) -> str | None:
    accession_text = str(accession or "").strip()
    entry_text = str(entry_name or "").strip()
    if accession_text.startswith(("NP_", "XP_", "YP_", "WP_")):
        return f"https://www.ncbi.nlm.nih.gov/protein/{quote(accession_text)}"
    if accession_text or entry_text:
        return f"https://www.uniprot.org/uniprotkb/{quote(accession_text or entry_text)}"
    return None


def _uniprot_resource_url(
    accession: str | None,
    entry_name: str | None,
    gene_symbol: str | None,
    species: str | None,
) -> str | None:
    accession_text = str(accession or "").strip()
    entry_text = str(entry_name or "").strip()
    if accession_text and not accession_text.startswith(("NP_", "XP_", "YP_", "WP_")):
        return f"https://www.uniprot.org/uniprotkb/{quote(accession_text)}"
    if entry_text and "_" in entry_text and not accession_text.startswith(("NP_", "XP_", "YP_", "WP_")):
        return f"https://www.uniprot.org/uniprotkb/{quote(entry_text)}"
    symbol = str(gene_symbol or "").strip()
    normalized = str(species or "").lower()
    taxon_map = {
        "human": "9606",
        "mouse": "10090",
        "macaca_fascicularis": "9541",
    }
    taxon_id = taxon_map.get(normalized)
    if symbol and taxon_id:
        query = f"(gene_exact:{symbol}) AND (organism_id:{taxon_id})"
        return f"https://www.uniprot.org/uniprotkb?query={quote(query)}"
    return None


def _protein_resource_label(accession: str | None) -> str:
    accession_text = str(accession or "").strip()
    if accession_text.startswith(("NP_", "XP_", "YP_", "WP_")):
        return "RefSeq"
    return "UniProt"


def _gene_resource_url(gene_symbol: str | None, species: str | None) -> str | None:
    symbol = str(gene_symbol or "").strip()
    if not symbol:
        return None
    normalized = str(species or "").lower()
    if normalized == "human":
        return f"https://www.genenames.org/tools/search/#!/?query={quote(symbol)}"
    taxon_map = {
        "mouse": "10090",
        "macaca_fascicularis": "9541",
    }
    taxon_id = taxon_map.get(normalized)
    if taxon_id:
        query = f"{symbol}[Gene Name] AND {taxon_id}[Taxonomy ID]"
        return f"https://www.ncbi.nlm.nih.gov/gene/?term={quote(query)}"
    return None


def _gene_resource_label(species: str | None) -> str:
    return "HGNC" if str(species or "").lower() == "human" else "Gene"


def _citeab_resource_url(gene_symbol: str | None, species: str | None) -> str | None:
    symbol = str(gene_symbol or "").strip()
    normalized_species = str(species or "").lower()
    if not symbol or normalized_species not in {"", "human", "homo sapiens"}:
        return None
    return f"https://www.citeab.com/antibodies/search?q={quote(symbol)}"


def _render_resource_links(
    *,
    accession: str | None,
    entry_name: str | None,
    gene_symbol: str | None,
    species: str | None,
) -> str:
    links: list[str] = []
    protein_url = _protein_resource_url(accession, entry_name)
    if protein_url:
        links.append(
            f'<a class="inline-link" href="{escape(protein_url)}" target="_blank" rel="noreferrer">{escape(_protein_resource_label(accession))}</a>'
        )
    uniprot_url = _uniprot_resource_url(accession, entry_name, gene_symbol, species)
    if uniprot_url and (not protein_url or "uniprot.org/uniprotkb/" not in protein_url):
        links.append(
            f'<a class="inline-link" href="{escape(uniprot_url)}" target="_blank" rel="noreferrer">UniProt</a>'
        )
    gene_url = _gene_resource_url(gene_symbol, species)
    if gene_url:
        links.append(
            f'<a class="inline-link" href="{escape(gene_url)}" target="_blank" rel="noreferrer">{escape(_gene_resource_label(species))}</a>'
        )
    citeab_url = _citeab_resource_url(gene_symbol, species)
    if citeab_url:
        links.append(
            f'<a class="inline-link" href="{escape(citeab_url)}" target="_blank" rel="noreferrer">CiteAb</a>'
        )
    return " ".join(links) if links else '<span class="muted">n/a</span>'


def _construct_display_group(construct: dict[str, Any]) -> tuple[int, int, int, str]:
    name = str(construct.get("name") or "")
    start = int(construct.get("start") or 0)
    end = int(construct.get("end") or 0)
    membrane_order = {
        "membrane_expression_full_length": 0,
        "membrane_expression_c_tail_trimmed": 1,
        "membrane_expression_n_tail_trimmed": 2,
        "membrane_expression_terminal_trimmed": 3,
    }
    if name in membrane_order:
        return (5, membrane_order[name], start, name)
    if name.startswith("membrane_expression_gpcr_c_tail_keep_"):
        return (5, 4, start, name)
    if name == "full_ectodomain":
        return (0, start, end, name)
    if name.startswith("pdb_"):
        return (1, start, end, name)
    if name.startswith("strict_") or name.startswith("structural_region_"):
        return (3, start, end, name)
    if name.startswith("lenient_"):
        return (4, start, end, name)
    return (2, start, end, name)


def _sorted_constructs_for_display(constructs: list[dict[str, Any]]) -> list[dict[str, Any]]:
    def summary_order(construct: dict[str, Any]) -> tuple[int, int, int, str]:
        group, *boundary = _construct_display_group(construct)
        # Match the builder while keeping the Construct Details tab categories unchanged.
        return ((0, 1, 4, 2, 3, 5)[group], *boundary)

    return sorted(constructs, key=summary_order)


def _partition_construct_dicts(constructs: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    soluble: list[dict[str, Any]] = []
    membrane: list[dict[str, Any]] = []
    for construct in constructs:
        if str(construct.get("classification") or "") == "membrane_expression":
            membrane.append(construct)
        else:
            soluble.append(construct)
    return soluble, membrane


def _construct_group_label(construct: dict[str, Any]) -> str:
    name = str(construct.get("name") or "")
    classification = str(construct.get("classification") or "")
    evidence_text = " ".join(str(item) for item in (construct.get("evidence") or []))
    evidence_lower = evidence_text.lower()
    name_lower = name.lower()
    if name == "membrane_expression_full_length":
        return "Full-length membrane construct"
    if name in {
        "membrane_expression_c_tail_trimmed",
        "membrane_expression_n_tail_trimmed",
        "membrane_expression_terminal_trimmed",
    } or name.startswith("membrane_expression_gpcr_c_tail_keep_"):
        return "Terminal-tail truncation"
    if name == "full_ectodomain":
        return "Topology"
    if name.startswith("pdb_"):
        return "PDB-backed"
    if "alphafold strict" in evidence_lower or name.startswith("strict_") or name.startswith("structural_region_"):
        return "AlphaFold strict"
    if "alphafold lenient" in evidence_lower or name.startswith("lenient_"):
        return "AlphaFold lenient"
    if classification == "membrane_expression":
        return "Membrane expression"
    if "topology-derived" in evidence_lower:
        return "Topology"
    if "interpro" in evidence_lower:
        if classification == "repeat_module" or "repeat" in name_lower or "repeat" in evidence_lower:
            return "InterPro repeat"
        return "InterPro domain"
    if "uniprot" in evidence_lower:
        return "UniProt domain"
    if classification == "repeat_module":
        return "Repeat module"
    if classification in {"complete_domain", "multi_domain_unit"}:
        return "Domain"
    return "Construct"


def _construct_support_summary(construct: dict[str, Any]) -> str:
    evidence = [str(item).strip() for item in (construct.get("evidence") or []) if str(item).strip()]
    rationale = str(construct.get("rationale") or "").strip()
    parts: list[str] = []
    if evidence:
        parts.append(_scientific_region_text(evidence[0]))
    if rationale and rationale not in parts:
        parts.append(_scientific_region_text(rationale))
    return " | ".join(parts[:2]) or "n/a"


def _scientific_region_text(value: Any) -> str:
    text = str(value or "")
    replacements = (
        (r"\bECD-focused\b", "extracellular-region"),
        (r"\bECD\b", "extracellular region"),
        (r"\b[Ee]ctodomain\b", "design region"),
    )
    for pattern, replacement in replacements:
        text = re.sub(pattern, replacement, text)
    return text


def _rcsb_structure_url(pdb_id: str) -> str:
    return f"https://www.rcsb.org/structure/{quote(pdb_id.strip().upper())}"


def _construct_name_html(construct: dict[str, Any]) -> str:
    """Escaped construct name, linked to its RCSB entry when PDB-backed."""
    name = str(construct.get("name") or "")
    pdb_id = str(construct.get("pdb_id") or "").strip()
    if not pdb_id:
        return escape(name)
    return (
        f'<a class="inline-link" href="{escape(_rcsb_structure_url(pdb_id))}" '
        f'target="_blank" rel="noreferrer">{escape(name)}</a>'
    )


# Isometric-cube glyph that reads as "3D" on the viewer call-to-action button.
_VIEW_3D_ICON = (
    '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" '
    'stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">'
    '<path d="M21 7.5 12 12 3 7.5 12 3z"/>'
    '<path d="M3 7.5v9L12 21l9-4.5v-9"/>'
    '<path d="M12 12v9"/></svg>'
)


def _viewer_focus_button(start: Any, end: Any, label: Any, *, compact: bool = False) -> str:
    """Prominent button that focuses a construct's residue range in the 3D viewer.

    Returns an empty string when boundaries are missing so callers can supply
    their own fallback. The ``construct-highlight`` class is the JS click hook.
    """
    if start is None or end is None or str(start).strip() == "" or str(end).strip() == "":
        return ""
    classes = "view-3d-btn construct-highlight" + (" view-3d-btn--compact" if compact else "")
    text = "View in 3D" if compact else "View in 3D structure"
    return (
        f'<button type="button" class="{classes}" '
        f'data-start="{escape(str(start))}" data-end="{escape(str(end))}" '
        f'data-label="{escape(str(label or "construct"))}">'
        f'{_VIEW_3D_ICON}<span>{escape(text)}</span></button>'
    )


def _render_precomputed_construct_summary_rows(
    constructs: list[dict[str, Any]],
    *,
    empty_message: str,
    allow_structure_actions: bool = True,
) -> str:
    if not constructs:
        return f"<tr><td colspan=\"6\">{escape(empty_message)}</td></tr>"
    rows: list[str] = []
    for item in constructs:
        start = item.get("start")
        end = item.get("end")
        boundary = f"{start}-{end}" if start is not None and end is not None else ""
        action = (
            _viewer_focus_button(start, end, item.get("name"), compact=True)
            if allow_structure_actions
            else ""
        )
        rows.append(
            "<tr><td><code>{name}</code></td><td>{group}</td><td>{boundary}</td><td>{length}</td><td>{support}</td><td>{action}</td></tr>".format(
                name=_construct_name_html(item),
                group=escape(_construct_group_label(item)),
                boundary=escape(boundary),
                length=escape(str(item.get("length") or "")),
                support=escape(_construct_support_summary(item)),
                action=action,
            )
        )
    return "".join(rows)


def _fmt_evalue(value: Any) -> str:
    if value is None:
        return ""
    try:
        return f"{float(value):.2e}"
    except Exception:
        return str(value)


def _alignment_metrics_from_strings(aligned_query: str, aligned_subject: str) -> dict[str, Any]:
    matches = 0
    aligned_positions = 0
    query_start = None
    query_end = None
    subject_start = None
    subject_end = None
    query_position = 0
    subject_position = 0
    for query_residue, subject_residue in zip(aligned_query, aligned_subject, strict=True):
        if query_residue != "-":
            query_position += 1
        if subject_residue != "-":
            subject_position += 1
        if query_residue == "-" or subject_residue == "-":
            continue
        aligned_positions += 1
        if query_residue == subject_residue:
            matches += 1
        if query_start is None:
            query_start = query_position
            subject_start = subject_position
        query_end = query_position
        subject_end = subject_position
    identity = (100.0 * matches / aligned_positions) if aligned_positions else None
    return {
        "matches": matches,
        "aligned_positions": aligned_positions,
        "identity": identity,
        "query_start": query_start,
        "query_end": query_end,
        "subject_start": subject_start,
        "subject_end": subject_end,
    }


def _format_alignment_block(aligned_query: str, aligned_subject: str, *, width: int = 70) -> str:
    if not aligned_query or not aligned_subject:
        return ""
    blocks: list[str] = []
    for index in range(0, len(aligned_query), width):
        query_block = aligned_query[index : index + width]
        subject_block = aligned_subject[index : index + width]
        match_block = "".join("|" if q == s and q != "-" else " " for q, s in zip(query_block, subject_block, strict=True))
        blocks.append(f"human  {query_block}")
        blocks.append(f"       {match_block}")
        blocks.append(f"hit    {subject_block}")
        blocks.append("")
    return "\n".join(blocks).rstrip()


def _format_alignment_block_html(
    aligned_query: str,
    aligned_subject: str,
    query_annotation: str | None,
    *,
    width: int = 70,
    query_label: str = "human",
) -> str:
    if not aligned_query or not aligned_subject:
        return ""
    annotation = query_annotation or "." * len(aligned_query)
    blocks: list[str] = []
    for index in range(0, len(aligned_query), width):
        query_block = aligned_query[index : index + width]
        subject_block = aligned_subject[index : index + width]
        annotation_block = annotation[index : index + width]
        match_block = "".join("|" if q == s and q != "-" else " " for q, s in zip(query_block, subject_block, strict=True))
        blocks.append(
            f"{query_label[:5]:<5}  "
            + _annotated_alignment_sequence_html(query_block, annotation_block)
        )
        blocks.append("       " + escape(match_block))
        blocks.append("hit    " + escape(subject_block))
        blocks.append("")
    return "\n".join(blocks).rstrip()


def _annotated_alignment_sequence_html(sequence: str, annotation: str) -> str:
    pieces: list[str] = []
    for residue, residue_class in zip(sequence, annotation, strict=False):
        text = escape(residue)
        if residue_class == "S" and residue != "-":
            pieces.append(f'<span class="alignment-residue alignment-surface-ecto">{text}</span>')
        else:
            pieces.append(text)
    return "".join(pieces)


def _render_cross_reactivity_alignment(
    hit: dict[str, Any],
    *,
    rank: int,
    query_label: str = "human",
) -> str:
    query_alignment = str(hit.get("query_alignment") or "")
    subject_alignment = str(hit.get("subject_alignment") or "")
    if query_alignment and subject_alignment:
        alignment_source = str(hit.get("alignment_source") or "BLAST")
        query_start = hit.get("query_start")
        query_end = hit.get("query_end")
        subject_start = hit.get("subject_start")
        subject_end = hit.get("subject_end")
        range_text = "n/a"
        if None not in (query_start, query_end, subject_start, subject_end):
            range_text = f"{query_label} {query_start}-{query_end} | hit {subject_start}-{subject_end}"
        return (
            f"<div class=\"alignment-meta\"><span><strong>Rank:</strong> {escape(str(rank))}</span>"
            f"<span><strong>Bit score:</strong> {escape(_fmt_number(hit.get('bitscore')))}</span>"
            f"<span><strong>Alignment source:</strong> {escape(alignment_source)}</span>"
            f"<span><strong>Aligned range:</strong> {escape(range_text)}</span></div>"
            f"<div class=\"alignment-legend\"><span class=\"alignment-residue alignment-surface-ecto\">A</span> {escape(query_label)} extracellular accessible residue</div>"
            f"<pre class=\"alignment-block\">{_format_alignment_block_html(query_alignment, subject_alignment, hit.get('query_alignment_annotation'), query_label=query_label)}</pre>"
        )

    return "<p class=\"viewer-info\">No stored BLAST alignment is available for this hit.</p>"


def _fetch_cross_reactivity_target(hit: dict[str, Any], *, batch_dir: Path) -> dict[str, Any] | None:
    subject_id = str(hit.get("subject_id") or "").strip()
    accession = extract_accession(subject_id)
    if accession is None and subject_id.startswith(("NP_", "XP_", "YP_", "WP_")):
        accession = subject_id
    entry_name = extract_entry_name(subject_id)
    if accession and str(accession).startswith(("NP_", "XP_", "YP_", "WP_")):
        return _fetch_refseq_target(str(accession), batch_dir=batch_dir)
    if accession or entry_name:
        return _fetch_uniprot_target(str(accession or entry_name), batch_dir=batch_dir)
    return None


def _wrap_fasta_sequence(sequence: str, *, width: int = 80) -> str:
    if not sequence:
        return ""
    return "\n".join(sequence[index : index + width] for index in range(0, len(sequence), width))


def _identity_header(label: str) -> str:
    normalized = re.sub(r"[^a-z0-9]+", "_", str(label or "identity").strip().lower()).strip("_")
    return normalized or "identity"


def _surface_identity_header(identity_label: str) -> str:
    return f"extracellular_surface_{_identity_header(identity_label)}"


def _export_tsv(
    rows: list[dict[str, str]],
    *,
    identity_header: str = "identity_to_human",
    surface_identity_header: str = "extracellular_surface_identity_to_human",
) -> str:
    header = ["name", "species", "boundary", "sequence", identity_header, surface_identity_header]
    lines = ["\t".join(header)]
    for row in rows:
        lines.append("\t".join([row.get("name", ""), row.get("species", ""), row.get("boundary", ""), row.get("sequence", ""), row.get("identity", ""), row.get("surface_identity", "")]))
    return "\n".join(lines)


def _export_fasta(
    rows: list[dict[str, str]],
    *,
    identity_header: str = "identity_to_human",
    surface_identity_header: str = "extracellular_surface_identity_to_human",
) -> str:
    records: list[str] = []
    for row in rows:
        header = (
            f">{row.get('name', '')} species={row.get('species', '')} boundary={row.get('boundary', '')} "
            f"{identity_header}={row.get('identity', '') or 'n/a'} "
            f"{surface_identity_header}={row.get('surface_identity', '') or 'n/a'}"
        )
        records.append(f"{header}\n{_wrap_fasta_sequence(row.get('sequence', ''))}".rstrip())
    return "\n".join(record for record in records if record)


def _resolve_existing_path(path_text: str | None, batch_dir: Path) -> Path | None:
    if not path_text:
        return None
    path = Path(path_text)
    if path.is_absolute():
        candidates = [
            path,
            batch_dir / path.name,
            Path.cwd() / path.name,
        ]
    else:
        candidates = [
            batch_dir / path,
            batch_dir / path.name,
            Path.cwd() / path,
        ]
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    return None


def _resolve_alphafold_artifact(
    accession: str | None,
    suffix: str,
    *,
    batch_dir: Path,
    allow_af3_structures: bool = False,
) -> Path | None:
    accession_text = str(accession or "").strip()
    if not accession_text:
        return None
    suffix_text = suffix if suffix.startswith(".") else f".{suffix}"
    file_name = f"{accession_text}{suffix_text}"
    bases: list[Path] = []
    for base in (batch_dir, batch_dir.parent, batch_dir.parent.parent, Path.cwd()):
        if base not in bases:
            bases.append(base)
    candidates: list[Path] = []
    for base in bases:
        candidates.extend(
            [
                base / "alphafold" / file_name,
                base / "data" / "alphafold" / file_name,
                base / ".agdesign2" / "data" / "alphafold" / file_name,
            ]
        )
        if suffix_text == ".pae.json":
            candidates.extend(
                [
                    base / "alphafold" / f"{file_name}.gz",
                    base / "data" / "alphafold" / f"{file_name}.gz",
                    base / ".agdesign2" / "data" / "alphafold" / f"{file_name}.gz",
                ]
            )
    primary = next((candidate.resolve() for candidate in candidates if candidate.exists()), None)
    if primary is not None:
        return primary
    if not allow_af3_structures:
        return None
    candidates = []
    for base in bases:
        candidates.extend(
            [
                base / "alphafold3" / file_name,
                base / "data" / "alphafold3" / file_name,
                base / ".agdesign2" / "data" / "alphafold3" / file_name,
            ]
        )
        if suffix_text == ".pae.json":
            candidates.extend(
                [
                    base / "alphafold3" / f"{file_name}.gz",
                    base / "data" / "alphafold3" / f"{file_name}.gz",
                    base / ".agdesign2" / "data" / "alphafold3" / f"{file_name}.gz",
                ]
            )
    return next((candidate.resolve() for candidate in candidates if candidate.exists()), None)


def _copy_portal_structure(
    *,
    report_data: dict[str, Any],
    structures_dir: Path,
    batch_dir: Path,
    allow_af3_structures: bool = False,
) -> str | None:
    target = report_data.get("target", {})
    accession = target.get("accession")
    entry_name = target.get("entry_name") or accession
    if not accession or not entry_name:
        return None
    destination = structures_dir / f"{_safe_filename(str(entry_name).lower())}.pdb"
    source = _resolve_alphafold_artifact(
        accession,
        ".pdb",
        batch_dir=batch_dir,
        allow_af3_structures=allow_af3_structures,
    )
    destination_matches = destination.exists() and _alphafold_pdb_matches_target_sequence(report_data, destination)
    if source is None:
        if destination_matches:
            return str(destination)
        return None
    if not _alphafold_pdb_matches_target_sequence(report_data, source):
        if destination_matches:
            return str(destination)
        return None
    try:
        if (
            not destination_matches
            or not destination.exists()
            or source.stat().st_mtime > destination.stat().st_mtime
            or source.stat().st_size != destination.stat().st_size
        ):
            temp_destination = destination.with_suffix(destination.suffix + ".tmp")
            shutil.copyfile(source, temp_destination)
            temp_destination.replace(destination)
        source_meta = source.with_suffix(".meta.json")
        if source_meta.exists():
            meta_destination = destination.with_suffix(".meta.json")
            temp_meta_destination = meta_destination.with_suffix(meta_destination.suffix + ".tmp")
            shutil.copyfile(source_meta, temp_meta_destination)
            temp_meta_destination.replace(meta_destination)
    except Exception:
        return None
    if not _alphafold_pdb_matches_target_sequence(report_data, destination):
        return None
    return str(destination)


def _alphafold_pdb_matches_target_sequence(report_data: dict[str, Any], pdb_path: Path) -> bool:
    target_sequence = str((report_data.get("target") or {}).get("sequence") or "").strip()
    if not target_sequence:
        return False
    try:
        stat = pdb_path.stat()
        cache_key = (
            str(pdb_path.resolve()),
            int(stat.st_mtime_ns),
            int(stat.st_size),
            hashlib.sha1(
                target_sequence.encode("ascii", "ignore"), usedforsecurity=False
            ).hexdigest(),
        )
    except Exception:
        cache_key = None
    if cache_key is not None and cache_key in _ALPHAFOLD_SEQUENCE_MATCH_CACHE:
        return _ALPHAFOLD_SEQUENCE_MATCH_CACHE[cache_key]
    try:
        structure = parse_alphafold_pdb(pdb_path)
    except Exception:
        return False
    residue_names = structure.get("residue_names") or {}
    if not residue_names:
        return False
    structure_sequence = "".join(
        _THREE_TO_ONE.get(str(residue_names[residue_id]).upper(), "X")
        for residue_id in sorted(residue_names)
    )
    matches = structure_sequence == target_sequence
    if cache_key is not None:
        _ALPHAFOLD_SEQUENCE_MATCH_CACHE[cache_key] = matches
    return matches


def _structure_source_label(pdb_path: Path) -> str:
    meta_path = pdb_path.with_suffix(".meta.json")
    if meta_path.exists():
        try:
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
        except Exception:
            meta = {}
        if meta.get("source") == "alphafold_server" or meta.get("model_type") == "AF3":
            return "AlphaFold 3 model"
    if is_af3_artifact(pdb_path):
        return "AlphaFold 3 model"
    return "AlphaFold"


def _report_uses_af3_structure(report_data: dict[str, Any]) -> bool:
    for note in report_data.get("notes") or []:
        if not isinstance(note, dict):
            continue
        source = str(note.get("source") or "").strip().casefold()
        message = str(note.get("message") or "").strip().casefold()
        if source in {"alphafold server", "alphafold 3", "af3"}:
            return True
        if "alphafold 3" in message:
            return True
    metadata = report_data.get("metadata")
    if isinstance(metadata, dict):
        source = str(metadata.get("structure_source") or "").strip().casefold()
        model_type = str(metadata.get("model_type") or "").strip().casefold()
        if source in {"alphafold server", "alphafold 3", "af3", "alphafold_server"}:
            return True
        if model_type == "af3":
            return True
    return False


def _copy_report_assets(
    *,
    report_data: dict[str, Any],
    report_json_path: Path | None,
    report_assets_dir: Path,
    batch_dir: Path,
) -> Path | None:
    target = report_data.get("target") or {}
    entry_name = str(target.get("entry_name") or "").strip().lower()
    if not entry_name:
        return None
    asset_dir_name = f"{entry_name}_report_assets"
    candidates = [
        report_json_path.parent / asset_dir_name if report_json_path is not None else None,
        batch_dir / asset_dir_name,
        Path.cwd() / asset_dir_name,
    ]
    source = next((candidate.resolve() for candidate in candidates if candidate is not None and candidate.exists()), None)
    if source is None or not source.is_dir():
        return None
    symlinks = [path.relative_to(source) for path in source.rglob("*") if path.is_symlink()]
    if symlinks:
        raise ValueError(f"Report asset directory contains symbolic links: {', '.join(map(str, symlinks))}")
    destination = report_assets_dir / source.name
    try:
        if destination.exists():
            shutil.rmtree(destination)
        shutil.copytree(source, destination, symlinks=True)
    except OSError as exc:
        raise RuntimeError(f"Could not copy report assets from {source}") from exc
    return destination


def _load_viewer_payload(
    *,
    report: dict[str, Any],
    batch_dir: Path,
    pdb_path: Path | None = None,
    allow_af3_structures: bool = False,
) -> dict[str, Any]:
    target = report.get("target", {})
    ectodomain = report.get("ectodomain") or {}
    topology = report.get("topology") or {}
    topology_class = str(topology.get("topology_class") or "")
    plot_scope_full_length = topology_class.startswith("multipass")
    accession = target.get("accession")
    target_sequence = str(target.get("sequence") or "")
    if not accession:
        return {}
    pdb_path = pdb_path or _resolve_alphafold_artifact(
        accession,
        ".pdb",
        batch_dir=batch_dir,
        allow_af3_structures=allow_af3_structures,
    )
    pae_path = _resolve_alphafold_artifact(
        accession,
        ".pae.json",
        batch_dir=batch_dir,
        allow_af3_structures=allow_af3_structures,
    )
    if pdb_path is not None and not _alphafold_pdb_matches_target_sequence(report, pdb_path):
        pdb_path = None
        pae_path = None
    payload: dict[str, Any] = {}
    if pdb_path is not None:
        try:
            structure = parse_alphafold_pdb(pdb_path)
            residue_plddt = structure.get("plddt") or {}
            max_residue = max(residue_plddt) if residue_plddt else len(target.get("sequence") or "")
            payload["plddt"] = [round(float(residue_plddt.get(index, 0.0)), 2) for index in range(1, max_residue + 1)]
        except Exception:
            pass
    if pdb_path is not None and target_sequence:
        try:
            payload.update(_build_sequence_structure_mapping(target_sequence, pdb_path))
        except Exception:
            pass
    if pae_path is not None:
        try:
            pae_matrix = load_pae_matrix(pae_path)
            if not pae_matches_sequence_length(pae_matrix, len(target_sequence)):
                raise ValueError(
                    "AlphaFold PAE dimensions do not match the canonical report sequence."
                )
            if plot_scope_full_length:
                submatrix = pae_matrix
            elif ectodomain.get("start") is not None and ectodomain.get("end") is not None:
                ecto_start = int(ectodomain["start"])
                ecto_end = int(ectodomain["end"])
                submatrix = [row[ecto_start - 1 : ecto_end] for row in pae_matrix[ecto_start - 1 : ecto_end]]
            else:
                submatrix = pae_matrix
            downsampled = _downsample_matrix(submatrix, max_size=220)
            payload["paeMatrix"] = downsampled
            payload["paeMax"] = round(
                max((max(row) for row in downsampled if row), default=30.0),
                2,
            )
        except Exception:
            pass
    payload["cysteineAnalysis"] = report.get("cysteine_analysis") or []
    payload["residueAnnotations"] = report.get("residue_annotations") or []
    payload["furinSites"] = report.get("furin_sites") or []
    payload["ptms"] = _report_ptms(report)
    payload["gpcrAnnotation"] = report.get("gpcr_annotation") or None
    payload["homologMappings"] = _load_homolog_mappings(report=report, batch_dir=batch_dir)
    return payload


def _downsample_matrix(matrix: list[list[float]], *, max_size: int) -> list[list[float]]:
    if not matrix:
        return []
    size = len(matrix)
    if size <= max_size:
        return [[round(float(value), 2) for value in row] for row in matrix]
    step = size / max_size
    reduced: list[list[float]] = []
    for row_index in range(max_size):
        row_start = int(row_index * step)
        row_end = max(row_start + 1, int((row_index + 1) * step))
        reduced_row: list[float] = []
        for col_index in range(max_size):
            col_start = int(col_index * step)
            col_end = max(col_start + 1, int((col_index + 1) * step))
            total = 0.0
            count = 0
            for source_row in range(row_start, min(row_end, size)):
                for source_col in range(col_start, min(col_end, size)):
                    total += float(matrix[source_row][source_col])
                    count += 1
            reduced_row.append(round(total / max(count, 1), 2))
        reduced.append(reduced_row)
    return reduced


def _load_homolog_mappings(*, report: dict[str, Any], batch_dir: Path) -> list[dict[str, Any]]:
    target = report.get("target", {})
    ectodomain = report.get("ectodomain") or {}
    topology = report.get("topology") or {}
    topology_class = str(topology.get("topology_class") or "")
    use_full_length = topology_class.startswith("multipass")
    human_sequence = str(target.get("sequence") or "")
    if not human_sequence:
        return []
    if use_full_length:
        human_region_start = 1
        human_region_end = len(human_sequence)
    elif ectodomain.get("start") is not None and ectodomain.get("end") is not None:
        human_region_start = int(ectodomain["start"])
        human_region_end = int(ectodomain["end"])
    elif any(record.get("available") and record.get("ectodomain_start") is not None and record.get("ectodomain_end") is not None for record in report.get("ectodomain_homology") or []):
        human_region_start = 1
        human_region_end = len(human_sequence)
    else:
        return []
    human_region_sequence = slice_sequence(human_sequence, human_region_start, human_region_end)
    mappings: list[dict[str, Any]] = []
    for record in report.get("ectodomain_homology") or []:
        if not record.get("available"):
            continue
        homolog_target = _fetch_homolog_target(record, batch_dir=batch_dir)
        if not homolog_target:
            continue
        homolog_sequence = str(homolog_target.get("sequence") or "")
        if not homolog_sequence:
            continue
        if use_full_length:
            homolog_region_start = 1
            homolog_region_end = len(homolog_sequence)
        else:
            homolog_ecto_start = record.get("ectodomain_start")
            homolog_ecto_end = record.get("ectodomain_end")
            if homolog_ecto_start is None or homolog_ecto_end is None:
                continue
            homolog_region_start = int(homolog_ecto_start)
            homolog_region_end = int(homolog_ecto_end)
        homolog_region_sequence = slice_sequence(homolog_sequence, homolog_region_start, homolog_region_end)
        if not homolog_region_sequence:
            continue
        alignment = _portal_align_sequences(human_region_sequence, homolog_region_sequence)
        mappings.append(
            {
                "species": record.get("species"),
                "geneSymbol": homolog_target.get("gene_symbol") or record.get("entry_name") or record.get("accession"),
                "entryName": homolog_target.get("entry_name") or record.get("entry_name") or record.get("accession"),
                "accession": homolog_target.get("accession") or record.get("accession"),
                "humanRegionStart": human_region_start,
                "humanRegionEnd": human_region_end,
                "homologRegionStart": homolog_region_start,
                "homologRegionEnd": homolog_region_end,
                "ectodomainStart": homolog_region_start,
                "ectodomainEnd": homolog_region_end,
                "ectodomainSequence": homolog_region_sequence,
                "alignedHumanEctodomain": alignment["aligned_query"],
                "alignedHomologEctodomain": alignment["aligned_subject"],
            }
        )
    existing_species = {str(item.get("species") or "").lower() for item in mappings}
    for fallback_mapping in _load_homolog_mappings_from_constructs(
        report=report,
        human_region_start=human_region_start,
        human_region_end=human_region_end,
        human_region_sequence=human_region_sequence,
    ):
        species_key = str(fallback_mapping.get("species") or "").lower()
        if species_key and species_key in existing_species:
            continue
        mappings.append(fallback_mapping)
        existing_species.add(species_key)
    return mappings


def _load_homolog_mappings_from_constructs(
    *,
    report: dict[str, Any],
    human_region_start: int,
    human_region_end: int,
    human_region_sequence: str,
) -> list[dict[str, Any]]:
    constructs = report.get("construct_details") or []
    candidates: list[dict[str, Any]] = []
    for construct in constructs:
        if int(construct.get("start") or 0) != human_region_start:
            continue
        if int(construct.get("end") or 0) != human_region_end:
            continue
        if not construct.get("homologs"):
            continue
        candidates.append(construct)
    if not candidates:
        for construct in constructs:
            name = str(construct.get("name") or "")
            classification = str(construct.get("classification") or "")
            if name == "full_ectodomain" or classification == "full_ectodomain":
                if construct.get("homologs"):
                    candidates.append(construct)
                    break
    if not candidates:
        return []

    construct = max(candidates, key=lambda item: len(str(item.get("sequence") or "")))
    construct_sequence = str(construct.get("sequence") or human_region_sequence)
    if construct_sequence and construct_sequence != human_region_sequence:
        human_alignment_reference = construct_sequence
    else:
        human_alignment_reference = human_region_sequence

    mappings: list[dict[str, Any]] = []
    for homolog in construct.get("homologs") or []:
        if not homolog.get("available"):
            continue
        homolog_sequence = str(homolog.get("sequence") or "")
        homolog_start = homolog.get("start")
        homolog_end = homolog.get("end")
        if not homolog_sequence or homolog_start is None or homolog_end is None:
            continue
        alignment = _portal_align_sequences(human_alignment_reference, homolog_sequence)
        mappings.append(
            {
                "species": homolog.get("species"),
                "geneSymbol": report.get("target", {}).get("gene_symbol") or homolog.get("entry_name") or homolog.get("accession"),
                "entryName": homolog.get("entry_name") or homolog.get("accession"),
                "accession": homolog.get("accession"),
                "identityToHuman": homolog.get("identity_to_human"),
                "humanRegionStart": int(construct.get("start") or human_region_start),
                "humanRegionEnd": int(construct.get("end") or human_region_end),
                "homologRegionStart": int(homolog_start),
                "homologRegionEnd": int(homolog_end),
                "ectodomainStart": int(homolog_start),
                "ectodomainEnd": int(homolog_end),
                "ectodomainSequence": homolog_sequence,
                "alignedHumanEctodomain": alignment["aligned_query"],
                "alignedHomologEctodomain": alignment["aligned_subject"],
            }
        )
    return mappings


def _build_sequence_structure_mapping(target_sequence: str, pdb_path: Path) -> dict[str, Any]:
    structure = parse_alphafold_pdb(pdb_path)
    residue_names = structure.get("residue_names") or {}
    residue_ids = sorted(int(residue_id) for residue_id in residue_names)
    if not residue_ids or not target_sequence:
        return {}
    structure_sequence = "".join(_THREE_TO_ONE.get(str(residue_names[residue_id]).upper(), "X") for residue_id in residue_ids)
    alignment = _portal_align_sequences(target_sequence, structure_sequence)
    sequence_to_structure: list[int | None] = [None] * len(target_sequence)
    structure_to_sequence: dict[str, int] = {}
    query_position = 0
    subject_position = 0
    for query_residue, subject_residue in zip(alignment["aligned_query"], alignment["aligned_subject"], strict=True):
        if query_residue != "-":
            query_position += 1
        if subject_residue != "-":
            subject_position += 1
        if query_residue == "-" or subject_residue == "-":
            continue
        structure_residue = residue_ids[subject_position - 1]
        sequence_to_structure[query_position - 1] = structure_residue
        structure_to_sequence[str(structure_residue)] = query_position
    return {
        "sequenceToStructureResidues": sequence_to_structure,
        "structureToSequenceResidues": structure_to_sequence,
        "structureResidueCount": len(residue_ids),
    }


def _portal_align_sequences(query: str, subject: str) -> dict[str, str]:
    alignment = global_align(query, subject)
    result = {
        "aligned_query": alignment.aligned_query,
        "aligned_subject": alignment.aligned_subject,
    }
    return result


def _fetch_homolog_target(record: dict[str, Any], *, batch_dir: Path) -> dict[str, Any] | None:
    accession = str(record.get("accession") or "")
    entry_name = str(record.get("entry_name") or "")
    if not (accession or entry_name):
        return None
    if accession.startswith(("NP_", "XP_")):
        result = _fetch_refseq_target(accession, batch_dir=batch_dir)
    else:
        result = _fetch_uniprot_target(accession or entry_name, batch_dir=batch_dir)
    # Fall back to the ortholog reference table (keyed by accession across human,
    # mouse and macaque). This is the only source for the human sequence in the
    # mouse portal, where the human entry was never fetched into the live cache
    # (human is the base, never a homolog, in the human portal).
    if result is None and accession:
        result = _fetch_refseq_target_from_ortholog_table(accession, batch_dir=batch_dir)
    return result


def _fetch_uniprot_target(accession_or_entry: str, *, batch_dir: Path) -> dict[str, Any] | None:
    try:
        url = f"https://rest.uniprot.org/uniprotkb/{quote(accession_or_entry)}.json"
        cache_path = _cached_response_path("uniprot_entry", url, ".json", batch_dir=batch_dir)
        if not cache_path.exists():
            return None
        entry = json.loads(cache_path.read_text(encoding="utf-8"))
        if not entry:
            return None
        gene_symbol = None
        genes = entry.get("genes") or []
        if genes:
            gene_name = genes[0].get("geneName")
            if isinstance(gene_name, dict):
                gene_symbol = gene_name.get("value")
        sequence = entry.get("sequence", {}).get("value", "")
        return {
            "accession": entry.get("primaryAccession"),
            "entry_name": entry.get("uniProtkbId"),
            "gene_symbol": gene_symbol,
            "sequence": sequence,
        }
    except Exception:
        return None


def _fetch_refseq_target(accession: str, *, batch_dir: Path) -> dict[str, Any] | None:
    try:
        fasta_url = (
            f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id={quote(accession)}&rettype=fasta&retmode=text"
        )
        cache_path = _cached_response_path("refseq_fasta", fasta_url, ".fasta", batch_dir=batch_dir)
        if not cache_path.exists():
            return _fetch_refseq_target_from_ortholog_table(accession, batch_dir=batch_dir)
        fasta_text = cache_path.read_text(encoding="utf-8")
        lines = [line.strip() for line in fasta_text.splitlines() if line.strip()]
        if not lines or not lines[0].startswith(">"):
            return _fetch_refseq_target_from_ortholog_table(accession, batch_dir=batch_dir)
        header = lines[0][1:]
        sequence = "".join(lines[1:])
        gene_symbol = accession
        for token in header.replace(",", " ").split():
            cleaned = token.strip("[]();")
            if cleaned.isupper() and 2 <= len(cleaned) <= 10 and "_" not in cleaned and not cleaned.startswith(("XP", "NP")):
                gene_symbol = cleaned
                break
        return {
            "accession": accession,
            "entry_name": accession,
            "gene_symbol": gene_symbol,
            "sequence": sequence,
        }
    except Exception:
        pass
    return _fetch_refseq_target_from_ortholog_table(accession, batch_dir=batch_dir)


def _cached_response_path(namespace: str, url: str, suffix: str, *, batch_dir: Path) -> Path:
    for cache_dir in (
        batch_dir / ".agdesign2" / "cache",
        batch_dir.parent / ".agdesign2" / "cache",
        Path.cwd() / ".agdesign2" / "cache",
    ):
        path = FileCache(cache_dir.resolve()).path_for(namespace, url, suffix)
        if path.exists():
            return path
    return FileCache((batch_dir / ".agdesign2" / "cache").resolve()).path_for(namespace, url, suffix)


_ORTHOLOG_REFSEQ_SPECIES_COLUMNS = (
    ("human", "human_gene_symbol", "human_refseq_accession", "human_refseq_sequence"),
    ("mouse", "mouse_gene_symbol", "mouse_refseq_accession", "mouse_refseq_sequence"),
    (
        "macaca_fascicularis",
        "macaca_fascicularis_gene_symbol",
        "macaca_fascicularis_refseq_accession",
        "macaca_fascicularis_refseq_sequence",
    ),
)


def _ortholog_refseq_index(table_path: Path) -> dict[str, dict[str, Any] | None]:
    """Map each refseq accession to its target dict (or None) from the ortholog
    table, cached per (path, mtime). Records the FIRST (row, species)-order
    occurrence of each accession — including the empty-sequence case, which the
    original per-call scan short-circuited to None — so lookups are byte-for-byte
    equivalent to the former linear scan.
    """
    try:
        mtime = table_path.stat().st_mtime_ns
    except OSError:
        return {}
    cache_key = (str(table_path), mtime)
    cached = _ORTHOLOG_REFSEQ_INDEX_CACHE.get(cache_key)
    if cached is not None:
        return cached
    index: dict[str, dict[str, Any] | None] = {}
    try:
        with table_path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                for _species, gene_column, accession_column, sequence_column in _ORTHOLOG_REFSEQ_SPECIES_COLUMNS:
                    acc = str(row.get(accession_column) or "")
                    if not acc or acc in index:
                        continue
                    sequence = str(row.get(sequence_column) or "")
                    if not sequence:
                        index[acc] = None
                    else:
                        index[acc] = {
                            "accession": acc,
                            "entry_name": acc,
                            "gene_symbol": row.get(gene_column) or acc,
                            "sequence": sequence,
                        }
    except Exception:
        return {}
    _ORTHOLOG_REFSEQ_INDEX_CACHE[cache_key] = index
    return index


def _fetch_refseq_target_from_ortholog_table(accession: str, *, batch_dir: Path) -> dict[str, Any] | None:
    table_path = batch_dir / "ortholog_reference_table.tsv"
    if not table_path.exists():
        return None
    result = _ortholog_refseq_index(table_path).get(accession)
    return dict(result) if result is not None else None


def _resolve_asset_for_construct(
    path_text: str | None,
    fallback: Path,
    *,
    base_dir: Path,
) -> Path | None:
    if path_text:
        resolved = _resolve_existing_path(path_text, base_dir)
        if resolved is not None:
            return resolved
    if fallback.exists():
        return fallback.resolve()
    return None


def _relative_link(target: Path, from_dir: Path) -> str:
    return os.path.relpath(Path(target).resolve(), Path(from_dir).resolve())


def _format_region(region: dict[str, Any]) -> str:
    start = region.get("start")
    end = region.get("end")
    if start is None or end is None:
        return ""
    return f"{start}-{end}"


def _fmt_pct(value: Any) -> str:
    if value is None:
        return ""
    try:
        return f"{float(value):.2f}%"
    except Exception:
        return str(value)


def _fmt_surface_identity(item: dict[str, Any]) -> str:
    value = item.get("extracellular_surface_identity")
    if value is None:
        value = item.get("extracellular_surface_identity_to_human")
    formatted = _fmt_pct(value)
    if not formatted:
        return ""
    aligned = item.get("extracellular_surface_aligned_positions")
    if aligned is None:
        return formatted
    return f"{formatted} ({aligned} aa)"


def _fmt_fraction_pct(value: Any) -> str:
    if value is None:
        return ""
    try:
        return f"{float(value) * 100:.2f}%"
    except Exception:
        return str(value)


def _fmt_number(value: Any) -> str:
    if value is None:
        return ""
    try:
        return f"{float(value):.2f}"
    except Exception:
        return str(value)


def _render_identity_matrix_cell(value: Any) -> str:
    if value is None:
        return '<td class="matrix-cell matrix-cell-empty"></td>'
    try:
        numeric = max(0.0, min(100.0, float(value)))
    except Exception:
        return f'<td class="matrix-cell">{escape(str(value))}</td>'
    ratio = numeric / 100.0
    lightness = 97.0 - (ratio * 52.0)
    text_color = "#f7fbff" if ratio >= 0.72 else "#17324a"
    border_color = "#0e8c88" if ratio >= 0.72 else "#d5e4ef"
    style = (
        f"background:hsl(184 58% {lightness:.1f}%);"
        f"color:{text_color};"
        f"font-weight:{700 if ratio >= 0.55 else 600};"
        f"box-shadow:inset 0 0 0 1px {border_color};"
    )
    label = _fmt_pct(numeric)
    return (
        f'<td class="matrix-cell" data-identity="{numeric:.2f}" '
        f'style="{style}" title="{escape(label)}">{escape(label)}</td>'
    )


def _truncate(text: str, limit: int) -> str:
    if len(text) <= limit:
        return text
    return text[: limit - 3] + "..."


def _safe_filename(text: str) -> str:
    return "".join(ch if ch.isalnum() or ch in {"-", "_"} else "_" for ch in text)


def _construct_workbench_css() -> str:
    """Return the shared styles for the integrated construct workbench."""

    return """
#interactive-construct-builder.construct-workbench {
  padding: 0;
  overflow: hidden;
  background: #fff;
  scroll-margin-top: 12px;
}
.construct-workbench-header {
  margin: 0;
  padding: 15px 20px;
  border-bottom: 1px solid var(--line);
  background: #f8fbfc;
}
.construct-workbench-header h2 {
  margin: 0;
  font-family: Inter, ui-sans-serif, system-ui, sans-serif;
  font-size: 1.55rem;
  font-weight: 700;
  letter-spacing: -0.025em;
}
.construct-workbench-header .construct-choice-help {
  margin: 6px 0 0;
}
.construct-workbench-shell {
  display: grid;
  grid-template-columns: 190px minmax(560px, 1fr) minmax(312px, 330px);
  gap: 0;
  align-items: stretch;
  min-width: 0;
}
.construct-workbench-picker {
  display: flex;
  flex-direction: column;
  min-width: 0;
  min-height: 0;
  padding: 12px;
  border-right: 1px solid var(--line);
}
.construct-workbench-picker-heading {
  display: grid;
  gap: 3px;
  margin-bottom: 10px;
}
.construct-workbench-picker-heading > strong,
.construct-workbench-toolbar-label,
.construct-workbench-current .selection-live-header h3,
.construct-workbench-evidence .plot-card h3 {
  margin: 0;
  color: var(--ink);
  font-family: Inter, ui-sans-serif, system-ui, sans-serif;
  font-size: 0.9rem;
  font-weight: 700;
  line-height: 1.2;
  letter-spacing: -0.015em;
}
.construct-workbench-picker-heading span {
  color: var(--muted);
  font-size: 0.76rem;
  line-height: 1.35;
}
.construct-workbench-picker > input {
  width: 100%;
  min-height: 36px;
  margin-bottom: 8px;
  padding: 8px 10px;
  border: 1px solid var(--line);
  border-radius: 9px;
  font: inherit;
}
.construct-workbench-list {
  display: grid;
  flex: 1 1 0;
  min-height: 0;
  gap: 2px;
  align-content: start;
  grid-auto-rows: max-content;
  overflow: auto;
  overscroll-behavior: contain;
}
.construct-workbench-option[hidden] {
  display: none;
}
.construct-workbench-option {
  display: grid;
  grid-template-columns: minmax(0, 1fr) auto;
  gap: 2px 8px;
  width: 100%;
  padding: 9px 8px;
  border: 0;
  border-bottom: 1px solid var(--line);
  background: transparent;
  color: var(--ink);
  text-align: left;
  cursor: pointer;
}
.construct-workbench-option:hover {
  background: #f3f8f8;
}
.construct-workbench-option[aria-pressed="true"] {
  background: #e5f4f2;
  box-shadow: inset 3px 0 #147d7e;
}
.construct-workbench-option strong {
  grid-column: 1 / -1;
  overflow-wrap: anywhere;
  font: 750 0.8rem/1.25 ui-monospace, SFMono-Regular, Menlo, monospace;
}
.construct-workbench-option small {
  grid-column: 1;
  color: var(--muted);
  font-size: 0.7rem;
}
.construct-workbench-boundary {
  color: #147d7e;
  font: 750 0.74rem/1.3 ui-monospace, SFMono-Regular, Menlo, monospace;
}
.construct-workbench-view-all {
  width: 100%;
  margin-top: 10px;
  padding: 9px;
  border: 1px solid var(--line);
  border-radius: 9px;
  background: #fff;
  color: #147d7e;
  font-weight: 750;
  cursor: pointer;
}
.construct-workbench-center {
  display: grid;
  min-width: 0;
  gap: 0;
  align-content: start;
  background: #fff;
}
.construct-workbench-center .viewer-panel {
  display: grid;
  min-width: 0;
  gap: 0;
}
.construct-workbench-center .viewer-toolbar {
  display: flex;
  flex-wrap: nowrap;
  align-items: center;
  min-height: 44px;
  gap: 8px;
  padding: 7px 9px;
  overflow-x: auto;
  border: 0;
  border-bottom: 1px solid rgba(16, 42, 67, 0.14);
  border-radius: 0;
  background: #f3f8f9;
  scrollbar-width: thin;
}
.construct-workbench-toolbar-label {
  flex: 0 0 auto;
  white-space: nowrap;
}
.construct-workbench-toolbar-controls {
  display: flex;
  flex: 1 0 auto;
  min-width: max-content;
  flex-wrap: nowrap;
  align-items: center;
  gap: 5px;
}
.construct-workbench-toolbar-group {
  display: inline-flex;
  flex: 0 0 auto;
  gap: 0;
}
.construct-workbench-toolbar-controls .mini-button {
  min-height: 30px;
  margin: 0 0 -1px -1px;
  padding: 5px 6px;
  border-radius: 0;
  font-family: Inter, ui-sans-serif, system-ui, sans-serif;
  font-size: 0.72rem;
  font-weight: 600;
  white-space: nowrap;
}
.construct-workbench-toolbar-controls .mini-button:first-child {
  border-radius: 7px 0 0 7px;
}
.construct-workbench-toolbar-controls .mini-button:last-child {
  border-radius: 0 7px 7px 0;
}
.construct-workbench-viewer .structure-viewer {
  position: relative;
  min-height: 410px;
  border: 0;
  border-radius: 0;
}
.construct-workbench-structure-legend,
.construct-workbench-structure-status {
  position: absolute;
  bottom: 8px;
  z-index: 3;
  padding: 4px 7px;
  border: 1px solid rgba(20, 125, 126, 0.22);
  border-radius: 6px;
  background: rgba(255, 255, 255, 0.9);
  color: #536f85;
  font-size: 0.68rem;
  line-height: 1.2;
  pointer-events: none;
}
.construct-workbench-structure-legend {
  right: 8px;
}
.construct-workbench-structure-status {
  left: 8px;
  max-width: calc(100% - 225px);
  overflow: hidden;
  font-weight: 700;
  text-overflow: ellipsis;
  white-space: nowrap;
}
.construct-workbench-evidence {
  display: grid;
  grid-template-columns: minmax(0, 1.2fr) minmax(260px, 0.8fr);
  gap: 0;
  align-items: start;
  border-top: 1px solid rgba(16, 42, 67, 0.14);
}
.construct-workbench-evidence .plot-card {
  display: grid;
  grid-template-rows: auto auto auto;
  min-width: 0;
  padding: 0;
  overflow: hidden;
  border: 0;
  border-radius: 0;
  background: #fff;
  box-shadow: none;
}
.construct-workbench-evidence .plot-card + .plot-card {
  border-left: 1px solid rgba(16, 42, 67, 0.14);
}
.construct-workbench-evidence .plot-card-header {
  display: flex;
  justify-content: space-between;
  align-items: center;
  flex-wrap: nowrap;
  min-height: 44px;
  gap: 6px;
  margin: 0;
  padding: 7px 9px;
  border-bottom: 1px solid rgba(16, 42, 67, 0.1);
  background: #f3f8f9;
}
.construct-workbench-evidence .plot-card h3 {
  white-space: nowrap;
}
.construct-workbench-evidence .plot-controls {
  display: flex;
  flex: 0 0 auto;
  flex-wrap: nowrap;
  width: max-content;
  gap: 2px;
  padding: 3px;
  border-radius: 9px;
  background: #e5eff1;
}
.construct-workbench-evidence .plot-controls .mini-button {
  min-height: 26px;
  margin: 0;
  padding: 4px 6px;
  border: 0;
  border-radius: 6px;
  background: transparent;
  font-family: Inter, ui-sans-serif, system-ui, sans-serif;
  font-size: 0.68rem;
  font-weight: 600;
  white-space: nowrap;
}
.construct-workbench-evidence .plot-controls .mini-button:hover {
  background: #fff;
  box-shadow: 0 1px 4px rgba(16, 42, 67, 0.12);
}
.construct-workbench-evidence .pae-canvas {
  display: block;
  width: 100%;
  max-width: none;
  margin: 0;
  border: 0;
  border-radius: 0;
  box-shadow: none;
}
.construct-workbench-evidence .viewer-info {
  display: flex;
  align-items: center;
  min-height: 42px;
  padding: 8px 10px;
  border-top: 1px solid rgba(16, 42, 67, 0.08);
  background: #f8fafb;
  color: #526979;
  font-family: Inter, ui-sans-serif, system-ui, sans-serif;
  font-size: 0.7rem;
  font-variant-numeric: tabular-nums;
  font-weight: 500;
  line-height: 1.4;
  overflow-wrap: anywhere;
}
.construct-workbench-evidence .topology-legend {
  display: none;
}
.construct-workbench-sequence-details {
  padding: 10px 12px;
  border-top: 1px solid var(--line);
}
.construct-workbench-sequence-details > summary {
  color: #147d7e;
  font-weight: 800;
  cursor: pointer;
}
.construct-workbench-sequence-details .sequence-panel {
  margin-top: 10px;
}
.construct-workbench-current {
  height: 100%;
  max-height: none;
  margin: 0;
  padding: 12px;
  overflow: auto;
  border: 0;
  border-left: 1px solid var(--line);
  border-radius: 0;
  box-shadow: none;
}
.construct-workbench-current .selection-live-header {
  align-items: center;
  flex-wrap: wrap;
  margin-bottom: 6px;
}
.construct-workbench-current .selection-live-header .button-row {
  flex-wrap: nowrap;
  gap: 0;
}
.construct-workbench-current .selection-live-header .mini-button {
  min-height: 32px;
  padding: 6px 8px;
  border-color: #147d7e;
  border-radius: 0;
  background: #147d7e;
  color: #fff;
  font-size: 0.76rem;
  font-weight: 800;
  white-space: nowrap;
}
.construct-workbench-current .selection-live-header h3 {
  flex: 0 0 auto;
  white-space: nowrap;
}
.construct-workbench-current .selection-live-header .mini-button:first-child {
  border-radius: 7px 0 0 7px;
}
.construct-workbench-current .selection-live-header .mini-button:last-child {
  border-left-color: rgba(255, 255, 255, 0.35);
  border-radius: 0 7px 7px 0;
}
.construct-workbench-boundary-editor {
  display: grid;
  grid-template-columns: minmax(68px, 0.7fr) minmax(68px, 0.7fr) minmax(124px, 1.2fr);
  gap: 6px;
  margin: 6px 0 8px;
  padding: 7px 0;
  border-top: 1px solid var(--line);
  border-bottom: 1px solid var(--line);
}
.construct-workbench-boundary-editor label {
  display: grid;
  gap: 4px;
  color: var(--muted);
  font-size: 0.76rem;
  font-weight: 750;
}
.construct-workbench-boundary-editor input {
  width: 100%;
  min-height: 32px;
  padding: 5px 7px;
  border: 1px solid var(--line);
  border-radius: 7px;
  font: 750 0.86rem/1 ui-monospace, SFMono-Regular, Menlo, monospace;
}
.construct-workbench-boundary-editor button {
  align-self: end;
  min-height: 32px;
  padding: 5px 8px;
  border: 0;
  border-radius: 7px;
  background: #147d7e;
  color: #fff;
  font-size: 0.76rem;
  font-weight: 800;
  cursor: pointer;
}
.construct-workbench-boundary-message {
  grid-column: 1 / -1;
  min-height: 0;
  color: var(--muted);
  font-size: 0.75rem;
}
.construct-workbench-boundary-message:empty {
  display: none;
}
.construct-workbench-species {
  display: grid;
  gap: 0;
  margin: 10px 0;
  border-top: 1px solid var(--line);
}
.construct-workbench-species-row {
  display: grid;
  gap: 6px;
  padding: 10px 0;
  border-bottom: 1px solid var(--line);
}
.construct-workbench-species-row > div {
  display: flex;
  justify-content: space-between;
  gap: 10px;
  color: var(--muted);
  font-size: 0.78rem;
}
.construct-workbench-species-row strong {
  color: var(--ink);
}
.construct-workbench-species-row code {
  max-height: 62px;
  overflow: auto;
  white-space: pre-wrap;
  word-break: break-all;
  font-size: 0.78rem;
}
.construct-workbench-species-empty {
  color: var(--muted);
  font-size: 0.86rem;
  line-height: 1.45;
}
.construct-workbench-original-table {
  display: none;
}
.construct-workbench-raw-exports {
  margin-top: 12px;
  color: #147d7e;
  font-weight: 750;
}
.construct-workbench-raw-exports textarea {
  min-height: 90px;
  font-weight: 500;
}
.construct-workbench-current .selection-warnings {
  margin: 8px 0;
}
.construct-workbench-current .selection-warning-block {
  padding: 8px 10px;
}
@media (max-width: 1160px) {
  .construct-workbench-shell {
    grid-template-columns: 190px minmax(470px, 1fr) minmax(300px, 340px);
  }
}
@media (max-width: 980px) {
  .construct-workbench-shell {
    grid-template-columns: 1fr;
  }
  .construct-workbench-picker {
    max-height: 420px;
    border-right: 0;
    border-bottom: 1px solid var(--line);
  }
  .construct-workbench-list {
    flex-basis: auto;
    max-height: 260px;
  }
  .construct-workbench-current {
    border-top: 1px solid var(--line);
    border-left: 0;
  }
}
@media (max-width: 720px) {
  .construct-workbench-evidence {
    grid-template-columns: 1fr;
  }
  .construct-workbench-evidence .plot-card + .plot-card {
    border-top: 1px solid rgba(16, 42, 67, 0.14);
    border-left: 0;
  }
  .construct-workbench-current .selection-live-header {
    align-items: flex-start;
    flex-direction: column;
  }
  .construct-workbench-boundary-editor {
    grid-template-columns: 1fr 1fr;
  }
  .construct-workbench-boundary-editor button {
    grid-column: 1 / -1;
  }
}
"""


def _portal_css(*, theme: str = "human") -> str:
    theme_overrides = ""
    if theme == "mouse":
        theme_overrides = """
body {
  background:
    radial-gradient(circle at top right, rgba(203, 58, 76, 0.07), transparent 34%),
    radial-gradient(circle at 15% 10%, rgba(185, 72, 90, 0.055), transparent 32%),
    linear-gradient(180deg, #fff7f8 0%, #fffdfd 34%, #faf6f7 100%);
}
.site-hero {
  background:
    radial-gradient(circle at top right, rgba(255, 189, 194, 0.18), transparent 30%),
    linear-gradient(135deg, #6e2434 0%, #b64a5d 54%, #d17682 100%);
  box-shadow: 0 28px 70px rgba(166, 62, 80, 0.14);
}
:root {
  --bg: #fff7f8;
  --bg-deep: #6e2434;
  --bg-panel: rgba(121, 42, 58, 0.58);
  --card-soft: #fff9fa;
  --ink: #33262b;
  --muted: #785f66;
  --line: #efdce1;
  --accent: #bd4d61;
  --accent-2: #96364c;
  --accent-soft: #f9e8ec;
  --shadow: 0 22px 60px rgba(174, 75, 94, 0.09);
}
"""
    return """
:root {
  --bg: #f5f8fc;
  --bg-deep: #071a2f;
  --bg-panel: rgba(8, 27, 48, 0.74);
  --card: #ffffff;
  --card-soft: #f8fbff;
  --ink: #102338;
  --muted: #587089;
  --line: #d8e3f0;
  --accent: #0e8c88;
  --accent-2: #1f4d8c;
  --accent-soft: #dff7f5;
  --hero-line: rgba(255, 255, 255, 0.12);
  --error: #b7424b;
  --ok: #1f8f66;
  --warn: #b97416;
  --shadow: 0 22px 60px rgba(10, 35, 66, 0.12);
}
* { box-sizing: border-box; }
body {
  margin: 0;
  font-family: "Avenir Next", "Segoe UI", "Helvetica Neue", Arial, sans-serif;
  color: var(--ink);
  background:
    radial-gradient(circle at top right, rgba(19, 140, 136, 0.14), transparent 32%),
    radial-gradient(circle at 15% 10%, rgba(31, 77, 140, 0.14), transparent 30%),
    linear-gradient(180deg, #eef4fb 0%, #f8fbff 28%, #f4f7fb 100%);
}
.page {
  width: min(1480px, calc(100vw - 40px));
  margin: 24px auto 72px;
}
.site-hero {
  position: relative;
  overflow: hidden;
  margin-bottom: 22px;
  padding: 24px 28px 28px;
  border-radius: 28px;
  background:
    radial-gradient(circle at top right, rgba(14, 140, 136, 0.25), transparent 28%),
    linear-gradient(135deg, #071a2f 0%, #102f55 45%, #0a5a66 100%);
  color: #f7fbff;
  box-shadow: 0 28px 70px rgba(8, 27, 48, 0.22);
}
.site-hero-index {
  padding-top: 22px;
  padding-bottom: 24px;
}
.site-hero::after {
  content: "";
  position: absolute;
  inset: auto -10% -45% 35%;
  height: 320px;
  background: radial-gradient(circle, rgba(255,255,255,0.12), transparent 58%);
  pointer-events: none;
}
.brand-bar {
  position: relative;
  z-index: 1;
  display: flex;
  justify-content: space-between;
  align-items: flex-start;
  gap: 24px;
  flex-wrap: nowrap;
  margin-bottom: 18px;
}
.site-hero-index .brand-bar {
  margin-bottom: 14px;
}
.brand-lockup {
  display: inline-flex;
  align-items: center;
  gap: 16px;
  color: inherit;
}
.brand-home-link,
.brand-kicker,
.brand-product {
  color: #fff !important;
  text-decoration: none !important;
}
.site-hero .brand-home-link,
.site-hero .brand-kicker,
.site-hero .brand-product,
.site-hero .brand-home-link:visited,
.site-hero .brand-kicker:visited,
.site-hero .brand-product:visited {
  color: #fff !important;
  text-decoration: none !important;
}
.brand-home-link {
  display: inline-flex;
  align-items: center;
}
.brand-home-link:hover,
.brand-kicker:hover,
.brand-product:hover {
  opacity: 0.86;
}
.brand-logo {
  display: block;
  height: 48px;
  width: auto;
}
.brand-copy {
  display: grid;
  gap: 2px;
}
.brand-actions {
  position: relative;
  z-index: 2;
  display: flex;
  flex-direction: column;
  align-items: center;
  gap: 10px;
  justify-content: flex-end;
  margin-left: auto;
  max-width: min(820px, 68vw);
}
.brand-secondary-actions {
  display: inline-flex;
  align-items: center;
  justify-content: flex-end;
  gap: 8px;
  width: 100%;
}
.site-nav {
  display: inline-flex;
  align-items: center;
  justify-content: flex-end;
  flex-wrap: nowrap;
  gap: 6px;
  padding: 5px;
  border-radius: 999px;
  border: 1px solid rgba(255,255,255,0.14);
  background: rgba(255,255,255,0.08);
  backdrop-filter: blur(14px);
}
.site-nav a {
  display: inline-flex;
  align-items: center;
  min-height: 30px;
  padding: 7px 10px;
  border-radius: 999px;
  color: rgba(246, 251, 255, 0.82);
  text-decoration: none;
  font-weight: 700;
  font-size: 0.82rem;
}
.site-nav a:hover,
.site-nav a.active {
  color: #071a2f;
  background: rgba(255,255,255,0.9);
}
.species-switcher {
  display: inline-flex;
  align-items: center;
  gap: 4px;
  padding: 4px;
  border-radius: 999px;
  border: 1px solid rgba(255,255,255,0.16);
  background: rgba(255,255,255,0.08);
  backdrop-filter: blur(14px);
  white-space: nowrap;
}
.species-current,
.species-link {
  display: inline-flex;
  align-items: center;
  min-height: 28px;
  padding: 5px 10px;
  border-radius: 999px;
  font-size: 0.8rem;
  font-weight: 800;
  text-decoration: none;
}
.species-current {
  color: #071a2f;
  background: rgba(255,255,255,0.92);
}
.species-link {
  color: rgba(246, 251, 255, 0.84);
}
.species-link:hover {
  color: #071a2f;
  background: rgba(255,255,255,0.9);
}
.brand-kicker {
  font-size: 0.82rem;
  letter-spacing: 0.12em;
  text-transform: uppercase;
  color: rgba(255,255,255,0.7);
}
.brand-product {
  font-size: 1.18rem;
  font-weight: 700;
  letter-spacing: 0.01em;
}
.hero-chip {
  display: inline-flex;
  align-items: center;
  padding: 8px 12px;
  border-radius: 999px;
  border: 1px solid rgba(255,255,255,0.16);
  background: rgba(255,255,255,0.08);
  color: #eef6ff;
  font-size: 0.84rem;
  white-space: nowrap;
  backdrop-filter: blur(16px);
}
.hero-chip-link {
  text-decoration: none;
}
.site-hero-doc .hero-copy h1 {
  max-width: 26ch;
}
.site-hero-doc .hero-grid {
  grid-template-columns: minmax(0, 1fr);
  gap: 0;
}
.hero-grid {
  position: relative;
  z-index: 1;
  display: grid;
  grid-template-columns: minmax(0, 1.45fr) minmax(300px, 0.55fr);
  gap: 24px;
  align-items: center;
}
.hero-copy h1 {
  margin: 0 0 12px;
  font-size: clamp(2rem, 3vw, 3.65rem);
  line-height: 1.02;
  letter-spacing: -0.03em;
  max-width: 16ch;
}
.site-hero-index .hero-copy h1 {
  max-width: 20ch;
}
.hero-copy p {
  margin: 0;
}
.eyebrow {
  margin: 0 0 10px;
  color: #91d7cf;
  font-size: 0.85rem;
  letter-spacing: 0.16em;
  text-transform: uppercase;
}
.hero-text {
  max-width: 62ch;
  color: rgba(241, 247, 255, 0.82);
  font-size: 1.02rem;
  line-height: 1.6;
}
.hero-panel {
  padding: 18px;
  border-radius: 22px;
  border: 1px solid var(--hero-line);
  background: rgba(255,255,255,0.08);
  backdrop-filter: blur(14px);
}
.site-hero-index .hero-panel {
  padding: 14px;
}
.hero-actions-panel {
  display: flex;
  align-items: center;
  justify-content: flex-end;
}
.stats {
  display: grid;
  grid-template-columns: repeat(2, minmax(120px, 1fr));
  gap: 14px;
}
.stat, .card, .construct-card {
  background: var(--card);
  border: 1px solid var(--line);
  border-radius: 22px;
  box-shadow: var(--shadow);
}
.partner-warning-card {
  margin: 0 0 18px;
  padding: 20px 22px;
  border-radius: 22px;
  border: 1px solid rgba(20, 120, 130, 0.28);
  background:
    linear-gradient(135deg, rgba(231, 249, 247, 0.96), rgba(240, 249, 255, 0.92)),
    #eefaf8;
  box-shadow: 0 18px 42px rgba(12, 84, 94, 0.1);
}
.partner-warning-card h2 {
  margin: 0 0 8px;
  color: #0b4f5a;
}
.partner-warning-card p {
  max-width: 100ch;
  margin: 0 0 10px;
  color: #153f48;
  line-height: 1.55;
}
.partner-warning-card p:last-child {
  margin-bottom: 0;
}
.warning-kicker {
  color: #0f7a7f !important;
  font-size: 0.78rem;
  font-weight: 800;
  letter-spacing: 0.14em;
  text-transform: uppercase;
}
.stat {
  padding: 14px 16px;
  text-align: center;
  background: rgba(255,255,255,0.1);
  border-color: rgba(255,255,255,0.1);
  box-shadow: none;
  color: #f7fbff;
}
.stat-primary {
  grid-column: 1 / -1;
  text-align: left;
  background: rgba(255,255,255,0.16);
}
.stat span {
  display: block;
  font-size: 1.65rem;
  font-weight: 700;
  color: #ffffff;
}
.stat label {
  color: rgba(241,247,255,0.72);
  font-size: 0.88rem;
}
.status-strip {
  display: flex;
  flex-wrap: nowrap;
  gap: 6px;
  margin-top: 10px;
  white-space: nowrap;
}
.status-strip span {
  display: inline-flex;
  align-items: center;
  justify-content: center;
  flex: 1 1 0;
  min-width: 0;
  min-height: 28px;
  padding: 6px 8px;
  border-radius: 999px;
  background: rgba(255,255,255,0.1);
  border: 1px solid rgba(255,255,255,0.12);
  color: rgba(241,247,255,0.82);
  font-size: clamp(0.74rem, 0.82vw, 0.85rem);
  font-weight: 700;
  overflow: hidden;
  text-overflow: ellipsis;
}
.controls {
  display: flex;
  gap: 12px;
  margin-bottom: 12px;
  flex-wrap: wrap;
}
.controls input, .controls select, .controls button {
  border: 1px solid var(--line);
  background: rgba(255,255,255,0.82);
  border-radius: 999px;
  padding: 12px 16px;
  font: inherit;
  color: var(--ink);
  box-shadow: 0 8px 18px rgba(16, 35, 56, 0.05);
}
.controls input { flex: 1; }
.controls button { cursor: pointer; }
.disease-cell {
  display: grid;
  gap: 3px;
  min-width: 180px;
}
.muted-small {
  color: var(--muted);
  font-size: 0.82rem;
}
.sort-summary {
  display: inline-flex;
  align-items: center;
  padding: 12px 16px;
  border: 1px solid var(--line);
  border-radius: 999px;
  background: rgba(255,255,255,0.82);
  color: var(--muted);
  font-weight: 600;
  box-shadow: 0 8px 18px rgba(16, 35, 56, 0.05);
}
.pagination-bar {
  display: flex;
  flex-wrap: wrap;
  gap: 10px;
  align-items: center;
  justify-content: flex-end;
  margin: -4px 0 16px;
}
.page-size-control {
  display: inline-flex;
  align-items: center;
  gap: 8px;
  margin-right: auto;
  padding: 6px 8px 6px 12px;
  border: 1px solid var(--line);
  border-radius: 999px;
  background: rgba(255,255,255,0.82);
  color: var(--muted);
  font-weight: 700;
  box-shadow: 0 8px 18px rgba(16, 35, 56, 0.05);
}
.page-size-control select {
  min-height: 32px;
  border: 1px solid var(--line);
  border-radius: 999px;
  background: #ffffff;
  color: var(--ink);
  font: inherit;
  font-weight: 700;
  padding: 4px 26px 4px 10px;
}
.pagination-summary {
  min-width: 220px;
  color: var(--muted);
  font-weight: 700;
  text-align: center;
}
.card {
  padding: 22px;
  margin-bottom: 18px;
  overflow: auto;
  background: linear-gradient(180deg, #ffffff 0%, #fbfdff 100%);
}
.table-card {
  padding: 8px 12px 12px;
}
.stack {
  display: grid;
  grid-template-columns: 1fr;
  gap: 18px;
}
.grid.two {
  display: grid;
  grid-template-columns: repeat(2, minmax(0, 1fr));
  gap: 18px;
}
.doc-layout {
  display: grid;
  grid-template-columns: minmax(0, 1fr);
  gap: 18px;
}
.doc-grid,
.download-grid {
  display: grid;
  grid-template-columns: repeat(2, minmax(0, 1fr));
  gap: 18px;
}
.calculator-layout {
  display: grid;
  grid-template-columns: minmax(0, 1.18fr) minmax(320px, 0.82fr);
  gap: 18px;
  align-items: start;
}
.calculator-hero-panel {
  display: grid;
  gap: 10px;
  align-self: stretch;
  align-content: center;
  color: rgba(241, 247, 255, 0.84);
  font-weight: 700;
}
.calculator-card {
  overflow: visible;
}
.calculator-card-header {
  display: flex;
  justify-content: space-between;
  align-items: center;
  gap: 14px;
  margin-bottom: 12px;
}
.calc-status {
  display: inline-flex;
  align-items: center;
  min-height: 32px;
  padding: 6px 11px;
  border-radius: 999px;
  font-size: 0.82rem;
  font-weight: 800;
  white-space: nowrap;
}
.calc-status.waiting {
  background: #fff4d7;
  color: #8a5b00;
}
.calc-status.ready {
  background: #e8f8f1;
  color: var(--ok);
}
.calculator-form {
  display: grid;
  border: 1px solid var(--line);
  border-radius: 16px;
  overflow: hidden;
  background: #fff;
}
.calculator-row {
  display: grid;
  grid-template-columns: minmax(190px, 0.95fr) minmax(160px, 1fr) minmax(118px, 0.5fr);
  gap: 12px;
  align-items: center;
  min-height: 78px;
  padding: 14px;
  border-bottom: 1px solid var(--line);
}
.calculator-row:last-child {
  border-bottom: 0;
}
.calculator-row label {
  display: grid;
  gap: 4px;
  color: var(--ink);
  font-weight: 800;
  line-height: 1.25;
}
.calculator-row label span {
  color: var(--muted);
  font-size: 0.78rem;
  font-weight: 600;
  line-height: 1.35;
}
.calculator-row label .required-tag {
  display: inline-flex;
  width: max-content;
  padding: 2px 7px;
  border-radius: 999px;
  background: var(--accent-soft);
  color: var(--accent);
  font-size: 0.7rem;
  font-weight: 900;
  text-transform: uppercase;
  letter-spacing: 0.06em;
}
.calculator-row input,
.calculator-row select {
  width: 100%;
  min-height: 42px;
  border: 1px solid var(--line);
  border-radius: 12px;
  background: #fff;
  color: var(--ink);
  font: inherit;
  font-size: 0.95rem;
}
.calculator-row input {
  padding: 9px 11px;
  font-variant-numeric: tabular-nums;
}
.calculator-row select {
  padding: 9px 10px;
  cursor: pointer;
}
.calculator-row select:disabled {
  color: var(--muted);
  cursor: default;
  background: #f5f8fc;
}
.calculator-row input.calculated {
  background: #f5fbff;
  color: #24435f;
}
.calculator-row input.user-set {
  border-color: var(--accent);
  box-shadow: inset 0 0 0 1px rgba(14, 140, 136, 0.22);
}
.calculator-actions {
  display: flex;
  justify-content: flex-end;
  gap: 10px;
  flex-wrap: wrap;
  margin-top: 14px;
}
.calc-message {
  min-height: 92px;
  margin-bottom: 16px;
  padding: 13px 14px;
  border: 1px solid var(--line);
  border-radius: 16px;
  background: #f6faff;
  color: var(--muted);
  line-height: 1.5;
}
.calc-message strong {
  color: var(--ink);
}
.calc-message.error {
  border-color: #e1a4aa;
  background: #fff2f3;
  color: #8d313a;
}
.calc-message.ready {
  border-color: #a5ddc8;
  background: #f1fbf6;
  color: #1f6d50;
}
.calculator-facts {
  display: grid;
  grid-template-columns: 1fr auto;
  gap: 0 14px;
  margin: 0 0 16px;
}
.calculator-facts dt,
.calculator-facts dd {
  margin: 0;
  padding: 10px 0;
  border-bottom: 1px solid var(--line);
}
.calculator-facts dt {
  color: var(--muted);
  font-weight: 700;
}
.calculator-facts dd {
  text-align: right;
  font-weight: 800;
  font-variant-numeric: tabular-nums;
}
.calculator-notes {
  font-size: 0.92rem;
}
.doc-card {
  overflow: visible;
}
.doc-card .inline-link {
  max-width: 100%;
  overflow-wrap: anywhere;
}
.doc-card p {
  margin: 0 0 12px;
  color: var(--muted);
  line-height: 1.65;
}
.doc-card p:last-child {
  margin-bottom: 0;
}
.doc-list {
  margin: 0;
  padding-left: 22px;
  color: var(--muted);
  line-height: 1.7;
}
.doc-list li + li {
  margin-top: 8px;
}
#choose-constructs,
#after-export,
#construct-classes,
#construct-feature-review {
  scroll-margin-top: 16px;
}
.construct-choice-workflow {
  display: grid;
  grid-template-columns: repeat(4, minmax(0, 1fr));
  gap: 28px;
  list-style: none;
  margin: 20px 0;
  padding: 0;
}
.construct-choice-workflow li {
  position: relative;
  padding: 14px;
  border: 1px solid #b9dbe4;
  border-radius: 12px;
  background: #f0f8fa;
}
.construct-choice-workflow strong {
  display: block;
  margin-bottom: 6px;
  color: #0d1f3c;
}
.construct-choice-workflow span {
  color: var(--muted);
  font-size: 0.9rem;
  line-height: 1.5;
}
.construct-choice-workflow li:not(:last-child)::after {
  content: "";
  position: absolute;
  top: calc(50% - 5px);
  right: -19px;
  width: 8px;
  height: 8px;
  border-top: 2px solid #007fa3;
  border-right: 2px solid #007fa3;
  transform: rotate(45deg);
}
.construct-choice-table {
  table-layout: fixed;
  overflow-wrap: anywhere;
}
.doc-card .construct-choice-more {
  margin-top: 14px;
}
@media (max-width: 760px) {
  .construct-choice-workflow {
    grid-template-columns: 1fr;
  }
  .construct-choice-workflow li:not(:last-child)::after {
    top: auto;
    right: calc(50% - 5px);
    bottom: -19px;
    transform: rotate(135deg);
  }
}
.code-block {
  overflow: auto;
  margin: 12px 0;
  padding: 14px;
  border-radius: 14px;
  background: #081b30;
  color: #eef7ff;
  font-size: 0.9rem;
  line-height: 1.5;
}
table {
  width: 100%;
  border-collapse: collapse;
  font-size: 0.95rem;
}
th, td {
  text-align: left;
  padding: 10px 12px;
  border-bottom: 1px solid var(--line);
  vertical-align: top;
}
th {
  position: sticky;
  top: 0;
  background: #f6faff;
  color: #26435f;
}
.sort-header {
  display: inline-flex;
  align-items: center;
  gap: 6px;
  padding: 0;
  border: 0;
  background: transparent;
  color: inherit;
  font: inherit;
  font-weight: 700;
  cursor: pointer;
}
.sort-header .sort-arrow {
  color: var(--muted);
  font-size: 0.9em;
}
.sort-header.active {
  color: var(--accent);
}
.sort-header.active .sort-arrow {
  color: var(--accent);
}
.status {
  display: inline-block;
  padding: 4px 10px;
  border-radius: 999px;
  font-size: 0.82rem;
  font-weight: 700;
}
.status.ok, .status.skipped_existing { background: #e8f8f1; color: var(--ok); }
.status.pending, .status.running { background: #fff4d7; color: #8a5b00; }
.status.error { background: #fdebed; color: var(--error); }
.structure-badge {
  display: inline-flex;
  align-items: center;
  justify-content: center;
  min-width: 46px;
  padding: 4px 10px;
  border-radius: 999px;
  font-size: 0.78rem;
  font-weight: 800;
  letter-spacing: 0.02em;
}
.structure-badge.yes {
  background: rgba(15, 154, 112, 0.14);
  color: #0f6c53;
}
.structure-badge.no {
  background: rgba(120, 137, 160, 0.14);
  color: #5a6c80;
}
.link-group { display: flex; gap: 10px; flex-wrap: wrap; align-items: center; }
.link-group-wrap {
  margin-top: 18px;
}
.button {
  display: inline-flex;
  align-items: center;
  justify-content: center;
  text-decoration: none;
  color: white;
  background: linear-gradient(135deg, var(--accent), var(--accent-2));
  border: 0;
  padding: 11px 16px;
  border-radius: 999px;
  font-weight: 700;
  font: inherit;
  box-shadow: 0 12px 28px rgba(15, 93, 124, 0.28);
  cursor: pointer;
}
.button:disabled,
.mini-button:disabled {
  opacity: 0.62;
  cursor: progress;
}
.mini-button {
  display: inline-flex;
  align-items: center;
  gap: 6px;
  border: 1px solid var(--line);
  background: #fff;
  color: var(--accent);
  padding: 6px 10px;
  border-radius: 999px;
  font: inherit;
  font-size: 0.88rem;
  cursor: pointer;
}
.mini-button.active {
  background: var(--accent);
  border-color: var(--accent);
  color: #fff;
}
.mini-button:disabled {
  cursor: not-allowed;
  opacity: 0.48;
}
.construct-actions {
  margin: 4px 0 14px;
}
.view-3d-btn {
  display: inline-flex;
  align-items: center;
  gap: 8px;
  border: 1px solid var(--accent);
  background: var(--accent);
  color: #fff;
  padding: 9px 16px;
  border-radius: 999px;
  font: inherit;
  font-size: 0.92rem;
  font-weight: 600;
  line-height: 1;
  cursor: pointer;
  box-shadow: 0 1px 3px rgba(0, 0, 0, 0.16);
  transition: background 0.15s ease, box-shadow 0.15s ease, transform 0.04s ease;
}
.view-3d-btn:hover {
  background: var(--accent-2);
  border-color: var(--accent-2);
  box-shadow: 0 3px 9px rgba(0, 0, 0, 0.22);
}
.view-3d-btn:active {
  transform: translateY(1px);
}
.view-3d-btn:focus-visible {
  outline: 2px solid var(--accent-2);
  outline-offset: 2px;
}
.view-3d-btn svg {
  width: 17px;
  height: 17px;
  flex: none;
}
.view-3d-btn--compact {
  padding: 6px 12px;
  font-size: 0.84rem;
  box-shadow: none;
}
.view-3d-btn--compact svg {
  width: 15px;
  height: 15px;
}
.button-row {
  display: inline-flex;
  flex-wrap: wrap;
  gap: 8px;
}
.inline-link {
  display: inline-flex;
  align-items: center;
  gap: 4px;
  margin-right: 8px;
  color: var(--accent);
  font-weight: 700;
  text-decoration: none;
}
.inline-link:hover {
  text-decoration: underline;
}
.muted {
  color: var(--muted);
}
.meta {
  display: grid;
  grid-template-columns: max-content 1fr;
  gap: 10px 14px;
  margin: 0;
}
.meta dt { font-weight: 700; }
.meta dd { margin: 0; color: var(--muted); }
.card h2 {
  margin: 0 0 16px;
  font-size: 1.08rem;
  letter-spacing: -0.01em;
}
.viewer-layout {
  display: grid;
  grid-template-columns: minmax(340px, 0.9fr) minmax(320px, 1.1fr);
  gap: 18px;
}
.viewer-panel, .sequence-panel {
  display: grid;
  gap: 12px;
}
.viewer-toolbar {
  display: flex;
  gap: 10px;
  align-items: center;
  flex-wrap: wrap;
}
.viewer-note, .viewer-info {
  color: var(--muted);
}
.structure-viewer {
  position: relative;
  width: 100%;
  min-height: 420px;
  border: 1px solid var(--line);
  border-radius: 16px;
  background: #fff;
  overflow: hidden;
}
.interactive-plots {
  display: grid;
  grid-template-columns: minmax(0, 1fr) minmax(0, 1fr);
  gap: 18px;
  margin-top: 8px;
}
.topology-legend {
  display: flex;
  flex-wrap: wrap;
  gap: 8px 12px;
  margin-top: 8px;
  color: var(--muted);
  font-size: 0.82rem;
}
.topology-legend-item {
  display: inline-flex;
  align-items: center;
  gap: 6px;
  white-space: nowrap;
}
.topology-swatch {
  display: inline-block;
  width: 12px;
  height: 12px;
  border-radius: 999px;
  border: 1px solid rgba(0, 0, 0, 0.12);
}
.gpcr-swatch {
  background: #1f1f1f;
  border-radius: 2px;
}
.selection-live-card {
  margin-top: 18px;
  padding: 16px;
  border: 1px solid var(--line);
  border-radius: 16px;
  background: #fff;
}
.selection-live-header {
  display: flex;
  justify-content: space-between;
  gap: 12px;
  align-items: center;
  flex-wrap: wrap;
  margin-bottom: 10px;
}
.selection-live-header h3 {
  margin: 0;
}
.selection-export-help {
  margin: 0 0 10px;
  font-size: 0.85rem;
}
.selection-table-wrap {
  overflow: auto;
}
.selection-warnings {
  margin: 10px 0 14px;
}
.selection-warning-block {
  padding: 10px 12px;
  border-radius: 12px;
  border: 1px solid var(--line);
  background: #fcfaf6;
  margin-bottom: 10px;
}
.selection-warning-block.cysteine {
  border-color: #d5b26b;
  background: #fff7e6;
}
.selection-warning-block.furin {
  border-color: #d39a9a;
  background: #fff1f1;
}
.selection-warning-block strong {
  display: block;
  margin-bottom: 6px;
}
.selection-warning-block ul {
  margin: 0;
  padding-left: 18px;
}
.mutation-controls {
  display: grid;
  gap: 8px;
  margin-top: 12px;
  padding: 12px;
  border: 1px solid var(--line);
  border-radius: 12px;
  background: #fcfdff;
}
.mutation-checkboxes {
  display: flex;
  flex-wrap: wrap;
  gap: 10px 14px;
}
.mutation-checkbox {
  display: inline-flex;
  align-items: center;
  gap: 8px;
  font: 600 0.88rem/1.2 ui-sans-serif, system-ui, sans-serif;
  color: #24435f;
}
.mutation-checkbox input {
  margin: 0;
}
.selection-tsv {
  width: 100%;
  min-height: 110px;
  margin-top: 12px;
  padding: 10px 12px;
  border: 1px solid var(--line);
  border-radius: 12px;
  font: 0.82rem/1.4 ui-monospace, SFMono-Regular, Menlo, monospace;
  background: #fcfaf6;
  color: var(--ink);
}
.plot-card {
  padding: 16px;
  border: 1px solid var(--line);
  border-radius: 16px;
  background: #fff;
}
.plot-card-header {
  display: flex;
  justify-content: space-between;
  align-items: center;
  gap: 12px;
  flex-wrap: wrap;
  margin-bottom: 12px;
}
.plot-card h3 {
  margin: 0;
}
.plot-controls {
  display: flex;
  gap: 8px;
  flex-wrap: wrap;
}
.plot-svg {
  width: 100%;
  height: 240px;
  display: block;
  border: 1px solid var(--line);
  border-radius: 12px;
  background: white;
  cursor: crosshair;
}
.pae-canvas {
  width: 100%;
  max-width: 520px;
  aspect-ratio: 1 / 1;
  display: block;
  border: 1px solid var(--line);
  border-radius: 12px;
  background: white;
  cursor: crosshair;
}
.sequence-header {
  display: flex;
  flex-wrap: wrap;
  justify-content: space-between;
  gap: 8px;
  color: var(--muted);
}
.sequence-viewer {
  border: 1px solid var(--line);
  border-radius: 16px;
  background: #fff;
  padding: 12px;
  max-height: 480px;
  overflow: auto;
}
.sequence-row {
  display: grid;
  grid-template-columns: 48px 1fr;
  gap: 10px;
  align-items: start;
  margin-bottom: 8px;
}
.sequence-index {
  color: var(--muted);
  font: 600 0.82rem/1.8 ui-monospace, SFMono-Regular, Menlo, monospace;
  text-align: right;
}
.sequence-residues {
  display: flex;
  flex-wrap: wrap;
  gap: 2px;
}
.residue {
  border: 0;
  background: transparent;
  color: #3e2f1c;
  padding: 2px 3px;
  border-radius: 6px;
  font: 600 0.8rem/1 ui-monospace, SFMono-Regular, Menlo, monospace;
  cursor: pointer;
}
.residue:hover {
  background: #f7efe3;
}
.residue.ecto-residue {
  background: #eef6ff;
}
.residue.active-range {
  background: #dbeadf;
  color: #1f5c3e;
}
.residue.active-residue {
  background: #fce7ea;
  color: #a61e2e;
}
.construct-list { display: grid; gap: 16px; }
.report-page { scroll-behavior: smooth; }
[data-report-section] { scroll-margin-top: 16px; }
.report-sections-button {
  position: fixed;
  right: calc(16px + env(safe-area-inset-right));
  bottom: calc(16px + env(safe-area-inset-bottom));
  z-index: 20;
  min-height: 44px;
  padding: 10px 18px;
  border: 1px solid #fff;
  border-radius: 24px;
  background: #0D1F3C;
  color: #fff;
  font: inherit;
  font-weight: 600;
  box-shadow: 0 4px 18px rgba(13, 31, 60, 0.2);
  cursor: pointer;
}
#report-sections-menu {
  position: fixed;
  inset: auto calc(16px + env(safe-area-inset-right)) calc(72px + env(safe-area-inset-bottom)) auto;
  margin: 0;
  width: min(300px, calc(100vw - 32px - env(safe-area-inset-left) - env(safe-area-inset-right)));
  max-height: min(70dvh, calc(100dvh - 96px - env(safe-area-inset-bottom)));
  overflow-y: auto;
  overscroll-behavior: contain;
  padding: 8px;
  border: 1px solid var(--line);
  border-radius: 16px;
  background: #fff;
  color: var(--ink);
  box-shadow: 0 8px 32px rgba(13, 31, 60, 0.2);
}
#report-sections-menu nav { display: grid; gap: 2px; }
#report-sections-menu a {
  display: flex;
  align-items: center;
  min-height: 44px;
  padding: 10px 12px;
  border-radius: 8px;
  color: inherit;
  text-decoration: none;
}
#report-sections-menu a:hover,
#report-sections-menu a[aria-current="location"] { background: #dff7f5; color: #075e5a; }
.report-sections-button:focus-visible,
#report-sections-menu a:focus-visible { outline: 2px solid #007FA3; outline-offset: 2px; }
@media (prefers-reduced-motion: reduce) {
  .report-page { scroll-behavior: auto; }
}
@media print {
  .report-sections-button, #report-sections-menu { display: none; }
}
.construct-tabs {
  display: flex;
  gap: 8px;
  overflow-x: auto;
  padding: 4px 4px 12px;
  margin-bottom: 12px;
}
.construct-tabs [role="tab"] {
  flex: 0 0 auto;
  padding: 10px 14px;
  border: 1px solid var(--line);
  border-radius: 10px;
  background: var(--card);
  color: var(--ink);
  font: inherit;
  cursor: pointer;
}
.construct-tabs [aria-selected="true"] {
  background: #0D1F3C;
  color: #fff;
}
.construct-tabs [role="tab"]:focus-visible {
  outline: 2px solid #007FA3;
  outline-offset: 2px;
}
#construct-details [role="tabpanel"][hidden] { display: none; }
#construct-details .construct-list,
#construct-details .construct-card,
#construct-details .construct-export-card {
  grid-template-columns: minmax(0, 1fr);
}
#construct-details .construct-copy > h3,
#construct-details .construct-copy > p,
#construct-details .construct-grid li {
  overflow-wrap: anywhere;
}
.construct-card {
  display: grid;
  grid-template-columns: 1fr;
  gap: 18px;
  padding: 18px;
  background: linear-gradient(180deg, #ffffff 0%, #f9fcff 100%);
}
.construct-card h3 { margin-top: 0; }
.construct-grid {
  display: grid;
  grid-template-columns: 1fr;
  gap: 12px 18px;
}
.construct-export-card {
  display: grid;
  gap: 12px;
  margin: 12px 0 16px;
  padding: 14px;
  border: 1px solid var(--line);
  border-radius: 14px;
  background: #fcfdff;
}
.construct-export-card h4 {
  margin: 0;
}
.export-box-grid {
  display: grid;
  grid-template-columns: 1fr;
  gap: 12px;
}
.export-label {
  display: inline-block;
  margin-bottom: 6px;
  font-weight: 700;
  color: var(--muted);
}
.construct-export-box {
  min-height: 132px;
  margin-top: 0;
}
.construct-sequence-table code.sequence {
  max-height: 96px;
  overflow: auto;
}
.construct-grid h4 { margin-bottom: 6px; }
.construct-card ul { margin: 0; padding-left: 18px; }
.construct-images {
  display: grid;
  grid-template-columns: repeat(2, minmax(0, 1fr));
  gap: 12px;
  align-content: start;
}
.construct-images img {
  width: 100%;
  border-radius: 12px;
  border: 1px solid var(--line);
  background: white;
}
.alignment-meta {
  display: flex;
  flex-wrap: wrap;
  gap: 10px 16px;
  margin: 8px 0 10px;
  color: var(--muted);
  font-size: 0.88rem;
}
.alignment-legend {
  display: inline-flex;
  align-items: center;
  gap: 7px;
  margin: 0 0 8px;
  color: var(--muted);
  font-size: 0.82rem;
}
.alignment-residue {
  display: inline-block;
  width: 1ch;
  min-width: 1ch;
  text-align: center;
  font: inherit;
  letter-spacing: 0;
  box-sizing: content-box;
}
.alignment-surface-ecto {
  border-radius: 3px;
  background: #f7c948;
  color: #2f2300;
  font-weight: 800;
}
.alignment-block {
  margin: 0;
  padding: 10px 12px;
  border: 1px solid var(--line);
  border-radius: 12px;
  background: #fcfaf6;
  font: 0.8rem/1.45 ui-monospace, SFMono-Regular, Menlo, monospace;
  letter-spacing: 0;
  tab-size: 4;
  white-space: pre;
  overflow: auto;
}
.alignment-block * {
  font: inherit;
  letter-spacing: inherit;
}
.hit-description {
  margin-top: 4px;
  max-width: 36ch;
}
.cross-reactivity-scroll details summary {
  cursor: pointer;
  font-weight: 700;
}
code.sequence {
  display: inline-block;
  white-space: pre-wrap;
  word-break: break-word;
  font-size: 0.86rem;
}
.matrix-wrap {
  overflow: auto;
  margin-top: 14px;
}
.matrix-wrap table {
  min-width: max-content;
}
.matrix-wrap tbody th {
  position: sticky;
  left: 0;
  z-index: 1;
  background: #f6faff;
}
.matrix-cell {
  text-align: center;
  white-space: nowrap;
  min-width: 84px;
}
.matrix-cell-empty {
  background: #f9fbfe;
  color: #94a6b8;
  box-shadow: inset 0 0 0 1px #e4edf5;
}
.table-scroll {
  overflow: auto;
}
.cross-reactivity-scroll {
  max-height: 420px;
}
.site-footer {
  width: 100%;
  max-width: 100%;
  margin: 34px 0 0;
  padding: 26px;
  border: 1px solid var(--line);
  border-radius: 24px;
  background:
    linear-gradient(180deg, rgba(255,255,255,0.92), rgba(248,251,255,0.96)),
    radial-gradient(circle at top right, rgba(14, 140, 136, 0.1), transparent 38%);
  box-shadow: 0 18px 50px rgba(10, 35, 66, 0.08);
  overflow: hidden;
}
.footer-grid {
  display: flex;
  align-items: flex-start;
  justify-content: space-between;
  gap: 24px;
  flex-wrap: wrap;
}
.footer-brand {
  display: flex;
  align-items: center;
  gap: 16px;
  max-width: 620px;
}
.footer-logo-link {
  display: inline-flex;
  align-items: center;
  flex: 0 0 auto;
}
.footer-logo-link:hover {
  opacity: 0.82;
}
.footer-logo {
  display: block;
  height: 42px;
  width: auto;
}
.footer-brand strong {
  display: block;
  font-size: 1.08rem;
  letter-spacing: -0.01em;
}
.footer-brand p {
  margin: 4px 0 0;
  color: var(--muted);
  line-height: 1.45;
}
.footer-links {
  display: flex;
  flex-wrap: wrap;
  gap: 8px;
  justify-content: flex-end;
}
.footer-links a {
  display: inline-flex;
  align-items: center;
  min-height: 34px;
  padding: 8px 12px;
  border: 1px solid var(--line);
  border-radius: 999px;
  background: #ffffff;
  color: var(--accent-2);
  font-weight: 700;
  text-decoration: none;
}
.footer-links a:hover {
  border-color: rgba(14, 140, 136, 0.32);
  color: var(--accent);
}
.footer-meta {
  display: flex;
  flex-wrap: wrap;
  gap: 10px 18px;
  margin-top: 20px;
  padding-top: 16px;
  border-top: 1px solid var(--line);
  color: var(--muted);
  font-size: 0.9rem;
}
.footer-meta a {
  color: var(--accent-2);
  text-decoration: none;
  font-weight: 700;
}
.section-note {
  margin: 0 0 10px;
  color: var(--muted);
}
.index-default-note {
  margin: -2px 0 18px;
  padding: 10px 14px;
  border: 1px solid rgba(14, 140, 136, 0.18);
  border-radius: 14px;
  background: rgba(14, 140, 136, 0.06);
}
a { color: var(--accent-2); }
input::placeholder {
  color: #7d93a8;
}
select, input, button, textarea {
  outline: none;
}
select:focus, input:focus, button:focus, textarea:focus, .residue:focus {
  box-shadow: 0 0 0 3px rgba(14, 140, 136, 0.18);
}
@media (max-width: 980px) {
  .hero-grid, .grid.two, .stack, .construct-card, .viewer-layout, .doc-grid, .download-grid, .calculator-layout { grid-template-columns: 1fr; display: grid; }
  .calculator-row { grid-template-columns: 1fr; min-height: 0; }
  .calculator-card-header { align-items: flex-start; flex-direction: column; }
  .calc-status { white-space: normal; }
  .interactive-plots { grid-template-columns: 1fr; }
  .construct-images { grid-template-columns: 1fr; }
  .stats { grid-template-columns: repeat(2, 1fr); }
  .stat-primary { grid-column: span 2; text-align: center; }
  .footer-brand { align-items: flex-start; }
  .footer-links { justify-content: flex-start; }
  .site-hero { padding: 22px 20px 24px; }
  .brand-bar { flex-wrap: wrap; align-items: center; }
  .brand-actions { align-items: flex-start; justify-content: flex-start; max-width: 100%; width: 100%; }
  .brand-secondary-actions { justify-content: flex-start; }
  .site-nav { justify-content: flex-start; flex-wrap: wrap; border-radius: 18px; }
  .hero-chip { display: none; }
  .brand-logo { height: 40px; }
  .page { width: min(1480px, calc(100vw - 24px)); }
}
.agent-banner {display:grid;grid-template-columns:1fr auto;align-items:center;gap:10px 28px;margin:0 0 16px;padding:20px 24px;border:1px solid #abd8d5;border-left:4px solid #0e8c88;border-radius:14px;background:#f4fcfb;}
.agent-banner .agent-eyebrow {margin:0 0 5px;color:#076b67;font-weight:700;font-size:11px;letter-spacing:.1em;text-transform:uppercase;}
.agent-banner h2 {margin:0 0 6px;font-size:19px;line-height:1.35;letter-spacing:-.25px;}
.agent-banner p {margin:0;color:#415d72;font-size:14px;line-height:1.6;}
.agent-actions {display:flex;align-items:center;gap:18px;flex-wrap:wrap;}
.agent-guide-link {color:#076b67;font-weight:650;font-size:14px;text-underline-offset:4px;}
.agent-copy {border:1px solid #08716d;border-radius:9px;background:#08716d;color:white;padding:11px 16px;font:inherit;font-size:13px;font-weight:650;cursor:pointer;min-height:44px;}
.agent-copy:hover {background:#075e5b;}
.agent-prompt {grid-column:1/-1;font-size:12px;color:#415d72;}
.agent-prompt summary {cursor:pointer;width:fit-content;text-decoration:underline;text-underline-offset:3px;padding:4px 0;}
.agent-prompt label {display:block;margin:10px 0;}
.agent-prompt textarea {display:block;width:100%;max-width:100%;margin:0 0 8px;padding:12px;border:1px solid #abd8d5;border-radius:8px;font-family:inherit;font-size:13px;line-height:1.6;background:white;color:#102338;}
.agent-copy-status {grid-column:1/-1;color:#076b67;font-size:12px;}
.agent-copy-status:empty {display:none;}
.agent-essentials {margin:0 0 16px;padding:0 4px;color:#415d72;font-size:12px;line-height:1.6;}
.agent-essentials a {color:#076b67;text-underline-offset:3px;}
.guide-layout {display:grid;grid-template-columns:220px minmax(0,1fr);gap:28px;align-items:start;}
.guide-toc {position:sticky;top:24px;display:grid;gap:16px;padding:20px 0;font-size:13px;line-height:1.5;}
.guide-toc strong {font-size:11px;letter-spacing:.08em;text-transform:uppercase;color:#587089;}
.guide-toc a {color:#1f4d8c;text-decoration:none;}
.guide-toc a:hover {text-decoration:underline;}
.guide-content {padding:30px 36px;font-size:15px;line-height:1.75;overflow-wrap:anywhere;}
.guide-content h2 {font-size:22px;margin:36px 0 14px;scroll-margin-top:24px;line-height:1.3;}
.guide-content li {margin:10px 0;}
.guide-content ul,.guide-content ol {padding-left:24px;}
.guide-content a {color:#076b67;text-underline-offset:3px;}
.guide-content blockquote {margin:24px 0;border-left:3px solid #abd8d5;padding:1px 20px;color:#415d72;}
.guide-content code {background:#eff4f7;padding:2px 4px;border-radius:3px;font-size:.9em;}
.agent-banner a:focus-visible,.agent-banner button:focus-visible,.agent-banner summary:focus-visible,.agent-banner textarea:focus-visible,.guide-toc a:focus-visible {outline:3px solid #1f4d8c;outline-offset:4px;}
@media(max-width:800px) {
 .agent-banner {grid-template-columns:minmax(0,1fr);padding:18px;gap:14px;}
 .agent-actions {gap:14px;}
 .agent-banner h2 {font-size:18px;}
 .guide-layout {grid-template-columns:minmax(0,1fr);gap:16px;}
 .guide-toc {position:static;padding:0;gap:10px;}
 .guide-content {padding:20px;font-size:14px;}
}
.footer-citation {margin:20px 0 0;color:var(--muted);line-height:1.7;font-size:.9rem;}
.footer-citation a {color:var(--accent-2);text-underline-offset:3px;}
.paper-citation {overflow-wrap:anywhere;}
#cite-openantigens {scroll-margin-top:24px;}

""" + theme_overrides + _construct_workbench_css()
