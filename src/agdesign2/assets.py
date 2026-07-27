from __future__ import annotations

import json
import os
import types
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import fields, is_dataclass
from pathlib import Path
from time import perf_counter
from typing import Any, get_args, get_origin, get_type_hints

from .config import AnalysisConfig
from .af3 import find_local_af3_artifacts
from .models import AnalysisReport
from .pipeline import AntigenAnalyzer
from .structure_utils import (
    load_pae_matrix,
    pae_matches_sequence_length,
    parse_alphafold_pdb,
    residue_names_to_sequence,
)


def generate_assets_for_summary(
    summary_path: str | Path,
    *,
    config: AnalysisConfig,
    jobs: int = 1,
    resume: bool = True,
    verbose: bool = False,
) -> list[dict[str, Any]]:
    summary_file = Path(summary_path)
    if not summary_file.exists():
        return []
    rows = json.loads(summary_file.read_text(encoding="utf-8"))
    tasks = [
        (index, row)
        for index, row in enumerate(rows, start=1)
        if row.get("status") in {"ok", "skipped_existing"} and row.get("json_report")
    ]
    normalized_jobs = _normalize_jobs(jobs)
    results: list[dict[str, Any] | None] = [None] * len(tasks)
    if normalized_jobs <= 1 or len(tasks) <= 1:
        for task_index, (row_index, row) in enumerate(tasks):
            results[task_index] = _generate_assets_worker(config, row, row_index, len(rows), resume, verbose)
        return [item for item in results if item is not None]

    try:
        with ProcessPoolExecutor(max_workers=min(normalized_jobs, len(tasks))) as executor:
            futures = {
                executor.submit(_generate_assets_worker, config, row, row_index, len(rows), resume, verbose): task_index
                for task_index, (row_index, row) in enumerate(tasks)
            }
            for future in as_completed(futures):
                results[futures[future]] = future.result()
    except (OSError, PermissionError):
        for task_index, (row_index, row) in enumerate(tasks):
            results[task_index] = _generate_assets_worker(config, row, row_index, len(rows), resume, verbose)
    return [item for item in results if item is not None]


def generate_assets_for_report_json(
    report_json_path: str | Path,
    *,
    config: AnalysisConfig,
    resume: bool = True,
) -> dict[str, Any]:
    report_path = Path(report_json_path)
    started = perf_counter()
    try:
        payload = json.loads(report_path.read_text(encoding="utf-8"))
        disabled_assets_removed = _remove_disabled_asset_outputs(
            payload=payload,
            report_path=report_path,
            config=config,
        )
        duplicate_assets_removed = _reuse_duplicate_boundary_quality_plots(
            payload=payload,
            report_path=report_path,
        )
        report = _hydrate_dataclass(payload, AnalysisReport)
        if disabled_assets_removed or duplicate_assets_removed:
            report.write_json(report_path)
            report.write_markdown(report_path.with_suffix(".md"))
        accession = report.target.accession
        pdb_path = config.data_dir / "alphafold" / f"{accession}.pdb"
        pae_path = config.data_dir / "alphafold" / f"{accession}.pae.json"
        if (not pdb_path.exists() or not pae_path.exists()) and config.enable_local_af3_fallback:
            fallback = find_local_af3_artifacts(
                accession,
                data_dir=config.data_dir,
                canonical_sequence=report.target.sequence,
            )
            if fallback is not None:
                pdb_path, pae_path = fallback
        if not report.construct_details:
            return {
                "json_report": str(report_path),
                "status": "skipped_no_constructs",
                "duration_seconds": round(perf_counter() - started, 2),
                "error": "",
            }
        if report.ectodomain is None:
            return {
                "json_report": str(report_path),
                "status": "skipped_no_ectodomain",
                "duration_seconds": round(perf_counter() - started, 2),
                "error": "",
            }
        if not pdb_path.exists():
            return {
                "json_report": str(report_path),
                "status": "skipped_missing_alphafold_pdb",
                "duration_seconds": round(perf_counter() - started, 2),
                "error": "",
            }
        if not pae_path.exists():
            return {
                "json_report": str(report_path),
                "status": "skipped_missing_alphafold_pae",
                "duration_seconds": round(perf_counter() - started, 2),
                "error": "",
            }
        structure = parse_alphafold_pdb(pdb_path)
        if residue_names_to_sequence(structure.get("residue_names") or {}) != report.target.sequence:
            _remove_model_dependent_asset_outputs(payload=payload, report_path=report_path)
            report = _hydrate_dataclass(payload, AnalysisReport)
            report.write_json(report_path)
            report.write_markdown(report_path.with_suffix(".md"))
            return {
                "json_report": str(report_path),
                "status": "skipped_incompatible_alphafold",
                "duration_seconds": round(perf_counter() - started, 2),
                "error": "Cached AlphaFold PDB sequence does not match the canonical report sequence.",
            }
        if resume and _report_assets_complete(payload=payload, report_path=report_path, config=config):
            return {
                "json_report": str(report_path),
                "status": "skipped_existing",
                "duration_seconds": round(perf_counter() - started, 2),
                "error": "",
            }
        pae_matrix = load_pae_matrix(pae_path)
        if not pae_matches_sequence_length(pae_matrix, len(report.target.sequence)):
            _remove_model_dependent_asset_outputs(payload=payload, report_path=report_path)
            report = _hydrate_dataclass(payload, AnalysisReport)
            report.write_json(report_path)
            report.write_markdown(report_path.with_suffix(".md"))
            return {
                "json_report": str(report_path),
                "status": "skipped_incompatible_alphafold",
                "duration_seconds": round(perf_counter() - started, 2),
                "error": "Cached AlphaFold PAE dimensions do not match the canonical report sequence.",
            }
        analyzer = AntigenAnalyzer(config=config)
        assets_dir = report_path.parent / f"{report.target.entry_name.lower()}_report_assets"
        analyzer._generate_construct_assets(
            construct_details=report.construct_details,
            assets_dir=assets_dir,
            alphafold_pdb_path=pdb_path,
            residue_plddt=structure["plddt"],
            pae_matrix=pae_matrix,
            ectodomain=report.ectodomain,
        )
        report.write_json(report_path)
        report.write_markdown(report_path.with_suffix(".md"))
        return {
            "json_report": str(report_path),
            "status": "ok",
            "duration_seconds": round(perf_counter() - started, 2),
            "error": "",
        }
    except Exception as exc:
        return {
            "json_report": str(report_path),
            "status": "error",
            "duration_seconds": round(perf_counter() - started, 2),
            "error": str(exc),
        }


def _generate_assets_worker(
    config: AnalysisConfig,
    row: dict[str, Any],
    row_index: int,
    total: int,
    resume: bool,
    verbose: bool,
) -> dict[str, Any]:
    if verbose:
        print(f"[agdesign2-assets] {row_index}/{total} {row.get('query') or row.get('resolved_entry_name')}", flush=True)
    return generate_assets_for_report_json(row["json_report"], config=config, resume=resume)


def _report_assets_complete(*, payload: dict[str, Any], report_path: Path, config: AnalysisConfig) -> bool:
    constructs = payload.get("construct_details") or []
    if not constructs:
        return True
    if not config.render_structure_images and not config.render_quality_plots:
        return True

    target = payload.get("target") if isinstance(payload.get("target"), dict) else {}
    entry_name = str(target.get("entry_name") or "").lower()
    assets_dir = report_path.parent / f"{entry_name}_report_assets" if entry_name else report_path.parent
    ectodomain = payload.get("ectodomain") if isinstance(payload.get("ectodomain"), dict) else None
    ectodomain_start = _coerce_int(ectodomain.get("start")) if ectodomain else None
    ectodomain_end = _coerce_int(ectodomain.get("end")) if ectodomain else None
    expected_paths = []
    for construct in constructs:
        if not isinstance(construct, dict):
            continue
        start = _coerce_int(construct.get("start"))
        end = _coerce_int(construct.get("end"))
        if (
            ectodomain_start is not None
            and ectodomain_end is not None
            and start is not None
            and end is not None
            and (start < ectodomain_start or end > ectodomain_end)
        ):
            continue
        safe_name = _safe_filename(str(construct.get("name") or "construct"))
        if config.render_structure_images:
            value = construct.get("structure_image")
            expected_paths.append(report_path.parent / str(value) if value else assets_dir / f"{safe_name}_structure.png")
        if config.render_quality_plots:
            value = construct.get("quality_plot")
            expected_paths.append(report_path.parent / str(value) if value else assets_dir / f"{safe_name}_quality.png")
    return bool(expected_paths) and all(path.exists() for path in expected_paths)


def _remove_disabled_asset_outputs(
    *,
    payload: dict[str, Any],
    report_path: Path,
    config: AnalysisConfig,
) -> bool:
    constructs = payload.get("construct_details") or []
    changed = False
    for construct in constructs:
        if not isinstance(construct, dict):
            continue
        if not config.render_structure_images and construct.get("structure_image") is not None:
            construct["structure_image"] = None
            changed = True
        if not config.render_quality_plots and construct.get("quality_plot") is not None:
            construct["quality_plot"] = None
            changed = True

    target = payload.get("target") if isinstance(payload.get("target"), dict) else {}
    entry_name = str(target.get("entry_name") or "").lower()
    if not entry_name:
        return changed
    assets_dir = report_path.parent / f"{entry_name}_report_assets"
    patterns = []
    if not config.render_structure_images:
        patterns.append("*_structure.png")
    if not config.render_quality_plots:
        patterns.append("*_quality.png")
    for pattern in patterns:
        for path in assets_dir.glob(pattern):
            if path.is_file():
                path.unlink()
                changed = True
    return changed


def _remove_model_dependent_asset_outputs(
    *,
    payload: dict[str, Any],
    report_path: Path,
) -> bool:
    constructs = payload.get("construct_details") or []
    changed = False
    for construct in constructs:
        if not isinstance(construct, dict):
            continue
        for field_name in ("structure_image", "quality_plot"):
            if construct.get(field_name) is not None:
                construct[field_name] = None
                changed = True

    target = payload.get("target") if isinstance(payload.get("target"), dict) else {}
    entry_name = str(target.get("entry_name") or "").lower()
    if not entry_name:
        return changed
    assets_dir = report_path.parent / f"{entry_name}_report_assets"
    for pattern in ("*_structure.png", "*_quality.png"):
        for path in assets_dir.glob(pattern):
            if path.is_file():
                path.unlink()
                changed = True
    return changed


def _reuse_duplicate_boundary_quality_plots(
    *,
    payload: dict[str, Any],
    report_path: Path,
) -> bool:
    constructs = [
        construct
        for construct in (payload.get("construct_details") or [])
        if isinstance(construct, dict)
    ]
    by_boundary: dict[tuple[int, int], list[dict[str, Any]]] = {}
    for construct in constructs:
        start = _coerce_int(construct.get("start"))
        end = _coerce_int(construct.get("end"))
        if start is None or end is None:
            continue
        by_boundary.setdefault((start, end), []).append(construct)

    old_paths: set[str] = set()
    final_paths: set[str] = set()
    changed = False
    for group in by_boundary.values():
        paths = [
            str(construct.get("quality_plot"))
            for construct in group
            if construct.get("quality_plot")
        ]
        old_paths.update(paths)
        if len(group) <= 1 or not paths:
            final_paths.update(paths)
            continue
        canonical = next(
            (
                value
                for value in paths
                if (report_path.parent / value).is_file()
            ),
            paths[0],
        )
        final_paths.add(canonical)
        for construct in group:
            if construct.get("quality_plot") != canonical:
                construct["quality_plot"] = canonical
                changed = True

    target = payload.get("target") if isinstance(payload.get("target"), dict) else {}
    entry_name = str(target.get("entry_name") or "").lower()
    if not entry_name:
        return changed
    assets_dir = (report_path.parent / f"{entry_name}_report_assets").resolve()
    for value in old_paths - final_paths:
        path = (report_path.parent / value).resolve()
        if (
            path.parent == assets_dir
            and path.name.endswith("_quality.png")
            and path.is_file()
        ):
            path.unlink()
            changed = True
    return changed


def _safe_filename(text: str) -> str:
    return "".join(ch if ch.isalnum() or ch in {"-", "_"} else "_" for ch in text)


def _coerce_int(value: Any) -> int | None:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _normalize_jobs(jobs: int | None) -> int:
    if jobs is None or jobs <= 0:
        return max(1, min(4, os.cpu_count() or 1))
    return max(1, jobs)


def _hydrate_dataclass(value: Any, annotation: Any) -> Any:
    if value is None:
        return None
    origin = get_origin(annotation)
    if origin in {list, tuple}:
        args = get_args(annotation)
        item_type = args[0] if args else Any
        return [_hydrate_dataclass(item, item_type) for item in value]
    if origin is dict:
        return value
    if origin in {types.UnionType, getattr(types, "UnionType", None)} or str(origin) == "typing.Union":
        for arg in get_args(annotation):
            if arg is type(None):
                continue
            try:
                return _hydrate_dataclass(value, arg)
            except Exception:
                continue
        return value
    if annotation is Any:
        return value
    if isinstance(annotation, type) and is_dataclass(annotation):
        hints = get_type_hints(annotation)
        kwargs = {}
        for field in fields(annotation):
            if isinstance(value, dict) and field.name in value:
                kwargs[field.name] = _hydrate_dataclass(value[field.name], hints.get(field.name, Any))
        return annotation(**kwargs)
    return value
