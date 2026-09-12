from __future__ import annotations

import csv
import json
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from time import perf_counter
from typing import Any

from .config import AnalysisConfig
from .http import _atomic_write
from .pipeline import AntigenAnalyzer


@dataclass(slots=True)
class BatchRow:
    query: str
    source_column: str
    record: dict[str, str]


def load_batch_rows(tsv_path: str | Path) -> list[BatchRow]:
    path = Path(tsv_path)
    rows: list[BatchRow] = []
    for record in _load_input_records(path):
        query, source_column = _select_query(record)
        if path.suffix.lower() == ".csv" and str(record.get("uniprot_name") or "").strip():
            query, source_column = str(record["uniprot_name"]).strip(), "uniprot_name"
        if query is None or source_column is None:
            continue
        rows.append(BatchRow(query=query, source_column=source_column, record=record))
    return rows


def _load_input_records(path: Path) -> list[dict[str, str]]:
    delimiter = "," if path.suffix.lower() == ".csv" else "\t"
    with path.open("r", encoding="utf-8", newline="") as handle:
        records = list(csv.DictReader(handle, delimiter=delimiter))
    if path.suffix.lower() == ".csv":
        for record in records:
            for field, alias in (
                ("uniprot_id", "accession"),
                ("uniprot_name", "uniprot_entry_name"),
                ("gene", "gene_symbol"),
            ):
                if not str(record.get(field) or "").strip():
                    record[field] = str(record.get(alias) or "")
    return records


def run_batch_from_tsv(
    *,
    analyzer: AntigenAnalyzer,
    tsv_path: str | Path,
    output_dir: str | Path,
    limit: int | None = None,
    resume: bool = True,
    jobs: int = 1,
) -> list[dict[str, str | int | float | None]]:
    rows = load_batch_rows(tsv_path)
    if limit is not None:
        rows = rows[:limit]

    output_root = Path(output_dir)
    output_root.mkdir(parents=True, exist_ok=True)
    summary_path = output_root / "batch_summary.json"
    dashboard_path = output_root / "batch_summary.md"
    results: list[dict[str, str | int | float | None] | None] = [None] * len(rows)

    total = len(rows)
    normalized_jobs = _normalize_jobs(jobs)
    if normalized_jobs <= 1 or total <= 1:
        for index, row in enumerate(rows, start=1):
            result = _run_batch_row(
                analyzer=analyzer,
                row=row,
                index=index,
                total=total,
                output_root=output_root,
                resume=resume,
            )
            results[index - 1] = result
            if index & (index - 1) == 0:
                _write_partial_batch_outputs(results, summary_path=summary_path, dashboard_path=dashboard_path)
        _write_partial_batch_outputs(results, summary_path=summary_path, dashboard_path=dashboard_path)
        return [item for item in results if item is not None]

    config = analyzer.config
    pending_jobs: list[tuple[int, BatchRow]] = []
    for index, row in enumerate(rows, start=1):
        expected_stem = _expected_stem(row)
        json_path = output_root / f"{expected_stem}.json"
        markdown_path = output_root / f"{expected_stem}.md"
        if resume and json_path.exists() and markdown_path.exists():
            result = _base_batch_result(row=row, index=index, total=total)
            result["status"] = "skipped_existing"
            result["json_report"] = str(json_path)
            result["markdown_report"] = str(markdown_path)
            result["duration_seconds"] = 0.0
            results[index - 1] = result
        else:
            pending_jobs.append((index, row))

    completed = total - len(pending_jobs)
    _write_partial_batch_outputs(results, summary_path=summary_path, dashboard_path=dashboard_path)
    try:
        with ProcessPoolExecutor(
            max_workers=min(normalized_jobs, len(pending_jobs) or 1),
            initializer=_init_batch_worker,
            initargs=(config,),
        ) as executor:
            futures = {
                executor.submit(
                    _run_batch_worker,
                    config,
                    row,
                    index,
                    total,
                    output_root,
                ): index
                for index, row in pending_jobs
            }
            for future in as_completed(futures):
                index = futures[future]
                results[index - 1] = future.result()
                completed += 1
                if completed & (completed - 1) == 0:
                    _write_partial_batch_outputs(results, summary_path=summary_path, dashboard_path=dashboard_path)
    except (OSError, PermissionError):
        for index, row in pending_jobs:
            if results[index - 1] is not None:
                continue
            results[index - 1] = _run_batch_worker(config, row, index, total, output_root)
            completed += 1
            if completed & (completed - 1) == 0:
                _write_partial_batch_outputs(results, summary_path=summary_path, dashboard_path=dashboard_path)
    _write_partial_batch_outputs(results, summary_path=summary_path, dashboard_path=dashboard_path)
    return [item for item in results if item is not None]


def _normalize_jobs(jobs: int | None) -> int:
    if jobs is None or jobs <= 0:
        return max(1, min(4, os.cpu_count() or 1))
    return max(1, jobs)


def _base_batch_result(*, row: BatchRow, index: int, total: int) -> dict[str, str | int | float | None]:
    base_record = {key: value for key, value in row.record.items()}
    return {
        **base_record,
        "query": row.query,
        "source_column": row.source_column,
        "batch_index": index,
        "batch_total": total,
        "status": "started",
        "resolved_entry_name": None,
        "json_report": None,
        "markdown_report": None,
        "duration_seconds": None,
        "error": None,
    }


def _run_batch_row(
    *,
    analyzer: AntigenAnalyzer,
    row: BatchRow,
    index: int,
    total: int,
    output_root: Path,
    resume: bool,
) -> dict[str, str | int | float | None]:
    result = _base_batch_result(row=row, index=index, total=total)
    print(f"[agdesign2-batch] {index}/{total} {row.query} ({row.source_column})", flush=True)
    started = perf_counter()
    try:
        expected_stem = _expected_stem(row)
        json_path = output_root / f"{expected_stem}.json"
        markdown_path = output_root / f"{expected_stem}.md"
        if resume and json_path.exists() and markdown_path.exists():
            result["status"] = "skipped_existing"
            result["json_report"] = str(json_path)
            result["markdown_report"] = str(markdown_path)
        else:
            report = analyzer.analyze_target(row.query, output_dir=output_root)
            stem = f"{report.target.entry_name.lower()}_report"
            result["status"] = "ok"
            result["resolved_entry_name"] = report.target.entry_name
            result["json_report"] = str(output_root / f"{stem}.json")
            result["markdown_report"] = str(output_root / f"{stem}.md")
    except Exception as exc:
        result["status"] = "error"
        result["error"] = str(exc)
    result["duration_seconds"] = round(perf_counter() - started, 2)
    return result


# One analyzer per worker process, built once by the pool initializer and
# reused across every row that worker handles. Building it per row re-parsed the
# ortholog/family/paralog/topology reference files for every target. The
# sequential path already reuses a single analyzer, so this matches that
# proven-safe pattern (the reference caches are read-only after load).
_WORKER_ANALYZER: "AntigenAnalyzer | None" = None


def _init_batch_worker(config: AnalysisConfig) -> None:
    global _WORKER_ANALYZER
    _WORKER_ANALYZER = AntigenAnalyzer(config=config)


def _run_batch_worker(
    config: AnalysisConfig,
    row: BatchRow,
    index: int,
    total: int,
    output_root: Path,
) -> dict[str, str | int | float | None]:
    # Reuse the worker-local analyzer when the pool initializer created one;
    # fall back to a fresh analyzer on the in-process (no-pool) error path.
    analyzer = _WORKER_ANALYZER if _WORKER_ANALYZER is not None else AntigenAnalyzer(config=config)
    return _run_batch_row(
        analyzer=analyzer,
        row=row,
        index=index,
        total=total,
        output_root=output_root,
        resume=False,
    )


def _write_partial_batch_outputs(
    results: list[dict[str, str | int | float | None] | None],
    *,
    summary_path: Path,
    dashboard_path: Path,
) -> None:
    written = [item for item in results if item is not None]
    _atomic_write(summary_path, (json.dumps(written, indent=2) + "\n").encode("utf-8"))
    _atomic_write(dashboard_path, render_batch_summary_markdown(written).encode("utf-8"))


def retry_batch_from_tsv(
    *,
    analyzer: AntigenAnalyzer,
    tsv_path: str | Path,
    output_dir: str | Path,
    summary_path: str | Path | None = None,
    require_alphafold: bool = True,
) -> list[dict[str, str | int | float | None]]:
    rows = load_batch_rows(tsv_path)
    output_root = Path(output_dir)
    output_root.mkdir(parents=True, exist_ok=True)
    summary_file = Path(summary_path) if summary_path is not None else output_root / "batch_summary.json"
    dashboard_path = output_root / "batch_summary.md"

    existing_results = []
    if summary_file.exists():
        existing_results = json.loads(summary_file.read_text(encoding="utf-8"))
    # Key by (query, source_column) — the same identity the resume path and
    # refresh job use. Keying by query alone collided when two distinct rows
    # produced the same query string, letting one row inherit the other's
    # status/report/AlphaFold-retry decision.
    existing_by_key = {
        (str(item.get("query")), str(item.get("source_column"))): item
        for item in existing_results
        if item.get("query")
    }

    total = len(rows)
    merged_results: list[dict[str, str | int | float | None]] = []
    retry_targets: list[tuple[int, BatchRow]] = []

    for index, row in enumerate(rows, start=1):
        base_record = {key: value for key, value in row.record.items()}
        current = dict(existing_by_key.get((str(row.query), str(row.source_column))) or {})
        merged = {
            **base_record,
            "query": row.query,
            "source_column": row.source_column,
            "batch_index": index,
            "batch_total": total,
            "status": current.get("status", "pending"),
            "resolved_entry_name": current.get("resolved_entry_name"),
            "json_report": current.get("json_report"),
            "markdown_report": current.get("markdown_report"),
            "duration_seconds": current.get("duration_seconds"),
            "error": current.get("error"),
        }
        merged_results.append(merged)
        if _needs_retry(merged, output_root=output_root, require_alphafold=require_alphafold):
            retry_targets.append((index, row))

    for index, row in retry_targets:
        print(f"[agdesign2-retry] {index}/{total} {row.query} ({row.source_column})", flush=True)
        started = perf_counter()
        result = merged_results[index - 1]
        try:
            report = analyzer.analyze_target(row.query, output_dir=output_root)
            stem = f"{report.target.entry_name.lower()}_report"
            result["status"] = "ok"
            result["resolved_entry_name"] = report.target.entry_name
            result["json_report"] = str(output_root / f"{stem}.json")
            result["markdown_report"] = str(output_root / f"{stem}.md")
            result["error"] = None
        except Exception as exc:
            result["status"] = "error"
            result["error"] = str(exc)
        result["duration_seconds"] = round(perf_counter() - started, 2)
        summary_file.write_text(json.dumps(merged_results, indent=2) + "\n", encoding="utf-8")
        dashboard_path.write_text(render_batch_summary_markdown(merged_results), encoding="utf-8")
    if not retry_targets:
        summary_file.write_text(json.dumps(merged_results, indent=2) + "\n", encoding="utf-8")
        dashboard_path.write_text(render_batch_summary_markdown(merged_results), encoding="utf-8")
    return merged_results


def rebuild_batch_summary_from_tsv(
    *,
    tsv_path: str | Path,
    output_dir: str | Path,
    summary_path: str | Path | None = None,
) -> list[dict[str, str | int | float | None]]:
    rows = load_batch_rows(tsv_path)
    output_root = Path(output_dir)
    output_root.mkdir(parents=True, exist_ok=True)
    summary_file = Path(summary_path) if summary_path is not None else output_root / "batch_summary.json"
    dashboard_path = output_root / "batch_summary.md"

    report_records = _index_existing_reports(output_root)
    state_records = _load_refresh_state_records(output_root)
    total = len(rows)
    rebuilt: list[dict[str, str | int | float | None]] = []

    for index, row in enumerate(rows, start=1):
        base_record = {key: value for key, value in row.record.items()}
        state_record = state_records.get((row.query, row.source_column))
        result: dict[str, str | int | float | None] = {
            **base_record,
            "query": row.query,
            "source_column": row.source_column,
            "batch_index": index,
            "batch_total": total,
            "status": state_record.get("status", "error") if state_record else "error",
            "resolved_entry_name": state_record.get("resolved_entry_name") if state_record else None,
            "json_report": state_record.get("json_report") if state_record else None,
            "markdown_report": state_record.get("markdown_report") if state_record else None,
            "duration_seconds": None,
            "error": state_record.get("error") if state_record else None,
        }
        report_record = _locate_existing_report(row=row, report_records=report_records, output_root=output_root)
        if report_record is not None:
            json_path = report_record["json_path"]
            markdown_path = json_path.with_suffix(".md")
            if markdown_path.exists():
                result["status"] = "ok"
                result["resolved_entry_name"] = report_record.get("entry_name")
                result["json_report"] = str(json_path)
                result["markdown_report"] = str(markdown_path)
                result["error"] = None
            else:
                result["error"] = f"Missing markdown report: {markdown_path.name}"
        else:
            if result["status"] not in {"pending", "running", "error"}:
                result["status"] = "error"
            if result["status"] == "error" and not result.get("error"):
                expected_stem = _expected_stem(row)
                result["error"] = f"Missing report files for {expected_stem}"
            if result["status"] in {"pending", "running"}:
                result["error"] = None
        rebuilt.append(result)

    summary_file.write_text(json.dumps(rebuilt, indent=2) + "\n", encoding="utf-8")
    dashboard_path.write_text(render_batch_summary_markdown(rebuilt), encoding="utf-8")
    return rebuilt


def _load_refresh_state_records(output_root: Path) -> dict[tuple[str, str], dict[str, Any]]:
    state_path = output_root / "accessible_portal_refresh_state.json"
    if not state_path.exists():
        return {}
    try:
        payload = json.loads(state_path.read_text(encoding="utf-8"))
    except Exception:
        return {}
    records: dict[tuple[str, str], dict[str, Any]] = {}
    for item in payload.get("items", []):
        query = item.get("query")
        source_column = item.get("source_column")
        if not query or not source_column:
            continue
        records[(str(query), str(source_column))] = item
    return records


def write_batch_summary_markdown(
    summary_path: str | Path,
    output_path: str | Path | None = None,
) -> Path:
    summary_file = Path(summary_path)
    results = json.loads(summary_file.read_text(encoding="utf-8"))
    destination = Path(output_path) if output_path is not None else summary_file.with_suffix(".md")
    destination.write_text(render_batch_summary_markdown(results), encoding="utf-8")
    return destination


def render_batch_summary_markdown(results: list[dict[str, str | int | float | None]]) -> str:
    total = max(
        (
            int(item.get("batch_total"))
            for item in results
            if item.get("batch_total") is not None
        ),
        default=len(results),
    )
    completed = len(results)
    ok = sum(1 for item in results if item.get("status") == "ok")
    skipped = sum(1 for item in results if item.get("status") == "skipped_existing")
    errors = sum(1 for item in results if item.get("status") == "error")
    avg_duration = _safe_average(
        [
            float(item["duration_seconds"])
            for item in results
            if item.get("duration_seconds") is not None
        ]
    )
    percent = (completed / total * 100.0) if total else 0.0

    lines = [
        "# Batch Summary",
        "",
        f"- Progress: {completed}/{total} ({percent:.1f}%)",
        f"- Successful: {ok}",
        f"- Skipped existing: {skipped}",
        f"- Errors: {errors}",
        f"- Average duration: {avg_duration:.2f}s" if avg_duration is not None else "- Average duration: n/a",
        "",
    ]

    if results:
        lines.extend(
            [
                "## Recent Results",
                "",
                "| # | Query | Status | Resolved Entry | Seconds | Markdown | Error |",
                "| --- | --- | --- | --- | --- | --- | --- |",
            ]
        )
        for item in results[-15:]:
            markdown_report = item.get("markdown_report")
            markdown_label = Path(str(markdown_report)).name if markdown_report else ""
            lines.append(
                "| {index} | `{query}` | {status} | `{resolved}` | {seconds} | `{markdown}` | {error} |".format(
                    index=item.get("batch_index", ""),
                    query=item.get("query", ""),
                    status=item.get("status", ""),
                    resolved=item.get("resolved_entry_name") or "",
                    seconds=item.get("duration_seconds") if item.get("duration_seconds") is not None else "",
                    markdown=markdown_label,
                    error=_truncate(str(item.get("error") or ""), 120),
                )
            )
        lines.append("")

    error_rows = [item for item in results if item.get("status") == "error"]
    if error_rows:
        lines.extend(
            [
                "## Recent Errors",
                "",
                "| # | Query | Error |",
                "| --- | --- | --- |",
            ]
        )
        for item in error_rows[-10:]:
            lines.append(
                f"| {item.get('batch_index', '')} | `{item.get('query', '')}` | {_truncate(str(item.get('error') or ''), 180)} |"
            )
        lines.append("")

    return "\n".join(lines) + "\n"


def _select_query(record: dict[str, str]) -> tuple[str | None, str | None]:
    gene = (record.get("gene") or "").strip()
    if gene:
        return gene, "gene"

    entry_name = (record.get("uniprot_name") or "").strip()
    if entry_name:
        return entry_name, "uniprot_name"

    accession = (record.get("uniprot_id") or "").strip()
    if accession:
        return accession, "uniprot_id"
    return None, None


def _expected_stem(row: BatchRow) -> str:
    if row.source_column == "uniprot_name":
        return f"{row.query.lower()}_report"
    if row.source_column == "gene":
        return f"{row.query.lower()}_human_report"
    return f"{row.query.lower()}_report"


def _safe_average(values: list[float]) -> float | None:
    if not values:
        return None
    return sum(values) / len(values)


def _truncate(text: str, limit: int) -> str:
    if len(text) <= limit:
        return text
    return text[: limit - 3] + "..."


def _needs_retry(
    result: dict[str, str | int | float | None],
    *,
    output_root: Path,
    require_alphafold: bool,
) -> bool:
    if result.get("status") == "error":
        return True
    if not require_alphafold:
        return False
    json_report = result.get("json_report")
    if not json_report:
        return True
    report_path = Path(str(json_report))
    if not report_path.exists():
        candidate = output_root / report_path.name
        if candidate.exists():
            report_path = candidate
        else:
            return True
    try:
        payload = json.loads(report_path.read_text(encoding="utf-8"))
    except Exception:
        return True
    accession = ((payload.get("target") or {}).get("accession") or "").strip()
    if not accession:
        return True
    return _find_alphafold_pdb(accession, output_root=output_root) is None


def _find_alphafold_pdb(accession: str, *, output_root: Path) -> Path | None:
    """Locate an AlphaFold PDB in the current batch or snapshot layout."""

    file_name = f"{accession}.pdb"
    bases: list[Path] = []
    for base in (output_root, output_root.parent, output_root.parent.parent, Path.cwd()):
        if base not in bases:
            bases.append(base)
    candidates: list[Path] = []
    for base in bases:
        candidates.extend(
            [
                base / "alphafold" / file_name,
                base / "data" / "alphafold" / file_name,
                base / ".agdesign2" / "data" / "alphafold" / file_name,
                base / "alphafold3" / file_name,
                base / "data" / "alphafold3" / file_name,
                base / ".agdesign2" / "data" / "alphafold3" / file_name,
            ]
        )
    return next((candidate.resolve() for candidate in candidates if candidate.exists()), None)


def _index_existing_reports(output_root: Path) -> dict[str, dict[str, Any]]:
    by_stem: dict[str, dict[str, Any]] = {}
    by_entry_name: dict[str, dict[str, Any]] = {}
    by_accession: dict[str, dict[str, Any]] = {}
    for json_path in sorted(output_root.glob("*_report.json")):
        try:
            payload = json.loads(json_path.read_text(encoding="utf-8"))
        except Exception:
            continue
        target = payload.get("target") or {}
        entry_name = str(target.get("entry_name") or "").strip()
        accession = str(target.get("accession") or "").strip()
        record = {
            "json_path": json_path,
            "entry_name": entry_name,
            "accession": accession,
        }
        by_stem[json_path.stem] = record
        if entry_name:
            by_entry_name[entry_name] = record
        if accession:
            by_accession[accession] = record
    return {
        "by_stem": by_stem,
        "by_entry_name": by_entry_name,
        "by_accession": by_accession,
    }


def _locate_existing_report(
    *,
    row: BatchRow,
    report_records: dict[str, dict[str, Any]],
    output_root: Path,
) -> dict[str, Any] | None:
    expected_stem = _expected_stem(row)
    direct_path = output_root / f"{expected_stem}.json"
    if direct_path.exists():
        return report_records["by_stem"].get(expected_stem) or {"json_path": direct_path, "entry_name": None, "accession": None}

    candidates = [
        report_records["by_stem"].get(expected_stem),
        report_records["by_entry_name"].get(row.query),
        report_records["by_entry_name"].get((row.record.get("uniprot_name") or "").strip()),
        report_records["by_accession"].get((row.record.get("uniprot_id") or "").strip()),
    ]
    for candidate in candidates:
        if candidate is not None:
            return candidate
    return None
