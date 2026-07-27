from __future__ import annotations

import json
import os
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .batch import load_batch_rows, rebuild_batch_summary_from_tsv
from .config import AnalysisConfig
from .pipeline import AntigenAnalyzer
from .portal import build_portal
from .target_sets import DEFAULT_ACCESSIBLE_BUCKETS, build_accessible_target_tsv


def _utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


@dataclass(slots=True)
class RefreshJobResult:
    state_path: Path
    target_tsv_path: Path
    summary_path: Path
    processed_this_run: int
    total_targets: int
    completed_targets: int
    portal_built: bool


def run_accessible_portal_refresh_job(
    *,
    source_csv: str | Path,
    output_tsv: str | Path,
    output_dir: str | Path,
    state_path: str | Path,
    allowed_buckets: tuple[str, ...] = DEFAULT_ACCESSIBLE_BUCKETS,
    portal_interval: int = 25,
    summary_interval: int = 10,
    max_targets: int | None = None,
    max_attempts: int | None = None,
    verbose: bool = False,
    enable_complex_portal: bool = False,
) -> RefreshJobResult:
    source_path = Path(source_csv)
    target_tsv_path = Path(output_tsv)
    output_root = Path(output_dir)
    output_root.mkdir(parents=True, exist_ok=True)
    state_file = Path(state_path)
    state_file.parent.mkdir(parents=True, exist_ok=True)
    summary_path = output_root / "batch_summary.json"
    lock_path = state_file.with_suffix(state_file.suffix + ".lock")

    build_accessible_target_tsv(
        source_csv=source_path,
        output_tsv=target_tsv_path,
        allowed_buckets=allowed_buckets,
    )
    rows = load_batch_rows(target_tsv_path)
    state = _load_or_initialize_state(
        state_file=state_file,
        source_path=source_path,
        target_tsv_path=target_tsv_path,
        rows=rows,
        output_root=output_root,
        summary_path=summary_path,
        portal_interval=portal_interval,
        summary_interval=summary_interval,
        allowed_buckets=allowed_buckets,
    )

    lock_fd = _acquire_lock(lock_path)
    try:
        analyzer = AntigenAnalyzer(
            config=AnalysisConfig(
                verbose_progress=verbose,
                enable_complex_portal=enable_complex_portal,
            )
        )
        processed_this_run = 0
        portal_built = False
        for item in state["items"]:
            if item.get("status") == "ok":
                continue
            # Optionally stop retrying targets that have already failed
            # max_attempts times so a permanently-failing target does not
            # re-consume the per-run budget on every invocation. Off by default
            # (None) to preserve existing retry-forever behavior.
            if (
                max_attempts is not None
                and item.get("status") == "error"
                and int(item.get("attempts") or 0) >= max_attempts
            ):
                continue
            if max_targets is not None and processed_this_run >= max_targets:
                break
            processed_this_run += 1
            item["status"] = "running"
            item["started_at"] = _utc_now()
            item["attempts"] = int(item.get("attempts") or 0) + 1
            _write_state(state_file, state)
            query = str(item.get("query") or "")
            if verbose:
                print(f"[agdesign2-refresh] {item.get('batch_index')}/{state['total_targets']} {query}", flush=True)
            try:
                report = analyzer.analyze_target(query, output_dir=output_root)
                stem = f"{report.target.entry_name.lower()}_report"
                item["status"] = "ok"
                item["resolved_entry_name"] = report.target.entry_name
                item["json_report"] = str(output_root / f"{stem}.json")
                item["markdown_report"] = str(output_root / f"{stem}.md")
                item["error"] = None
            except Exception as exc:
                item["status"] = "error"
                item["error"] = str(exc)
            item["finished_at"] = _utc_now()
            state["updated_at"] = item["finished_at"]
            _write_state(state_file, state)

            if _should_sync_summary(state, processed_this_run):
                rebuild_batch_summary_from_tsv(
                    tsv_path=target_tsv_path,
                    output_dir=output_root,
                    summary_path=summary_path,
                )
                state["last_summary_sync_at"] = _utc_now()
                _write_state(state_file, state)
            if _should_sync_portal(state, processed_this_run):
                rebuild_batch_summary_from_tsv(
                    tsv_path=target_tsv_path,
                    output_dir=output_root,
                    summary_path=summary_path,
                )
                build_portal(summary_path, refresh_reports=False, refresh_stale_reports=False, verbose=verbose)
                state["last_portal_sync_at"] = _utc_now()
                _write_state(state_file, state)
                portal_built = True

        rebuild_batch_summary_from_tsv(
            tsv_path=target_tsv_path,
            output_dir=output_root,
            summary_path=summary_path,
        )
        if _all_targets_complete(state):
            build_portal(summary_path, refresh_reports=False, refresh_stale_reports=False, verbose=verbose)
            portal_built = True
            state["completed_at"] = _utc_now()
        state["updated_at"] = _utc_now()
        _write_state(state_file, state)
        completed_targets = sum(1 for item in state["items"] if item.get("status") == "ok")
        return RefreshJobResult(
            state_path=state_file,
            target_tsv_path=target_tsv_path,
            summary_path=summary_path,
            processed_this_run=processed_this_run,
            total_targets=state["total_targets"],
            completed_targets=completed_targets,
            portal_built=portal_built,
        )
    finally:
        _release_lock(lock_fd, lock_path)


def _load_or_initialize_state(
    *,
    state_file: Path,
    source_path: Path,
    target_tsv_path: Path,
    rows: list[Any],
    output_root: Path,
    summary_path: Path,
    portal_interval: int,
    summary_interval: int,
    allowed_buckets: tuple[str, ...],
) -> dict[str, Any]:
    if state_file.exists():
        state = json.loads(state_file.read_text(encoding="utf-8"))
        existing_signature = [(item.get("query"), item.get("source_column")) for item in state.get("items", [])]
        new_signature = [(row.query, row.source_column) for row in rows]
        if existing_signature == new_signature:
            _reconcile_state_with_existing_reports(state, rows=rows, output_root=output_root)
            state["updated_at"] = _utc_now()
            _write_state(state_file, state)
            return state
        existing_by_key = {
            (item.get("query"), item.get("source_column")): item
            for item in state.get("items", [])
        }
        existing_reports = _index_existing_report_stems(output_root)
        migrated_items: list[dict[str, Any]] = []
        total_targets = len(rows)
        for index, row in enumerate(rows, start=1):
            key = (row.query, row.source_column)
            existing = existing_by_key.get(key)
            if existing is not None:
                item = {
                    **{name: value for name, value in row.record.items()},
                    **existing,
                    "query": row.query,
                    "source_column": row.source_column,
                    "batch_index": index,
                    "batch_total": total_targets,
                }
                if item.get("status") == "running":
                    item["status"] = "pending"
                report_paths = _matching_existing_report(row, existing_reports)
                if report_paths is not None:
                    item["status"] = "ok"
                    item["resolved_entry_name"] = report_paths["entry_name"]
                    item["json_report"] = str(report_paths["json"])
                    item["markdown_report"] = str(report_paths["markdown"])
                    item["error"] = None
                migrated_items.append(item)
                continue
            report_paths = _matching_existing_report(row, existing_reports)
            if report_paths is not None:
                migrated_items.append(
                    {
                        **{key_name: value for key_name, value in row.record.items()},
                        "query": row.query,
                        "source_column": row.source_column,
                        "batch_index": index,
                        "batch_total": total_targets,
                        "status": "ok",
                        "resolved_entry_name": report_paths["entry_name"],
                        "json_report": str(report_paths["json"]),
                        "markdown_report": str(report_paths["markdown"]),
                        "error": None,
                        "attempts": 0,
                        "started_at": None,
                        "finished_at": None,
                    }
                )
                continue
            migrated_items.append(
                {
                    **{key_name: value for key_name, value in row.record.items()},
                    "query": row.query,
                    "source_column": row.source_column,
                    "batch_index": index,
                    "batch_total": total_targets,
                    "status": "pending",
                    "resolved_entry_name": None,
                    "json_report": None,
                    "markdown_report": None,
                    "error": None,
                    "attempts": 0,
                    "started_at": None,
                    "finished_at": None,
                }
            )
        state.update(
            {
                "version": 2,
                "source_csv": str(source_path),
                "target_tsv": str(target_tsv_path),
                "output_dir": str(output_root),
                "summary_path": str(summary_path),
                "portal_interval": int(max(portal_interval, 1)),
                "summary_interval": int(max(summary_interval, 1)),
                "allowed_buckets": list(allowed_buckets),
                "updated_at": _utc_now(),
                "completed_at": None,
                "total_targets": total_targets,
                "items": migrated_items,
            }
        )
        _write_state(state_file, state)
        return state

    items: list[dict[str, Any]] = []
    total_targets = len(rows)
    for index, row in enumerate(rows, start=1):
        items.append(
            {
                **{key: value for key, value in row.record.items()},
                "query": row.query,
                "source_column": row.source_column,
                "batch_index": index,
                "batch_total": total_targets,
                "status": "pending",
                "resolved_entry_name": None,
                "json_report": None,
                "markdown_report": None,
                "error": None,
                "attempts": 0,
                "started_at": None,
                "finished_at": None,
            }
        )
    state = {
        "version": 1,
        "source_csv": str(source_path),
        "target_tsv": str(target_tsv_path),
        "output_dir": str(output_root),
        "summary_path": str(summary_path),
        "portal_interval": int(max(portal_interval, 1)),
        "summary_interval": int(max(summary_interval, 1)),
        "allowed_buckets": list(allowed_buckets),
        "created_at": _utc_now(),
        "updated_at": _utc_now(),
        "completed_at": None,
        "last_summary_sync_at": None,
        "last_portal_sync_at": None,
        "total_targets": total_targets,
        "items": items,
    }
    _write_state(state_file, state)
    return state


def _reconcile_state_with_existing_reports(state: dict[str, Any], *, rows: list[Any], output_root: Path) -> None:
    existing_reports = _index_existing_report_stems(output_root)
    state_items = {
        (item.get("query"), item.get("source_column")): item
        for item in state.get("items", [])
    }
    for row in rows:
        item = state_items.get((row.query, row.source_column))
        if item is None:
            continue
        report_paths = _matching_existing_report(row, existing_reports)
        if report_paths is None:
            # Reset both stale "ok" items whose report vanished and "running"
            # items left over from a crashed/killed prior run; otherwise a
            # leftover "running" keeps _all_targets_complete False and the final
            # portal build is skipped. (The migration path already does this.)
            if item.get("status") in {"ok", "running"}:
                item["status"] = "pending"
                item["resolved_entry_name"] = None
                item["json_report"] = None
                item["markdown_report"] = None
                item["error"] = None
            continue
        item["status"] = "ok"
        item["resolved_entry_name"] = report_paths["entry_name"]
        item["json_report"] = str(report_paths["json"])
        item["markdown_report"] = str(report_paths["markdown"])
        item["error"] = None


def _index_existing_report_stems(output_root: Path) -> dict[str, dict[str, Any]]:
    records: dict[str, dict[str, Any]] = {}
    for json_path in output_root.glob("*_report.json"):
        markdown_path = json_path.with_suffix(".md")
        if not markdown_path.exists():
            continue
        stem = json_path.name.removesuffix("_report.json")
        accession = None
        entry_name = stem.upper()
        gene_symbol = None
        try:
            payload = json.loads(json_path.read_text(encoding="utf-8"))
            target = payload.get("target") or {}
            accession = target.get("accession")
            entry_name = target.get("entry_name") or entry_name
            gene_symbol = target.get("gene_symbol")
        except Exception:
            pass
        records[stem] = {
            "accession": accession,
            "entry_name": entry_name,
            "gene_symbol": gene_symbol,
            "json": json_path,
            "markdown": markdown_path,
        }
    return records


def _matching_existing_report(row: Any, records: dict[str, dict[str, Any]]) -> dict[str, Any] | None:
    expected_accession = str(row.record.get("uniprot_id") or row.record.get("accession") or "").strip()
    expected_gene = str(row.record.get("gene_symbol") or row.record.get("gene") or "").strip().upper()
    candidates = [
        row.record.get("uniprot_name"),
        row.record.get("gene_symbol"),
        row.record.get("gene"),
        row.query,
    ]
    for candidate in candidates:
        text = str(candidate or "").strip()
        if not text:
            continue
        stem = text.lower()
        if stem.endswith("_human"):
            key = stem
        else:
            key = f"{stem}_human"
        record = records.get(key)
        if record is None:
            continue
        record_accession = str(record.get("accession") or "").strip()
        record_gene = str(record.get("gene_symbol") or "").strip().upper()
        if expected_accession and record_accession and expected_accession != record_accession:
            continue
        if expected_gene and record_gene and expected_gene != record_gene:
            continue
        return record
    return None


def _should_sync_summary(state: dict[str, Any], processed_this_run: int) -> bool:
    return processed_this_run > 0 and processed_this_run % int(state.get("summary_interval") or 1) == 0


def _should_sync_portal(state: dict[str, Any], processed_this_run: int) -> bool:
    return processed_this_run > 0 and processed_this_run % int(state.get("portal_interval") or 1) == 0


def _all_targets_complete(state: dict[str, Any]) -> bool:
    return all(item.get("status") in {"ok", "error"} for item in state.get("items", []))


def _write_state(state_file: Path, state: dict[str, Any]) -> None:
    state_file.write_text(json.dumps(state, indent=2) + "\n", encoding="utf-8")


def _load_lock_metadata(lock_path: Path) -> dict[str, Any] | None:
    if not lock_path.exists():
        return None
    try:
        return json.loads(lock_path.read_text(encoding="utf-8"))
    except Exception:
        return None


def _acquire_lock(lock_path: Path) -> int:
    metadata = _load_lock_metadata(lock_path)
    if metadata:
        existing_pid = int(metadata.get("pid") or 0)
        if existing_pid > 0:
            try:
                os.kill(existing_pid, 0)
            except OSError:
                lock_path.unlink(missing_ok=True)
            else:
                raise RuntimeError(f"Refresh job already running with PID {existing_pid}")
    fd = os.open(lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    with os.fdopen(fd, "w", encoding="utf-8") as handle:
        json.dump({"pid": os.getpid(), "created_at": _utc_now()}, handle)
    return os.open(lock_path, os.O_RDONLY)


def _release_lock(lock_fd: int, lock_path: Path) -> None:
    try:
        os.close(lock_fd)
    except OSError:
        pass
    lock_path.unlink(missing_ok=True)
