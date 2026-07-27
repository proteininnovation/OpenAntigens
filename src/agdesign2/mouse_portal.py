from __future__ import annotations

import csv
import json
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .af3 import copy_af3_catalog_artifacts
from .batch import run_batch_from_tsv
from .cache import FileCache
from .config import AnalysisConfig
from .http import HttpClient
from .module_refresh import refresh_report_modules
from .pipeline import AntigenAnalyzer
from .portal import build_portal
from .sequence_utils import global_align, slice_sequence
from .uniprot import UniProtClient
from .uniprot_prefetch import prefetch_targets_from_tsv


@dataclass(slots=True)
class MousePortalBuildResult:
    summary_path: Path
    output_dir: Path
    portal_index_path: Path | None
    included_count: int
    skipped_count: int
    skipped_manifest_path: Path


# Bump when the mouse analysis pipeline changes in a way that makes previously
# written reports stale. The build records this in reports/.mouse_pipeline_version;
# on a later run a mismatch (or a missing sentinel, e.g. reports left by the old
# annotation-transfer code) disables resume so every target is re-analyzed rather
# than silently reused.
MOUSE_PIPELINE_VERSION = "2026-07-no-default-af3"


def build_mouse_portal_from_human_orthologs(
    summary_path: str | Path,
    *,
    output_dir: str | Path,
    build_site: bool = True,
    bundle_vendor_assets: bool = False,
    fetch_mouse_alphafold: bool = True,
    af3_catalog: str | Path | None = None,
    jobs: int = 1,
    resume: bool = True,
    verbose: bool = False,
) -> MousePortalBuildResult:
    human_summary_path = Path(summary_path)
    human_batch_dir = human_summary_path.parent
    mouse_output_dir = Path(output_dir)
    mouse_output_dir.mkdir(parents=True, exist_ok=True)
    reports_dir = mouse_output_dir / "reports"
    reports_dir.mkdir(parents=True, exist_ok=True)
    mouse_summary_path = mouse_output_dir / "batch_summary.json"
    skipped_manifest_path = mouse_output_dir / "mouse_skipped_orthologs.json"
    mouse_targets_tsv = mouse_output_dir / "mouse_targets.tsv"

    ortholog_table_path = human_batch_dir / "ortholog_reference_table.tsv"
    ortholog_rows = _load_ortholog_rows(ortholog_table_path)
    if ortholog_table_path.exists():
        shutil.copy2(ortholog_table_path, mouse_output_dir / ortholog_table_path.name)
    human_summary = json.loads(human_summary_path.read_text(encoding="utf-8"))
    skipped: list[dict[str, str]] = []
    mouse_rows: list[dict[str, str]] = []
    provenance_by_entry: dict[str, dict[str, Any]] = {}
    human_to_mouse_entry: dict[str, str] = {}
    seen_mouse_targets: set[str] = set()

    uniprot_client = _mouse_uniprot_client(mouse_output_dir)
    for index, item in enumerate(human_summary, start=1):
        report_path = _resolve_path(item.get("json_report"), base=human_batch_dir)
        if report_path is None:
            skipped.append(_skip_row(item, "missing_human_report"))
            continue
        try:
            human_report = json.loads(report_path.read_text(encoding="utf-8"))
        except Exception as exc:
            skipped.append(_skip_row(item, f"invalid_human_report: {exc}"))
            continue

        mouse_target = _mouse_target_for_human_report(
            human_report=human_report,
            ortholog_rows=ortholog_rows,
            uniprot_client=uniprot_client,
        )
        if mouse_target is None:
            skipped.append(_skip_row(item, "missing_available_mouse_ortholog"))
            continue
        mouse_entry = str(mouse_target.get("entry_name") or mouse_target.get("accession") or "").strip()
        if not mouse_entry:
            skipped.append(_skip_row(item, "missing_mouse_identifier"))
            continue
        mouse_accession = str(mouse_target.get("accession") or "").strip()
        mouse_key = (mouse_accession or mouse_entry).upper()
        if mouse_key in seen_mouse_targets:
            skipped.append(_skip_row(item, f"duplicate_mouse_target:{mouse_entry}"))
            continue
        seen_mouse_targets.add(mouse_key)
        human_target = human_report.get("target") or {}
        human_entry = str(human_target.get("entry_name") or "").strip()
        provenance_by_entry[mouse_key] = {
            "human_source_target": {
                "accession": human_target.get("accession"),
                "entry_name": human_target.get("entry_name"),
                "gene_symbol": human_target.get("gene_symbol"),
                "protein_name": human_target.get("protein_name"),
            },
            "human_family_context": human_report.get("family_context"),
            "human_full_length_family_context": human_report.get("full_length_family_context"),
            "mouse_sequence_provenance": _mouse_sequence_provenance(mouse_target),
        }
        if human_entry:
            human_to_mouse_entry[human_entry.upper()] = mouse_accession or mouse_entry
        mouse_rows.append(
            {
                "uniprot_id": str(mouse_target.get("accession") or ""),
                "uniprot_name": mouse_entry,
                "gene_symbol": str(mouse_target.get("gene_symbol") or ""),
                "topology_bucket": item.get("topology_bucket") or _topology_bucket(human_report),
                "primary_topology_class": item.get("primary_topology_class") or _topology_class(human_report),
                "mouse_ortholog_accession": str(mouse_target.get("accession") or ""),
                "mouse_ortholog_entry_name": mouse_entry,
                "mouse_ortholog_gene_symbol": str(mouse_target.get("gene_symbol") or ""),
                "mouse_refseq_accession": str(mouse_target.get("refseq_accession") or ""),
                "mouse_sequence_source": str(mouse_target.get("sequence_source") or ""),
                "human_source_entry_name": (human_report.get("target") or {}).get("entry_name") or "",
                "human_source_accession": (human_report.get("target") or {}).get("accession") or "",
                "human_source_gene_symbol": (human_report.get("target") or {}).get("gene_symbol") or "",
            }
        )
        if verbose:
            print(f"[openantigens-mouse] included {mouse_entry}", flush=True)

    _write_mouse_target_tsv(mouse_targets_tsv, mouse_rows)
    af3_catalog_path = _resolve_af3_catalog(af3_catalog)
    snapshot_data_dir = _snapshot_data_dir(human_batch_dir)
    mouse_config = AnalysisConfig(
        target_species="mouse",
        fetch_alphafold=fetch_mouse_alphafold,
        enable_local_af3_fallback=af3_catalog_path is not None,
        cache_dir=mouse_output_dir / ".agdesign2" / "cache",
        data_dir=mouse_output_dir / ".agdesign2" / "data",
        blast_db_dir=snapshot_data_dir / "blastdb",
        ortholog_fasta_dir=snapshot_data_dir / "proteomes",
        ortholog_blast_db_dir=snapshot_data_dir / "ortholog_blastdb",
        prefer_precomputed_references=False,
        # Per-target BLAST during analyze is slow (one blastp per target, per DB).
        # Cross-reactivity is instead computed once in bulk below via search_many,
        # the same fast path the human snapshot uses.
        run_cross_reactivity=False,
        generate_assets=False,
        enable_complex_portal=False,
    )
    if mouse_rows:
        prefetch_targets_from_tsv(
            tsv_path=mouse_targets_tsv,
            cache_dir=mouse_config.cache_dir,
        )

    if af3_catalog_path is not None and mouse_rows:
        copy_af3_catalog_artifacts(af3_catalog_path, data_dir=mouse_config.data_dir, target_tsv=mouse_targets_tsv)
    if mouse_rows:
        _require_blast_tools_for_mouse_cross_reactivity()
    # Resume only when existing reports were produced by the current pipeline
    # version. A missing/mismatched sentinel means the reports are stale (e.g.
    # left by an older pipeline), so re-analyze everything rather than silently
    # reusing them.
    effective_resume, recorded_version = _mouse_resume_decision(
        reports_dir,
        resume,
        fetch_mouse_alphafold=fetch_mouse_alphafold,
    )
    if resume and not effective_resume and any(reports_dir.glob("*.json")):
        print(
            "[openantigens-mouse] existing reports are from pipeline version "
            f"{recorded_version or 'unknown (pre-versioning)'}; current is "
            f"{MOUSE_PIPELINE_VERSION} -- re-analyzing all targets instead of resuming.",
            flush=True,
        )
    mouse_analyzer = AntigenAnalyzer(config=mouse_config)
    batch_results = run_batch_from_tsv(
        analyzer=mouse_analyzer,
        tsv_path=mouse_targets_tsv,
        output_dir=reports_dir,
        resume=effective_resume,
        jobs=jobs,
    )
    _record_mouse_pipeline_version(
        reports_dir,
        fetch_mouse_alphafold=fetch_mouse_alphafold,
    )
    mouse_summary = _finalize_mouse_batch_results(
        batch_results,
        provenance_by_entry=provenance_by_entry,
        human_to_mouse_entry=human_to_mouse_entry,
        reports_dir=reports_dir,
    )
    mouse_summary_path.write_text(json.dumps(mouse_summary, indent=2) + "\n", encoding="utf-8")
    skipped_manifest_path.write_text(json.dumps(skipped, indent=2) + "\n", encoding="utf-8")
    # Compute cross-reactivity once in bulk -- a single multi-FASTA blastp per
    # database via search_many -- instead of one blastp per target. Re-reads the
    # reports written above and fills in the hits before the portal is rendered.
    if mouse_rows:
        refresh_report_modules(
            mouse_summary_path, modules=["cross-reactivity"], analyzer=mouse_analyzer, verbose=verbose
        )
    _warn_if_cross_reactivity_missing(
        mouse_summary, reports_dir=reports_dir, run_cross_reactivity=bool(mouse_rows)
    )

    portal_index_path = None
    if build_site:
        portal_index_path = build_portal(
            mouse_summary_path,
            output_dir=mouse_output_dir / "portal",
            refresh_reports=False,
            refresh_stale_reports=False,
            bundle_vendor_assets=bundle_vendor_assets,
            portal_title="OpenAntigens Mouse",
            portal_subtitle="Mouse orthologs for OpenAntigens human targets.",
            portal_intro="OpenAntigens Mouse applies the same structure-aware construct pipeline to resolved mouse orthologs of human targets, using compatible mouse AlphaFold structure evidence when available, running mouse BLAST searches, and recording the originating human target on every report.",
            theme="mouse",
            include_disease_context=False,
            sibling_link=("OpenAntigens Human", "../index.html"),
            detail_sibling_link=("OpenAntigens Human", "../../index.html"),
        )

    return MousePortalBuildResult(
        summary_path=mouse_summary_path,
        output_dir=mouse_output_dir,
        portal_index_path=portal_index_path,
        included_count=len(mouse_summary),
        skipped_count=len(skipped),
        skipped_manifest_path=skipped_manifest_path,
    )


def _mouse_resume_decision(
    reports_dir: Path,
    requested_resume: bool,
    *,
    fetch_mouse_alphafold: bool = True,
) -> tuple[bool, str | None]:
    """Resume the mouse batch only if existing reports carry the current pipeline
    version sentinel. Returns (effective_resume, recorded_version)."""
    version_file = reports_dir / ".mouse_pipeline_version"
    recorded = version_file.read_text(encoding="utf-8").strip() if version_file.exists() else None
    expected = _mouse_pipeline_stamp(fetch_mouse_alphafold)
    return (requested_resume and recorded == expected), recorded


def _record_mouse_pipeline_version(
    reports_dir: Path,
    *,
    fetch_mouse_alphafold: bool = True,
) -> None:
    reports_dir.mkdir(parents=True, exist_ok=True)
    (reports_dir / ".mouse_pipeline_version").write_text(
        _mouse_pipeline_stamp(fetch_mouse_alphafold) + "\n",
        encoding="utf-8",
    )


def _mouse_pipeline_stamp(fetch_mouse_alphafold: bool) -> str:
    return f"{MOUSE_PIPELINE_VERSION};fetch_alphafold={int(fetch_mouse_alphafold)}"


def _warn_if_cross_reactivity_missing(
    summary: list[dict[str, Any]],
    *,
    reports_dir: Path,
    run_cross_reactivity: bool,
    sample: int = 300,
) -> None:
    """Loudly warn if cross-reactivity is empty across a sample of reports despite
    being enabled -- catches BLAST/database misconfiguration that would otherwise
    ship a silently empty mouse portal."""
    if not run_cross_reactivity:
        return
    ok = [row for row in summary if row.get("status") == "ok"]
    if len(ok) < 20:
        return
    checked = empty = 0
    for row in ok[:sample]:
        path = _resolve_path(row.get("json_report"), base=reports_dir)
        if path is None or not path.exists():
            continue
        try:
            report = json.loads(path.read_text(encoding="utf-8"))
        except Exception:
            continue
        checked += 1
        if not (report.get("cross_reactivity_hits") or report.get("full_length_cross_reactivity_hits")):
            empty += 1
    if checked and empty / checked > 0.9:
        print(
            f"[openantigens-mouse] WARNING: {empty}/{checked} sampled reports have NO cross-reactivity "
            "hits despite run_cross_reactivity=True. Check that BLAST+ and the human/mouse/cynomolgus "
            "BLAST databases are available to the mouse build.",
            flush=True,
        )


def _mouse_target_for_human_report(
    *,
    human_report: dict[str, Any],
    ortholog_rows: dict[str, dict[str, str]],
    uniprot_client: UniProtClient,
) -> dict[str, Any] | None:
    human_target = human_report.get("target") or {}
    human_entry = str(human_target.get("entry_name") or "").strip()
    ortholog = ortholog_rows.get(human_entry.upper())
    mouse_homology = _first_available_mouse_homology(human_report)
    return _mouse_target_from_sources(
        human_report=human_report,
        human_target=human_target,
        ortholog=ortholog,
        mouse_homology=mouse_homology,
        uniprot_client=uniprot_client,
    )


def _write_mouse_target_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    fieldnames = [
        "uniprot_id",
        "uniprot_name",
        "gene_symbol",
        "topology_bucket",
        "primary_topology_class",
        "mouse_ortholog_accession",
        "mouse_ortholog_entry_name",
        "mouse_ortholog_gene_symbol",
        "mouse_refseq_accession",
        "mouse_sequence_source",
        "human_source_entry_name",
        "human_source_accession",
        "human_source_gene_symbol",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _resolve_af3_catalog(value: str | Path | None) -> Path | None:
    if value is not None:
        candidate = Path(value)
        if candidate.exists():
            return candidate
    return None


def _require_blast_tools_for_mouse_cross_reactivity() -> None:
    missing = [tool for tool in ("blastp", "makeblastdb") if shutil.which(tool) is None]
    if missing:
        raise RuntimeError(
            "Mouse portal cross-reactivity requires BLAST+ on PATH; missing: "
            + ", ".join(missing)
            + ". Run `agdesign2 setup-blast` or install BLAST+ before building the mouse portal."
        )


def _snapshot_data_dir(human_batch_dir: Path) -> Path:
    snapshot_dir = human_batch_dir.parent.parent
    data_dir = snapshot_dir / "data"
    if data_dir.exists():
        return data_dir
    return human_batch_dir / ".agdesign2" / "data"


def _finalize_mouse_batch_results(
    batch_results: list[dict[str, str | int | float | None]],
    *,
    provenance_by_entry: dict[str, dict[str, Any]],
    human_to_mouse_entry: dict[str, str],
    reports_dir: Path,
) -> list[dict[str, str | int | float | None]]:
    finalized: list[dict[str, str | int | float | None]] = []
    reports_by_entry: dict[str, dict[str, Any]] = {}
    for result in batch_results:
        if result.get("status") not in {"ok", "skipped_existing"}:
            continue
        report_path = _resolve_path(result.get("json_report"), base=reports_dir)
        if report_path is None:
            continue
        report = json.loads(report_path.read_text(encoding="utf-8"))
        target = report.get("target") or {}
        mouse_entry = str(target.get("entry_name") or result.get("resolved_entry_name") or result.get("query") or "").strip()
        if mouse_entry:
            reports_by_entry[mouse_entry.upper()] = report
        mouse_accession = str(target.get("accession") or "").strip()
        if mouse_accession:
            reports_by_entry[mouse_accession.upper()] = report

    for result in batch_results:
        if result.get("status") not in {"ok", "skipped_existing"}:
            continue
        report_path = _resolve_path(result.get("json_report"), base=reports_dir)
        if report_path is None:
            continue
        report = json.loads(report_path.read_text(encoding="utf-8"))
        target = report.get("target") or {}
        mouse_entry = str(target.get("entry_name") or result.get("resolved_entry_name") or result.get("query") or "").strip()
        mouse_accession = str(target.get("accession") or "").strip()
        provenance = provenance_by_entry.get((mouse_accession or mouse_entry).upper(), {})
        metadata = report.get("metadata") if isinstance(report.get("metadata"), dict) else {}
        report["metadata"] = {
            **metadata,
            "portal_species": "mouse",
            "source_database": "mouse_openantigens_structure_analysis",
            "human_source_target": provenance.get("human_source_target") or {},
            "construct_identity_label": "Identity to mouse source",
            "mouse_sequence_provenance": provenance.get("mouse_sequence_provenance") or {},
        }
        report["family_context"] = _mouse_family_context_from_human_context(
            provenance.get("human_family_context"),
            reports_by_entry=reports_by_entry,
            human_to_mouse_entry=human_to_mouse_entry,
            target_entry=mouse_entry,
            matrix_scope="ectodomain",
        )
        report["full_length_family_context"] = _mouse_family_context_from_human_context(
            provenance.get("human_full_length_family_context") or provenance.get("human_family_context"),
            reports_by_entry=reports_by_entry,
            human_to_mouse_entry=human_to_mouse_entry,
            target_entry=mouse_entry,
            matrix_scope="full_length",
        )
        _append_note(
            report,
            "info",
            "This mouse report was generated by running the OpenAntigens structure-aware construct pipeline on the mouse ortholog target.",
            "OpenAntigens Mouse",
        )
        source = report["metadata"]["human_source_target"]
        if source:
            _append_note(
                report,
                "info",
                f"Source human OpenAntigens target: {source.get('entry_name') or 'n/a'} {source.get('accession') or ''}.",
                "OpenAntigens Mouse",
            )
        report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        finalized_row = dict(result)
        finalized_row["status"] = "ok"
        finalized_row["json_report"] = str(report_path)
        markdown_report = result.get("markdown_report")
        if markdown_report:
            finalized_row["markdown_report"] = str(_resolve_path(markdown_report, base=reports_dir) or markdown_report)
        finalized.append(finalized_row)
    return finalized


def _mouse_family_context_from_human_context(
    human_context: Any,
    *,
    reports_by_entry: dict[str, dict[str, Any]],
    human_to_mouse_entry: dict[str, str],
    target_entry: str,
    matrix_scope: str,
) -> dict[str, Any] | None:
    if not isinstance(human_context, dict):
        return None
    members = human_context.get("members") or []
    if not isinstance(members, list):
        return None
    resolved: list[dict[str, Any]] = []
    seen: set[str] = set()
    for human_member in members:
        if not isinstance(human_member, dict):
            continue
        human_entry = str(human_member.get("entry_name") or "").strip().upper()
        mouse_entry = human_to_mouse_entry.get(human_entry)
        if not mouse_entry or mouse_entry.upper() in seen:
            continue
        mouse_report = reports_by_entry.get(mouse_entry.upper())
        if mouse_report is None:
            continue
        mouse_target = mouse_report.get("target") or {}
        sequence = str(mouse_target.get("sequence") or "")
        if matrix_scope == "ectodomain":
            ectodomain = mouse_report.get("ectodomain") or {}
            start = _coerce_int(ectodomain.get("start"))
            end = _coerce_int(ectodomain.get("end"))
            if start is None or end is None:
                continue
            matrix_sequence = slice_sequence(sequence, start, end)
        else:
            start = None
            end = None
            matrix_sequence = sequence
        if not matrix_sequence:
            continue
        seen.add(mouse_entry.upper())
        resolved.append(
            {
                "gene_symbol": str(mouse_target.get("gene_symbol") or ""),
                "gene_name": str(mouse_target.get("protein_name") or mouse_target.get("gene_symbol") or ""),
                "accession": str(mouse_target.get("accession") or ""),
                "entry_name": str(mouse_target.get("entry_name") or mouse_entry),
                "sequence_length": len(sequence) if sequence else None,
                "ectodomain_start": start,
                "ectodomain_end": end,
                "ectodomain_length": len(matrix_sequence) if matrix_scope == "ectodomain" else None,
                "notes": [
                    "Mouse ortholog substituted for the corresponding human family member; only mouse proteins with human orthologs in this OpenAntigens Mouse build are included."
                ],
                "_matrix_sequence": matrix_sequence,
            }
        )
    if len(resolved) < 2:
        return None
    target_upper = target_entry.upper()
    resolved.sort(key=lambda item: 0 if str(item.get("entry_name") or "").upper() == target_upper else 1)
    labels = [str(item.get("entry_name") or item.get("accession") or item.get("gene_symbol") or "") for item in resolved]
    matrix, coverage = _identity_and_coverage_matrices([str(item["_matrix_sequence"]) for item in resolved])
    for item in resolved:
        item.pop("_matrix_sequence", None)
    family_names = [str(name) for name in (human_context.get("family_names") or []) if str(name).strip()]
    if not family_names and human_context.get("family_name"):
        family_names = [str(human_context.get("family_name"))]
    source = "Mouse ortholog family matrix recalculated from human OpenAntigens family membership"
    if matrix_scope == "full_length":
        source += "; full-length matrix recomputed locally"
    return {
        "gene_symbol": str((reports_by_entry.get(target_upper, {}).get("target") or {}).get("gene_symbol") or target_entry),
        "source": source,
        "family_names": family_names,
        "members": resolved,
        "identity_matrix_labels": labels,
        "identity_matrix": matrix,
        "coverage_matrix_labels": labels,
        "coverage_matrix": coverage,
        "metadata": {
            **(human_context.get("metadata") if isinstance(human_context.get("metadata"), dict) else {}),
            "matrix_scope": matrix_scope,
            "identity_denominator": f"row mouse {matrix_scope.replace('_', '-')} sequence length",
            "coverage_denominator": f"row mouse {matrix_scope.replace('_', '-')} sequence length",
            "mouse_ortholog_matrix_caveat": (
                "Only mouse proteins whose human OpenAntigens family members have resolved mouse orthologs in this mouse build are included."
            ),
        },
    }


def _identity_and_coverage_matrices(sequences: list[str]) -> tuple[list[list[float | None]], list[list[float | None]]]:
    identity_matrix: list[list[float | None]] = []
    coverage_matrix: list[list[float | None]] = []
    for row_sequence in sequences:
        identity_row: list[float | None] = []
        coverage_row: list[float | None] = []
        row_length = len(row_sequence)
        for column_sequence in sequences:
            if not row_sequence or not column_sequence:
                identity_row.append(None)
                coverage_row.append(None)
                continue
            alignment = global_align(row_sequence, column_sequence)
            identity_row.append(round(100.0 * alignment.matches / row_length, 2) if row_length else None)
            coverage_row.append(round(100.0 * alignment.aligned_positions / row_length, 2) if row_length else None)
        identity_matrix.append(identity_row)
        coverage_matrix.append(coverage_row)
    return identity_matrix, coverage_matrix


def _coerce_int(value: Any) -> int | None:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _mouse_uniprot_client(output_dir: Path) -> UniProtClient:
    return UniProtClient(HttpClient(FileCache(output_dir / ".agdesign2" / "cache")))


def _mouse_target_from_sources(
    *,
    human_report: dict[str, Any],
    human_target: dict[str, Any],
    ortholog: dict[str, str] | None,
    mouse_homology: dict[str, Any] | None,
    uniprot_client: UniProtClient,
) -> dict[str, Any] | None:
    homolog_accession = str((mouse_homology or {}).get("accession") or "").strip()
    homolog_entry = str((mouse_homology or {}).get("entry_name") or "").strip()
    refseq_accession = str((ortholog or {}).get("mouse_refseq_accession") or "").strip()
    refseq_sequence = str((ortholog or {}).get("mouse_refseq_sequence") or "").strip()
    accession = homolog_accession or refseq_accession
    entry_name = homolog_entry or accession
    sequence = refseq_sequence
    sequence_source = "ortholog_reference_table.mouse_refseq_sequence" if refseq_sequence else ""
    if not accession and not sequence:
        return None
    gene_symbol = _mouse_gene_symbol(human_target, ortholog)
    uniprot_target = _resolve_mouse_uniprot_target(gene_symbol, uniprot_client=uniprot_client)
    if uniprot_target is not None:
        accession = uniprot_target["accession"]
        entry_name = uniprot_target["entry_name"]
        sequence = uniprot_target["sequence"]
        sequence_source = "uniprot_reviewed_mouse_gene_symbol"
    protein_name = _mouse_protein_name(
        human_report=human_report,
        human_target=human_target,
        accession=accession,
        entry_name=entry_name,
        gene_symbol=gene_symbol,
    )
    if uniprot_target is not None and uniprot_target.get("protein_name"):
        protein_name = uniprot_target["protein_name"]
    return {
        "accession": accession,
        "entry_name": entry_name,
        "gene_symbol": gene_symbol,
        "protein_name": protein_name,
        "organism": "Mus musculus",
        "taxon_id": 10090,
        "sequence": sequence,
        "identifier_source": (
            "uniprot_reviewed_mouse_gene_symbol"
            if uniprot_target is not None
            else "human_report_mouse_homolog"
            if homolog_accession or homolog_entry
            else "ortholog_reference_table.mouse_refseq_accession"
        ),
        "sequence_source": sequence_source,
        "refseq_accession": refseq_accession,
        "uniprot_resolution_note": (
            "Reviewed mouse UniProt record resolved from mouse ortholog gene symbol."
            if uniprot_target is not None
            else "No reviewed mouse UniProt record was resolved from the mouse ortholog gene symbol; using explicit RefSeq ortholog sequence."
        ),
    }


def _resolve_mouse_uniprot_target(gene_symbol: str, *, uniprot_client: UniProtClient) -> dict[str, str] | None:
    symbol = str(gene_symbol or "").strip()
    if not symbol:
        return None
    try:
        resolution = uniprot_client.resolve_target(symbol, target_species_taxon=10090)
    except Exception:
        return None
    target = resolution.target
    if not target.accession or not target.entry_name or not str(target.entry_name).endswith("_MOUSE"):
        return None
    if not target.sequence:
        return None
    return {
        "accession": target.accession,
        "entry_name": target.entry_name,
        "sequence": target.sequence,
        "protein_name": target.protein_name,
    }


def _mouse_protein_name(
    *,
    human_report: dict[str, Any],
    human_target: dict[str, Any],
    accession: str,
    entry_name: str,
    gene_symbol: str,
) -> str:
    for hit in (human_report.get("full_length_cross_reactivity_hits") or []) + (human_report.get("cross_reactivity_hits") or []):
        subject_id = str(hit.get("subject_id") or "")
        description = str(hit.get("description") or "")
        if not _hit_matches_mouse_identifier(subject_id, description, accession=accession, entry_name=entry_name):
            continue
        name = _protein_name_from_blast_description(description)
        if name:
            return name
    human_name = str(human_target.get("protein_name") or "").strip()
    if human_name:
        return human_name
    return f"{gene_symbol or entry_name or accession} mouse ortholog"


def _hit_matches_mouse_identifier(description_id: str, description: str, *, accession: str, entry_name: str) -> bool:
    tokens = {token for token in (accession, entry_name) if token}
    haystack = f"{description_id} {description}"
    return any(token in haystack for token in tokens)


def _protein_name_from_blast_description(description: str) -> str:
    text = description.strip()
    if not text:
        return ""
    if "|" in text:
        parts = text.split("|", 2)
        if len(parts) == 3:
            text = parts[2].split(" ", 1)[1] if " " in parts[2] else ""
    for marker in (" OS=", " OX=", " GN=", " PE=", " SV="):
        if marker in text:
            text = text.split(marker, 1)[0]
            break
    return text.strip()


def _mouse_sequence_provenance(mouse_target: dict[str, Any]) -> dict[str, str]:
    return {
        "identifier_source": str(mouse_target.get("identifier_source") or ""),
        "sequence_source": str(mouse_target.get("sequence_source") or ""),
        "refseq_accession": str(mouse_target.get("refseq_accession") or ""),
        "alphafold_accession": str(mouse_target.get("alphafold_accession") or ""),
        "alphafold_status": str(mouse_target.get("alphafold_status") or ""),
    }


def _first_available_mouse_homology(report: dict[str, Any]) -> dict[str, Any] | None:
    for section in ("ectodomain_homology", "species_name_matches"):
        for item in report.get(section) or []:
            if str(item.get("species") or "").lower() == "mouse" and item.get("available"):
                return item
    for construct in report.get("construct_details") or []:
        for item in construct.get("homologs") or []:
            if str(item.get("species") or "").lower() == "mouse" and item.get("available"):
                return item
    return None


def _load_ortholog_rows(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return {
        str(row.get("input_uniprot_name") or "").strip().upper(): row
        for row in rows
        if str(row.get("input_uniprot_name") or "").strip()
    }


def _resolve_path(value: Any, *, base: Path) -> Path | None:
    if not value:
        return None
    path = Path(str(value))
    candidates = [path] if path.is_absolute() else [base / path, path]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def _skip_row(item: dict[str, Any], reason: str) -> dict[str, str]:
    return {
        "query": str(item.get("query") or ""),
        "resolved_entry_name": str(item.get("resolved_entry_name") or item.get("entry_name") or ""),
        "reason": reason,
    }


def _mouse_gene_symbol(human_target: dict[str, Any], ortholog: dict[str, str] | None) -> str:
    if ortholog is not None:
        value = str(ortholog.get("mouse_gene_symbol") or "").strip()
        if value:
            return value
    return str(human_target.get("gene_symbol") or "").strip().capitalize()


def _topology_bucket(report: dict[str, Any]) -> str:
    topology = report.get("topology") or {}
    topology_class = str(topology.get("topology_class") or "").lower()
    if "secret" in topology_class:
        return "Secreted"
    if "gpi" in topology_class:
        return "GPI"
    if "single" in topology_class:
        return "Single-pass"
    if "multi" in topology_class:
        return "Multipass"
    return ""


def _topology_class(report: dict[str, Any]) -> str:
    topology = report.get("topology") or {}
    return str(topology.get("topology_class") or "")


def _append_note(report: dict[str, Any], severity: str, message: str, source: str) -> None:
    notes = report.setdefault("notes", [])
    notes.append({"severity": severity, "message": message, "source": source})
