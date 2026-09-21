from __future__ import annotations

import csv
import gzip
import json
import os
import random
import re
import shutil
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from time import perf_counter

from .af3 import copy_af3_catalog_artifacts
from .assets import generate_assets_for_summary
from .batch import run_batch_from_tsv
from .config import AnalysisConfig
from .family_alignments import build_family_alignments_from_tsv
from .module_refresh import refresh_report_modules
from .mouse_portal import build_mouse_portal_from_human_orthologs
from .open_targets import build_open_targets_associations_from_bulk_downloads
from .ortholog_table import build_ortholog_table_from_tsv, load_ortholog_rows
from .paralogs import build_paralog_reference_from_tsv
from .pipeline import AntigenAnalyzer
from .portal import _write_portal_htaccess, build_portal
from .release_validation import validate_release_data
from .target_sets import DEFAULT_ACCESSIBLE_BUCKETS, build_accessible_target_tsv
from .uniprot_prefetch import prefetch_reviewed_taxa, prefetch_targets_from_tsv


@dataclass(slots=True)
class DeploymentPrecomputeResult:
    steps_planned: list[str]
    steps_completed: list[str]
    summary_path: Path
    output_dir: Path


@dataclass(slots=True)
class FreshSnapshotResult:
    snapshot_dir: Path
    target_tsv: Path
    output_dir: Path
    summary_path: Path
    portal_dir: Path
    public_site_dir: Path
    benchmark_path: Path
    manifest_path: Path
    timings: list[dict[str, object]]
    selected_smoke_tsv: Path | None = None
    full_smoke_manifest_path: Path | None = None
    mouse_output_dir: Path | None = None
    mouse_public_site_dir: Path | None = None


FRESH_SNAPSHOT_STEPS = (
    "build-target-tsv",
    "prefetch-uniprot",
    "setup-blast",
    "build-ortholog-table",
    "build-family-alignments",
    "build-paralog-reference",
    "analyze-targets",
    "generate-assets",
    "build-open-targets",
    "refresh-cross-reactivity",
    "build-portal",
    "package-public-site",
)


def build_fresh_snapshot(
    *,
    source_csv: str | Path = "data/inputs/accessible_human_topology_refined_accessible_only.csv",
    snapshot_root: str | Path = "snapshots",
    snapshot_name: str | None = None,
    allowed_buckets: tuple[str, ...] = DEFAULT_ACCESSIBLE_BUCKETS,
    limit: int | None = None,
    sample_size: int | None = None,
    sample_seed: int = 20260530,
    sample_mode: str | None = None,
    include_mouse: bool = True,
    fetch_mouse_alphafold: bool = True,
    force: bool = False,
    resume: bool = False,
    reanalyze_reports: bool = False,
    taxa: tuple[int, ...] = (9606, 10090),
    uniprot_prefetch_mode: str = "targets",
    uniprot_page_size: int = 500,
    open_targets_page_size: int = 100,
    open_targets_max_pages: int = 5,
    jobs: int = 0,
    asset_jobs: int = 0,
    ortholog_jobs: int = 1,
    af3_catalog: str | Path | None = None,
    render_structure_images: bool = False,
    render_quality_plots: bool = True,
    dry_run: bool = False,
    verbose: bool = False,
) -> FreshSnapshotResult:
    created_at = datetime.now(timezone.utc).replace(microsecond=0)
    name = snapshot_name or f"openantigen-fresh-{created_at.date().isoformat()}"
    snapshot_dir = Path(snapshot_root) / name
    target_tsv = snapshot_dir / "inputs" / "accessible_targets.tsv"
    selected_smoke_tsv = snapshot_dir / "inputs" / "selected_smoke_targets.tsv"
    output_dir = snapshot_dir / "outputs" / "surfy_batch"
    summary_path = output_dir / "batch_summary.json"
    portal_dir = output_dir / "portal"
    public_site_dir = snapshot_dir / "public_site"
    mouse_output_dir = snapshot_dir / "mouse_data"
    mouse_public_site_dir = public_site_dir / "mouse"
    benchmark_path = snapshot_dir / "snapshot_benchmark.json"
    manifest_path = snapshot_dir / "snapshot_manifest.json"
    full_smoke_manifest_path = snapshot_dir / "full_smoke_manifest.json"
    af3_catalog_path = _resolve_default_af3_catalog(af3_catalog)
    timings: list[dict[str, object]] = []
    sample_mode_text = str(sample_mode or "").strip() or None
    if sample_size is not None and sample_size <= 0:
        raise ValueError("--sample-size must be greater than zero.")
    if sample_size is not None and limit is not None:
        raise ValueError("--limit cannot be combined with --sample-size; choose one target selection mode.")
    if sample_mode_text not in {None, "stratified-random"}:
        raise ValueError("--sample-mode must be stratified-random when provided.")
    if sample_mode_text is not None and sample_size is None:
        raise ValueError("--sample-mode requires --sample-size.")
    if sample_size is not None and sample_mode_text is None:
        sample_mode_text = "stratified-random"

    if dry_run:
        return FreshSnapshotResult(
            snapshot_dir=snapshot_dir,
            target_tsv=target_tsv,
            output_dir=output_dir,
            summary_path=summary_path,
            portal_dir=portal_dir,
            public_site_dir=public_site_dir,
            benchmark_path=benchmark_path,
            manifest_path=manifest_path,
            selected_smoke_tsv=selected_smoke_tsv if sample_size is not None else None,
            full_smoke_manifest_path=full_smoke_manifest_path if sample_size is not None else None,
            mouse_output_dir=mouse_output_dir if include_mouse else None,
            mouse_public_site_dir=mouse_public_site_dir if include_mouse else None,
            timings=[
                {"step": step, "status": "planned"}
                for step in (
                    FRESH_SNAPSHOT_STEPS[:1]
                    + (("prepare-af3-catalog-artifacts",) if af3_catalog_path is not None else ())
                    + FRESH_SNAPSHOT_STEPS[1:]
                    + (("build-mouse-portal",) if include_mouse else ())
                    + (("validate-release-data",) if sample_size is None and limit is None else ())
                )
            ],
        )

    if force and resume:
        raise ValueError("--force and --resume cannot be used together.")

    if snapshot_dir.exists():
        if force:
            shutil.rmtree(snapshot_dir)
        elif not resume:
            raise FileExistsError(f"Snapshot already exists: {snapshot_dir}. Use --force, --resume, or a new --snapshot-name.")
    snapshot_dir.mkdir(parents=True, exist_ok=True)
    (snapshot_dir / "inputs").mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    if resume and benchmark_path.exists():
        try:
            existing_progress = json.loads(benchmark_path.read_text(encoding="utf-8"))
            existing_timings = existing_progress.get("timings") if isinstance(existing_progress, dict) else None
            if isinstance(existing_timings, list):
                timings.extend(row for row in existing_timings if isinstance(row, dict))
        except Exception:
            pass

    cache_dir = snapshot_dir / "cache"
    data_dir = snapshot_dir / "data"
    config = AnalysisConfig(
        cache_dir=cache_dir,
        data_dir=data_dir,
        enable_local_af3_fallback=af3_catalog_path is not None,
        blast_db_dir=data_dir / "blastdb",
        ortholog_fasta_dir=data_dir / "proteomes",
        ortholog_blast_db_dir=data_dir / "ortholog_blastdb",
        precomputed_ortholog_table_path=output_dir / "ortholog_reference_table.tsv",
        precomputed_family_alignment_index_path=output_dir / "family_alignments" / "family_alignment_index.tsv",
        precomputed_paralog_index_path=output_dir / "paralogs" / "paralog_index.tsv",
        generate_assets=False,
        render_structure_images=render_structure_images,
        render_quality_plots=render_quality_plots,
        run_cross_reactivity=False,
        enable_complex_portal=True,
        verbose_progress=verbose,
    )
    analyzer = AntigenAnalyzer(config=config)

    def run_step(name: str, callback, *, skip_if=None) -> object:
        _log(verbose, name)
        started = perf_counter()
        if resume and skip_if is not None and skip_if():
            row = {
                "step": name,
                "status": "skipped_existing",
                "duration_seconds": 0,
                "finished_utc": datetime.now(timezone.utc).isoformat(),
            }
            timings.append(row)
            _write_snapshot_progress(benchmark_path, timings)
            return None
        row: dict[str, object] = {"step": name, "status": "started", "started_utc": datetime.now(timezone.utc).isoformat()}
        timings.append(row)
        _write_snapshot_progress(benchmark_path, timings)
        try:
            result = callback()
        except Exception as exc:
            row.update(
                {
                    "status": "error",
                    "duration_seconds": round(perf_counter() - started, 3),
                    "error": str(exc),
                    "finished_utc": datetime.now(timezone.utc).isoformat(),
                }
            )
            _write_snapshot_progress(benchmark_path, timings)
            raise
        row.update(
            {
                "status": "ok",
                "duration_seconds": round(perf_counter() - started, 3),
                "finished_utc": datetime.now(timezone.utc).isoformat(),
            }
        )
        _write_snapshot_progress(benchmark_path, timings)
        return result

    def make_target_tsv() -> list[dict[str, str]]:
        rows = build_accessible_target_tsv(
            source_csv=source_csv,
            output_tsv=target_tsv,
            allowed_buckets=allowed_buckets,
        )
        if sample_size is not None:
            rows, sample_manifest = _select_stratified_smoke_rows(
                rows,
                sample_size=sample_size,
                sample_seed=sample_seed,
            )
            _write_target_tsv(target_tsv, rows)
            _write_target_tsv(selected_smoke_tsv, rows)
            _write_full_smoke_manifest(
                full_smoke_manifest_path,
                sample_manifest=sample_manifest,
                selected_tsv=selected_smoke_tsv,
                include_mouse=include_mouse,
                fetch_mouse_alphafold=fetch_mouse_alphafold,
            )
        elif limit is not None:
            rows = rows[:limit]
            _write_target_tsv(target_tsv, rows)
        return rows

    run_step("build-target-tsv", make_target_tsv, skip_if=lambda: target_tsv.exists())
    if af3_catalog_path is not None:
        run_step(
            "prepare-af3-catalog-artifacts",
            lambda: copy_af3_catalog_artifacts(af3_catalog_path, data_dir=data_dir, target_tsv=target_tsv),
            skip_if=lambda: (data_dir / "alphafold3" / "af3_catalog_copy_manifest.json").exists(),
        )
    run_step(
        "prefetch-uniprot",
        lambda: (
            prefetch_targets_from_tsv(tsv_path=target_tsv, cache_dir=cache_dir, chunk_size=uniprot_page_size)
            if uniprot_prefetch_mode == "targets"
            else prefetch_reviewed_taxa(cache_dir=cache_dir, taxa=taxa, page_size=uniprot_page_size)
        ),
        # Re-enter on resume: completed pages are cached, while a partial cache
        # must not make an interrupted prefetch look complete.
        skip_if=None,
    )
    run_step(
        "setup-blast",
        lambda: (analyzer.blast_client.ensure_databases(), analyzer.blast_client.ensure_ortholog_databases()),
        skip_if=lambda: _has_any_file(config.blast_db_dir) and _has_any_file(config.ortholog_blast_db_dir),
    )
    run_step(
        "build-ortholog-table",
        lambda: build_ortholog_table_from_tsv(
            analyzer=analyzer,
            tsv_path=target_tsv,
            output_path=output_dir / "ortholog_reference_table.tsv",
            resume=resume,
            verbose=verbose,
            jobs=ortholog_jobs,
        ),
        skip_if=lambda: _ortholog_table_complete(output_dir / "ortholog_reference_table.tsv", target_tsv),
    )
    run_step(
        "build-family-alignments",
        lambda: build_family_alignments_from_tsv(
            analyzer=analyzer,
            tsv_path=target_tsv,
            output_dir=output_dir / "family_alignments",
            verbose=verbose,
        ),
        skip_if=lambda: (output_dir / "family_alignments" / "family_alignment_index.tsv").exists(),
    )
    run_step(
        "build-paralog-reference",
        lambda: build_paralog_reference_from_tsv(
            analyzer=analyzer,
            tsv_path=target_tsv,
            output_dir=output_dir / "paralogs",
            verbose=verbose,
        ),
        skip_if=lambda: (output_dir / "paralogs" / "paralog_index.tsv").exists(),
    )
    run_step(
        "analyze-targets",
        lambda: run_batch_from_tsv(
            analyzer=analyzer,
            tsv_path=target_tsv,
            output_dir=output_dir,
            # --reanalyze-reports forces a full re-derivation of every report so a
            # report-logic change lands while the cached upstream artifacts above
            # (UniProt/AlphaFold/BLAST, ortholog/family/paralog tables) are reused.
            resume=resume and not reanalyze_reports,
            jobs=jobs,
        ),
    )
    _validate_batch_completion(target_tsv, summary_path)
    if resume:
        run_step(
            "refresh-complex-portal",
            lambda: refresh_report_modules(summary_path, modules=["complex-portal"], verbose=verbose, analyzer=analyzer),
        )
    run_step(
        "generate-assets",
        lambda: generate_assets_for_summary(
            summary_path,
            config=config,
            jobs=asset_jobs or jobs,
            resume=resume,
            verbose=verbose,
        ),
        skip_if=lambda: (not render_structure_images and not render_quality_plots),
    )
    run_step(
        "build-open-targets",
        lambda: build_open_targets_associations_from_bulk_downloads(
            summary_path,
            output_dir=output_dir,
            data_dir=data_dir / "open_targets",
            verbose=verbose,
        ),
        skip_if=lambda: (output_dir / "open_targets_disease_associations.tsv").exists(),
    )
    run_step(
        "refresh-cross-reactivity",
        lambda: refresh_report_modules(summary_path, modules=["cross-reactivity"], verbose=verbose, analyzer=analyzer),
    )
    run_step(
        "build-portal",
        lambda: build_portal(
            summary_path,
            refresh_reports=False,
            refresh_stale_reports=False,
            verbose=verbose,
            enable_complex_portal=False,
            bundle_vendor_assets=True,
        ),
    )
    run_step("package-public-site", lambda: _package_public_site(portal_dir, public_site_dir))
    if include_mouse:
        run_step(
            "build-mouse-portal",
            lambda: _build_mouse_site_for_snapshot(
                summary_path=summary_path,
                mouse_output_dir=mouse_output_dir,
                mouse_public_site_dir=mouse_public_site_dir,
                resume=resume and not reanalyze_reports,
                fetch_mouse_alphafold=fetch_mouse_alphafold,
                af3_catalog=af3_catalog_path,
                jobs=jobs,
                verbose=verbose,
            ),
            skip_if=lambda: resume and not reanalyze_reports and (mouse_public_site_dir / "index.html").exists(),
        )
    if sample_size is None and limit is None:
        run_step(
            "validate-release-data",
            lambda: validate_release_data(public_site_dir, require_full_catalog=True),
        )

    _write_snapshot_manifest(
        manifest_path=manifest_path,
        snapshot_dir=snapshot_dir,
        source_csv=Path(source_csv),
        target_tsv=target_tsv,
        output_dir=output_dir,
        portal_dir=portal_dir,
        public_site_dir=public_site_dir,
        benchmark_path=benchmark_path,
        created_at=created_at,
        limit=limit,
        sample_size=sample_size,
        sample_seed=sample_seed,
        sample_mode=sample_mode_text,
        include_mouse=include_mouse,
        fetch_mouse_alphafold=fetch_mouse_alphafold,
        selected_smoke_tsv=selected_smoke_tsv if sample_size is not None else None,
        full_smoke_manifest_path=full_smoke_manifest_path if sample_size is not None else None,
        mouse_output_dir=mouse_output_dir if include_mouse else None,
        mouse_public_site_dir=mouse_public_site_dir if include_mouse else None,
        allowed_buckets=allowed_buckets,
        uniprot_prefetch_mode=uniprot_prefetch_mode,
        resume=resume,
        jobs=jobs,
        asset_jobs=asset_jobs or jobs,
        af3_catalog=af3_catalog_path,
        render_structure_images=render_structure_images,
        render_quality_plots=render_quality_plots,
        timings=timings,
    )
    return FreshSnapshotResult(
        snapshot_dir=snapshot_dir,
        target_tsv=target_tsv,
        output_dir=output_dir,
        summary_path=summary_path,
        portal_dir=portal_dir,
        public_site_dir=public_site_dir,
        benchmark_path=benchmark_path,
        manifest_path=manifest_path,
        selected_smoke_tsv=selected_smoke_tsv if sample_size is not None else None,
        full_smoke_manifest_path=full_smoke_manifest_path if sample_size is not None else None,
        mouse_output_dir=mouse_output_dir if include_mouse else None,
        mouse_public_site_dir=mouse_public_site_dir if include_mouse else None,
        timings=timings,
    )


def prepare_deployment_data(
    *,
    target_tsv: str | Path = "accessible_secreted_gpi_singlepass.tsv",
    summary_path: str | Path = "outputs/surfy_batch/batch_summary.json",
    output_dir: str | Path = "outputs/surfy_batch",
    cache_dir: str | Path = ".agdesign2/cache",
    data_dir: str | Path = ".agdesign2/data",
    taxa: tuple[int, ...] = (9606, 10090),
    uniprot_page_size: int = 500,
    open_targets_page_size: int = 100,
    open_targets_max_pages: int = 5,
    skip_uniprot: bool = False,
    skip_blast_setup: bool = False,
    skip_orthologs: bool = False,
    skip_family_alignments: bool = False,
    skip_paralogs: bool = False,
    skip_open_targets: bool = False,
    skip_cross_reactivity: bool = False,
    skip_portal: bool = False,
    dry_run: bool = False,
    verbose: bool = False,
) -> DeploymentPrecomputeResult:
    target_tsv_path = Path(target_tsv)
    summary_file = Path(summary_path)
    output_root = Path(output_dir)
    data_root = Path(data_dir)
    config = AnalysisConfig(
        cache_dir=Path(cache_dir),
        data_dir=data_root,
        blast_db_dir=data_root / "blastdb",
        ortholog_fasta_dir=data_root / "proteomes",
        ortholog_blast_db_dir=data_root / "ortholog_blastdb",
        precomputed_ortholog_table_path=output_root / "ortholog_reference_table.tsv",
        precomputed_family_alignment_index_path=output_root / "family_alignments" / "family_alignment_index.tsv",
        precomputed_paralog_index_path=output_root / "paralogs" / "paralog_index.tsv",
        generate_assets=False,
        render_structure_images=False,
        render_quality_plots=False,
        enable_complex_portal=False,
        verbose_progress=verbose,
    )
    steps = _planned_steps(
        skip_uniprot=skip_uniprot,
        skip_blast_setup=skip_blast_setup,
        skip_orthologs=skip_orthologs,
        skip_family_alignments=skip_family_alignments,
        skip_paralogs=skip_paralogs,
        skip_open_targets=skip_open_targets,
        skip_cross_reactivity=skip_cross_reactivity,
        skip_portal=skip_portal,
    )
    completed: list[str] = []
    if dry_run:
        return DeploymentPrecomputeResult(
            steps_planned=steps,
            steps_completed=completed,
            summary_path=summary_file,
            output_dir=output_root,
        )

    output_root.mkdir(parents=True, exist_ok=True)
    analyzer: AntigenAnalyzer | None = None

    if not skip_uniprot:
        _log(verbose, "bulk-prefetching reviewed UniProt entries")
        prefetch_reviewed_taxa(
            cache_dir=Path(cache_dir),
            taxa=taxa,
            page_size=uniprot_page_size,
        )
        completed.append("prefetch-uniprot")

    if not skip_blast_setup:
        analyzer = analyzer or AntigenAnalyzer(config=config)
        _log(verbose, "preparing BLAST databases")
        analyzer.blast_client.ensure_databases()
        analyzer.blast_client.ensure_ortholog_databases()
        completed.append("setup-blast")

    if not skip_orthologs:
        analyzer = analyzer or AntigenAnalyzer(config=config)
        _log(verbose, "building ortholog reference table")
        build_ortholog_table_from_tsv(
            analyzer=analyzer,
            tsv_path=target_tsv_path,
            output_path=output_root / "ortholog_reference_table.tsv",
            resume=True,
            verbose=verbose,
        )
        completed.append("build-ortholog-table")

    if not skip_family_alignments:
        analyzer = analyzer or AntigenAnalyzer(config=config)
        _log(verbose, "building family alignments")
        build_family_alignments_from_tsv(
            analyzer=analyzer,
            tsv_path=target_tsv_path,
            output_dir=output_root / "family_alignments",
            verbose=verbose,
        )
        completed.append("build-family-alignments")

    if not skip_paralogs:
        analyzer = analyzer or AntigenAnalyzer(config=config)
        _log(verbose, "building paralog reference sets")
        build_paralog_reference_from_tsv(
            analyzer=analyzer,
            tsv_path=target_tsv_path,
            output_dir=output_root / "paralogs",
            verbose=verbose,
        )
        completed.append("build-paralog-reference")

    if not skip_open_targets:
        _log(verbose, "building Open Targets disease associations")
        build_open_targets_associations_from_bulk_downloads(
            summary_file,
            output_dir=output_root,
            data_dir=analyzer.config.data_dir / "open_targets",
            verbose=verbose,
        )
        completed.append("build-open-targets")

    if not skip_cross_reactivity:
        _log(verbose, "refreshing cross-reactivity with bulk BLAST")
        refresh_report_modules(
            summary_file,
            modules=["cross-reactivity"],
            verbose=verbose,
        )
        completed.append("refresh-cross-reactivity")

    if not skip_portal:
        _log(verbose, "building static portal")
        build_portal(
            summary_file,
            refresh_reports=False,
            refresh_stale_reports=False,
            verbose=verbose,
            enable_complex_portal=False,
            bundle_vendor_assets=True,
        )
        completed.append("build-portal")

    return DeploymentPrecomputeResult(
        steps_planned=steps,
        steps_completed=completed,
        summary_path=summary_file,
        output_dir=output_root,
    )


def _planned_steps(
    *,
    skip_uniprot: bool,
    skip_blast_setup: bool,
    skip_orthologs: bool,
    skip_family_alignments: bool,
    skip_paralogs: bool,
    skip_open_targets: bool,
    skip_cross_reactivity: bool,
    skip_portal: bool,
) -> list[str]:
    return [
        name
        for name, skip in (
            ("prefetch-uniprot", skip_uniprot),
            ("setup-blast", skip_blast_setup),
            ("build-ortholog-table", skip_orthologs),
            ("build-family-alignments", skip_family_alignments),
            ("build-paralog-reference", skip_paralogs),
            ("build-open-targets", skip_open_targets),
            ("refresh-cross-reactivity", skip_cross_reactivity),
            ("build-portal", skip_portal),
        )
        if not skip
    ]


def _write_target_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    fieldnames = [
        "uniprot_id",
        "uniprot_name",
        "gene_symbol",
        "topology_bucket",
        "primary_topology_class",
        "protein_name",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _select_stratified_smoke_rows(
    rows: list[dict[str, str]],
    *,
    sample_size: int,
    sample_seed: int,
) -> tuple[list[dict[str, str]], dict[str, object]]:
    if sample_size <= 0:
        raise ValueError("sample_size must be greater than zero.")

    keyed_rows: list[dict[str, str]] = []
    seen: set[str] = set()
    for row in rows:
        key = _target_unique_key(row)
        if not key or key in seen:
            continue
        seen.add(key)
        keyed_rows.append(row)

    ordered_rows = sorted(
        keyed_rows,
        key=lambda row: (
            str(row.get("topology_bucket") or ""),
            str(row.get("uniprot_name") or ""),
            str(row.get("uniprot_id") or ""),
            str(row.get("gene_symbol") or ""),
        ),
    )
    if len(ordered_rows) <= sample_size:
        selected = ordered_rows
        return selected, _sample_manifest(
            selected,
            sample_size=sample_size,
            sample_seed=sample_seed,
            available_counts=_stratum_counts(ordered_rows),
            selected_by_stratum={},
            fill_count=len(selected),
        )

    rng = random.Random(sample_seed)
    strata = _smoke_strata(ordered_rows)
    available_strata = [(name, candidates) for name, candidates in strata.items() if candidates]
    selected: list[dict[str, str]] = []
    selected_keys: set[str] = set()
    selected_by_stratum: dict[str, int] = {}

    if available_strata:
        base_quota = sample_size // len(available_strata)
        remainder = sample_size % len(available_strata)
        for index, (name, candidates) in enumerate(available_strata):
            quota = base_quota + (1 if index < remainder else 0)
            pool = list(candidates)
            rng.shuffle(pool)
            before = len(selected)
            for row in pool:
                if len(selected) >= sample_size or len(selected) - before >= quota:
                    break
                key = _target_unique_key(row)
                if key in selected_keys:
                    continue
                selected.append(row)
                selected_keys.add(key)
            selected_by_stratum[name] = len(selected) - before

    fill_pool = list(ordered_rows)
    rng.shuffle(fill_pool)
    fill_count = 0
    for row in fill_pool:
        if len(selected) >= sample_size:
            break
        key = _target_unique_key(row)
        if key in selected_keys:
            continue
        selected.append(row)
        selected_keys.add(key)
        fill_count += 1

    return selected, _sample_manifest(
        selected,
        sample_size=sample_size,
        sample_seed=sample_seed,
        available_counts={name: len(candidates) for name, candidates in strata.items()},
        selected_by_stratum=selected_by_stratum,
        fill_count=fill_count,
    )


def _target_unique_key(row: dict[str, str]) -> str:
    return (
        str(row.get("uniprot_name") or "").strip().upper()
        or str(row.get("uniprot_id") or "").strip().upper()
        or str(row.get("gene_symbol") or "").strip().upper()
    )


def _smoke_strata(rows: list[dict[str, str]]) -> dict[str, list[dict[str, str]]]:
    return {
        "Secreted": [row for row in rows if _bucket(row) == "secreted"],
        "GPI": [row for row in rows if _bucket(row) == "gpi"],
        "Single-pass": [row for row in rows if _bucket(row) == "single-pass"],
        "Multipass": [row for row in rows if _bucket(row) == "multipass"],
        "GPCR-like": [row for row in rows if _row_text_contains(row, ("gpcr", "g protein-coupled"))],
        "Obligatory-complex": [row for row in rows if _row_text_contains(row, ("obligatory", "complex"))],
    }


def _bucket(row: dict[str, str]) -> str:
    return str(row.get("topology_bucket") or "").strip().lower()


def _row_text_contains(row: dict[str, str], needles: tuple[str, ...]) -> bool:
    text = " ".join(
        str(row.get(field) or "")
        for field in ("primary_topology_class", "protein_name", "gene_symbol", "uniprot_name", "topology_bucket")
    ).lower()
    return any(needle in text for needle in needles)


def _stratum_counts(rows: list[dict[str, str]]) -> dict[str, int]:
    return {name: len(candidates) for name, candidates in _smoke_strata(rows).items()}


def _sample_manifest(
    selected: list[dict[str, str]],
    *,
    sample_size: int,
    sample_seed: int,
    available_counts: dict[str, int],
    selected_by_stratum: dict[str, int],
    fill_count: int,
) -> dict[str, object]:
    return {
        "sample_mode": "stratified-random",
        "sample_size_requested": sample_size,
        "sample_seed": sample_seed,
        "selected_count": len(selected),
        "available_counts": available_counts,
        "selected_by_stratum": selected_by_stratum,
        "random_fill_count": fill_count,
        "selected_entries": [
            {
                "uniprot_id": row.get("uniprot_id") or "",
                "uniprot_name": row.get("uniprot_name") or "",
                "gene_symbol": row.get("gene_symbol") or "",
                "topology_bucket": row.get("topology_bucket") or "",
                "primary_topology_class": row.get("primary_topology_class") or "",
            }
            for row in selected
        ],
    }


def _write_full_smoke_manifest(
    path: Path,
    *,
    sample_manifest: dict[str, object],
    selected_tsv: Path,
    include_mouse: bool,
    fetch_mouse_alphafold: bool,
) -> None:
    manifest = {
        "database": "OpenAntigens",
        "snapshot_mode": "full-smoke",
        "selected_smoke_targets_tsv": str(selected_tsv),
        "sample": sample_manifest,
        "mouse": {
            "include_mouse": include_mouse,
            "fetch_mouse_alphafold": fetch_mouse_alphafold,
            "source": "human OpenAntigens reports with annotated mouse orthologs",
        },
    }
    path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def _build_mouse_site_for_snapshot(
    *,
    summary_path: Path,
    mouse_output_dir: Path,
    mouse_public_site_dir: Path,
    resume: bool,
    fetch_mouse_alphafold: bool,
    af3_catalog: Path | None,
    jobs: int,
    verbose: bool,
) -> Path:
    result = build_mouse_portal_from_human_orthologs(
        summary_path,
        output_dir=mouse_output_dir,
        resume=resume,
        build_site=True,
        bundle_vendor_assets=True,
        fetch_mouse_alphafold=fetch_mouse_alphafold,
        af3_catalog=af3_catalog,
        jobs=jobs,
        verbose=verbose,
    )
    source_portal = mouse_output_dir / "portal"
    if result.portal_index_path is None or not source_portal.exists():
        raise FileNotFoundError(f"Mouse portal build did not create a portal directory: {source_portal}")
    if mouse_public_site_dir.exists():
        shutil.rmtree(mouse_public_site_dir)
    shutil.copytree(
        source_portal,
        mouse_public_site_dir,
        ignore=lambda _directory, names: {
            name
            for name in names
            if name in {".DS_Store", "Thumbs.db", "__pycache__", "public_site"}
            or (name.startswith(".") and name not in {".well-known", ".htaccess"})
        },
    )
    return mouse_public_site_dir


def _write_snapshot_progress(path: Path, timings: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"timings": timings}, indent=2) + "\n", encoding="utf-8")


def _has_any_file(path: Path) -> bool:
    if not path.exists():
        return False
    return any(candidate.is_file() for candidate in path.rglob("*"))


def _ortholog_table_complete(output_tsv: Path, target_tsv: Path) -> bool:
    """True only when the ortholog reference table covers every target.

    The table is written incrementally (one row per target), so a mere
    existence check would treat a partially-built table as done and skip the
    step on resume — leaving downstream analysis with missing cross-species
    homolog data. Compare the written row count against the expected target
    count instead, so an interrupted build resumes the ortholog step and runs
    it to completion.
    """
    json_path = output_tsv.with_suffix(".json")
    if not json_path.exists():
        return False
    try:
        rows = json.loads(json_path.read_text(encoding="utf-8"))
        expected = len(load_ortholog_rows(target_tsv))
    except Exception:
        return False
    return isinstance(rows, list) and expected > 0 and len(rows) >= expected


def _package_public_site(portal_dir: Path, public_site_dir: Path) -> Path:
    if public_site_dir.exists():
        shutil.rmtree(public_site_dir)
    shutil.copytree(
        portal_dir,
        public_site_dir,
        ignore=lambda _directory, names: {
            name
            for name in names
            if name in {".DS_Store", "Thumbs.db", "__pycache__", "public_site"}
            or (name.startswith(".") and name not in {".well-known", ".htaccess"})
        },
    )
    for name in ("LICENSE", "DATA_LICENSE.md"):
        source = Path.cwd() / name
        if source.exists():
            shutil.copy2(source, public_site_dir / name)
    return public_site_dir


# --- public_site disk-size reduction for the disk-quota-limited GoDaddy host ---
#
# The packaged public_site/ is ~26 GB, over the host's ~25 GB quota. We store the
# large browser-fetched text assets gzip-compressed under their public filenames.
# The precompressed portal .htaccess supplies Content-Encoding directly, avoiding
# the mod_rewrite dependency that failed on the live GoDaddy host in August 2026.

_COMPRESSIBLE_SUFFIXES = (".js", ".css", ".svg")
_GZIP_LEVEL = 6  # matches the ratio measured on real assets (report .js ~7.5x)


def _is_gzip_file(path: Path) -> bool:
    with path.open("rb") as handle:
        return handle.read(2) == b"\x1f\x8b"


def _gzip_file_atomic(src: Path, level: int = _GZIP_LEVEL) -> int | None:
    """Replace ``src`` atomically with gzip bytes while retaining its filename.

    Returns the compressed size, or ``None`` when ``src`` is already gzip data.
    An interrupted write leaves the original intact.
    """
    if _is_gzip_file(src):
        return None
    tmp = src.with_name(src.name + ".gzip.tmp")
    try:
        with src.open("rb") as fin, gzip.open(tmp, "wb", compresslevel=level) as fout:
            shutil.copyfileobj(fin, fout, length=1024 * 1024)
        os.replace(tmp, src)
    finally:
        if tmp.exists():
            tmp.unlink()
    return src.stat().st_size


def _dir_size(directory: Path) -> int:
    return sum(p.stat().st_size for p in directory.rglob("*") if p.is_file())


def _strip_open_targets_json_button(downloads_html: Path) -> bool:
    """Remove the now-dead 'Open Targets JSON' button from downloads.html.

    Returns True if the file changed. The .json flat file is dropped (the sibling
    .tsv holds the same data), so its button would otherwise 404.
    """
    if not downloads_html.exists():
        return False
    text = downloads_html.read_text(encoding="utf-8")
    pattern = re.compile(
        r'[ \t]*<a class="button" href="open_targets_disease_associations\.json">'
        r"[^<]*</a>\n?"
    )
    updated = pattern.sub("", text)
    if updated == text:
        return False
    downloads_html.write_text(updated, encoding="utf-8")
    return True


def compress_public_site(
    public_site_dir: Path | str,
    *,
    drop_unreferenced_structures: bool = False,
    drop_redundant_downloads: bool = False,
    verbose: bool = False,
) -> dict[str, int]:
    """Shrink a packaged public_site/ to fit the GoDaddy disk quota. Idempotent.

    Gzips browser-fetched text assets (.js/.css/.svg) in place, retaining their
    public filenames. The compressed portal .htaccess supplies the encoding and
    MIME headers without relying on URL rewriting.

    ``drop_unreferenced_structures`` removes ``structures/`` and
    ``mouse/structures/``: the 3D viewer inlines PDB text into each report's JS,
    so those standalone copies are referenced by nothing (the source structures
    remain under the snapshot's outputs/ and are regenerable).

    ``drop_redundant_downloads`` removes
    ``open_targets_disease_associations.json`` (the sibling .tsv holds the same
    data and stays a plain, wget-friendly download) and strips its now-dead button
    from ``downloads.html``.

    Returns a summary of file/byte counts.
    """
    public_site_dir = Path(public_site_dir)
    if not (public_site_dir / "index.html").exists():
        raise FileNotFoundError(
            f"not a packaged public_site (no index.html): {public_site_dir}"
        )

    suffixes = set(_COMPRESSIBLE_SUFFIXES)

    # Materialize the target list before replacing files in place.
    targets = [
        path
        for path in public_site_dir.rglob("*")
        if path.suffix in suffixes
        and path.is_file()
        and not path.is_symlink()
    ]

    summary = {
        "files_compressed": 0,
        "bytes_before": 0,
        "bytes_after": 0,
        "structures_bytes_removed": 0,
        "downloads_bytes_removed": 0,
    }
    for path in targets:
        before = path.stat().st_size
        after = _gzip_file_atomic(path)
        if after is None:
            continue
        summary["files_compressed"] += 1
        summary["bytes_before"] += before
        summary["bytes_after"] += after
        if verbose and summary["files_compressed"] % 2000 == 0:
            print(f"[compress-public-site] gzipped {summary['files_compressed']} files")

    if drop_unreferenced_structures:
        for rel in ("structures", "mouse/structures"):
            directory = public_site_dir / rel
            if directory.is_dir():
                summary["structures_bytes_removed"] += _dir_size(directory)
                shutil.rmtree(directory)
                if verbose:
                    print(f"[compress-public-site] removed {rel}/")

    if drop_redundant_downloads:
        json_flat = public_site_dir / "open_targets_disease_associations.json"
        if json_flat.exists():
            summary["downloads_bytes_removed"] += json_flat.stat().st_size
            json_flat.unlink()
        if _strip_open_targets_json_button(public_site_dir / "downloads.html") and verbose:
            print("[compress-public-site] stripped Open Targets JSON download button")

    # Refresh every .htaccess so the precompressed-serving rules are present. A
    # snapshot packaged before these rules existed would otherwise 404 on .js,
    # and mod_rewrite rules in the root .htaccess do not inherit into mouse/,
    # which ships its own .htaccess — so update each one in place.
    htaccess_dirs = {public_site_dir}
    htaccess_dirs.update(path.parent for path in public_site_dir.rglob(".htaccess"))
    for directory in htaccess_dirs:
        _write_portal_htaccess(directory, precompressed=True)

    return summary


def _write_snapshot_manifest(
    *,
    manifest_path: Path,
    snapshot_dir: Path,
    source_csv: Path,
    target_tsv: Path,
    output_dir: Path,
    portal_dir: Path,
    public_site_dir: Path,
    benchmark_path: Path,
    created_at: datetime,
    limit: int | None,
    sample_size: int | None,
    sample_seed: int,
    sample_mode: str | None,
    include_mouse: bool,
    fetch_mouse_alphafold: bool,
    selected_smoke_tsv: Path | None,
    full_smoke_manifest_path: Path | None,
    mouse_output_dir: Path | None,
    mouse_public_site_dir: Path | None,
    allowed_buckets: tuple[str, ...],
    uniprot_prefetch_mode: str,
    resume: bool,
    jobs: int,
    asset_jobs: int,
    af3_catalog: Path | None,
    render_structure_images: bool,
    render_quality_plots: bool,
    timings: list[dict[str, object]],
) -> None:
    summary_path = output_dir / "batch_summary.json"
    target_count = _count_tsv_rows(target_tsv)
    summary_count = _count_json_list(summary_path)
    manifest = {
        "database": "OpenAntigens",
        "snapshot_created_utc": created_at.isoformat(),
        "snapshot_dir": str(snapshot_dir),
        "source_csv": str(source_csv),
        "target_tsv": str(target_tsv),
        "output_dir": str(output_dir),
        "portal_dir": str(portal_dir),
        "public_site_dir": str(public_site_dir),
        "benchmark_path": str(benchmark_path),
        "limit": limit,
        "sample": {
            "sample_size": sample_size,
            "sample_seed": sample_seed,
            "sample_mode": sample_mode,
            "selected_smoke_tsv": str(selected_smoke_tsv) if selected_smoke_tsv is not None else None,
            "full_smoke_manifest_path": str(full_smoke_manifest_path) if full_smoke_manifest_path is not None else None,
        },
        "mouse_portal": {
            "include_mouse": include_mouse,
            "fetch_mouse_alphafold": fetch_mouse_alphafold,
            "mouse_output_dir": str(mouse_output_dir) if mouse_output_dir is not None else None,
            "mouse_public_site_dir": str(mouse_public_site_dir) if mouse_public_site_dir is not None else None,
        },
        "allowed_buckets": list(allowed_buckets),
        "target_count": target_count,
        "summary_count": summary_count,
        "timings": timings,
        "freshness_policy": {
            "local_inputs_reused": ["surface/secreted protein universe CSV"],
            "local_caches_reused": [str(snapshot_dir / "cache"), str(snapshot_dir / "data")] if resume else [],
            "snapshot_local_cache": str(snapshot_dir / "cache"),
            "snapshot_local_data": str(snapshot_dir / "data"),
            "uniprot_prefetch_mode": uniprot_prefetch_mode,
            "resume": resume,
            "analyze_jobs": jobs,
            "asset_jobs": asset_jobs,
            "af3_catalog": str(af3_catalog) if af3_catalog is not None else None,
            "render_structure_images": render_structure_images,
            "render_quality_plots": render_quality_plots,
            "notes": "Retrieval and enrichment use snapshot-local cache/data directories. Resume reuses those artifacts and validates target completion before packaging.",
        },
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def _validate_batch_completion(target_tsv: Path, summary_path: Path) -> None:
    from collections import Counter
    from .batch import load_batch_rows
    from .portal import _resolve_existing_path

    expected = Counter((row.query, row.source_column) for row in load_batch_rows(target_tsv))
    rows = json.loads(summary_path.read_text())
    actual = Counter((row.get("query"), row.get("source_column")) for row in rows)
    if actual != expected:
        raise ValueError("Batch summary does not account for every requested target")
    for row in rows:
        if row.get("status") not in {"ok", "skipped_existing"}:
            raise ValueError(f"Target {row.get('query')} is incomplete: {row.get('error') or row.get('status')}")
        report_path = _resolve_existing_path(row.get("json_report"), summary_path.parent)
        if report_path is None or not report_path.is_file():
            raise ValueError(f"Completed target has no report: {row.get('query')}")
        json.loads(report_path.read_text())


def _count_tsv_rows(path: Path) -> int:
    if not path.exists():
        return 0
    with path.open("r", encoding="utf-8", newline="") as handle:
        return max(0, sum(1 for _ in handle) - 1)


def _count_json_list(path: Path) -> int:
    if not path.exists():
        return 0
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return 0
    return len(payload) if isinstance(payload, list) else 0


def _resolve_default_af3_catalog(value: str | Path | None) -> Path | None:
    if value is not None:
        candidate = Path(value)
        if candidate.exists():
            return candidate
    return None


def _log(verbose: bool, message: str) -> None:
    if verbose:
        print(f"[openantigen-deploy] {message}", flush=True)
