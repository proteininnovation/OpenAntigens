from __future__ import annotations

import argparse
import os
from pathlib import Path

from .assets import generate_assets_for_summary
from .af3 import catalogue_af3_structures, prepare_af3_structures
from .batch import rebuild_batch_summary_from_tsv, retry_batch_from_tsv, run_batch_from_tsv, write_batch_summary_markdown
from .config import AnalysisConfig
from .deployment import build_fresh_snapshot, compress_public_site, prepare_deployment_data
from .family_alignments import build_family_alignments_from_tsv
from .module_refresh import refresh_report_modules
from .mouse_portal import build_mouse_portal_from_human_orthologs
from .open_targets import (
    OPEN_TARGETS_DEFAULT_RELEASE,
    build_open_targets_associations_from_bulk_downloads,
    build_open_targets_associations_from_summary,
)
from .ortholog_table import build_ortholog_table_from_tsv
from .paralogs import build_paralog_reference_from_tsv
from .pipeline import AntigenAnalyzer
from .portal import build_portal
from .refresh_job import run_accessible_portal_refresh_job
from .target_sets import DEFAULT_ACCESSIBLE_BUCKETS, build_accessible_target_tsv
from .test_subset import build_test_portal_subset
from .uniprot_prefetch import prefetch_reviewed_taxa


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="AgDesign2 portal and antigen-design data workflow")
    subparsers = parser.add_subparsers(dest="command")

    analyze = subparsers.add_parser("analyze", help="Analyze a target protein")
    analyze.add_argument("query", help="UniProt accession, entry name, or official human gene symbol")
    analyze.add_argument("--output-dir", default="outputs", help="Directory for report JSON output")
    analyze.add_argument("--verbose", action="store_true", help="Print stage-by-stage progress updates")
    analyze.add_argument("--complex-portal", action="store_true", help="Enable Complex Portal lookups")

    analyze_tsv = subparsers.add_parser("analyze-tsv", help="Analyze many targets listed in a CSV or TSV file")
    analyze_tsv.add_argument("tsv_path", help="CSV or TSV with a target identifier column")
    analyze_tsv.add_argument("--output-dir", default="outputs/batch", help="Directory for batch reports")
    analyze_tsv.add_argument("--limit", type=int, default=None, help="Only process the first N rows")
    analyze_tsv.add_argument("--no-resume", action="store_true", help="Re-run targets even if report files already exist")
    analyze_tsv.add_argument("--jobs", type=int, default=1, help="Parallel target-analysis workers (0 = auto)")
    analyze_tsv.add_argument("--verbose", action="store_true", help="Print stage-by-stage progress updates")
    analyze_tsv.add_argument("--complex-portal", action="store_true", help="Enable Complex Portal lookups")

    retry_batch = subparsers.add_parser(
        "retry-batch",
        help="Re-run batch rows that errored or are missing a local AlphaFold structure",
    )
    retry_batch.add_argument("tsv_path", help="CSV or TSV with a target identifier column")
    retry_batch.add_argument("--output-dir", default="outputs/batch", help="Directory for batch reports")
    retry_batch.add_argument(
        "--summary-path",
        default=None,
        help="Existing batch_summary.json to update in place (defaults to OUTPUT_DIR/batch_summary.json)",
    )
    retry_batch.add_argument(
        "--no-alphafold-check",
        action="store_true",
        help="Only retry explicit errors and ignore rows that are merely missing a local AlphaFold structure",
    )
    retry_batch.add_argument("--verbose", action="store_true", help="Print stage-by-stage progress updates")
    retry_batch.add_argument("--complex-portal", action="store_true", help="Enable Complex Portal lookups")

    summarize_batch = subparsers.add_parser("summarize-batch", help="Render a Markdown dashboard from batch_summary.json")
    summarize_batch.add_argument("summary_path", help="Path to batch_summary.json")
    summarize_batch.add_argument("--output-path", default=None, help="Optional destination for the Markdown dashboard")

    generate_assets = subparsers.add_parser("generate-assets", help="Generate or refresh construct images/plots for report JSON files")
    generate_assets.add_argument("summary_path", help="Path to batch_summary.json")
    generate_assets.add_argument("--data-dir", default=".agdesign2/data", help="Directory containing AlphaFold artifacts")
    generate_assets.add_argument("--jobs", type=int, default=1, help="Parallel asset-generation workers (0 = auto)")
    generate_assets.add_argument("--no-resume", action="store_true", help="Regenerate assets even when paths already exist")
    generate_assets.add_argument("--no-structure-images", action="store_true", help="Skip PyMOL structure PNG rendering")
    generate_assets.add_argument("--no-quality-plots", action="store_true", help="Skip pLDDT/PAE quality plot rendering")
    generate_assets.add_argument("--verbose", action="store_true", help="Print progress")

    prepare_af3 = subparsers.add_parser(
        "prepare-af3-structures",
        help="Unpack AlphaFold Server / AF3 structure exports into snapshot-local fallback artifacts",
    )
    prepare_af3.add_argument("zip_path", help="AlphaFold Server export ZIP")
    prepare_af3.add_argument(
        "--snapshot-dir",
        default="snapshots/openantigen-2026-05-02",
        help="Snapshot directory containing data/batch_summary.json",
    )
    prepare_af3.add_argument(
        "--output-dir",
        default=None,
        help="Destination directory for converted AF3 artifacts (default: SNAPSHOT_DIR/data/alphafold3)",
    )
    prepare_af3.add_argument("--verbose", action="store_true", help="Print AF3 target conversion progress")

    catalogue_af3 = subparsers.add_parser(
        "catalogue-af3-structures",
        help="Catalogue AlphaFold Server / AF3 ZIP exports into SQLite plus reusable fallback artifacts",
    )
    catalogue_af3.add_argument("input_path", help="AF3 ZIP file or directory containing AF3 ZIP exports")
    catalogue_af3.add_argument(
        "--output-dir",
        default="data/af3_catalog",
        help="Destination catalogue directory (default: data/af3_catalog)",
    )
    catalogue_af3.add_argument("--verbose", action="store_true", help="Print archive-level catalogue progress")

    rebuild_batch_summary = subparsers.add_parser(
        "rebuild-batch-summary",
        help="Reconstruct batch_summary.json and batch_summary.md from a target list plus existing report files",
    )
    rebuild_batch_summary.add_argument(
        "tsv_path",
        nargs="?",
        default="data/inputs/accessible_human_topology_refined_accessible_only.csv",
        help="CSV or TSV describing the original batch input order (default: canonical accessible-human CSV)",
    )
    rebuild_batch_summary.add_argument("--output-dir", default="outputs/surfy_batch", help="Directory containing batch report files")
    rebuild_batch_summary.add_argument(
        "--summary-path",
        default=None,
        help="Destination batch_summary.json path (defaults to OUTPUT_DIR/batch_summary.json)",
    )

    build_accessible_targets = subparsers.add_parser(
        "build-accessible-targets",
        help="Build a TSV of accessible human targets from the refined topology CSV",
    )
    build_accessible_targets.add_argument(
        "--source-csv",
        default="data/inputs/accessible_human_topology_refined_accessible_only.csv",
        help="Refined accessible-topology CSV",
    )
    build_accessible_targets.add_argument(
        "--output-tsv",
        default="accessible_secreted_gpi_singlepass.tsv",
        help="Destination TSV for filtered targets",
    )
    build_accessible_targets.add_argument(
        "--buckets",
        nargs="+",
        default=list(DEFAULT_ACCESSIBLE_BUCKETS),
        help="Topology buckets to include (default: Secreted GPI Single-pass Multipass)",
    )
    mouse_from_human = subparsers.add_parser(
        "mouse-from-human-orthologs",
        help="Build OpenAntigens Mouse from annotated mouse orthologs in a human OpenAntigens batch",
    )
    mouse_from_human.add_argument("summary_path", help="Human batch_summary.json")
    mouse_from_human.add_argument("--output-dir", required=True, help="Destination mouse batch/portal directory")
    mouse_from_human.add_argument(
        "--no-build-site",
        action="store_true",
        help="Only write mouse batch_summary.json/report files and skip portal HTML",
    )
    mouse_from_human.add_argument(
        "--bundle-vendor-assets",
        action="store_true",
        help="Copy/download browser libraries into the mouse portal directory",
    )
    mouse_from_human.add_argument(
        "--fetch-mouse-alphafold",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Fetch mouse AlphaFold DB structures when missing locally (default: on; use --no-fetch-mouse-alphafold to use local compatible artifacts only)",
    )
    mouse_from_human.add_argument(
        "--af3-catalog",
        default=None,
        help="Optional AF3 catalogue SQLite DB; AF3 fallback is disabled unless this option is supplied.",
    )
    mouse_from_human.add_argument("--jobs", type=int, default=1, help="Parallel mouse target-analysis workers")
    mouse_from_human.add_argument(
        "--no-resume",
        action="store_true",
        help="Re-analyze every mouse target even if reports exist (re-analysis is also forced automatically when the pipeline version changes)",
    )
    mouse_from_human.add_argument("--verbose", action="store_true", help="Print included target progress")

    refresh_portal_job = subparsers.add_parser(
        "refresh-accessible-portal",
        help="Run a resumable full refresh job for the accessible Secreted/GPI/Single-pass/Multipass portal dataset",
    )
    refresh_portal_job.add_argument(
        "--source-csv",
        default="data/inputs/accessible_human_topology_refined_accessible_only.csv",
        help="Refined accessible-topology CSV",
    )
    refresh_portal_job.add_argument(
        "--output-tsv",
        default="snapshots/local-refresh/accessible_secreted_gpi_singlepass.tsv",
        help="Filtered target TSV path",
    )
    refresh_portal_job.add_argument(
        "--output-dir",
        default="snapshots/local-refresh/outputs/surfy_batch",
        help="Directory containing refreshed reports and the portal",
    )
    refresh_portal_job.add_argument(
        "--state-path",
        default=None,
        help="Persistent refresh state JSON",
    )
    refresh_portal_job.add_argument(
        "--buckets",
        nargs="+",
        default=list(DEFAULT_ACCESSIBLE_BUCKETS),
        help="Topology buckets to include (default: Secreted GPI Single-pass Multipass)",
    )
    refresh_portal_job.add_argument(
        "--portal-interval",
        type=int,
        default=25,
        help="Rebuild the static portal every N processed targets",
    )
    refresh_portal_job.add_argument(
        "--summary-interval",
        type=int,
        default=10,
        help="Rebuild batch_summary.json every N processed targets",
    )
    refresh_portal_job.add_argument(
        "--max-targets",
        type=int,
        default=None,
        help="Optional cap on the number of targets processed in this invocation",
    )
    refresh_portal_job.add_argument(
        "--max-attempts",
        type=int,
        default=None,
        help="Optional cap on retries of a failed target across invocations; "
        "exhausted targets are skipped so they stop consuming the per-run budget",
    )
    refresh_portal_job.add_argument("--verbose", action="store_true", help="Print row-by-row progress updates")
    refresh_portal_job.add_argument("--complex-portal", action="store_true", help="Enable Complex Portal lookups")

    build_portal_parser = subparsers.add_parser("build-portal", help="Build a static HTML portal from batch outputs")
    build_portal_parser.add_argument("summary_path", help="Path to batch_summary.json")
    build_portal_parser.add_argument("--output-dir", default=None, help="Optional portal output directory")
    build_portal_parser.add_argument(
        "--refresh-reports",
        action="store_true",
        help="Re-run all reports from batch_summary.json before building the portal",
    )
    build_portal_parser.add_argument(
        "--no-refresh-stale-reports",
        action="store_true",
        help="Do not automatically refresh reports that are older than the precomputed ortholog/family data",
    )
    build_portal_parser.add_argument("--verbose", action="store_true", help="Print refresh progress updates")
    build_portal_parser.add_argument("--complex-portal", action="store_true", help="Enable Complex Portal during report refreshes")
    build_portal_parser.add_argument(
        "--bundle-vendor-assets",
        action="store_true",
        help="Copy/download browser libraries into the portal directory for a self-contained static site",
    )

    build_portal_db_parser = subparsers.add_parser(
        "build-portal-db",
        help="Build a SQLite portal database from existing batch JSON/report files",
    )
    build_portal_db_parser.add_argument("summary_path", help="Path to batch_summary.json")
    build_portal_db_parser.add_argument("--db", required=True, help="Destination SQLite database path")
    build_portal_db_parser.add_argument(
        "--assets-dir",
        default=None,
        help="Directory for copied portal assets (default: DB parent/assets)",
    )
    build_portal_db_parser.add_argument(
        "--bundle-vendor-assets",
        action="store_true",
        help="Copy/download browser libraries into the SQLite portal assets directory",
    )

    serve_portal_db_parser = subparsers.add_parser(
        "serve-portal-db",
        help="Serve a dynamic portal from a SQLite portal database",
    )
    serve_portal_db_parser.add_argument("db", help="Path to openantigen.sqlite")
    serve_portal_db_parser.add_argument(
        "--assets-dir",
        default=None,
        help="Directory containing portal assets (default: assets_dir recorded in the DB)",
    )
    serve_portal_db_parser.add_argument("--host", default="127.0.0.1", help="Host interface for the local web server")
    serve_portal_db_parser.add_argument("--port", type=int, default=8000, help="Port for the local web server")

    build_portal_from_db_parser = subparsers.add_parser(
        "build-portal-from-db",
        help="Build a static HTML portal from a SQLite portal database",
    )
    build_portal_from_db_parser.add_argument("db", help="Path to openantigen.sqlite")
    build_portal_from_db_parser.add_argument("--output-dir", required=True, help="Destination portal output directory")
    build_portal_from_db_parser.add_argument(
        "--assets-dir",
        default=None,
        help="Directory containing portal assets (default: assets_dir recorded in the DB)",
    )
    build_portal_from_db_parser.add_argument(
        "--bundle-vendor-assets",
        action="store_true",
        help="Copy/download browser libraries into the portal directory for a self-contained static site",
    )

    build_test_portal = subparsers.add_parser(
        "build-test-portal",
        help="Build a compact portal snapshot with representative target classes for upload testing",
    )
    build_test_portal.add_argument(
        "summary_path",
        nargs="?",
        default="snapshots/openantigen-2026-05-02/data/batch_summary.json",
        help="Source batch_summary.json (default: snapshots/openantigen-2026-05-02/data/batch_summary.json)",
    )
    build_test_portal.add_argument(
        "--output-dir",
        default="snapshots/openantigen-test-10",
        help="Destination snapshot root containing data/ and public_site/ (default: snapshots/openantigen-test-10)",
    )
    build_test_portal.add_argument(
        "--per-kind",
        type=int,
        default=10,
        help="Targets to select per kind: multipass, GPCR, single-pass, GPI, secreted, obligatory partner",
    )
    build_test_portal.add_argument(
        "--no-build-portal",
        action="store_true",
        help="Only write the subset data files and skip rendering public_site/",
    )
    build_test_portal.add_argument(
        "--bundle-vendor-assets",
        action="store_true",
        help="Copy/download browser libraries into the test portal for a self-contained static site",
    )
    build_test_portal.add_argument("--verbose", action="store_true", help="Print copied support files")

    serve_portal_parser = subparsers.add_parser(
        "serve-portal",
        help="Serve a dynamic portal directly from batch JSON/report files",
    )
    serve_portal_parser.add_argument("summary_path", help="Path to batch_summary.json")
    serve_portal_parser.add_argument("--host", default="127.0.0.1", help="Host interface for the local web server")
    serve_portal_parser.add_argument("--port", type=int, default=8000, help="Port for the local web server")
    serve_portal_parser.add_argument(
        "--refresh-reports",
        action="store_true",
        help="Re-run all reports from batch_summary.json before serving the portal",
    )
    serve_portal_parser.add_argument(
        "--refresh-stale-reports",
        action="store_true",
        help="Automatically refresh reports that are older than the precomputed ortholog/family data before serving",
    )
    serve_portal_parser.add_argument(
        "--no-refresh-stale-reports",
        action="store_true",
        help=argparse.SUPPRESS,
    )
    serve_portal_parser.add_argument("--verbose", action="store_true", help="Print refresh progress updates")
    serve_portal_parser.add_argument("--complex-portal", action="store_true", help="Enable Complex Portal during report refreshes")

    refresh_module_parser = subparsers.add_parser(
        "refresh-report-module",
        help="Recompute one or more report modules in place without rerunning the full analysis",
    )
    refresh_module_parser.add_argument("summary_path", help="Path to batch_summary.json")
    refresh_module_parser.add_argument(
        "--module",
        dest="modules",
        action="append",
        required=True,
        help="Module to recompute. Repeat for multiple modules. Supported: complex-portal, cross-reactivity, family-context, target-metadata",
    )
    refresh_module_parser.add_argument(
        "--query",
        dest="queries",
        action="append",
        default=[],
        help="Optional query or resolved entry name to limit the refresh to specific targets",
    )
    refresh_module_parser.add_argument(
        "--only-zero-cross-reactivity",
        action="store_true",
        help="When refreshing cross-reactivity, skip reports that already have relevant BLAST hits",
    )
    refresh_module_parser.add_argument("--verbose", action="store_true", help="Print row-by-row progress updates")

    open_targets_parser = subparsers.add_parser(
        "build-open-targets",
        help="Precompute Open Targets disease association scores for portal targets",
    )
    open_targets_parser.add_argument(
        "summary_path",
        nargs="?",
        default="outputs/surfy_batch/batch_summary.json",
        help="Path to batch_summary.json (default: outputs/surfy_batch/batch_summary.json)",
    )
    open_targets_parser.add_argument(
        "--output-dir",
        default=None,
        help="Destination directory for open_targets_disease_associations TSV/JSON (defaults to summary directory)",
    )
    open_targets_parser.add_argument("--page-size", type=int, default=100, help="Open Targets association page size")
    open_targets_parser.add_argument("--max-pages", type=int, default=5, help="Maximum association pages per target")
    open_targets_parser.add_argument("--limit", type=int, default=None, help="Optional number of uncached targets to process")
    open_targets_parser.add_argument("--no-resume", action="store_true", help="Rebuild from scratch instead of resuming")
    open_targets_parser.add_argument(
        "--source",
        choices=["bulk", "api"],
        default="bulk",
        help="Use Open Targets bulk release downloads or per-target GraphQL API calls (default: bulk)",
    )
    open_targets_parser.add_argument(
        "--release",
        default=OPEN_TARGETS_DEFAULT_RELEASE,
        help=f"Open Targets bulk release version, or 'latest' to auto-resolve the newest release (default: {OPEN_TARGETS_DEFAULT_RELEASE})",
    )
    open_targets_parser.add_argument(
        "--bulk-data-dir",
        default=".agdesign2/data/open_targets",
        help="Directory for downloaded Open Targets bulk Parquet datasets",
    )
    open_targets_parser.add_argument(
        "--no-download",
        action="store_true",
        help="Use existing local Open Targets bulk datasets without downloading",
    )
    open_targets_parser.add_argument("--verbose", action="store_true", help="Print target-by-target progress")

    ortholog_table = subparsers.add_parser(
        "build-ortholog-table",
        help="Precompute a human/mouse/cynomolgus monkey ortholog and RefSeq reference table from a target list",
    )
    ortholog_table.add_argument(
        "tsv_path",
        nargs="?",
        default="data/inputs/accessible_human_topology_refined_accessible_only.csv",
        help="CSV or TSV with a target identifier column (default: canonical accessible-human CSV)",
    )
    ortholog_table.add_argument(
        "--output-path",
        default="outputs/surfy_batch/ortholog_reference_table.tsv",
        help="Destination TSV path for the ortholog reference table",
    )
    ortholog_table.add_argument(
        "--no-resume",
        action="store_true",
        help="Rebuild all rows from scratch instead of resuming from an existing partial table",
    )
    ortholog_table.add_argument(
        "--force",
        action="store_true",
        help="Alias for --no-resume",
    )
    ortholog_table.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="Number of concurrent worker threads for the network-bound per-target "
        "lookups (default: 1 = sequential). Output is identical regardless of value; "
        "keep modest (e.g. 4-8) to avoid NCBI/UniProt rate limiting.",
    )
    ortholog_table.add_argument("--verbose", action="store_true", help="Print row-by-row progress updates")

    family_alignments = subparsers.add_parser(
        "build-family-alignments",
        help="Precompute ectodomain-only human family alignments and directional identity matrices from a target list",
    )
    family_alignments.add_argument(
        "tsv_path",
        nargs="?",
        default="data/inputs/accessible_human_topology_refined_accessible_only.csv",
        help="CSV or TSV with human targets (default: canonical accessible-human CSV)",
    )
    family_alignments.add_argument(
        "--output-dir",
        default="outputs/surfy_batch/family_alignments",
        help="Destination directory for per-family alignment outputs",
    )
    family_alignments.add_argument("--verbose", action="store_true", help="Print row-by-row progress updates")

    paralog_reference = subparsers.add_parser(
        "build-paralog-reference",
        help="Precompute broader human paralog sets and directional ectodomain matrices from a target list",
    )
    paralog_reference.add_argument(
        "tsv_path",
        nargs="?",
        default="data/inputs/accessible_human_topology_refined_accessible_only.csv",
        help="CSV or TSV with human targets (default: canonical accessible-human CSV)",
    )
    paralog_reference.add_argument(
        "--output-dir",
        default="outputs/surfy_batch/paralogs",
        help="Destination directory for per-target paralog outputs",
    )
    paralog_reference.add_argument("--verbose", action="store_true", help="Print row-by-row progress updates")

    setup = subparsers.add_parser("setup-blast", help="Download Swiss-Prot FASTA files and build BLAST DBs")
    setup.add_argument("--data-dir", default=".agdesign2/data", help="Location for cached data and BLAST DBs")

    prefetch_uniprot = subparsers.add_parser(
        "prefetch-uniprot",
        help="Bulk-prefetch reviewed UniProt entries by taxonomy into local cache",
    )
    prefetch_uniprot.add_argument(
        "--taxa",
        nargs="+",
        type=int,
        default=[9606, 10090],
        help="Taxonomy IDs to prefetch (default: 9606 10090)",
    )
    prefetch_uniprot.add_argument(
        "--cache-dir",
        default=".agdesign2/cache",
        help="Cache directory used by the pipeline",
    )
    prefetch_uniprot.add_argument(
        "--page-size",
        type=int,
        default=200,
        help="UniProt page size for bulk search",
    )
    prefetch_uniprot.add_argument(
        "--max-pages",
        type=int,
        default=None,
        help="Optional safety cap on pages per taxonomy",
    )
    prepare_deployment = subparsers.add_parser(
        "prepare-deployment-data",
        help="Run the bulk/precomputed data steps needed before publishing the portal",
    )
    prepare_deployment.add_argument(
        "--target-tsv",
        default="accessible_secreted_gpi_singlepass.tsv",
        help="Target TSV used for ortholog/family/paralog precomputes",
    )
    prepare_deployment.add_argument(
        "--summary-path",
        default="outputs/surfy_batch/batch_summary.json",
        help="Batch summary used for Open Targets, cross-reactivity refresh, and portal build",
    )
    prepare_deployment.add_argument("--output-dir", default="outputs/surfy_batch", help="Portal output/data directory")
    prepare_deployment.add_argument("--cache-dir", default=".agdesign2/cache", help="Pipeline cache directory")
    prepare_deployment.add_argument("--data-dir", default=".agdesign2/data", help="Pipeline data directory")
    prepare_deployment.add_argument("--taxa", nargs="+", type=int, default=[9606, 10090], help="UniProt taxa to bulk-prefetch")
    prepare_deployment.add_argument("--uniprot-page-size", type=int, default=500, help="UniProt bulk prefetch page size")
    prepare_deployment.add_argument("--open-targets-page-size", type=int, default=100, help="Open Targets page size")
    prepare_deployment.add_argument("--open-targets-max-pages", type=int, default=5, help="Open Targets max pages per target")
    prepare_deployment.add_argument("--skip-uniprot", action="store_true", help="Skip reviewed UniProt bulk prefetch")
    prepare_deployment.add_argument("--skip-blast-setup", action="store_true", help="Skip BLAST database setup")
    prepare_deployment.add_argument("--skip-orthologs", action="store_true", help="Skip ortholog reference table")
    prepare_deployment.add_argument("--skip-family-alignments", action="store_true", help="Skip family alignment precompute")
    prepare_deployment.add_argument("--skip-paralogs", action="store_true", help="Skip paralog reference precompute")
    prepare_deployment.add_argument("--skip-open-targets", action="store_true", help="Skip Open Targets precompute")
    prepare_deployment.add_argument("--skip-cross-reactivity", action="store_true", help="Skip bulk BLAST cross-reactivity refresh")
    prepare_deployment.add_argument("--skip-portal", action="store_true", help="Skip static portal rebuild")
    prepare_deployment.add_argument("--dry-run", action="store_true", help="Print planned steps without running them")
    prepare_deployment.add_argument("--verbose", action="store_true", help="Print progress")
    fresh_snapshot = subparsers.add_parser(
        "build-fresh-snapshot",
        help="Build a complete fresh OpenAntigens website snapshot using clean snapshot-local cache/data directories",
    )
    fresh_snapshot.add_argument(
        "--source-csv",
        default="data/inputs/accessible_human_topology_refined_accessible_only.csv",
        help="Surface/secreted protein universe CSV; this is the only local biological input reused",
    )
    fresh_snapshot.add_argument("--snapshot-root", default="snapshots", help="Directory where fresh snapshots are created")
    fresh_snapshot.add_argument("--snapshot-name", default=None, help="Optional snapshot folder name")
    fresh_snapshot.add_argument(
        "--buckets",
        nargs="+",
        default=list(DEFAULT_ACCESSIBLE_BUCKETS),
        help="Topology buckets to include",
    )
    fresh_snapshot.add_argument("--limit", type=int, default=None, help="Optional target count for benchmarking/subsets")
    fresh_snapshot.add_argument("--sample-size", type=int, default=None, help="Pinned smoke target count; activates stratified random sampling")
    fresh_snapshot.add_argument("--sample-seed", type=int, default=20260530, help="Random seed for pinned smoke sampling")
    fresh_snapshot.add_argument(
        "--sample-mode",
        choices=["stratified-random"],
        default=None,
        help="Sampling mode for pinned smoke snapshots",
    )
    fresh_snapshot.add_argument(
        "--include-mouse",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Build structure-analyzed OpenAntigens Mouse under public_site/mouse (default: on; use --no-include-mouse to skip)",
    )
    fresh_snapshot.add_argument(
        "--fetch-mouse-alphafold",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Fetch mouse AlphaFold DB structures when missing locally (default: on; use --no-fetch-mouse-alphafold to use local compatible artifacts only)",
    )
    fresh_snapshot.add_argument("--force", action="store_true", help="Replace an existing snapshot directory")
    fresh_snapshot.add_argument("--resume", action="store_true", help="Resume an existing incomplete snapshot directory")
    fresh_snapshot.add_argument(
        "--reanalyze-reports",
        action="store_true",
        help="With --resume, force re-analysis of all target reports (reusing the cached "
        "UniProt/AlphaFold/BLAST and precomputed ortholog/family/paralog artifacts) so a "
        "report-logic change is picked up without a cold rebuild. Pair with "
        "--no-structure-images/--no-quality-plots to avoid asset regrowth.",
    )
    fresh_snapshot.add_argument("--taxa", nargs="+", type=int, default=[9606, 10090], help="UniProt taxa to bulk-prefetch")
    fresh_snapshot.add_argument(
        "--uniprot-prefetch-mode",
        choices=["targets", "taxa"],
        default="targets",
        help="UniProt fresh prefetch mode. 'targets' bulk-fetches only target TSV entries; 'taxa' mirrors reviewed taxa.",
    )
    fresh_snapshot.add_argument(
        "--uniprot-page-size",
        type=int,
        default=500,
        help="UniProt bulk page size (target-list queries are capped at the API limit of 100 identifiers)",
    )
    fresh_snapshot.add_argument("--open-targets-page-size", type=int, default=100, help="Open Targets page size")
    fresh_snapshot.add_argument("--open-targets-max-pages", type=int, default=5, help="Open Targets max pages per target")
    fresh_snapshot.add_argument("--jobs", type=int, default=0, help="Parallel target-analysis workers (0 = auto)")
    fresh_snapshot.add_argument("--asset-jobs", type=int, default=0, help="Parallel asset-generation workers (0 = use --jobs/auto)")
    fresh_snapshot.add_argument(
        "--ortholog-jobs",
        type=int,
        default=1,
        help="Concurrent worker threads for the network-bound ortholog-table build "
        "(default: 1 = sequential). Output is identical regardless; keep modest "
        "(e.g. 4-8) to avoid NCBI/UniProt rate limiting.",
    )
    fresh_snapshot.add_argument(
        "--af3-catalog",
        default=None,
        help="Optional AF3 catalogue SQLite DB to copy snapshot-matching fallback artifacts from before analysis",
    )
    structure_images = fresh_snapshot.add_mutually_exclusive_group()
    structure_images.add_argument(
        "--structure-images",
        dest="render_structure_images",
        action="store_true",
        help="Opt in to static PyMOL construct structure PNGs",
    )
    structure_images.add_argument(
        "--no-structure-images",
        dest="render_structure_images",
        action="store_false",
        help="Do not render static PyMOL construct structure PNGs (default)",
    )
    fresh_snapshot.set_defaults(render_structure_images=False)
    fresh_snapshot.add_argument("--no-quality-plots", action="store_true", help="Skip static pLDDT/PAE construct quality plots")
    fresh_snapshot.add_argument("--dry-run", action="store_true", help="Print planned paths/steps without running work")
    fresh_snapshot.add_argument("--verbose", action="store_true", help="Print progress")

    compress_site = subparsers.add_parser(
        "compress-public-site",
        help="Gzip browser-fetched assets (.js/.css/.svg) on disk and optionally "
        "prune dead weight, to fit a disk-quota-limited static host. Idempotent.",
    )
    compress_site.add_argument("public_site_dir", help="Path to a packaged public_site/ directory")
    compress_site.add_argument(
        "--include-html",
        action="store_true",
        help="Also gzip report HTML on disk (saves more, but a .gz-only report page "
        "404s for the rare non-gzip crawler/monitor; directory index.html is always "
        "left intact). Off by default.",
    )
    compress_site.add_argument(
        "--drop-unreferenced-structures",
        action="store_true",
        help="Remove structures/ and mouse/structures/ — the 3D viewer inlines PDB "
        "text into each report's JS, so those copies are linked from nothing.",
    )
    compress_site.add_argument(
        "--drop-redundant-downloads",
        action="store_true",
        help="Remove the Open Targets .json flat file (the .tsv holds the same data) "
        "and strip its dead button from downloads.html.",
    )
    compress_site.add_argument("--verbose", action="store_true", help="Print progress")
    return parser


def main() -> None:
    os.environ.setdefault("MPLCONFIGDIR", str(Path(".agdesign2/mplconfig").resolve()))
    parser = build_parser()
    args = parser.parse_args()
    if args.command == "setup-blast":
        data_dir = Path(args.data_dir)
        config = AnalysisConfig(
            data_dir=data_dir,
            blast_db_dir=data_dir / "blastdb",
            ortholog_fasta_dir=data_dir / "proteomes",
            ortholog_blast_db_dir=data_dir / "ortholog_blastdb",
        )
        analyzer = AntigenAnalyzer(config=config)
        analyzer.blast_client.ensure_databases()
        analyzer.blast_client.ensure_ortholog_databases()
        print(f"Cross-reactivity BLAST databases prepared in {config.blast_db_dir}")
        print(f"Ortholog proteome BLAST databases prepared in {config.ortholog_blast_db_dir}")
        return
    if args.command == "prefetch-uniprot":
        summary = prefetch_reviewed_taxa(
            cache_dir=Path(args.cache_dir),
            taxa=tuple(args.taxa),
            page_size=args.page_size,
            max_pages=args.max_pages,
        )
        print(f"UniProt pages fetched: {summary['pages']}")
        print(f"UniProt entries processed: {summary['entries']}")
        print(f"UniProt entry cache files written: {summary['written_entry_cache']}")
        return
    if args.command == "prepare-deployment-data":
        result = prepare_deployment_data(
            target_tsv=args.target_tsv,
            summary_path=args.summary_path,
            output_dir=args.output_dir,
            cache_dir=args.cache_dir,
            data_dir=args.data_dir,
            taxa=tuple(args.taxa),
            uniprot_page_size=args.uniprot_page_size,
            open_targets_page_size=args.open_targets_page_size,
            open_targets_max_pages=args.open_targets_max_pages,
            skip_uniprot=args.skip_uniprot,
            skip_blast_setup=args.skip_blast_setup,
            skip_orthologs=args.skip_orthologs,
            skip_family_alignments=args.skip_family_alignments,
            skip_paralogs=args.skip_paralogs,
            skip_open_targets=args.skip_open_targets,
            skip_cross_reactivity=args.skip_cross_reactivity,
            skip_portal=args.skip_portal,
            dry_run=args.dry_run,
            verbose=args.verbose,
        )
        print(result.summary_path)
        print(result.output_dir)
        print("Planned steps:")
        for step in result.steps_planned:
            print(f"- {step}")
        if result.steps_completed:
            print("Completed steps:")
            for step in result.steps_completed:
                print(f"- {step}")
        return
    if args.command == "build-fresh-snapshot":
        result = build_fresh_snapshot(
            source_csv=args.source_csv,
            snapshot_root=args.snapshot_root,
            snapshot_name=args.snapshot_name,
            allowed_buckets=tuple(args.buckets),
            limit=args.limit,
            sample_size=args.sample_size,
            sample_seed=args.sample_seed,
            sample_mode=args.sample_mode,
            include_mouse=args.include_mouse,
            fetch_mouse_alphafold=args.fetch_mouse_alphafold,
            force=args.force,
            resume=args.resume,
            reanalyze_reports=args.reanalyze_reports,
            taxa=tuple(args.taxa),
            uniprot_prefetch_mode=args.uniprot_prefetch_mode,
            uniprot_page_size=args.uniprot_page_size,
            open_targets_page_size=args.open_targets_page_size,
            open_targets_max_pages=args.open_targets_max_pages,
            jobs=args.jobs,
            asset_jobs=args.asset_jobs,
            ortholog_jobs=args.ortholog_jobs,
            af3_catalog=args.af3_catalog,
            render_structure_images=args.render_structure_images,
            render_quality_plots=not args.no_quality_plots,
            dry_run=args.dry_run,
            verbose=args.verbose,
        )
        print(result.snapshot_dir)
        print(result.target_tsv)
        print(result.summary_path)
        print(result.portal_dir)
        print(result.public_site_dir)
        print(result.benchmark_path)
        print(result.manifest_path)
        if result.selected_smoke_tsv is not None:
            print(result.selected_smoke_tsv)
        if result.full_smoke_manifest_path is not None:
            print(result.full_smoke_manifest_path)
        if result.mouse_public_site_dir is not None:
            print(result.mouse_public_site_dir)
        for row in result.timings:
            status = row.get("status")
            duration = row.get("duration_seconds")
            suffix = f" {duration}s" if duration is not None else ""
            print(f"{row.get('step')}: {status}{suffix}")
        return
    if args.command == "compress-public-site":
        summary = compress_public_site(
            Path(args.public_site_dir),
            gzip_html=args.include_html,
            drop_unreferenced_structures=args.drop_unreferenced_structures,
            drop_redundant_downloads=args.drop_redundant_downloads,
            verbose=args.verbose,
        )
        saved = summary["bytes_before"] - summary["bytes_after"]
        total = saved + summary["structures_bytes_removed"] + summary["downloads_bytes_removed"]
        print(f"files_compressed: {summary['files_compressed']}")
        print(f"compression_bytes_saved: {saved}")
        print(f"structures_bytes_removed: {summary['structures_bytes_removed']}")
        print(f"downloads_bytes_removed: {summary['downloads_bytes_removed']}")
        print(f"total_bytes_reclaimed: {total}")
        return
    if args.command == "analyze":
        analyzer = AntigenAnalyzer(
            config=AnalysisConfig(
                verbose_progress=args.verbose,
                enable_complex_portal=args.complex_portal,
            )
        )
        report = analyzer.analyze_target(args.query, output_dir=args.output_dir)
        output_root = Path(args.output_dir)
        stem = f"{report.target.entry_name.lower()}_report"
        print(output_root / f"{stem}.json")
        print(output_root / f"{stem}.md")
        return
    if args.command == "analyze-tsv":
        analyzer = AntigenAnalyzer(
            config=AnalysisConfig(
                verbose_progress=args.verbose,
                enable_complex_portal=args.complex_portal,
            )
        )
        results = run_batch_from_tsv(
            analyzer=analyzer,
            tsv_path=args.tsv_path,
            output_dir=args.output_dir,
            limit=args.limit,
            resume=not args.no_resume,
            jobs=args.jobs,
        )
        print(Path(args.output_dir) / "batch_summary.json")
        print(f"Processed {len(results)} rows")
        return
    if args.command == "retry-batch":
        analyzer = AntigenAnalyzer(
            config=AnalysisConfig(
                verbose_progress=args.verbose,
                enable_complex_portal=args.complex_portal,
            )
        )
        results = retry_batch_from_tsv(
            analyzer=analyzer,
            tsv_path=args.tsv_path,
            output_dir=args.output_dir,
            summary_path=args.summary_path,
            require_alphafold=not args.no_alphafold_check,
        )
        summary_path = Path(args.summary_path) if args.summary_path is not None else Path(args.output_dir) / "batch_summary.json"
        print(summary_path)
        print(f"Tracked {len(results)} rows")
        return
    if args.command == "summarize-batch":
        output_path = write_batch_summary_markdown(args.summary_path, args.output_path)
        print(output_path)
        return
    if args.command == "generate-assets":
        data_dir = Path(args.data_dir)
        results = generate_assets_for_summary(
            args.summary_path,
            config=AnalysisConfig(
                data_dir=data_dir,
                blast_db_dir=data_dir / "blastdb",
                ortholog_fasta_dir=data_dir / "proteomes",
                ortholog_blast_db_dir=data_dir / "ortholog_blastdb",
                render_structure_images=not args.no_structure_images,
                render_quality_plots=not args.no_quality_plots,
                verbose_progress=args.verbose,
            ),
            jobs=args.jobs,
            resume=not args.no_resume,
            verbose=args.verbose,
        )
        ok = sum(1 for item in results if item.get("status") == "ok" or str(item.get("status") or "").startswith("skipped"))
        errors = sum(1 for item in results if item.get("status") == "error")
        print(f"Asset jobs completed: {ok}; errors: {errors}")
        return
    if args.command == "prepare-af3-structures":
        result = prepare_af3_structures(
            args.zip_path,
            snapshot_dir=args.snapshot_dir,
            output_dir=args.output_dir,
            verbose=args.verbose,
        )
        ok = sum(1 for row in result.rows if row.get("status") == "ok")
        errors = sum(1 for row in result.rows if row.get("status") != "ok")
        print(result.manifest_tsv)
        print(result.manifest_json)
        print(f"AF3 targets prepared: {ok}; errors: {errors}")
        return
    if args.command == "catalogue-af3-structures":
        result = catalogue_af3_structures(
            args.input_path,
            output_dir=args.output_dir,
            verbose=args.verbose,
        )
        ok = sum(1 for row in result.rows if row.get("status") == "ok")
        selected = sum(1 for row in result.rows if row.get("selected_for_accession"))
        errors = sum(1 for row in result.rows if row.get("status") != "ok")
        print(result.db_path)
        print(result.manifest_tsv)
        print(result.manifest_json)
        print(f"AF3 targets catalogued: {ok}; selected artifacts: {selected}; errors: {errors}")
        return
    if args.command == "rebuild-batch-summary":
        results = rebuild_batch_summary_from_tsv(
            tsv_path=args.tsv_path,
            output_dir=args.output_dir,
            summary_path=args.summary_path,
        )
        summary_path = Path(args.summary_path) if args.summary_path is not None else Path(args.output_dir) / "batch_summary.json"
        print(summary_path)
        print(f"Tracked {len(results)} rows")
        return
    if args.command == "build-accessible-targets":
        rows = build_accessible_target_tsv(
            source_csv=args.source_csv,
            output_tsv=args.output_tsv,
            allowed_buckets=tuple(args.buckets),
        )
        print(Path(args.output_tsv))
        print(f"Targets written: {len(rows)}")
        return
    if args.command == "mouse-from-human-orthologs":
        result = build_mouse_portal_from_human_orthologs(
            args.summary_path,
            output_dir=args.output_dir,
            build_site=not args.no_build_site,
            bundle_vendor_assets=args.bundle_vendor_assets,
            fetch_mouse_alphafold=args.fetch_mouse_alphafold,
            af3_catalog=args.af3_catalog,
            jobs=args.jobs,
            resume=not args.no_resume,
            verbose=args.verbose,
        )
        print(result.summary_path)
        if result.portal_index_path is not None:
            print(result.portal_index_path)
        print(f"Mouse targets included: {result.included_count}")
        print(f"Human targets skipped: {result.skipped_count}")
        print(result.skipped_manifest_path)
        return
    if args.command == "refresh-accessible-portal":
        result = run_accessible_portal_refresh_job(
            source_csv=args.source_csv,
            output_tsv=args.output_tsv,
            output_dir=args.output_dir,
            state_path=args.state_path or str(Path(args.output_dir) / "accessible_portal_refresh_state.json"),
            allowed_buckets=tuple(args.buckets),
            portal_interval=args.portal_interval,
            summary_interval=args.summary_interval,
            max_targets=args.max_targets,
            max_attempts=args.max_attempts,
            verbose=args.verbose,
            enable_complex_portal=args.complex_portal,
        )
        print(result.state_path)
        print(result.summary_path)
        print(result.target_tsv_path)
        print(
            f"Processed {result.processed_this_run} target(s) this run; "
            f"completed {result.completed_targets}/{result.total_targets}; "
            f"portal_built={'yes' if result.portal_built else 'no'}"
        )
        return
    if args.command == "build-portal":
        index_path = build_portal(
            args.summary_path,
            args.output_dir,
            refresh_reports=args.refresh_reports,
            refresh_stale_reports=not args.no_refresh_stale_reports,
            verbose=args.verbose,
            enable_complex_portal=args.complex_portal,
            bundle_vendor_assets=args.bundle_vendor_assets,
        )
        print(index_path)
        return
    if args.command == "build-portal-db":
        from .portal_db import build_portal_database

        db_path = build_portal_database(
            args.summary_path,
            args.db,
            assets_dir=args.assets_dir,
            bundle_vendor_assets=args.bundle_vendor_assets,
        )
        print(db_path)
        return
    if args.command == "serve-portal-db":
        from .portal_db import run_db_portal

        run_db_portal(args.db, assets_dir=args.assets_dir, host=args.host, port=args.port)
        return
    if args.command == "build-portal-from-db":
        from .portal_db import build_portal_from_database

        index_path = build_portal_from_database(
            args.db,
            args.output_dir,
            assets_dir=args.assets_dir,
            bundle_vendor_assets=args.bundle_vendor_assets,
        )
        print(index_path)
        return
    if args.command == "build-test-portal":
        result = build_test_portal_subset(
            args.summary_path,
            output_dir=args.output_dir,
            per_kind=args.per_kind,
            build_portal_site=not args.no_build_portal,
            bundle_vendor_assets=args.bundle_vendor_assets,
            verbose=args.verbose,
        )
        print(result.summary_path)
        if result.portal_index_path is not None:
            print(result.portal_index_path)
        print(f"Selected targets: {result.selected_count}")
        for kind, entries in result.selected_by_kind.items():
            missing = result.missing_by_kind.get(kind, 0)
            suffix = f" ({missing} missing)" if missing else ""
            print(f"{kind}: {len(entries)}{suffix}")
        return
    if args.command == "serve-portal":
        from .dynamic_portal import run_dynamic_portal

        run_dynamic_portal(
            args.summary_path,
            host=args.host,
            port=args.port,
            refresh_reports=args.refresh_reports,
            refresh_stale_reports=args.refresh_stale_reports and not args.no_refresh_stale_reports,
            verbose=args.verbose,
            enable_complex_portal=args.complex_portal,
        )
        return
    if args.command == "refresh-report-module":
        summary_path, _, updated = refresh_report_modules(
            args.summary_path,
            modules=args.modules,
            queries=args.queries,
            only_zero_cross_reactivity=args.only_zero_cross_reactivity,
            verbose=args.verbose,
        )
        print(summary_path)
        print(f"Reports updated: {updated}")
        return
    if args.command == "build-open-targets":
        if args.source == "bulk":
            result = build_open_targets_associations_from_bulk_downloads(
                args.summary_path,
                output_dir=args.output_dir,
                release=args.release,
                data_dir=args.bulk_data_dir,
                download=not args.no_download,
                verbose=args.verbose,
            )
        else:
            result = build_open_targets_associations_from_summary(
                args.summary_path,
                output_dir=args.output_dir,
                page_size=args.page_size,
                max_pages=args.max_pages,
                limit=args.limit,
                resume=not args.no_resume,
                verbose=args.verbose,
            )
        print(result.output_tsv)
        print(result.output_json)
        print(f"Targets written: {len(result.targets)}")
        print(f"Associations written: {len(result.associations)}")
        return
    if args.command == "build-ortholog-table":
        analyzer = AntigenAnalyzer(
            config=AnalysisConfig(
                verbose_progress=args.verbose,
                generate_assets=False,
                render_structure_images=False,
                render_quality_plots=False,
                enable_complex_portal=False,
            )
        )
        output_tsv, output_json, rows = build_ortholog_table_from_tsv(
            analyzer=analyzer,
            tsv_path=args.tsv_path,
            output_path=args.output_path,
            resume=not (args.no_resume or args.force),
            verbose=args.verbose,
            jobs=args.jobs,
        )
        print(output_tsv)
        print(output_json)
        print(f"Rows written: {len(rows)}")
        return
    if args.command == "build-family-alignments":
        analyzer = AntigenAnalyzer(
            config=AnalysisConfig(
                verbose_progress=args.verbose,
                generate_assets=False,
                render_structure_images=False,
                render_quality_plots=False,
                enable_complex_portal=False,
            )
        )
        summary_tsv, summary_json, summaries = build_family_alignments_from_tsv(
            analyzer=analyzer,
            tsv_path=args.tsv_path,
            output_dir=args.output_dir,
            verbose=args.verbose,
        )
        print(summary_tsv)
        print(summary_json)
        print(f"Families written: {len(summaries)}")
        return
    if args.command == "build-paralog-reference":
        analyzer = AntigenAnalyzer(
            config=AnalysisConfig(
                verbose_progress=args.verbose,
                generate_assets=False,
                render_structure_images=False,
                render_quality_plots=False,
                enable_complex_portal=False,
            )
        )
        summary_tsv, summary_json, summaries = build_paralog_reference_from_tsv(
            analyzer=analyzer,
            tsv_path=args.tsv_path,
            output_dir=args.output_dir,
            verbose=args.verbose,
        )
        print(summary_tsv)
        print(summary_json)
        print(f"Targets written: {len(summaries)}")
        return
    parser.print_help()


if __name__ == "__main__":
    main()
