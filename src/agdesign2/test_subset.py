from __future__ import annotations

import csv
import json
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any


TEST_PORTAL_KIND_ORDER = (
    "multipass",
    "gpcr",
    "single_pass",
    "gpi",
    "secreted",
    "obligatory_partner",
)

TEST_PORTAL_KIND_LABELS = {
    "multipass": "Multipass",
    "gpcr": "GPCR",
    "single_pass": "Single-pass",
    "gpi": "GPI",
    "secreted": "Secreted",
    "obligatory_partner": "Obligatory partner",
}


@dataclass(frozen=True)
class TestPortalSubsetResult:
    summary_path: Path
    manifest_path: Path
    portal_index_path: Path | None
    selected_count: int
    selected_by_kind: dict[str, list[str]]
    missing_by_kind: dict[str, int]


@dataclass
class _Candidate:
    index: int
    row: dict[str, Any]
    report: dict[str, Any] | None
    report_json_path: Path | None
    markdown_path: Path | None
    entry_name: str
    kinds: set[str]


def build_test_portal_subset(
    summary_path: str | Path,
    *,
    output_dir: str | Path = "snapshots/openantigen-test-10",
    per_kind: int = 10,
    build_portal_site: bool = True,
    bundle_vendor_assets: bool = False,
    verbose: bool = False,
) -> TestPortalSubsetResult:
    """Create a small static-portal snapshot with representative target kinds."""
    summary_file = Path(summary_path)
    output_root = Path(output_dir)
    subset_batch_dir = output_root / "data"
    portal_dir = output_root / "public_site"
    subset_batch_dir.mkdir(parents=True, exist_ok=True)

    rows = json.loads(summary_file.read_text(encoding="utf-8"))
    candidates = _load_candidates(rows, summary_file.parent)
    selected, selected_by_kind = _select_candidates(candidates, per_kind=per_kind)
    selected_keys = {candidate.entry_name.upper() for candidate in selected}

    subset_rows = _write_subset_reports(selected, subset_batch_dir=subset_batch_dir)
    summary_out = subset_batch_dir / "batch_summary.json"
    summary_out.write_text(json.dumps(subset_rows, indent=2) + "\n", encoding="utf-8")

    _copy_support_files(summary_file.parent, subset_batch_dir, selected_keys=selected_keys, verbose=verbose)
    manifest = _build_manifest(
        source_summary=summary_file,
        output_root=output_root,
        per_kind=per_kind,
        selected=selected,
        selected_by_kind=selected_by_kind,
        candidates=candidates,
    )
    manifest_path = subset_batch_dir / "test_subset_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    portal_index_path: Path | None = None
    if build_portal_site:
        from .portal import build_portal

        portal_index_path = build_portal(
            summary_out,
            portal_dir,
            refresh_reports=False,
            refresh_stale_reports=False,
            verbose=verbose,
            bundle_vendor_assets=bundle_vendor_assets,
        )
        shutil.copy2(manifest_path, portal_dir / manifest_path.name)

    missing_by_kind = {
        kind: max(0, per_kind - len(selected_by_kind.get(kind, []))) for kind in TEST_PORTAL_KIND_ORDER
    }
    return TestPortalSubsetResult(
        summary_path=summary_out,
        manifest_path=manifest_path,
        portal_index_path=portal_index_path,
        selected_count=len(selected),
        selected_by_kind={kind: [candidate.entry_name for candidate in values] for kind, values in selected_by_kind.items()},
        missing_by_kind=missing_by_kind,
    )


def _load_candidates(rows: list[dict[str, Any]], batch_dir: Path) -> list[_Candidate]:
    candidates: list[_Candidate] = []
    for index, row in enumerate(rows):
        if row.get("status") not in (None, "ok"):
            continue
        report_json_path = _resolve_path(row.get("json_report"), batch_dir)
        report = _load_json(report_json_path)
        entry_name = _entry_name(row, report)
        if not entry_name:
            continue
        markdown_path = _resolve_path(row.get("markdown_report"), batch_dir)
        kinds = _candidate_kinds(row, report)
        if not kinds:
            continue
        candidates.append(
            _Candidate(
                index=index,
                row=row,
                report=report,
                report_json_path=report_json_path,
                markdown_path=markdown_path,
                entry_name=entry_name,
                kinds=kinds,
            )
        )
    return candidates


def _select_candidates(
    candidates: list[_Candidate],
    *,
    per_kind: int,
) -> tuple[list[_Candidate], dict[str, list[_Candidate]]]:
    selected: list[_Candidate] = []
    selected_names: set[str] = set()
    selected_by_kind: dict[str, list[_Candidate]] = {}

    for kind in TEST_PORTAL_KIND_ORDER:
        pool = [candidate for candidate in candidates if kind in candidate.kinds]
        if kind == "multipass":
            non_gpcr = [candidate for candidate in pool if "gpcr" not in candidate.kinds]
            pool = non_gpcr + [candidate for candidate in pool if candidate not in non_gpcr]

        chosen: list[_Candidate] = []
        for candidate in pool:
            key = candidate.entry_name.upper()
            if key in selected_names:
                continue
            chosen.append(candidate)
            if len(chosen) >= per_kind:
                break
        if len(chosen) < per_kind:
            chosen_keys = {candidate.entry_name.upper() for candidate in chosen}
            for candidate in pool:
                key = candidate.entry_name.upper()
                if key in chosen_keys:
                    continue
                chosen.append(candidate)
                if len(chosen) >= per_kind:
                    break

        selected_by_kind[kind] = chosen
        for candidate in chosen:
            key = candidate.entry_name.upper()
            if key not in selected_names:
                selected.append(candidate)
                selected_names.add(key)

    selected.sort(key=lambda candidate: candidate.index)
    return selected, selected_by_kind


def _candidate_kinds(row: dict[str, Any], report: dict[str, Any] | None) -> set[str]:
    kinds: set[str] = set()
    bucket = str(row.get("topology_bucket") or "").lower().replace("_", "-")
    topology = report.get("topology") if isinstance(report, dict) else {}
    topology_class = str((topology or {}).get("topology_class") or row.get("primary_topology_class") or "").lower()

    if "multipass" in bucket or topology_class.startswith("multipass"):
        kinds.add("multipass")
    if "single-pass" in bucket or "single-pass" in topology_class.replace("_", "-"):
        kinds.add("single_pass")
    if "gpi" in bucket or "gpi" in topology_class:
        kinds.add("gpi")
    if "secreted" in bucket or topology_class.startswith("secreted") or "secreted" in topology_class:
        kinds.add("secreted")

    if _is_gpcr(report):
        kinds.add("gpcr")
    if _has_obligatory_partner(report):
        kinds.add("obligatory_partner")
    return kinds


def _is_gpcr(report: dict[str, Any] | None) -> bool:
    if not isinstance(report, dict):
        return False
    if report.get("gpcr_annotation"):
        return True
    canonical_family = report.get("canonical_family") or {}
    family_name = str(canonical_family.get("name") or "").lower()
    if "g protein-coupled receptor" in family_name or "gpcr" in family_name:
        return True
    for annotation in report.get("interpro_annotations") or []:
        text = " ".join(str(annotation.get(key) or "") for key in ("name", "description", "accession")).lower()
        if "g protein-coupled receptor" in text or "gpcr" in text:
            return True
    return False


def _has_obligatory_partner(report: dict[str, Any] | None) -> bool:
    if not isinstance(report, dict):
        return False
    return any(bool(item.get("obligatory")) for item in report.get("assembly_requirements") or [])


def _write_subset_reports(candidates: list[_Candidate], *, subset_batch_dir: Path) -> list[dict[str, Any]]:
    subset_rows: list[dict[str, Any]] = []
    for batch_index, candidate in enumerate(candidates, start=1):
        entry_slug = candidate.entry_name.lower()
        report_name = f"{entry_slug}_report.json"
        markdown_name = f"{entry_slug}_report.md"
        if candidate.report_json_path and candidate.report_json_path.exists():
            shutil.copy2(candidate.report_json_path, subset_batch_dir / report_name)
            _copy_report_asset_dir(candidate.report_json_path.parent, subset_batch_dir, entry_slug)
        if candidate.markdown_path and candidate.markdown_path.exists():
            shutil.copy2(candidate.markdown_path, subset_batch_dir / markdown_name)

        row = dict(candidate.row)
        row["batch_index"] = batch_index
        row["batch_total"] = len(candidates)
        row["json_report"] = report_name
        if (subset_batch_dir / markdown_name).exists():
            row["markdown_report"] = markdown_name
        subset_rows.append(row)
    return subset_rows


def _copy_support_files(source_batch_dir: Path, subset_batch_dir: Path, *, selected_keys: set[str], verbose: bool) -> None:
    copied_names: set[str] = set()
    for source in _support_file_candidates(source_batch_dir):
        if not source.exists() or source.is_dir():
            continue
        if source.name in copied_names:
            continue
        destination = subset_batch_dir / source.name
        if source.name == "open_targets_disease_associations.json":
            _write_filtered_open_targets_json(source, destination, selected_keys=selected_keys)
        elif source.name == "open_targets_disease_associations.tsv":
            _write_filtered_open_targets_tsv(source, destination, selected_keys=selected_keys)
        else:
            shutil.copy2(source, destination)
        copied_names.add(source.name)
        if verbose:
            print(f"Copied support file: {destination}")
    _copy_vendor_assets(source_batch_dir, subset_batch_dir, verbose=verbose)


def _support_file_candidates(source_batch_dir: Path) -> list[Path]:
    bases = [source_batch_dir, source_batch_dir.parent, source_batch_dir.parent / "outputs" / "surfy_batch"]
    candidates: list[Path] = []
    for base in bases:
        for name in (
            "ortholog_reference_table.tsv",
            "family_alignment_index.tsv",
            "paralog_reference_table.tsv",
            "open_targets_disease_associations.json",
            "open_targets_disease_associations.tsv",
        ):
            path = base / name
            if path not in candidates:
                candidates.append(path)
    return candidates


def _copy_vendor_assets(source_batch_dir: Path, subset_batch_dir: Path, *, verbose: bool) -> None:
    destination_dir = subset_batch_dir / "vendor"
    for source in (
        Path.cwd() / "assets" / "vendor" / "3Dmol-min.js",
        source_batch_dir / "vendor" / "3Dmol-min.js",
        source_batch_dir.parent / "vendor" / "3Dmol-min.js",
        source_batch_dir.parent / "public_site" / "3Dmol-min.js",
        source_batch_dir.parent / "outputs" / "surfy_batch" / "portal" / "3Dmol-min.js",
    ):
        if not source.exists() or source.stat().st_size <= 0:
            continue
        destination_dir.mkdir(parents=True, exist_ok=True)
        destination = destination_dir / "3Dmol-min.js"
        shutil.copy2(source, destination)
        if verbose:
            print(f"Copied vendor asset: {destination}")
        return


def _write_filtered_open_targets_json(source: Path, destination: Path, *, selected_keys: set[str]) -> None:
    try:
        payload = json.loads(source.read_text(encoding="utf-8"))
    except Exception:
        return
    targets = [
        row
        for row in payload.get("targets", [])
        if str(row.get("entry_name") or "").upper() in selected_keys
        or str(row.get("gene_symbol") or "").upper() in selected_keys
        or str(row.get("uniprot_accession") or "").upper() in selected_keys
    ]
    target_entries = {str(row.get("entry_name") or "").upper() for row in targets}
    associations = [
        row for row in payload.get("associations", []) if str(row.get("entry_name") or "").upper() in target_entries
    ]
    filtered = dict(payload)
    filtered["target_count"] = len(targets)
    filtered["association_count"] = len(associations)
    filtered["targets"] = targets
    filtered["associations"] = associations
    destination.write_text(json.dumps(filtered, indent=2) + "\n", encoding="utf-8")


def _write_filtered_open_targets_tsv(source: Path, destination: Path, *, selected_keys: set[str]) -> None:
    with source.open("r", encoding="utf-8", newline="") as input_handle, destination.open(
        "w", encoding="utf-8", newline=""
    ) as output_handle:
        reader = csv.DictReader(input_handle, delimiter="\t")
        if reader.fieldnames is None:
            return
        writer = csv.DictWriter(output_handle, fieldnames=reader.fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in reader:
            if str(row.get("entry_name") or "").upper() in selected_keys:
                writer.writerow(row)


def _build_manifest(
    *,
    source_summary: Path,
    output_root: Path,
    per_kind: int,
    selected: list[_Candidate],
    selected_by_kind: dict[str, list[_Candidate]],
    candidates: list[_Candidate],
) -> dict[str, Any]:
    available_by_kind = {
        kind: sum(1 for candidate in candidates if kind in candidate.kinds) for kind in TEST_PORTAL_KIND_ORDER
    }
    return {
        "source_summary": str(source_summary),
        "output_root": str(output_root),
        "per_kind": per_kind,
        "selected_count": len(selected),
        "selected_entries": [candidate.entry_name for candidate in selected],
        "selected_by_kind": {
            kind: [candidate.entry_name for candidate in selected_by_kind.get(kind, [])]
            for kind in TEST_PORTAL_KIND_ORDER
        },
        "available_by_kind": available_by_kind,
        "missing_by_kind": {
            kind: max(0, per_kind - len(selected_by_kind.get(kind, []))) for kind in TEST_PORTAL_KIND_ORDER
        },
        "kind_labels": TEST_PORTAL_KIND_LABELS,
    }


def _copy_report_asset_dir(source_dir: Path, subset_batch_dir: Path, entry_slug: str) -> None:
    asset_dir = source_dir / f"{entry_slug}_report_assets"
    if not asset_dir.is_dir():
        return
    shutil.copytree(asset_dir, subset_batch_dir / asset_dir.name, dirs_exist_ok=True)


def _entry_name(row: dict[str, Any], report: dict[str, Any] | None) -> str:
    if isinstance(report, dict):
        target = report.get("target") or {}
        entry_name = str(target.get("entry_name") or "").strip()
        if entry_name:
            return entry_name
    return str(row.get("resolved_entry_name") or row.get("uniprot_name") or row.get("query") or "").strip()


def _resolve_path(path_text: str | None, batch_dir: Path) -> Path | None:
    if not path_text:
        return None
    path = Path(path_text)
    candidates = [path] if path.is_absolute() else [path, batch_dir / path, batch_dir / path.name, Path.cwd() / path]
    return next((candidate.resolve() for candidate in candidates if candidate.exists()), None)


def _load_json(path: Path | None) -> dict[str, Any] | None:
    if path is None:
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None
