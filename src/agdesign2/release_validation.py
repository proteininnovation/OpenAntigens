from __future__ import annotations

import gzip
import json
from pathlib import Path
from typing import Any


MIN_HUMAN_TARGETS = 5_000
MIN_MOUSE_TARGETS = 4_500
MIN_BLAST_COVERAGE = 0.80
MIN_OPEN_TARGETS_COVERAGE = 0.80
MIN_LOCAL_STRUCTURE_COVERAGE = 0.80
MIN_CATALOG_SIZE_FOR_STRUCTURE_COVERAGE = 1_000


def validate_structure_coverage(portal_dir: str | Path) -> None:
    portal_dir = Path(portal_dir)
    for relative_path in (
        Path("downloads/agdesign2_portal_index.json"),
        Path("mouse/downloads/agdesign2_portal_index.json"),
    ):
        path = portal_dir / relative_path
        if not path.is_file():
            continue
        rows = _load_rows(path)
        if len(rows) < MIN_CATALOG_SIZE_FOR_STRUCTURE_COVERAGE:
            continue
        available = sum(int(row.get("has_alphafold_structure") or 0) for row in rows)
        coverage = available / len(rows)
        if coverage < MIN_LOCAL_STRUCTURE_COVERAGE:
            raise ValueError(
                f"local AlphaFold coverage is {available}/{len(rows)} ({coverage:.1%}) in "
                f"{relative_path}; minimum is {MIN_LOCAL_STRUCTURE_COVERAGE:.0%}"
            )


def validate_release_data(portal_dir: str | Path, *, require_full_catalog: bool = False) -> None:
    portal_dir = Path(portal_dir)
    human = _load_rows(portal_dir / "downloads/agdesign2_portal_index.json")
    full_release = require_full_catalog or len(human) >= 1_000
    _validate_catalog(
        human,
        label="human",
        catalog_root=portal_dir,
        metadata_path=portal_dir / "portal_metadata.json",
        minimum_targets=MIN_HUMAN_TARGETS if require_full_catalog else 0,
        strict=full_release,
    )

    if full_release:
        _require_coverage(human, "cross_reactivity_count", MIN_BLAST_COVERAGE, "human BLAST")
        _require_coverage(human, "top_disease_name", MIN_OPEN_TARGETS_COVERAGE, "Open Targets")
        disease_tsv = portal_dir / "open_targets_disease_associations.tsv"
        if not disease_tsv.is_file() or disease_tsv.stat().st_size == 0:
            raise ValueError("Open Targets association TSV is missing or empty")
        with disease_tsv.open("rb") as handle:
            header = handle.readline().decode("utf-8", errors="replace")
            association_rows = sum(bool(line.strip()) for line in handle)
        if "disease_name" not in header or "open_targets_url" not in header:
            raise ValueError("Open Targets association TSV has an invalid header")
        expected_associations = sum(int(row.get("disease_count") or 0) for row in human)
        if association_rows != expected_associations:
            raise ValueError(
                "Open Targets association TSV row count does not match the portal index: "
                f"{association_rows} != {expected_associations}"
            )

    mouse_path = portal_dir / "mouse/downloads/agdesign2_portal_index.json"
    if require_full_catalog and not mouse_path.is_file():
        raise ValueError("mouse portal index is missing")
    if mouse_path.is_file():
        mouse = _load_rows(mouse_path)
        _validate_catalog(
            mouse,
            label="mouse",
            catalog_root=portal_dir / "mouse",
            metadata_path=portal_dir / "mouse/portal_metadata.json",
            minimum_targets=MIN_MOUSE_TARGETS if require_full_catalog else 0,
            strict=require_full_catalog or len(mouse) >= 1_000,
        )
        if require_full_catalog or len(mouse) >= 1_000:
            _require_coverage(mouse, "cross_reactivity_count", MIN_BLAST_COVERAGE, "mouse BLAST")


def _load_rows(path: Path) -> list[dict[str, Any]]:
    if not path.is_file():
        raise ValueError(f"portal index is missing: {path}")
    data = path.read_bytes()
    if data.startswith(b"\x1f\x8b"):
        data = gzip.decompress(data)
    payload = json.loads(data)
    if not isinstance(payload, list) or not payload:
        raise ValueError(f"portal index is not a non-empty list: {path}")
    return payload


def _validate_catalog(
    rows: list[dict[str, Any]],
    *,
    label: str,
    catalog_root: Path,
    metadata_path: Path,
    minimum_targets: int,
    strict: bool,
) -> None:
    if minimum_targets and len(rows) < minimum_targets:
        raise ValueError(f"{label} catalog has only {len(rows)} targets; minimum is {minimum_targets}")
    if not strict:
        return
    bad_statuses = [row for row in rows if row.get("status") not in {"ok", "skipped_existing"}]
    if bad_statuses:
        raise ValueError(f"{label} catalog contains {len(bad_statuses)} failed or incomplete targets")
    missing_pages = sum(not row.get("detail_page") for row in rows)
    if missing_pages:
        raise ValueError(f"{label} catalog contains {missing_pages} targets without detail pages")
    detail_pages = {Path(str(row["detail_page"])) for row in rows}
    missing_files = [path for path in detail_pages if not (catalog_root / path).is_file()]
    if missing_files:
        raise ValueError(f"{label} catalog is missing {len(missing_files)} rendered detail pages")
    missing_scripts = [
        path
        for path in detail_pages
        if not (catalog_root / "report_scripts" / path.with_suffix(".js").name).is_file()
    ]
    if missing_scripts:
        raise ValueError(f"{label} catalog is missing {len(missing_scripts)} report scripts")
    if not metadata_path.is_file():
        raise ValueError(f"{label} portal metadata is missing: {metadata_path}")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if int(metadata.get("target_count") or 0) != len(rows):
        raise ValueError(
            f"{label} metadata target_count does not match the index: "
            f"{metadata.get('target_count')!r} != {len(rows)}"
        )


def _require_coverage(
    rows: list[dict[str, Any]],
    field: str,
    minimum: float,
    label: str,
) -> None:
    populated = sum(bool(row.get(field)) for row in rows)
    coverage = populated / len(rows)
    if coverage < minimum:
        raise ValueError(
            f"{label} coverage is {populated}/{len(rows)} ({coverage:.1%}); minimum is {minimum:.0%}"
        )
