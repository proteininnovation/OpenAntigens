from __future__ import annotations

import json
from pathlib import Path

from agdesign2.test_subset import build_test_portal_subset


def test_build_test_portal_subset_selects_requested_kinds(tmp_path: Path) -> None:
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    rows = []
    specs = [
        ("ABC1_HUMAN", "Multipass", "multipass_surface", None, False),
        ("GPCR1_HUMAN", "Multipass", "multipass_compact", "G protein-coupled receptor, rhodopsin-like", False),
        ("SP1_HUMAN", "Single-pass", "single_pass_surface", None, False),
        ("GPI1_HUMAN", "GPI", "gpi_anchored_surface", None, False),
        ("SEC1_HUMAN", "Secreted", "secreted_soluble", None, False),
        ("OB1_HUMAN", "Single-pass", "single_pass_surface", None, True),
    ]
    for index, (entry, bucket, topology_class, family_name, obligatory) in enumerate(specs, start=1):
        slug = entry.lower()
        report_path = source_dir / f"{slug}_report.json"
        markdown_path = source_dir / f"{slug}_report.md"
        report = {
            "target": {"entry_name": entry, "accession": f"P{index:05d}", "gene_symbol": entry.split("_", 1)[0]},
            "topology": {"topology_class": topology_class},
            "canonical_family": {"name": family_name} if family_name else None,
            "assembly_requirements": [{"obligatory": obligatory}] if obligatory else [],
        }
        report_path.write_text(json.dumps(report), encoding="utf-8")
        markdown_path.write_text(f"# {entry}\n", encoding="utf-8")
        rows.append(
            {
                "status": "ok",
                "resolved_entry_name": entry,
                "topology_bucket": bucket,
                "primary_topology_class": topology_class,
                "json_report": report_path.name,
                "markdown_report": markdown_path.name,
            }
        )
    summary_path = source_dir / "batch_summary.json"
    summary_path.write_text(json.dumps(rows), encoding="utf-8")

    result = build_test_portal_subset(
        summary_path,
        output_dir=tmp_path / "test-snapshot",
        per_kind=1,
        build_portal_site=False,
    )

    subset_rows = json.loads(result.summary_path.read_text(encoding="utf-8"))
    manifest = json.loads(result.manifest_path.read_text(encoding="utf-8"))
    assert result.selected_count == 6
    assert len(subset_rows) == 6
    assert all((result.summary_path.parent / row["json_report"]).exists() for row in subset_rows)
    assert manifest["selected_by_kind"]["multipass"] == ["ABC1_HUMAN"]
    assert manifest["selected_by_kind"]["gpcr"] == ["GPCR1_HUMAN"]
    assert manifest["selected_by_kind"]["single_pass"] == ["SP1_HUMAN"]
    assert manifest["selected_by_kind"]["gpi"] == ["GPI1_HUMAN"]
    assert manifest["selected_by_kind"]["secreted"] == ["SEC1_HUMAN"]
    assert manifest["selected_by_kind"]["obligatory_partner"] == ["OB1_HUMAN"]
