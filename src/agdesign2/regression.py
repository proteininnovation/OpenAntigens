"""Regression panel: a pinned set of ~120 targets that covers the portal's use cases.

The panel is built with the production command, ``build-fresh-snapshot``, plus ``--target-panel``,
which narrows the target TSV after the normal source-CSV and bucket selection; every later step,
and the publish gate with its catalog minimums scaled to the panel, is the production code path.
``check-regression-panel`` then checks each target's report against the expectations written in
the panel and, given a baseline build, lists what changed.

Expectation tokens (space-separated in the ``expect`` column, added to the category defaults):
  class=A|B        topology class is one of the listed values
  region|no_region a design region is / is not assigned
  constructs>=N    at least N construct recommendations; constructs=0 for none
  region_start=N   design region starts / ends at residue N
  region_end=N
  gpcr             a GPCR annotation is present
  assembly         an assembly requirement (obligate partner) is listed
  mouse            the target has a page in the mouse portal
  note~TEXT        a report note contains TEXT (case-insensitive; use _ for spaces)
"""

from __future__ import annotations

import csv
import gzip
import json
import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


DEFAULT_PANEL_PATH = Path("data/regression/panel.tsv")

CATEGORY_DEFAULTS: dict[str, str] = {
    "gpcr": "class=multipass_compact|multipass_mixed gpcr constructs>=1",
    "single_pass": "class=single_pass region constructs>=1",
    "single_pass_inferred": "class=single_pass region note~inferred",
    "very_large": "region constructs>=1",
    "multipass_compact": "class=multipass_compact no_region constructs>=1",
    "multipass_mixed": "class=multipass_mixed region",
    "gpi": "class=gpi_anchored region constructs>=1",
    "gpi_with_tm": "region constructs>=1",
    "gpi_untrimmed": "class=gpi_anchored_untrimmed region note~no_processed_chain",
    "secreted": "class=secreted region constructs>=1",
    "secreted_propeptide": "class=secreted region constructs>=1",
    "secreted_short": "class=secreted region constructs=0",
    "secreted_fallback": "class=secreted_external_evidence region",
    "secreted_processed_chain": "class=secreted region constructs>=1 note~processed_chain",
    "curated_override": "region constructs>=1 note~curated_design_region",
    "ig_constant": "region",
    "assembly": "assembly",
    "complex_partner": "class=single_pass region",
    "beta_stranded": "class=secreted|non_membrane|single_pass",
    "no_design_region": "no_region constructs=0",
    "no_alphafold": "region",
}

_NUMERIC = re.compile(r"^(constructs)(>=|=)(\d+)$")


@dataclass(frozen=True)
class PanelRow:
    gene_symbol: str
    category: str
    expect: tuple[str, ...]
    why: str


@dataclass
class TargetResult:
    gene_symbol: str
    category: str
    failures: list[str] = field(default_factory=list)
    summary: dict[str, Any] = field(default_factory=dict)


def load_panel(path: str | Path = DEFAULT_PANEL_PATH) -> list[PanelRow]:
    rows: list[PanelRow] = []
    seen: set[str] = set()
    with Path(path).open(encoding="utf-8", newline="") as handle:
        for line_number, row in enumerate(csv.DictReader(handle, delimiter="\t"), start=2):
            gene = (row.get("gene_symbol") or "").strip()
            category = (row.get("category") or "").strip()
            if not gene or gene.startswith("#"):
                continue
            if category not in CATEGORY_DEFAULTS:
                raise ValueError(f"{path}:{line_number}: unknown category {category!r}")
            if gene.upper() in seen:
                raise ValueError(f"{path}:{line_number}: {gene} is listed twice")
            seen.add(gene.upper())
            tokens = tuple(CATEGORY_DEFAULTS[category].split()) + tuple((row.get("expect") or "").split())
            for token in tokens:
                _parse_token(token)  # fail on a malformed panel, not at check time
            rows.append(PanelRow(gene, category, tokens, (row.get("why") or "").strip()))
    if not rows:
        raise ValueError(f"{path} lists no targets")
    return rows


def select_panel_rows(rows: list[dict[str, str]], panel: list[PanelRow]) -> list[dict[str, str]]:
    """Narrow the production target TSV rows to the panel, in panel order. Every panel target must
    be in the production selection: a target the production build would not include cannot test it."""
    by_gene = {str(row.get("gene_symbol") or "").upper(): row for row in rows}
    missing = [item.gene_symbol for item in panel if item.gene_symbol.upper() not in by_gene]
    if missing:
        raise ValueError(
            "regression panel targets missing from the production target selection: " + ", ".join(missing)
        )
    return [by_gene[item.gene_symbol.upper()] for item in panel]


def scaled_mouse_minimum(human_targets: int, *, human_minimum: int, mouse_minimum: int) -> int:
    """Production expects at least mouse_minimum/human_minimum mouse pages per human target."""
    return math.floor(human_targets * mouse_minimum / human_minimum)


def _parse_token(token: str) -> tuple[str, str]:
    if token in {"region", "no_region", "gpcr", "assembly", "mouse"}:
        return token, ""
    match = _NUMERIC.match(token)
    if match:
        return match.group(1) + match.group(2), match.group(3)
    for prefix in ("class=", "region_start=", "region_end=", "note~"):
        if token.startswith(prefix) and len(token) > len(prefix):
            return prefix, token[len(prefix):]
    raise ValueError(f"unknown expectation token {token!r}")


def _load_reports(report_dir: Path, pattern: str) -> dict[str, dict[str, Any]]:
    reports: dict[str, dict[str, Any]] = {}
    for path in sorted(report_dir.glob(pattern)):
        report = json.loads(path.read_text(encoding="utf-8"))
        gene = str(((report.get("target") or {}).get("gene_symbol")) or "").upper()
        if gene:
            reports[gene] = report
    return reports


def _load_index(path: Path) -> list[dict[str, Any]]:
    if not path.is_file():
        return []
    data = path.read_bytes()
    if data.startswith(b"\x1f\x8b"):
        data = gzip.decompress(data)
    return json.loads(data)


def summarize_report(report: dict[str, Any]) -> dict[str, Any]:
    region = report.get("ectodomain") or None
    return {
        "topology_class": (report.get("topology") or {}).get("topology_class"),
        "region": f"{region['start']}-{region['end']}" if region else "",
        "constructs": len(report.get("construct_recommendations") or []),
        "gpcr": bool(report.get("gpcr_annotation")),
        "assembly": len(report.get("assembly_requirements") or []),
    }


def snapshot_summaries(snapshot_dir: str | Path) -> dict[str, dict[str, Any]]:
    """Per human gene symbol: the report summary plus whether a mouse page exists."""
    snapshot = Path(snapshot_dir)
    human = _load_reports(snapshot / "outputs" / "surfy_batch", "*_human_report.json")
    mouse_index = _load_index(snapshot / "public_site" / "mouse" / "downloads" / "agdesign2_portal_index.json")
    with_mouse = {str(row.get("human_source_gene_symbol") or "").upper() for row in mouse_index}
    summaries = {}
    for gene, report in human.items():
        summary = summarize_report(report)
        summary["mouse"] = gene in with_mouse
        summary["notes"] = [str(note.get("message") or "") for note in report.get("notes") or []]
        summary["_report"] = report
        summaries[gene] = summary
    return summaries


def check_target(row: PanelRow, summary: dict[str, Any] | None) -> TargetResult:
    result = TargetResult(row.gene_symbol, row.category)
    if summary is None:
        result.failures.append("no human report in this build")
        return result
    result.summary = {key: value for key, value in summary.items() if not key.startswith("_") and key != "notes"}
    region = (summary["_report"].get("ectodomain") or {})
    for token in row.expect:
        kind, value = _parse_token(token)
        ok = True
        if kind == "class=":
            ok = summary["topology_class"] in value.split("|")
        elif kind == "region":
            ok = bool(summary["region"])
        elif kind == "no_region":
            ok = not summary["region"]
        elif kind == "constructs>=":
            ok = summary["constructs"] >= int(value)
        elif kind == "constructs=":
            ok = summary["constructs"] == int(value)
        elif kind == "region_start=":
            ok = region.get("start") == int(value)
        elif kind == "region_end=":
            ok = region.get("end") == int(value)
        elif kind in {"gpcr", "mouse"}:
            ok = bool(summary[kind])
        elif kind == "assembly":
            ok = summary["assembly"] > 0
        elif kind == "note~":
            needle = value.replace("_", " ").lower()
            ok = any(needle in note.lower() for note in summary["notes"])
        if not ok:
            observed = {
                "class=": summary["topology_class"], "region": summary["region"] or "none",
                "no_region": summary["region"], "constructs>=": summary["constructs"],
                "constructs=": summary["constructs"], "region_start=": region.get("start"),
                "region_end=": region.get("end"),
            }.get(kind, "absent")
            result.failures.append(f"{token} (observed {observed})")
    return result


def diff_against_baseline(
    current: dict[str, dict[str, Any]], baseline: dict[str, dict[str, Any]], genes: list[str]
) -> list[tuple[str, str, Any, Any]]:
    changes = []
    for gene in genes:
        now, before = current.get(gene.upper()), baseline.get(gene.upper())
        if now is None or before is None:
            if (now is None) != (before is None):
                changes.append((gene, "report", "present" if before else "absent", "present" if now else "absent"))
            continue
        for key in ("topology_class", "region", "constructs", "gpcr", "assembly", "mouse"):
            if now.get(key) != before.get(key):
                changes.append((gene, key, before.get(key), now.get(key)))
    return changes


def check_regression_snapshot(
    snapshot_dir: str | Path,
    *,
    panel_path: str | Path = DEFAULT_PANEL_PATH,
    baseline_dir: str | Path | None = None,
) -> tuple[list[TargetResult], list[tuple[str, str, Any, Any]], str]:
    panel = load_panel(panel_path)
    summaries = snapshot_summaries(snapshot_dir)
    results = [check_target(row, summaries.get(row.gene_symbol.upper())) for row in panel]
    changes = (
        diff_against_baseline(summaries, snapshot_summaries(baseline_dir), [row.gene_symbol for row in panel])
        if baseline_dir is not None
        else []
    )
    failed = [result for result in results if result.failures]
    lines = [
        f"# Regression panel check: {Path(snapshot_dir).name}",
        "",
        f"{len(results) - len(failed)} of {len(results)} panel targets meet their expectations.",
        "",
    ]
    if failed:
        lines += ["## Failures", "", "| Target | Category | Failed expectation |", "|---|---|---|"]
        lines += [f"| {r.gene_symbol} | {r.category} | {'; '.join(r.failures)} |" for r in failed]
        lines.append("")
    if baseline_dir is not None:
        lines += [f"## Changes against {Path(baseline_dir).name}", ""]
        if changes:
            lines += ["| Target | Field | Baseline | This build |", "|---|---|---|---|"]
            lines += [f"| {gene} | {key} | {before} | {now} |" for gene, key, before, now in changes]
        else:
            lines.append("No panel target changed topology class, design region, construct count, flags or mouse page.")
        lines.append("")
    by_category: dict[str, list[TargetResult]] = {}
    for result in results:
        by_category.setdefault(result.category, []).append(result)
    lines += ["## Coverage by category", "", "| Category | Targets | Passing |", "|---|---|---|"]
    lines += [
        f"| {category} | {len(items)} | {sum(not item.failures for item in items)} |"
        for category, items in by_category.items()
    ]
    return results, changes, "\n".join(lines) + "\n"
