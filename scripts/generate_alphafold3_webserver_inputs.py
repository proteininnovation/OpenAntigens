#!/usr/bin/env python3
"""Generate AlphaFold Server JSON batches for OpenAntigens targets missing AFDB PDBs."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path
from typing import Any


STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
SAFE_NAME_RE = re.compile(r"[^A-Za-z0-9_.-]+")


def safe_job_name(gene_symbol: str, accession: str) -> str:
    raw = f"OpenAntigen_AF3_{gene_symbol}_{accession}"
    return SAFE_NAME_RE.sub("_", raw).strip("_")[:80]


def load_reports(reports_dir: Path) -> list[dict[str, Any]]:
    reports: list[dict[str, Any]] = []
    for report_path in sorted(reports_dir.glob("*_report.json")):
        with report_path.open() as handle:
            report = json.load(handle)
        target = report.get("target") or {}
        accession = str(target.get("accession") or "").strip()
        gene_symbol = str(target.get("gene_symbol") or accession or report_path.stem).strip()
        sequence = str(target.get("sequence") or "").replace(" ", "").replace("\n", "").upper()
        if not accession:
            continue
        reports.append(
            {
                "accession": accession,
                "gene_symbol": gene_symbol,
                "entry_name": str(target.get("entry_name") or "").strip(),
                "protein_name": str(target.get("protein_name") or "").strip(),
                "sequence": sequence,
                "sequence_length": len(sequence),
                "report_path": str(report_path.relative_to(reports_dir)),
            }
        )
    return reports


def classify_missing(
    reports: list[dict[str, Any]], alphafold_dir: Path, max_tokens: int
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    ready: list[dict[str, Any]] = []
    skipped: list[dict[str, Any]] = []
    present: list[dict[str, Any]] = []
    seen_accessions: set[str] = set()

    for record in reports:
        accession = record["accession"]
        if accession in seen_accessions:
            skipped.append({**record, "skip_reason": "duplicate_accession"})
            continue
        seen_accessions.add(accession)

        pdb_path = alphafold_dir / f"{accession}.pdb"
        if pdb_path.exists() and pdb_path.stat().st_size > 0:
            present.append(record)
            continue

        sequence = record["sequence"]
        invalid = sorted(set(sequence) - STANDARD_AA)
        if not sequence:
            skipped.append({**record, "skip_reason": "missing_sequence"})
        elif invalid:
            skipped.append({**record, "skip_reason": f"non_standard_amino_acids:{''.join(invalid)}"})
        elif len(sequence) < 4:
            skipped.append({**record, "skip_reason": "sequence_shorter_than_4_aa"})
        elif len(sequence) > max_tokens:
            skipped.append({**record, "skip_reason": f"sequence_longer_than_{max_tokens}_aa"})
        else:
            ready.append(record)

    return ready, skipped, present


def make_server_job(record: dict[str, Any]) -> dict[str, Any]:
    return {
        "name": safe_job_name(record["gene_symbol"], record["accession"]),
        "modelSeeds": [],
        "sequences": [
            {
                "proteinChain": {
                    "sequence": record["sequence"],
                    "count": 1,
                }
            }
        ],
        "dialect": "alphafoldserver",
        "version": 1,
    }


def write_json(path: Path, payload: Any) -> None:
    with path.open("w") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")


def write_manifest(path: Path, rows: list[dict[str, Any]]) -> None:
    fieldnames = [
        "batch_file",
        "batch_index",
        "job_index",
        "job_name",
        "accession",
        "gene_symbol",
        "entry_name",
        "protein_name",
        "sequence_length",
        "report_path",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_skipped(path: Path, rows: list[dict[str, Any]]) -> None:
    fieldnames = [
        "skip_reason",
        "accession",
        "gene_symbol",
        "entry_name",
        "protein_name",
        "sequence_length",
        "report_path",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--snapshot-dir",
        type=Path,
        default=Path("snapshots/openantigen-2026-05-02"),
        help="OpenAntigens snapshot directory.",
    )
    parser.add_argument(
        "--reports-dir",
        type=Path,
        help="Directory containing *_report.json files. Defaults to <snapshot>/outputs/surfy_batch.",
    )
    parser.add_argument(
        "--alphafold-dir",
        type=Path,
        help="Directory containing local AlphaFold PDB files. Defaults to <snapshot>/data/alphafold.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        help="Output directory for AlphaFold Server JSON batches.",
    )
    parser.add_argument("--batch-size", type=int, default=100)
    parser.add_argument("--max-tokens", type=int, default=5000)
    args = parser.parse_args()

    snapshot_dir = args.snapshot_dir
    reports_dir = args.reports_dir or snapshot_dir / "outputs" / "surfy_batch"
    alphafold_dir = args.alphafold_dir or snapshot_dir / "data" / "alphafold"
    out_dir = args.out_dir or snapshot_dir / "alphafold3_webserver_inputs"
    out_dir.mkdir(parents=True, exist_ok=True)

    reports = load_reports(reports_dir)
    ready, skipped, present = classify_missing(reports, alphafold_dir, args.max_tokens)
    ready.sort(key=lambda row: (row["gene_symbol"], row["accession"]))

    manifest_rows: list[dict[str, Any]] = []
    batch_count = math.ceil(len(ready) / args.batch_size) if ready else 0
    for batch_index in range(batch_count):
        start = batch_index * args.batch_size
        batch_records = ready[start : start + args.batch_size]
        jobs = [make_server_job(record) for record in batch_records]
        batch_file = f"openantigen_af3_missing_{batch_index + 1:03d}.json"
        write_json(out_dir / batch_file, jobs)
        for job_index, (record, job) in enumerate(zip(batch_records, jobs, strict=True), start=1):
            manifest_rows.append(
                {
                    "batch_file": batch_file,
                    "batch_index": batch_index + 1,
                    "job_index": job_index,
                    "job_name": job["name"],
                    "accession": record["accession"],
                    "gene_symbol": record["gene_symbol"],
                    "entry_name": record["entry_name"],
                    "protein_name": record["protein_name"],
                    "sequence_length": record["sequence_length"],
                    "report_path": record["report_path"],
                }
            )

    write_manifest(out_dir / "manifest.tsv", manifest_rows)
    write_skipped(out_dir / "skipped.tsv", skipped)
    summary = {
        "report_count": len(reports),
        "alphafold_pdb_present_count": len(present),
        "missing_and_ready_count": len(ready),
        "skipped_count": len(skipped),
        "batch_size": args.batch_size,
        "batch_count": batch_count,
        "max_tokens": args.max_tokens,
        "json_format": "AlphaFold Server top-level list of alphafoldserver jobs",
    }
    write_json(out_dir / "summary.json", summary)
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
