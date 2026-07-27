from __future__ import annotations

import csv
from pathlib import Path


DEFAULT_ACCESSIBLE_BUCKETS = ("Secreted", "GPI", "Single-pass", "Multipass")


def build_accessible_target_tsv(
    *,
    source_csv: str | Path,
    output_tsv: str | Path,
    allowed_buckets: tuple[str, ...] = DEFAULT_ACCESSIBLE_BUCKETS,
) -> list[dict[str, str]]:
    source_path = Path(source_csv)
    output_path = Path(output_tsv)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with source_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)

    allowed = {value.strip() for value in allowed_buckets}
    filtered: list[dict[str, str]] = []
    seen: set[str] = set()
    for row in rows:
        bucket = str(row.get("simplified_topology_bucket") or "").strip()
        if bucket not in allowed:
            continue
        entry_name = str(row.get("uniprot_entry_name") or "").strip()
        accession = str(row.get("accession") or "").strip()
        gene = str(row.get("gene_symbol") or "").strip()
        key = entry_name or accession or gene
        if not key or key in seen:
            continue
        seen.add(key)
        filtered.append(
            {
                "uniprot_id": accession,
                "uniprot_name": entry_name,
                "gene_symbol": gene,
                "topology_bucket": bucket,
                "primary_topology_class": str(row.get("primary_topology_class") or "").strip(),
                "protein_name": str(row.get("protein_name") or "").strip(),
            }
        )

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "uniprot_id",
                "uniprot_name",
                "gene_symbol",
                "topology_bucket",
                "primary_topology_class",
                "protein_name",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(filtered)
    return filtered
