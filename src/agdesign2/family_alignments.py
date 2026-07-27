from __future__ import annotations

import csv
import json
import re
from dataclasses import asdict, dataclass
from pathlib import Path

from .ortholog_table import load_ortholog_rows
from .pipeline import AntigenAnalyzer
from .sequence_utils import global_align, slice_sequence
from .topology import derive_ectodomain


@dataclass(slots=True)
class FamilyAlignmentMember:
    input_index: int
    input_uniprot_id: str | None
    input_uniprot_name: str | None
    input_gene: str | None
    input_prot_family: str | None
    accession: str
    entry_name: str
    gene_symbol: str | None
    protein_name: str
    family_accession: str
    family_name: str
    ectodomain_start: int
    ectodomain_end: int
    ectodomain_sequence: str
    matrix_label: str


@dataclass(slots=True)
class PairwiseAlignmentRecord:
    query_label: str
    subject_label: str
    query_accession: str
    subject_accession: str
    matches: int
    aligned_positions: int
    query_to_subject_identity: float | None
    subject_to_query_identity: float | None
    aligned_query: str
    aligned_subject: str


def build_family_alignments_from_tsv(
    *,
    analyzer: AntigenAnalyzer,
    tsv_path: str | Path,
    output_dir: str | Path,
    verbose: bool = False,
) -> tuple[Path, Path, list[dict[str, str | int]]]:
    batch_rows = load_ortholog_rows(tsv_path)
    output_root = Path(output_dir)
    output_root.mkdir(parents=True, exist_ok=True)
    family_dir = output_root / "families"
    family_dir.mkdir(parents=True, exist_ok=True)
    summary_tsv = output_root / "family_alignment_index.tsv"
    summary_json = output_root / "family_alignment_index.json"

    grouped: dict[str, list[FamilyAlignmentMember]] = {}
    family_names: dict[str, str] = {}

    total = len(batch_rows)
    for index, batch_row in enumerate(batch_rows, start=1):
        if verbose:
            print(f"[agdesign2-families] {index}/{total} {batch_row.query}", flush=True)
        member = _build_family_member(analyzer=analyzer, batch_row=batch_row, input_index=index)
        if member is None:
            continue
        grouped.setdefault(member.family_accession, []).append(member)
        family_names[member.family_accession] = member.family_name

    summaries: list[dict[str, str | int]] = []
    for family_accession in sorted(grouped):
        members = grouped[family_accession]
        family_name = family_names[family_accession]
        _assign_unique_labels(members)
        payload = _build_family_payload(
            family_accession=family_accession,
            family_name=family_name,
            members=members,
        )
        family_stem = _safe_filename(family_accession.lower())
        family_json = family_dir / f"{family_stem}.json"
        family_fasta = family_dir / f"{family_stem}.fasta"
        family_json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        family_fasta.write_text(_render_family_fasta(members), encoding="utf-8")
        summaries.append(
            {
                "family_accession": family_accession,
                "family_name": family_name,
                "member_count": len(members),
                "json_path": str(family_json),
                "fasta_path": str(family_fasta),
            }
        )

    _write_summary(summary_tsv=summary_tsv, summary_json=summary_json, summaries=summaries)
    return summary_tsv, summary_json, summaries


def _build_family_member(
    *,
    analyzer: AntigenAnalyzer,
    batch_row,
    input_index: int,
) -> FamilyAlignmentMember | None:
    try:
        resolution = analyzer.uniprot_client.resolve_target(batch_row.query)
    except Exception:
        return None
    target = resolution.target
    entry = resolution.entry
    features = analyzer.uniprot_client.get_features(entry)
    topology = derive_ectodomain(features, len(target.sequence), analyzer.config)
    # Mirror the report pipeline: curated-secreted, no-transmembrane targets
    # (e.g. IL1B, FGF2, immunoglobulins) derive their design region from the
    # mature chain/peptide/propeptide fallback rather than topology alone.
    # Without this they get no ectodomain here and are silently dropped from
    # their own family alignment, so they vanish from their own family context.
    analyzer._apply_secreted_universe_fallback(target=target, topology=topology, features=features)
    ectodomain = topology.ectodomain
    if ectodomain is None:
        return None
    try:
        annotations = analyzer.interpro_client.fetch_annotations(target.accession)
    except Exception:
        return None
    canonical_family = analyzer._select_canonical_family(
        annotations,
        sequence_length=len(target.sequence),
        protein_name=target.protein_name,
        gene_symbol=target.gene_symbol,
    )
    if canonical_family is None or not canonical_family.accession.startswith("IPR"):
        return None
    ectodomain_sequence = slice_sequence(target.sequence, ectodomain.start, ectodomain.end)
    if not ectodomain_sequence:
        return None
    return FamilyAlignmentMember(
        input_index=input_index,
        input_uniprot_id=batch_row.record.get("uniprot_id"),
        input_uniprot_name=batch_row.record.get("uniprot_name"),
        input_gene=batch_row.record.get("gene"),
        input_prot_family=batch_row.record.get("prot_family"),
        accession=target.accession,
        entry_name=target.entry_name,
        gene_symbol=target.gene_symbol,
        protein_name=target.protein_name,
        family_accession=canonical_family.accession,
        family_name=canonical_family.name,
        ectodomain_start=ectodomain.start,
        ectodomain_end=ectodomain.end,
        ectodomain_sequence=ectodomain_sequence,
        matrix_label=target.entry_name,
    )


def _assign_unique_labels(members: list[FamilyAlignmentMember]) -> None:
    counts: dict[str, int] = {}
    for member in members:
        label = member.entry_name or member.matrix_label
        counts[label] = counts.get(label, 0) + 1
    duplicate_offsets: dict[str, int] = {}
    for member in members:
        label = member.entry_name or member.matrix_label
        if counts[label] > 1:
            duplicate_offsets[label] = duplicate_offsets.get(label, 0) + 1
            suffix = (
                getattr(member, "input_index", None)
                or getattr(member, "accession", None)
                or getattr(member, "gene_symbol", None)
                or duplicate_offsets[label]
            )
            member.matrix_label = f"{label}__{suffix}"


def _build_family_payload(
    *,
    family_accession: str,
    family_name: str,
    members: list[FamilyAlignmentMember],
) -> dict[str, object]:
    labels = [member.matrix_label for member in members]
    matrix = _build_directional_identity_matrix(members)
    pairwise_alignments = _build_pairwise_alignments(members)
    return {
        "family_accession": family_accession,
        "family_name": family_name,
        "member_count": len(members),
        "identity_matrix_note": "Directional matrix: each cell reports row protein -> column protein identity, normalized by the row ectodomain length.",
        "identity_matrix_labels": labels,
        "identity_matrix": matrix,
        "members": [asdict(member) for member in members],
        "pairwise_alignments": [asdict(record) for record in pairwise_alignments],
    }


def _build_directional_identity_matrix(members: list[FamilyAlignmentMember]) -> list[list[float | None]]:
    matrix: list[list[float | None]] = []
    for row_member in members:
        row: list[float | None] = []
        for col_member in members:
            if row_member.matrix_label == col_member.matrix_label:
                row.append(100.0)
                continue
            alignment = global_align(row_member.ectodomain_sequence, col_member.ectodomain_sequence)
            if not row_member.ectodomain_sequence:
                row.append(None)
                continue
            row.append(round(100.0 * alignment.matches / len(row_member.ectodomain_sequence), 2))
        matrix.append(row)
    return matrix


def _build_pairwise_alignments(members: list[FamilyAlignmentMember]) -> list[PairwiseAlignmentRecord]:
    records: list[PairwiseAlignmentRecord] = []
    for index, query_member in enumerate(members):
        for subject_member in members[index + 1 :]:
            alignment = global_align(query_member.ectodomain_sequence, subject_member.ectodomain_sequence)
            query_identity = (
                round(100.0 * alignment.matches / len(query_member.ectodomain_sequence), 2)
                if query_member.ectodomain_sequence
                else None
            )
            subject_identity = (
                round(100.0 * alignment.matches / len(subject_member.ectodomain_sequence), 2)
                if subject_member.ectodomain_sequence
                else None
            )
            records.append(
                PairwiseAlignmentRecord(
                    query_label=query_member.matrix_label,
                    subject_label=subject_member.matrix_label,
                    query_accession=query_member.accession,
                    subject_accession=subject_member.accession,
                    matches=alignment.matches,
                    aligned_positions=alignment.aligned_positions,
                    query_to_subject_identity=query_identity,
                    subject_to_query_identity=subject_identity,
                    aligned_query=alignment.aligned_query,
                    aligned_subject=alignment.aligned_subject,
                )
            )
    return records


def _render_family_fasta(members: list[FamilyAlignmentMember]) -> str:
    lines: list[str] = []
    for member in members:
        header = (
            f">{member.matrix_label} accession={member.accession} gene={member.gene_symbol or ''} "
            f"family={member.family_accession} ectodomain={member.ectodomain_start}-{member.ectodomain_end}"
        ).strip()
        lines.append(header)
        for start in range(0, len(member.ectodomain_sequence), 80):
            lines.append(member.ectodomain_sequence[start : start + 80])
    return "\n".join(lines) + ("\n" if lines else "")


def _write_summary(
    *,
    summary_tsv: Path,
    summary_json: Path,
    summaries: list[dict[str, str | int]],
) -> None:
    fieldnames = ["family_accession", "family_name", "member_count", "json_path", "fasta_path"]
    with summary_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(summaries)
    summary_json.write_text(json.dumps(summaries, indent=2) + "\n", encoding="utf-8")


def _safe_filename(value: str) -> str:
    return re.sub(r"[^a-zA-Z0-9._-]+", "_", value).strip("_") or "family"
