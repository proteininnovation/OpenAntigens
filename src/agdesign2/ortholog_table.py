from __future__ import annotations

import csv
import io
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from .batch import BatchRow, _load_input_records
from .http import _atomic_write
from .ortholog_proteomes import OrthologProteomeResolver
from .pipeline import AntigenAnalyzer
from .sequence_utils import global_align


@dataclass(slots=True)
class OrthologTableRow:
    input_index: int
    input_uniprot_id: str | None
    input_uniprot_name: str | None
    input_gene: str | None
    input_prot_family: str | None
    query: str
    source_column: str
    status: str
    error: str | None
    human_gene_symbol: str
    mouse_gene_symbol: str | None
    macaca_fascicularis_gene_symbol: str | None
    canonical_family_accession: str | None
    canonical_family_name: str | None
    human_refseq_accession: str | None
    mouse_refseq_accession: str | None
    macaca_fascicularis_refseq_accession: str | None
    human_uniprot_accession: str | None
    human_refseq_sequence: str | None
    mouse_refseq_sequence: str | None
    macaca_fascicularis_refseq_sequence: str | None
    mouse_identity_to_human: float | None
    macaca_fascicularis_identity_to_human: float | None
    human_notes: str
    mouse_notes: str
    macaca_fascicularis_notes: str


def build_ortholog_table_from_tsv(
    *,
    analyzer: AntigenAnalyzer,
    tsv_path: str | Path,
    output_path: str | Path,
    resume: bool = True,
    verbose: bool = False,
    jobs: int = 1,
) -> tuple[Path, Path, list[OrthologTableRow]]:
    batch_rows = load_ortholog_rows(tsv_path)
    output_tsv = Path(output_path)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    output_json = output_tsv.with_suffix(".json")
    if not resume:
        _atomic_write(output_json, b"[]\n")
        for checkpoint in output_json.with_suffix(".checkpoints").glob("*.json"):
            checkpoint.unlink()
    existing_rows = _load_existing_rows(output_json=output_json, batch_rows=batch_rows) if resume else {}
    # Build (download + index) the mouse UniProt and macaque RefSeq proteomes once,
    # serially, so worker threads only do lock-free local lookups — no per-target
    # NCBI calls. Skip the work entirely when every row is already resumed.
    resolver = OrthologProteomeResolver(
        http=analyzer.refseq_client.http,
        data_dir=analyzer.config.ortholog_fasta_dir,
        species_taxonomy=analyzer.config.species_taxonomy,
        verbose=verbose,
    )
    if len(existing_rows) < len(batch_rows):
        resolver.prewarm(["mouse", "macaca_fascicularis"])
    if max(1, jobs) <= 1:
        return _build_ortholog_table_sequential(
            analyzer=analyzer,
            resolver=resolver,
            batch_rows=batch_rows,
            existing_rows=existing_rows,
            output_tsv=output_tsv,
            output_json=output_json,
            verbose=verbose,
        )
    return _build_ortholog_table_parallel(
        analyzer=analyzer,
        resolver=resolver,
        batch_rows=batch_rows,
        existing_rows=existing_rows,
        output_tsv=output_tsv,
        output_json=output_json,
        verbose=verbose,
        jobs=max(1, jobs),
    )


def _build_ortholog_table_sequential(
    *,
    analyzer: AntigenAnalyzer,
    resolver: OrthologProteomeResolver,
    batch_rows: list[BatchRow],
    existing_rows: dict[int, OrthologTableRow],
    output_tsv: Path,
    output_json: Path,
    verbose: bool,
) -> tuple[Path, Path, list[OrthologTableRow]]:
    table_rows: list[OrthologTableRow] = []
    total = len(batch_rows)
    for index, batch_row in enumerate(batch_rows, start=1):
        existing_row = existing_rows.get(index)
        if existing_row is not None:
            table_rows.append(existing_row)
            if verbose:
                print(
                    f"[agdesign2-orthologs] {index}/{total} {batch_row.query} (resume {existing_row.status})",
                    flush=True,
                )
            continue
        if verbose:
            print(f"[agdesign2-orthologs] {index}/{total} {batch_row.query}", flush=True)
        table_rows.append(_build_table_row(analyzer=analyzer, resolver=resolver, batch_row=batch_row, input_index=index))
        _write_row_checkpoint(output_json, table_rows[-1])
    _write_outputs(output_tsv=output_tsv, output_json=output_json, rows=table_rows)
    return output_tsv, output_json, table_rows


def _build_ortholog_table_parallel(
    *,
    analyzer: AntigenAnalyzer,
    resolver: OrthologProteomeResolver,
    batch_rows: list[BatchRow],
    existing_rows: dict[int, OrthologTableRow],
    output_tsv: Path,
    output_json: Path,
    verbose: bool,
    jobs: int,
) -> tuple[Path, Path, list[OrthologTableRow]]:
    # Results are placed by input index, so the final table is byte-identical to
    # the sequential build regardless of completion order. Per-row work is
    # independent network I/O; the shared HTTP cache is written atomically
    # (http._atomic_write) so concurrent threads cannot corrupt a cache file.
    total = len(batch_rows)
    results: list[OrthologTableRow | None] = [None] * total
    pending: list[tuple[int, BatchRow]] = []
    for index, batch_row in enumerate(batch_rows, start=1):
        existing_row = existing_rows.get(index)
        if existing_row is not None:
            results[index - 1] = existing_row
        else:
            pending.append((index, batch_row))

    if pending:
        with ThreadPoolExecutor(max_workers=min(jobs, len(pending))) as executor:
            futures = {
                executor.submit(_build_table_row, analyzer=analyzer, resolver=resolver, batch_row=batch_row, input_index=index): (index, batch_row)
                for index, batch_row in pending
            }
            for future in as_completed(futures):
                index, batch_row = futures[future]
                results[index - 1] = future.result()
                if verbose:
                    print(f"[agdesign2-orthologs] {index}/{total} {batch_row.query}", flush=True)
                _write_row_checkpoint(output_json, results[index - 1])

    table_rows = [r for r in results if r is not None]
    _write_outputs(output_tsv=output_tsv, output_json=output_json, rows=table_rows)
    return output_tsv, output_json, table_rows


def _build_table_row(
    *, analyzer: AntigenAnalyzer, resolver: OrthologProteomeResolver, batch_row: BatchRow, input_index: int
) -> OrthologTableRow:
    query = batch_row.query
    record = batch_row.record
    row_data = _base_row_data(batch_row=batch_row, input_index=input_index)
    row_data["status"] = "ok"
    row_data["error"] = None

    human_gene_symbol = str(record.get("gene") or query)
    canonical_family = None

    ortholog_lookup = None
    human_notes: list[str] = []
    mouse_notes: list[str] = []
    macaca_notes: list[str] = []
    try:
        resolution = analyzer.uniprot_client.resolve_target(query)
        human_target = resolution.target
        human_gene_symbol = human_target.gene_symbol or human_gene_symbol
        human_notes.extend(note.message for note in resolution.notes)
    except Exception as exc:
        row_data["status"] = "error"
        row_data["error"] = str(exc)
        human_notes.append(f"Target resolution failed: {exc}")
        row_data["human_gene_symbol"] = human_gene_symbol
        row_data["human_notes"] = " ".join(human_notes)
        row_data["mouse_notes"] = " ".join(mouse_notes)
        row_data["macaca_fascicularis_notes"] = " ".join(macaca_notes)
        return OrthologTableRow(**row_data)

    try:
        ortholog_lookup = analyzer.hgnc_client.fetch_hcop_orthologs(human_gene_symbol)
    except Exception as exc:
        human_notes.append(f"HCOP ortholog lookup failed: {exc}")

    try:
        canonical_family = analyzer._select_canonical_family(
            analyzer.interpro_client.fetch_annotations(human_target.accession),
            sequence_length=len(human_target.sequence),
            protein_name=human_target.protein_name,
            gene_symbol=human_target.gene_symbol,
        )
    except Exception as exc:
        human_notes.append(f"Canonical family lookup failed: {exc}")

    # Human reference protein: the resolved canonical UniProt (SwissProt) target
    # itself — no extra lookup needed.
    human_sequence = human_target.sequence or None
    human_protein_accession = human_target.accession
    if not human_sequence:
        human_notes.append("No canonical human protein sequence resolved.")

    mouse_prediction = ortholog_lookup.predictions.get("mouse") if ortholog_lookup else None
    if mouse_prediction is None:
        mouse_notes.append("No HCOP mouse ortholog prediction resolved.")
        mouse_protein = None
    else:
        mouse_notes.append(
            f"HCOP ortholog symbol {mouse_prediction.gene_symbol or 'n/a'} from {mouse_prediction.gene_symbol_source or 'unknown'} with {mouse_prediction.evidence_count} sources."
        )
        mouse_protein = resolver.resolve(
            "mouse",
            gene_symbol=mouse_prediction.gene_symbol,
            gene_id=mouse_prediction.ncbi_gene_id,
            fallback_symbol=human_gene_symbol,
        )
        if mouse_protein is not None:
            mouse_notes.append(f"Resolved from {_source_label(mouse_protein.source)} {mouse_protein.accession}.")
        else:
            mouse_notes.append("No mouse UniProt protein resolved (reviewed or reference proteome).")

    macaca_prediction = ortholog_lookup.predictions.get("macaca_fascicularis") if ortholog_lookup else None
    if macaca_prediction is None:
        macaca_notes.append("No HCOP cynomolgus monkey ortholog prediction resolved.")
        macaca_protein = None
    else:
        macaca_notes.append(
            f"HCOP ortholog symbol {macaca_prediction.gene_symbol or 'n/a'} from {macaca_prediction.gene_symbol_source or 'unknown'} with {macaca_prediction.evidence_count} sources."
        )
        macaca_protein = resolver.resolve(
            "macaca_fascicularis",
            gene_symbol=macaca_prediction.gene_symbol,
            gene_id=macaca_prediction.ncbi_gene_id,
            fallback_symbol=human_gene_symbol,
        )
        if macaca_protein is not None:
            macaca_notes.append(f"Resolved from {_source_label(macaca_protein.source)} {macaca_protein.accession}.")
        else:
            macaca_notes.append("No cynomolgus monkey RefSeq protein resolved.")

    mouse_sequence = mouse_protein.sequence if mouse_protein is not None else None
    macaca_sequence = macaca_protein.sequence if macaca_protein is not None else None

    row_data.update(
        {
            "human_gene_symbol": human_gene_symbol,
            "mouse_gene_symbol": mouse_prediction.gene_symbol if mouse_prediction is not None else None,
            "macaca_fascicularis_gene_symbol": macaca_prediction.gene_symbol if macaca_prediction is not None else None,
            "canonical_family_accession": canonical_family.accession if canonical_family is not None else None,
            "canonical_family_name": canonical_family.name if canonical_family is not None else None,
            "human_refseq_accession": human_protein_accession,
            "mouse_refseq_accession": mouse_protein.accession if mouse_protein is not None else None,
            "macaca_fascicularis_refseq_accession": macaca_protein.accession if macaca_protein is not None else None,
            "human_uniprot_accession": human_target.accession,
            "human_refseq_sequence": human_sequence,
            "mouse_refseq_sequence": mouse_sequence,
            "macaca_fascicularis_refseq_sequence": macaca_sequence,
            "mouse_identity_to_human": _sequence_identity(human_sequence, mouse_sequence),
            "macaca_fascicularis_identity_to_human": _sequence_identity(human_sequence, macaca_sequence),
            "human_notes": " ".join(human_notes),
            "mouse_notes": " ".join(mouse_notes),
            "macaca_fascicularis_notes": " ".join(macaca_notes),
        }
    )
    return OrthologTableRow(**row_data)


_SOURCE_LABELS = {
    "uniprot": "UniProt SwissProt (reviewed)",
    "uniprot-trembl": "UniProt TrEMBL (unreviewed)",
    "refseq": "RefSeq",
}


def _source_label(source: str) -> str:
    return _SOURCE_LABELS.get(source, source)


def _sequence_identity(query_sequence: str | None, subject_sequence: str | None) -> float | None:
    if not query_sequence or not subject_sequence:
        return None
    return round(global_align(query_sequence, subject_sequence).identity, 2)


def _write_row_checkpoint(output_json: Path, row: OrthologTableRow) -> None:
    path = output_json.with_suffix(".checkpoints") / f"{row.input_index}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    _atomic_write(path, (json.dumps(asdict(row)) + "\n").encode("utf-8"))


def _write_outputs(*, output_tsv: Path, output_json: Path, rows: list[OrthologTableRow]) -> None:
    row_dicts = [asdict(row) for row in rows]
    fieldnames = list(OrthologTableRow.__dataclass_fields__)
    handle = io.StringIO(newline="")
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(row_dicts)
    _atomic_write(output_tsv, handle.getvalue().encode("utf-8"))
    _atomic_write(output_json, (json.dumps(row_dicts, indent=2) + "\n").encode("utf-8"))
    for checkpoint in output_json.with_suffix(".checkpoints").glob("*.json"):
        checkpoint.unlink()


def _load_existing_rows(
    *,
    output_json: Path,
    batch_rows: list[BatchRow],
) -> dict[int, OrthologTableRow]:
    payload = json.loads(output_json.read_text(encoding="utf-8")) if output_json.exists() else []
    if not isinstance(payload, list):
        raise ValueError(f"Expected an ortholog row list in {output_json}")
    for checkpoint in sorted(output_json.with_suffix(".checkpoints").glob("*.json")):
        payload.append(json.loads(checkpoint.read_text(encoding="utf-8")))
    existing_rows: dict[int, OrthologTableRow] = {}
    for fallback_index, item in enumerate(payload, start=1):
        if not isinstance(item, dict):
            continue
        input_index = _coerce_int(item.get("input_index")) or fallback_index
        if not (1 <= input_index <= len(batch_rows)):
            continue
        batch_row = batch_rows[input_index - 1]
        if not _existing_row_matches_batch(item=item, batch_row=batch_row):
            continue
        normalized = _normalize_existing_row(item=item, batch_row=batch_row, input_index=input_index)
        if normalized.status not in {"ok", "error"}:
            continue
        existing_rows[input_index] = normalized
    return existing_rows


def _normalize_existing_row(*, item: dict[str, Any], batch_row: BatchRow, input_index: int) -> OrthologTableRow:
    data = _base_row_data(batch_row=batch_row, input_index=input_index)
    for field_name in OrthologTableRow.__dataclass_fields__:
        if field_name in {
            "input_index",
            "input_uniprot_id",
            "input_uniprot_name",
            "input_gene",
            "input_prot_family",
            "query",
            "source_column",
        }:
            continue
        if field_name in item:
            data[field_name] = item.get(field_name)
    if data["status"] in {None, "", "pending"}:
        data["status"] = "ok"
    if data["human_gene_symbol"] is None:
        data["human_gene_symbol"] = str(batch_row.record.get("gene") or batch_row.query)
    for field_name in ("human_notes", "mouse_notes", "macaca_fascicularis_notes"):
        if data[field_name] is None:
            data[field_name] = ""
    return OrthologTableRow(**data)


def _base_row_data(*, batch_row: BatchRow, input_index: int) -> dict[str, Any]:
    record = batch_row.record
    return {
        "input_index": input_index,
        "input_uniprot_id": record.get("uniprot_id"),
        "input_uniprot_name": record.get("uniprot_name"),
        "input_gene": record.get("gene"),
        "input_prot_family": record.get("prot_family"),
        "query": batch_row.query,
        "source_column": batch_row.source_column,
        "status": "pending",
        "error": None,
        "human_gene_symbol": str(record.get("gene") or batch_row.query),
        "mouse_gene_symbol": None,
        "macaca_fascicularis_gene_symbol": None,
        "canonical_family_accession": None,
        "canonical_family_name": None,
        "human_refseq_accession": None,
        "mouse_refseq_accession": None,
        "macaca_fascicularis_refseq_accession": None,
        "human_uniprot_accession": None,
        "human_refseq_sequence": None,
        "mouse_refseq_sequence": None,
        "macaca_fascicularis_refseq_sequence": None,
        "mouse_identity_to_human": None,
        "macaca_fascicularis_identity_to_human": None,
        "human_notes": "",
        "mouse_notes": "",
        "macaca_fascicularis_notes": "",
    }


def _coerce_int(value: Any) -> int | None:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _existing_row_matches_batch(*, item: dict[str, Any], batch_row: BatchRow) -> bool:
    record = batch_row.record
    comparable_pairs = [
        (item.get("input_uniprot_name"), record.get("uniprot_name")),
        (item.get("query"), batch_row.query),
        (item.get("input_gene"), record.get("gene")),
        (item.get("input_uniprot_id"), record.get("uniprot_id")),
        (item.get("human_gene_symbol"), record.get("gene")),
    ]
    matched = False
    for left, right in comparable_pairs:
        left_text = str(left or "").strip()
        right_text = str(right or "").strip()
        if not left_text or not right_text:
            continue
        if left_text != right_text:
            return False
        matched = True
    return matched


def load_ortholog_rows(tsv_path: str | Path) -> list[BatchRow]:
    path = Path(tsv_path)
    rows: list[BatchRow] = []
    for record in _load_input_records(path):
        query, source_column = _select_ortholog_query(record)
        if query is None or source_column is None:
            continue
        rows.append(BatchRow(query=query, source_column=source_column, record=record))
    return rows


def _select_ortholog_query(record: dict[str, str]) -> tuple[str | None, str | None]:
    for column in ("uniprot_name", "gene", "uniprot_id"):
        value = (record.get(column) or "").strip()
        if value:
            return value, column
    return None, None
