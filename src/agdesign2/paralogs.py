from __future__ import annotations

import csv
import json
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any
from urllib.parse import quote

from .family_alignments import _assign_unique_labels
from .http import HttpClient
from .ortholog_table import load_ortholog_rows
from .pipeline import AntigenAnalyzer
from .sequence_utils import global_align, slice_sequence
from .topology import derive_ectodomain


@dataclass(slots=True)
class EnsemblParalogPrediction:
    gene_symbol: str | None
    ensembl_gene_id: str | None
    ensembl_protein_id: str | None
    homology_type: str | None
    source_to_target_identity: float | None
    target_to_source_identity: float | None
    taxonomy_level: str | None
    dn: float | None
    ds: float | None
    notes: list[str] = field(default_factory=list)


@dataclass(slots=True)
class ParalogMember:
    gene_symbol: str
    gene_name: str
    accession: str | None
    entry_name: str | None
    sequence_length: int | None
    ectodomain_start: int | None
    ectodomain_end: int | None
    ectodomain_length: int | None
    ectodomain_sequence: str
    matrix_label: str
    sources: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    identity_to_target: float | None = None
    coverage_to_target: float | None = None
    evidence_score: float = 0.0


@dataclass(slots=True)
class ResolvedParalogTarget:
    input_index: int
    batch_row: Any
    base_payload: dict[str, Any]
    target: Any
    canonical_family: Any
    ectodomain_start: int | None
    ectodomain_end: int | None
    ectodomain_sequence: str


class EnsemblParalogClient:
    def __init__(self, http: HttpClient) -> None:
        self.http = http

    def fetch_human_paralogs(self, gene_symbol: str | None) -> list[EnsemblParalogPrediction]:
        if not gene_symbol:
            return []
        url = (
            "https://rest.ensembl.org/homology/symbol/human/"
            f"{quote(gene_symbol)}?content-type=application/json;type=paralogues;format=full;sequence=none;target_species=human"
        )
        payload = self.http.fetch_json(url, cache_namespace="ensembl_paralogs")
        data = payload.get("data") if isinstance(payload, dict) else None
        if not isinstance(data, list) or not data:
            return []
        homologies = data[0].get("homologies") if isinstance(data[0], dict) else None
        if not isinstance(homologies, list):
            return []
        predictions: list[EnsemblParalogPrediction] = []
        pending_lookups: list[tuple[int, str]] = []
        for item in homologies:
            if not isinstance(item, dict):
                continue
            target = item.get("target") if isinstance(item.get("target"), dict) else {}
            ensembl_gene_id = _text(item.get("id")) or _text(target.get("id"))
            gene_symbol_text = (
                _text(target.get("display_id"))
                or _text(target.get("gene_symbol"))
                or _text(target.get("external_id"))
            )
            notes: list[str] = []
            if gene_symbol_text is None and ensembl_gene_id:
                pending_lookups.append((len(predictions), ensembl_gene_id))
            predictions.append(
                EnsemblParalogPrediction(
                    gene_symbol=gene_symbol_text,
                    ensembl_gene_id=ensembl_gene_id,
                    ensembl_protein_id=_text(target.get("protein_id")),
                    homology_type=_text(item.get("type")),
                    source_to_target_identity=_coerce_float(
                        target.get("perc_id") if "perc_id" in target else item.get("source_perc_id")
                    ),
                    target_to_source_identity=_coerce_float(
                        item.get("target_perc_id") if "target_perc_id" in item else target.get("target_perc_id")
                    ),
                    taxonomy_level=_text(item.get("taxonomy_level")),
                    dn=_coerce_float(item.get("dn")),
                    ds=_coerce_float(item.get("ds")),
                    notes=notes,
                )
            )
        if pending_lookups:
            with ThreadPoolExecutor(max_workers=min(4, len(pending_lookups))) as executor:
                futures = [
                    (prediction_index, executor.submit(self.lookup_gene_symbol, ensembl_gene_id))
                    for prediction_index, ensembl_gene_id in pending_lookups
                ]
                for prediction_index, future in futures:
                    gene_symbol_text = future.result()
                    if gene_symbol_text:
                        predictions[prediction_index].gene_symbol = gene_symbol_text
                        predictions[prediction_index].notes.append("Gene symbol resolved through Ensembl lookup.")
        return predictions

    def lookup_gene_symbol(self, ensembl_gene_id: str) -> str | None:
        url = f"https://rest.ensembl.org/lookup/id/{quote(ensembl_gene_id)}?content-type=application/json"
        payload = self.http.fetch_json(url, cache_namespace="ensembl_lookup")
        if not isinstance(payload, dict):
            return None
        return _text(payload.get("display_name")) or _text(payload.get("description"))


def build_paralog_reference_from_tsv(
    *,
    analyzer: AntigenAnalyzer,
    tsv_path: str | Path,
    output_dir: str | Path,
    verbose: bool = False,
    ensembl_client: EnsemblParalogClient | None = None,
) -> tuple[Path, Path, list[dict[str, Any]]]:
    batch_rows = load_ortholog_rows(tsv_path)
    output_root = Path(output_dir)
    output_root.mkdir(parents=True, exist_ok=True)
    paralog_dir = output_root / "targets"
    paralog_dir.mkdir(parents=True, exist_ok=True)
    summary_tsv = output_root / "paralog_index.tsv"
    summary_json = output_root / "paralog_index.json"
    default_http = getattr(analyzer.hgnc_client, "http", None) or getattr(analyzer.uniprot_client, "http", None) or HttpClient()
    ensembl = ensembl_client or EnsemblParalogClient(default_http)
    summaries: list[dict[str, Any]] = []
    payloads_by_index: dict[int, dict[str, Any]] = {}
    family_groups: dict[str, list[ResolvedParalogTarget]] = {}
    seed_cache: dict[str, ResolvedParalogTarget | None] = {}
    member_cache: dict[str, ParalogMember | None] = {}
    family_alignment_cache: dict[str, list[dict[str, Any]]] = {}
    ensembl_cache: dict[str, list[EnsemblParalogPrediction]] = {}
    hgnc_cache: dict[str, Any] = {}
    total = len(batch_rows)

    for index, batch_row in enumerate(batch_rows, start=1):
        if verbose:
            print(f"[agdesign2-paralogs] {index}/{total} {batch_row.query}", flush=True)
        cache_key = str(batch_row.query or "").strip().upper()
        cached_seed = seed_cache.get(cache_key)
        if cached_seed is not None:
            seed = _clone_resolved_seed(cached_seed, batch_row=batch_row, input_index=index)
        else:
            seed = _resolve_paralog_target(analyzer=analyzer, batch_row=batch_row, input_index=index)
            seed_cache[cache_key] = seed
        if seed is None:
            payloads_by_index[index] = _error_payload(batch_row=batch_row, input_index=index, message="Target resolution failed.")
            continue
        if seed.base_payload.get("status") != "ok":
            payloads_by_index[index] = seed.base_payload
            continue
        canonical_family = seed.canonical_family
        if canonical_family is not None and str(canonical_family.accession).startswith("IPR"):
            family_groups.setdefault(str(canonical_family.accession), []).append(seed)
        else:
            payloads_by_index[index] = _build_target_specific_paralog_payload(
                analyzer=analyzer,
                seed=seed,
                ensembl_client=ensembl,
                member_cache=member_cache,
                family_alignment_cache=family_alignment_cache,
                ensembl_cache=ensembl_cache,
                hgnc_cache=hgnc_cache,
            )

    for family_accession, seeds in family_groups.items():
        shared_payload = _build_shared_family_payload(
            analyzer=analyzer,
            seeds=seeds,
            ensembl_client=ensembl,
            member_cache=member_cache,
            family_alignment_cache=family_alignment_cache,
            ensembl_cache=ensembl_cache,
            hgnc_cache=hgnc_cache,
        )
        for seed in seeds:
            payloads_by_index[seed.input_index] = _assemble_target_payload_from_shared(seed=seed, shared_payload=shared_payload)

    for index, batch_row in enumerate(batch_rows, start=1):
        payload = payloads_by_index.get(index) or _error_payload(
            batch_row=batch_row,
            input_index=index,
            message="No paralog payload generated.",
        )
        entry_name = str(payload.get("target_entry_name") or batch_row.query)
        target_json = paralog_dir / f"{_safe_filename(entry_name.lower())}.json"
        target_json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        summaries.append(
            {
                "input_index": index,
                "query": batch_row.query,
                "target_entry_name": payload.get("target_entry_name") or "",
                "target_accession": payload.get("target_accession") or "",
                "target_gene_symbol": payload.get("target_gene_symbol") or "",
                "canonical_family_accession": payload.get("canonical_family_accession") or "",
                "canonical_family_name": payload.get("canonical_family_name") or "",
                "member_count": int(payload.get("member_count") or 0),
                "status": payload.get("status") or "",
                "error": payload.get("error") or "",
                "json_path": str(target_json),
            }
        )
        _write_paralog_summary(summary_tsv=summary_tsv, summary_json=summary_json, summaries=summaries)
    return summary_tsv, summary_json, summaries


def _initialize_base_payload(*, batch_row, input_index: int) -> dict[str, Any]:
    return {
        "input_index": input_index,
        "query": batch_row.query,
        "source_column": batch_row.source_column,
        "input_uniprot_id": batch_row.record.get("uniprot_id"),
        "input_uniprot_name": batch_row.record.get("uniprot_name"),
        "input_gene": batch_row.record.get("gene"),
        "input_prot_family": batch_row.record.get("prot_family"),
        "status": "ok",
        "error": "",
        "target_entry_name": "",
        "target_accession": "",
        "target_gene_symbol": batch_row.record.get("gene") or batch_row.query,
        "canonical_family_accession": "",
        "canonical_family_name": "",
        "family_names": [],
        "identity_matrix_note": "Directional matrix: each cell reports row protein -> column protein identity, normalized by the row ectodomain length.",
        "coverage_matrix_note": "Directional matrix: each cell reports the percent of the row ectodomain aligned to the column ectodomain.",
        "identity_matrix_labels": [],
        "identity_matrix": [],
        "coverage_matrix_labels": [],
        "coverage_matrix": [],
        "members": [],
        "pairwise_alignments": [],
        "metadata": {},
        "member_count": 0,
    }


def _error_payload(*, batch_row, input_index: int, message: str) -> dict[str, Any]:
    payload = _initialize_base_payload(batch_row=batch_row, input_index=input_index)
    payload["status"] = "error"
    payload["error"] = message
    return payload


def _clone_resolved_seed(
    seed: ResolvedParalogTarget,
    *,
    batch_row,
    input_index: int,
) -> ResolvedParalogTarget:
    base_payload = dict(seed.base_payload)
    base_payload.update(
        {
            "input_index": input_index,
            "query": batch_row.query,
            "source_column": batch_row.source_column,
            "input_uniprot_id": batch_row.record.get("uniprot_id"),
            "input_uniprot_name": batch_row.record.get("uniprot_name"),
            "input_gene": batch_row.record.get("gene"),
            "input_prot_family": batch_row.record.get("prot_family"),
        }
    )
    return ResolvedParalogTarget(
        input_index=input_index,
        batch_row=batch_row,
        base_payload=base_payload,
        target=seed.target,
        canonical_family=seed.canonical_family,
        ectodomain_start=seed.ectodomain_start,
        ectodomain_end=seed.ectodomain_end,
        ectodomain_sequence=seed.ectodomain_sequence,
    )


def _resolve_paralog_target(*, analyzer: AntigenAnalyzer, batch_row, input_index: int) -> ResolvedParalogTarget | None:
    base_payload = _initialize_base_payload(batch_row=batch_row, input_index=input_index)
    try:
        resolution = analyzer.uniprot_client.resolve_target(batch_row.query)
        target = resolution.target
        features = analyzer.uniprot_client.get_features(resolution.entry)
        topology = derive_ectodomain(features, len(target.sequence), analyzer.config)
        analyzer._apply_secreted_universe_fallback(target=target, topology=topology, features=features)
        ectodomain = topology.ectodomain
        if ectodomain is None:
            raise ValueError("No ectodomain / designable extracellular region available for paralog analysis.")
        ectodomain_sequence = slice_sequence(target.sequence, ectodomain.start, ectodomain.end)
        if not ectodomain_sequence:
            raise ValueError("Empty ectodomain sequence.")
        base_payload.update(
            {
                "target_entry_name": target.entry_name,
                "target_accession": target.accession,
                "target_gene_symbol": target.gene_symbol or base_payload["target_gene_symbol"],
            }
        )
    except Exception as exc:
        base_payload["status"] = "error"
        base_payload["error"] = str(exc)
        return ResolvedParalogTarget(
            input_index=input_index,
            batch_row=batch_row,
            base_payload=base_payload,
            target=None,
            canonical_family=None,
            ectodomain_start=None,
            ectodomain_end=None,
            ectodomain_sequence="",
        )
    precomputed_ortholog = analyzer.precomputed_store.find_ortholog_record(target=target)
    canonical_family = analyzer.precomputed_store.canonical_family_from_record(precomputed_ortholog)
    if canonical_family is None:
        try:
            annotations = analyzer.interpro_client.fetch_annotations(target.accession)
            canonical_family = analyzer._select_canonical_family(
                annotations,
                sequence_length=len(target.sequence),
                protein_name=target.protein_name,
                gene_symbol=target.gene_symbol,
            )
        except Exception:
            canonical_family = None
    if canonical_family is not None:
        base_payload["canonical_family_accession"] = canonical_family.accession
        base_payload["canonical_family_name"] = canonical_family.name
    return ResolvedParalogTarget(
        input_index=input_index,
        batch_row=batch_row,
        base_payload=base_payload,
        target=target,
        canonical_family=canonical_family,
        ectodomain_start=ectodomain.start,
        ectodomain_end=ectodomain.end,
        ectodomain_sequence=ectodomain_sequence,
    )


def _resolve_member_by_symbol(
    *,
    analyzer: AntigenAnalyzer,
    gene_symbol: str,
    member_cache: dict[str, ParalogMember | None],
) -> ParalogMember | None:
    cache_key = str(gene_symbol or "").strip().upper()
    if cache_key in member_cache:
        cached = member_cache[cache_key]
        return _clone_member(cached) if cached is not None else None
    try:
        resolution = analyzer.uniprot_client.resolve_target(gene_symbol)
    except Exception:
        member_cache[cache_key] = None
        return None
    target = resolution.target
    features = analyzer.uniprot_client.get_features(resolution.entry)
    topology = derive_ectodomain(features, len(target.sequence), analyzer.config)
    analyzer._apply_secreted_universe_fallback(target=target, topology=topology, features=features)
    ectodomain = topology.ectodomain
    if ectodomain is None:
        member_cache[cache_key] = None
        return None
    ectodomain_sequence = slice_sequence(target.sequence, ectodomain.start, ectodomain.end)
    if not ectodomain_sequence:
        member_cache[cache_key] = None
        return None
    member = ParalogMember(
        gene_symbol=target.gene_symbol or gene_symbol,
        gene_name=target.protein_name,
        accession=target.accession,
        entry_name=target.entry_name,
        sequence_length=len(target.sequence),
        ectodomain_start=ectodomain.start,
        ectodomain_end=ectodomain.end,
        ectodomain_length=len(ectodomain_sequence),
        ectodomain_sequence=ectodomain_sequence,
        matrix_label=target.entry_name,
        notes=[],
    )
    member_cache[cache_key] = member
    return _clone_member(member)


def _clone_member(member: ParalogMember) -> ParalogMember:
    return ParalogMember(
        gene_symbol=member.gene_symbol,
        gene_name=member.gene_name,
        accession=member.accession,
        entry_name=member.entry_name,
        sequence_length=member.sequence_length,
        ectodomain_start=member.ectodomain_start,
        ectodomain_end=member.ectodomain_end,
        ectodomain_length=member.ectodomain_length,
        ectodomain_sequence=member.ectodomain_sequence,
        matrix_label=member.matrix_label,
        sources=list(member.sources),
        notes=list(member.notes),
        identity_to_target=member.identity_to_target,
        coverage_to_target=member.coverage_to_target,
        evidence_score=member.evidence_score,
    )


def _build_target_specific_paralog_payload(
    *,
    analyzer: AntigenAnalyzer,
    seed: ResolvedParalogTarget,
    ensembl_client: EnsemblParalogClient,
    member_cache: dict[str, ParalogMember | None],
    family_alignment_cache: dict[str, list[dict[str, Any]]],
    ensembl_cache: dict[str, list[EnsemblParalogPrediction]],
    hgnc_cache: dict[str, Any],
) -> dict[str, Any]:
    shared_payload = _build_shared_member_payload(
        analyzer=analyzer,
        seeds=[seed],
        ensembl_client=ensembl_client,
        member_cache=member_cache,
        family_alignment_cache=family_alignment_cache,
        ensembl_cache=ensembl_cache,
        hgnc_cache=hgnc_cache,
    )
    return _assemble_target_payload_from_shared(seed=seed, shared_payload=shared_payload)


def _build_shared_family_payload(
    *,
    analyzer: AntigenAnalyzer,
    seeds: list[ResolvedParalogTarget],
    ensembl_client: EnsemblParalogClient,
    member_cache: dict[str, ParalogMember | None],
    family_alignment_cache: dict[str, list[dict[str, Any]]],
    ensembl_cache: dict[str, list[EnsemblParalogPrediction]],
    hgnc_cache: dict[str, Any],
) -> dict[str, Any]:
    shared_payload = _build_shared_member_payload(
        analyzer=analyzer,
        seeds=seeds,
        ensembl_client=ensembl_client,
        member_cache=member_cache,
        family_alignment_cache=family_alignment_cache,
        ensembl_cache=ensembl_cache,
        hgnc_cache=hgnc_cache,
    )
    family_accession = str(seeds[0].canonical_family.accession) if seeds and seeds[0].canonical_family is not None else ""
    shared_payload["metadata"] = {
        **(shared_payload.get("metadata") or {}),
        "shared_family_payload": True,
        "shared_family_accession": family_accession,
        "shared_target_entries": [seed.target.entry_name for seed in seeds if seed.target is not None],
    }
    return shared_payload


def _build_shared_member_payload(
    *,
    analyzer: AntigenAnalyzer,
    seeds: list[ResolvedParalogTarget],
    ensembl_client: EnsemblParalogClient,
    member_cache: dict[str, ParalogMember | None],
    family_alignment_cache: dict[str, list[dict[str, Any]]],
    ensembl_cache: dict[str, list[EnsemblParalogPrediction]],
    hgnc_cache: dict[str, Any],
) -> dict[str, Any]:
    members_by_key: dict[str, ParalogMember] = {}
    family_names: list[str] = []
    source_counts = {"ensembl": 0, "interpro_family": 0, "hgnc": 0}
    metadata: dict[str, Any] = {}
    seed_sequences = [seed.ectodomain_sequence for seed in seeds if seed.ectodomain_sequence]
    if not seeds:
        return {}
    canonical_family = seeds[0].canonical_family
    if canonical_family is not None and canonical_family.name:
        family_names.append(canonical_family.name)

    for seed in seeds:
        target = seed.target
        if target is None:
            continue
        target_member = ParalogMember(
            gene_symbol=target.gene_symbol or target.entry_name,
            gene_name=target.protein_name,
            accession=target.accession,
            entry_name=target.entry_name,
            sequence_length=len(target.sequence),
            ectodomain_start=seed.ectodomain_start,
            ectodomain_end=seed.ectodomain_end,
            ectodomain_length=len(seed.ectodomain_sequence),
            ectodomain_sequence=seed.ectodomain_sequence,
            matrix_label=target.entry_name,
            sources=["self"],
            notes=["Target protein in this shared paralog family set."],
            identity_to_target=100.0,
            coverage_to_target=100.0,
            evidence_score=1000.0,
        )
        members_by_key[_member_key(accession=target.accession, entry_name=target.entry_name, gene_symbol=target.gene_symbol)] = target_member

    for member in _load_family_alignment_members(
        analyzer=analyzer,
        canonical_family=canonical_family,
        family_alignment_cache=family_alignment_cache,
    ):
        candidate = ParalogMember(
            gene_symbol=member.get("gene_symbol") or member.get("entry_name") or "unknown",
            gene_name=member.get("protein_name") or member.get("gene_symbol") or "unknown",
            accession=_text(member.get("accession")),
            entry_name=_text(member.get("entry_name")),
            sequence_length=None,
            ectodomain_start=_coerce_int(member.get("ectodomain_start")),
            ectodomain_end=_coerce_int(member.get("ectodomain_end")),
            ectodomain_length=len(str(member.get("ectodomain_sequence") or "")) or None,
            ectodomain_sequence=str(member.get("ectodomain_sequence") or ""),
            matrix_label=_text(member.get("entry_name")) or _text(member.get("accession")) or str(member.get("gene_symbol") or "candidate"),
            sources=["interpro_family"],
            notes=["Included through the precomputed InterPro family alignment set."],
            evidence_score=40.0,
        )
        _update_identity_to_seed_set(seed_sequences=seed_sequences, member=candidate)
        _merge_member(members_by_key, candidate)
        source_counts["interpro_family"] += 1

    use_live_family_sources = source_counts["interpro_family"] == 0
    if not use_live_family_sources:
        metadata["live_family_sources_skipped"] = (
            "Skipped Ensembl/HGNC fallback because precomputed InterPro family members were available."
        )

    seen_gene_symbols: set[str] = set()
    for seed in seeds if use_live_family_sources else []:
        gene_symbol = seed.target.gene_symbol if seed.target is not None else None
        if not gene_symbol or gene_symbol in seen_gene_symbols:
            continue
        seen_gene_symbols.add(gene_symbol)
        if gene_symbol in ensembl_cache:
            paralog_predictions = ensembl_cache[gene_symbol]
        else:
            try:
                paralog_predictions = ensembl_client.fetch_human_paralogs(gene_symbol)
            except Exception as exc:
                paralog_predictions = []
                metadata.setdefault("ensembl_errors", {})[gene_symbol] = str(exc)
            ensembl_cache[gene_symbol] = paralog_predictions
        for prediction in paralog_predictions:
            if not prediction.gene_symbol:
                continue
            candidate = _resolve_member_by_symbol(analyzer=analyzer, gene_symbol=prediction.gene_symbol, member_cache=member_cache)
            if candidate is None:
                continue
            candidate.sources.append("ensembl")
            if prediction.homology_type:
                candidate.notes.append(f"Ensembl homology type: {prediction.homology_type}.")
            if prediction.source_to_target_identity is not None or prediction.target_to_source_identity is not None:
                candidate.notes.append(
                    "Ensembl protein identity source->target {source}%, target->source {target}%.".format(
                        source=_fmt_pct(prediction.source_to_target_identity),
                        target=_fmt_pct(prediction.target_to_source_identity),
                    )
                )
            candidate.evidence_score += 100.0
            _update_identity_to_seed_set(seed_sequences=seed_sequences, member=candidate)
            _merge_member(members_by_key, candidate)
            source_counts["ensembl"] += 1

        if gene_symbol in hgnc_cache:
            hgnc_context = hgnc_cache[gene_symbol]
        else:
            try:
                hgnc_context = analyzer.hgnc_client.fetch_family_context(gene_symbol)
            except Exception as exc:
                hgnc_context = None
                metadata.setdefault("hgnc_errors", {})[gene_symbol] = str(exc)
            hgnc_cache[gene_symbol] = hgnc_context
        if hgnc_context is not None:
            for family_name in hgnc_context.family_names:
                if family_name and family_name not in family_names:
                    family_names.append(family_name)
            for family_member in hgnc_context.members:
                candidate = _resolve_member_by_symbol(analyzer=analyzer, gene_symbol=family_member.gene_symbol, member_cache=member_cache)
                if candidate is None:
                    continue
                candidate.sources.append("hgnc")
                candidate.notes.append("Included through HGNC family membership.")
                candidate.evidence_score += 20.0
                _update_identity_to_seed_set(seed_sequences=seed_sequences, member=candidate)
                _merge_member(members_by_key, candidate)
                source_counts["hgnc"] += 1

    selected = _select_shared_members(
        members=list(members_by_key.values()),
        target_entry_names={seed.target.entry_name for seed in seeds if seed.target is not None},
        max_members=max(analyzer.config.max_paralog_context_members, len(seeds)),
    )
    _assign_unique_labels(selected)
    identity_matrix, coverage_matrix, pairwise_alignments = _build_matrices(selected)
    return {
        "family_names": family_names,
        "identity_matrix_labels": [member.matrix_label for member in selected],
        "identity_matrix": identity_matrix,
        "coverage_matrix_labels": [member.matrix_label for member in selected],
        "coverage_matrix": coverage_matrix,
        "members": [asdict(member) for member in selected],
        "pairwise_alignments": pairwise_alignments,
        "member_count": len(selected),
        "metadata": {
            **metadata,
            "candidate_source_counts": source_counts,
            "selected_member_cap": max(analyzer.config.max_paralog_context_members, len(seeds)),
        },
    }


def _assemble_target_payload_from_shared(*, seed: ResolvedParalogTarget, shared_payload: dict[str, Any]) -> dict[str, Any]:
    payload = dict(seed.base_payload)
    payload.update(
        {
            "family_names": list(shared_payload.get("family_names") or []),
            "identity_matrix_labels": list(shared_payload.get("identity_matrix_labels") or []),
            "identity_matrix": shared_payload.get("identity_matrix") or [],
            "coverage_matrix_labels": list(shared_payload.get("coverage_matrix_labels") or []),
            "coverage_matrix": shared_payload.get("coverage_matrix") or [],
            "members": shared_payload.get("members") or [],
            "pairwise_alignments": shared_payload.get("pairwise_alignments") or [],
            "member_count": int(shared_payload.get("member_count") or 0),
            "metadata": {
                **(shared_payload.get("metadata") or {}),
                "shared_matrix_reused": bool((shared_payload.get("metadata") or {}).get("shared_family_payload")),
            },
        }
    )
    return payload


def _update_identity_to_seed_set(*, seed_sequences: list[str], member: ParalogMember) -> None:
    if not seed_sequences or not member.ectodomain_sequence:
        return
    best_identity: float | None = None
    best_coverage: float | None = None
    for sequence in seed_sequences:
        alignment = global_align(sequence, member.ectodomain_sequence)
        identity = round(100.0 * alignment.matches / len(sequence), 2) if sequence else None
        coverage = round(alignment.coverage, 2)
        if best_identity is None or (identity or 0.0) > (best_identity or 0.0):
            best_identity = identity
        if best_coverage is None or (coverage or 0.0) > (best_coverage or 0.0):
            best_coverage = coverage
    member.identity_to_target = best_identity
    member.coverage_to_target = best_coverage
    member.notes.append(
        "Best directional identity to any family seed {identity}%, coverage {coverage}%.".format(
            identity=_fmt_pct(best_identity),
            coverage=_fmt_pct(best_coverage),
        )
    )
    member.evidence_score += (best_identity or 0.0) / 10.0 + (best_coverage or 0.0) / 20.0


def _select_shared_members(*, members: list[ParalogMember], target_entry_names: set[str], max_members: int) -> list[ParalogMember]:
    target_members = [member for member in members if member.entry_name in target_entry_names]
    other_members = [member for member in members if member.entry_name not in target_entry_names]
    target_members.sort(key=lambda member: member.entry_name or "")
    other_members.sort(
        key=lambda member: (
            -member.evidence_score,
            -(member.identity_to_target or 0.0),
            -(member.coverage_to_target or 0.0),
            member.gene_symbol,
        )
    )
    limit = max(max_members, len(target_members))
    return [*target_members, *other_members[: max(0, limit - len(target_members))]]


def _build_matrices(members: list[ParalogMember]) -> tuple[list[list[float | None]], list[list[float | None]], list[dict[str, Any]]]:
    identity_matrix: list[list[float | None]] = []
    coverage_matrix: list[list[float | None]] = []
    pairwise_alignments: list[dict[str, Any]] = []
    for row_member in members:
        identity_row: list[float | None] = []
        coverage_row: list[float | None] = []
        for col_member in members:
            if row_member.matrix_label == col_member.matrix_label:
                identity_row.append(100.0)
                coverage_row.append(100.0)
                continue
            alignment = global_align(row_member.ectodomain_sequence, col_member.ectodomain_sequence)
            if not row_member.ectodomain_sequence:
                identity_row.append(None)
                coverage_row.append(None)
                continue
            identity_row.append(round(100.0 * alignment.matches / len(row_member.ectodomain_sequence), 2))
            coverage_row.append(round(alignment.coverage, 2))
        identity_matrix.append(identity_row)
        coverage_matrix.append(coverage_row)
    for index, query_member in enumerate(members):
        for subject_member in members[index + 1 :]:
            alignment = global_align(query_member.ectodomain_sequence, subject_member.ectodomain_sequence)
            pairwise_alignments.append(
                {
                    "query_label": query_member.matrix_label,
                    "subject_label": subject_member.matrix_label,
                    "query_accession": query_member.accession,
                    "subject_accession": subject_member.accession,
                    "matches": alignment.matches,
                    "aligned_positions": alignment.aligned_positions,
                    "query_to_subject_identity": round(100.0 * alignment.matches / len(query_member.ectodomain_sequence), 2)
                    if query_member.ectodomain_sequence
                    else None,
                    "subject_to_query_identity": round(100.0 * alignment.matches / len(subject_member.ectodomain_sequence), 2)
                    if subject_member.ectodomain_sequence
                    else None,
                    "query_to_subject_coverage": round(alignment.coverage, 2),
                    "aligned_query": alignment.aligned_query,
                    "aligned_subject": alignment.aligned_subject,
                }
            )
    return identity_matrix, coverage_matrix, pairwise_alignments


def _load_family_alignment_members(
    *,
    analyzer: AntigenAnalyzer,
    canonical_family,
    family_alignment_cache: dict[str, list[dict[str, Any]]],
) -> list[dict[str, Any]]:
    if canonical_family is None or not str(canonical_family.accession).startswith("IPR"):
        return []
    family_accession = str(canonical_family.accession)
    if family_accession in family_alignment_cache:
        return family_alignment_cache[family_accession]
    index_path = analyzer.config.precomputed_family_alignment_index_path
    if not index_path.exists():
        return []
    try:
        with index_path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            family_json_path = next(
                (
                    str(row.get("json_path") or "").strip()
                    for row in reader
                    if str(row.get("family_accession") or "").strip() == family_accession
                ),
                "",
            )
    except Exception:
        return []
    if not family_json_path:
        return []
    resolved = _resolve_data_path(family_json_path, base=index_path.parent)
    if resolved is None or not resolved.exists():
        return []
    try:
        payload = json.loads(resolved.read_text(encoding="utf-8"))
    except Exception:
        return []
    members = payload.get("members")
    loaded_members = [item for item in members if isinstance(item, dict)] if isinstance(members, list) else []
    family_alignment_cache[family_accession] = loaded_members
    return loaded_members


def _merge_member(members_by_key: dict[str, ParalogMember], candidate: ParalogMember) -> None:
    key = _member_key(accession=candidate.accession, entry_name=candidate.entry_name, gene_symbol=candidate.gene_symbol)
    existing = members_by_key.get(key)
    if existing is None:
        candidate.sources = sorted(set(candidate.sources))
        candidate.notes = _dedupe_notes(candidate.notes)
        members_by_key[key] = candidate
        return
    existing.sources = sorted(set([*existing.sources, *candidate.sources]))
    existing.notes = _dedupe_notes([*existing.notes, *candidate.notes])
    existing.evidence_score = max(existing.evidence_score, candidate.evidence_score)
    if existing.identity_to_target is None or (candidate.identity_to_target or 0.0) > (existing.identity_to_target or 0.0):
        existing.identity_to_target = candidate.identity_to_target
    if existing.coverage_to_target is None or (candidate.coverage_to_target or 0.0) > (existing.coverage_to_target or 0.0):
        existing.coverage_to_target = candidate.coverage_to_target
    if existing.ectodomain_length is None and candidate.ectodomain_length is not None:
        existing.ectodomain_length = candidate.ectodomain_length
    if (not existing.ectodomain_sequence) and candidate.ectodomain_sequence:
        existing.ectodomain_sequence = candidate.ectodomain_sequence
    if existing.ectodomain_start is None and candidate.ectodomain_start is not None:
        existing.ectodomain_start = candidate.ectodomain_start
    if existing.ectodomain_end is None and candidate.ectodomain_end is not None:
        existing.ectodomain_end = candidate.ectodomain_end
    if existing.accession is None and candidate.accession is not None:
        existing.accession = candidate.accession
    if existing.entry_name is None and candidate.entry_name is not None:
        existing.entry_name = candidate.entry_name


def _update_identity_to_target(*, target_sequence: str, member: ParalogMember) -> None:
    if not target_sequence or not member.ectodomain_sequence:
        return
    alignment = global_align(target_sequence, member.ectodomain_sequence)
    member.identity_to_target = round(100.0 * alignment.matches / len(target_sequence), 2) if target_sequence else None
    member.coverage_to_target = round(alignment.coverage, 2)
    member.notes.append(
        "Target ectodomain directional identity {identity}%, coverage {coverage}%.".format(
            identity=_fmt_pct(member.identity_to_target),
            coverage=_fmt_pct(member.coverage_to_target),
        )
    )
    member.evidence_score += (member.identity_to_target or 0.0) / 10.0 + (member.coverage_to_target or 0.0) / 20.0


def _member_key(*, accession: str | None, entry_name: str | None, gene_symbol: str | None) -> str:
    return "|".join(
        [
            (accession or "").strip().upper(),
            (entry_name or "").strip().upper(),
            (gene_symbol or "").strip().upper(),
        ]
    )


def _is_target_member(*, target, member: dict[str, Any]) -> bool:
    accession = _text(member.get("accession"))
    entry_name = _text(member.get("entry_name"))
    gene_symbol = _text(member.get("gene_symbol"))
    return bool(
        (accession and accession == target.accession)
        or (entry_name and entry_name == target.entry_name)
        or (gene_symbol and target.gene_symbol and gene_symbol == target.gene_symbol)
    )


def _resolve_data_path(value: str, *, base: Path) -> Path | None:
    path = Path(value)
    candidates = [path, Path.cwd() / path, base / path, base.parent / path]
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    return None


def _write_paralog_summary(*, summary_tsv: Path, summary_json: Path, summaries: list[dict[str, Any]]) -> None:
    fieldnames = [
        "input_index",
        "query",
        "target_entry_name",
        "target_accession",
        "target_gene_symbol",
        "canonical_family_accession",
        "canonical_family_name",
        "member_count",
        "status",
        "error",
        "json_path",
    ]
    with summary_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(summaries)
    summary_json.write_text(json.dumps(summaries, indent=2) + "\n", encoding="utf-8")


def _safe_filename(value: str) -> str:
    return "".join(ch if ch.isalnum() or ch in {"_", "-"} else "_" for ch in value)


def _text(value: Any) -> str | None:
    text = str(value or "").strip()
    return text or None


def _coerce_int(value: Any) -> int | None:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _coerce_float(value: Any) -> float | None:
    try:
        return round(float(value), 2)
    except (TypeError, ValueError):
        return None


def _fmt_pct(value: float | None) -> str:
    return f"{value:.2f}" if value is not None else "n/a"


def _dedupe_notes(notes: list[str]) -> list[str]:
    seen: set[str] = set()
    ordered: list[str] = []
    for note in notes:
        cleaned = str(note or "").strip()
        if not cleaned or cleaned in seen:
            continue
        seen.add(cleaned)
        ordered.append(cleaned)
    return ordered
