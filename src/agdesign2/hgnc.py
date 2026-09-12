from __future__ import annotations

import csv
from dataclasses import dataclass
import json
from pathlib import Path
from urllib.parse import quote
from typing import Any

from .config import AnalysisConfig
from .exceptions import ExternalServiceError
from .http import HttpClient
from .models import FamilyContext, FamilyMember
from .sequence_utils import global_align, slice_sequence
from .topology import derive_ectodomain
from .uniprot import UniProtClient


INTEGRIN_ALPHA_GENES = (
    "ITGA1",
    "ITGA2",
    "ITGA2B",
    "ITGA3",
    "ITGA4",
    "ITGA5",
    "ITGA6",
    "ITGA7",
    "ITGA8",
    "ITGA9",
    "ITGA10",
    "ITGA11",
    "ITGAD",
    "ITGAE",
    "ITGAL",
    "ITGAM",
    "ITGAV",
    "ITGAX",
)

INTEGRIN_BETA_GENES = (
    "ITGB1",
    "ITGB2",
    "ITGB3",
    "ITGB4",
    "ITGB5",
    "ITGB6",
    "ITGB7",
    "ITGB8",
)


@dataclass(slots=True)
class OrthologPrediction:
    species: str
    gene_symbol: str | None
    gene_name: str | None
    ncbi_gene_id: str | None
    ensembl_id: str | None
    mod_id: str | None
    gene_symbol_source: str | None
    evidence_count: int


@dataclass(slots=True)
class OrthologLookup:
    human_gene_symbol: str
    human_ncbi_gene_id: str | None
    predictions: dict[str, OrthologPrediction]


class HGNCClient:
    def __init__(self, http: HttpClient, uniprot_client: UniProtClient, config: AnalysisConfig) -> None:
        self.http = http
        self.uniprot_client = uniprot_client
        self.config = config
        self._local_report_index: dict[str, str] | None = None
        self._local_report_cache: dict[str, dict[str, Any] | None] = {}
        self._identity_cache: dict[tuple[str, str], tuple[float | None, float | None]] = {}

    def fetch_family_context(self, gene_symbol: str | None) -> FamilyContext | None:
        if not gene_symbol:
            return None
        url = f"https://rest.genenames.org/fetch/symbol/{gene_symbol}"
        try:
            payload = self.http.fetch_json(
                url,
                headers={"Accept": "application/json"},
                cache_namespace="hgnc",
            )
        except ExternalServiceError as exc:
            fallback = self._fallback_family_context(gene_symbol, [], exc)
            if fallback is not None:
                return fallback
            raise
        response = payload.get("response", {})
        docs = response.get("docs", [])
        if not docs:
            return self._fallback_family_context(gene_symbol, [], ExternalServiceError("HGNC symbol lookup returned no documents."))
        doc = docs[0]
        families = doc.get("gene_group") or []
        if isinstance(families, str):
            families = [families]
        family_ids = doc.get("gene_group_id") or []
        if isinstance(family_ids, int):
            family_ids = [family_ids]
        try:
            selected_family = self._select_preferred_family(families, family_ids)
        except ExternalServiceError as exc:
            fallback = self._fallback_family_context(gene_symbol, families, exc)
            if fallback is not None:
                return fallback
            raise
        if selected_family is None:
            fallback = self._fallback_family_context(
                gene_symbol,
                families,
                ExternalServiceError("HGNC family lookup returned no members."),
            )
            if fallback is not None:
                return fallback
            return None
        selected_family_name, selected_family_id, selected_docs = selected_family
        if len(selected_docs) > self.config.max_family_context_members:
            return None
        members = self._build_family_members_from_docs(selected_docs)
        labels = [member.gene_symbol for member in members]
        matrix = self._build_identity_matrix(members)
        return FamilyContext(
            gene_symbol=gene_symbol,
            source="HGNC",
            family_names=[selected_family_name],
            members=members,
            identity_matrix_labels=labels,
            identity_matrix=matrix,
            metadata={
                "hgnc_id": doc.get("hgnc_id"),
                "gene_group_id": family_ids,
                "selected_gene_group_id": selected_family_id,
                "all_gene_groups": families,
            },
        )

    def fetch_hcop_orthologs(self, gene_symbol: str | None) -> OrthologLookup | None:
        if not gene_symbol:
            return None
        params = [
            "q_type=symbol",
            f"q={quote(gene_symbol)}",
            "q_tax_id=9606",
            "t_tax_id=10090",
            "t_tax_id=9541",
        ]
        url = f"https://www.genenames.org/cgi-bin/orthologs/hcop?{'&'.join(params)}"
        payload = self.http.fetch_json(
            url,
            headers={"Accept": "application/json"},
            cache_namespace="hcop",
        )
        predictions = payload.get("predictions", {})
        selected: dict[str, OrthologPrediction] = {}
        species_aliases = {
            "mouse": ("mouse",),
            "macaca_fascicularis": ("macaque",),
        }
        for species, aliases in species_aliases.items():
            candidates: list[dict] = []
            for alias in aliases:
                values = predictions.get(alias)
                if isinstance(values, list):
                    candidates.extend(item for item in values if isinstance(item, dict))
            best = self._select_best_hcop_prediction(species, candidates)
            if best is not None:
                selected[species] = best
        return OrthologLookup(
            human_gene_symbol=str(payload.get("gene_symbol") or gene_symbol),
            human_ncbi_gene_id=str(payload.get("ncbi_gene_id") or "") or None,
            predictions=selected,
        )

    def _select_preferred_family(
        self,
        family_names: list[str],
        family_ids: list[int],
    ) -> tuple[str, int, list[dict]] | None:
        if not family_ids:
            return None
        candidates: list[tuple[str, int, list[dict]]] = []
        last_error: ExternalServiceError | None = None
        for index, family_id in enumerate(family_ids):
            family_name = family_names[index] if index < len(family_names) else f"HGNC family {family_id}"
            try:
                docs = self._fetch_family_docs_for_id(family_id)
            except ExternalServiceError as exc:
                last_error = exc
                continue
            if docs:
                candidates.append((family_name, family_id, docs))
        if candidates:
            return min(candidates, key=lambda item: (len(item[2]), item[0]))
        if last_error is not None:
            raise last_error
        return None

    def _fallback_family_context(
        self,
        gene_symbol: str,
        family_names: list[str],
        error: ExternalServiceError,
    ) -> FamilyContext | None:
        if gene_symbol in INTEGRIN_ALPHA_GENES:
            family_label = family_names or ["Integrin alpha chains"]
            members = self._build_symbol_fallback_members(INTEGRIN_ALPHA_GENES, resolve_from_uniprot=False)
        elif gene_symbol in INTEGRIN_BETA_GENES:
            family_label = family_names or ["Integrin beta chains"]
            members = self._build_symbol_fallback_members(INTEGRIN_BETA_GENES, resolve_from_uniprot=False)
        else:
            return None
        if len(members) > self.config.max_family_context_members:
            return None
        labels = [member.gene_symbol for member in members]
        matrix = self._build_identity_matrix(members)
        return FamilyContext(
            gene_symbol=gene_symbol,
            source="HGNC fallback",
            family_names=family_label,
            members=members,
            identity_matrix_labels=labels,
            identity_matrix=matrix,
            metadata={
                "fallback": True,
                "fallback_reason": str(error),
                "identity_matrix_note": "Matrix may contain n/a values because fallback members were kept without live UniProt resolution.",
            },
        )

    def _build_symbol_fallback_members(
        self,
        symbols: tuple[str, ...],
        *,
        resolve_from_uniprot: bool,
    ) -> list[FamilyMember]:
        members: list[FamilyMember] = []
        for symbol in symbols:
            member = self._resolve_member_by_symbol(symbol) if resolve_from_uniprot else None
            if member is None:
                members.append(
                    FamilyMember(
                        gene_symbol=symbol,
                        gene_name=symbol,
                        accession=None,
                        entry_name=None,
                        notes=["Fallback family membership inferred without UniProt resolution."],
                    )
                )
            else:
                members.append(member)
        return members

    def _resolve_member_by_symbol(self, symbol: str) -> FamilyMember | None:
        try:
            payload = self.uniprot_client.search(
                f"(gene_exact:{symbol}) AND (organism_id:9606) AND (reviewed:true)",
                size=3,
            )
        except Exception:
            return None
        results = payload.get("results", [])
        if not results:
            return None
        accession = results[0].get("primaryAccession")
        if not accession:
            return None
        try:
            entry = self.uniprot_client._fetch_entry(accession)
        except Exception:
            return FamilyMember(
                gene_symbol=symbol,
                gene_name=symbol,
                accession=accession,
                entry_name=None,
                notes=["UniProt accession resolved, but entry fetch failed during fallback."],
            )
        sequence = entry.get("sequence", {}).get("value", "")
        features = self.uniprot_client.get_features(entry)
        topology = derive_ectodomain(features, len(sequence), self.config)
        ectodomain = topology.ectodomain
        return FamilyMember(
            gene_symbol=symbol,
            gene_name=entry.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", symbol),
            accession=accession,
            entry_name=entry.get("uniProtkbId"),
            sequence_length=len(sequence) if sequence else None,
            ectodomain_start=ectodomain.start if ectodomain else None,
            ectodomain_end=ectodomain.end if ectodomain else None,
            ectodomain_length=(ectodomain.end - ectodomain.start + 1) if ectodomain else None,
            notes=[note.message for note in topology.notes] or ["Fallback family membership inferred from gene symbol pattern."],
        )

    def _fetch_family_docs_for_id(self, family_id: int) -> list[dict]:
        payload = self.http.fetch_json(
            f"https://rest.genenames.org/fetch/gene_group_id/{family_id}",
            headers={"Accept": "application/json"},
            cache_namespace="hgnc_family",
        )
        return payload.get("response", {}).get("docs", [])

    def _select_best_hcop_prediction(
        self,
        species: str,
        candidates: list[dict],
    ) -> OrthologPrediction | None:
        if not candidates:
            return None
        preferred_sources = {
            "mouse": {"MGI": 0, "VGNC": 1, "NCBI": 2},
            "macaca_fascicularis": {"VGNC": 0, "NCBI": 1},
        }.get(species, {})

        def sort_key(candidate: dict) -> tuple[int, int, str]:
            source = str(candidate.get("gene_symbol_source") or "")
            source_rank = preferred_sources.get(source, 5)
            evidence_count = len(candidate.get("evidence") or [])
            gene_symbol = str(candidate.get("gene_symbol") or "")
            return (source_rank, -evidence_count, gene_symbol)

        chosen = sorted(candidates, key=sort_key)[0]
        return OrthologPrediction(
            species=species,
            gene_symbol=str(chosen.get("gene_symbol") or "") or None,
            gene_name=str(chosen.get("gene_name") or "") or None,
            ncbi_gene_id=str(chosen.get("ncbi_gene_id") or "") or None,
            ensembl_id=str(chosen.get("ensembl_id") or "") or None,
            mod_id=str(chosen.get("mod_id") or "") or None,
            gene_symbol_source=str(chosen.get("gene_symbol_source") or "") or None,
            evidence_count=len(chosen.get("evidence") or []),
        )

    def _build_family_members_from_docs(self, docs: list[dict]) -> list[FamilyMember]:
        members: dict[str, FamilyMember] = {}
        for doc in docs:
            symbol = doc.get("symbol")
            if not symbol or symbol in members:
                continue
            accession, entry_name, ectodomain_start, ectodomain_end, ectodomain_length, notes = (
                self._resolve_member_uniprot(doc)
            )
            sequence_length = None
            if accession:
                try:
                    entry = self.uniprot_client.fetch_entry_by_name(entry_name) if entry_name else None
                    if not entry and accession:
                        entry = self.uniprot_client._fetch_entry(accession)
                    sequence_length = len(entry.get("sequence", {}).get("value", "")) if entry else None
                except Exception:
                    sequence_length = None
            members[symbol] = FamilyMember(
                gene_symbol=symbol,
                gene_name=doc.get("name", symbol),
                accession=accession,
                entry_name=entry_name,
                sequence_length=sequence_length,
                ectodomain_start=ectodomain_start,
                ectodomain_end=ectodomain_end,
                ectodomain_length=ectodomain_length,
                notes=notes,
            )
        return sorted(members.values(), key=lambda item: item.gene_symbol)

    def _resolve_member_uniprot(
        self,
        doc: dict,
    ) -> tuple[str | None, str | None, int | None, int | None, int | None, list[str]]:
        local = self._resolve_member_from_local_report(doc)
        if local is not None:
            return local
        notes: list[str] = []
        accessions = doc.get("uniprot_ids") or []
        if isinstance(accessions, str):
            accessions = [accessions]
        for accession in accessions:
            try:
                entry = self.uniprot_client._fetch_entry(accession)
            except Exception as exc:
                notes.append(f"UniProt fetch failed for {accession}: {exc}")
                continue
            sequence = entry.get("sequence", {}).get("value", "")
            features = self.uniprot_client.get_features(entry)
            topology = derive_ectodomain(features, len(sequence), self.config)
            ectodomain = topology.ectodomain
            member_notes = notes + [note.message for note in topology.notes]
            if ectodomain:
                return (
                    accession,
                    entry.get("uniProtkbId"),
                    ectodomain.start,
                    ectodomain.end,
                    ectodomain.end - ectodomain.start + 1,
                    member_notes,
                )
            return (
                accession,
                entry.get("uniProtkbId"),
                None,
                None,
                None,
                member_notes,
            )
        return None, None, None, None, None, notes or ["No UniProt accession available from HGNC."]

    def _build_identity_matrix(self, members: list[FamilyMember]) -> list[list[float | None]]:
        sequences: dict[str, str | None] = {}
        for member in members:
            if not member.accession or member.ectodomain_start is None or member.ectodomain_end is None:
                sequences[member.gene_symbol] = None
                continue
            local_sequence = self._resolve_member_ectodomain_sequence_from_local_report(member)
            if local_sequence is not None:
                sequences[member.gene_symbol] = local_sequence
                continue
            try:
                entry = self.uniprot_client._fetch_entry(member.accession)
                sequence = entry.get("sequence", {}).get("value", "")
                sequences[member.gene_symbol] = slice_sequence(
                    sequence, member.ectodomain_start, member.ectodomain_end
                )
            except Exception:
                sequences[member.gene_symbol] = None

        matrix: list[list[float | None]] = []
        for row_member in members:
            row: list[float | None] = []
            row_sequence = sequences.get(row_member.gene_symbol)
            for col_member in members:
                col_sequence = sequences.get(col_member.gene_symbol)
                if row_sequence is None or col_sequence is None:
                    row.append(None)
                    continue
                if row_member.gene_symbol == col_member.gene_symbol:
                    row.append(100.0)
                    continue
                row_identity, _ = self._pairwise_directional_identity(row_sequence, col_sequence)
                row.append(row_identity)
            matrix.append(row)
        return matrix

    def _pairwise_directional_identity(
        self,
        row_sequence: str,
        col_sequence: str,
    ) -> tuple[float | None, float | None]:
        cache_key = (row_sequence, col_sequence)
        cached = self._identity_cache.get(cache_key)
        if cached is not None:
            return cached

        reciprocal_key = (col_sequence, row_sequence)
        reciprocal = self._identity_cache.get(reciprocal_key)
        if reciprocal is not None:
            result = (reciprocal[1], reciprocal[0])
            self._identity_cache[cache_key] = result
            return result

        if not row_sequence or not col_sequence:
            result = (None, None)
            self._identity_cache[cache_key] = result
            return result

        matches = global_align(row_sequence, col_sequence).matches

        row_identity = round(100.0 * matches / len(row_sequence), 2) if row_sequence else None
        col_identity = round(100.0 * matches / len(col_sequence), 2) if col_sequence else None
        result = (row_identity, col_identity)
        self._identity_cache[cache_key] = result
        self._identity_cache[reciprocal_key] = (col_identity, row_identity)
        return result

    def _resolve_member_from_local_report(
        self,
        doc: dict[str, Any],
    ) -> tuple[str | None, str | None, int | None, int | None, int | None, list[str]] | None:
        report = self._load_local_report_for_doc(doc)
        if report is None:
            return None
        target = report.get("target") or {}
        ectodomain = report.get("ectodomain") or {}
        start = ectodomain.get("start")
        end = ectodomain.get("end")
        ectodomain_length = None
        if start is not None and end is not None:
            ectodomain_length = int(end) - int(start) + 1
        return (
            str(target.get("accession") or "") or None,
            str(target.get("entry_name") or "") or None,
            int(start) if start is not None else None,
            int(end) if end is not None else None,
            ectodomain_length,
            ["Resolved from local batch report."],
        )

    def _resolve_member_ectodomain_sequence_from_local_report(self, member: FamilyMember) -> str | None:
        report = self._load_local_report_for_member(member)
        if report is None:
            return None
        target = report.get("target") or {}
        ectodomain = report.get("ectodomain") or {}
        start = member.ectodomain_start if member.ectodomain_start is not None else ectodomain.get("start")
        end = member.ectodomain_end if member.ectodomain_end is not None else ectodomain.get("end")
        sequence = str(target.get("sequence") or "")
        if not sequence or start is None or end is None:
            return None
        return slice_sequence(sequence, int(start), int(end))

    def _load_local_report_for_member(self, member: FamilyMember) -> dict[str, Any] | None:
        index = self._load_local_report_index()
        candidates = [
            member.entry_name,
            index.get(str(member.accession or "").strip()),
            index.get(str(member.gene_symbol or "").strip().upper()),
        ]
        for entry_name in candidates:
            payload = self._load_local_report_payload(entry_name)
            if payload is not None:
                return payload
        return None

    def _load_local_report_for_doc(self, doc: dict[str, Any]) -> dict[str, Any] | None:
        index = self._load_local_report_index()
        accessions = doc.get("uniprot_ids") or []
        if isinstance(accessions, str):
            accessions = [accessions]
        symbol = str(doc.get("symbol") or "").strip().upper()
        candidates = [index.get(symbol)]
        candidates.extend(index.get(str(accession).strip()) for accession in accessions)
        for entry_name in candidates:
            payload = self._load_local_report_payload(entry_name)
            if payload is not None:
                return payload
        return None

    def _load_local_report_index(self) -> dict[str, str]:
        if self._local_report_index is not None:
            return self._local_report_index
        index: dict[str, str] = {}
        path = self.config.precomputed_ortholog_table_path
        if not path.exists():
            self._local_report_index = index
            return index
        try:
            with path.open(encoding="utf-8") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                for row in reader:
                    entry_name = str(row.get("input_uniprot_name") or "").strip()
                    if not entry_name:
                        continue
                    accession = str(row.get("human_uniprot_accession") or row.get("input_uniprot_id") or "").strip()
                    gene_symbol = str(row.get("human_gene_symbol") or row.get("input_gene") or "").strip().upper()
                    if accession and accession not in index:
                        index[accession] = entry_name
                    if gene_symbol and gene_symbol not in index:
                        index[gene_symbol] = entry_name
        except Exception:
            index = {}
        self._local_report_index = index
        return index

    def _load_local_report_payload(self, entry_name: str | None) -> dict[str, Any] | None:
        if not entry_name:
            return None
        key = str(entry_name).strip()
        if not key:
            return None
        if key in self._local_report_cache:
            return self._local_report_cache[key]
        report_dir = self.config.precomputed_ortholog_table_path.parent
        report_path = report_dir / f"{key.lower()}_report.json"
        if not report_path.exists():
            self._local_report_cache[key] = None
            return None
        try:
            payload = json.loads(report_path.read_text(encoding="utf-8"))
        except Exception:
            payload = None
        self._local_report_cache[key] = payload
        return payload
