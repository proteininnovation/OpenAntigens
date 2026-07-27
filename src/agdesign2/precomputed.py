from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .config import AnalysisConfig
from .models import CanonicalFamily, FamilyContext, FamilyMember, HomologyRecord, Region, TargetRecord
from .sequence_utils import global_align, map_query_region_to_subject, slice_sequence


_MATCH_OK_STATUSES = {"ok", "completed", ""}


class _MatchIndex:
    """O(1) replacement for the linear scans that resolve a target to a
    precomputed ortholog/paralog row.

    Both scans matched the first in-file-order row where either the entry name
    equals the target's, or (accession, gene_symbol) both equal the target's.
    The index preserves that exact tie-break by storing each key's first
    (position, item) and, on lookup, returning the lower-positioned of the two
    candidate matches — identical results, without the O(rows) scan per target.
    """

    __slots__ = ("_by_entry", "_by_pair")

    def __init__(self, items, *, status_of, entry_of, accession_of, gene_of) -> None:
        self._by_entry: dict[str, tuple[int, Any]] = {}
        self._by_pair: dict[tuple[str, str], tuple[int, Any]] = {}
        for pos, item in enumerate(items):
            if status_of(item) not in _MATCH_OK_STATUSES:
                continue
            entry = entry_of(item)
            if entry and entry not in self._by_entry:
                self._by_entry[entry] = (pos, item)
            accession = accession_of(item)
            gene = gene_of(item)
            if accession and gene:
                key = (accession, gene)
                if key not in self._by_pair:
                    self._by_pair[key] = (pos, item)

    def lookup(self, *, entry_name, accession, gene_symbol):
        candidates: list[tuple[int, Any]] = []
        if entry_name:
            hit = self._by_entry.get(entry_name)
            if hit is not None:
                candidates.append(hit)
        if accession and gene_symbol:
            hit = self._by_pair.get((accession, gene_symbol))
            if hit is not None:
                candidates.append(hit)
        if not candidates:
            return None
        return min(candidates, key=lambda item: item[0])[1]


@dataclass(slots=True)
class PrecomputedOrthologRecord:
    row: dict[str, str]

    @property
    def status(self) -> str:
        return str(self.row.get("status") or "").strip().lower()

    @property
    def entry_name(self) -> str | None:
        return _clean(self.row.get("input_uniprot_name"))

    @property
    def accession(self) -> str | None:
        return _clean(self.row.get("human_uniprot_accession")) or _clean(self.row.get("input_uniprot_id"))

    @property
    def gene_symbol(self) -> str | None:
        return _clean(self.row.get("human_gene_symbol")) or _clean(self.row.get("input_gene"))

    @property
    def canonical_family_accession(self) -> str | None:
        return _clean(self.row.get("canonical_family_accession"))

    @property
    def canonical_family_name(self) -> str | None:
        return _clean(self.row.get("canonical_family_name"))


class PrecomputedReferenceStore:
    def __init__(self, config: AnalysisConfig) -> None:
        self.config = config
        self._ortholog_rows: list[PrecomputedOrthologRecord] | None = None
        self._family_index: dict[str, str] | None = None
        self._paralog_index: list[dict[str, str]] | None = None
        self._ortholog_match_index: _MatchIndex | None = None
        self._paralog_match_index: _MatchIndex | None = None

    def find_ortholog_record(self, *, target: TargetRecord) -> PrecomputedOrthologRecord | None:
        if not self.config.prefer_precomputed_references:
            return None
        if self._ortholog_match_index is None:
            self._ortholog_match_index = _MatchIndex(
                self._load_ortholog_rows(),
                status_of=lambda record: record.status,
                entry_of=lambda record: record.entry_name,
                accession_of=lambda record: record.accession,
                gene_of=lambda record: record.gene_symbol,
            )
        return self._ortholog_match_index.lookup(
            entry_name=target.entry_name,
            accession=target.accession,
            gene_symbol=target.gene_symbol,
        )

    def canonical_family_from_record(self, record: PrecomputedOrthologRecord | None) -> CanonicalFamily | None:
        if record is None:
            return None
        accession = record.canonical_family_accession
        if not accession or not accession.startswith("IPR"):
            return None
        return CanonicalFamily(
            accession=accession,
            name=record.canonical_family_name or accession,
            source_database="INTERPRO",
            fragment_count=1,
            notes=["Loaded from precomputed ortholog reference table."],
        )

    def load_family_context(
        self,
        *,
        target: TargetRecord,
        canonical_family: CanonicalFamily | None,
        max_members: int,
    ) -> FamilyContext | None:
        if not self.config.prefer_precomputed_references or canonical_family is None:
            return None
        family_path = self._load_family_index().get(canonical_family.accession)
        if not family_path:
            return None
        resolved = self._resolve_data_path(family_path, base=self.config.precomputed_family_alignment_index_path.parent)
        if resolved is None or not resolved.exists():
            return None
        try:
            payload = json.loads(resolved.read_text(encoding="utf-8"))
        except Exception:
            return None
        members_payload = payload.get("members") or []
        if not isinstance(members_payload, list):
            return None
        if len(members_payload) > max_members:
            return None
        members: list[FamilyMember] = []
        for item in members_payload:
            if not isinstance(item, dict):
                continue
            members.append(
                FamilyMember(
                    gene_symbol=str(item.get("gene_symbol") or ""),
                    gene_name=str(item.get("protein_name") or item.get("gene_symbol") or ""),
                    accession=_clean(item.get("accession")),
                    entry_name=_clean(item.get("entry_name")),
                    ectodomain_start=_coerce_int(item.get("ectodomain_start")),
                    ectodomain_end=_coerce_int(item.get("ectodomain_end")),
                    ectodomain_length=len(str(item.get("ectodomain_sequence") or "")) or None,
                    sequence_length=None,
                    notes=["Loaded from precomputed family alignment set."],
                )
            )
        labels = payload.get("identity_matrix_labels") or []
        matrix = payload.get("identity_matrix") or []
        if not isinstance(labels, list) or not isinstance(matrix, list):
            return None
        return FamilyContext(
            gene_symbol=target.gene_symbol or target.entry_name,
            source="Precomputed InterPro family alignments",
            family_names=[str(payload.get("family_name") or canonical_family.name)],
            members=members,
            identity_matrix_labels=[str(label) for label in labels],
            identity_matrix=[[(_coerce_float(value)) for value in row] for row in matrix if isinstance(row, list)],
            coverage_matrix_labels=[str(label) for label in (payload.get("coverage_matrix_labels") or labels)],
            coverage_matrix=[
                [(_coerce_float(value)) for value in row]
                for row in (payload.get("coverage_matrix") or [])
                if isinstance(row, list)
            ],
            metadata={
                "family_accession": canonical_family.accession,
                "precomputed_json_path": str(resolved),
                "pairwise_alignments": payload.get("pairwise_alignments") or [],
            },
        )

    def load_paralog_context(
        self,
        *,
        target: TargetRecord,
        max_members: int,
    ) -> FamilyContext | None:
        if not self.config.prefer_precomputed_references:
            return None
        index_row = self._find_paralog_row(target=target)
        if index_row is None:
            return None
        resolved = self._resolve_data_path(
            str(index_row.get("json_path") or ""),
            base=self.config.precomputed_paralog_index_path.parent,
        )
        if resolved is None or not resolved.exists():
            return None
        try:
            payload = json.loads(resolved.read_text(encoding="utf-8"))
        except Exception:
            return None
        members_payload = payload.get("members") or []
        if not isinstance(members_payload, list) or len(members_payload) > max_members:
            return None
        members: list[FamilyMember] = []
        for item in members_payload:
            if not isinstance(item, dict):
                continue
            notes = [str(note) for note in (item.get("notes") or []) if str(note).strip()]
            sources = [str(source) for source in (item.get("sources") or []) if str(source).strip()]
            if sources:
                notes = [f"Sources: {', '.join(sources)}.", *notes]
            members.append(
                FamilyMember(
                    gene_symbol=str(item.get("gene_symbol") or ""),
                    gene_name=str(item.get("gene_name") or item.get("gene_symbol") or ""),
                    accession=_clean(item.get("accession")),
                    entry_name=_clean(item.get("entry_name")),
                    ectodomain_start=_coerce_int(item.get("ectodomain_start")),
                    ectodomain_end=_coerce_int(item.get("ectodomain_end")),
                    ectodomain_length=_coerce_int(item.get("ectodomain_length")),
                    sequence_length=_coerce_int(item.get("sequence_length")),
                    notes=notes,
                )
            )
        labels = payload.get("identity_matrix_labels") or []
        matrix = payload.get("identity_matrix") or []
        coverage_labels = payload.get("coverage_matrix_labels") or labels
        coverage_matrix = payload.get("coverage_matrix") or []
        if not isinstance(labels, list) or not isinstance(matrix, list):
            return None
        family_names = [str(name) for name in (payload.get("family_names") or []) if str(name).strip()]
        canonical_family_name = _clean(payload.get("canonical_family_name"))
        if canonical_family_name and canonical_family_name not in family_names:
            family_names.insert(0, canonical_family_name)
        return FamilyContext(
            gene_symbol=target.gene_symbol or target.entry_name,
            source="Precomputed paralog reference set",
            family_names=family_names,
            members=members,
            identity_matrix_labels=[str(label) for label in labels],
            identity_matrix=[[(_coerce_float(value)) for value in row] for row in matrix if isinstance(row, list)],
            coverage_matrix_labels=[str(label) for label in coverage_labels] if isinstance(coverage_labels, list) else [],
            coverage_matrix=[
                [(_coerce_float(value)) for value in row]
                for row in coverage_matrix
                if isinstance(row, list)
            ],
            metadata={
                "precomputed_json_path": str(resolved),
                "identity_matrix_note": str(payload.get("identity_matrix_note") or ""),
                "coverage_matrix_note": str(payload.get("coverage_matrix_note") or ""),
                "pairwise_alignments": payload.get("pairwise_alignments") or [],
                "canonical_family_accession": _clean(payload.get("canonical_family_accession")),
                "candidate_source_counts": payload.get("metadata", {}).get("candidate_source_counts")
                if isinstance(payload.get("metadata"), dict)
                else None,
            },
        )

    def build_species_matches(
        self,
        *,
        target: TargetRecord,
        ectodomain: Region | None,
        record: PrecomputedOrthologRecord,
    ) -> tuple[list[HomologyRecord], list[HomologyRecord], dict[str, dict[str, Any]]]:
        matches: list[HomologyRecord] = []
        homology: list[HomologyRecord] = []
        context: dict[str, dict[str, Any]] = {}
        for species, prefix, organism, taxon_id in (
            ("mouse", "mouse", "Mus musculus", 10090),
            ("macaca_fascicularis", "macaca_fascicularis", "Macaca fascicularis", 9541),
        ):
            accession = _clean(record.row.get(f"{prefix}_refseq_accession"))
            sequence = _clean(record.row.get(f"{prefix}_refseq_sequence"))
            gene_symbol = _clean(record.row.get(f"{prefix}_gene_symbol"))
            species_notes = [
                "Loaded from precomputed ortholog reference table.",
            ]
            if not accession or not sequence:
                match = HomologyRecord(
                    species=species,
                    accession=accession,
                    entry_name=accession,
                    ectodomain_start=None,
                    ectodomain_end=None,
                    available=False,
                    notes=species_notes + ["No precomputed sequence available for this species."],
                )
                matches.append(match)
                homology.append(match)
                context[species] = {"available": False, "notes": match.notes}
                continue
            homolog_target = TargetRecord(
                accession=accession,
                entry_name=accession,
                gene_symbol=gene_symbol,
                protein_name=gene_symbol or accession,
                organism=organism,
                taxon_id=taxon_id,
                sequence=sequence,
            )
            if ectodomain is None:
                homolog_ectodomain = Region(
                    start=1,
                    end=len(homolog_target.sequence),
                    label="Transferred full-length protein",
                    source="precomputed-ortholog-alignment",
                    confidence=0.8,
                )
                homolog_ecto_sequence = homolog_target.sequence
                query_to_subject = global_align(target.sequence, homolog_target.sequence)
            else:
                homolog_ectodomain, homolog_ecto_sequence, query_to_subject = self._transfer_ectodomain(
                    target=target,
                    homolog_target=homolog_target,
                    ectodomain=ectodomain,
                )
            match = HomologyRecord(
                species=species,
                accession=homolog_target.accession,
                entry_name=homolog_target.entry_name,
                ectodomain_start=homolog_ectodomain.start if homolog_ectodomain else None,
                ectodomain_end=homolog_ectodomain.end if homolog_ectodomain else None,
                available=True,
                notes=species_notes,
            )
            matches.append(match)
            if homolog_ectodomain and homolog_ecto_sequence and query_to_subject:
                homology.append(
                    HomologyRecord(
                        species=species,
                        accession=homolog_target.accession,
                        entry_name=homolog_target.entry_name,
                        ectodomain_start=homolog_ectodomain.start,
                        ectodomain_end=homolog_ectodomain.end,
                        identity=round(query_to_subject.identity, 2),
                        coverage=round(query_to_subject.coverage, 2),
                        available=True,
                        notes=species_notes,
                    )
                )
                context[species] = {
                    "available": True,
                    "target": homolog_target,
                    "ectodomain": homolog_ectodomain,
                    "ectodomain_sequence": homolog_ecto_sequence,
                    "alignment": query_to_subject,
                    "notes": species_notes,
                }
            else:
                homology.append(match)
                context[species] = {
                    "available": True,
                    "target": homolog_target,
                    "ectodomain": homolog_ectodomain,
                    "ectodomain_sequence": homolog_ecto_sequence,
                    "alignment": query_to_subject,
                    "notes": species_notes + ["Unable to transfer ectodomain boundaries from the precomputed sequences."],
                }
        return matches, homology, context

    def _transfer_ectodomain(
        self,
        *,
        target: TargetRecord,
        homolog_target: TargetRecord,
        ectodomain: Region | None,
    ) -> tuple[Region | None, str | None, Any | None]:
        if ectodomain is None:
            return None, None, None
        full_alignment = global_align(target.sequence, homolog_target.sequence)
        mapped_start, mapped_end, mapped_sequence, _ = map_query_region_to_subject(
            full_alignment,
            query_start=ectodomain.start,
            query_end=ectodomain.end,
            subject_sequence=homolog_target.sequence,
        )
        if mapped_start is None or mapped_end is None or mapped_sequence is None:
            return None, None, None
        homolog_ectodomain = Region(
            start=mapped_start,
            end=mapped_end,
            label="Transferred ectodomain",
            source="precomputed-ortholog-alignment",
            confidence=0.8,
        )
        query_ecto = slice_sequence(target.sequence, ectodomain.start, ectodomain.end)
        query_to_subject = global_align(query_ecto, mapped_sequence)
        return homolog_ectodomain, mapped_sequence, query_to_subject

    def _load_ortholog_rows(self) -> list[PrecomputedOrthologRecord]:
        if self._ortholog_rows is not None:
            return self._ortholog_rows
        path = self.config.precomputed_ortholog_table_path
        if not path.exists():
            self._ortholog_rows = []
            return self._ortholog_rows
        try:
            with path.open("r", encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                self._ortholog_rows = [PrecomputedOrthologRecord(dict(row)) for row in reader]
        except Exception:
            self._ortholog_rows = []
        return self._ortholog_rows

    def _load_family_index(self) -> dict[str, str]:
        if self._family_index is not None:
            return self._family_index
        path = self.config.precomputed_family_alignment_index_path
        if not path.exists():
            self._family_index = {}
            return self._family_index
        try:
            with path.open("r", encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                self._family_index = {
                    str(row.get("family_accession") or "").strip(): str(row.get("json_path") or "").strip()
                    for row in reader
                    if str(row.get("family_accession") or "").strip() and str(row.get("json_path") or "").strip()
                }
        except Exception:
            self._family_index = {}
        return self._family_index

    def _find_paralog_row(self, *, target: TargetRecord) -> dict[str, str] | None:
        if self._paralog_match_index is None:
            self._paralog_match_index = _MatchIndex(
                self._load_paralog_index(),
                status_of=lambda row: str(row.get("status") or "").strip().lower(),
                entry_of=lambda row: _clean(row.get("target_entry_name")),
                accession_of=lambda row: _clean(row.get("target_accession")),
                gene_of=lambda row: _clean(row.get("target_gene_symbol")),
            )
        return self._paralog_match_index.lookup(
            entry_name=target.entry_name,
            accession=target.accession,
            gene_symbol=target.gene_symbol,
        )

    def _load_paralog_index(self) -> list[dict[str, str]]:
        if self._paralog_index is not None:
            return self._paralog_index
        path = self.config.precomputed_paralog_index_path
        if not path.exists():
            self._paralog_index = []
            return self._paralog_index
        try:
            with path.open("r", encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                self._paralog_index = [dict(row) for row in reader]
        except Exception:
            self._paralog_index = []
        return self._paralog_index

    def _resolve_data_path(self, value: str, *, base: Path) -> Path | None:
        path = Path(value)
        candidates = [
            path,
            Path.cwd() / path,
            base / path,
            base.parent / path,
        ]
        for candidate in candidates:
            if candidate.exists():
                return candidate.resolve()
        return None


def _clean(value: Any) -> str | None:
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
