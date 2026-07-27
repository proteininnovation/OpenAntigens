from __future__ import annotations

import html
from typing import Any

from .http import HttpClient
from .models import GPCRdbAnnotation, GPCRdbMotif, GPCRdbResidue, GPCRdbSegment


GPCRDB_BASE = "https://gpcrdb.org"
GPCRDB_SERVICES = f"{GPCRDB_BASE}/services"


class GPCRdbClient:
    def __init__(self, http: HttpClient) -> None:
        self.http = http

    def fetch_annotation(
        self,
        *,
        entry_name: str,
        accession: str,
        sequence: str,
    ) -> GPCRdbAnnotation | None:
        protein_payload: dict[str, Any] | None = None
        gpcr_entry_name = str(entry_name or "").lower()
        if gpcr_entry_name:
            protein_payload = self._fetch_protein(gpcr_entry_name)
        if protein_payload is None and accession:
            protein_payload = self._fetch_protein_by_accession(accession)
        if not protein_payload:
            return None

        resolved_entry_name = str(protein_payload.get("entry_name") or gpcr_entry_name).lower()
        residues = self._fetch_residues(resolved_entry_name)
        if residues and sequence:
            residues = [
                residue
                for residue in residues
                if residue.sequence_number <= len(sequence)
                and sequence[residue.sequence_number - 1] == residue.amino_acid
            ]
        family_slug = _clean_text(protein_payload.get("family"))
        family_path = self._fetch_family_path(family_slug) if family_slug else []
        family_name = family_path[-1] if family_path else None
        gpcr_class = family_path[0] if family_path else None
        return GPCRdbAnnotation(
            entry_name=resolved_entry_name,
            accession=_clean_text(protein_payload.get("accession")),
            name=_clean_text(protein_payload.get("name")),
            family_slug=family_slug,
            family_name=family_name,
            family_path=family_path,
            gpcr_class=gpcr_class,
            species=_clean_text(protein_payload.get("species")),
            source=_clean_text(protein_payload.get("source")),
            residue_numbering_scheme=_clean_text(protein_payload.get("residue_numbering_scheme")),
            url=f"{GPCRDB_BASE}/protein/{resolved_entry_name}",
            segments=_segments_from_residues(residues),
            residues=residues,
            conserved_motifs=_conserved_motifs(residues),
            warnings=[] if residues else ["GPCRdb protein resolved, but residue/segment mapping was unavailable."],
        )

    def _fetch_protein(self, entry_name: str) -> dict[str, Any] | None:
        try:
            payload = self.http.fetch_json(
                f"{GPCRDB_SERVICES}/protein/{entry_name}/",
                cache_namespace="gpcrdb",
                cache_key=f"protein:{entry_name}",
            )
        except Exception:
            return None
        return payload if isinstance(payload, dict) and payload.get("entry_name") else None

    def _fetch_protein_by_accession(self, accession: str) -> dict[str, Any] | None:
        try:
            payload = self.http.fetch_json(
                f"{GPCRDB_SERVICES}/protein/accession/{accession}/",
                cache_namespace="gpcrdb",
                cache_key=f"protein_accession:{accession}",
            )
        except Exception:
            return None
        return payload if isinstance(payload, dict) and payload.get("entry_name") else None

    def _fetch_residues(self, entry_name: str) -> list[GPCRdbResidue]:
        try:
            payload = self.http.fetch_json(
                f"{GPCRDB_SERVICES}/residues/{entry_name}/",
                cache_namespace="gpcrdb",
                cache_key=f"residues:{entry_name}",
            )
        except Exception:
            return []
        if not isinstance(payload, list):
            return []
        residues: list[GPCRdbResidue] = []
        for item in payload:
            if not isinstance(item, dict):
                continue
            sequence_number = item.get("sequence_number")
            amino_acid = item.get("amino_acid")
            if not isinstance(sequence_number, int) or not isinstance(amino_acid, str) or not amino_acid:
                continue
            residues.append(
                GPCRdbResidue(
                    sequence_number=sequence_number,
                    amino_acid=amino_acid,
                    protein_segment=_clean_text(item.get("protein_segment")),
                    display_generic_number=_clean_text(item.get("display_generic_number")),
                )
            )
        return residues

    def _fetch_family_path(self, family_slug: str) -> list[str]:
        path: list[str] = []
        seen: set[str] = set()
        slug: str | None = family_slug
        while slug and slug not in seen:
            seen.add(slug)
            try:
                payload = self.http.fetch_json(
                    f"{GPCRDB_SERVICES}/proteinfamily/{slug}/",
                    cache_namespace="gpcrdb",
                    cache_key=f"family:{slug}",
                )
            except Exception:
                break
            if not isinstance(payload, dict):
                break
            name = _clean_text(payload.get("name"))
            if name and name.lower() != "parent family":
                path.append(name)
            parent = payload.get("parent")
            slug = parent.get("slug") if isinstance(parent, dict) else None
        return list(reversed(path))


def _segments_from_residues(residues: list[GPCRdbResidue]) -> list[GPCRdbSegment]:
    segments: list[GPCRdbSegment] = []
    current_name: str | None = None
    current_residues: list[GPCRdbResidue] = []
    for residue in sorted(residues, key=lambda item: item.sequence_number):
        name = residue.protein_segment or "Unknown"
        if current_name is None:
            current_name = name
        if name != current_name:
            segments.append(_segment_from_group(current_name, current_residues))
            current_name = name
            current_residues = []
        current_residues.append(residue)
    if current_name is not None and current_residues:
        segments.append(_segment_from_group(current_name, current_residues))
    return segments


def _segment_from_group(name: str, residues: list[GPCRdbResidue]) -> GPCRdbSegment:
    generic_numbers = [residue.display_generic_number for residue in residues if residue.display_generic_number]
    return GPCRdbSegment(
        name=name,
        start=min(residue.sequence_number for residue in residues),
        end=max(residue.sequence_number for residue in residues),
        generic_start=generic_numbers[0] if generic_numbers else None,
        generic_end=generic_numbers[-1] if generic_numbers else None,
    )


def _conserved_motifs(residues: list[GPCRdbResidue]) -> list[GPCRdbMotif]:
    by_generic: dict[str, GPCRdbResidue] = {}
    for residue in residues:
        generic = _canonical_generic(residue.display_generic_number)
        if generic:
            by_generic[generic] = residue
    motifs = [
        ("DRY / E/DRY activation motif", ["3.49", "3.50", "3.51"]),
        ("CWxP / toggle-switch motif", ["6.47", "6.48", "6.49", "6.50"]),
        ("NPxxY / helix-7 motif", ["7.49", "7.50", "7.51", "7.52", "7.53"]),
        ("Na+ pocket / class-A D2.50 position", ["2.50"]),
        ("Class-A N1.50 position", ["1.50"]),
    ]
    detected: list[GPCRdbMotif] = []
    for name, generics in motifs:
        members = [by_generic[generic] for generic in generics if generic in by_generic]
        if not members:
            continue
        detected.append(
            GPCRdbMotif(
                name=name,
                positions=[item.sequence_number for item in members],
                generic_numbers=[item.display_generic_number or "" for item in members],
                sequence="".join(item.amino_acid for item in members),
                status="detected" if len(members) == len(generics) else "partial",
            )
        )
    return detected


def _canonical_generic(value: str | None) -> str | None:
    if not value:
        return None
    text = value.strip()
    if "x" in text:
        text = text.split("x", 1)[0]
    return text


def _clean_text(value: Any) -> str | None:
    if value is None:
        return None
    text = html.unescape(str(value)).strip()
    return text or None
