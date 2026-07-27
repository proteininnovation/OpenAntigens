from __future__ import annotations

from typing import Any

from .http import HttpClient
from .models import DomainAnnotation


INTERPRO_API_BASE = "https://www.ebi.ac.uk/interpro/api"


class InterProClient:
    def __init__(self, http: HttpClient) -> None:
        self.http = http

    def fetch_annotations(self, accession: str) -> list[DomainAnnotation]:
        annotations: list[DomainAnnotation] = []
        for database in ("interpro", "pfam"):
            annotations.extend(self._fetch_database_annotations(database, accession))
        return self._dedupe_annotations(annotations)

    def _fetch_database_annotations(self, database: str, accession: str) -> list[DomainAnnotation]:
        url = f"{INTERPRO_API_BASE}/entry/{database}/protein/uniprot/{accession}?page_size=200"
        annotations: list[DomainAnnotation] = []
        while url:
            payload = self.http.fetch_json(
                url,
                headers={"Accept": "application/json"},
                cache_namespace="interpro",
            )
            for result in payload.get("results", []):
                annotations.extend(self._parse_result(database, result))
            next_url = payload.get("next")
            url = next_url if isinstance(next_url, str) and next_url else ""
        return annotations

    def _parse_result(self, database: str, result: dict[str, Any]) -> list[DomainAnnotation]:
        metadata = result.get("metadata", {})
        accession = str(metadata.get("accession") or "")
        if not accession:
            return []
        name = str(metadata.get("name") or accession)
        entry_type = self._normalize_type(str(metadata.get("type") or "unknown"))
        source_database = str(
            metadata.get("source_database")
            or metadata.get("member_database")
            or database
        ).upper()
        integrated = metadata.get("integrated")
        integrated_accession = None
        integrated_name = None
        if isinstance(integrated, dict):
            integrated_accession = integrated.get("accession")
            integrated_name = integrated.get("name")

        annotations: list[DomainAnnotation] = []
        proteins = result.get("proteins", [])
        if not proteins:
            proteins = [result]
        for protein in proteins:
            signature = protein.get("signature", {})
            locations = protein.get("entry_protein_locations") or protein.get("protein_locations") or []
            if not locations and result.get("entry_protein_locations"):
                locations = result["entry_protein_locations"]
            local_integrated_accession = integrated_accession
            local_integrated_name = integrated_name
            if isinstance(signature, dict):
                signature_entry = signature.get("entry")
                if isinstance(signature_entry, dict):
                    local_integrated_accession = local_integrated_accession or signature_entry.get("accession")
                    local_integrated_name = local_integrated_name or signature_entry.get("name")
            for location in locations:
                representative = bool(location.get("representative"))
                for fragment in location.get("fragments", []):
                    start = self._location_value(fragment.get("start"))
                    end = self._location_value(fragment.get("end"))
                    if start is None or end is None:
                        continue
                    annotations.append(
                        DomainAnnotation(
                            accession=accession,
                            name=name,
                            type=entry_type,
                            source_database=source_database,
                            start=start,
                            end=end,
                            representative=representative,
                            integrated_accession=local_integrated_accession,
                            integrated_name=local_integrated_name,
                            metadata={
                                "source_database_raw": database,
                                "member_database": metadata.get("member_database"),
                                "signature_accession": signature.get("accession") if isinstance(signature, dict) else None,
                            },
                        )
                    )
        return annotations

    def _location_value(self, value: Any) -> int | None:
        if isinstance(value, int):
            return value
        if isinstance(value, str) and value.isdigit():
            return int(value)
        if isinstance(value, dict):
            raw = value.get("value")
            if isinstance(raw, int):
                return raw
            if isinstance(raw, str) and raw.isdigit():
                return int(raw)
        return None

    def _normalize_type(self, value: str) -> str:
        return value.strip().lower().replace(" ", "_")

    def _dedupe_annotations(self, annotations: list[DomainAnnotation]) -> list[DomainAnnotation]:
        deduped: dict[tuple[str, str, int, int, str], DomainAnnotation] = {}
        for annotation in annotations:
            key = (
                annotation.source_database,
                annotation.accession,
                annotation.start,
                annotation.end,
                annotation.type,
            )
            current = deduped.get(key)
            if current is None or (annotation.representative and not current.representative):
                deduped[key] = annotation
        return sorted(
            deduped.values(),
            key=lambda item: (item.start, item.end, item.source_database, item.accession),
        )
