from __future__ import annotations

import csv
from io import StringIO
from urllib.parse import quote

from .exceptions import ExternalServiceError
from .http import HttpClient
from .models import ComplexPortalComplex, ComplexPortalParticipant


COMPLEX_PORTAL_BASE = "https://www.ebi.ac.uk/intact/complex-ws"
COMPLEX_PORTAL_HUMAN_COMPLEXTAB = "https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/9606.tsv"
COMPLEX_PORTAL_HUMAN_PREDICTED_COMPLEXTAB = "https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/9606_predicted.tsv"
COMPLEX_PORTAL_MOUSE_COMPLEXTAB = "https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/10090.tsv"
COMPLEX_PORTAL_BULK_CATALOGUES = {
    9606: (COMPLEX_PORTAL_HUMAN_COMPLEXTAB, "Homo sapiens"),
    10090: (COMPLEX_PORTAL_MOUSE_COMPLEXTAB, "Mus musculus"),
}
COMPLEXTAB_COLUMNS = (
    "#Complex ac",
    "Recommended name",
    "Aliases for complex",
    "Taxonomy identifier",
    "Identifiers (and stoichiometry) of molecules in complex",
    "Evidence Code",
    "Experimental evidence",
    "Go Annotations",
    "Cross references",
    "Description",
    "Complex properties",
    "Complex assembly",
    "Ligand",
    "Disease",
    "Agonist",
    "Antagonist",
    "Comment",
    "Source",
    "Expanded participant list",
)


class ComplexPortalClient:
    def __init__(self, http: HttpClient) -> None:
        self.http = http
        self._bulk_complex_indexes: dict[tuple[int, bool], dict[str, list[ComplexPortalComplex]]] = {}
        self._bulk_complex_index_errors: dict[tuple[int, bool], Exception] = {}

    def fetch_human_complexes_for_target(
        self,
        *,
        accession: str,
        include_predicted: bool = False,
    ) -> list[ComplexPortalComplex]:
        """Return curated human Complex Portal records from one cached ComplexTab download."""
        return self._bulk_complexes_for_accession(
            taxon_id=9606,
            accession=accession,
            include_predicted=include_predicted,
        )

    def _bulk_complexes_for_accession(
        self,
        *,
        taxon_id: int,
        accession: str,
        include_predicted: bool = False,
    ) -> list[ComplexPortalComplex]:
        if include_predicted and taxon_id != 9606:
            raise ValueError("Predicted ComplexTab lookup is only available for human targets")
        normalized_accession = self._normalize_identifier(accession) or accession.upper()
        complexes = list(self._bulk_complex_index(taxon_id=taxon_id, predicted_complex=False).get(normalized_accession, []))
        if include_predicted:
            complexes.extend(self._bulk_complex_index(taxon_id=taxon_id, predicted_complex=True).get(normalized_accession, []))
        return sorted(complexes, key=lambda item: (item.predicted_complex, -(item.confidence_score or 0), item.name))

    def _bulk_complex_index(
        self,
        *,
        taxon_id: int,
        predicted_complex: bool,
    ) -> dict[str, list[ComplexPortalComplex]]:
        key = (taxon_id, predicted_complex)
        if key in self._bulk_complex_index_errors:
            raise self._bulk_complex_index_errors[key]
        if key not in self._bulk_complex_indexes:
            url = COMPLEX_PORTAL_HUMAN_PREDICTED_COMPLEXTAB if predicted_complex else COMPLEX_PORTAL_BULK_CATALOGUES[taxon_id][0]
            try:
                text = self.http.fetch_text(
                    url,
                    cache_namespace="complex_portal_complextab",
                    suffix=".tsv",
                )
                self._bulk_complex_indexes[key] = self._parse_complextab(
                    text,
                    taxon_id=taxon_id,
                    predicted_complex=predicted_complex,
                )
            except Exception as exc:
                self._bulk_complex_index_errors[key] = exc
                raise
        return self._bulk_complex_indexes[key]

    def _parse_complextab(
        self,
        text: str,
        *,
        taxon_id: int,
        predicted_complex: bool,
    ) -> dict[str, list[ComplexPortalComplex]]:
        index: dict[str, list[ComplexPortalComplex]] = {}
        reader = csv.DictReader(StringIO(text), delimiter="\t")
        if tuple(reader.fieldnames or ()) != COMPLEXTAB_COLUMNS:
            raise ExternalServiceError("Unexpected ComplexTab header")
        row_count = 0
        for row in reader:
            row_count += 1
            complex_ac = self._complextab_value(row, "#Complex ac")
            if not complex_ac:
                continue
            row_taxon_id = self._complextab_value(row, "Taxonomy identifier")
            if row_taxon_id != str(taxon_id):
                raise ExternalServiceError("Unexpected ComplexTab taxonomy")
            evidence_code, evidence_description = self._parse_complextab_evidence(
                self._complextab_value(row, "Evidence Code")
            )
            participants = self._parse_complextab_participants(
                self._complextab_value(row, "Expanded participant list")
            )
            if not participants:
                continue
            complex_item = ComplexPortalComplex(
                complex_ac=complex_ac,
                name=self._complextab_value(row, "Recommended name") or complex_ac,
                species=f"{COMPLEX_PORTAL_BULK_CATALOGUES[taxon_id][1]}; {row_taxon_id}",
                predicted_complex=predicted_complex,
                evidence_code=evidence_code,
                evidence_description=evidence_description,
                complex_assemblies=self._split_complextab_value(self._complextab_value(row, "Complex assembly")),
                properties=self._split_complextab_value(self._complextab_value(row, "Complex properties")),
                participants=participants,
            )
            for participant in participants:
                accession = self._normalize_identifier(participant.identifier)
                if accession:
                    index.setdefault(accession, []).append(complex_item)
        if not row_count:
            raise ExternalServiceError("ComplexTab catalogue contains no data rows")
        if not index:
            raise ExternalServiceError("ComplexTab catalogue contains no participants")
        return index

    @staticmethod
    def _complextab_value(row: dict[str | None, str | list[str] | None], column: str) -> str:
        value = row.get(column)
        if not isinstance(value, str):
            return ""
        value = value.strip()
        return "" if value == "-" else value

    @staticmethod
    def _split_complextab_value(value: str) -> list[str]:
        return [item.strip() for item in value.split("|") if item.strip() and item.strip() != "-"]

    @staticmethod
    def _parse_complextab_evidence(value: str) -> tuple[str | None, str | None]:
        if not value:
            return None, None
        if value.endswith(")") and "(" in value:
            code, description = value.split("(", 1)
            return code.strip() or None, description[:-1].strip() or None
        return value, None

    def _parse_complextab_participants(self, value: str) -> list[ComplexPortalParticipant]:
        participants: list[ComplexPortalParticipant] = []
        for token in self._split_complextab_value(value):
            identifier = token
            stoichiometry = None
            if token.endswith(")") and "(" in token:
                identifier, stoichiometry = token.rsplit("(", 1)
                identifier = identifier.strip()
                stoichiometry = stoichiometry[:-1].strip() or None
            if identifier:
                participants.append(
                    ComplexPortalParticipant(
                        identifier=identifier,
                        name=identifier,
                        interactor_type="protein",
                        stoichiometry=stoichiometry,
                    )
                )
        return participants

    def fetch_complexes_for_target(
        self,
        *,
        accession: str,
        gene_symbol: str | None,
        taxon_id: int | None,
        entry_name: str | None = None,
        max_pages: int = 2,
        include_predicted: bool = False,
    ) -> list[ComplexPortalComplex]:
        if taxon_id in COMPLEX_PORTAL_BULK_CATALOGUES:
            return self._bulk_complexes_for_accession(
                taxon_id=taxon_id,
                accession=accession,
                include_predicted=include_predicted,
            )

        queries = [accession]
        if gene_symbol and gene_symbol not in queries:
            queries.append(gene_symbol)
        if entry_name:
            base_name = entry_name.split("_", 1)[0]
            if base_name and base_name not in queries:
                queries.append(base_name)

        complexes: dict[str, ComplexPortalComplex] = {}
        for query in queries:
            for page in range(max_pages):
                payload = self._search(query=query, page=page)
                elements = payload.get("elements", [])
                if not elements:
                    break
                for element in elements:
                    if not self._search_hit_matches_target(element, accession=accession, gene_symbol=gene_symbol, taxon_id=taxon_id):
                        continue
                    complex_ac = element.get("complexAC")
                    if not complex_ac or complex_ac in complexes:
                        continue
                    detail = self._fetch_complex_detail(complex_ac)
                    if not self._detail_matches_target(detail, accession=accession, gene_symbol=gene_symbol, taxon_id=taxon_id):
                        continue
                    complexes[complex_ac] = self._build_complex(detail)
        return sorted(complexes.values(), key=lambda item: (item.predicted_complex, -(item.confidence_score or 0), item.name))

    def _search(self, *, query: str, page: int) -> dict:
        url = f"{COMPLEX_PORTAL_BASE}/search/*?query={quote(query)}&format=json&page={page}"
        return self.http.fetch_json(url, cache_namespace="complex_portal_search")

    def _fetch_complex_detail(self, complex_ac: str) -> dict:
        url = f"{COMPLEX_PORTAL_BASE}/complex/{quote(complex_ac)}"
        return self.http.fetch_json(url, cache_namespace="complex_portal_detail")

    def _search_hit_matches_target(
        self,
        element: dict,
        *,
        accession: str,
        gene_symbol: str | None,
        taxon_id: int | None,
    ) -> bool:
        if taxon_id is not None and not self._species_matches(element.get("organismName"), taxon_id):
            return False
        for interactor in element.get("interactors", []):
            if self._participant_matches_target(interactor.get("identifier"), interactor.get("name"), accession, gene_symbol):
                return True
        return False

    def _detail_matches_target(
        self,
        detail: dict,
        *,
        accession: str,
        gene_symbol: str | None,
        taxon_id: int | None,
    ) -> bool:
        if taxon_id is not None and not self._species_matches(detail.get("species"), taxon_id):
            return False
        for participant in detail.get("participants", []):
            if self._participant_matches_target(participant.get("identifier"), participant.get("name"), accession, gene_symbol):
                return True
        return False

    def _participant_matches_target(
        self,
        identifier: str | None,
        name: str | None,
        accession: str,
        gene_symbol: str | None,
    ) -> bool:
        candidate_ids = {self._normalize_identifier(identifier)}
        if name:
            candidate_ids.add(name.upper())
        candidate_ids.discard(None)
        return accession.upper() in candidate_ids or (gene_symbol and gene_symbol.upper() in candidate_ids)

    def _normalize_identifier(self, identifier: str | None) -> str | None:
        if not identifier:
            return None
        token = identifier.split(":", 1)[-1]
        token = token.split("-", 1)[0]
        return token.upper()

    def _species_matches(self, species_text: str | None, taxon_id: int) -> bool:
        if not species_text:
            return False
        return str(taxon_id) in str(species_text)

    def _build_complex(self, detail: dict) -> ComplexPortalComplex:
        evidence = detail.get("evidenceType") or {}
        participants: list[ComplexPortalParticipant] = []
        for participant in detail.get("participants", []):
            binding_regions: list[str] = []
            for feature in participant.get("linkedFeatures", []) or []:
                if str(feature.get("featureType", "")).lower() != "binding region":
                    continue
                ranges = feature.get("ranges") or []
                for range_text in ranges:
                    if range_text not in binding_regions:
                        binding_regions.append(range_text)
            participants.append(
                ComplexPortalParticipant(
                    identifier=participant.get("identifier", ""),
                    name=participant.get("name", participant.get("identifier", "")),
                    description=participant.get("description"),
                    stoichiometry=participant.get("stochiometry"),
                    bio_role=participant.get("bioRole"),
                    interactor_type=participant.get("interactorType"),
                    binding_regions=binding_regions,
                )
            )
        return ComplexPortalComplex(
            complex_ac=detail.get("complexAc", ""),
            name=detail.get("name", ""),
            species=detail.get("species", ""),
            predicted_complex=bool(detail.get("predictedComplex", False)),
            evidence_code=evidence.get("identifier"),
            evidence_description=evidence.get("description"),
            confidence_score=evidence.get("confidenceScore"),
            complex_assemblies=list(detail.get("complexAssemblies") or []),
            properties=list(detail.get("properties") or []),
            participants=participants,
        )
