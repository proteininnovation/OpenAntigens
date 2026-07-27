from __future__ import annotations

from dataclasses import dataclass
import re
from typing import Any
from urllib.parse import quote

from .exceptions import ResolutionError
from .http import HttpClient
from .models import AnalysisNote, Feature, ExperimentalConstruct, Region, TargetRecord


UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb"


@dataclass(slots=True)
class TargetResolution:
    target: TargetRecord
    entry: dict[str, Any]
    notes: list[AnalysisNote]
    query_type: str


class UniProtClient:
    def __init__(self, http: HttpClient) -> None:
        self.http = http

    def resolve_target(self, query: str, target_species_taxon: int = 9606) -> TargetResolution:
        cleaned = query.strip()
        if "_" in cleaned:
            query_type = "entry_name"
            entry = self._fetch_entry(cleaned)
            notes = [AnalysisNote(severity="info", message="Resolved from UniProt entry name.", source="UniProt")]
            return TargetResolution(target=self._target_from_entry(entry), entry=entry, notes=notes, query_type=query_type)
        if self._looks_like_accession(cleaned):
            query_type = "accession"
            entry = self._fetch_entry(cleaned)
            notes = [AnalysisNote(severity="info", message="Resolved from UniProt accession.", source="UniProt")]
            return TargetResolution(target=self._target_from_entry(entry), entry=entry, notes=notes, query_type=query_type)
        query_type = "gene_symbol"
        search_query = f'(gene_exact:{cleaned}) AND (organism_id:{target_species_taxon}) AND (reviewed:true)'
        results = self.search(search_query, size=5)
        entries = results.get("results", [])
        if not entries:
            raise ResolutionError(f"No reviewed UniProt record found for gene symbol '{query}'.")
        if len(entries) > 1:
            accessions = ", ".join(item.get("primaryAccession", "?") for item in entries[:5])
            raise ResolutionError(
                f"Ambiguous gene symbol '{query}'. Candidate reviewed accessions: {accessions}"
            )
        accession = entries[0]["primaryAccession"]
        entry = self._fetch_entry(accession)
        notes = [AnalysisNote(severity="info", message="Resolved from official gene symbol.", source="UniProt")]
        return TargetResolution(target=self._target_from_entry(entry), entry=entry, notes=notes, query_type=query_type)

    def search(self, query: str, *, size: int = 10) -> dict[str, Any]:
        url = f"{UNIPROT_BASE}/search?query={quote(query)}&format=json&size={size}"
        return self.http.fetch_json(url, cache_namespace="uniprot_search")

    def fetch_entry_by_name(self, entry_name: str) -> dict[str, Any] | None:
        try:
            return self._fetch_entry(entry_name)
        except Exception:
            return None

    def get_features(self, entry: dict[str, Any]) -> list[Feature]:
        features: list[Feature] = []
        for feature in entry.get("features", []):
            location = feature.get("location", {})
            start = self._location_value(location.get("start"))
            end = self._location_value(location.get("end"))
            if not start or not end:
                continue
            description = feature.get("description")
            raw_type = feature.get("type", "UNKNOWN")
            metadata = {"raw_type": raw_type}
            if "ligand" in feature:
                metadata["ligand"] = feature["ligand"]
            if "featureId" in feature:
                metadata["feature_id"] = feature["featureId"]
            features.append(
                Feature(
                    type=self._normalize_feature_type(raw_type),
                    start=start,
                    end=end,
                    description=description,
                    metadata=metadata,
                )
            )
        return sorted(features, key=lambda item: (item.start, item.end, item.type))

    def get_comment_texts(self, entry: dict[str, Any], comment_type: str) -> list[str]:
        texts: list[str] = []
        for comment in entry.get("comments", []):
            if str(comment.get("commentType", "")).upper() != comment_type.upper():
                continue
            texts.extend(self._extract_comment_texts(comment))
        return [text for text in texts if text]

    def get_experimental_constructs(self, entry: dict[str, Any]) -> list[ExperimentalConstruct]:
        constructs: list[ExperimentalConstruct] = []
        for reference in entry.get("uniProtKBCrossReferences", []):
            if reference.get("database") != "PDB":
                continue
            properties = {item.get("key"): item.get("value") for item in reference.get("properties", [])}
            chains = self._parse_chain_ranges(properties.get("Chains", ""))
            if not chains:
                continue
            constructs.append(
                ExperimentalConstruct(
                    pdb_id=reference.get("id", ""),
                    method=properties.get("Method"),
                    resolution=properties.get("Resolution"),
                    chains=chains,
                )
            )
        return constructs

    def fetch_same_name_species_match(self, base_name: str, species_suffix: str) -> dict[str, Any] | None:
        return self.fetch_entry_by_name(f"{base_name}_{species_suffix}")

    def _fetch_entry(self, accession_or_id: str) -> dict[str, Any]:
        url = f"{UNIPROT_BASE}/{quote(accession_or_id)}.json"
        return self.http.fetch_json(url, cache_namespace="uniprot_entry")

    def _target_from_entry(self, entry: dict[str, Any]) -> TargetRecord:
        sequence_block = entry.get("sequence", {})
        genes = entry.get("genes", [])
        gene_symbol = None
        if genes:
            gene_name = genes[0].get("geneName")
            if isinstance(gene_name, dict):
                gene_symbol = gene_name.get("value")
        protein_description = entry.get("proteinDescription", {})
        protein_name = (
            protein_description.get("recommendedName", {})
            .get("fullName", {})
            .get("value")
            or entry.get("uniProtkbId", entry.get("primaryAccession", "Unknown protein"))
        )
        return TargetRecord(
            accession=entry.get("primaryAccession", ""),
            entry_name=entry.get("uniProtkbId", ""),
            gene_symbol=gene_symbol,
            protein_name=protein_name,
            organism=entry.get("organism", {}).get("scientificName", ""),
            taxon_id=entry.get("organism", {}).get("taxonId"),
            sequence=sequence_block.get("value", ""),
            canonical_isoform_id=self._canonical_isoform_id(entry),
            alternative_names=self._alternative_protein_names(protein_description),
            gene_synonyms=self._gene_synonyms(genes),
        )

    def _alternative_protein_names(self, protein_description: dict[str, Any]) -> list[str]:
        names: list[str] = []

        def add_name(value: Any) -> None:
            if isinstance(value, dict):
                text = value.get("value")
            else:
                text = value
            cleaned = str(text or "").strip()
            if cleaned and cleaned not in names:
                names.append(cleaned)

        def collect_name_block(block: Any) -> None:
            if not isinstance(block, dict):
                return
            add_name(block.get("fullName"))
            for item in block.get("shortNames") or []:
                add_name(item)

        for block in protein_description.get("alternativeNames") or []:
            collect_name_block(block)
        for container_key in ("contains", "includes"):
            for component in protein_description.get(container_key) or []:
                collect_name_block(component.get("recommendedName"))
                for block in component.get("alternativeNames") or []:
                    collect_name_block(block)
        for item in protein_description.get("cdAntigenNames") or []:
            cd_name = str(item.get("value") if isinstance(item, dict) else item).strip()
            if not cd_name:
                continue
            add_name(cd_name)
            add_name(f"CD antigen {cd_name}")
            all_labels = " ".join(names).lower()
            if "heavy chain" in all_labels and not cd_name.lower().endswith("hc"):
                add_name(f"{cd_name}hc")
        return names

    def _gene_synonyms(self, genes: list[dict[str, Any]]) -> list[str]:
        synonyms: list[str] = []
        for gene in genes:
            for key in ("synonyms", "orderedLocusNames", "orfNames"):
                for item in gene.get(key) or []:
                    if isinstance(item, dict):
                        value = item.get("value")
                    else:
                        value = item
                    cleaned = str(value or "").strip()
                    if cleaned and cleaned not in synonyms:
                        synonyms.append(cleaned)
        return synonyms

    def _canonical_isoform_id(self, entry: dict[str, Any]) -> str | None:
        accession = str(entry.get("primaryAccession", "")).strip()
        if not accession:
            return None
        for comment in entry.get("comments", []):
            if str(comment.get("commentType", "")).upper() != "ALTERNATIVE PRODUCTS":
                continue
            for isoform in comment.get("isoforms", []):
                if str(isoform.get("isoformSequenceStatus", "")).lower() != "displayed":
                    continue
                isoform_ids = isoform.get("isoformIds") or []
                for isoform_id in isoform_ids:
                    cleaned = str(isoform_id or "").strip()
                    if cleaned:
                        return cleaned
        return accession

    def _parse_chain_ranges(self, chains: str) -> list[Region]:
        regions: list[Region] = []
        for part in chains.split(","):
            block = part.strip()
            if "=" not in block:
                continue
            chain_name, raw_ranges = block.split("=", 1)
            for raw_range in raw_ranges.split("/"):
                text = raw_range.strip()
                if "-" not in text:
                    continue
                start_text, end_text = text.split("-", 1)
                if not start_text.isdigit() or not end_text.isdigit():
                    continue
                regions.append(
                    Region(
                        start=int(start_text),
                        end=int(end_text),
                        label=f"Chain {chain_name.strip()}",
                        source="PDB",
                        metadata={"chain": chain_name.strip()},
                    )
                )
        return regions

    def _looks_like_accession(self, value: str) -> bool:
        patterns = (
            r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$",
            r"^[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]$",
            r"^[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9][A-Z][A-Z0-9]{2}[0-9]$",
        )
        return any(re.fullmatch(pattern, value) for pattern in patterns)

    def _location_value(self, block: Any) -> int | None:
        if isinstance(block, dict):
            value = block.get("value")
            if isinstance(value, int):
                return value
            if isinstance(value, str) and value.isdigit():
                return int(value)
        if isinstance(block, str) and block.isdigit():
            return int(block)
        return None

    def _normalize_feature_type(self, raw_type: str) -> str:
        normalized = raw_type.strip().lower().replace(" ", "_")
        aliases = {
            "transmembrane": "TRANSMEM",
            "topological_domain": "TOPO_DOM",
            "signal": "SIGNAL",
            "domain": "DOMAIN",
            "region": "REGION",
            "binding_site": "BINDING",
            "site": "SITE",
            "disulfide_bond": "DISULFID",
            "glycosylation": "CARBOHYD",
            "modified_residue": "MOD_RES",
        }
        return aliases.get(normalized, raw_type.upper().replace(" ", "_"))

    def _extract_comment_texts(self, payload: Any) -> list[str]:
        collected: list[str] = []
        if isinstance(payload, dict):
            value = payload.get("value")
            if isinstance(value, str) and value.strip():
                collected.append(value.strip())
            for key in ("texts", "subcellularLocations", "molecule", "note", "interactions"):
                nested = payload.get(key)
                if nested is not None:
                    collected.extend(self._extract_comment_texts(nested))
        elif isinstance(payload, list):
            for item in payload:
                collected.extend(self._extract_comment_texts(item))
        elif isinstance(payload, str) and payload.strip():
            collected.append(payload.strip())
        return collected
