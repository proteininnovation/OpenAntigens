from __future__ import annotations

import json
import re
import time
from dataclasses import dataclass
from typing import Any
from urllib.parse import quote

from .exceptions import ExternalServiceError
from .http import HttpClient
from .models import AnalysisNote, TargetRecord


EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# HTTP statuses worth retrying: rate limiting (429) and transient server errors.
# HttpClient formats these as "HTTP error for {url}: {code}" and connection
# failures as "Connection error for {url}: {reason}".
_RETRYABLE_HTTP_CODES = ("429", "500", "502", "503", "504")


def _is_retryable_external_error(exc: ExternalServiceError) -> bool:
    text = str(exc)
    if text.startswith("Connection error"):
        return True
    return any(text.endswith(f": {code}") for code in _RETRYABLE_HTTP_CODES)


@dataclass(slots=True)
class RefSeqProtein:
    target: TargetRecord
    notes: list[AnalysisNote]
    title: str


class RefSeqClient:
    def __init__(self, http: HttpClient) -> None:
        self.http = http

    def fetch_canonical_protein(
        self,
        *,
        gene_symbol: str,
        organism: str,
        gene_id: str | int | None = None,
    ) -> RefSeqProtein | None:
        ids: list[str] = []
        for term in self._candidate_terms(gene_symbol=gene_symbol, organism=organism, gene_id=gene_id):
            search_url = (
                f"{EUTILS_BASE}/esearch.fcgi?db=protein&retmode=json&retmax=20&term={quote(term)}"
            )
            search_payload = self._fetch_json_with_retry(search_url, cache_namespace="refseq_esearch")
            ids = search_payload.get("esearchresult", {}).get("idlist", [])
            if ids:
                break
        if not ids:
            return None

        summary_url = (
            f"{EUTILS_BASE}/esummary.fcgi?db=protein&retmode=json&id={quote(','.join(ids))}"
        )
        summary_payload = self._fetch_json_with_retry(summary_url, cache_namespace="refseq_esummary")
        candidates = []
        for uid in summary_payload.get("result", {}).get("uids", []):
            item = summary_payload["result"].get(uid, {})
            accession = item.get("caption")
            title = item.get("title", "")
            if accession:
                candidates.append((self._rank_candidate(accession, title), accession, title))
        if not candidates:
            return None

        _, accession, title = sorted(candidates, key=lambda item: item[0])[0]
        fasta_url = (
            f"{EUTILS_BASE}/efetch.fcgi?db=protein&id={quote(accession)}&rettype=fasta&retmode=text"
        )
        fasta_text = self._fetch_text_with_retry(
            fasta_url,
            cache_namespace="refseq_fasta",
            suffix=".fasta",
            headers={"Accept": "text/plain"},
        )
        header, sequence = self._parse_fasta(fasta_text)
        notes = [
            AnalysisNote(
                severity="info",
                message=f"Resolved from RefSeq canonical protein candidate {accession}.",
                source="RefSeq",
            )
        ]
        taxon_id = self._taxon_id_for_organism(organism)
        return RefSeqProtein(
            target=TargetRecord(
                accession=accession,
                entry_name=accession,
                gene_symbol=gene_symbol,
                protein_name=title or accession,
                organism=organism,
                taxon_id=taxon_id,
                sequence=sequence,
            ),
            notes=notes,
            title=title,
        )

    def _rank_candidate(self, accession: str, title: str) -> tuple[int, int, str]:
        lowered = title.lower()
        select_rank = 0 if "mane select" in lowered or "refseq select" in lowered else 1
        prefix_rank = 0 if accession.startswith("NP_") else 1 if accession.startswith("XP_") else 2
        isoform_rank = self._isoform_rank(title)
        return (select_rank, prefix_rank, isoform_rank, accession)

    def _isoform_rank(self, title: str) -> int:
        match = re.search(r"isoform\s+X?(\d+)", title, flags=re.IGNORECASE)
        if match:
            return int(match.group(1))
        return 0

    def _parse_fasta(self, fasta_text: str) -> tuple[str, str]:
        lines = [line.strip() for line in fasta_text.splitlines() if line.strip()]
        if not lines or not lines[0].startswith(">"):
            raise ValueError("Unexpected RefSeq FASTA response.")
        return lines[0][1:], "".join(lines[1:])

    def _candidate_terms(
        self,
        *,
        gene_symbol: str,
        organism: str,
        gene_id: str | int | None,
    ) -> list[str]:
        terms: list[str] = []
        gene_id_text = str(gene_id) if gene_id is not None else None
        if gene_id_text:
            if organism == "Homo sapiens":
                terms.append(f"{gene_id_text}[Gene ID] AND srcdb_refseq[prop] AND MANE Select[keyword]")
                terms.append(f"{gene_id_text}[Gene ID] AND srcdb_refseq[prop] AND refseq_select[filter]")
            elif organism == "Mus musculus":
                terms.append(f"{gene_id_text}[Gene ID] AND srcdb_refseq[prop] AND refseq_select[filter]")
            terms.append(f"{gene_id_text}[Gene ID] AND srcdb_refseq[prop]")
            terms.append(f"{gene_id_text}[Gene ID] AND {organism}[Organism] AND srcdb_refseq[prop]")
        if organism == "Homo sapiens":
            terms.append(
                f"{gene_symbol}[Gene Name] AND {organism}[Organism] AND srcdb_refseq[prop] AND MANE Select[keyword]"
            )
        if organism in {"Homo sapiens", "Mus musculus"}:
            terms.append(
                f"{gene_symbol}[Gene Name] AND {organism}[Organism] AND srcdb_refseq[prop] AND refseq_select[filter]"
            )
        terms.append(f"{gene_symbol}[Gene Name] AND {organism}[Organism] AND srcdb_refseq[prop]")
        return list(dict.fromkeys(terms))

    def _fetch_json_with_retry(self, url: str, *, cache_namespace: str) -> Any:
        delay_seconds = 1.0
        for attempt in range(3):
            try:
                return self.http.fetch_json(url, cache_namespace=cache_namespace)
            except ExternalServiceError as exc:
                if not _is_retryable_external_error(exc) or attempt == 2:
                    raise
                time.sleep(delay_seconds)
                delay_seconds *= 2.0
        raise RuntimeError("Unreachable retry loop.")

    def _fetch_text_with_retry(
        self,
        url: str,
        *,
        cache_namespace: str,
        suffix: str,
        headers: dict[str, str] | None = None,
    ) -> str:
        delay_seconds = 1.0
        for attempt in range(3):
            try:
                return self.http.fetch_text(
                    url,
                    cache_namespace=cache_namespace,
                    suffix=suffix,
                    headers=headers,
                )
            except ExternalServiceError as exc:
                if not _is_retryable_external_error(exc) or attempt == 2:
                    raise
                time.sleep(delay_seconds)
                delay_seconds *= 2.0
        raise RuntimeError("Unreachable retry loop.")

    def _taxon_id_for_organism(self, organism: str) -> int | None:
        return {
            "Homo sapiens": 9606,
            "Mus musculus": 10090,
            "Macaca fascicularis": 9541,
        }.get(organism)
