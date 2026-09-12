from __future__ import annotations

import csv
import json
from pathlib import Path
from urllib.parse import quote

from .cache import FileCache
from .http import HttpClient


UNIPROT_SEARCH_BASE = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_ENTRY_BASE = "https://rest.uniprot.org/uniprotkb"


def prefetch_reviewed_taxa(
    *,
    cache_dir: Path,
    taxa: tuple[int, ...],
    page_size: int = 200,
    max_pages: int | None = None,
) -> dict[str, int]:
    cache = FileCache(cache_dir)
    http = HttpClient(cache=cache, timeout=60)
    summary = {
        "pages": 0,
        "entries": 0,
        "written_entry_cache": 0,
    }
    for taxon_id in taxa:
        taxon_stats = _prefetch_taxon(
            http=http,
            cache=cache,
            taxon_id=taxon_id,
            page_size=page_size,
            max_pages=max_pages,
        )
        summary["pages"] += taxon_stats["pages"]
        summary["entries"] += taxon_stats["entries"]
        summary["written_entry_cache"] += taxon_stats["written_entry_cache"]
    return summary


def prefetch_targets_from_tsv(
    *,
    tsv_path: str | Path,
    cache_dir: Path,
    chunk_size: int = 50,
) -> dict[str, int]:
    """Bulk-prefetch UniProt entries listed in a target TSV.

    This is the deployment-friendly path for fresh snapshots: it uses the
    locally curated target universe as the only local biological input, then
    retrieves just those current UniProt records in chunks rather than mirroring
    entire reviewed proteomes.
    """
    cache = FileCache(cache_dir)
    http = HttpClient(cache=cache, timeout=60)
    identifiers = _target_identifiers_from_tsv(tsv_path)
    summary = {
        "chunks": 0,
        "requested_identifiers": len(identifiers),
        "entries": 0,
        "written_entry_cache": 0,
    }
    for chunk in _chunks(identifiers, min(100, max(1, chunk_size))):
        query = " OR ".join(_uniprot_identifier_query(identifier) for identifier in chunk)
        url = f"{UNIPROT_SEARCH_BASE}?query={quote(query)}&format=json&size={len(chunk)}"
        payload = _fetch_page_with_retry(
            http=http,
            url=url,
            cache_namespace="uniprot_search_targets",
            cache_key=url,
        )
        summary["chunks"] += 1
        for entry in payload.get("results", []):
            summary["entries"] += 1
            summary["written_entry_cache"] += _cache_entry(cache, entry)
    return summary


def _prefetch_taxon(
    *,
    http: HttpClient,
    cache: FileCache,
    taxon_id: int,
    page_size: int,
    max_pages: int | None,
) -> dict[str, int]:
    pages = entries = written_entry_cache = 0
    query = f"(organism_id:{taxon_id}) AND (reviewed:true)"
    url = f"{UNIPROT_SEARCH_BASE}?query={quote(query)}&format=json&size={min(500, max(25, page_size))}"
    seen_urls: set[str] = set()
    seen_accessions: set[str] = set()
    while url and (max_pages is None or pages < max_pages):
        if url in seen_urls:
            raise ValueError("UniProt returned a repeated pagination URL")
        seen_urls.add(url)
        payload = _fetch_page_with_retry(
            http=http, url=url, cache_namespace="uniprot_search_bulk", cache_key=url, pagination=True,
        )
        results = payload.get("results", [])
        accessions = {entry["primaryAccession"] for entry in results}
        if accessions & seen_accessions:
            raise ValueError("UniProt repeated entries across pagination pages")
        seen_accessions.update(accessions)
        pages += 1
        for entry in results:
            entries += 1
            written_entry_cache += _cache_entry(cache, entry)
        url = payload.get("_next_url")
    return {
        "pages": pages,
        "entries": entries,
        "written_entry_cache": written_entry_cache,
    }


def _target_identifiers_from_tsv(tsv_path: str | Path) -> list[str]:
    path = Path(tsv_path)
    identifiers: list[str] = []
    seen: set[str] = set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            identifier = str(row.get("uniprot_id") or row.get("input_uniprot_id") or "").strip()
            if not identifier:
                identifier = str(row.get("uniprot_name") or row.get("input_uniprot_name") or "").strip()
            if not identifier or identifier in seen:
                continue
            seen.add(identifier)
            identifiers.append(identifier)
    return identifiers


def _uniprot_identifier_query(identifier: str) -> str:
    if "_" in identifier:
        return f"id:{identifier}"
    return f"accession:{identifier}"


def _chunks(values: list[str], size: int) -> list[list[str]]:
    return [values[index : index + size] for index in range(0, len(values), size)]


def _cache_entry(cache: FileCache, entry: dict[str, object]) -> int:
    written = 0
    identifiers: list[str] = []
    primary_accession = str(entry.get("primaryAccession") or "").strip()
    entry_name = str(entry.get("uniProtkbId") or "").strip()
    if primary_accession:
        identifiers.append(primary_accession)
    if entry_name:
        identifiers.append(entry_name)
    payload = json.dumps(entry)
    for identifier in identifiers:
        key_url = f"{UNIPROT_ENTRY_BASE}/{quote(identifier)}.json"
        path = cache.path_for("uniprot_entry", key_url, ".json")
        if not path.exists():
            path.write_text(payload, encoding="utf-8")
            written += 1
    return written


def _fetch_page_with_retry(
    *,
    http: HttpClient,
    url: str,
    cache_namespace: str,
    cache_key: str,
    attempts: int = 4,
    pagination: bool = False,
) -> dict[str, object]:
    last_error: Exception | None = None
    for _ in range(attempts):
        try:
            if pagination:
                payload, next_url = http.fetch_json_page(url)
                if isinstance(payload, dict):
                    payload = {**payload, "_next_url": next_url}
            else:
                payload = http.fetch_json(url, cache_namespace=cache_namespace, cache_key=cache_key)
            if isinstance(payload, dict):
                return payload
            raise TypeError(f"UniProt returned {type(payload).__name__}, expected a JSON object")
        except Exception as exc:
            last_error = exc
    raise RuntimeError(f"UniProt prefetch failed after {attempts} attempts: {url}") from last_error
