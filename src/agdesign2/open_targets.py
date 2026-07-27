from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
import csv
import json
import re
from pathlib import Path
from typing import Any, Protocol
from urllib.parse import urljoin, urlsplit
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

from .exceptions import ExternalServiceError


OPEN_TARGETS_GRAPHQL_URL = "https://api.platform.opentargets.org/api/v4/graphql"
OPEN_TARGETS_FTP_BASE_URL = "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform"
OPEN_TARGETS_DEFAULT_RELEASE = "latest"
OPEN_TARGETS_BULK_DATASETS = (
    "target",
    "disease",
    "association_overall_direct",
    "association_by_datasource_indirect",
)


def _require_http_url(url: str) -> None:
    parsed = urlsplit(url)
    if parsed.scheme not in {"http", "https"} or not parsed.netloc:
        raise ExternalServiceError(f"Unsupported Open Targets URL: {url}")


@dataclass(slots=True)
class OpenTargetsBuildResult:
    output_tsv: Path
    output_json: Path
    targets: list[dict[str, Any]]
    associations: list[dict[str, Any]]


class OpenTargetsApi(Protocol):
    def resolve_target(self, gene_symbol: str) -> dict[str, Any] | None:
        ...

    def fetch_associated_diseases(
        self,
        ensembl_id: str,
        *,
        page_size: int,
        max_pages: int,
    ) -> tuple[int | None, list[dict[str, Any]]]:
        ...


class OpenTargetsClient:
    def __init__(self, *, endpoint: str = OPEN_TARGETS_GRAPHQL_URL, timeout: int = 60) -> None:
        self.endpoint = endpoint
        self.timeout = timeout

    def graphql(self, query: str, variables: dict[str, Any] | None = None) -> dict[str, Any]:
        _require_http_url(self.endpoint)
        payload = json.dumps({"query": query, "variables": variables or {}}).encode("utf-8")
        request = Request(
            self.endpoint,
            data=payload,
            headers={
                "User-Agent": "openantigen/0.1",
                "Accept": "application/json",
                "Content-Type": "application/json",
            },
            method="POST",
        )
        try:
            with urlopen(request, timeout=self.timeout) as response:  # nosec B310
                data = json.loads(response.read().decode("utf-8"))
        except HTTPError as exc:
            raise ExternalServiceError(f"Open Targets HTTP error: {exc.code}") from exc
        except URLError as exc:
            raise ExternalServiceError(f"Open Targets connection error: {exc.reason}") from exc
        if data.get("errors"):
            message = "; ".join(str(item.get("message") or item) for item in data["errors"])
            raise ExternalServiceError(f"Open Targets GraphQL error: {message}")
        return data.get("data") or {}

    def resolve_target(self, gene_symbol: str) -> dict[str, Any] | None:
        cleaned = str(gene_symbol or "").strip()
        if not cleaned:
            return None
        query = """
        query SearchTarget($q: String!) {
          search(queryString: $q) {
            hits {
              object {
                __typename
                ... on Target {
                  id
                  approvedSymbol
                  approvedName
                }
              }
              score
            }
          }
        }
        """
        data = self.graphql(query, {"q": cleaned})
        hits = ((data.get("search") or {}).get("hits") or [])
        exact: list[dict[str, Any]] = []
        fallback: list[dict[str, Any]] = []
        for hit in hits:
            target = hit.get("object") or {}
            if target.get("__typename") != "Target":
                continue
            row = {
                "ensembl_id": target.get("id"),
                "approved_symbol": target.get("approvedSymbol"),
                "approved_name": target.get("approvedName"),
                "search_score": hit.get("score"),
            }
            if str(target.get("approvedSymbol") or "").upper() == cleaned.upper():
                exact.append(row)
            else:
                fallback.append(row)
        return (exact or fallback or [None])[0]

    def fetch_associated_diseases(
        self,
        ensembl_id: str,
        *,
        page_size: int,
        max_pages: int,
    ) -> tuple[int | None, list[dict[str, Any]]]:
        query = """
        query TargetDiseases($id: String!, $page: Pagination) {
          target(ensemblId: $id) {
            id
            approvedSymbol
            associatedDiseases(page: $page) {
              count
              rows {
                score
                disease {
                  id
                  name
                }
                datasourceScores {
                  id
                  score
                }
              }
            }
          }
        }
        """
        all_rows: list[dict[str, Any]] = []
        total_count: int | None = None
        for page_index in range(max(0, max_pages)):
            data = self.graphql(query, {"id": ensembl_id, "page": {"index": page_index, "size": page_size}})
            target = data.get("target") or {}
            associated = target.get("associatedDiseases") or {}
            if total_count is None:
                count = associated.get("count")
                total_count = int(count) if isinstance(count, int) else None
            rows = associated.get("rows") or []
            if not rows:
                break
            all_rows.extend(rows)
            # Stop once we have the reported total. When count is null/missing
            # (total_count is None) we must NOT treat it as 0 — that would break
            # after the first page and silently drop later associations. Instead
            # keep paging until a short/empty page or max_pages is reached.
            if total_count is not None and len(all_rows) >= total_count:
                break
            if len(rows) < page_size:
                break
        return total_count, all_rows


def build_open_targets_associations_from_summary(
    summary_path: str | Path,
    *,
    output_dir: str | Path | None = None,
    client: OpenTargetsApi | None = None,
    page_size: int = 100,
    max_pages: int = 5,
    limit: int | None = None,
    resume: bool = True,
    verbose: bool = False,
) -> OpenTargetsBuildResult:
    summary_file = Path(summary_path)
    batch_dir = summary_file.parent
    destination = Path(output_dir) if output_dir is not None else batch_dir
    destination.mkdir(parents=True, exist_ok=True)
    output_tsv = destination / "open_targets_disease_associations.tsv"
    output_json = destination / "open_targets_disease_associations.json"
    api = client or OpenTargetsClient()
    existing_targets, existing_associations = _load_existing_open_targets(output_json) if resume else ({}, [])
    results = json.loads(summary_file.read_text(encoding="utf-8"))
    target_rows: list[dict[str, Any]] = []
    association_rows: list[dict[str, Any]] = []
    processed = 0

    for index, item in enumerate(results, start=1):
        if limit is not None and processed >= limit:
            break
        report_data = _load_report_data(item, batch_dir)
        target = (report_data or {}).get("target") or {}
        entry_name = str(target.get("entry_name") or item.get("resolved_entry_name") or item.get("query") or "").strip()
        gene_symbol = str(target.get("gene_symbol") or item.get("gene_symbol") or item.get("query") or "").split("_", 1)[0].strip()
        accession = str(target.get("accession") or "").strip()
        cache_key = entry_name.upper() or f"ROW_{index}"
        if resume and cache_key in existing_targets:
            target_row = existing_targets[cache_key]
            target_rows.append(target_row)
            association_rows.extend([row for row in existing_associations if str(row.get("entry_name") or "").upper() == cache_key])
            continue

        target_row, rows_for_target = _build_one_open_targets_record(
            api=api,
            input_index=index,
            item=item,
            entry_name=entry_name,
            gene_symbol=gene_symbol,
            accession=accession,
            page_size=page_size,
            max_pages=max_pages,
        )
        target_rows.append(target_row)
        association_rows.extend(rows_for_target)
        processed += 1
        _write_open_targets_outputs(output_tsv, output_json, target_rows, association_rows)
        if verbose:
            status = target_row.get("status")
            top = target_row.get("top_disease_name") or "-"
            print(f"[open-targets] {index}/{len(results)} {entry_name or gene_symbol}: {status} {top}")

    _write_open_targets_outputs(output_tsv, output_json, target_rows, association_rows)
    return OpenTargetsBuildResult(
        output_tsv=output_tsv,
        output_json=output_json,
        targets=target_rows,
        associations=association_rows,
    )


def build_open_targets_associations_from_bulk_downloads(
    summary_path: str | Path,
    *,
    output_dir: str | Path | None = None,
    release: str = OPEN_TARGETS_DEFAULT_RELEASE,
    ftp_base_url: str = OPEN_TARGETS_FTP_BASE_URL,
    data_dir: str | Path = Path(".agdesign2/data/open_targets"),
    download: bool = True,
    verbose: bool = False,
) -> OpenTargetsBuildResult:
    """Build Open Targets associations from release Parquet downloads.

    This is the scalable path for portal builds. It downloads/uses the bulk
    Open Targets Platform datasets and performs local joins instead of calling
    GraphQL once per target.
    """

    summary_file = Path(summary_path)
    batch_dir = summary_file.parent
    destination = Path(output_dir) if output_dir is not None else batch_dir
    destination.mkdir(parents=True, exist_ok=True)
    output_tsv = destination / "open_targets_disease_associations.tsv"
    output_json = destination / "open_targets_disease_associations.json"
    requested_release = release
    release = resolve_open_targets_release(
        release=release,
        ftp_base_url=ftp_base_url,
        data_dir=data_dir,
        allow_network=download,
    )
    release_dir = Path(data_dir) / release
    if verbose:
        print(f"[open-targets] using Open Targets release {release}", flush=True)
    if download:
        download_open_targets_bulk_datasets(
            release=release,
            ftp_base_url=ftp_base_url,
            data_dir=release_dir,
            verbose=verbose,
        )

    target_rows_from_summary = _summary_target_rows(summary_file, batch_dir)
    target_lookup = _load_open_targets_bulk_targets(release_dir / "target")
    disease_lookup = _load_open_targets_bulk_diseases(release_dir / "disease")
    ensembl_ids = {
        target_lookup[str(row.get("gene_symbol") or "").upper()]["ensembl_id"]
        for row in target_rows_from_summary
        if str(row.get("gene_symbol") or "").upper() in target_lookup
    }
    direct_rows_raw = _load_open_targets_bulk_associations(
        release_dir / "association_overall_direct",
        target_ids=ensembl_ids,
    )
    indirect_rows_raw = _load_open_targets_bulk_associations(
        release_dir / "association_by_datasource_indirect",
        target_ids=ensembl_ids,
        include_datasource=True,
    )
    direct_by_pair: dict[tuple[str, str], float] = {}
    indirect_by_pair: dict[tuple[str, str], dict[str, Any]] = {}
    disease_ids_by_target: dict[str, set[str]] = {}
    for item in direct_rows_raw:
        target_id = str(item.get("targetId") or item.get("target_id") or "")
        disease_id = str(item.get("diseaseId") or item.get("disease_id") or "")
        if not target_id or not disease_id:
            continue
        direct_by_pair[(target_id, disease_id)] = _safe_float(item.get("score"))
        disease_ids_by_target.setdefault(target_id, set()).add(disease_id)
    for item in indirect_rows_raw:
        target_id = str(item.get("targetId") or item.get("target_id") or "")
        disease_id = str(item.get("diseaseId") or item.get("disease_id") or "")
        datasource_id = str(item.get("datasourceId") or item.get("datasource_id") or "")
        if not target_id or not disease_id:
            continue
        score = _safe_float(item.get("score"))
        pair = (target_id, disease_id)
        bucket = indirect_by_pair.setdefault(pair, {"score": 0.0, "datasource_scores": []})
        bucket["score"] = max(float(bucket["score"]), score)
        if datasource_id:
            bucket["datasource_scores"].append({"id": datasource_id, "score": score})
        disease_ids_by_target.setdefault(target_id, set()).add(disease_id)
    associations_by_target: dict[str, list[dict[str, Any]]] = {}
    for target_id, disease_ids in disease_ids_by_target.items():
        rows_for_target: list[dict[str, Any]] = []
        for disease_id in disease_ids:
            disease = disease_lookup.get(disease_id)
            if disease is None:
                continue
            indirect = indirect_by_pair.get((target_id, disease_id), {})
            indirect_datasource_scores = sorted(
                indirect.get("datasource_scores") or [],
                key=lambda row: _safe_float(row.get("score")),
                reverse=True,
            )
            direct_score = direct_by_pair.get((target_id, disease_id))
            indirect_score = indirect.get("score")
            rows_for_target.append(
                {
                    "score": indirect_score if indirect_score is not None else direct_score,
                    "direct_score": direct_score if direct_score is not None else "",
                    "indirect_score": indirect_score if indirect_score is not None else "",
                    "disease": {"id": disease_id, "name": disease.get("name") or disease_id},
                    "datasourceScores": [],
                    "indirectDatasourceScores": indirect_datasource_scores,
                }
            )
        associations_by_target[target_id] = rows_for_target

    target_rows: list[dict[str, Any]] = []
    association_rows: list[dict[str, Any]] = []
    for row in target_rows_from_summary:
        gene_symbol = str(row.get("gene_symbol") or "").upper()
        target_record = target_lookup.get(gene_symbol)
        if target_record is None:
            target_row = _open_targets_error_row(row, error=f"No Open Targets target found for {gene_symbol}")
            target_rows.append(target_row)
            continue
        ensembl_id = target_record["ensembl_id"]
        disease_rows = associations_by_target.get(ensembl_id, [])
        associations = _normalize_open_targets_associations(
            input_index=int(row.get("input_index") or 0),
            entry_name=str(row.get("entry_name") or ""),
            gene_symbol=str(row.get("gene_symbol") or ""),
            uniprot_accession=str(row.get("uniprot_accession") or ""),
            ensembl_id=ensembl_id,
            disease_rows=disease_rows,
        )
        _sort_and_rank_associations(associations)
        top = associations[0] if associations else {}
        disease_names = sorted({str(item.get("disease_name") or "") for item in associations if item.get("disease_name")}, key=str.lower)
        target_rows.append(
            {
                "input_index": row.get("input_index"),
                "query": row.get("query") or "",
                "entry_name": row.get("entry_name") or "",
                "gene_symbol": row.get("gene_symbol") or "",
                "uniprot_accession": row.get("uniprot_accession") or "",
                "ensembl_id": ensembl_id,
                "status": "ok",
                "error": "",
                "association_count_total": len(associations),
                "association_count_downloaded": len(associations),
                "top_disease_id": top.get("disease_id") or "",
                "top_disease_name": top.get("disease_name") or "",
                "top_disease_score": _blank_if_none(top.get("score")),
                "all_disease_names": "|".join(disease_names),
                "open_targets_url": f"https://platform.opentargets.org/target/{ensembl_id}/associations",
            }
        )
        association_rows.extend(associations)

    _write_open_targets_outputs(
        output_tsv,
        output_json,
        target_rows,
        association_rows,
        release=release,
        requested_release=requested_release,
    )
    return OpenTargetsBuildResult(
        output_tsv=output_tsv,
        output_json=output_json,
        targets=target_rows,
        associations=association_rows,
    )


def download_open_targets_bulk_datasets(
    *,
    release: str = OPEN_TARGETS_DEFAULT_RELEASE,
    ftp_base_url: str = OPEN_TARGETS_FTP_BASE_URL,
    data_dir: str | Path,
    jobs: int = 4,
    verbose: bool = False,
) -> None:
    root = Path(data_dir)
    for dataset in OPEN_TARGETS_BULK_DATASETS:
        dataset_url = f"{ftp_base_url.rstrip('/')}/{release}/output/{dataset}/"
        destination = root / dataset
        destination.mkdir(parents=True, exist_ok=True)
        urls = _list_open_targets_parquet_urls(dataset_url)
        if not urls:
            raise ExternalServiceError(f"No Open Targets Parquet files found at {dataset_url}")
        downloads: list[tuple[str, Path]] = []
        for url in urls:
            target_path = destination / Path(url).name
            if target_path.exists() and target_path.stat().st_size > 0:
                continue
            if verbose:
                print(f"[open-targets] downloading {url}", flush=True)
            downloads.append((url, target_path))
        if not downloads:
            continue
        with ThreadPoolExecutor(max_workers=min(max(1, jobs), len(downloads))) as executor:
            futures = {
                executor.submit(_download_file, url, target_path): url
                for url, target_path in downloads
            }
            for future in as_completed(futures):
                future.result()


def resolve_open_targets_release(
    *,
    release: str = OPEN_TARGETS_DEFAULT_RELEASE,
    ftp_base_url: str = OPEN_TARGETS_FTP_BASE_URL,
    data_dir: str | Path = Path(".agdesign2/data/open_targets"),
    allow_network: bool = True,
) -> str:
    requested = str(release or OPEN_TARGETS_DEFAULT_RELEASE).strip()
    if requested.lower() not in {"latest", "auto"}:
        return requested
    if allow_network:
        try:
            remote_releases = _list_open_targets_releases(ftp_base_url.rstrip("/") + "/")
        except ExternalServiceError:
            remote_releases = []
        if remote_releases:
            return remote_releases[-1]
    local_releases = _list_local_open_targets_releases(Path(data_dir))
    if local_releases:
        return local_releases[-1]
    raise ExternalServiceError(
        "Unable to resolve latest Open Targets release. Use --release VERSION or "
        "download once without --no-download."
    )


def _list_open_targets_releases(url: str) -> list[str]:
    _require_http_url(url)
    try:
        with urlopen(url, timeout=60) as response:  # nosec B310
            html = response.read().decode("utf-8", "replace")
    except URLError as exc:
        raise ExternalServiceError(f"Unable to list Open Targets releases at {url}: {exc.reason}") from exc
    releases: list[str] = []
    for href in re.findall(r'href=["\']([^"\']+)["\']', html):
        cleaned = href.strip("/").strip()
        if _open_targets_release_sort_key(cleaned) is not None:
            releases.append(cleaned)
    return sorted(set(releases), key=lambda item: _open_targets_release_sort_key(item) or ())


def _list_local_open_targets_releases(data_dir: Path) -> list[str]:
    if not data_dir.exists():
        return []
    releases = [
        path.name
        for path in data_dir.iterdir()
        if path.is_dir() and _open_targets_release_sort_key(path.name) is not None
    ]
    return sorted(set(releases), key=lambda item: _open_targets_release_sort_key(item) or ())


def _open_targets_release_sort_key(value: str) -> tuple[int, ...] | None:
    if not re.fullmatch(r"\d+(?:\.\d+)*", value):
        return None
    return tuple(int(part) for part in value.split("."))


def _list_open_targets_parquet_urls(url: str) -> list[str]:
    _require_http_url(url)
    try:
        with urlopen(url, timeout=60) as response:  # nosec B310
            html = response.read().decode("utf-8", "replace")
    except URLError as exc:
        raise ExternalServiceError(f"Unable to list Open Targets dataset {url}: {exc.reason}") from exc
    urls: list[str] = []
    for href in re.findall(r'href=["\']([^"\']+)["\']', html):
        if href.startswith("?") or href.startswith("../"):
            continue
        absolute = urljoin(url, href)
        if href.endswith(".parquet"):
            urls.append(absolute)
    return sorted(set(urls))


def _download_file(url: str, destination: Path) -> None:
    _require_http_url(url)
    temp_path = destination.with_suffix(destination.suffix + ".tmp")
    try:
        with urlopen(url, timeout=300) as response, temp_path.open("wb") as handle:  # nosec B310
            while True:
                chunk = response.read(1024 * 1024)
                if not chunk:
                    break
                handle.write(chunk)
    except URLError as exc:
        raise ExternalServiceError(f"Unable to download {url}: {exc.reason}") from exc
    temp_path.replace(destination)


def _summary_target_rows(summary_file: Path, batch_dir: Path) -> list[dict[str, Any]]:
    results = json.loads(summary_file.read_text(encoding="utf-8"))
    rows: list[dict[str, Any]] = []
    for index, item in enumerate(results, start=1):
        report_data = _load_report_data(item, batch_dir)
        target = (report_data or {}).get("target") or {}
        entry_name = str(target.get("entry_name") or item.get("resolved_entry_name") or item.get("query") or "").strip()
        gene_symbol = str(target.get("gene_symbol") or item.get("gene_symbol") or item.get("query") or "").split("_", 1)[0].strip()
        accession = str(target.get("accession") or "").strip()
        rows.append(
            {
                "input_index": index,
                "query": item.get("query") or "",
                "entry_name": entry_name,
                "gene_symbol": gene_symbol,
                "uniprot_accession": accession,
            }
        )
    return rows


def _load_open_targets_bulk_targets(path: Path) -> dict[str, dict[str, Any]]:
    table = _read_parquet_dataset(path, columns=["id", "approvedSymbol", "approvedName"])
    rows: dict[str, dict[str, Any]] = {}
    for item in table:
        symbol = str(item.get("approvedSymbol") or "").upper()
        ensembl_id = str(item.get("id") or "")
        if not symbol or not ensembl_id:
            continue
        rows[symbol] = {
            "ensembl_id": ensembl_id,
            "approved_symbol": item.get("approvedSymbol") or symbol,
            "approved_name": item.get("approvedName") or "",
        }
    return rows


def _load_open_targets_bulk_diseases(path: Path) -> dict[str, dict[str, Any]]:
    table = _read_open_targets_disease_dataset(path)
    rows: dict[str, dict[str, Any]] = {}
    for item in table:
        disease_id = str(item.get("id") or "")
        if disease_id and _is_open_targets_disease_like(item):
            rows[disease_id] = {
                "id": disease_id,
                "name": item.get("name") or disease_id,
                "ancestors": item.get("ancestors") or [],
                "therapeuticAreas": item.get("therapeuticAreas") or [],
            }
    return rows


def _read_open_targets_disease_dataset(path: Path) -> list[dict[str, Any]]:
    """Read Open Targets disease metadata with optional ontology context."""

    try:
        import pyarrow.dataset as ds
    except Exception as exc:
        raise ExternalServiceError(
            "Bulk Open Targets mode requires pyarrow. Install project dependencies with "
            "`.venv/bin/python -m pip install -e .` or install `pyarrow`."
        ) from exc
    dataset = ds.dataset(path, format="parquet")
    available = set(dataset.schema.names)
    columns = [column for column in ("id", "name", "ancestors", "therapeuticAreas") if column in available]
    missing = [column for column in ("id", "name") if column not in available]
    if missing:
        raise ExternalServiceError(f"Open Targets disease dataset {path} is missing expected columns: {', '.join(missing)}")
    return dataset.to_table(columns=columns).to_pylist()


_OPEN_TARGETS_DISEASE_ANCESTOR_IDS = {
    "EFO_0000408",  # disease
}


_OPEN_TARGETS_NON_DISEASE_IDS = {
    "EFO_0001444",  # measurement
    "EFO_0000651",  # phenotype
    "GO_0008150",  # biological_process
}


_OPEN_TARGETS_NON_DISEASE_ID_PREFIXES = (
    "GO_",
    "HP_",
    "MP_",
    "OBA_",
)


_OPEN_TARGETS_NON_DISEASE_NAME_EXACT = {
    "amount",
    "anatomical entity attribute",
    "anatomical entity morphology",
    "anatomical structure size",
    "attribute of cell",
    "behavior",
    "body height at birth",
    "body weight",
    "body weight gain",
    "body weights and measures",
    "biological process attribute",
    "biological_process",
    "birth weight",
    "cell size",
    "cell volume",
    "cellular component attribute",
    "cellular component size",
    "coffee consumption",
    "diastolic blood pressure",
    "drinking behavior",
    "duration",
    "hematocrit",
    "height-adjusted body mass index",
    "infant body height",
    "measurement",
    "multicellular organismal process",
    "overweight body mass index status",
    "pathologic process",
    "phenotype",
    "phenotypic abnormality",
    "platelet crit",
    "protein amount",
    "protein attribute",
    "sitting height ratio",
    "smoking behavior",
    "smoking cessation",
    "smoking initiation",
    "systolic blood pressure",
    "response to stimulus",
    "vital signs",
    "alcohol drinking",
    "age at initiation of smoking",
    "breastfeeding duration",
    "p wave duration",
    "sum of basophil and neutrophil counts",
    "sum of eosinophil and basophil counts",
    "sum of neutrophil and eosinophil counts",
    "weight-to-muscle ratio",
}


_OPEN_TARGETS_NON_DISEASE_NAME_SUFFIXES = (
    " amount",
    " attribute",
    " biomarker",
    " count",
    " density",
    " exposure measurement",
    " level",
    " measurement",
    " morphology trait",
    " quality",
    " size trait",
    " volume",
)


_OPEN_TARGETS_NON_DISEASE_NAME_PREFIXES = (
    "level of ",
)


def _as_string_set(value: Any) -> set[str]:
    if value is None:
        return set()
    if isinstance(value, str):
        return {value}
    if isinstance(value, dict):
        return {str(item) for item in value.values() if item is not None}
    try:
        return {str(item) for item in value if item is not None}
    except TypeError:
        return {str(value)}


def _is_open_targets_non_disease_name(name: str) -> bool:
    text = " ".join(str(name or "").strip().lower().split())
    if not text:
        return False
    if text in _OPEN_TARGETS_NON_DISEASE_NAME_EXACT:
        return True
    if any(text.startswith(prefix) for prefix in _OPEN_TARGETS_NON_DISEASE_NAME_PREFIXES):
        return True
    return any(text.endswith(suffix) for suffix in _OPEN_TARGETS_NON_DISEASE_NAME_SUFFIXES)


def _is_open_targets_disease_like(item: dict[str, Any]) -> bool:
    disease_id = str(item.get("id") or item.get("disease_id") or "").strip()
    name = str(item.get("name") or item.get("disease_name") or "").strip()
    if (
        disease_id in _OPEN_TARGETS_NON_DISEASE_IDS
        or disease_id.startswith(_OPEN_TARGETS_NON_DISEASE_ID_PREFIXES)
        or _is_open_targets_non_disease_name(name)
    ):
        return False
    if disease_id.startswith(("MONDO_", "Orphanet_", "ORPHA:")):
        return True
    ancestors = _as_string_set(item.get("ancestors"))
    if ancestors & _OPEN_TARGETS_DISEASE_ANCESTOR_IDS:
        return True
    # GraphQL/resumed rows often only carry ID/name. In that case, apply only
    # explicit exclusions and keep terms unless they are clearly non-disease.
    if not ancestors:
        return not _is_open_targets_non_disease_name(name)
    return False


def _load_open_targets_bulk_associations(
    path: Path,
    *,
    target_ids: set[str],
    include_datasource: bool = False,
) -> list[dict[str, Any]]:
    if not target_ids:
        return []
    columns = ["targetId", "diseaseId", "score"]
    column_aliases = {"score": ["score", "associationScore"]}
    if include_datasource:
        columns.append("datasourceId")
        column_aliases["datasourceId"] = ["datasourceId", "aggregationValue"]
    return _read_parquet_dataset(
        path,
        columns=columns,
        filter_target_ids=target_ids,
        column_aliases=column_aliases,
    )


def _read_parquet_dataset(
    path: Path,
    *,
    columns: list[str],
    filter_target_ids: set[str] | None = None,
    column_aliases: dict[str, list[str]] | None = None,
) -> list[dict[str, Any]]:
    try:
        import pyarrow.compute as pc
        import pyarrow.dataset as ds
    except Exception as exc:
        raise ExternalServiceError(
            "Bulk Open Targets mode requires pyarrow. Install project dependencies with "
            "`.venv/bin/python -m pip install -e .` or install `pyarrow`."
        ) from exc
    dataset = ds.dataset(path, format="parquet")
    schema_names = set(dataset.schema.names)
    aliases = column_aliases or {}
    column_map: dict[str, str] = {}
    missing: list[str] = []
    for column in columns:
        candidates = aliases.get(column, [column])
        source_column = next((candidate for candidate in candidates if candidate in schema_names), None)
        if source_column is None:
            missing.append(column)
        else:
            column_map[column] = source_column
    if missing:
        raise ExternalServiceError(f"Open Targets dataset {path} is missing expected columns: {', '.join(missing)}")
    filter_expr = None
    if filter_target_ids is not None:
        if "targetId" not in schema_names:
            raise ExternalServiceError(f"Open Targets association dataset {path} is missing targetId.")
        filter_expr = pc.field("targetId").isin(sorted(filter_target_ids))
    source_columns = list(dict.fromkeys(column_map.values()))
    table = dataset.to_table(columns=source_columns, filter=filter_expr)
    rows = table.to_pylist()
    if all(column_map[column] == column for column in columns):
        return rows
    normalized_rows: list[dict[str, Any]] = []
    for row in rows:
        normalized = dict(row)
        for output_column, source_column in column_map.items():
            if output_column != source_column:
                normalized[output_column] = row.get(source_column)
        normalized_rows.append(normalized)
    return normalized_rows


def _open_targets_error_row(row: dict[str, Any], *, error: str) -> dict[str, Any]:
    return {
        "input_index": row.get("input_index"),
        "query": row.get("query") or "",
        "entry_name": row.get("entry_name") or "",
        "gene_symbol": row.get("gene_symbol") or "",
        "uniprot_accession": row.get("uniprot_accession") or "",
        "ensembl_id": "",
        "status": "error",
        "error": error,
        "association_count_total": 0,
        "association_count_downloaded": 0,
        "top_disease_id": "",
        "top_disease_name": "",
        "top_disease_score": "",
        "all_disease_names": "",
        "open_targets_url": "",
    }


def _build_one_open_targets_record(
    *,
    api: OpenTargetsApi,
    input_index: int,
    item: dict[str, Any],
    entry_name: str,
    gene_symbol: str,
    accession: str,
    page_size: int,
    max_pages: int,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    base = {
        "input_index": input_index,
        "query": item.get("query") or "",
        "entry_name": entry_name,
        "gene_symbol": gene_symbol,
        "uniprot_accession": accession,
        "ensembl_id": "",
        "status": "error",
        "error": "",
        "association_count_total": 0,
        "association_count_downloaded": 0,
        "top_disease_id": "",
        "top_disease_name": "",
        "top_disease_score": "",
        "all_disease_names": "",
        "open_targets_url": "",
    }
    if not gene_symbol:
        base["error"] = "Missing gene symbol"
        return base, []
    try:
        resolved = api.resolve_target(gene_symbol)
        if not resolved or not resolved.get("ensembl_id"):
            base["error"] = f"No Open Targets target found for {gene_symbol}"
            return base, []
        ensembl_id = str(resolved["ensembl_id"])
        total_count, disease_rows = api.fetch_associated_diseases(
            ensembl_id,
            page_size=page_size,
            max_pages=max_pages,
        )
    except Exception as exc:
        base["error"] = str(exc)
        return base, []

    associations = _normalize_open_targets_associations(
        input_index=input_index,
        entry_name=entry_name,
        gene_symbol=gene_symbol,
        uniprot_accession=accession,
        ensembl_id=ensembl_id,
        disease_rows=disease_rows,
    )
    associations = [
        row
        for row in associations
        if _is_open_targets_disease_like({"id": row.get("disease_id"), "name": row.get("disease_name")})
    ]
    _sort_and_rank_associations(associations)
    top = associations[0] if associations else {}
    disease_names = sorted({str(row.get("disease_name") or "") for row in associations if row.get("disease_name")}, key=str.lower)
    base.update(
        {
            "ensembl_id": ensembl_id,
            "status": "ok",
            "error": "",
            "association_count_total": total_count if total_count is not None else len(associations),
            "association_count_downloaded": len(associations),
            "top_disease_id": top.get("disease_id") or "",
            "top_disease_name": top.get("disease_name") or "",
            "top_disease_score": _blank_if_none(top.get("score")),
            "all_disease_names": "|".join(disease_names),
            "open_targets_url": f"https://platform.opentargets.org/target/{ensembl_id}/associations",
        }
    )
    return base, associations


def _normalize_open_targets_associations(
    *,
    input_index: int,
    entry_name: str,
    gene_symbol: str,
    uniprot_accession: str,
    ensembl_id: str,
    disease_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for rank, item in enumerate(disease_rows, start=1):
        disease = item.get("disease") or {}
        datasource_scores = item.get("datasourceScores") or []
        indirect_datasource_scores = item.get("indirectDatasourceScores") or []
        direct_score = item.get("direct_score")
        indirect_score = item.get("indirect_score")
        score = item.get("score") if item.get("score") is not None else ""
        if direct_score is None:
            direct_score = score
        if indirect_score is None:
            indirect_score = ""
        rows.append(
            {
                "input_index": input_index,
                "entry_name": entry_name,
                "gene_symbol": gene_symbol,
                "uniprot_accession": uniprot_accession,
                "ensembl_id": ensembl_id,
                "rank": rank,
                "disease_id": disease.get("id") or "",
                "disease_name": disease.get("name") or "",
                "score": score,
                "direct_score": direct_score,
                "indirect_score": indirect_score,
                "datasource_scores": ";".join(
                    f"{score.get('id')}:{score.get('score')}"
                    for score in datasource_scores
                    if score.get("id") is not None and score.get("score") is not None
                ),
                "indirect_datasource_scores": ";".join(
                    f"{score.get('id')}:{score.get('score')}"
                    for score in indirect_datasource_scores
                    if score.get("id") is not None and score.get("score") is not None
                ),
                "open_targets_url": f"https://platform.opentargets.org/evidence/{ensembl_id}/{disease.get('id')}",
            }
        )
    return rows


def _load_report_data(item: dict[str, Any], batch_dir: Path) -> dict[str, Any] | None:
    raw_path = item.get("json_report")
    if not raw_path:
        return None
    path = Path(str(raw_path))
    if not path.is_absolute():
        path = batch_dir / path
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None


def _load_existing_open_targets(output_json: Path) -> tuple[dict[str, dict[str, Any]], list[dict[str, Any]]]:
    if not output_json.exists():
        return {}, []
    try:
        payload = json.loads(output_json.read_text(encoding="utf-8"))
    except Exception:
        return {}, []
    targets = {
        str(item.get("entry_name") or "").upper(): item
        for item in payload.get("targets", [])
        if item.get("entry_name")
    }
    associations = [dict(row) for row in payload.get("associations", [])]
    return targets, associations


def _write_open_targets_outputs(
    output_tsv: Path,
    output_json: Path,
    targets: list[dict[str, Any]],
    associations: list[dict[str, Any]],
    *,
    release: str | None = None,
    requested_release: str | None = None,
) -> None:
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "input_index",
        "entry_name",
        "gene_symbol",
        "uniprot_accession",
        "ensembl_id",
        "rank",
        "disease_id",
        "disease_name",
        "score",
        "direct_score",
        "indirect_score",
        "datasource_scores",
        "indirect_datasource_scores",
        "open_targets_url",
    ]
    tsv_tmp = output_tsv.with_suffix(output_tsv.suffix + ".tmp")
    json_tmp = output_json.with_suffix(output_json.suffix + ".tmp")
    with tsv_tmp.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in associations:
            writer.writerow({field: row.get(field, "") for field in fields})
    payload = {
        "generated_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "open_targets_release": release or "",
        "open_targets_release_requested": requested_release or "",
        "target_count": len(targets),
        "association_count": len(associations),
        "targets": targets,
        "associations": associations,
    }
    json_tmp.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    tsv_tmp.replace(output_tsv)
    json_tmp.replace(output_json)


def _safe_float(value: Any) -> float:
    try:
        return float(value)
    except Exception:
        return 0.0


def _blank_if_none(value: Any) -> Any:
    return "" if value is None else value


def _sort_and_rank_associations(rows: list[dict[str, Any]]) -> None:
    rows.sort(key=lambda row: _safe_float(row.get("score")), reverse=True)
    for rank, row in enumerate(rows, start=1):
        row["rank"] = rank
