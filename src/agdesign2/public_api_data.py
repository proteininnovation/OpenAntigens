"""Compile immutable public data from frozen reports and finalized portal indexes.

This module reads no ambient caches and performs no network or scientific analysis.
Only explicitly projected fields enter the service bundle.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import math
import os
import re
import shutil
import struct
import tempfile
from html.parser import HTMLParser
from collections import Counter
from pathlib import Path, PurePosixPath
from typing import Any, Iterable
from urllib.parse import quote, urlsplit
from jsonschema import Draft202012Validator

from .public_api_schema import (
    CATEGORIES, DEFS, PUBLIC_SCHEMA, SCHEMA_VERSION, SECTIONS, SECTION_TYPES,
    SPECIES, STATUSES, PublicDataError, validate_schema,
)

INDEX_LIMIT = 32 * 1024 * 1024
INDEX_TOTAL_LIMIT = 64 * 1024 * 1024
MANIFEST_LIMIT = 64 * 1024 * 1024
TARGET_LIMIT = 8 * 1024 * 1024
SHARD_LIMIT = 1024 * 1024
PAGE_SIZE = 25
LICENSE = "CC BY 4.0 for OpenAntigens annotations; source data retain original terms."
LICENSES = {"software": "Apache-2.0", "annotations": "CC BY 4.0", "third_party": "Third-party sequences, structures and annotations retain original source terms; see DATA_LICENSE.md."}
SOURCE_DOMAINS = frozenset({"www.uniprot.org", "rest.uniprot.org", "www.rcsb.org", "www.ebi.ac.uk", "alphafold.ebi.ac.uk", "www.ncbi.nlm.nih.gov", "platform.opentargets.org", "identifiers.org", "doi.org", "orcid.org", "www.genenames.org", "www.ensembl.org", "gpcrdb.org"})
_PRIVATE_MARKERS = ("/Users/", "/private/", "/home/", "file://", "-----BEGIN PRIVATE KEY", "Authorization: Bearer", "PRIVATE_CANARY", "sk-proj-")


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")


def _hash(value: Any) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def semantic_hash(value: Any) -> str:
    """Portable payload digest: JSON numbers use their IEEE-754 bits, not float spelling."""
    checksum = hashlib.sha256()
    _semantic_update(checksum, value)
    return checksum.hexdigest()


def _semantic_update(checksum: Any, item: Any) -> None:
    if isinstance(item, (int, float)) and not isinstance(item, bool):
        if not math.isfinite(item) or (isinstance(item, int) and abs(item) > 2 ** 53 - 1):
            raise PublicDataError("Public numeric value is outside the portable contract")
        checksum.update(canonical_bytes({"$number": struct.pack(">d", float(item)).hex()}))
    elif isinstance(item, list):
        checksum.update(b"[")
        for position, child in enumerate(item):
            if position:
                checksum.update(b",")
            _semantic_update(checksum, child)
        checksum.update(b"]")
    elif isinstance(item, dict):
        checksum.update(b"{")
        for position, key in enumerate(sorted(item)):
            if position:
                checksum.update(b",")
            checksum.update(canonical_bytes(key) + b":")
            _semantic_update(checksum, item[key])
        checksum.update(b"}")
    else:
        checksum.update(canonical_bytes(item))


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _sha_file(path: Path) -> str:
    checksum = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            checksum.update(block)
    return checksum.hexdigest()


def _read_json(path: Path) -> Any:
    return _read_json_hashed(path)[0]


def _read_json_hashed(path: Path, *, maximum: int | None = None) -> tuple[Any, str]:
    try:
        if maximum is not None:
            with path.open("rb") as handle:
                raw = handle.read(maximum + 1)
            if len(raw) > maximum:
                raise PublicDataError("Declared frozen JSON exceeds its bounded input budget")
        else:
            raw = path.read_bytes()
        content = raw
        if content.startswith(b"\x1f\x8b"):
            if maximum is None:
                content = gzip.decompress(content)
            else:
                with gzip.GzipFile(fileobj=io.BytesIO(raw)) as handle:
                    content = handle.read(maximum + 1)
                if len(content) > maximum:
                    raise PublicDataError("Declared frozen JSON exceeds its decoded input budget")
        return json.loads(content), _sha(raw)
    except (OSError, ValueError) as exc:
        raise PublicDataError("Cannot read declared public-build input") from exc


def _safe_text(value: Any) -> str | None:
    if value is None or value == "":
        return None
    if not isinstance(value, str):
        raise PublicDataError("Expected recorded text")
    if len(value) > 16384 or any(marker in value for marker in _PRIVATE_MARKERS):
        raise PublicDataError("Unsafe or oversized text in public projection")
    if any(ord(char) < 32 and char not in "\n\t\r" for char in value):
        raise PublicDataError("Control character in public text")
    return value


def _number(value: Any, *, integer: bool = False) -> int | float | None:
    if value is None or value == "":
        return None
    if isinstance(value, bool):
        raise PublicDataError("Boolean is not a recorded numeric value")
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise PublicDataError("Invalid recorded numeric value") from exc
    if not math.isfinite(parsed) or (integer and not parsed.is_integer()):
        raise PublicDataError("Invalid recorded numeric value")
    return int(parsed) if integer else parsed


def _project(schema_name: str, record: Any) -> dict[str, Any]:
    """Copy only named, recursively typed fields; arbitrary metadata never passes."""
    if not isinstance(record, dict):
        raise PublicDataError("Expected recorded object")

    def project(schema: dict[str, Any], value: Any) -> Any:
        if "$ref" in schema:
            return project(DEFS[schema["$ref"].rsplit("/", 1)[-1]], value)
        if "anyOf" in schema:
            return None if value is None or value == "" else project(schema["anyOf"][0], value)
        kind = schema.get("type")
        if kind == "object":
            if not isinstance(value, dict):
                raise PublicDataError("Expected recorded nested object")
            return {key: project(child, value.get(key)) for key, child in schema["properties"].items()}
        if kind == "array":
            if value is None:
                return []
            if not isinstance(value, list):
                raise PublicDataError("Expected recorded list")
            return [project(schema["items"], child) for child in value]
        if kind in {"number", "integer"}:
            return _number(value, integer=kind == "integer")
        if kind == "string":
            return _safe_text(value)
        if kind == "boolean" and value is not None and not isinstance(value, bool):
            raise PublicDataError("Expected recorded Boolean")
        return value

    result = project(DEFS[schema_name], record)
    validate_schema(schema_name, result)
    return result


def _url(value: Any, *, static_origin: str) -> str | None:
    text = _safe_text(value)
    if text is None:
        return None
    parsed = urlsplit(text)
    allowed = SOURCE_DOMAINS | {urlsplit(static_origin).hostname}
    if parsed.scheme != "https" or parsed.hostname not in allowed or parsed.username or parsed.password or parsed.port not in (None, 443):
        return None
    return text


def _scan_public(value: Any) -> None:
    if isinstance(value, str):
        if len(value) > 1000000 or any(marker in value for marker in _PRIVATE_MARKERS):
            raise PublicDataError("Unsafe or oversized value in public artifact")
        if any(ord(char) < 32 and char not in "\n\t\r" for char in value):
            raise PublicDataError("Control character in public artifact")
    elif isinstance(value, list):
        for child in value:
            _scan_public(child)
    elif isinstance(value, dict):
        for key, child in value.items():
            _safe_text(key)
            _scan_public(child)


def _safe_path(root: Path, relative: str, *, exists: bool = True) -> Path:
    if not isinstance(relative, str) or not relative or "\\" in relative or "\x00" in relative:
        raise PublicDataError("Invalid public artifact path")
    parts = PurePosixPath(relative)
    if parts.is_absolute() or any(part in {".", ".."} for part in relative.split("/")) or "//" in relative:
        raise PublicDataError("Invalid public artifact path")
    path = root.joinpath(*parts.parts)
    for parent in (path, *path.parents):
        if parent == root.parent:
            break
        if parent.is_symlink():
            raise PublicDataError("Public artifacts cannot use symbolic links")
    if not path.resolve().is_relative_to(root.resolve()) or (exists and not path.is_file()):
        raise PublicDataError("Missing or out-of-root public artifact")
    return path


def _verified_final_model(report: dict[str, Any], model: Path, source: Any) -> str:
    from .portal import _alphafold_pdb_matches_target_sequence, _report_uses_af3_structure, _structure_source_label
    if source != "AlphaFold" or _report_uses_af3_structure(report):
        raise PublicDataError("Final predicted model has unsupported source provenance")
    metadata = model.with_suffix(".meta.json")
    if metadata.exists():
        _safe_path(model.parent, metadata.name)
        recorded = _read_json(metadata)
        if not isinstance(recorded, dict):
            raise PublicDataError("Invalid final model provenance")
        if recorded.get("source") not in (None, "alphafold", "alphafold_db", "alphafold_server") or recorded.get("model_type") not in (None, "AF2", "AF3"):
            raise PublicDataError("Final predicted model has unsupported source provenance")
        target = report.get("target") or {}
        if recorded.get("canonical_sequence") not in (None, target.get("sequence")) or recorded.get("artifact_validation_version") not in (None, 1):
            raise PublicDataError("Final model provenance disagrees with its recorded reference")
        if recorded.get("model_uniprot_accession") and recorded["model_uniprot_accession"] not in {target.get("accession"), target.get("canonical_isoform_id")}:
            raise PublicDataError("Final model accession provenance disagrees with its target")
    if _structure_source_label(model) != "AlphaFold":
        raise PublicDataError("Final predicted model has excluded source provenance")
    maximum = max(1024 * 1024, len(str((report.get("target") or {}).get("sequence") or "")) * 10000)
    with model.open("rb") as handle:
        compressed = handle.read(2) == b"\x1f\x8b"
    opener = gzip.open if compressed else open
    with opener(model, "rb") as handle:
        content = handle.read(maximum + 1)
    if len(content) > maximum:
        raise PublicDataError("Final predicted model exceeds its validation bound")
    # A fresh file avoids trusting a match cached before final packaging replaced the asset.
    with tempfile.TemporaryDirectory(prefix="openantigens-model-check-") as temporary:
        plain = Path(temporary) / "final.pdb"
        plain.write_bytes(content)
        if not _alphafold_pdb_matches_target_sequence(report, plain):
            raise PublicDataError("Final predicted model sequence disagrees with exported reference")
    return _sha_file(model)


def _frozen_pubtator_links(portal_root: Path, index: list[dict[str, Any]], *, selected_entries_only: bool = False) -> tuple[dict[str, str], str | None]:
    path = portal_root / "portal-index-data.js"
    if not path.exists():
        return {}, None
    path = _safe_path(portal_root, "portal-index-data.js")
    with path.open("rb") as handle:
        compressed = handle.read(2) == b"\x1f\x8b"
    opener = gzip.open if compressed else open
    with opener(path, "rb") as handle:
        content = handle.read(INDEX_LIMIT + 1)
    if len(content) > INDEX_LIMIT:
        raise PublicDataError("Frozen portal index script exceeds public index budget")
    text = content.decode("utf-8").strip()
    assignment = "window.OpenAntigenIndexRows = "
    if not text.startswith(assignment) or not text.endswith(";"):
        raise PublicDataError("Unrecognized frozen portal index script format")
    payload = json.loads(text[len(assignment):-1])
    expected = {str(row.get("entry_name") or "").lower(): row for row in index}
    result, seen = {}, set()

    class Links(HTMLParser):
        def __init__(self) -> None:
            super().__init__(convert_charrefs=True)
            self.links: set[str] = set()

        def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
            href = dict(attrs).get("href")
            if tag == "a" and href:
                parsed = urlsplit(href)
                if parsed.scheme == "https" and parsed.hostname == "www.ncbi.nlm.nih.gov" and parsed.path == "/research/pubtator3/docsum" and not parsed.username and not parsed.password:
                    self.links.add(href)

    if not isinstance(payload, list):
        raise PublicDataError("Invalid frozen public index rows")
    for row in payload:
        if not isinstance(row, dict) or not isinstance(row.get("dataset"), dict) or not isinstance(row.get("cells"), list):
            raise PublicDataError("Invalid frozen public index row")
        key = row["dataset"].get("entry")
        if selected_entries_only and key not in expected:
            continue
        if key not in expected or key in seen:
            raise PublicDataError("Frozen public index script and compact index identities disagree")
        seen.add(key)
        count = _number(expected[key].get("pubtator_count"), integer=True)
        if count is not None and count != _number(row["dataset"].get("pubtator"), integer=True):
            raise PublicDataError("Frozen PubTator count disagrees with its published query row")
        parser = Links()
        for cell in row["cells"]:
            if not isinstance(cell, str):
                raise PublicDataError("Invalid frozen public index cell")
            parser.feed(cell)
        if len(parser.links) > 1:
            raise PublicDataError("Ambiguous frozen PubTator query")
        if parser.links:
            result[key] = next(iter(parser.links))
    if seen != set(expected):
        raise PublicDataError("Frozen public index script is incomplete")
    return result, _sha_file(path)


def _identifier_source(accession: str | None) -> str | None:
    if accession is None:
        return None
    if re.fullmatch(r"(?:NP|XP|YP|WP|AP)_\d+(?:\.\d+)?", accession):
        return "refseq"
    if re.fullmatch(r"(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-\d+)?", accession):
        return "uniprot"
    raise PublicDataError("Unrecognized recorded target accession namespace")


def _sequence_provenance(report: dict[str, Any], entry: dict[str, Any], species: str) -> tuple[str | None, str | None]:
    target = report.get("target") or {}
    metadata = report.get("metadata") or {}
    provenance = metadata.get("mouse_sequence_provenance") or {} if species == "mouse" else {}
    source = _safe_text(provenance.get("sequence_source") or target.get("sequence_source") or (entry.get("mouse_sequence_source") if species == "mouse" else None))
    accession = _safe_text(target.get("sequence_accession"))
    if source and ("refseq" in source.lower() or "ortholog_reference_table" in source.lower()):
        accession = accession or _safe_text(provenance.get("refseq_accession") or target.get("refseq_accession") or (entry.get("mouse_refseq_accession") if species == "mouse" else None))
    elif source and "uniprot" in source.lower():
        accession = accession or _safe_text(target.get("accession"))
    elif species == "human" and target.get("accession"):
        accession = accession or _safe_text(target.get("accession"))
        source = _identifier_source(accession)
    elif _identifier_source(_safe_text(target.get("accession"))) == "refseq":
        accession = accession or _safe_text(target.get("accession"))
        source = source or "refseq"
    return accession, source


def _state(items: list[Any], *, missing: bool = False, complete: bool = False) -> tuple[str, str | None]:
    if items:
        return "available", None
    if complete:
        return "no_results", "The recorded query completed with no released results."
    return "unknown", "The frozen source does not distinguish absent results from an unrecorded query."


def _track(value: Any) -> str:
    return {"secreted": "secreted", "gpi": "gpi", "single-pass": "single_pass", "single_pass": "single_pass", "multipass": "multipass", "other": "other"}.get(str(value or "").strip().lower(), "unknown")


def _sequence_record(summary: dict[str, Any], sequence: str, *, construct_id: str | None = None, start: int = 1, end: int | None = None, basis: str = "reference", warnings: list[str] | None = None, evidence: list[str] | None = None, accession: str | None = None, source: str | None = None) -> dict[str, Any]:
    record = {key: summary[key] for key in ("target_id", "species", "taxon_id", "identifier_accession", "identifier_source")}
    record.update({"construct_id": construct_id, "sequence_accession": accession if accession is not None else summary["sequence_accession"], "sequence_source": source if source is not None else summary["sequence_source"], "sequence_basis": basis, "start": start, "end": len(sequence) if end is None else end, "length": len(sequence), "sequence": sequence, "sequence_sha256": _sha(sequence.encode("ascii")), "warnings": warnings or [], "evidence": evidence or []})
    validate_schema("sequence", record)
    return record


def _hazards(report: dict[str, Any], construct: dict[str, Any]) -> dict[str, Any]:
    start, end = _number(construct.get("start"), integer=True), _number(construct.get("end"), integer=True)

    def selected(records: Any, *, position: bool = False) -> list[Any]:
        if not isinstance(records, list):
            raise PublicDataError("Invalid recorded hazard list")
        if start is None or end is None:
            return []
        result = []
        for record in records:
            low = _number(record.get("position" if position else "start"), integer=True)
            high = low if position else _number(record.get("end"), integer=True)
            if low is not None and high is not None and low <= end and high >= start:
                result.append(record)
        return result

    cysteines = construct.get("cysteine_analysis") if "cysteine_analysis" in construct else selected(report.get("cysteine_analysis") or [], position=True)
    analysis_start, analysis_end = (start, end) if "cysteine_analysis" in construct else _report_cysteine_interval(report)
    cysteines = [{**record, "analysis_start": analysis_start, "analysis_end": analysis_end} for record in cysteines or []]
    return _project("hazards", {
        "cysteines": cysteines,
        "furin_sites": selected(report.get("furin_sites") or []),
        "ptms": construct.get("ptms") if "ptms" in construct else selected(report.get("ptms") or []),
        "ligand_interactions": construct.get("ligand_interactions") if "ligand_interactions" in construct else selected(report.get("ligand_interactions") or []),
        "complete_domains": construct.get("complete_domain_annotations") or [],
        "clipped_domains": construct.get("clipped_domain_annotations") or [],
    })


def _report_cysteine_interval(report: dict[str, Any]) -> tuple[int | None, int | None]:
    region = report.get("ectodomain") or {}
    if region:
        return _number(region.get("start"), integer=True), _number(region.get("end"), integer=True)
    if str((report.get("topology") or {}).get("topology_class") or "").startswith("multipass"):
        return 1, len(str((report.get("target") or {}).get("sequence") or ""))
    return None, None


def _hazard_warnings(hazards: dict[str, Any], *, start: int | None, end: int | None) -> list[str]:
    warnings = []
    for row in hazards["cysteines"]:
        interval = f"reference interval {row['analysis_start']}-{row['analysis_end']}" if row["analysis_start"] is not None and row["analysis_end"] is not None else "an unrecorded source interval"
        if row["warning"]:
            warnings.append(f"Recorded cysteine analysis for {interval}: {row['warning']}")
        if row["partner_outside_range"] is not None:
            partner = row["partner_outside_range"]
            if start is not None and end is not None and not start <= partner <= end:
                warnings.append(f"Recorded cysteine {row['position']} has annotated partner {partner} outside exported interval {start}-{end}.")
            else:
                warnings.append(f"Recorded cysteine {row['position']} has annotated partner {partner}; the outside-range flag belongs to the source analysis for {interval}.")
    for row in hazards["furin_sites"]:
        warnings.append(f"Recorded furin motif {row['motif'] or 'unknown'} at {row['start']}-{row['end']}." + (f" {row['curated_annotation']}" if row["curated_annotation"] else ""))
    for kind in ("ptms", "ligand_interactions"):
        for row in hazards[kind]:
            warnings.append(f"Recorded {row['type'] or kind} at {row['start']}-{row['end']}." + (f" {row['description']}" if row["description"] else ""))
    for row in hazards["clipped_domains"]:
        warnings.append(f"Recorded clipped domain {row['accession'] or row['name'] or 'unknown'} at {row['start']}-{row['end']}.")
    return warnings


def _constructs(report: dict[str, Any], summary: dict[str, Any]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    from .portal import _construct_category, _construct_size_warning

    reference = str((report.get("target") or {}).get("sequence") or "")
    suggestions = report.get("construct_details") or report.get("construct_recommendations") or []
    if not isinstance(suggestions, list):
        raise PublicDataError("Invalid recorded construct list")
    records: dict[str, dict[str, Any]] = {}
    sequences: dict[str, Any] = {}
    for construct in suggestions:
        if not isinstance(construct, dict):
            raise PublicDataError("Invalid recorded construct")
        start, end = _number(construct.get("start"), integer=True), _number(construct.get("end"), integer=True)
        sequence = str(construct.get("sequence") or "")
        length = _number(construct.get("length"), integer=True)
        category = _construct_category(construct)
        hazards = _hazards(report, construct)
        warnings = [_safe_text(x) for x in (construct.get("warnings") or [])]
        warnings.extend(_hazard_warnings(hazards, start=start, end=end))
        size_warning = _construct_size_warning(construct)
        if size_warning:
            warnings.insert(0, size_warning)
        evidence = [_safe_text(x) for x in (construct.get("evidence") or [])]
        state, reason, basis = "available", None, "reference_slice"
        if not sequence:
            state, reason, basis = "unavailable", "No exact sequence was recorded for this construct.", "unavailable"
        else:
            if not re.fullmatch(r"[A-Z*]+", sequence):
                raise PublicDataError("Invalid recorded protein sequence")
            if length is not None and length != len(sequence):
                raise PublicDataError("Construct sequence length disagrees with declared length")
            if not isinstance(start, int) or not isinstance(end, int) or not 1 <= start <= end <= len(reference):
                raise PublicDataError("Construct coordinates do not fit its reference")
            if sequence != reference[start - 1:end]:
                # The report's PDB windows normally contain native sequence, not deposited engineering.
                state, reason, basis = "unavailable", "Recorded sequence differs from its reference interval; explicit variant mapping is required.", "unavailable"
        record = {"construct_id": "", "target_id": summary["target_id"], "name": _safe_text(construct.get("name")), "category": category, "classification": _safe_text(construct.get("classification")), "start": start, "end": end, "length": len(sequence) if sequence else length, "sequence_basis": basis, "sequence_sha256": _sha(sequence.encode("ascii")) if sequence else None, "sequence_accession": summary["sequence_accession"], "sequence_source": summary["sequence_source"], "export_state": state, "export_reason": reason, "score": _number(construct.get("score")), "rationale": _safe_text(construct.get("rationale")), "warnings": warnings, "evidence": evidence, "pdb_id": _safe_text(construct.get("pdb_id")), "structural_metrics": _project("metrics", construct["structural_metrics"]) if construct.get("structural_metrics") else None}
        record["hazards"] = hazards
        identity = [summary["target_id"], category, record["classification"], {"name": record["name"], "pdb_id": record["pdb_id"]}, {"start": start, "end": end, "sequence_basis": basis}, record["sequence_sha256"]]
        cid = _hash(identity)
        record["construct_id"] = cid
        validate_schema("construct", record)
        if cid in records:
            if {k: v for k, v in records[cid].items() if k not in {"warnings", "evidence"}} != {k: v for k, v in record.items() if k not in {"warnings", "evidence"}}:
                raise PublicDataError("Duplicate construct identity has conflicting values")
            for key in ("warnings", "evidence"):
                records[cid][key] = sorted(set(records[cid][key] + record[key]))
            continue
        records[cid] = record
        sequences[cid] = _sequence_record(summary, sequence, construct_id=cid, start=start, end=end, basis=basis, warnings=warnings, evidence=evidence) if state == "available" else None
    for cid, sequence_record in sequences.items():
        if sequence_record is not None:
            sequence_record["warnings"] = records[cid]["warnings"]
            sequence_record["evidence"] = records[cid]["evidence"]
    return [records[key] for key in sorted(records)], {key: sequences[key] for key in sorted(sequences)}


def project_public_target(report: dict[str, Any] | None, entry: dict[str, Any], species: str, *, public_site_dir: str | Path | None = None, static_origin: str = "https://openantigens.org") -> dict[str, Any]:
    """Project one frozen target; absent evidence stays explicitly unknown."""
    if species not in SPECIES or not isinstance(entry, dict) or (report is not None and not isinstance(report, dict)):
        raise PublicDataError("Invalid public target inputs")
    report = report or {}
    target = report.get("target") or {}
    if not isinstance(target, dict):
        raise PublicDataError("Invalid report target")
    if target.get("taxon_id") not in (None, SPECIES[species]):
        raise PublicDataError("Target species disagrees with portal")
    accession = _safe_text(target.get("accession"))
    identifier_source = _identifier_source(accession)
    target_id = f"{SPECIES[species]}:{identifier_source}:{accession}" if accession else None
    entry_name = _safe_text(target.get("entry_name") or entry.get("entry_name"))
    if target.get("entry_name") and entry.get("entry_name") != target["entry_name"]:
        raise PublicDataError("Final portal entry disagrees with report identity")
    sequence = target.get("sequence") or ""
    if not isinstance(sequence, str) or (sequence and not re.fullmatch(r"[A-Z*]+", sequence)):
        raise PublicDataError("Invalid canonical protein sequence")
    sequence_accession, sequence_source = _sequence_provenance(report, entry, species)
    aliases = sorted(set(_safe_text(x) for x in (target.get("alternative_names") or []) + (target.get("gene_synonyms") or [])))
    raw_status = _safe_text(entry.get("status"))
    topology = report.get("topology") or {}
    design_region = None
    if report.get("ectodomain"):
        region = report["ectodomain"]
        start, end = _number(region.get("start"), integer=True), _number(region.get("end"), integer=True)
        if not isinstance(start, int) or not isinstance(end, int) or not 1 <= start <= end <= len(sequence):
            raise PublicDataError("Design region does not fit recorded reference")
        design_region = {"start": start, "end": end, "label": _safe_text(region.get("label")), "source": _safe_text(region.get("source"))}
    relative_page = entry.get("detail_page")
    report_url = None
    prefix = "mouse/" if species == "mouse" else ""
    if relative_page:
        if not isinstance(relative_page, str) or not relative_page.startswith("reports/"):
            raise PublicDataError("Final portal report URL is not species-relative")
        if public_site_dir is not None:
            _safe_path(Path(public_site_dir), prefix + relative_page)
        else:
            _safe_path(Path("/private/tmp"), prefix + relative_page, exists=False)
        report_url = static_origin.rstrip("/") + "/" + prefix + relative_page
    summary = {"target_id": target_id, "record_id": _hash([species, accession, entry_name, _safe_text(entry.get("gene_symbol")), _safe_text(entry.get("protein_name"))]), "species": species, "taxon_id": SPECIES[species], "identifier_source": identifier_source, "identifier_accession": accession, "entry_name": entry_name, "gene_symbol": _safe_text(target.get("gene_symbol") or entry.get("gene_symbol")), "protein_name": _safe_text(target.get("protein_name") or entry.get("protein_name")), "aliases": aliases, "sequence_accession": sequence_accession, "sequence_source": sequence_source, "sequence_sha256": _sha(sequence.encode("ascii")) if sequence else None, "sequence_length": len(sequence) if sequence else None, "status": raw_status if raw_status in STATUSES else "unknown", "status_raw": raw_status, "track": _track(entry.get("topology_bucket")), "track_raw": _safe_text(entry.get("topology_bucket")), "analyzed_topology": _safe_text(topology.get("topology_class")), "analyzed_topology_source": _safe_text(topology.get("source")), "design_region": design_region, "construct_count": None, "structure_state": "available" if entry.get("has_any_structure") in (True, 1) else "unavailable" if entry.get("has_any_structure") in (False, 0) else "unknown", "construct_state": "unknown", "section_availability": {}, "disease_ids": [], "family_names": [], "live_report_url": report_url, "search_text": ""}
    sections: dict[str, dict[str, Any]] = {}
    summary["canonical_isoform_id"] = _safe_text(target.get("canonical_isoform_id"))
    summary["unsupported_variants"] = [{"name": _safe_text(row.get("name")), "classification": _safe_text(row.get("category")), "state": "unavailable", "reason": "Linker, cassette and replacement variants require an explicit variant export contract beyond v1."} for row in report.get("gpcr_engineering_variants") or []]
    raw_sections = {"topology": [topology] if topology else [], "domains": report.get("interpro_annotations") or [], "ptms": [], "orthologs": report.get("ectodomain_homology") or [], "paralogs": (report.get("family_context") or {}).get("members") or [], "family": [report["family_context"]] if report.get("family_context") else [], "assembly": [], "structures": [], "similarity": [], "literature": [], "disease": entry.get("open_targets_diseases") or [] if species == "human" else []}
    for key, kind in (("ptms", "feature"), ("furin_sites", "furin"), ("cysteine_analysis", "cysteine")):
        for record in report.get(key) or []:
            if kind == "cysteine":
                low, high = _report_cysteine_interval(report)
                record = {**record, "analysis_start": low, "analysis_end": high}
            raw_sections["ptms"].append({"kind": kind, "feature": None, "furin": None, "cysteine": None, kind: record})
    for requirement in report.get("assembly_requirements") or []:
        raw_sections["assembly"].append({"kind": "requirement", "requirement": requirement, "complex": None, "lookup": None})
    for complex_record in report.get("complex_portal_complexes") or []:
        raw_sections["assembly"].append({"kind": "complex", "requirement": None, "complex": complex_record, "lookup": None})
    lookup_states = []
    lookup_records = report.get("complex_portal_lookup") or {}
    if not isinstance(lookup_records, dict):
        raise PublicDataError("Invalid recorded Complex Portal lookup states")
    for lookup_species, record in lookup_records.items():
        if not isinstance(record, dict):
            raise PublicDataError("Invalid recorded Complex Portal lookup")
        status = _safe_text(record.get("status"))
        lookup_states.append(status)
        lookup = {"species": _safe_text(lookup_species), "status": status, "source": "Complex Portal", "reason": "The recorded source lookup failed." if status == "error" else "Lookup was disabled for this frozen source." if status == "disabled" else None}
        raw_sections["assembly"].append({"kind": "lookup", "requirement": None, "complex": None, "lookup": lookup})
    for row in report.get("experimental_constructs") or []:
        raw_sections["structures"].append({"kind": "pdb", "pdb_id": row.get("pdb_id"), "method": row.get("method"), "resolution": row.get("resolution"), "chains": row.get("chains") or [], "source": "PDB", "live_asset_url": f"https://www.rcsb.org/structure/{quote(str(row['pdb_id']))}" if row.get("pdb_id") else None, "sha256": None})
    if entry.get("has_alphafold_structure") in (True, 1):
        if entry.get("structure_source") != "AlphaFold":
            raise PublicDataError("Final predicted model has unsupported source provenance")
        model_url, model_sha = None, None
        if public_site_dir is not None and relative_page:
            model_relative = prefix + f"structures/{Path(relative_page).stem}.pdb"
            candidate = Path(public_site_dir) / model_relative
            if candidate.exists():
                model_file = _safe_path(Path(public_site_dir), model_relative)
                model_url = static_origin.rstrip("/") + "/" + model_relative
                model_sha = _verified_final_model(report, model_file, entry.get("structure_source"))
        raw_sections["structures"].append({"kind": "predicted_model", "pdb_id": None, "method": None, "resolution": None, "chains": [], "source": entry.get("structure_source") or None, "live_asset_url": model_url, "sha256": model_sha})
    for key, scope in (("cross_reactivity_hits", "design_region"), ("full_length_cross_reactivity_hits", "full_length")):
        for row in report.get(key) or []:
            raw_sections["similarity"].append({**row, "scope": scope})
    literature = report.get("literature_context") or {}
    literature_count = _number(entry.get("pubtator_count", literature.get("count")), integer=True)
    literature_url = entry.get("pubtator_query_url") or literature.get("query_url")
    if not literature_url and literature_count is not None and literature_count > 0 and public_site_dir is not None:
        links, _ = _frozen_pubtator_links(Path(public_site_dir) / ("mouse" if species == "mouse" else ""), [entry], selected_entries_only=True)
        literature_url = links.get(str(entry.get("entry_name") or "").lower())
    if literature_count is not None or literature_url:
        raw_sections["literature"].append({"count": literature_count, "query_url": literature_url, "source": "PubTator"})
    for section, rows in raw_sections.items():
        if not isinstance(rows, list):
            raise PublicDataError("Invalid recorded evidence section")
        if section == "disease":
            rows = [{**row, "source": "Open Targets"} for row in rows]
        projected = [_project(SECTION_TYPES[section], row) for row in rows]
        for row in projected:
            for key in ("query_url", "open_targets_url", "live_asset_url"):
                if key in row:
                    row[key] = _url(row[key], static_origin=static_origin)
        state, reason = _state(projected)
        if section == "assembly" and not (report.get("assembly_requirements") or report.get("complex_portal_complexes")):
            if "error" in lookup_states:
                state, reason = "error", "At least one recorded Complex Portal lookup failed."
            elif lookup_states and all(value == "no_hits" for value in lookup_states):
                state, reason = "no_results", "The recorded Complex Portal lookups completed without hits."
            elif lookup_states and all(value == "disabled" for value in lookup_states):
                state, reason = "unavailable", "Complex Portal lookups were disabled in this frozen source."
        if section == "disease" and species == "mouse":
            state, reason = "not_applicable", "Mouse disease associations are outside this release capability."
        sections[section] = {"state": state, "reason": reason, "items": projected}
        summary["section_availability"][section] = state
    summary["disease_ids"] = sorted({row["disease_id"] for row in sections["disease"]["items"] if row["disease_id"]})
    summary["family_names"] = sorted({name for row in sections["family"]["items"] for name in row["family_names"]})
    if entry.get("family_name") and entry["family_name"] not in summary["family_names"]:
        summary["family_names"].append(_safe_text(entry["family_name"]))
    constructs, construct_sequences = _constructs(report, summary) if target_id else ([], {})
    declared_count = _number(entry.get("construct_count"), integer=True)
    native_count = len(report.get("construct_details") or report.get("construct_recommendations") or [])
    if declared_count is not None and declared_count != native_count:
        raise PublicDataError("Final portal construct count disagrees with projected native records")
    summary["construct_count"] = len(constructs) if target_id and ("construct_details" in report or "construct_recommendations" in report) else declared_count
    summary["portal_construct_count"] = declared_count
    summary["construct_state"] = "available" if summary["construct_count"] else "unavailable" if summary["construct_count"] == 0 else "unknown"
    search_parts = [summary[key] or "" for key in ("identifier_accession", "entry_name", "gene_symbol", "protein_name", "track_raw", "analyzed_topology")]
    summary["search_text"] = " ".join(search_parts + aliases + summary["family_names"]).lower()
    validate_schema("summary", summary)
    reference_hazards = _hazards(report, {"start": 1, "end": len(sequence)}) if sequence else None
    reference = _sequence_record(summary, sequence, warnings=_hazard_warnings(reference_hazards, start=1, end=len(sequence))) if target_id and sequence else None
    projected = {"summary": summary, "sections": sections, "constructs": constructs, "reference_sequence": reference, "construct_sequences": construct_sequences}
    _scan_public(projected)
    return projected


def _envelope(release_id: str, species: str, data: Any) -> dict[str, Any]:
    return {"release_id": release_id, "schema_version": SCHEMA_VERSION, "species": species, "data": data}


def _file_meta(content: bytes, media_type: str, schema: str | None = None, encoding: str | None = None) -> dict[str, Any]:
    return {"sha256": _sha(content), "bytes": len(content), "media_type": media_type, "encoding": encoding, "schema": schema}


def _write(root: Path, relative: str, content: bytes) -> None:
    path = _safe_path(root, relative, exists=False)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(content)


def _write_json(root: Path, relative: str, value: Any) -> bytes:
    content = canonical_bytes(value) + b"\n"
    _write(root, relative, content)
    return content


def _atomic_public_replace(root: Path, relative: str, content: bytes) -> None:
    path = _safe_path(root, relative)
    descriptor, temporary = tempfile.mkstemp(prefix=".public-metadata-", suffix=".tmp", dir=path.parent)
    temp = Path(temporary)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        _safe_path(root, relative)
        os.replace(temp, path)
    finally:
        if temp.exists():
            temp.unlink()


def _update_portal_metadata(public_site: Path, prepared: dict[str, tuple[bytes, dict[str, Any]]], release_id: str, static_origin: str) -> None:
    updates = {}
    for species, (old_bytes, metadata) in prepared.items():
        relative = ("mouse/" if species == "mouse" else "") + "portal_metadata.json"
        if _safe_path(public_site, relative).read_bytes() != old_bytes:
            raise PublicDataError("Portal metadata changed while compiling the public bundle")
        updated = {**metadata, "public_api": {"schema_version": SCHEMA_VERSION, "release_id": release_id, "state": "packaged", "manifest_url": static_origin.rstrip("/") + f"/data/releases/{release_id}/manifest.json", "downloads_url": static_origin.rstrip("/") + f"/data/releases/{release_id}/{species}/index.jsonl"}}
        updates[relative] = (old_bytes, canonical_bytes(updated) + b"\n")
    completed = []
    try:
        for relative, (_, content) in updates.items():
            _atomic_public_replace(public_site, relative, content)
            completed.append(relative)
    except Exception:
        for relative in reversed(completed):
            previous, content = updates[relative]
            if _safe_path(public_site, relative).read_bytes() != content:
                raise PublicDataError("Portal metadata changed during atomic rollback")
            _atomic_public_replace(public_site, relative, previous)
        raise


def _page_records(records: list[Any], make_envelope: Any) -> Iterable[list[Any]]:
    page: list[Any] = []
    for record in records:
        if len(canonical_bytes(make_envelope([record]))) + 1 > SHARD_LIMIT:
            raise PublicDataError("A public record exceeds the evidence shard budget")
        if page and (len(page) == PAGE_SIZE or len(canonical_bytes(make_envelope(page + [record]))) + 1 > SHARD_LIMIT):
            yield page
            page = []
        page.append(record)
    if page:
        yield page


def _tsv_chunks(rows: Iterable[dict[str, Any]], fields: tuple[str, ...], *, include_header: bool = True) -> Iterable[bytes]:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
    if include_header:
        writer.writeheader()
        yield output.getvalue().encode("utf-8")
        output.seek(0)
        output.truncate()
    for row in rows:
        values = {}
        for key in fields:
            value = row.get(key)
            if isinstance(value, (list, dict)):
                value = canonical_bytes(value).decode("utf-8")
            if isinstance(value, str) and value.lstrip().startswith(("=", "+", "-", "@")):
                value = "'" + value
            values[key] = value
        writer.writerow(values)
        yield output.getvalue().encode("utf-8")
        output.seek(0)
        output.truncate()


def sequence_fasta(record: dict[str, Any], *, release_id: str, provenance_url: str, records_url: str) -> str:
    """One exact released sequence, with machine-safe identity fields and LF wrapping."""
    validate_schema("sequence", record)
    if not re.fullmatch(r"[a-f0-9]{64}", release_id):
        raise PublicDataError("FASTA requires a literal public release ID")
    provenance = urlsplit(provenance_url)
    records = urlsplit(records_url)
    if provenance.scheme != "https" or records.scheme != "https" or provenance.netloc != records.netloc or provenance.username or provenance.password or records.username or records.password or provenance.path != f"/data/releases/{release_id}/manifest.json" or records.path != f"/data/releases/{release_id}/{record['species']}/records.jsonl.gz":
        raise PublicDataError("FASTA requires immutable provenance and complete-record references")
    fields = {"release_id": release_id, **{key: record[key] for key in ("target_id", "construct_id", "species", "taxon_id", "identifier_source", "identifier_accession", "sequence_accession", "sequence_source", "sequence_basis", "start", "end", "length")}, "sha256": record["sequence_sha256"], "warnings": canonical_bytes(record["warnings"]).decode("utf-8"), "evidence": canonical_bytes(record["evidence"]).decode("utf-8"), "provenance_url": provenance_url, "records_url": records_url, "record_selector": canonical_bytes({"target_id": record["target_id"], "construct_id": record["construct_id"]}).decode("utf-8")}
    header = " ".join(f"{key}={quote('null' if value is None else str(value), safe='')}" for key, value in fields.items())
    sequence = record["sequence"]
    primary = quote(f"{release_id}:{record['target_id']}:{record['construct_id'] or 'reference'}", safe="")
    return f">{primary} {header}\n" + "\n".join(sequence[i:i + 80] for i in range(0, len(sequence), 80)) + "\n"


def _load_sources(summary_path: Path, species: str, public_site: Path, static_origin: str) -> tuple[list[dict[str, Any]], dict[str, str]]:
    if summary_path.is_symlink() or not summary_path.is_file():
        raise PublicDataError("Missing or symlinked frozen summary")
    summary_rows, summary_hash = _read_json_hashed(summary_path)
    if not isinstance(summary_rows, list):
        raise PublicDataError("Frozen summary must contain target rows")
    portal_root = public_site / ("mouse" if species == "mouse" else "")
    index, index_hash = _read_json_hashed(_safe_path(portal_root, "downloads/agdesign2_portal_index.json"))
    if not isinstance(index, list) or not index:
        raise PublicDataError("Final species index must be nonempty")
    index_by_entry = {str(row.get("entry_name") or "").upper(): row for row in index}
    if len(index_by_entry) != len(index):
        raise PublicDataError("Duplicate final portal entry")
    hashes = {f"{species}:summary": summary_hash, f"{species}:portal_index": index_hash}
    literature_links, literature_hash = _frozen_pubtator_links(portal_root, index)
    if literature_hash:
        hashes[f"{species}:portal_index_script"] = literature_hash
    disease_by_entry: dict[str, list[dict[str, Any]]] = {}
    disease_tsv = portal_root / "open_targets_disease_associations.tsv"
    if species == "human" and disease_tsv.is_file():
        if disease_tsv.is_symlink():
            raise PublicDataError("Disease artifact cannot be a symlink")
        hashes["human:disease_tsv"] = _sha(disease_tsv.read_bytes())
        with disease_tsv.open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                disease_by_entry.setdefault(str(row.get("entry_name") or "").upper(), []).append(row)
    projected = []
    seen = set()
    root = summary_path.parent.resolve()
    for position, row in enumerate(summary_rows):
        if not isinstance(row, dict):
            raise PublicDataError("Invalid frozen summary row")
        report = None
        if row.get("json_report"):
            path = Path(str(row["json_report"]))
            path = path if path.is_absolute() else root / path
            if not path.absolute().is_relative_to(root):
                raise PublicDataError("Frozen report escapes its declared snapshot")
            path = _safe_path(root, path.absolute().relative_to(root).as_posix())
            report, report_hash = _read_json_hashed(path, maximum=64 * 1024 * 1024)
            if not isinstance(report, dict):
                raise PublicDataError("Invalid frozen report")
            hashes[f"{species}:report:{position}"] = report_hash
        key = str(((report or {}).get("target") or {}).get("entry_name") or row.get("resolved_entry_name") or row.get("query") or "").upper()
        if key in seen or key not in index_by_entry:
            raise PublicDataError("Duplicate or unmatched frozen source identity")
        seen.add(key)
        entry = index_by_entry[key]
        if species == "human":
            entry = {**entry, "open_targets_diseases": disease_by_entry.get(key, [])}
        if key.lower() in literature_links:
            entry = {**entry, "pubtator_query_url": literature_links[key.lower()]}
        projected.append(project_public_target(report, entry, species, public_site_dir=public_site, static_origin=static_origin))
        del report
    if seen != set(index_by_entry):
        raise PublicDataError("Frozen summary and finalized portal catalogue differ")
    projected.sort(key=lambda row: (row["summary"]["target_id"] or "", row["summary"]["record_id"]))
    ids = [row["summary"]["target_id"] for row in projected if row["summary"]["target_id"]]
    if len(ids) != len(set(ids)):
        raise PublicDataError("Duplicate species/accession public target identity")
    return projected, hashes


def build_public_api_bundle(human_summary_path: str | Path, mouse_summary_path: str | Path, output_dir: str | Path, *, public_site_dir: str | Path, software_commit: str, snapshot_id: str, previous_bundle_dir: str | Path | None = None, static_origin: str = "https://openantigens.org") -> Path:
    """Build current plus optional previous release; refuse existing output content."""
    if Path(output_dir).is_symlink() or Path(public_site_dir).is_symlink():
        raise PublicDataError("Public build roots cannot be symbolic links")
    output, public_site = Path(output_dir).resolve(), Path(public_site_dir).resolve()
    origin = urlsplit(static_origin)
    if origin.scheme != "https" or not origin.hostname or origin.path not in {"", "/"} or origin.query or origin.fragment or origin.username or origin.password:
        raise PublicDataError("Static origin must be a plain HTTPS origin")
    if output.exists() and (output.is_symlink() or any(output.iterdir())):
        raise PublicDataError("Public bundle destination must be absent or empty")
    if not public_site.is_dir():
        raise PublicDataError("Finalized public site is required")
    prepared_metadata = {}
    for species in SPECIES:
        relative = ("mouse/" if species == "mouse" else "") + "portal_metadata.json"
        metadata_path = _safe_path(public_site, relative)
        raw = metadata_path.read_bytes()
        metadata = _read_json(metadata_path)
        if not isinstance(metadata, dict):
            raise PublicDataError("Final portal metadata must be an object")
        prepared_metadata[species] = (raw, metadata)
    _safe_text(software_commit)
    _safe_text(snapshot_id)
    projected, inputs = {}, {}
    for species, summary in (("human", human_summary_path), ("mouse", mouse_summary_path)):
        projected[species], hashes = _load_sources(Path(summary).resolve(), species, public_site, static_origin)
        inputs.update(hashes)
    citation = _project("citation", _read_json(Path(__file__).parent / "data/openantigens_citation.json"))
    provenance = {"software_commit": software_commit, "snapshot_id": snapshot_id, "source_versions": {}, "source_dates": {}, "input_sha256": inputs, "notes": ["Historical source versions and retrieval dates are unknown unless recorded in the frozen input."]}
    validate_schema("provenance", provenance)
    asset_hashes = {}
    for rows in projected.values():
        for record in rows:
            for asset in record["sections"]["structures"]["items"]:
                if asset["sha256"] and asset["live_asset_url"]:
                    asset_hashes[urlsplit(asset["live_asset_url"]).path.lstrip("/")] = asset["sha256"]
    core = {"schema_version": SCHEMA_VERSION, "provenance": provenance, "payload_sha256": {**{species: semantic_hash(rows) for species, rows in projected.items()}, "schema": semantic_hash(PUBLIC_SCHEMA), "citation": semantic_hash(citation)}, "asset_sha256": asset_hashes}
    release_id = _hash(core)
    prefix = f"releases/{release_id}"
    files: dict[str, Any] = {}
    species_manifests = {}
    output.mkdir(parents=True, exist_ok=True)

    def write_json(relative: str, value: Any, schema: str, limit: int = SHARD_LIMIT) -> None:
        validate_schema(schema, value)
        _scan_public(value)
        content = canonical_bytes(value) + b"\n"
        if len(content) > limit:
            raise PublicDataError("Public shard exceeds its storage budget")
        _write(output, relative, content)
        files[relative] = _file_meta(content, "application/json", schema)

    for species, rows in projected.items():
        species_prefix = f"{prefix}/{species}"
        shard_paths: list[tuple[str, str]] = []
        target_paths: dict[str, str] = {}
        all_constructs, all_sequences = [], []
        summaries = [row["summary"] for row in rows]
        index_path = f"{species_prefix}/index.json"
        write_json(index_path, _envelope(release_id, species, summaries), "index", INDEX_LIMIT)
        shard_paths.append(("index", index_path))
        for projected_target in rows:
            summary = projected_target["summary"]
            target_id = summary["target_id"]
            key = _hash(target_id or summary["record_id"])
            target_path = f"{species_prefix}/targets/{key}.json"
            descriptors = {}
            for section, evidence in projected_target["sections"].items():
                descriptor = {"state": evidence["state"], "reason": evidence["reason"], "total": len(evidence["items"]), "pages": []}
                make = lambda items, section=section, evidence=evidence: _envelope(release_id, species, {"target_id": target_id, "section": section, "state": evidence["state"], "reason": evidence["reason"], "items": items})
                for page_number, page in enumerate(_page_records(evidence["items"], make)):
                    path = f"{species_prefix}/sections/{key}/{section}-{page_number}.json"
                    write_json(path, make(page), f"section_{section}")
                    descriptor["pages"].append(path)
                    shard_paths.append(("section", path))
                descriptors[section] = descriptor
            constructs = projected_target["constructs"]
            cstate, creason = _state(constructs)
            if not constructs and summary["construct_count"] == 0:
                cstate, creason = "no_results", "The source recorded zero released constructs."
            construct_descriptor = {"state": cstate, "reason": creason, "total": len(constructs), "pages": []}
            for page_number, page in enumerate(_page_records(constructs, lambda items: _envelope(release_id, species, items))):
                path = f"{species_prefix}/constructs/{key}-{page_number}.json"
                write_json(path, _envelope(release_id, species, page), "construct_page")
                construct_descriptor["pages"].append(path)
                shard_paths.append(("constructs", path))
            all_constructs.extend(constructs)
            sequence_paths = {}
            reference_path = None
            reference = projected_target["reference_sequence"]
            if reference:
                reference_path = f"{species_prefix}/sequences/{key}-reference.json"
                write_json(reference_path, _envelope(release_id, species, reference), "sequence_shard")
                shard_paths.append(("sequence", reference_path))
                all_sequences.append(reference)
            for cid, sequence_record in projected_target["construct_sequences"].items():
                sequence_paths[cid] = None
                if sequence_record is not None:
                    path = f"{species_prefix}/sequences/{cid}.json"
                    write_json(path, _envelope(release_id, species, sequence_record), "sequence_shard")
                    sequence_paths[cid] = path
                    shard_paths.append(("sequence", path))
                    all_sequences.append(sequence_record)
            value = {"summary": summary, "sections": descriptors, "constructs": construct_descriptor, "reference_sequence_path": reference_path, "construct_sequence_paths": sequence_paths}
            write_json(target_path, _envelope(release_id, species, value), "target_shard", TARGET_LIMIT)
            if target_id is not None:
                target_paths[target_id] = target_path
            shard_paths.append(("target", target_path))
        downloads = []

        def download(name: str, content: bytes | Iterable[bytes], media: str, count: int | None, schema: str | None = None, encoding: str | None = None) -> dict[str, Any]:
            path = f"{species_prefix}/downloads/{name}"
            destination = _safe_path(output, path, exists=False)
            destination.parent.mkdir(parents=True, exist_ok=True)
            chunks = [content] if isinstance(content, bytes) else content
            with destination.open("wb") as handle:
                if encoding == "gzip":
                    with gzip.GzipFile(filename="", mode="wb", fileobj=handle, compresslevel=6, mtime=0) as compressor:
                        for chunk in chunks:
                            compressor.write(chunk)
                else:
                    for chunk in chunks:
                        handle.write(chunk)
            metadata = {"sha256": _sha_file(destination), "bytes": destination.stat().st_size, "media_type": media, "encoding": encoding, "schema": schema}
            files[path] = metadata
            result = {"name": name, "url": static_origin.rstrip("/") + f"/data/releases/{release_id}/{species}/{name}", "path": path, **metadata, "row_count": count, "license": LICENSE}
            validate_schema("download", result)
            downloads.append(result)
            return result

        jsonl = lambda records: (canonical_bytes(row) + b"\n" for row in records)
        download("index.jsonl", jsonl(summaries), "application/x-ndjson", len(summaries), "summary")
        download("index.tsv", _tsv_chunks(summaries, tuple(DEFS["summary"]["properties"])), "text/tab-separated-values", len(summaries))
        download("constructs.jsonl", jsonl(all_constructs), "application/x-ndjson", len(all_constructs), "construct")
        download("constructs.tsv", _tsv_chunks(all_constructs, tuple(DEFS["construct"]["properties"])), "text/tab-separated-values", len(all_constructs))
        for name, selected in (("constructs.fasta", [row for row in all_sequences if row["construct_id"]]), ("reference.fasta", [row for row in all_sequences if row["construct_id"] is None])):
            download(name, (sequence_fasta(record, release_id=release_id, provenance_url=static_origin.rstrip("/") + f"/data/releases/{release_id}/manifest.json", records_url=static_origin.rstrip("/") + f"/data/releases/{release_id}/{species}/records.jsonl.gz").encode("utf-8") for record in selected), "text/x-fasta", len(selected))
        download("schema.json", canonical_bytes(PUBLIC_SCHEMA) + b"\n", "application/schema+json", None)
        download("citation.json", canonical_bytes(citation) + b"\n", "application/json", None, "citation")
        from .site_docs import citation_downloads
        for name, content in citation_downloads().items():
            download(name, content.encode("utf-8"), "application/x-bibtex" if name.endswith(".bib") else "application/x-research-info-systems", None)
        records = ({"kind": kind, "path": path, "payload": _read_json(output / path)} for kind, path in shard_paths)
        records_download = download("records.jsonl.gz", jsonl(records), "application/x-ndjson", len(shard_paths), None, "gzip")
        human_disease = species == "human" and (public_site / "open_targets_disease_associations.tsv").is_file()
        species_manifests[species] = {"taxon_id": SPECIES[species], "index_path": index_path, "target_paths": target_paths, "downloads": downloads, "counts": {"targets": len(rows), "resolved_targets": len(target_paths), "constructs": len(all_constructs), "status": dict(sorted(Counter(row["status"] for row in summaries).items()))}, "capabilities": {"sections": [section for section in SECTIONS if section != "disease" or human_disease], "disease": human_disease}, "records_download": records_download}
    manifest = {"release_id": release_id, "schema_version": SCHEMA_VERSION, "manifest_url": static_origin.rstrip("/") + f"/data/releases/{release_id}/manifest.json", "provenance": provenance, "citation": citation, "licenses": LICENSES, "species": species_manifests, "files": files, "core": core}
    validate_schema("release_manifest", manifest)
    manifest_path = f"{prefix}/manifest.json"
    manifest_bytes = _write_json(output, manifest_path, manifest)
    if len(manifest_bytes) > MANIFEST_LIMIT:
        raise PublicDataError("Release manifest exceeds the validated service manifest budget")
    release_pointer = {"release_id": release_id, "manifest_path": manifest_path, "sha256": _sha(manifest_bytes)}
    releases = [release_pointer]
    retired_releases = []
    if previous_bundle_dir is not None:
        previous = Path(previous_bundle_dir).resolve()
        old_root = validate_public_api_bundle(previous)
        newly_retired = []
        for pointer in old_root["releases"]:
            if pointer["release_id"] in {item["release_id"] for item in releases}:
                continue
            if len(releases) == 2:
                newly_retired.append({"release_id": pointer["release_id"], "retired_by_release_id": release_id, "archive_url": None})
                continue
            old_manifest = _read_json(_safe_path(previous, pointer["manifest_path"]))
            for path in [pointer["manifest_path"], *old_manifest["files"]]:
                _write(output, path, _safe_path(previous, path).read_bytes())
            releases.append(pointer)
        retained_ids = {item["release_id"] for item in releases}
        seen_retired = set()
        for tombstone in [*newly_retired, *old_root.get("retired_releases", [])]:
            if tombstone["release_id"] not in retained_ids | seen_retired:
                retired_releases.append(tombstone)
                seen_retired.add(tombstone["release_id"])
        retired_releases = retired_releases[:1000]
    root_manifest = {"schema_version": SCHEMA_VERSION, "current_release_id": release_id, "releases": releases, "retired_releases": retired_releases}
    if len(canonical_bytes(root_manifest)) + 1 > SHARD_LIMIT:
        raise PublicDataError("Public root manifest exceeds its bounded service budget")
    root_bytes = _write_json(output, "manifest.json", root_manifest)
    _write(output, "manifest.sha256", f"{_sha(root_bytes)}  manifest.json\n".encode("ascii"))
    validate_public_api_bundle(output)
    _copy_static_downloads(output, public_site)
    _update_portal_metadata(public_site, prepared_metadata, release_id, static_origin)
    return output / "manifest.json"


def _copy_static_downloads(bundle: Path, public_site: Path) -> None:
    root = _read_json(bundle / "manifest.json")
    for release in root["releases"]:
        manifest = _read_json(bundle / release["manifest_path"])
        rid = release["release_id"]
        for species_manifest in manifest["species"].values():
            for entry in species_manifest["downloads"]:
                path = urlsplit(entry["url"]).path.lstrip("/")
                source = _safe_path(bundle, entry["path"])
                destination = _safe_path(public_site, path, exists=False)
                if destination.exists() and _sha_file(destination) != entry["sha256"]:
                    raise PublicDataError("Immutable public download path already contains different bytes")
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(source, destination)
        content = _safe_path(bundle, release["manifest_path"]).read_bytes()
        relative = f"data/releases/{rid}/manifest.json"
        destination = _safe_path(public_site, relative, exists=False)
        if destination.exists() and destination.read_bytes() != content:
            raise PublicDataError("Immutable release manifest path already contains different bytes")
        _write(public_site, relative, content)
        _write(public_site, f"data/releases/{rid}/manifest.sha256", f"{_sha(content)}  manifest.json\n".encode("ascii"))


def _approved_release_schema(schema: Any) -> int:
    """Accept only offline, closed schemas compatible with the trusted public allowlist."""
    if not isinstance(schema, dict) or set(schema) - {"$schema", "$id", "$defs", "x-openantigens-fasta-version"} or schema.get("$schema") != PUBLIC_SCHEMA["$schema"] or schema.get("$id") != PUBLIC_SCHEMA["$id"] or not isinstance(schema.get("$defs"), dict):
        raise PublicDataError("Unsupported historical public schema")
    version = schema.get("x-openantigens-fasta-version", 1)
    if version not in (1, 2):
        raise PublicDataError("Unsupported historical FASTA contract")

    def compatible(old: Any, current: Any) -> None:
        if isinstance(old, dict):
            if not isinstance(current, dict) or set(old) != set(current):
                raise PublicDataError("Unsupported historical schema shape")
            for keyword, value in old.items():
                if keyword in {"$dynamicRef", "$recursiveRef"}:
                    raise PublicDataError("Dynamic public schema references are prohibited")
                if keyword == "$ref":
                    if not isinstance(value, str) or not value.startswith("#/$defs/") or value[len("#/$defs/"):] not in schema["$defs"] or value != current[keyword]:
                        raise PublicDataError("Nonlocal or incompatible public schema reference")
                elif keyword == "properties":
                    if not isinstance(value, dict) or set(value) - set(current[keyword]):
                        raise PublicDataError("Historical schema allows unknown public fields")
                    for field, child in value.items():
                        compatible(child, current[keyword][field])
                elif keyword == "required":
                    if set(value) != set(current[keyword]) & set(old.get("properties", {})):
                        raise PublicDataError("Historical schema weakens required public fields")
                else:
                    compatible(value, current[keyword])
        elif isinstance(old, list):
            if not isinstance(current, list) or len(old) != len(current):
                raise PublicDataError("Unsupported historical schema alternatives")
            for a, b in zip(old, current):
                compatible(a, b)
        elif type(old) is not type(current) or old != current:
            raise PublicDataError("Unsupported historical public schema value")

    for name, definition in schema["$defs"].items():
        if name not in DEFS:
            raise PublicDataError("Unknown historical public schema definition")
        compatible(definition, DEFS[name])
    for required in ("summary", "target_shard", "index", "construct", "construct_page", "sequence", "sequence_shard", "citation", "release_manifest", *[f"section_{key}" for key in SECTIONS]):
        if required not in schema["$defs"]:
            raise PublicDataError("Incomplete historical public schema")
    return version


def _validate_release_value(schema: dict[str, Any], name: str, value: Any) -> None:
    if name not in schema["$defs"]:
        raise PublicDataError("Unknown release-local public schema")
    error = next(Draft202012Validator({**schema, "$ref": f"#/$defs/{name}"}).iter_errors(value), None)
    if error is not None:
        raise PublicDataError("Public data does not match its approved release-local schema")


def _legacy_fasta(record: dict[str, Any]) -> str:
    fields = {key: record[key] for key in ("target_id", "construct_id", "species", "sequence_accession", "start", "end")}
    header = " ".join(f"{key}={quote(str(value), safe=':._-')}" for key, value in fields.items() if value is not None)
    sequence = record["sequence"]
    return f">{header} sha256={record['sequence_sha256']}\n" + "\n".join(sequence[i:i + 80] for i in range(0, len(sequence), 80)) + "\n"


class _ShardFiles:
    """Validated file references; only compact indexes stay resident across targets."""

    def __init__(self, root: Path, files: dict[str, Any], release_id: str):
        self.root, self.files, self.release_id = root, files, release_id
        self.paths: set[str] = set()
        self.indexes: dict[str, Any] = {}

    def add(self, relative: str, data: Any) -> None:
        self.paths.add(relative)
        if self.files[relative]["schema"] == "index":
            self.indexes[relative] = data

    def __len__(self) -> int:
        return len(self.paths)

    def __contains__(self, relative: str) -> bool:
        return relative in self.paths

    def __getitem__(self, relative: str) -> Any:
        if relative not in self.paths:
            raise KeyError(relative)
        if relative in self.indexes:
            return self.indexes[relative]
        info = self.files[relative]
        bound = TARGET_LIMIT if info["schema"] == "target_shard" else SHARD_LIMIT
        file = _safe_path(self.root, relative)
        with file.open("rb") as handle:
            content = handle.read(bound + 1)
        if len(content) > bound or len(content) != info["bytes"] or _sha(content) != info["sha256"]:
            raise PublicDataError("Public shard changed after manifest verification")
        data = json.loads(content)
        # These exact bytes were schema/privacy checked in the manifest pass.
        if data["release_id"] != self.release_id or f"/{data['species']}/" not in relative:
            raise PublicDataError("Public shard release/species mismatch")
        return data

    def get(self, relative: str, default: Any = None) -> Any:
        return self[relative] if relative in self.paths else default

def validate_public_api_bundle(bundle_dir: str | Path) -> dict[str, Any]:
    """Verify schemas, hashes, all references and complete overflow-download parity."""
    root = Path(bundle_dir).resolve()
    manifest_bytes = _safe_path(root, "manifest.json").read_bytes()
    if len(manifest_bytes) > SHARD_LIMIT:
        raise PublicDataError("Public root manifest exceeds its bounded service budget")
    if manifest_bytes.startswith(b"\x1f\x8b"):
        raise PublicDataError("Public root manifest must be plain bounded JSON")
    checksum = _safe_path(root, "manifest.sha256").read_text(encoding="ascii")
    if checksum != f"{_sha(manifest_bytes)}  manifest.json\n":
        raise PublicDataError("Root public manifest checksum mismatch")
    manifest = _read_json(root / "manifest.json")
    validate_schema("root_manifest", manifest)
    release_ids = [row["release_id"] for row in manifest["releases"]]
    if not 1 <= len(release_ids) <= 2 or len(set(release_ids)) != len(release_ids) or manifest["current_release_id"] not in release_ids:
        raise PublicDataError("Invalid current/retained public releases")
    retired_ids = set()
    for tombstone in manifest.get("retired_releases", []):
        retired_id = tombstone["release_id"]
        if retired_id in retired_ids or retired_id in release_ids or retired_id == tombstone["retired_by_release_id"]:
            raise PublicDataError("Invalid or duplicate retired public release identity")
        retired_ids.add(retired_id)
        if tombstone["archive_url"] is not None:
            archive = urlsplit(tombstone["archive_url"])
            if archive.scheme != "https" or not archive.hostname or archive.username or archive.password:
                raise PublicDataError("Unsafe recorded retired-release archive URL")
    allowed = {"manifest.json", "manifest.sha256"}
    total_index_bytes = 0
    for pointer in manifest["releases"]:
        path = _safe_path(root, pointer["manifest_path"])
        if path.stat().st_size > MANIFEST_LIMIT:
            raise PublicDataError("Release manifest exceeds its validated service budget")
        release_bytes = path.read_bytes()
        if release_bytes.startswith(b"\x1f\x8b"):
            raise PublicDataError("Public release manifest must be plain bounded JSON")
        if _sha(release_bytes) != pointer["sha256"]:
            raise PublicDataError("Release public manifest checksum mismatch")
        release = json.loads(release_bytes)
        validate_schema("release_manifest", release)
        _scan_public(release)
        rid = pointer["release_id"]
        if release["release_id"] != rid or _hash(release["core"]) != rid or release["provenance"] != release["core"]["provenance"]:
            raise PublicDataError("Public release identity mismatch")
        manifest_url = urlsplit(release["manifest_url"])
        if manifest_url.scheme != "https" or manifest_url.username or manifest_url.password or manifest_url.query or manifest_url.fragment or manifest_url.path != f"/data/releases/{rid}/manifest.json":
            raise PublicDataError("Invalid immutable public manifest URL")
        expected_prefix = f"releases/{rid}/"
        if pointer["manifest_path"] != expected_prefix + "manifest.json":
            raise PublicDataError("Release manifest is outside its immutable namespace")
        allowed.add(pointer["manifest_path"])
        schema_paths = [relative for relative, info in release["files"].items() if info["media_type"] == "application/schema+json"]
        if not schema_paths:
            raise PublicDataError("Release-local public schema is missing")
        schema_path = _safe_path(root, schema_paths[0])
        if schema_path.stat().st_size > SHARD_LIMIT:
            raise PublicDataError("Release schema exceeds service budget")
        release_schema = json.loads(schema_path.read_bytes())
        fasta_version = _approved_release_schema(release_schema)
        _validate_release_value(release_schema, "release_manifest", release)
        stored = _ShardFiles(root, release["files"], rid)
        schema_payloads = []
        for relative, info in release["files"].items():
            if not relative.startswith(expected_prefix):
                raise PublicDataError("Public file belongs to another release")
            file = _safe_path(root, relative)
            maximum = INDEX_LIMIT if info["schema"] == "index" else TARGET_LIMIT if info["schema"] == "target_shard" else 128 * 1024 * 1024 if "/downloads/" in relative and info["media_type"] != "application/schema+json" else SHARD_LIMIT
            if file.stat().st_size > maximum:
                raise PublicDataError("Public file exceeds the validated service budget")
            if info["schema"] == "index":
                total_index_bytes += file.stat().st_size
                if total_index_bytes > INDEX_TOTAL_LIMIT:
                    raise PublicDataError("Combined retained public indexes exceed service memory budget")
            parse_json = info["media_type"] in {"application/json", "application/schema+json"}
            content = file.read_bytes() if parse_json else None
            checksum = _sha(content) if content is not None else _sha_file(file)
            if file.stat().st_size != info["bytes"] or checksum != info["sha256"]:
                raise PublicDataError("Public file bytes or checksum mismatch")
            allowed.add(relative)
            if info["media_type"] == "application/schema+json":
                schema_payloads.append(json.loads(content))
                if content != schema_path.read_bytes():
                    raise PublicDataError("Public release schema copies disagree")
            if info["schema"]:
                if info["media_type"] not in {"application/json", "application/x-ndjson"}:
                    raise PublicDataError("Unsupported typed public file media type")
                if info["media_type"] == "application/json" and (info["encoding"] is not None or content.startswith(b"\x1f\x8b")):
                    raise PublicDataError("Public JSON shards require plain bounded JSON bytes")
                data = json.loads(content) if info["media_type"] == "application/json" else None
                if data is not None:
                    _validate_release_value(release_schema, info["schema"], data)
                    _scan_public(data)
                    if info["schema"] in {"index", "target_shard", "construct_page", "sequence_shard"} or info["schema"].startswith("section_"):
                        if data["release_id"] != rid or f"/{data['species']}/" not in relative:
                            raise PublicDataError("Shard release/species mismatch")
                        stored.add(relative, data)
                        maximum = INDEX_LIMIT if info["schema"] == "index" else TARGET_LIMIT if info["schema"] == "target_shard" else SHARD_LIMIT
                        if len(content) > maximum:
                            raise PublicDataError("Public storage budget exceeded")
                        if info["schema"] == "sequence_shard":
                            sequence = data["data"]
                            if sequence["length"] != len(sequence["sequence"]) or sequence["sequence_sha256"] != _sha(sequence["sequence"].encode("ascii")):
                                raise PublicDataError("Public sequence identity mismatch")
        if not schema_payloads or any(semantic_hash(value) != release["core"]["payload_sha256"].get("schema") for value in schema_payloads):
            raise PublicDataError("Public schema differs from its canonical payload identity")
        if semantic_hash(release["citation"]) != release["core"]["payload_sha256"].get("citation"):
            raise PublicDataError("Public citation differs from its canonical payload identity")
        asset_hashes = {}
        for species, species_manifest in release["species"].items():
            if species_manifest["taxon_id"] != SPECIES[species]:
                raise PublicDataError("Public species taxon mismatch")
            index = stored.get(species_manifest["index_path"])
            if index is None or index["species"] != species:
                raise PublicDataError("Public species index is missing")
            rows = index["data"]
            ids = {row["target_id"] for row in rows if row["target_id"]}
            if len(rows) != species_manifest["counts"]["targets"] or ids != set(species_manifest["target_paths"]) or len(ids) != sum(row["target_id"] is not None for row in rows):
                raise PublicDataError("Public catalogue and manifest disagree")
            if rows != sorted(rows, key=lambda row: (row["target_id"] or "", row["record_id"])):
                raise PublicDataError("Public index stable order mismatch")
            if len(ids) != species_manifest["counts"]["resolved_targets"] or dict(Counter(row["status"] for row in rows)) != species_manifest["counts"]["status"]:
                raise PublicDataError("Public species aggregate counts disagree")
            for row in rows:
                if row["taxon_id"] != SPECIES[species] or row["species"] != species:
                    raise PublicDataError("Public target species mismatch")
                if row["target_id"] is not None and row["target_id"] != f"{SPECIES[species]}:{row['identifier_source']}:{row['identifier_accession']}":
                    raise PublicDataError("Public target ID mismatch")
            # Target/section/construct/sequence cross-references are checked together
            # in the per-target semantic pass below, without loading a species twice.
            for download in species_manifest["downloads"]:
                metadata = release["files"].get(download["path"])
                if metadata is None or any(download[key] != metadata[key] for key in ("sha256", "bytes", "media_type", "encoding", "schema")):
                    raise PublicDataError("Public download metadata mismatch")
                url = urlsplit(download["url"])
                if url.scheme != "https" or url.netloc != manifest_url.netloc or url.username or url.password or url.query or url.fragment or url.path != f"/data/releases/{rid}/{species}/{download['name']}":
                    raise PublicDataError("Public download URL is outside immutable release")
            fallback = species_manifest["records_download"]
            if fallback not in species_manifest["downloads"] or fallback["name"] != "records.jsonl.gz":
                raise PublicDataError("Complete overflow download is missing")
            seen = set()
            with gzip.open(_safe_path(root, fallback["path"]), "rb") as handle:
                while True:
                    line = handle.readline(INDEX_LIMIT + 4096)
                    if not line:
                        break
                    if not line.endswith(b"\n") or len(seen) >= len(stored):
                        raise PublicDataError("Overflow download exceeds its declared record bounds")
                    value = json.loads(line)
                    if not isinstance(value, dict) or set(value) != {"kind", "path", "payload"} or value["path"] in seen or value["path"] not in stored or value["payload"] != stored[value["path"]]:
                        raise PublicDataError("Overflow download differs from its public shards")
                    if value["payload"]["species"] != species:
                        raise PublicDataError("Overflow download mixes species")
                    seen.add(value["path"])
            expected = {relative for relative in stored.paths if f"/{species}/" in relative}
            if seen != expected or len(seen) != fallback["row_count"]:
                raise PublicDataError("Overflow download is incomplete")
            target_records = {}
            for relative in stored.paths:
                if release["files"][relative]["schema"] == "target_shard" and f"/{species}/" in relative:
                    value = stored[relative]
                    record_id = value["data"]["summary"]["record_id"]
                    if record_id in target_records:
                        raise PublicDataError("Duplicate public target record")
                    target_records[record_id] = relative
            if set(target_records) != {row["record_id"] for row in rows}:
                raise PublicDataError("Public target shards do not cover the complete index")
            payload_digest = hashlib.sha256()
            payload_digest.update(b"[")
            construct_total = 0
            export_hashes = {name: hashlib.sha256() for name in ("index.jsonl", "index.tsv", "constructs.jsonl", "constructs.tsv", "reference.fasta", "constructs.fasta", "citation.json")}
            for row in rows:
                export_hashes["index.jsonl"].update(canonical_bytes(row) + b"\n")
            for chunk in _tsv_chunks(rows, tuple(DEFS["summary"]["properties"])):
                export_hashes["index.tsv"].update(chunk)
            for chunk in _tsv_chunks([], tuple(DEFS["construct"]["properties"])):
                export_hashes["constructs.tsv"].update(chunk)
            export_hashes["citation.json"].update(canonical_bytes(release["citation"]) + b"\n")

            def fasta_bytes(sequence: dict[str, Any]) -> bytes:
                return (sequence_fasta(sequence, release_id=rid, provenance_url=release["manifest_url"], records_url=species_manifest["records_download"]["url"]) if fasta_version == 2 else _legacy_fasta(sequence)).encode("utf-8")

            for target_position, summary in enumerate(rows):
                target_path = target_records[summary["record_id"]]
                target = stored[target_path]["data"]
                if target["summary"] != summary:
                    raise PublicDataError("Public target summary differs from index")
                if summary["target_id"] is not None and species_manifest["target_paths"].get(summary["target_id"]) != target_path:
                    raise PublicDataError("Public resolved target mapping differs from its indexed record")
                sections = {}
                for section, descriptor in target["sections"].items():
                    values = []
                    for relative in descriptor["pages"]:
                        value = stored.get(relative)
                        if value is None or release["files"][relative]["schema"] != f"section_{section}" or value["species"] != species or value["data"]["target_id"] != summary["target_id"] or value["data"]["section"] != section or len(value["data"]["items"]) > PAGE_SIZE:
                            raise PublicDataError("Invalid public evidence page reference")
                        if value["data"]["state"] != descriptor["state"] or value["data"]["reason"] != descriptor["reason"]:
                            raise PublicDataError("Public evidence availability mismatch")
                        values.extend(value["data"]["items"])
                    if len(values) != descriptor["total"] or descriptor["state"] != summary["section_availability"][section]:
                        raise PublicDataError("Public evidence counts or state differ from summary")
                    sections[section] = {"state": descriptor["state"], "reason": descriptor["reason"], "items": values}
                constructs = []
                for relative in target["constructs"]["pages"]:
                    value = stored.get(relative)
                    if value is None or release["files"][relative]["schema"] != "construct_page" or value["species"] != species or len(value["data"]) > PAGE_SIZE:
                        raise PublicDataError("Invalid public construct page reference")
                    constructs.extend(value["data"])
                if len(constructs) != target["constructs"]["total"] or constructs != sorted(constructs, key=lambda row: row["construct_id"]):
                    raise PublicDataError("Public construct count or ordering mismatch")
                construct_ids = {row["construct_id"] for row in constructs}
                if len(construct_ids) != len(constructs) or construct_ids != set(target["construct_sequence_paths"]):
                    raise PublicDataError("Public construct and sequence-map references disagree")
                if summary["construct_count"] is not None and summary["construct_count"] != len(constructs):
                    raise PublicDataError("Public summary construct count mismatch")
                sequences = {}
                for construct in constructs:
                    if construct["target_id"] != summary["target_id"]:
                        raise PublicDataError("Public construct belongs to another target")
                    expected_id = _hash([summary["target_id"], construct["category"], construct["classification"], {"name": construct["name"], "pdb_id": construct["pdb_id"]}, {"start": construct["start"], "end": construct["end"], "sequence_basis": construct["sequence_basis"]}, construct["sequence_sha256"]])
                    if construct["construct_id"] != expected_id:
                        raise PublicDataError("Public construct ID mismatch")
                    relative = target["construct_sequence_paths"].get(construct["construct_id"])
                    if relative is not None and relative not in stored:
                        raise PublicDataError("Missing public construct sequence reference")
                    value = stored.get(relative) if relative else None
                    sequence = value["data"] if value else None
                    if sequence is not None and (sequence["target_id"] != summary["target_id"] or sequence["construct_id"] != construct["construct_id"] or sequence["species"] != species or sequence["sequence_sha256"] != construct["sequence_sha256"] or sequence["length"] != construct["length"] or sequence["start"] != construct["start"] or sequence["end"] != construct["end"] or sequence["warnings"] != construct["warnings"] or sequence["evidence"] != construct["evidence"]):
                        raise PublicDataError("Public construct and sequence provenance disagree")
                    if (sequence is not None) != (construct["export_state"] == "available"):
                        raise PublicDataError("Public construct export state mismatch")
                    sequences[construct["construct_id"]] = sequence
                relative = target["reference_sequence_path"]
                if relative is not None and relative not in stored:
                    raise PublicDataError("Missing public reference sequence")
                reference = stored[relative]["data"] if relative else None
                if reference is not None:
                    if reference["target_id"] != summary["target_id"] or reference["construct_id"] is not None or reference["species"] != species or reference["sequence_sha256"] != summary["sequence_sha256"] or reference["length"] != summary["sequence_length"] or reference["start"] != 1 or reference["end"] != reference["length"]:
                        raise PublicDataError("Public reference sequence identity mismatch")
                    export_hashes["reference.fasta"].update(fasta_bytes(reference))
                for sequence in sequences.values():
                    if sequence:
                        if reference is None or sequence["sequence_basis"] != "reference_slice" or sequence["sequence"] != reference["sequence"][sequence["start"] - 1:sequence["end"]]:
                            raise PublicDataError("Public construct does not match its declared reference slice")
                        export_hashes["constructs.fasta"].update(fasta_bytes(sequence))
                if target_position:
                    payload_digest.update(b",")
                _semantic_update(payload_digest, {"summary": summary, "sections": sections, "constructs": constructs, "reference_sequence": reference, "construct_sequences": sequences})
                construct_total += len(constructs)
                for construct in constructs:
                    export_hashes["constructs.jsonl"].update(canonical_bytes(construct) + b"\n")
                for chunk in _tsv_chunks(constructs, tuple(DEFS["construct"]["properties"]), include_header=False):
                    export_hashes["constructs.tsv"].update(chunk)
                for asset in sections["structures"]["items"]:
                    if asset["sha256"] and asset["live_asset_url"]:
                        asset_hashes[urlsplit(asset["live_asset_url"]).path.lstrip("/")] = asset["sha256"]
            payload_digest.update(b"]")
            if construct_total != species_manifest["counts"]["constructs"]:
                raise PublicDataError("Public aggregate construct count mismatch")
            if payload_digest.hexdigest() != release["core"]["payload_sha256"].get(species):
                raise PublicDataError("Public records differ from immutable canonical payload identity")
            downloads_by_name = {row["name"]: row for row in species_manifest["downloads"]}
            if len(downloads_by_name) != len(species_manifest["downloads"]):
                raise PublicDataError("Duplicate public download name")
            for name, checksum in export_hashes.items():
                if name not in downloads_by_name or checksum.hexdigest() != downloads_by_name[name]["sha256"]:
                    raise PublicDataError("Public tabular or sequence export differs from projected records")
        if asset_hashes != release["core"]["asset_sha256"]:
            raise PublicDataError("Public asset references differ from canonical asset identity")
    actual = {str(path.relative_to(root)) for path in root.rglob("*") if path.is_file()}
    if actual != allowed or any(path.is_symlink() for path in root.rglob("*")):
        raise PublicDataError("Unexpected file in public API bundle")
    return manifest
