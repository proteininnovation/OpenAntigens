from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
import shutil
import sqlite3
import zipfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from Bio.PDB import MMCIFParser, PDBIO

from .structure_utils import load_pae_matrix, parse_alphafold_pdb


_THREE_TO_ONE = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "ASX": "B",
    "GLX": "Z",
    "SEC": "U",
    "PYL": "O",
    "UNK": "X",
}


@dataclass(slots=True)
class AF3PrepareResult:
    manifest_tsv: Path
    manifest_json: Path
    output_dir: Path
    rows: list[dict[str, Any]]


@dataclass(slots=True)
class AF3CatalogueResult:
    db_path: Path
    manifest_tsv: Path
    manifest_json: Path
    output_dir: Path
    rows: list[dict[str, Any]]


def prepare_af3_structures(
    zip_path: str | Path,
    *,
    snapshot_dir: str | Path,
    output_dir: str | Path | None = None,
    verbose: bool = False,
) -> AF3PrepareResult:
    archive = Path(zip_path)
    snapshot = Path(snapshot_dir)
    target_output = Path(output_dir) if output_dir is not None else snapshot / "data" / "alphafold3"
    target_output.mkdir(parents=True, exist_ok=True)
    reports = _index_snapshot_reports(snapshot)
    by_accession = {str(row.get("accession") or "").upper(): row for row in reports if row.get("accession")}
    by_sequence = {str(row.get("sequence") or ""): row for row in reports if row.get("sequence")}
    rows: list[dict[str, Any]] = []

    with zipfile.ZipFile(archive) as zf:
        names = zf.namelist()
        target_dirs = sorted({name.split("/", 1)[0] for name in names if name.startswith("openantigen_af3_") and "/" in name})
        terms_present = "terms_of_use.md" in names
        for index, target_dir in enumerate(target_dirs, start=1):
            row = _prepare_one_af3_target(
                zf,
                archive=archive,
                target_dir=target_dir,
                output_dir=target_output,
                reports_by_accession=by_accession,
                reports_by_sequence=by_sequence,
                terms_present=terms_present,
            )
            rows.append(row)
            if verbose:
                print(
                    f"[agdesign2-af3] {index}/{len(target_dirs)} {target_dir}: {row.get('status')} {row.get('entry_name') or row.get('accession') or ''}",
                    flush=True,
                )

    manifest_tsv = target_output / "af3_manifest.tsv"
    manifest_json = target_output / "af3_manifest.json"
    _write_manifest_tsv(manifest_tsv, rows)
    manifest_json.write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
    return AF3PrepareResult(manifest_tsv=manifest_tsv, manifest_json=manifest_json, output_dir=target_output, rows=rows)


def catalogue_af3_structures(
    input_path: str | Path,
    *,
    output_dir: str | Path = "data/af3_catalog",
    verbose: bool = False,
) -> AF3CatalogueResult:
    """Catalogue AlphaFold Server exports without expanding large raw payloads."""

    source = Path(input_path)
    archives = [source] if source.is_file() else sorted(source.glob("*.zip"))
    if not archives:
        raise FileNotFoundError(f"No AF3 ZIP files found: {source}")

    target_output = Path(output_dir)
    artifacts_dir = target_output / "alphafold3"
    artifacts_dir.mkdir(parents=True, exist_ok=True)
    db_path = target_output / "af3_catalog.sqlite"
    manifest_tsv = target_output / "manifest.tsv"
    manifest_json = target_output / "manifest.json"
    terms_dir = target_output / "terms_of_use"

    archive_rows: list[dict[str, Any]] = []
    member_rows: list[dict[str, Any]] = []
    model_rows: list[dict[str, Any]] = []
    rows: list[dict[str, Any]] = []
    for archive_index, archive in enumerate(archives, start=1):
        archive_result = _scan_af3_archive(archive)
        archive_rows.append(archive_result["archive"])
        member_rows.extend(archive_result["members"])
        model_rows.extend(archive_result["models"])
        rows.extend(archive_result["targets"])
        if verbose:
            print(
                f"[agdesign2-af3-catalog] {archive_index}/{len(archives)} {archive}: "
                f"{len(archive_result['targets'])} targets",
                flush=True,
            )

    winners = _select_catalogue_winners(rows)
    for row in rows:
        key = _catalogue_target_key(row)
        row["selected_for_accession"] = bool(row.get("status") == "ok" and winners.get(key) is row)
    for row in winners.values():
        _materialize_catalogue_artifact(row, artifacts_dir)

    _write_catalogue_sqlite(db_path, archive_rows=archive_rows, member_rows=member_rows, target_rows=rows, model_rows=model_rows)
    _write_catalogue_manifest_tsv(manifest_tsv, rows)
    manifest_json.write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
    _write_catalogue_terms_of_use(archives, terms_dir)
    return AF3CatalogueResult(
        db_path=db_path,
        manifest_tsv=manifest_tsv,
        manifest_json=manifest_json,
        output_dir=target_output,
        rows=rows,
    )


def copy_af3_catalog_artifacts(
    catalog_db: str | Path,
    *,
    data_dir: str | Path,
    target_tsv: str | Path | None = None,
) -> dict[str, Any]:
    database = Path(catalog_db)
    if not database.exists():
        raise FileNotFoundError(f"AF3 catalogue database not found: {database}")
    accessions = _target_tsv_accessions(Path(target_tsv)) if target_tsv is not None else None
    target_dir = Path(data_dir) / "alphafold3"
    target_dir.mkdir(parents=True, exist_ok=True)
    copied: list[str] = []
    skipped: list[str] = []
    with sqlite3.connect(database) as conn:
        conn.row_factory = sqlite3.Row
        rows = conn.execute(
            """
            SELECT accession, sequence_sha256, pdb_path, pae_path, meta_path
            FROM targets
            WHERE status = 'ok' AND selected_for_accession = 1
            ORDER BY accession
            """
        ).fetchall()
    for item in rows:
        accession = str(item["accession"] or "").strip()
        if not accession:
            skipped.append(accession)
            continue
        if accessions is not None and accession.upper() not in accessions:
            continue
        source_paths = [_catalogue_artifact_path(database.parent, str(item[key])) for key in ("pdb_path", "pae_path", "meta_path")]
        if not all(path.exists() for path in source_paths):
            raise FileNotFoundError(f"AF3 catalogue artifact triplet is incomplete for {accession}")
        for path in source_paths:
            shutil.copy2(path, target_dir / path.name)
        copied.append(accession)
    manifest = {
        "source_catalog": str(database),
        "data_dir": str(Path(data_dir)),
        "target_tsv": str(target_tsv) if target_tsv is not None else None,
        "copied_count": len(copied),
        "copied_accessions": copied,
        "generated_at_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
    }
    (target_dir / "af3_catalog_copy_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def _catalogue_artifact_path(catalogue_dir: Path, path_text: str) -> Path:
    path = Path(path_text)
    return path if path.is_absolute() else catalogue_dir / path


def find_local_af3_artifacts(
    accession: str,
    *,
    data_dir: Path,
    canonical_sequence: str | None = None,
) -> tuple[Path, Path] | None:
    accession_text = str(accession or "").strip()
    if not accession_text:
        return None
    directory = data_dir / "alphafold3"
    pdb_path = directory / f"{accession_text}.pdb"
    pae_path = directory / f"{accession_text}.pae.json"
    if not pae_path.exists():
        pae_path = directory / f"{accession_text}.pae.json.gz"
    meta_path = directory / f"{accession_text}.meta.json"
    if not pdb_path.exists() or not pae_path.exists() or not meta_path.exists():
        return None
    try:
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
    except Exception:
        return None
    if meta.get("source") != "alphafold_server" or meta.get("model_type") != "AF3":
        return None
    sequence = str(canonical_sequence or "").strip()
    if sequence:
        expected = hashlib.sha256(sequence.encode("utf-8")).hexdigest()
        if str(meta.get("sequence_sha256") or "") != expected:
            return None
        if _pdb_sequence(pdb_path) != sequence:
            return None
    try:
        load_pae_matrix(pae_path)
    except Exception:
        return None
    return pdb_path, pae_path


def is_af3_artifact(path: str | Path | None) -> bool:
    if path is None:
        return False
    return "alphafold3" in Path(path).parts


def _prepare_one_af3_target(
    zf: zipfile.ZipFile,
    *,
    archive: Path,
    target_dir: str,
    output_dir: Path,
    reports_by_accession: dict[str, dict[str, Any]],
    reports_by_sequence: dict[str, dict[str, Any]],
    terms_present: bool,
) -> dict[str, Any]:
    row: dict[str, Any] = {
        "target_dir": target_dir,
        "source_zip": archive.name,
        "status": "error",
        "error": "",
        "terms_present": terms_present,
        "distribution_policy": "internal_fallback",
    }
    try:
        accession_hint = _accession_from_target_dir(target_dir)
        job_request_name = _single_name(zf, target_dir, "_job_request.json")
        job_request = json.loads(zf.read(job_request_name))
        sequence = _sequence_from_job_request(job_request)
        row.update(
            {
                "accession_hint": accession_hint,
                "input_sequence_length": len(sequence),
                "input_sequence_sha256": hashlib.sha256(sequence.encode("utf-8")).hexdigest() if sequence else "",
            }
        )
        if not sequence:
            row["error"] = "No single protein sequence found in AF3 job request."
            return row
        report = reports_by_accession.get(accession_hint.upper()) if accession_hint else None
        if report is None:
            report = reports_by_sequence.get(sequence)
        if report is None:
            row["error"] = "Could not correlate AF3 job to a snapshot report by accession or exact sequence."
            return row
        accession = str(report.get("accession") or accession_hint or "").strip()
        entry_name = str(report.get("entry_name") or "").strip()
        canonical_sequence = str(report.get("sequence") or "")
        row.update({"accession": accession, "entry_name": entry_name, "gene_symbol": report.get("gene_symbol") or ""})
        if not accession:
            row["error"] = "Resolved report does not contain a UniProt accession."
            return row
        if canonical_sequence != sequence:
            row["error"] = "AF3 job sequence does not match the canonical snapshot report sequence."
            return row
        model = _select_best_model(zf, target_dir)
        if model is None:
            row["error"] = "No complete AF3 model/confidence/full-data triplet found."
            return row
        model_index, cif_name, summary_name, full_data_name, summary, full_data = model
        pdb_path = output_dir / f"{accession}.pdb"
        pae_path = output_dir / f"{accession}.pae.json"
        meta_path = output_dir / f"{accession}.meta.json"
        temp_cif = output_dir / f".{accession}.model_{model_index}.cif"
        try:
            temp_cif.write_bytes(zf.read(cif_name))
            _convert_cif_to_pdb(temp_cif, pdb_path)
        finally:
            temp_cif.unlink(missing_ok=True)
        structure_sequence = _pdb_sequence(pdb_path)
        if structure_sequence != canonical_sequence:
            pdb_path.unlink(missing_ok=True)
            row["error"] = "Converted AF3 PDB sequence does not match the canonical snapshot report sequence."
            return row
        pae = full_data.get("pae")
        if not isinstance(pae, list) or len(pae) != len(canonical_sequence):
            pdb_path.unlink(missing_ok=True)
            row["error"] = "AF3 full-data PAE matrix is missing or does not match sequence length."
            return row
        pae_path.write_text(json.dumps({"pae": pae}, separators=(",", ":")) + "\n", encoding="utf-8")
        meta = {
            "source": "alphafold_server",
            "model_type": "AF3",
            "distribution_policy": "internal_fallback",
            "source_zip": str(archive),
            "target_dir": target_dir,
            "model_index": model_index,
            "model_cif": cif_name,
            "summary_confidences": summary_name,
            "full_data": full_data_name,
            "accession": accession,
            "entry_name": entry_name,
            "sequence_sha256": hashlib.sha256(canonical_sequence.encode("utf-8")).hexdigest(),
            "sequence_length": len(canonical_sequence),
            "ranking_score": summary.get("ranking_score"),
            "ptm": summary.get("ptm"),
            "fraction_disordered": summary.get("fraction_disordered"),
            "terms_present": terms_present,
        }
        meta_path.write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
        row.update(
            {
                "status": "ok",
                "selected_model_index": model_index,
                "ranking_score": summary.get("ranking_score"),
                "ptm": summary.get("ptm"),
                "fraction_disordered": summary.get("fraction_disordered"),
                "pdb_path": str(pdb_path),
                "pae_path": str(pae_path),
                "meta_path": str(meta_path),
            }
        )
        return row
    except Exception as exc:
        row["error"] = str(exc)
        return row


def _scan_af3_archive(archive: Path) -> dict[str, list[dict[str, Any]] | dict[str, Any]]:
    archive_row = {
        "path": str(archive),
        "name": archive.name,
        "size_bytes": archive.stat().st_size,
        "sha256": _sha256_file(archive),
        "terms_present": False,
        "target_count": 0,
    }
    members: list[dict[str, Any]] = []
    targets: list[dict[str, Any]] = []
    models: list[dict[str, Any]] = []
    with zipfile.ZipFile(archive) as zf:
        names = zf.namelist()
        name_set = set(names)
        target_dirs = sorted({name.split("/", 1)[0] for name in names if name.startswith("openantigen_af3_") and "/" in name})
        archive_row["terms_present"] = "terms_of_use.md" in name_set
        archive_row["target_count"] = len(target_dirs)
        for info in zf.infolist():
            members.append(
                {
                    "archive_path": str(archive),
                    "member_name": info.filename,
                    "file_size": info.file_size,
                    "compress_size": info.compress_size,
                    "crc": f"{info.CRC:08x}",
                    "modified": _zip_modified_iso(info),
                }
            )
        for target_dir in target_dirs:
            target_row, target_models = _scan_af3_target(
                zf,
                archive=archive,
                target_dir=target_dir,
                names=name_set,
                terms_present=bool(archive_row["terms_present"]),
            )
            targets.append(target_row)
            models.extend(target_models)
    return {"archive": archive_row, "members": members, "targets": targets, "models": models}


def _scan_af3_target(
    zf: zipfile.ZipFile,
    *,
    archive: Path,
    target_dir: str,
    names: set[str],
    terms_present: bool,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    row: dict[str, Any] = {
        "target_dir": target_dir,
        "source_zip": str(archive),
        "status": "error",
        "error": "",
        "terms_present": terms_present,
        "distribution_policy": "public_with_notice",
        "selected_for_accession": False,
    }
    model_rows: list[dict[str, Any]] = []
    try:
        accession_hint = _accession_from_target_dir(target_dir)
        job_request_name = _single_name(zf, target_dir, "_job_request.json")
        job_request = json.loads(zf.read(job_request_name))
        sequence = _sequence_from_job_request(job_request)
        sequence_sha256 = hashlib.sha256(sequence.encode("utf-8")).hexdigest() if sequence else ""
        row.update(
            {
                "accession": accession_hint,
                "accession_hint": accession_hint,
                "job_request": job_request_name,
                "input_sequence_length": len(sequence),
                "input_sequence_sha256": sequence_sha256,
                "sequence_sha256": sequence_sha256,
            }
        )
        if not accession_hint:
            row["error"] = "Could not parse UniProt accession from AF3 target directory."
            return row, model_rows
        if not sequence:
            row["error"] = "No single protein sequence found in AF3 job request."
            return row, model_rows
        models = _scan_af3_models(zf, target_dir, names=names)
        for model_row in models:
            model_row["source_zip"] = str(archive)
        model_rows.extend(models)
        complete_models = [item for item in models if item.get("has_complete_triplet")]
        if not complete_models:
            row["error"] = "No complete AF3 model/confidence/full-data triplet found."
            return row, model_rows
        winner = max(
            complete_models,
            key=lambda item: (
                float(item.get("ranking_score") or 0.0),
                float(item.get("ptm") or 0.0),
                -int(item.get("model_index") or 0),
            ),
        )
        row.update(
            {
                "status": "ok",
                "selected_model_index": winner.get("model_index"),
                "ranking_score": winner.get("ranking_score"),
                "ptm": winner.get("ptm"),
                "fraction_disordered": winner.get("fraction_disordered"),
                "model_cif": winner.get("model_cif"),
                "summary_confidences": winner.get("summary_confidences"),
                "full_data": winner.get("full_data"),
            }
        )
        return row, model_rows
    except Exception as exc:
        row["error"] = str(exc)
        return row, model_rows


def _scan_af3_models(zf: zipfile.ZipFile, target_dir: str, *, names: set[str]) -> list[dict[str, Any]]:
    model_rows: list[dict[str, Any]] = []
    pattern = re.compile(r"_model_(\d+)\.cif$")
    for cif_name in zf.namelist():
        if not cif_name.startswith(f"{target_dir}/"):
            continue
        match = pattern.search(cif_name)
        if not match:
            continue
        index = int(match.group(1))
        summary_name = cif_name.replace(f"_model_{index}.cif", f"_summary_confidences_{index}.json")
        full_data_name = cif_name.replace(f"_model_{index}.cif", f"_full_data_{index}.json")
        summary: dict[str, Any] = {}
        if summary_name in names:
            summary = json.loads(zf.read(summary_name))
        model_rows.append(
            {
                "target_dir": target_dir,
                "model_index": index,
                "model_cif": cif_name,
                "summary_confidences": summary_name,
                "full_data": full_data_name,
                "has_complete_triplet": summary_name in names and full_data_name in names,
                "ranking_score": summary.get("ranking_score"),
                "ptm": summary.get("ptm"),
                "fraction_disordered": summary.get("fraction_disordered"),
            }
        )
    return model_rows


def _select_catalogue_winners(rows: list[dict[str, Any]]) -> dict[tuple[str, str], dict[str, Any]]:
    sequence_by_accession: dict[str, str] = {}
    winners: dict[tuple[str, str], dict[str, Any]] = {}
    for row in rows:
        if row.get("status") != "ok":
            continue
        accession = str(row.get("accession") or "").upper()
        sequence_sha256 = str(row.get("sequence_sha256") or "")
        if not accession or not sequence_sha256:
            continue
        existing_sequence = sequence_by_accession.get(accession)
        if existing_sequence is not None and existing_sequence != sequence_sha256:
            raise ValueError(f"Conflicting AF3 sequence hashes for {accession}")
        sequence_by_accession[accession] = sequence_sha256
        key = (accession, sequence_sha256)
        existing = winners.get(key)
        if existing is None or _catalogue_rank(row) > _catalogue_rank(existing):
            winners[key] = row
    return winners


def _catalogue_target_key(row: dict[str, Any]) -> tuple[str, str]:
    return (str(row.get("accession") or "").upper(), str(row.get("sequence_sha256") or ""))


def _catalogue_rank(row: dict[str, Any]) -> tuple[float, float, int, str]:
    return (
        float(row.get("ranking_score") or 0.0),
        float(row.get("ptm") or 0.0),
        -int(row.get("selected_model_index") or 0),
        str(row.get("source_zip") or ""),
    )


def _materialize_catalogue_artifact(row: dict[str, Any], artifacts_dir: Path) -> None:
    archive = Path(str(row.get("source_zip") or ""))
    accession = str(row.get("accession") or "").strip()
    sequence_sha256 = str(row.get("sequence_sha256") or "")
    sequence_length = int(row.get("input_sequence_length") or 0)
    target_dir = str(row.get("target_dir") or "")
    if not archive.exists() or not accession or not target_dir:
        raise ValueError(f"Cannot materialize incomplete AF3 catalogue row: {row}")
    artifacts_dir.mkdir(parents=True, exist_ok=True)
    pdb_path = artifacts_dir / f"{accession}.pdb"
    pae_path = artifacts_dir / f"{accession}.pae.json.gz"
    meta_path = artifacts_dir / f"{accession}.meta.json"
    temp_cif = artifacts_dir / f".{accession}.model_{row.get('selected_model_index')}.cif"
    with zipfile.ZipFile(archive) as zf:
        try:
            temp_cif.write_bytes(zf.read(str(row.get("model_cif") or "")))
            _convert_cif_to_pdb(temp_cif, pdb_path)
        finally:
            temp_cif.unlink(missing_ok=True)
        structure_sequence = _pdb_sequence(pdb_path)
        if hashlib.sha256(structure_sequence.encode("utf-8")).hexdigest() != sequence_sha256:
            pdb_path.unlink(missing_ok=True)
            raise ValueError(f"Converted AF3 PDB sequence does not match the AF3 job sequence for {accession}")
        full_data = json.loads(zf.read(str(row.get("full_data") or "")))
    pae = full_data.get("pae")
    if not isinstance(pae, list) or len(pae) != sequence_length:
        pdb_path.unlink(missing_ok=True)
        raise ValueError(f"AF3 full-data PAE matrix is missing or does not match sequence length for {accession}")
    _write_json_gz(pae_path, {"pae": pae})
    meta = {
        "source": "alphafold_server",
        "model_type": "AF3",
        "distribution_policy": row.get("distribution_policy") or "public_with_notice",
        "source_zip": str(archive),
        "target_dir": target_dir,
        "model_index": row.get("selected_model_index"),
        "model_cif": row.get("model_cif"),
        "summary_confidences": row.get("summary_confidences"),
        "full_data": row.get("full_data"),
        "accession": accession,
        "sequence_sha256": sequence_sha256,
        "sequence_length": sequence_length,
        "ranking_score": row.get("ranking_score"),
        "ptm": row.get("ptm"),
        "fraction_disordered": row.get("fraction_disordered"),
        "terms_present": row.get("terms_present"),
    }
    meta_path.write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    catalogue_root = artifacts_dir.parent
    row.update(
        {
            "pdb_path": str(pdb_path.relative_to(catalogue_root)),
            "pae_path": str(pae_path.relative_to(catalogue_root)),
            "meta_path": str(meta_path.relative_to(catalogue_root)),
        }
    )


def _write_json_gz(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8", compresslevel=9) as handle:
        json.dump(payload, handle, separators=(",", ":"))
        handle.write("\n")


def _write_catalogue_sqlite(
    path: Path,
    *,
    archive_rows: list[dict[str, Any]],
    member_rows: list[dict[str, Any]],
    target_rows: list[dict[str, Any]],
    model_rows: list[dict[str, Any]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with sqlite3.connect(path) as conn:
        conn.executescript(
            """
            DROP TABLE IF EXISTS metadata;
            DROP TABLE IF EXISTS archives;
            DROP TABLE IF EXISTS archive_members;
            DROP TABLE IF EXISTS targets;
            DROP TABLE IF EXISTS models;
            CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL);
            CREATE TABLE archives (
                path TEXT PRIMARY KEY,
                name TEXT NOT NULL,
                size_bytes INTEGER NOT NULL,
                sha256 TEXT NOT NULL,
                terms_present INTEGER NOT NULL,
                target_count INTEGER NOT NULL
            );
            CREATE TABLE archive_members (
                archive_path TEXT NOT NULL,
                member_name TEXT NOT NULL,
                file_size INTEGER NOT NULL,
                compress_size INTEGER NOT NULL,
                crc TEXT NOT NULL,
                modified TEXT NOT NULL,
                PRIMARY KEY (archive_path, member_name)
            );
            CREATE TABLE targets (
                source_zip TEXT NOT NULL,
                target_dir TEXT NOT NULL,
                accession TEXT NOT NULL,
                accession_hint TEXT NOT NULL,
                status TEXT NOT NULL,
                error TEXT NOT NULL,
                sequence_length INTEGER NOT NULL,
                sequence_sha256 TEXT NOT NULL,
                selected_model_index INTEGER,
                ranking_score REAL,
                ptm REAL,
                fraction_disordered REAL,
                distribution_policy TEXT NOT NULL,
                terms_present INTEGER NOT NULL,
                selected_for_accession INTEGER NOT NULL,
                pdb_path TEXT NOT NULL,
                pae_path TEXT NOT NULL,
                meta_path TEXT NOT NULL,
                PRIMARY KEY (source_zip, target_dir)
            );
            CREATE TABLE models (
                source_zip TEXT NOT NULL,
                target_dir TEXT NOT NULL,
                model_index INTEGER NOT NULL,
                model_cif TEXT NOT NULL,
                summary_confidences TEXT NOT NULL,
                full_data TEXT NOT NULL,
                has_complete_triplet INTEGER NOT NULL,
                ranking_score REAL,
                ptm REAL,
                fraction_disordered REAL,
                PRIMARY KEY (source_zip, target_dir, model_index, model_cif)
            );
            CREATE INDEX idx_af3_targets_accession ON targets(accession);
            CREATE INDEX idx_af3_targets_sequence ON targets(sequence_sha256);
            """
        )
        conn.executemany(
            "INSERT INTO metadata(key, value) VALUES (?, ?)",
            [
                ("schema_version", "1"),
                ("generated_at_utc", datetime.now(timezone.utc).replace(microsecond=0).isoformat()),
                ("target_count", str(len(target_rows))),
            ],
        )
        conn.executemany(
            """
            INSERT INTO archives(path, name, size_bytes, sha256, terms_present, target_count)
            VALUES (:path, :name, :size_bytes, :sha256, :terms_present, :target_count)
            """,
            [{**row, "terms_present": int(bool(row.get("terms_present")))} for row in archive_rows],
        )
        conn.executemany(
            """
            INSERT INTO archive_members(archive_path, member_name, file_size, compress_size, crc, modified)
            VALUES (:archive_path, :member_name, :file_size, :compress_size, :crc, :modified)
            """,
            member_rows,
        )
        conn.executemany(
            """
            INSERT INTO targets(
                source_zip, target_dir, accession, accession_hint, status, error, sequence_length, sequence_sha256,
                selected_model_index, ranking_score, ptm, fraction_disordered, distribution_policy, terms_present,
                selected_for_accession, pdb_path, pae_path, meta_path
            )
            VALUES (
                :source_zip, :target_dir, :accession, :accession_hint, :status, :error, :sequence_length, :sequence_sha256,
                :selected_model_index, :ranking_score, :ptm, :fraction_disordered, :distribution_policy, :terms_present,
                :selected_for_accession, :pdb_path, :pae_path, :meta_path
            )
            """,
            [_catalogue_sql_target_row(row) for row in target_rows],
        )
        conn.executemany(
            """
            INSERT INTO models(
                source_zip, target_dir, model_index, model_cif, summary_confidences, full_data, has_complete_triplet,
                ranking_score, ptm, fraction_disordered
            )
            VALUES (
                :source_zip, :target_dir, :model_index, :model_cif, :summary_confidences, :full_data, :has_complete_triplet,
                :ranking_score, :ptm, :fraction_disordered
            )
            """,
            [{**row, "has_complete_triplet": int(bool(row.get("has_complete_triplet")))} for row in model_rows],
        )


def _catalogue_sql_target_row(row: dict[str, Any]) -> dict[str, Any]:
    return {
        "source_zip": str(row.get("source_zip") or ""),
        "target_dir": str(row.get("target_dir") or ""),
        "accession": str(row.get("accession") or ""),
        "accession_hint": str(row.get("accession_hint") or ""),
        "status": str(row.get("status") or ""),
        "error": str(row.get("error") or ""),
        "sequence_length": int(row.get("input_sequence_length") or 0),
        "sequence_sha256": str(row.get("sequence_sha256") or row.get("input_sequence_sha256") or ""),
        "selected_model_index": row.get("selected_model_index"),
        "ranking_score": row.get("ranking_score"),
        "ptm": row.get("ptm"),
        "fraction_disordered": row.get("fraction_disordered"),
        "distribution_policy": str(row.get("distribution_policy") or ""),
        "terms_present": int(bool(row.get("terms_present"))),
        "selected_for_accession": int(bool(row.get("selected_for_accession"))),
        "pdb_path": str(row.get("pdb_path") or ""),
        "pae_path": str(row.get("pae_path") or ""),
        "meta_path": str(row.get("meta_path") or ""),
    }


def _write_catalogue_manifest_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "target_dir",
        "status",
        "error",
        "accession",
        "accession_hint",
        "input_sequence_length",
        "sequence_sha256",
        "selected_model_index",
        "ranking_score",
        "ptm",
        "fraction_disordered",
        "distribution_policy",
        "selected_for_accession",
        "pdb_path",
        "pae_path",
        "meta_path",
        "source_zip",
        "terms_present",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def _write_catalogue_terms_of_use(archives: list[Path], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest: list[dict[str, Any]] = []
    seen_hashes: dict[str, str] = {}
    for archive in archives:
        with zipfile.ZipFile(archive) as zf:
            if "terms_of_use.md" not in zf.namelist():
                manifest.append({"source_zip": str(archive), "status": "missing", "terms_file": "", "sha256": ""})
                continue
            content = zf.read("terms_of_use.md")
        digest = hashlib.sha256(content).hexdigest()
        destination_name = seen_hashes.get(digest)
        if destination_name is None:
            destination_name = "terms_of_use.md" if not seen_hashes else f"terms_of_use_{len(seen_hashes) + 1}.md"
            (output_dir / destination_name).write_bytes(content)
            seen_hashes[digest] = destination_name
        manifest.append(
            {
                "source_zip": str(archive),
                "status": "ok",
                "terms_file": destination_name,
                "sha256": digest,
            }
        )
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def _target_tsv_accessions(path: Path) -> set[str]:
    if not path.exists():
        raise FileNotFoundError(f"Target TSV not found: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        accessions = {
            str(row.get("uniprot_id") or row.get("accession") or "").strip().upper()
            for row in reader
        }
    return {accession for accession in accessions if accession}


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _zip_modified_iso(info: zipfile.ZipInfo) -> str:
    year, month, day, hour, minute, second = info.date_time
    return f"{year:04d}-{month:02d}-{day:02d}T{hour:02d}:{minute:02d}:{second:02d}"


def _index_snapshot_reports(snapshot: Path) -> list[dict[str, Any]]:
    candidates = [
        snapshot / "data" / "batch_summary.json",
        snapshot / "outputs" / "surfy_batch" / "batch_summary.json",
    ]
    summary_path = next((path for path in candidates if path.exists()), None)
    reports: list[dict[str, Any]] = []
    if summary_path is None:
        return reports
    try:
        rows = json.loads(summary_path.read_text(encoding="utf-8"))
    except Exception:
        return reports
    for row in rows:
        report_path = Path(str(row.get("json_report") or ""))
        if not report_path.is_absolute():
            report_path = summary_path.parent / report_path
        if not report_path.exists():
            fallback = summary_path.parent / f"{str(row.get('entry_name') or row.get('query') or '').lower()}_report.json"
            report_path = fallback if fallback.exists() else report_path
        if not report_path.exists():
            continue
        try:
            payload = json.loads(report_path.read_text(encoding="utf-8"))
        except Exception:
            continue
        target = payload.get("target") or {}
        reports.append(
            {
                "accession": str(target.get("accession") or "").strip(),
                "entry_name": str(target.get("entry_name") or "").strip(),
                "gene_symbol": str(target.get("gene_symbol") or "").strip(),
                "sequence": str(target.get("sequence") or "").strip(),
                "report_path": str(report_path),
            }
        )
    return reports


def _single_name(zf: zipfile.ZipFile, target_dir: str, suffix: str) -> str:
    matches = [name for name in zf.namelist() if name.startswith(f"{target_dir}/") and name.endswith(suffix)]
    if len(matches) != 1:
        raise ValueError(f"Expected one {suffix} in {target_dir}, found {len(matches)}.")
    return matches[0]


def _accession_from_target_dir(target_dir: str) -> str:
    match = re.search(r"_([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9])$", target_dir.upper())
    return match.group(1) if match else ""


def _sequence_from_job_request(job_request: Any) -> str:
    jobs = job_request if isinstance(job_request, list) else [job_request]
    sequences: list[str] = []
    for job in jobs:
        for item in (job or {}).get("sequences") or []:
            protein = item.get("proteinChain") if isinstance(item, dict) else None
            if protein and isinstance(protein.get("sequence"), str):
                count = int(protein.get("count") or 1)
                if count != 1:
                    return ""
                sequences.append(protein["sequence"].strip())
    return sequences[0] if len(sequences) == 1 else ""


def _select_best_model(zf: zipfile.ZipFile, target_dir: str) -> tuple[int, str, str, str, dict[str, Any], dict[str, Any]] | None:
    candidates: list[tuple[float, float, int, str, str, str, dict[str, Any], dict[str, Any]]] = []
    pattern = re.compile(r"_model_(\d+)\.cif$")
    for cif_name in zf.namelist():
        if not cif_name.startswith(f"{target_dir}/"):
            continue
        match = pattern.search(cif_name)
        if not match:
            continue
        index = int(match.group(1))
        summary_name = cif_name.replace(f"_model_{index}.cif", f"_summary_confidences_{index}.json")
        full_data_name = cif_name.replace(f"_model_{index}.cif", f"_full_data_{index}.json")
        if summary_name not in zf.namelist() or full_data_name not in zf.namelist():
            continue
        summary = json.loads(zf.read(summary_name))
        full_data = json.loads(zf.read(full_data_name))
        ranking_score = float(summary.get("ranking_score") or 0.0)
        ptm = float(summary.get("ptm") or 0.0)
        candidates.append((ranking_score, ptm, -index, cif_name, summary_name, full_data_name, summary, full_data))
    if not candidates:
        return None
    ranking_score, ptm, negative_index, cif_name, summary_name, full_data_name, summary, full_data = max(candidates)
    return -negative_index, cif_name, summary_name, full_data_name, summary, full_data


def _convert_cif_to_pdb(cif_path: Path, pdb_path: Path) -> None:
    structure = MMCIFParser(QUIET=True).get_structure(cif_path.stem, str(cif_path))
    pdb_path.parent.mkdir(parents=True, exist_ok=True)
    writer = PDBIO()
    writer.set_structure(structure)
    writer.save(str(pdb_path))


def _pdb_sequence(pdb_path: Path) -> str:
    structure = parse_alphafold_pdb(pdb_path)
    residue_names = structure.get("residue_names") or {}
    return "".join(_THREE_TO_ONE.get(str(residue_names[residue_id]).upper(), "X") for residue_id in sorted(residue_names))


def _write_manifest_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "target_dir",
        "status",
        "error",
        "entry_name",
        "gene_symbol",
        "accession",
        "accession_hint",
        "input_sequence_length",
        "selected_model_index",
        "ranking_score",
        "ptm",
        "fraction_disordered",
        "distribution_policy",
        "pdb_path",
        "pae_path",
        "meta_path",
        "source_zip",
        "terms_present",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})
