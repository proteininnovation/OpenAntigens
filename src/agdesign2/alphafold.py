from __future__ import annotations

import gzip
import json
from pathlib import Path

from .config import AnalysisConfig
from .exceptions import ExternalServiceError
from .http import HttpClient, _atomic_write
from .structure_utils import parse_alphafold_pdb, load_pae_matrix, pae_matches_sequence_length
from .af3 import find_local_af3_artifacts


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


class AlphaFoldClient:
    METADATA_HOSTS = (
        "https://alphafold.ebi.ac.uk",
        "https://alphafold.com",
        "https://www.alphafold.com",
    )

    def __init__(self, http: HttpClient, config: AnalysisConfig) -> None:
        self.http = http
        self.config = config

    def ensure_artifacts(
        self,
        accession: str,
        *,
        canonical_sequence: str | None = None,
        canonical_isoform_id: str | None = None,
    ) -> tuple[Path, Path]:
        if not self.config.fetch_alphafold:
            fallback = (
                find_local_af3_artifacts(
                    accession,
                    data_dir=self.config.data_dir,
                    canonical_sequence=canonical_sequence,
                )
                if self.config.enable_local_af3_fallback
                else None
            )
            if fallback is not None:
                return fallback
            raise ExternalServiceError(
                f"No compatible local AlphaFold artifact is available for {accession}; remote retrieval is disabled."
            )
        try:
            return self._ensure_alphafold_db_artifacts(
                accession,
                canonical_sequence=canonical_sequence,
                canonical_isoform_id=canonical_isoform_id,
            )
        except Exception:
            if self.config.enable_local_af3_fallback:
                fallback = find_local_af3_artifacts(
                    accession,
                    data_dir=self.config.data_dir,
                    canonical_sequence=canonical_sequence,
                )
                if fallback is not None:
                    return fallback
            raise

    def _ensure_alphafold_db_artifacts(
        self,
        accession: str,
        *,
        canonical_sequence: str | None = None,
        canonical_isoform_id: str | None = None,
    ) -> tuple[Path, Path]:
        metadata = self._fetch_prediction_metadata(accession, canonical_isoform_id=canonical_isoform_id)
        selected_model = self._select_model(
            metadata,
            accession=accession,
            canonical_sequence=canonical_sequence,
            canonical_isoform_id=canonical_isoform_id,
        )

        pdb_path = self.config.data_dir / "alphafold" / f"{accession}.pdb"
        pae_path = self.config.data_dir / "alphafold" / f"{accession}.pae.json"
        meta_path = self.config.data_dir / "alphafold" / f"{accession}.meta.json"

        if not self._is_cached_model_current(
            meta_path,
            selected_model=selected_model,
            canonical_sequence=canonical_sequence,
        ):
            pdb_path.unlink(missing_ok=True)
            pae_path.unlink(missing_ok=True)
        pdb_urls = self._candidate_pdb_urls(accession, metadata, selected_model)
        pae_urls = self._candidate_pae_urls(accession, metadata, selected_model)
        if selected_model is not None:
            self._download_first_available(pdb_urls, pdb_path)
            self._download_first_available(pae_urls, pae_path)
        elif not (pdb_path.exists() and pae_path.exists()):
            errors = []
            for version in ("v6", "v4"):
                pdb_path.unlink(missing_ok=True)
                pae_path.unlink(missing_ok=True)
                try:
                    self._download_first_available([url for url in pdb_urls if f"model_{version}." in url], pdb_path)
                    self._download_first_available([url for url in pae_urls if f"error_{version}." in url], pae_path)
                    break
                except ExternalServiceError as exc:
                    errors.append(str(exc))
            else:
                pdb_path.unlink(missing_ok=True)
                pae_path.unlink(missing_ok=True)
                raise ExternalServiceError("No complete AlphaFold PDB/PAE version: " + " | ".join(errors))
        self._validate_downloaded_model(pdb_path, canonical_sequence)
        length = len(canonical_sequence) if canonical_sequence else len(parse_alphafold_pdb(pdb_path)["residue_names"])
        if not pae_matches_sequence_length(load_pae_matrix(pae_path), length):
            pae_path.unlink(missing_ok=True)
            meta_path.unlink(missing_ok=True)
            raise ExternalServiceError("AlphaFold PAE dimensions do not match the structure sequence")
        self._write_metadata_sidecar(
            meta_path,
            selected_model=selected_model,
            canonical_sequence=canonical_sequence,
            canonical_isoform_id=canonical_isoform_id,
        )
        return pdb_path, pae_path

    def _fetch_prediction_metadata(self, accession: str, canonical_isoform_id: str | None = None) -> list[dict]:
        keys = [accession]
        isoform = str(canonical_isoform_id or "").strip()
        if isoform and isoform not in keys:
            keys.append(isoform)
        collected: list[dict] = []
        seen: set[str] = set()
        for key in keys:
            for host in self.METADATA_HOSTS:
                url = f"{host}/api/prediction/{key}"
                try:
                    payload = self.http.fetch_json(url, cache_namespace="alphafold_meta")
                except Exception:
                    continue
                items = payload if isinstance(payload, list) else [payload] if isinstance(payload, dict) else []
                for item in items:
                    if not isinstance(item, dict):
                        continue
                    model_id = str(item.get("modelEntityId") or item.get("entryId") or "").strip()
                    fingerprint = model_id or json.dumps(item, sort_keys=True)
                    if fingerprint in seen:
                        continue
                    seen.add(fingerprint)
                    collected.append(item)
        return collected

    def _candidate_pdb_urls(self, accession: str, metadata: list[dict], selected_model: dict | None) -> list[str]:
        urls: list[str] = []
        if selected_model is not None:
            for key in ("pdbUrl", "pdb_url", "modelUrl", "model_url"):
                value = selected_model.get(key)
                if isinstance(value, str) and value.strip():
                    urls.append(value.strip())
        else:
            for model in metadata:
                for key in ("pdbUrl", "pdb_url", "modelUrl", "model_url"):
                    value = model.get(key)
                    if isinstance(value, str) and value.strip():
                        urls.append(value.strip())
        if selected_model is not None:
            return _dedupe(urls)
        for host in self.METADATA_HOSTS:
            urls.append(f"{host}/files/AF-{accession}-F1-model_v6.pdb")
            urls.append(f"{host}/files/AF-{accession}-F1-model_v6.pdb.gz")
            urls.append(f"{host}/files/AF-{accession}-F1-model_v4.pdb")
            urls.append(f"{host}/files/AF-{accession}-F1-model_v4.pdb.gz")
        return _dedupe(urls)

    def _candidate_pae_urls(self, accession: str, metadata: list[dict], selected_model: dict | None) -> list[str]:
        urls: list[str] = []
        if selected_model is not None:
            for key in ("paeDocUrl", "pae_doc_url", "paeUrl", "pae_url", "predictedAlignedErrorUrl"):
                value = selected_model.get(key)
                if isinstance(value, str) and value.strip():
                    urls.append(value.strip())
        else:
            for model in metadata:
                for key in ("paeDocUrl", "pae_doc_url", "paeUrl", "pae_url", "predictedAlignedErrorUrl"):
                    value = model.get(key)
                    if isinstance(value, str) and value.strip():
                        urls.append(value.strip())
        if selected_model is not None:
            return _dedupe(urls)
        for host in self.METADATA_HOSTS:
            urls.append(f"{host}/files/AF-{accession}-F1-predicted_aligned_error_v6.json")
            urls.append(f"{host}/files/AF-{accession}-F1-predicted_aligned_error_v6.json.gz")
            urls.append(f"{host}/files/AF-{accession}-F1-predicted_aligned_error_v4.json")
            urls.append(f"{host}/files/AF-{accession}-F1-predicted_aligned_error_v4.json.gz")
        return _dedupe(urls)

    def _select_model(
        self,
        metadata: list[dict],
        *,
        accession: str,
        canonical_sequence: str | None,
        canonical_isoform_id: str | None,
    ) -> dict | None:
        sequence = str(canonical_sequence or "").strip()
        isoform_id = str(canonical_isoform_id or "").strip()
        if not metadata:
            return None

        def normalized_model_isoform(model: dict) -> str:
            return str(model.get("uniprotAccession") or "").strip()

        if sequence:
            exact_sequence_matches = [model for model in metadata if str(model.get("sequence") or "") == sequence]
            if isoform_id:
                isoform_exact = [model for model in exact_sequence_matches if normalized_model_isoform(model) == isoform_id]
                if isoform_exact:
                    return isoform_exact[0]
            accession_exact = [model for model in exact_sequence_matches if normalized_model_isoform(model) in {accession, isoform_id}]
            if accession_exact:
                return accession_exact[0]
            if exact_sequence_matches:
                return exact_sequence_matches[0]
            raise ExternalServiceError(
                f"No AlphaFold model matches the canonical UniProt sequence for {accession}."
            )

        if isoform_id:
            isoform_matches = [model for model in metadata if normalized_model_isoform(model) == isoform_id]
            if isoform_matches:
                return isoform_matches[0]
        accession_matches = [model for model in metadata if normalized_model_isoform(model) == accession]
        if accession_matches:
            return accession_matches[0]
        return metadata[0]

    def _download_first_available(self, urls: list[str], destination: Path) -> Path:
        if destination.exists():
            return destination
        errors: list[str] = []
        for url in urls:
            try:
                return self._download_url(url, destination)
            except Exception as exc:
                errors.append(f"{url}: {exc}")
        raise ExternalServiceError(
            "All AlphaFold download URLs failed for {destination}: {errors}".format(
                destination=destination.name,
                errors=" | ".join(errors[:6]),
            )
        )

    def _download_url(self, url: str, destination: Path) -> Path:
        if url.endswith(".gz"):
            temp_path = destination.with_name(f"{destination.name}.gz")
            self.http.download(url, temp_path)
            destination.parent.mkdir(parents=True, exist_ok=True)
            with gzip.open(temp_path, "rb") as handle:
                _atomic_write(destination, handle.read())
            temp_path.unlink(missing_ok=True)
            return destination
        return self.http.download(url, destination)

    def _validate_downloaded_model(self, pdb_path: Path, canonical_sequence: str | None) -> None:
        sequence = str(canonical_sequence or "").strip()
        if not sequence:
            return
        structure = parse_alphafold_pdb(pdb_path)
        residue_names = structure.get("residue_names") or {}
        structure_sequence = "".join(
            _THREE_TO_ONE.get(str(residue_names[residue_id]).upper(), "X")
            for residue_id in sorted(residue_names)
        )
        if structure_sequence != sequence:
            pdb_path.unlink(missing_ok=True)
            raise ExternalServiceError(
                "Downloaded AlphaFold structure does not match the canonical UniProt sequence."
            )

    def _is_cached_model_current(
        self,
        meta_path: Path,
        *,
        selected_model: dict | None,
        canonical_sequence: str | None,
    ) -> bool:
        if not meta_path.exists():
            return False
        try:
            cached = json.loads(meta_path.read_text(encoding="utf-8"))
        except Exception:
            return False
        selected_model_id = str((selected_model or {}).get("modelEntityId") or (selected_model or {}).get("entryId") or "")
        if cached.get("artifact_validation_version") != 1:
            return False
        if cached.get("model_entity_id") != selected_model_id:
            return False
        if (canonical_sequence or "") != str(cached.get("canonical_sequence") or ""):
            return False
        return True

    def _write_metadata_sidecar(
        self,
        meta_path: Path,
        *,
        selected_model: dict | None,
        canonical_sequence: str | None,
        canonical_isoform_id: str | None,
    ) -> None:
        payload = {
            "artifact_validation_version": 1,
            "model_entity_id": (selected_model or {}).get("modelEntityId") or (selected_model or {}).get("entryId"),
            "model_uniprot_accession": (selected_model or {}).get("uniprotAccession"),
            "canonical_isoform_id": canonical_isoform_id,
            "canonical_sequence": canonical_sequence,
        }
        meta_path.parent.mkdir(parents=True, exist_ok=True)
        meta_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def _dedupe(urls: list[str]) -> list[str]:
    seen: set[str] = set()
    ordered: list[str] = []
    for url in urls:
        cleaned = str(url or "").strip()
        if not cleaned or cleaned in seen:
            continue
        seen.add(cleaned)
        ordered.append(cleaned)
    return ordered
