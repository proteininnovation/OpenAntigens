from __future__ import annotations

import hashlib
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any

from .config import AnalysisConfig
from .exceptions import BlastDatabaseError
from .http import HttpClient
from .models import BlastHit


class BlastClient:
    def __init__(self, http: HttpClient, config: AnalysisConfig) -> None:
        self.http = http
        self.config = config

    def ensure_databases(self) -> dict[str, Path]:
        self.config.ensure_directories()
        databases: dict[str, Path] = {}
        for species in self.config.blast_species:
            taxon_id = self.config.species_taxonomy[species]
            fasta_path = self.config.data_dir / "swissprot" / f"{species}.fasta"
            db_prefix = self.config.blast_db_dir / species
            if not fasta_path.exists():
                url = (
                    "https://rest.uniprot.org/uniprotkb/stream?"
                    f"format=fasta&query=%28reviewed%3Atrue%20AND%20organism_id%3A{taxon_id}%29"
                )
                self.http.download(url, fasta_path)
            if not (db_prefix.parent / f"{db_prefix.name}.pin").exists():
                blast_fasta = self._cli_compatible_path(fasta_path)
                blast_prefix = self._cli_compatible_path(db_prefix)
                self._run(
                    [
                        self._blast_tool("makeblastdb"),
                        "-in",
                        str(blast_fasta),
                        "-dbtype",
                        "prot",
                        "-out",
                        str(blast_prefix),
                    ]
                )
            databases[species] = db_prefix
        return databases

    def ensure_ortholog_databases(self) -> dict[str, Path]:
        self.config.ensure_directories()
        databases: dict[str, Path] = {}
        for species in self.config.ortholog_blast_species:
            taxon_id = self.config.species_taxonomy[species]
            fasta_path = self.config.ortholog_fasta_dir / f"{species}.fasta"
            db_prefix = self.config.ortholog_blast_db_dir / species
            if not fasta_path.exists():
                url = (
                    "https://rest.uniprot.org/uniprotkb/stream?"
                    f"format=fasta&query=%28organism_id%3A{taxon_id}%29"
                )
                self.http.download(url, fasta_path)
            if not (db_prefix.parent / f"{db_prefix.name}.pin").exists():
                blast_fasta = self._cli_compatible_path(fasta_path)
                blast_prefix = self._cli_compatible_path(db_prefix)
                self._run(
                    [
                        self._blast_tool("makeblastdb"),
                        "-in",
                        str(blast_fasta),
                        "-dbtype",
                        "prot",
                        "-out",
                        str(blast_prefix),
                    ]
                )
            databases[species] = db_prefix
        return databases

    def search(
        self,
        query_sequence: str,
        *,
        target_accession: str,
        target_entry_name: str,
    ) -> list[BlastHit]:
        databases = self.ensure_databases()
        ortholog_databases: dict[str, Path] = {}
        all_hits: list[BlastHit] = []
        with tempfile.TemporaryDirectory() as tmpdir:
            query_path = Path(tmpdir) / "query.fasta"
            query_path.write_text(f">query\n{query_sequence}\n", encoding="utf-8")
            for species, db_prefix in databases.items():
                species_hits = self._search_database(
                    query_sequence=query_sequence,
                    query_path=query_path,
                    db_prefix=db_prefix,
                    species=species,
                    query_length=len(query_sequence),
                    target_accession=target_accession,
                    target_entry_name=target_entry_name,
                    output_path=Path(tmpdir) / f"{species}.tsv",
                )
                if not species_hits and species in self.config.cross_reactivity_fallback_species:
                    if not ortholog_databases:
                        ortholog_databases = self.ensure_ortholog_databases()
                    ortholog_prefix = ortholog_databases.get(species)
                    if ortholog_prefix is not None:
                        species_hits = self._search_database(
                            query_sequence=query_sequence,
                            query_path=query_path,
                            db_prefix=ortholog_prefix,
                            species=species,
                            query_length=len(query_sequence),
                            target_accession=target_accession,
                            target_entry_name=target_entry_name,
                            output_path=Path(tmpdir) / f"{species}_ortholog.tsv",
                        )
                all_hits.extend(species_hits)
        deduped_hits = self._dedupe_hits(all_hits)
        return sorted(deduped_hits, key=lambda hit: (-hit.bitscore, hit.evalue))

    def search_many(self, queries: list[dict[str, Any]]) -> dict[str, list[BlastHit]]:
        """Run cross-reactivity BLAST for many query sequences in bulk.

        Each query dict must contain `id` and `sequence`; `target_accession`
        and `target_entry_name` are used for self-hit filtering.
        """
        cleaned_queries = [
            {
                "id": str(item.get("id") or "").strip(),
                "sequence": str(item.get("sequence") or "").strip(),
                "target_accession": str(item.get("target_accession") or "").strip(),
                "target_entry_name": str(item.get("target_entry_name") or "").strip(),
            }
            for item in queries
            if str(item.get("id") or "").strip() and str(item.get("sequence") or "").strip()
        ]
        results: dict[str, list[BlastHit]] = {item["id"]: [] for item in cleaned_queries}
        if not cleaned_queries:
            return results

        databases = self.ensure_databases()
        ortholog_databases: dict[str, Path] = {}
        query_lengths = {item["id"]: len(item["sequence"]) for item in cleaned_queries}
        target_accessions = {item["id"]: item["target_accession"] for item in cleaned_queries}
        target_entry_names = {item["id"]: item["target_entry_name"] for item in cleaned_queries}

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            query_path = tmpdir_path / "queries.fasta"
            self._write_query_fasta(cleaned_queries, query_path)
            for species, db_prefix in databases.items():
                species_hits = self._search_database_many(
                    query_path=query_path,
                    db_prefix=db_prefix,
                    species=species,
                    query_lengths=query_lengths,
                    target_accessions=target_accessions,
                    target_entry_names=target_entry_names,
                    output_path=tmpdir_path / f"{species}.tsv",
                )
                for query_id, hits in species_hits.items():
                    results.setdefault(query_id, []).extend(hits)
                if species in self.config.cross_reactivity_fallback_species:
                    missing_queries = [
                        item
                        for item in cleaned_queries
                        if not species_hits.get(item["id"])
                    ]
                    if missing_queries:
                        if not ortholog_databases:
                            ortholog_databases = self.ensure_ortholog_databases()
                        ortholog_prefix = ortholog_databases.get(species)
                        if ortholog_prefix is not None:
                            fallback_path = tmpdir_path / f"{species}_fallback_queries.fasta"
                            self._write_query_fasta(missing_queries, fallback_path)
                            fallback_hits = self._search_database_many(
                                query_path=fallback_path,
                                db_prefix=ortholog_prefix,
                                species=species,
                                query_lengths=query_lengths,
                                target_accessions=target_accessions,
                                target_entry_names=target_entry_names,
                                output_path=tmpdir_path / f"{species}_ortholog.tsv",
                            )
                            for query_id, hits in fallback_hits.items():
                                results.setdefault(query_id, []).extend(hits)

        return {
            query_id: sorted(self._dedupe_hits(hits), key=lambda hit: (-hit.bitscore, hit.evalue))
            for query_id, hits in results.items()
        }

    def search_species(
        self,
        query_sequence: str,
        *,
        species: str,
        target_accession: str = "",
        target_entry_name: str = "",
        ortholog_search: bool = False,
    ) -> list[BlastHit]:
        databases = self.ensure_ortholog_databases() if ortholog_search else self.ensure_databases()
        db_prefix = databases.get(species)
        if db_prefix is None:
            raise BlastDatabaseError(f"No BLAST database configured for species '{species}'.")
        with tempfile.TemporaryDirectory() as tmpdir:
            query_path = Path(tmpdir) / "query.fasta"
            query_path.write_text(f">query\n{query_sequence}\n", encoding="utf-8")
            return self._search_database(
                query_sequence=query_sequence,
                query_path=query_path,
                db_prefix=db_prefix,
                species=species,
                query_length=len(query_sequence),
                target_accession=target_accession,
                target_entry_name=target_entry_name,
                output_path=Path(tmpdir) / f"{species}.tsv",
            )

    def _parse_hits(
        self,
        path: Path,
        *,
        species: str,
        query_length: int,
        target_accession: str,
        target_entry_name: str,
    ) -> list[BlastHit]:
        hits: list[BlastHit] = []
        if not path.exists():
            return hits
        for line in path.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            parts = line.split("\t", maxsplit=11)
            if len(parts) != 12:
                continue
            hit = self._hit_from_fields(
                parts,
                species=species,
                query_length=query_length,
                target_accession=target_accession,
                target_entry_name=target_entry_name,
            )
            if hit is not None:
                hits.append(hit)
        return hits

    def _parse_batch_hits(
        self,
        path: Path,
        *,
        species: str,
        query_lengths: dict[str, int],
        target_accessions: dict[str, str],
        target_entry_names: dict[str, str],
    ) -> dict[str, list[BlastHit]]:
        hits_by_query: dict[str, list[BlastHit]] = {}
        if not path.exists():
            return hits_by_query
        for line in path.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            parts = line.split("\t", maxsplit=12)
            if len(parts) != 13:
                continue
            query_id = parts[0]
            hit = self._hit_from_fields(
                parts[1:],
                species=species,
                query_length=query_lengths.get(query_id, 0),
                target_accession=target_accessions.get(query_id, ""),
                target_entry_name=target_entry_names.get(query_id, ""),
            )
            if hit is not None:
                hits_by_query.setdefault(query_id, []).append(hit)
        return hits_by_query

    def _hit_from_fields(
        self,
        parts: list[str],
        *,
        species: str,
        query_length: int,
        target_accession: str,
        target_entry_name: str,
    ) -> BlastHit | None:
        (
            subject_id,
            pident,
            length,
            qstart,
            qend,
            sstart,
            send,
            evalue,
            bitscore,
            qseq,
            sseq,
            title,
        ) = parts
        if (target_accession and target_accession in subject_id) or (
            target_entry_name and target_entry_name in subject_id
        ):
            return None
        alignment_length = int(length)
        aligned_query_span = abs(int(qend) - int(qstart)) + 1
        coverage = 100.0 * aligned_query_span / query_length if query_length else 0.0
        coverage = min(100.0, coverage)
        return BlastHit(
            subject_id=subject_id,
            description=title,
            species=_extract_species(title) or species,
            identity=float(pident),
            coverage=coverage,
            alignment_length=alignment_length,
            evalue=float(evalue),
            bitscore=float(bitscore),
            query_start=int(qstart),
            query_end=int(qend),
            subject_start=int(sstart),
            subject_end=int(send),
            query_alignment=qseq,
            subject_alignment=sseq,
            alignment_source="BLAST",
        )

    def _run(self, command: list[str]) -> None:
        try:
            subprocess.run(command, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as exc:
            raise BlastDatabaseError(exc.stderr or exc.stdout or "BLAST command failed.") from exc
        except FileNotFoundError as exc:
            raise BlastDatabaseError(f"BLAST command not found: {command[0]}") from exc

    def _blast_tool(self, name: str) -> str:
        resolved = shutil.which(name)
        if resolved:
            return resolved
        bin_dir = os.environ.get("AGDESIGN2_BLAST_BIN_DIR")
        candidates = []
        if bin_dir:
            candidates.append(Path(bin_dir) / name)
        home = Path.home()
        candidates.extend(
            [
                home / "miniforge3-openantigen" / "bin" / name,
                home / "miniforge3" / "bin" / name,
                home / "mambaforge" / "bin" / name,
                Path("/usr/local/ncbi/blast/bin") / name,
                Path("/opt/ncbi/blast/bin") / name,
            ]
        )
        for candidate in candidates:
            if candidate.is_file() and os.access(candidate, os.X_OK):
                return str(candidate)
        raise BlastDatabaseError(
            f"{name} was not found on PATH or in known BLAST locations. "
            "Set AGDESIGN2_BLAST_BIN_DIR to the directory containing blastp and makeblastdb."
        )

    def _dedupe_hits(self, hits: list[BlastHit]) -> list[BlastHit]:
        best_by_subject: dict[str, BlastHit] = {}
        for hit in hits:
            current = best_by_subject.get(hit.subject_id)
            if current is None or (hit.bitscore, hit.coverage, hit.identity) > (
                current.bitscore,
                current.coverage,
                current.identity,
            ):
                best_by_subject[hit.subject_id] = hit
        return list(best_by_subject.values())

    def _search_database(
        self,
        *,
        query_sequence: str,
        query_path: Path,
        db_prefix: Path,
        species: str,
        query_length: int,
        target_accession: str,
        target_entry_name: str,
        output_path: Path,
    ) -> list[BlastHit]:
        blast_db_prefix = self._cli_compatible_path(db_prefix)
        self._run(
            [
                self._blast_tool("blastp"),
                "-query",
                str(query_path),
                "-db",
                str(blast_db_prefix),
                "-evalue",
                str(self.config.blast_evalue),
                "-max_target_seqs",
                str(self.config.blast_max_hits),
                "-num_threads",
                str(self.config.blast_threads),
                "-outfmt",
                "6 sseqid pident length qstart qend sstart send evalue bitscore qseq sseq stitle",
                "-out",
                str(output_path),
            ]
        )
        return self._parse_hits(
            output_path,
            species=species,
            query_length=query_length,
            target_accession=target_accession,
            target_entry_name=target_entry_name,
        )

    def _search_database_many(
        self,
        *,
        query_path: Path,
        db_prefix: Path,
        species: str,
        query_lengths: dict[str, int],
        target_accessions: dict[str, str],
        target_entry_names: dict[str, str],
        output_path: Path,
    ) -> dict[str, list[BlastHit]]:
        blast_db_prefix = self._cli_compatible_path(db_prefix)
        self._run(
            [
                self._blast_tool("blastp"),
                "-query",
                str(query_path),
                "-db",
                str(blast_db_prefix),
                "-evalue",
                str(self.config.blast_evalue),
                "-max_target_seqs",
                str(self.config.blast_max_hits),
                "-num_threads",
                str(self.config.blast_threads),
                "-outfmt",
                "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qseq sseq stitle",
                "-out",
                str(output_path),
            ]
        )
        return self._parse_batch_hits(
            output_path,
            species=species,
            query_lengths=query_lengths,
            target_accessions=target_accessions,
            target_entry_names=target_entry_names,
        )

    def _write_query_fasta(self, queries: list[dict[str, str]], path: Path) -> None:
        lines: list[str] = []
        for item in queries:
            lines.append(f">{item['id']}")
            lines.append(item["sequence"])
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    def _cli_compatible_path(self, path: Path) -> Path:
        resolved = path if path.exists() else path.parent.resolve() / path.name
        text = str(resolved)
        if " " not in text:
            return resolved
        alias_root = Path(tempfile.gettempdir()) / f"agdesign2_blast_aliases_{os.getuid()}"
        alias_root.mkdir(mode=0o700, parents=True, exist_ok=True)
        if (
            alias_root.is_symlink()
            or not alias_root.is_dir()
            or alias_root.stat().st_uid != os.getuid()
        ):
            return resolved
        alias_root.chmod(0o700)
        parent = resolved.parent
        digest = hashlib.sha1(
            text.encode("utf-8"), usedforsecurity=False
        ).hexdigest()[:12]
        parent_alias = alias_root / f"{parent.name}_{digest}"
        if parent_alias.exists() and not parent_alias.is_symlink():
            return resolved
        if not parent_alias.exists():
            os.symlink(parent, parent_alias, target_is_directory=True)
        return parent_alias / resolved.name


def _extract_species(title: str) -> str | None:
    marker = " OS="
    if marker not in title:
        return None
    tail = title.split(marker, 1)[1]
    end = tail.find(" OX=")
    if end == -1:
        return tail.strip()
    return tail[:end].strip()


def extract_accession(subject_id: str) -> str | None:
    parts = subject_id.split("|")
    if len(parts) >= 2 and parts[1]:
        return parts[1]
    return None


def extract_entry_name(subject_id: str) -> str | None:
    parts = subject_id.split("|")
    if len(parts) >= 3 and parts[2]:
        return parts[2]
    return None
