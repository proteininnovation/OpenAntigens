"""Local, bulk-download ortholog protein resolution.

Replaces the per-target NCBI E-utilities queries (esearch/esummary/efetch) that
the ortholog table previously issued for every mouse/macaque ortholog — those
trip NCBI rate limiting under any concurrency and silently lose ~a third of the
orthologs. Instead we download a small number of bulk proteomes once and resolve
every ortholog by a local dictionary lookup:

- human, mouse: UniProt SwissProt (reviewed) proteome, gene-symbol indexed.
- macaca_fascicularis: the current RefSeq assembly's protein set, indexed by
  gene symbol and NCBI GeneID via the assembly feature table (one bulk download).

The only NCBI traffic is the single macaque RefSeq assembly download; human and
mouse are pure UniProt bulk streams.
"""

from __future__ import annotations

import gzip
import io
import json
import urllib.parse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

from .http import HttpClient

_UNIPROT_STREAM = "https://rest.uniprot.org/uniprotkb/stream"
_DATASETS_API = "https://api.ncbi.nlm.nih.gov/datasets/v2"
_NCBI_FTP_GENOMES = "https://ftp.ncbi.nlm.nih.gov/genomes/all"

# Species that resolve from UniProt SwissProt; everything else uses RefSeq.
_UNIPROT_SPECIES = {"human", "mouse"}

# Species for which we add an unreviewed (TrEMBL) fallback from the UniProt
# reference proteome. Many real genes (e.g. mouse Abcc8, Adamts9) have no
# reviewed SwissProt entry; reviewed entries stay authoritative and TrEMBL only
# fills genes with no reviewed protein. Maps species -> reference proteome ID.
_UNIPROT_REFERENCE_PROTEOME = {"mouse": "UP000000589"}


@dataclass(slots=True)
class ResolvedOrthologProtein:
    accession: str
    sequence: str
    source: str  # "uniprot" (reviewed) | "uniprot-trembl" (unreviewed) | "refseq"
    organism: str


def _open_maybe_gz(path: Path) -> io.TextIOBase:
    if path.suffix == ".gz" or path.name.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8")
    return path.open("r", encoding="utf-8")


def _iter_fasta(path: Path) -> Iterator[tuple[str, str]]:
    header: str | None = None
    chunks: list[str] = []
    with _open_maybe_gz(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def _parse_uniprot_header(header: str) -> tuple[str | None, str | None]:
    """``sp|P0C7M3|SFTA3_HUMAN ... GN=SFTA3 PE=1 SV=1`` -> (accession, gene)."""
    accession = None
    parts = header.split("|")
    if len(parts) >= 3:
        accession = parts[1].strip()
    gene = None
    for token in header.split():
        if token.startswith("GN="):
            gene = token[3:].strip()
            break
    return accession, gene


class OrthologProteomeResolver:
    def __init__(
        self,
        *,
        http: HttpClient,
        data_dir: Path,
        species_taxonomy: dict[str, int],
        verbose: bool = False,
    ) -> None:
        self.http = http
        self.data_dir = Path(data_dir)
        self.species_taxonomy = species_taxonomy
        self.verbose = verbose
        self.data_dir.mkdir(parents=True, exist_ok=True)
        # species -> {"by_symbol": {GENE: protein}, "by_gene_id": {id: protein}}
        self._indices: dict[str, dict[str, dict[str, ResolvedOrthologProtein]]] = {}

    def prewarm(self, species: list[str]) -> None:
        """Build (download + index) the given species' proteomes up front.

        Called once before the parallel ortholog build so the heavy download/parse
        happens serially and worker threads only do lock-free dict reads.
        """
        for name in species:
            if name != "human":  # human reference is the resolved target itself
                self._species_index(name)

    def resolve(
        self,
        species: str,
        *,
        gene_symbol: str | None = None,
        gene_id: str | int | None = None,
        fallback_symbol: str | None = None,
    ) -> ResolvedOrthologProtein | None:
        """Resolve an ortholog by GeneID, then species symbol, then a fallback.

        ``fallback_symbol`` (typically the human gene symbol) recovers orthologs
        where HCOP reports a species-specific symbol (e.g. mouse ``Ank``,
        ``Car10``) but the proteome indexes the protein under the human-style
        symbol (``ANKH``, ``CA10``).
        """
        index = self._species_index(species)
        candidates: list[tuple[int, ResolvedOrthologProtein]] = []
        if gene_id is not None:
            hit = index["by_gene_id"].get(str(gene_id).strip())
            if hit is not None:
                candidates.append((0, hit))
        for position, symbol in enumerate((gene_symbol, fallback_symbol), start=1):
            if symbol:
                hit = index["by_symbol"].get(symbol.strip().upper())
                if hit is not None:
                    candidates.append((position, hit))
        if not candidates:
            return None
        return min(candidates, key=lambda item: (item[1].source == "uniprot-trembl", item[0]))[1]

    def _log(self, message: str) -> None:
        if self.verbose:
            print(f"[ortholog-proteomes] {message}", flush=True)

    def _species_index(self, species: str) -> dict[str, dict[str, ResolvedOrthologProtein]]:
        if species not in self._indices:
            if species in _UNIPROT_SPECIES:
                self._indices[species] = self._build_uniprot_index(species)
            else:
                self._indices[species] = self._build_refseq_index(species)
        return self._indices[species]

    # ---- UniProt (human, mouse) ---------------------------------------------
    def _build_uniprot_index(self, species: str) -> dict[str, dict[str, ResolvedOrthologProtein]]:
        taxon = self.species_taxonomy[species]
        organism = species.replace("_", " ").title()

        # Reviewed SwissProt: authoritative, indexed first so it always wins.
        fasta_path = self.data_dir / f"ortholog_{species}_swissprot.fasta"
        if not fasta_path.exists():
            self._log(f"downloading UniProt SwissProt proteome for {species} (taxon {taxon})")
            self.http.download(self._uniprot_stream_url(f"(reviewed:true) AND (organism_id:{taxon})"), fasta_path)
        by_symbol: dict[str, ResolvedOrthologProtein] = {}
        self._index_uniprot_fasta(fasta_path, by_symbol, organism=organism, source="uniprot")
        reviewed_count = len(by_symbol)

        # Unreviewed TrEMBL fallback (reference proteome): fills genes that have
        # no reviewed entry without ever overriding a reviewed one.
        proteome = _UNIPROT_REFERENCE_PROTEOME.get(species)
        if proteome is not None:
            trembl_path = self.data_dir / f"ortholog_{species}_trembl.fasta"
            if not trembl_path.exists():
                self._log(f"downloading UniProt TrEMBL fallback for {species} (reference proteome {proteome})")
                self.http.download(
                    self._uniprot_stream_url(
                        f"(organism_id:{taxon}) AND (proteome:{proteome}) AND (reviewed:false)"
                    ),
                    trembl_path,
                )
            self._index_uniprot_fasta(
                trembl_path, by_symbol, organism=organism, source="uniprot-trembl", only_new=True
            )
            self._log(
                f"{species}: indexed {reviewed_count} SwissProt genes "
                f"+ {len(by_symbol) - reviewed_count} TrEMBL-only fallback genes"
            )
        else:
            self._log(f"{species}: indexed {reviewed_count} SwissProt genes")
        return {"by_symbol": by_symbol, "by_gene_id": {}}

    @staticmethod
    def _uniprot_stream_url(query: str) -> str:
        return f"{_UNIPROT_STREAM}?format=fasta&query={urllib.parse.quote(query)}"

    @staticmethod
    def _index_uniprot_fasta(
        fasta_path: Path,
        by_symbol: dict[str, ResolvedOrthologProtein],
        *,
        organism: str,
        source: str,
        only_new: bool = False,
    ) -> None:
        for header, sequence in _iter_fasta(fasta_path):
            accession, gene = _parse_uniprot_header(header)
            if not gene or not sequence:
                continue
            key = gene.upper()
            existing = by_symbol.get(key)
            if only_new and existing is not None and existing.source != source:
                continue  # never override a higher-priority (reviewed) entry
            # Prefer the longer (canonical) sequence when a gene maps to several.
            if existing is None or len(sequence) > len(existing.sequence):
                by_symbol[key] = ResolvedOrthologProtein(
                    accession=accession or gene, sequence=sequence, source=source, organism=organism
                )

    # ---- RefSeq (macaca_fascicularis) ---------------------------------------
    def _build_refseq_index(self, species: str) -> dict[str, dict[str, ResolvedOrthologProtein]]:
        taxon = self.species_taxonomy[species]
        faa = self.data_dir / f"ortholog_{species}_protein.faa.gz"
        feature_table = self.data_dir / f"ortholog_{species}_feature_table.txt.gz"
        if not faa.exists() or not feature_table.exists():
            accession, name = self._resolve_refseq_assembly(taxon)
            base = self._assembly_ftp_base(accession, name)
            self._log(f"downloading {species} RefSeq proteins ({accession} {name})")
            if not faa.exists():
                self.http.download(f"{base}/{accession}_{name}_protein.faa.gz", faa)
            if not feature_table.exists():
                self.http.download(f"{base}/{accession}_{name}_feature_table.txt.gz", feature_table)

        sequences: dict[str, str] = {}
        for header, sequence in _iter_fasta(faa):
            acc = header.split()[0] if header else ""
            if acc:
                sequences[acc] = sequence
                sequences[acc.split(".")[0]] = sequence  # version-agnostic

        by_symbol: dict[str, ResolvedOrthologProtein] = {}
        by_gene_id: dict[str, ResolvedOrthologProtein] = {}
        organism = "Macaca fascicularis"
        with _open_maybe_gz(feature_table) as handle:
            header_cols = handle.readline().lstrip("#").strip().split("\t")
            idx = {col.strip(): i for i, col in enumerate(header_cols)}
            f_feature = idx.get("feature", 0)
            f_product = idx.get("product_accession")
            f_symbol = idx.get("symbol")
            f_geneid = idx.get("GeneID")
            for line in handle:
                cols = line.rstrip("\n").split("\t")
                if len(cols) <= f_feature or cols[f_feature] != "CDS":
                    continue
                product = cols[f_product].strip() if f_product is not None and f_product < len(cols) else ""
                if not product:
                    continue
                sequence = sequences.get(product) or sequences.get(product.split(".")[0])
                if not sequence:
                    continue
                symbol = cols[f_symbol].strip() if f_symbol is not None and f_symbol < len(cols) else ""
                gene_id = cols[f_geneid].strip() if f_geneid is not None and f_geneid < len(cols) else ""
                protein = ResolvedOrthologProtein(
                    accession=product, sequence=sequence, source="refseq", organism=organism
                )
                # Prefer curated NP_ over predicted XP_, then the longest isoform.
                def _better(new: ResolvedOrthologProtein, old: ResolvedOrthologProtein | None) -> bool:
                    if old is None:
                        return True
                    new_np = new.accession.startswith("NP_")
                    old_np = old.accession.startswith("NP_")
                    if new_np != old_np:
                        return new_np
                    return len(new.sequence) > len(old.sequence)

                if symbol and _better(protein, by_symbol.get(symbol.upper())):
                    by_symbol[symbol.upper()] = protein
                if gene_id and _better(protein, by_gene_id.get(gene_id)):
                    by_gene_id[gene_id] = protein
        self._log(f"{species}: indexed {len(by_symbol)} RefSeq genes ({len(by_gene_id)} by GeneID)")
        return {"by_symbol": by_symbol, "by_gene_id": by_gene_id}

    def _resolve_refseq_assembly(self, taxon: int) -> tuple[str, str]:
        url = (
            f"{_DATASETS_API}/genome/taxon/{taxon}/dataset_report"
            "?filters.assembly_source=refseq&filters.has_annotation=true&page_size=5"
        )
        payload = self.http.fetch_json(url, cache_namespace="ncbi_datasets_report")
        for report in payload.get("reports", []) or []:
            accession = report.get("accession")
            name = (report.get("assembly_info") or {}).get("assembly_name")
            if accession and name:
                return accession, name.replace(" ", "_")
        raise RuntimeError(f"No annotated RefSeq assembly found for taxon {taxon}")

    @staticmethod
    def _assembly_ftp_base(accession: str, name: str) -> str:
        # GCF_037993035.2 -> GCF/037/993/035/GCF_037993035.2_<name>
        prefix, digits = accession.split("_", 1)
        digits = digits.split(".")[0]
        parts = [digits[i : i + 3] for i in range(0, 9, 3)]
        return f"{_NCBI_FTP_GENOMES}/{prefix}/{'/'.join(parts)}/{accession}_{name}"
