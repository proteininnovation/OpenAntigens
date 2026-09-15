"""Local ortholog proteome resolution: UniProt (human/mouse) + RefSeq (macaca),
no per-target NCBI. Uses on-disk fixture proteomes so it's fully hermetic."""

from __future__ import annotations

import gzip
import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from agdesign2.ortholog_proteomes import OrthologProteomeResolver, _parse_uniprot_header

_TAXA = {"human": 9606, "mouse": 10090, "macaca_fascicularis": 9541}


class UniProtHeaderTests(unittest.TestCase):
    def test_parse(self) -> None:
        acc, gene = _parse_uniprot_header("sp|P0C7M3|SFTA3_HUMAN Surfactant ... GN=SFTA3 PE=1 SV=1")
        self.assertEqual(acc, "P0C7M3")
        self.assertEqual(gene, "SFTA3")


class ResolverTests(unittest.TestCase):
    def _resolver(self, data_dir: Path) -> OrthologProteomeResolver:
        # Pre-place fixture proteomes so no download happens.
        # Ankh is indexed under the human-style symbol so the HCOP mouse symbol
        # "Ank" only resolves via the fallback_symbol path.
        (data_dir / "ortholog_mouse_swissprot.fasta").write_text(
            ">sp|P11111|CADM1_MOUSE Cell adhesion molecule 1 OS=Mus musculus OX=10090 GN=Cadm1 PE=1 SV=1\nMOUSECADM1SEQ\n"
            ">sp|P22222|OPRM_MOUSE Mu opioid receptor OS=Mus musculus OX=10090 GN=Oprm1 PE=1 SV=1\nMOUSEOPRM1SEQ\n"
            ">sp|P33333|ANKH_MOUSE Progressive ankylosis OS=Mus musculus OX=10090 GN=Ankh PE=1 SV=1\nMOUSEANKHSEQ\n"
            ">sp|Q64444|CAH4_MOUSE Carbonic anhydrase 4 OS=Mus musculus OX=10090 GN=Ca4 PE=1 SV=1\nMOUSECA4REVIEWED\n",
            encoding="utf-8",
        )
        # Unreviewed (TrEMBL) reference-proteome fallback: Abcc8 has no reviewed
        # entry; the Cadm1 tr entry must never override the reviewed one.
        # A0A003 is a longer Abcc8 isoform than A0A001; longest must win within TrEMBL.
        (data_dir / "ortholog_mouse_trembl.fasta").write_text(
            ">tr|A0A001|A0A001_MOUSE ATP-binding cassette OS=Mus musculus OX=10090 GN=Abcc8 PE=4 SV=1\nMOUSEABCC8SEQ\n"
            ">tr|A0A003|A0A003_MOUSE ATP-binding cassette isoform OS=Mus musculus OX=10090 GN=Abcc8 PE=4 SV=1\nMOUSEABCC8SEQLONGER\n"
            ">tr|A0A002|A0A002_MOUSE Cell adhesion molecule 1 OS=Mus musculus OX=10090 GN=Cadm1 PE=4 SV=1\nTREMBLCADM1DECOY\n"
            ">tr|F6ST32|F6ST32_MOUSE Carbonic anhydrase 4 (Fragment) OS=Mus musculus OX=10090 GN=Car4 PE=1 SV=3\nCA4FRAGMENT\n",
            encoding="utf-8",
        )
        # macaca RefSeq fixture: protein FASTA + feature table (CDS rows).
        with gzip.open(data_dir / "ortholog_macaca_fascicularis_protein.faa.gz", "wt", encoding="utf-8") as fh:
            fh.write(">XP_0001.1 cell adhesion molecule 1 [Macaca fascicularis]\nMACCADM1SEQLONGER\n")
            fh.write(">XP_0002.1 mu opioid receptor [Macaca fascicularis]\nMACOPRM1SEQ\n")
        cols = ["# feature", "class", "assembly", "assembly_unit", "seq_type", "chromosome",
                "genomic_accession", "start", "end", "strand", "product_accession",
                "non-redundant_refseq", "related_accession", "name", "symbol", "GeneID"]
        with gzip.open(data_dir / "ortholog_macaca_fascicularis_feature_table.txt.gz", "wt", encoding="utf-8") as fh:
            fh.write("\t".join(cols) + "\n")
            row = lambda product, symbol, gid: "\t".join(
                ["CDS", "with_protein", "", "", "", "1", "NC_1", "1", "9", "+", product, "", "", "", symbol, gid]
            )
            fh.write(row("XP_0001.1", "CADM1", "111") + "\n")
            fh.write(row("XP_0002.1", "OPRM1", "222") + "\n")
        return OrthologProteomeResolver(
            http=None, data_dir=data_dir, species_taxonomy=_TAXA
        )

    def test_mouse_uniprot_by_symbol(self) -> None:
        with TemporaryDirectory() as tmp:
            r = self._resolver(Path(tmp))
            hit = r.resolve("mouse", gene_symbol="Cadm1")
            self.assertIsNotNone(hit)
            self.assertEqual(hit.accession, "P11111")
            self.assertEqual(hit.sequence, "MOUSECADM1SEQ")
            self.assertEqual(hit.source, "uniprot")
            # case-insensitive
            self.assertIsNotNone(r.resolve("mouse", gene_symbol="oprm1"))
            # missing gene -> None (not an error)
            self.assertIsNone(r.resolve("mouse", gene_symbol="NOSUCHGENE"))

    def test_human_symbol_fallback(self) -> None:
        with TemporaryDirectory() as tmp:
            r = self._resolver(Path(tmp))
            # HCOP mouse symbol "Ank" misses; human symbol "ANKH" recovers it.
            self.assertIsNone(r.resolve("mouse", gene_symbol="Ank"))
            hit = r.resolve("mouse", gene_symbol="Ank", fallback_symbol="ANKH")
            self.assertIsNotNone(hit)
            self.assertEqual(hit.accession, "P33333")
            self.assertEqual(hit.source, "uniprot")

    def test_trembl_fallback_and_reviewed_precedence(self) -> None:
        with TemporaryDirectory() as tmp:
            r = self._resolver(Path(tmp))
            # Abcc8 has no reviewed entry -> resolves to the longest unreviewed TrEMBL isoform.
            tre = r.resolve("mouse", gene_symbol="Abcc8")
            self.assertIsNotNone(tre)
            self.assertEqual(tre.accession, "A0A003")
            self.assertEqual(tre.sequence, "MOUSEABCC8SEQLONGER")
            self.assertEqual(tre.source, "uniprot-trembl")
            # A TrEMBL entry never overrides the reviewed SwissProt protein.
            cadm1 = r.resolve("mouse", gene_symbol="Cadm1")
            self.assertEqual(cadm1.accession, "P11111")
            self.assertEqual(cadm1.source, "uniprot")

    def test_reviewed_record_wins_across_candidate_symbols(self) -> None:
        with TemporaryDirectory() as tmp:
            r = self._resolver(Path(tmp))
            hit = r.resolve("mouse", gene_symbol="Car4", fallback_symbol="CA4")
            self.assertEqual(hit.accession, "Q64444")
            self.assertEqual(hit.source, "uniprot")

    def test_macaca_refseq_by_symbol_and_gene_id(self) -> None:
        with TemporaryDirectory() as tmp:
            r = self._resolver(Path(tmp))
            by_sym = r.resolve("macaca_fascicularis", gene_symbol="CADM1")
            self.assertEqual(by_sym.accession, "XP_0001.1")
            self.assertEqual(by_sym.sequence, "MACCADM1SEQLONGER")
            self.assertEqual(by_sym.source, "refseq")
            # GeneID lookup also works (robust to symbol synonyms)
            by_id = r.resolve("macaca_fascicularis", gene_id="222")
            self.assertEqual(by_id.accession, "XP_0002.1")


class AssemblyPathTests(unittest.TestCase):
    def test_ftp_base_construction(self) -> None:
        base = OrthologProteomeResolver._assembly_ftp_base("GCF_037993035.2", "T2T-MFA8v1.1")
        self.assertEqual(
            base,
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/037/993/035/GCF_037993035.2_T2T-MFA8v1.1",
        )


if __name__ == "__main__":
    unittest.main()
