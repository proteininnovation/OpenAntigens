"""The cached ortholog-table RefSeq lookup must match the former
per-call linear scan exactly, including the empty-sequence short-circuit."""

from __future__ import annotations

import csv
import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from agdesign2.portal import _fetch_homolog_target, _fetch_refseq_target_from_ortholog_table

_SPECIES = (
    ("human_gene_symbol", "human_refseq_accession", "human_refseq_sequence"),
    ("mouse_gene_symbol", "mouse_refseq_accession", "mouse_refseq_sequence"),
    ("macaca_fascicularis_gene_symbol", "macaca_fascicularis_refseq_accession", "macaca_fascicularis_refseq_sequence"),
)


def _ref_scan(table_path: Path, accession: str):
    with table_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            for gene_column, accession_column, sequence_column in _SPECIES:
                if str(row.get(accession_column) or "") != accession:
                    continue
                sequence = str(row.get(sequence_column) or "")
                if not sequence:
                    return None
                return {
                    "accession": accession,
                    "entry_name": accession,
                    "gene_symbol": row.get(gene_column) or accession,
                    "sequence": sequence,
                }
    return None


class OrthologRefseqIndexParityTests(unittest.TestCase):
    def _write_table(self, batch_dir: Path) -> None:
        fields = [
            "human_gene_symbol", "human_refseq_accession", "human_refseq_sequence",
            "mouse_gene_symbol", "mouse_refseq_accession", "mouse_refseq_sequence",
            "macaca_fascicularis_gene_symbol", "macaca_fascicularis_refseq_accession", "macaca_fascicularis_refseq_sequence",
        ]
        with (batch_dir / "ortholog_reference_table.tsv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            # NP_1 found in human with sequence; XP_9 found in macaca with sequence.
            writer.writerow({"human_gene_symbol": "AAA", "human_refseq_accession": "NP_1", "human_refseq_sequence": "MAAA",
                             "mouse_refseq_accession": "NP_2", "mouse_refseq_sequence": "MBBB",
                             "macaca_fascicularis_refseq_accession": "XP_9", "macaca_fascicularis_refseq_sequence": "MCCC",
                             "macaca_fascicularis_gene_symbol": "CCC"})
            # NP_3 first appears with EMPTY sequence (short-circuit to None)...
            writer.writerow({"human_refseq_accession": "NP_3", "human_refseq_sequence": ""})
            # ...even though a later row has NP_3 with a sequence.
            writer.writerow({"human_gene_symbol": "DDD", "human_refseq_accession": "NP_3", "human_refseq_sequence": "MDDD"})

    def test_parity(self) -> None:
        with TemporaryDirectory() as tmp:
            batch_dir = Path(tmp)
            self._write_table(batch_dir)
            table_path = batch_dir / "ortholog_reference_table.tsv"
            for acc in ["NP_1", "NP_2", "XP_9", "NP_3", "NP_404", ""]:
                self.assertEqual(
                    _fetch_refseq_target_from_ortholog_table(acc, batch_dir=batch_dir),
                    _ref_scan(table_path, acc),
                    f"mismatch for {acc!r}",
                )

    def test_returned_dict_is_a_copy(self) -> None:
        with TemporaryDirectory() as tmp:
            batch_dir = Path(tmp)
            self._write_table(batch_dir)
            first = _fetch_refseq_target_from_ortholog_table("NP_1", batch_dir=batch_dir)
            first["sequence"] = "MUTATED"
            second = _fetch_refseq_target_from_ortholog_table("NP_1", batch_dir=batch_dir)
            self.assertEqual(second["sequence"], "MAAA")


class HomologTargetTableFallbackTests(unittest.TestCase):
    """A UniProt (non-RefSeq) homolog accession must resolve from the ortholog
    table when it is not in the live cache — the human sequence in the mouse
    portal's Live Selected Region viewer depends on this."""

    def test_human_uniprot_accession_resolves_via_table(self) -> None:
        with TemporaryDirectory() as tmp:
            batch_dir = Path(tmp)
            fields = [
                "human_gene_symbol", "human_refseq_accession", "human_refseq_sequence",
                "mouse_gene_symbol", "mouse_refseq_accession", "mouse_refseq_sequence",
                "macaca_fascicularis_gene_symbol", "macaca_fascicularis_refseq_accession", "macaca_fascicularis_refseq_sequence",
            ]
            with (batch_dir / "ortholog_reference_table.tsv").open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
                writer.writeheader()
                writer.writerow({
                    "human_gene_symbol": "CADM1", "human_refseq_accession": "Q9ZTEST1", "human_refseq_sequence": "MHUMANSEQ",
                    "mouse_gene_symbol": "Cadm1", "mouse_refseq_accession": "Q9ZTEST2", "mouse_refseq_sequence": "MMOUSESEQ",
                })
            record = {"species": "human", "accession": "Q9ZTEST1", "entry_name": "CADM1_HUMAN", "available": True}
            target = _fetch_homolog_target(record, batch_dir=batch_dir)
            self.assertIsNotNone(target)
            self.assertEqual(target["sequence"], "MHUMANSEQ")
            self.assertEqual(target["gene_symbol"], "CADM1")


if __name__ == "__main__":
    unittest.main()
