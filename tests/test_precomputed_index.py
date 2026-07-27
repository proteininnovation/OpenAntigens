"""The indexed ortholog/paralog lookups must return exactly what the
former O(n) linear scans returned, including the first-in-file-order tie-break.

Each test pairs the indexed store method against an inline reference
implementation of the original linear scan and asserts identical results
(by object identity) across edge cases.
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from agdesign2.config import AnalysisConfig
from agdesign2.models import TargetRecord
from agdesign2.precomputed import PrecomputedOrthologRecord, PrecomputedReferenceStore, _clean

_OK = {"ok", "completed", ""}


def _ref_ortholog(rows, target):
    for record in rows:
        if record.status not in _OK:
            continue
        if record.entry_name and record.entry_name == target.entry_name:
            return record
        if (
            record.accession
            and target.accession
            and record.accession == target.accession
            and record.gene_symbol
            and target.gene_symbol
            and record.gene_symbol == target.gene_symbol
        ):
            return record
    return None


def _ref_paralog(rows, target):
    for row in rows:
        if str(row.get("status") or "").strip().lower() not in _OK:
            continue
        ten = _clean(row.get("target_entry_name"))
        if ten and ten == target.entry_name:
            return row
        tacc = _clean(row.get("target_accession"))
        tgene = _clean(row.get("target_gene_symbol"))
        if tacc and tacc == target.accession and tgene and target.gene_symbol and tgene == target.gene_symbol:
            return row
    return None


def _t(entry=None, accession=None, gene=None):
    return TargetRecord(
        accession=accession, entry_name=entry, gene_symbol=gene,
        protein_name="x", organism="Homo sapiens", taxon_id=9606, sequence="M" * 10,
    )


class OrthologIndexParityTests(unittest.TestCase):
    def setUp(self) -> None:
        def rec(entry, acc, gene, status):
            return PrecomputedOrthologRecord(row={
                "input_uniprot_name": entry, "human_uniprot_accession": acc,
                "human_gene_symbol": gene, "status": status,
            })

        self.rows = [
            rec("A_HUMAN", "P1", "A", "ok"),
            rec("B_HUMAN", "P2", "B", "error"),    # skipped by status
            rec("B_HUMAN", "P2", "B", "completed"),
            rec("", "P5", "E", "ok"),              # no entry name; pair-only
            rec("F_HUMAN", "P1", "A", "ok"),       # entry F but pair (P1,A) duplicates row 0
        ]
        self.store = PrecomputedReferenceStore(
            AnalysisConfig(prefer_precomputed_references=True)
        )
        self.store._ortholog_rows = self.rows

    def _check(self, target):
        self.assertIs(
            self.store.find_ortholog_record(target=target),
            _ref_ortholog(self.rows, target),
        )

    def test_parity_across_cases(self) -> None:
        self._check(_t(entry="A_HUMAN", accession="P1", gene="A"))   # entry match
        self._check(_t(entry="ZZ_HUMAN", accession="P5", gene="E"))  # pair-only match
        self._check(_t(entry="B_HUMAN"))                             # skips errored dup
        self._check(_t(entry="NOPE"))                                # no match -> None
        # tie-break: entry F matches row4, but pair (P1,A) matches earlier row0
        self._check(_t(entry="F_HUMAN", accession="P1", gene="A"))
        self._check(_t(accession="P1", gene="A"))                    # pair without entry


class ParalogIndexParityTests(unittest.TestCase):
    def setUp(self) -> None:
        self.rows = [
            {"target_entry_name": "A_HUMAN", "target_accession": "P1", "target_gene_symbol": "A", "status": "ok"},
            {"target_entry_name": "B_HUMAN", "target_accession": "P2", "target_gene_symbol": "B", "status": "failed"},
            {"target_entry_name": "B_HUMAN", "target_accession": "P2", "target_gene_symbol": "B", "status": "ok"},
            {"target_entry_name": "F_HUMAN", "target_accession": "P1", "target_gene_symbol": "A", "status": "ok"},
        ]
        self.store = PrecomputedReferenceStore(AnalysisConfig(prefer_precomputed_references=True))
        self.store._paralog_index = self.rows

    def _check(self, target):
        self.assertIs(self.store._find_paralog_row(target=target), _ref_paralog(self.rows, target))

    def test_parity_across_cases(self) -> None:
        self._check(_t(entry="A_HUMAN", accession="P1", gene="A"))
        self._check(_t(entry="B_HUMAN"))
        self._check(_t(entry="F_HUMAN", accession="P1", gene="A"))  # tie-break
        self._check(_t(entry="NOPE"))


if __name__ == "__main__":
    unittest.main()
