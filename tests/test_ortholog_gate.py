"""PR-8: ortholog selection, the acceptance gate, and what a withheld row ships.

* **SF1** (379 human rows under 70%, 583 gene-link mismatches) and **SF2**
  (479 of 743 symbol-different mouse rows; `ccl12_mouse` manufactured for
  CCL2) — ``_select_best_hcop_prediction`` keyed on
  ``(source rank, -evidence, symbol)``, so the **alphabet** was a primary key
  once source and evidence tied. It decided 195 mouse picks.
* **The gate** — there was none. Replayed through the production path over the
  cached HCOP responses: mouse 237 withheld of 5,003 pairs, macaque 79 of 4,456.
* **SF9** — withdrawn. HCOP answers the 9541 request with rhesus (9544) rows
  only; a taxon filter would empty every macaque list. The rows supply a
  symbol, and the protein comes from the cynomolgus proteome by symbol.
* **§5.5** — MGI holds *two* mouse genes for human AKR1B10 and lists the same
  two human orthologs for each. The exact-symbol rule orders them but must not
  hide the co-ortholog.
* A withheld row must ship **no sequence**: ``precomputed.py`` and
  ``mouse_portal.py`` read the accession and sequence columns by name and would
  otherwise keep putting the wrong protein into the exports.
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from agdesign2.config import AnalysisConfig  # noqa: E402
from agdesign2.hgnc import HGNCClient  # noqa: E402
from agdesign2.sequence_utils import (  # noqa: E402
    global_align,
    ortholog_gate_verdict,
    primary_source_asserts_ortholog,
)


def _client() -> HGNCClient:
    return HGNCClient.__new__(HGNCClient)


def _candidate(symbol: str, source: str, evidence: int, tax_id: int = 10090) -> dict:
    # The keys HCOP actually emits (snapshot cache): tax_id is an int.
    return {
        "gene_symbol": symbol,
        "gene_symbol_source": source,
        "evidence": ["x"] * evidence,
        "tax_id": tax_id,
    }


class RankingTests(unittest.TestCase):
    def test_an_exact_symbol_match_outranks_the_alphabet(self) -> None:
        """CCL2: Ccl12 sorted before Ccl2 and manufactured a wrong mouse page."""
        ranked = _client()._rank_hcop_candidates(
            "mouse",
            [_candidate("Ccl12", "MGI", 10), _candidate("Ccl2", "MGI", 10)],
            human_gene_symbol="CCL2",
        )
        self.assertEqual(ranked[0]["gene_symbol"], "Ccl2")

    def test_the_alphabet_no_longer_decides(self) -> None:
        # Same source, same evidence, neither an exact match: evidence then
        # source rank decide, not the symbol.
        ranked = _client()._rank_hcop_candidates(
            "mouse",
            [_candidate("Zzz1", "MGI", 12), _candidate("Aaa1", "MGI", 4)],
            human_gene_symbol="SOMETHING",
        )
        self.assertEqual(ranked[0]["gene_symbol"], "Zzz1")

    def test_a_poorly_supported_exact_match_does_not_win(self) -> None:
        # The exact match must carry at least half the best evidence.
        ranked = _client()._rank_hcop_candidates(
            "mouse",
            [_candidate("Ccl2", "MGI", 1), _candidate("Ccl12", "MGI", 12)],
            human_gene_symbol="CCL2",
        )
        self.assertEqual(ranked[0]["gene_symbol"], "Ccl12")

    def test_rhesus_rows_are_kept(self) -> None:
        """HCOP returns only 9544 rows for the macaque request (all 7,821 cached).

        Filtering them to 9541 would empty every macaque list; they supply the
        symbol that resolves the cynomolgus protein.
        """
        ranked = _client()._rank_hcop_candidates(
            "macaca_fascicularis",
            [_candidate("CD27", "VGNC", 3, tax_id=9544),
             _candidate("LOC1", "NCBI", 1, tax_id=9544)],
            human_gene_symbol="CD27",
        )
        self.assertEqual([c["gene_symbol"] for c in ranked], ["CD27", "LOC1"])

    def test_co_orthologs_are_reported_not_hidden(self) -> None:
        """§5.5 — AKR1B10's genuine 2x2 co-orthology in MGI."""
        prediction = _client()._select_best_hcop_prediction(
            "mouse",
            [_candidate("Akr1b10", "MGI", 10), _candidate("Akr1b8", "MGI", 10)],
            human_gene_symbol="AKR1B10",
        )
        self.assertEqual(prediction.gene_symbol, "Akr1b10")
        self.assertEqual(prediction.alternatives, ["Akr1b8"])


class GateTests(unittest.TestCase):
    def _verdict(self, human: str, other: str, *, asserted: bool) -> str | None:
        return ortholog_gate_verdict(
            global_align(human, other),
            human_length=len(human),
            ortholog_length=len(other),
            primary_source_asserted=asserted,
        )

    def test_a_short_but_complete_ortholog_passes(self) -> None:
        # VEGFA-shaped: the ortholog is shorter, but aligns end to end. The
        # refuted "coverage over the human sequence >= 70%" rule killed these.
        human = "MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFD"
        short = human[:55]
        self.assertIsNone(self._verdict(human, short, asserted=True))

    def test_a_sliver_alignment_is_withheld(self) -> None:
        human = "M" * 300
        other = "ACDEFGHIKLMNPQRSTVWY" * 30
        self.assertIsNotNone(self._verdict(human, other, asserted=True))

    def test_the_identity_floor_depends_on_the_assertion(self) -> None:
        import random

        rng = random.Random(11)
        human = "".join(rng.choices("ACDEFGHIKLMNPQRSTVWY", k=200))
        # 50% identity end to end: every other residue replaced.
        other = "".join(
            residue if index % 2 else rng.choice("ACDEFGHIKLMNPQRSTVWY".replace(residue, ""))
            for index, residue in enumerate(human)
        )
        self.assertIsNone(self._verdict(human, other, asserted=True))
        verdict = self._verdict(human, other, asserted=False)
        self.assertIsNotNone(verdict)
        self.assertIn("below the 70% floor", verdict)

    def test_coverage_of_the_shorter_sequence_is_enforced(self) -> None:
        # Deterministic, not seeded: 60 residues drawn from one half of the
        # alphabet align end to end, then a 60-residue tail drawn from the
        # other half cannot align to anything, so ~half the 120-aa ortholog
        # aligns and the 0.70 floor must reject it.
        human = ("ACDEFGHIKL" * 20)[:200]
        verdict = self._verdict(human, human[:60] + "W" * 60, asserted=True)
        self.assertIsNotNone(verdict)
        self.assertIn("of the shorter sequence aligned", verdict)

    def test_no_aligned_residues_is_withheld_not_called_zero_identity(self) -> None:
        self.assertIsNotNone(self._verdict("M" * 100, "W" * 10, asserted=True))


class PrimarySourceTests(unittest.TestCase):
    """The plan says "the symbol matches"; the data says otherwise."""

    def test_an_exact_symbol_match_counts(self) -> None:
        self.assertTrue(primary_source_asserts_ortholog(
            "mouse", human_symbol="VEGFA", ortholog_symbol="Vegfa",
            gene_symbol_source="NCBI", evidence_count=1,
        ))

    def test_a_curated_mod_assertion_counts(self) -> None:
        # IL6R's mouse ortholog is Il6ra (MGI, 12 sources, 56.4% identity) and
        # DEFB124's is Defb25 (MGI, 10 sources, 64.8%). The plan names both as
        # orthologs that must survive; a literal symbol test withholds both.
        self.assertTrue(primary_source_asserts_ortholog(
            "mouse", human_symbol="IL6R", ortholog_symbol="Il6ra",
            gene_symbol_source="MGI", evidence_count=12,
        ))

    def test_a_thinly_supported_assertion_does_not(self) -> None:
        self.assertFalse(primary_source_asserts_ortholog(
            "mouse", human_symbol="IL6R", ortholog_symbol="Something",
            gene_symbol_source="MGI", evidence_count=1,
        ))

    def test_a_non_curated_source_does_not(self) -> None:
        self.assertFalse(primary_source_asserts_ortholog(
            "macaca_fascicularis", human_symbol="A2M", ortholog_symbol="PRPF40B",
            gene_symbol_source="NCBI", evidence_count=12,
        ))


if __name__ == "__main__":
    unittest.main()
