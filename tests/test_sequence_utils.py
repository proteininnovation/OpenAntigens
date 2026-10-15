from __future__ import annotations

import random
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from agdesign2.sequence_utils import global_align, map_query_region_to_subject

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


class SequenceMappingTests(unittest.TestCase):
    def test_global_align_reports_identity_and_coverage(self) -> None:
        alignment = global_align("ABCDEFG", "ABXDEFG")
        self.assertEqual(alignment.matches, 6)
        self.assertEqual(alignment.aligned_positions, 7)
        self.assertAlmostEqual(alignment.identity, 85.71428571428571)
        self.assertEqual(alignment.coverage, 100.0)

    def test_maps_query_region_to_subject_region(self) -> None:
        alignment = global_align("ABCDEFG", "ABXDEFG")
        start, end, sequence, notes = map_query_region_to_subject(
            alignment,
            query_start=3,
            query_end=5,
            subject_sequence="ABXDEFG",
        )
        self.assertEqual((start, end, sequence), (3, 5, "XDE"))
        self.assertEqual(notes, [])

    def test_rejects_severely_incomplete_boundary_projection(self) -> None:
        alignment = global_align("ABCDEFGHIJ", "ABCDEF")
        start, end, sequence, notes = map_query_region_to_subject(
            alignment,
            query_start=1,
            query_end=10,
            subject_sequence="ABCDEF",
        )
        self.assertEqual((start, end, sequence), (None, None, None))
        self.assertEqual(
            notes,
            ["Boundary projection rejected: 6/10 query residues aligned to the homolog (60.0% < 70.0% minimum)."],
        )

    def test_accepts_boundary_projection_at_seventy_percent_coverage(self) -> None:
        alignment = global_align("ABCDEFGHIJ", "ABCDEFG")
        start, end, sequence, notes = map_query_region_to_subject(
            alignment,
            query_start=1,
            query_end=10,
            subject_sequence="ABCDEFG",
        )
        self.assertEqual((start, end, sequence), (1, 7, "ABCDEFG"))
        self.assertEqual(notes, ["3 construct residues align to gaps in the homolog ectodomain."])

    def test_there_is_one_scoring_scheme_and_no_backend_switch(self) -> None:
        """F33: three backends with three tie-breaks made pinned numbers fragile
        and let a silent fallback hide the scoring defect. One backend now."""
        source = (ROOT / "src" / "agdesign2" / "sequence_utils.py").read_text(encoding="utf-8")
        self.assertNotIn("AGDESIGN2_ALIGNMENT_BACKEND", source)
        self.assertNotIn("import parasail", source)

    def test_an_internal_insertion_is_aligned_as_a_gap(self) -> None:
        # Long enough that paying -11 for the internal gap beats sliding the
        # whole query off a free end gap; on a 7-residue query the two are tied.
        query = "MKWVTFISLLFLFSSAYSRG"
        subject = "MKWVTFISLLFPLFSSAYSRG"
        alignment = global_align(query, subject)
        self.assertEqual(alignment.aligned_query, "MKWVTFISLLF-LFSSAYSRG")
        self.assertEqual(alignment.aligned_subject, subject)
        self.assertEqual(alignment.matches, len(query))


class IdentityMetricGuardTests(unittest.TestCase):
    """F33: the scheme that was shipped made unrelated sequences look homologous.

    match +1 / **mismatch 0** / gap -1 let a short query cherry-pick matches out
    of a long subject. Measured on the shipped scheme: 20 seeded random pairs of
    83 vs 550 aa averaged **75.6%** identity at 100% coverage, and 200 vs 1,866
    aa averaged **86.9%**. Those numbers were published: a wrong-gene macaque
    row read 79.7%, and the mouse family matrix put GPR63 at 79.31% to LGR4.
    """

    def test_unrelated_sequences_of_very_different_length_do_not_look_homologous(self) -> None:
        rng = random.Random(0)
        for short_length, long_length in ((83, 550), (200, 1866)):
            for _ in range(20):
                query = "".join(rng.choices(AMINO_ACIDS, k=short_length))
                subject = "".join(rng.choices(AMINO_ACIDS, k=long_length))
                alignment = global_align(query, subject)
                with self.subTest(lengths=(short_length, long_length)):
                    # The directional-matrix value (matches / row length).
                    self.assertLess(100.0 * alignment.matches / len(query), 30.0)

    def test_equal_length_unrelated_pairs_stay_low(self) -> None:
        rng = random.Random(1)
        for _ in range(10):
            query = "".join(rng.choices(AMINO_ACIDS, k=300))
            subject = "".join(rng.choices(AMINO_ACIDS, k=300))
            alignment = global_align(query, subject)
            self.assertLess(100.0 * alignment.matches / len(query), 25.0)

    def test_a_contained_sequence_scores_100_percent_over_its_own_length(self) -> None:
        """The semantics tests/test_hgnc.py's directional matrix relies on."""
        rng = random.Random(2)
        core = "".join(rng.choices(AMINO_ACIDS, k=200))
        flank_left = "".join(rng.choices(AMINO_ACIDS, k=30))
        flank_right = "".join(rng.choices(AMINO_ACIDS, k=30))
        alignment = global_align(core, flank_left + core + flank_right)
        self.assertEqual(alignment.matches, len(core))
        self.assertEqual(100.0 * alignment.matches / len(core), 100.0)

    def test_gpr63_to_lgr4_is_not_a_family_relationship(self) -> None:
        """The cell that motivated F33, on the real snapshot sequences.

        The mouse family matrix published GPR63 -> LGR4 at **79.31%**: an 87-aa
        ectodomain against a 520-aa one, two unrelated GPCR/LRR-GPCR
        ectodomains. The honest value is 17.24%.
        """
        import json

        reports = Path(
            "/Users/Andre.Teixeira/Code/AgDesign2/snapshots/openantigen-2026-07-23"
            "/mouse_data/reports"
        )
        if not reports.is_dir():  # the snapshot is not part of the repo
            self.skipTest("mouse snapshot not available")

        def ectodomain(entry: str) -> str:
            report = json.loads((reports / f"{entry}_mouse_report.json").read_text(encoding="utf-8"))
            region = report["ectodomain"]
            return str(report["target"]["sequence"])[int(region["start"]) - 1 : int(region["end"])]

        row, column = ectodomain("gpr63"), ectodomain("lgr4")
        self.assertEqual((len(row), len(column)), (87, 520), "fixture drifted")
        alignment = global_align(row, column)
        self.assertLess(100.0 * alignment.matches / len(row), 50.0)

    def test_characters_outside_blosum62_are_scored_as_unknown_not_crashed(self) -> None:
        # Selenocysteine (U), pyrrolysine (O) and J are outside BLOSUM62's
        # alphabet. They must score, and must still match themselves.
        alignment = global_align("MKUOJWVT", "MKUOJWVT")
        self.assertEqual(alignment.matches, 8)
        self.assertEqual(alignment.identity, 100.0)


if __name__ == "__main__":
    unittest.main()
