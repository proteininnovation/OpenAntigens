from __future__ import annotations

import os
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from agdesign2.sequence_utils import global_align, map_query_region_to_subject


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

    def test_python_alignment_backend_can_be_forced(self) -> None:
        previous = os.environ.get("AGDESIGN2_ALIGNMENT_BACKEND")
        os.environ["AGDESIGN2_ALIGNMENT_BACKEND"] = "python"
        try:
            alignment = global_align("ABCDEFG", "ABXCDEFG")
        finally:
            if previous is None:
                os.environ.pop("AGDESIGN2_ALIGNMENT_BACKEND", None)
            else:
                os.environ["AGDESIGN2_ALIGNMENT_BACKEND"] = previous
        self.assertEqual(alignment.aligned_query, "AB-CDEFG")
        self.assertEqual(alignment.aligned_subject, "ABXCDEFG")


if __name__ == "__main__":
    unittest.main()
