from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from agdesign2.gpcrdb import GPCRdbClient


class FakeHttp:
    def fetch_json(self, url: str, **kwargs):
        if "/protein/oprm_human/" in url:
            return {
                "entry_name": "oprm_human",
                "name": "&mu; receptor",
                "accession": "P35372",
                "family": "001_002_022_003",
                "species": "Homo sapiens",
                "source": "SWISSPROT",
                "residue_numbering_scheme": "GPCRdb(A)",
            }
        if "/residues/oprm_human/" in url:
            return [
                {"sequence_number": 1, "amino_acid": "M", "protein_segment": "N-term", "display_generic_number": None},
                {"sequence_number": 2, "amino_acid": "D", "protein_segment": "TM1", "display_generic_number": "1.50x50"},
                {"sequence_number": 3, "amino_acid": "R", "protein_segment": "TM3", "display_generic_number": "3.50x50"},
                {"sequence_number": 4, "amino_acid": "Y", "protein_segment": "TM3", "display_generic_number": "3.51x51"},
            ]
        if "/proteinfamily/001_002_022_003/" in url:
            return {"slug": "001_002_022_003", "name": "&mu; receptor", "parent": {"slug": "001_002_022"}}
        if "/proteinfamily/001_002_022/" in url:
            return {"slug": "001_002_022", "name": "Opioid receptors", "parent": {"slug": "001"}}
        if "/proteinfamily/001/" in url:
            return {"slug": "001", "name": "Class A (Rhodopsin)", "parent": {"slug": "000"}}
        if "/proteinfamily/000/" in url:
            return {"slug": "000", "name": "Parent family", "parent": None}
        raise AssertionError(url)


class GPCRdbClientTests(unittest.TestCase):
    def test_fetch_annotation_builds_segments_family_path_and_motifs(self) -> None:
        annotation = GPCRdbClient(FakeHttp()).fetch_annotation(
            entry_name="OPRM_HUMAN",
            accession="P35372",
            sequence="MDRY",
        )

        self.assertIsNotNone(annotation)
        assert annotation is not None
        self.assertEqual(annotation.entry_name, "oprm_human")
        self.assertEqual(annotation.name, "μ receptor")
        self.assertEqual(annotation.family_path, ["Class A (Rhodopsin)", "Opioid receptors", "μ receptor"])
        self.assertEqual([segment.name for segment in annotation.segments], ["N-term", "TM1", "TM3"])
        self.assertTrue(any(motif.name == "DRY / E/DRY activation motif" for motif in annotation.conserved_motifs))


if __name__ == "__main__":
    unittest.main()
