import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from agdesign2.cache import FileCache
from agdesign2.portal import _fetch_homolog_target, _portal_align_sequences
from agdesign2.sequence_utils import AlignmentResult, global_align, map_query_region_to_subject


class PortalMappingRegressionTests(unittest.TestCase):
    def test_homolog_sequence_follows_snapshot_and_file_updates(self):
        with tempfile.TemporaryDirectory() as directory:
            for snapshot, sequence in (("first", "ACDE"), ("second", "ACDF"), ("second", "ACDG")):
                batch_dir = Path(directory) / snapshot
                cache = FileCache(batch_dir / ".agdesign2" / "cache")
                cache.path_for("uniprot_entry", "https://rest.uniprot.org/uniprotkb/P99998.json", ".json").write_text(
                    json.dumps({"primaryAccession": "P99998", "sequence": {"value": sequence}})
                )
                with self.subTest(snapshot=snapshot, sequence=sequence):
                    self.assertEqual(_fetch_homolog_target({"accession": "P99998"}, batch_dir=batch_dir)["sequence"], sequence)

    def test_portal_alignment_observes_backend_changes(self):
        with mock.patch.dict(os.environ, {"AGDESIGN2_ALIGNMENT_BACKEND": "python"}):
            _portal_align_sequences("A", "AA")
        with mock.patch.dict(os.environ, {"AGDESIGN2_ALIGNMENT_BACKEND": "auto"}):
            expected = global_align("A", "AA")
            self.assertEqual(_portal_align_sequences("A", "AA"), {
                "aligned_query": expected.aligned_query,
                "aligned_subject": expected.aligned_subject,
            })

    def test_shared_alignment_for_unrelated_sequences(self):
        import random
        rng = random.Random(47)
        for length in (40, 200, 500):
            query = "".join(rng.choices("ACDEFGHIKLMNPQRSTVWY", k=length))
            subject = query[:10] + "PPPP" + query[15:]
            with self.subTest(length=length):
                expected = global_align(query, subject)
                actual = _portal_align_sequences(query, subject)
                self.assertEqual(actual["aligned_query"], expected.aligned_query)
                self.assertEqual(actual["aligned_subject"], expected.aligned_subject)

    def test_ptprc_strict_domain_2_preserves_complete_mouse_mapping(self):
        fixture = json.loads((Path(__file__).parent / "fixtures/ptprc_mapping.json").read_text())
        mapping = _portal_align_sequences(fixture["human_sequence"], fixture["mouse_sequence"])
        alignment = AlignmentResult(**mapping, matches=0, aligned_positions=0, query_covered=0)
        start, end, sequence, _ = map_query_region_to_subject(
            alignment,
            query_start=390 - fixture["human_start"] + 1,
            query_end=481 - fixture["human_start"] + 1,
            subject_sequence=fixture["mouse_sequence"],
        )
        expected = fixture["construct"]
        self.assertEqual(start + fixture["mouse_start"] - 1, expected["start"])
        self.assertEqual(end + fixture["mouse_start"] - 1, expected["end"])
        self.assertEqual(sequence, expected["sequence"])
        self.assertEqual(len(sequence), 96)
