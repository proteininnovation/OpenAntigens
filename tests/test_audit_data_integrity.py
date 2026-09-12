import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

from agdesign2.alphafold import AlphaFoldClient
from agdesign2.batch import BatchRow
from agdesign2.config import AnalysisConfig
from agdesign2.exceptions import ExternalServiceError
from agdesign2.http import HttpClient
from agdesign2.models import Feature, TargetRecord
from agdesign2.paralogs import _resolve_paralog_target
from agdesign2.pipeline import AntigenAnalyzer
from agdesign2.sequence_utils import find_furin_sites
from agdesign2.structure_utils import analyze_cysteines, summarize_construct_quality, build_pae_stats
from agdesign2.uniprot_prefetch import _prefetch_taxon
from tests.test_portal_db import _minimal_methionine_pdb


class AuditDataIntegrityTests(unittest.TestCase):
    def test_cursor_pagination_and_repeated_page_rejection(self):
        http = MagicMock()
        http.fetch_json_page.side_effect = [
            ({"results": [{"primaryAccession": "A"}]}, "https://example.test/?cursor=next"),
            ({"results": [{"primaryAccession": "B"}]}, None),
        ]
        with patch("agdesign2.uniprot_prefetch._cache_entry", return_value=1):
            result = _prefetch_taxon(http=http, cache=None, taxon_id=9606, page_size=25, max_pages=None)
        self.assertEqual(result["entries"], 2)
        self.assertEqual(http.fetch_json_page.call_args_list[1].args[0], "https://example.test/?cursor=next")
        http.fetch_json_page.side_effect = [
            ({"results": [{"primaryAccession": "A"}]}, "https://example.test/?cursor=next"),
            ({"results": [{"primaryAccession": "A"}]}, None),
        ]
        with patch("agdesign2.uniprot_prefetch._cache_entry", return_value=1):
            with self.assertRaisesRegex(ValueError, "repeated entries"):
                _prefetch_taxon(http=http, cache=None, taxon_id=9606, page_size=25, max_pages=None)

    def test_page_link_is_read_from_response_headers(self):
        client = HttpClient()
        with patch.object(client, "_request_response", return_value=(b'{"results": []}', {"Link": '<https://example.test/?cursor=x>; rel="next"'})):
            self.assertEqual(client.fetch_json_page("https://example.test"), ({"results": []}, "https://example.test/?cursor=x"))

    def test_download_replaces_unverified_and_corrupt_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            destination = Path(tmp) / "model"
            destination.write_bytes(b"partial")
            client = HttpClient()
            with patch.object(client, "_request", return_value=b"complete") as request:
                client.download("https://example.test/model", destination)
                client.download("https://example.test/model", destination)
                self.assertEqual(request.call_count, 1)
                destination.write_bytes(b"corrupt")
                client.download("https://example.test/model", destination)
                self.assertEqual(request.call_count, 2)
            self.assertEqual(destination.read_bytes(), b"complete")

    def test_mismatched_pae_is_rejected_before_metadata_is_written(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            http = MagicMock()
            http.fetch_json.return_value = [{"modelEntityId": "test", "sequence": "M" * 10, "pdbUrl": "https://example.test/model.pdb", "paeDocUrl": "https://example.test/pae.json"}]
            def download(url, destination):
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination.write_text(_minimal_methionine_pdb(10) if destination.suffix == ".pdb" else json.dumps({"pae": [[1] * 3 for _ in range(3)]}))
                return destination
            http.download.side_effect = download
            with self.assertRaisesRegex(ExternalServiceError, "PAE dimensions"):
                AlphaFoldClient(http, AnalysisConfig(data_dir=root)).ensure_artifacts("P00001", canonical_sequence="M" * 10)
            self.assertFalse((root / "alphafold/P00001.meta.json").exists())
            with self.assertRaisesRegex(ValueError, "boundaries exceed"):
                summarize_construct_quality(1, 10, {}, [[1] * 3 for _ in range(3)], AnalysisConfig())

    def test_overlapping_motifs_and_ambiguous_cysteines(self):
        self.assertEqual(find_furin_sites("RRRRR", "R.[KR]R"), [(1, 4, "RRRR"), (2, 5, "RRRR")])
        with self.assertRaisesRegex(ValueError, "empty"):
            find_furin_sites("RR", "R*")
        atoms = {1: {"SG": (0., 0., 0.)}, 2: {"SG": (2., 0., 0.)}, 3: {"SG": (4., 0., 0.)}}
        findings = analyze_cysteines(atoms, dict.fromkeys(atoms, "CYS"), 1, 3)
        self.assertTrue(all(f.paired_with is None and "Ambiguous" in f.warning for f in findings))
        pair = analyze_cysteines(atoms, {1: "CYS", 2: "CYS"}, 1, 2)
        self.assertEqual([f.paired_with for f in pair], [2, 1])

    def test_invalid_pae_values_and_shape(self):
        for matrix in ([[1, 2]], [[float("nan")]], [[-1]]):
            with self.subTest(matrix=matrix), self.assertRaises(ValueError):
                build_pae_stats(matrix)

    def test_secreted_paralog_uses_pipeline_region(self):
        analyzer = AntigenAnalyzer.__new__(AntigenAnalyzer)
        analyzer.config = AnalysisConfig()
        analyzer._secreted_universe_index = {"SYNTH_HUMAN"}
        target = TargetRecord(accession="SYNTH", entry_name="SYNTH_HUMAN", gene_symbol="SYNTH", protein_name="synthetic", organism="synthetic", taxon_id=9606, sequence="M" * 100)
        analyzer.precomputed_store = MagicMock()
        analyzer.precomputed_store.canonical_family_from_record.return_value = None
        analyzer.uniprot_client = MagicMock()
        analyzer.uniprot_client.resolve_target.return_value = SimpleNamespace(target=target, entry={})
        analyzer.uniprot_client.get_features.return_value = [Feature(type="CHAIN", start=20, end=100)]
        result = _resolve_paralog_target(analyzer=analyzer, batch_row=BatchRow("SYNTH_HUMAN", "uniprot_name", {}), input_index=1)
        self.assertEqual((result.ectodomain_start, result.ectodomain_end), (20, 100))
