from __future__ import annotations

import json
import sys
import tempfile
import threading
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from agdesign2.config import AnalysisConfig
from agdesign2.family_alignments import _assign_unique_labels
from agdesign2.models import TargetRecord
from agdesign2.models import CanonicalFamily, DomainAnnotation, Feature
from agdesign2.paralogs import EnsemblParalogClient, ParalogMember, build_paralog_reference_from_tsv
from agdesign2.precomputed import PrecomputedReferenceStore
from agdesign2.uniprot import TargetResolution


class FakeHttpClient:
    def fetch_json(self, url: str, *, headers=None, cache_namespace="json", cache_key=None):
        if "/homology/symbol/human/EGFR" in url:
            return {
                "data": [
                    {
                        "id": "ENSG00000146648",
                        "homologies": [
                            {
                                "id": "ENSG00000141736",
                                "type": "within_species_paralog",
                                "taxonomy_level": "Vertebrata",
                                "target": {
                                    "id": "ENSG00000141736",
                                    "protein_id": "ENSP00000275493",
                                    "perc_id": 62.5,
                                },
                                "target_perc_id": 58.1,
                            }
                        ],
                    }
                ]
            }
        if "/lookup/id/ENSG00000141736" in url:
            return {"display_name": "ERBB2"}
        raise AssertionError(f"Unexpected URL: {url}")


class ParalogClientTests(unittest.TestCase):
    def test_fetch_human_paralogs_resolves_gene_symbol_via_lookup(self) -> None:
        client = EnsemblParalogClient(FakeHttpClient())
        predictions = client.fetch_human_paralogs("EGFR")
        self.assertEqual(len(predictions), 1)
        prediction = predictions[0]
        self.assertEqual(prediction.gene_symbol, "ERBB2")
        self.assertEqual(prediction.ensembl_gene_id, "ENSG00000141736")
        self.assertEqual(prediction.ensembl_protein_id, "ENSP00000275493")
        self.assertEqual(prediction.homology_type, "within_species_paralog")
        self.assertEqual(prediction.source_to_target_identity, 62.5)
        self.assertEqual(prediction.target_to_source_identity, 58.1)

    def test_fetch_human_paralogs_resolves_missing_symbols_concurrently_in_input_order(self) -> None:
        gene_ids = [f"ENSG{index}" for index in range(6)]
        four_started = threading.Event()
        lock = threading.Lock()
        active = 0
        max_active = 0

        class ConcurrentLookupHttp:
            def fetch_json(self, url: str, *, headers=None, cache_namespace="json", cache_key=None):
                nonlocal active, max_active
                if "/homology/symbol/human/TEST" in url:
                    return {
                        "data": [
                            {
                                "homologies": [
                                    {
                                        "id": gene_id,
                                        "type": "within_species_paralog",
                                        "target": {"id": gene_id, "protein_id": f"ENSP{index}"},
                                    }
                                    for index, gene_id in enumerate(gene_ids)
                                ]
                            }
                        ]
                    }
                for index, gene_id in enumerate(gene_ids):
                    if f"/lookup/id/{gene_id}" in url:
                        with lock:
                            active += 1
                            max_active = max(max_active, active)
                            if active == 4:
                                four_started.set()
                        try:
                            if not four_started.wait(timeout=2):
                                raise AssertionError("Missing-symbol lookups did not run concurrently.")
                            return {"display_name": f"GENE{index}"}
                        finally:
                            with lock:
                                active -= 1
                raise AssertionError(f"Unexpected URL: {url}")

        predictions = EnsemblParalogClient(ConcurrentLookupHttp()).fetch_human_paralogs("TEST")

        self.assertEqual([prediction.ensembl_gene_id for prediction in predictions], gene_ids)
        self.assertEqual([prediction.gene_symbol for prediction in predictions], [f"GENE{index}" for index in range(6)])
        self.assertEqual(
            [prediction.notes for prediction in predictions],
            [["Gene symbol resolved through Ensembl lookup."]] * 6,
        )
        self.assertEqual(max_active, 4)

    def test_fetch_human_paralogs_raises_first_lookup_error_in_input_order(self) -> None:
        class FailingLookupHttp:
            def fetch_json(self, url: str, *, headers=None, cache_namespace="json", cache_key=None):
                if "/homology/symbol/human/TEST" in url:
                    return {
                        "data": [
                            {
                                "homologies": [
                                    {"id": "ENSG_FIRST", "target": {"id": "ENSG_FIRST"}},
                                    {"id": "ENSG_SECOND", "target": {"id": "ENSG_SECOND"}},
                                ]
                            }
                        ]
                    }
                if "/lookup/id/ENSG_FIRST" in url:
                    raise RuntimeError("first lookup failed")
                if "/lookup/id/ENSG_SECOND" in url:
                    raise RuntimeError("second lookup failed")
                raise AssertionError(f"Unexpected URL: {url}")

        with self.assertRaisesRegex(RuntimeError, "first lookup failed"):
            EnsemblParalogClient(FailingLookupHttp()).fetch_human_paralogs("TEST")


class PrecomputedParalogStoreTests(unittest.TestCase):
    def test_load_paralog_context_reads_directional_identity_and_coverage(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            paralog_dir = root / "paralogs"
            paralog_dir.mkdir()
            payload_path = paralog_dir / "test_human.json"
            payload_path.write_text(
                json.dumps(
                    {
                        "target_entry_name": "TEST_HUMAN",
                        "target_accession": "P00001",
                        "target_gene_symbol": "TEST",
                        "canonical_family_accession": "IPR000001",
                        "canonical_family_name": "Test receptors",
                        "family_names": ["Test receptors", "HGNC helper family"],
                        "identity_matrix_labels": ["TEST_HUMAN", "PARA1_HUMAN"],
                        "identity_matrix": [[100.0, 75.0], [60.0, 100.0]],
                        "coverage_matrix_labels": ["TEST_HUMAN", "PARA1_HUMAN"],
                        "coverage_matrix": [[100.0, 92.0], [88.0, 100.0]],
                        "members": [
                            {
                                "gene_symbol": "TEST",
                                "gene_name": "Test receptor",
                                "accession": "P00001",
                                "entry_name": "TEST_HUMAN",
                                "sequence_length": 500,
                                "ectodomain_start": 25,
                                "ectodomain_end": 300,
                                "ectodomain_length": 276,
                                "sources": ["self"],
                                "notes": ["Target protein."],
                            },
                            {
                                "gene_symbol": "PARA1",
                                "gene_name": "Paralog one",
                                "accession": "Q00002",
                                "entry_name": "PARA1_HUMAN",
                                "sequence_length": 430,
                                "ectodomain_start": 40,
                                "ectodomain_end": 280,
                                "ectodomain_length": 241,
                                "sources": ["ensembl", "interpro_family"],
                                "notes": ["High-confidence paralog."],
                            },
                        ],
                        "metadata": {"candidate_source_counts": {"ensembl": 1, "interpro_family": 1, "hgnc": 0}},
                    }
                )
                + "\n",
                encoding="utf-8",
            )
            index_path = paralog_dir / "paralog_index.tsv"
            index_path.write_text(
                "input_index\tquery\ttarget_entry_name\ttarget_accession\ttarget_gene_symbol\tcanonical_family_accession\tcanonical_family_name\tmember_count\tstatus\terror\tjson_path\n"
                f"1\tTEST_HUMAN\tTEST_HUMAN\tP00001\tTEST\tIPR000001\tTest receptors\t2\tok\t\t{payload_path}\n",
                encoding="utf-8",
            )
            config = AnalysisConfig(precomputed_paralog_index_path=index_path)
            store = PrecomputedReferenceStore(config)
            target = TargetRecord(
                accession="P00001",
                entry_name="TEST_HUMAN",
                gene_symbol="TEST",
                protein_name="Test receptor",
                organism="Homo sapiens",
                taxon_id=9606,
                sequence="M" * 500,
            )
            context = store.load_paralog_context(target=target, max_members=5)
        self.assertIsNotNone(context)
        assert context is not None
        self.assertEqual(context.source, "Precomputed paralog reference set")
        self.assertEqual(context.family_names, ["Test receptors", "HGNC helper family"])
        self.assertEqual(context.identity_matrix[0][1], 75.0)
        self.assertEqual(context.coverage_matrix[1][0], 88.0)
        self.assertIn("Sources: ensembl, interpro_family.", context.members[1].notes[0])


class FakeUniProtClient:
    def __init__(self) -> None:
        self.entries = {
            "A_HUMAN": self._entry("P00001", "A_HUMAN", "GENEA", "MMMMMMMMMMAAAAAATTTTCC"),
            "B_HUMAN": self._entry("P00002", "B_HUMAN", "GENEB", "MMMMMMMMMMAAATTTTCC"),
        }

    def _entry(self, accession: str, entry_name: str, gene_symbol: str, sequence: str) -> dict:
        return {
            "primaryAccession": accession,
            "uniProtkbId": entry_name,
            "sequence": {"value": sequence},
            "organism": {"scientificName": "Homo sapiens", "taxonId": 9606},
            "genes": [{"geneName": {"value": gene_symbol}}],
            "proteinDescription": {"recommendedName": {"fullName": {"value": f"{gene_symbol} protein"}}},
        }

    def resolve_target(self, query: str, target_species_taxon: int = 9606) -> TargetResolution:
        entry = self.entries[query]
        return TargetResolution(
            target=TargetRecord(
                accession=entry["primaryAccession"],
                entry_name=entry["uniProtkbId"],
                gene_symbol=entry["genes"][0]["geneName"]["value"],
                protein_name=entry["proteinDescription"]["recommendedName"]["fullName"]["value"],
                organism="Homo sapiens",
                taxon_id=9606,
                sequence=entry["sequence"]["value"],
            ),
            entry=entry,
            notes=[],
            query_type="entry_name",
        )

    def get_features(self, entry: dict) -> list[Feature]:
        seq = entry["sequence"]["value"]
        ecto_end = len(seq) - 6
        return [
            Feature(type="SIGNAL", start=1, end=10, description="Signal peptide"),
            Feature(type="TOPO_DOM", start=11, end=ecto_end, description="Extracellular"),
            Feature(type="TRANSMEM", start=ecto_end + 1, end=ecto_end + 4, description="Helical"),
            Feature(type="TOPO_DOM", start=ecto_end + 5, end=len(seq), description="Cytoplasmic"),
        ]


class FakeInterProClient:
    def fetch_annotations(self, accession: str):
        return [
            DomainAnnotation(
                accession="IPR123456",
                name="Shared test family",
                type="family",
                source_database="INTERPRO",
                start=1,
                end=10,
            )
        ]


class FakePrecomputedStore:
    def find_ortholog_record(self, *, target):
        return None

    def canonical_family_from_record(self, record):
        return None


class FakeHGNCClient:
    def fetch_family_context(self, gene_symbol: str | None):
        return None


class FakeAnalyzer:
    def _apply_secreted_universe_fallback(self, **kwargs):
        return None

    def __init__(self) -> None:
        self.config = AnalysisConfig(max_paralog_context_members=10)
        self.uniprot_client = FakeUniProtClient()
        self.interpro_client = FakeInterProClient()
        self.precomputed_store = FakePrecomputedStore()
        self.hgnc_client = FakeHGNCClient()

    def _select_canonical_family(self, annotations, *, sequence_length: int, protein_name: str, gene_symbol: str | None):
        for annotation in annotations:
            if annotation.type == "family":
                return CanonicalFamily(
                    accession=annotation.accession,
                    name=annotation.name,
                    source_database=annotation.source_database,
                    fragment_count=1,
                    start=annotation.start,
                    end=annotation.end,
                )
        return None


class FakeEnsemblClient:
    def fetch_human_paralogs(self, gene_symbol: str | None):
        return []


class ParalogBuildTests(unittest.TestCase):
    def test_unique_label_assignment_supports_paralog_members_without_input_index(self) -> None:
        members = [
            ParalogMember(
                gene_symbol="GENEA",
                gene_name="Gene A",
                accession="P00001",
                entry_name="DUP_HUMAN",
                sequence_length=10,
                ectodomain_start=1,
                ectodomain_end=10,
                ectodomain_length=10,
                ectodomain_sequence="AAAAAAAAAA",
                matrix_label="DUP_HUMAN",
            ),
            ParalogMember(
                gene_symbol="GENEB",
                gene_name="Gene B",
                accession="P00002",
                entry_name="DUP_HUMAN",
                sequence_length=10,
                ectodomain_start=1,
                ectodomain_end=10,
                ectodomain_length=10,
                ectodomain_sequence="BBBBBBBBBB",
                matrix_label="DUP_HUMAN",
            ),
        ]

        _assign_unique_labels(members)

        self.assertEqual(members[0].matrix_label, "DUP_HUMAN__P00001")
        self.assertEqual(members[1].matrix_label, "DUP_HUMAN__P00002")

    def test_same_interpro_family_reuses_same_matrix_for_multiple_targets(self) -> None:
        analyzer = FakeAnalyzer()
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            tsv_path = root / "membrane.tsv"
            tsv_path.write_text(
                "\n".join(
                    [
                        "uniprot_id\tuniprot_name\tprot_family\tgene",
                        "P00001\tA_HUMAN\tTest family\tGENEA",
                        "P00002\tB_HUMAN\tTest family\tGENEB",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            _, _, summaries = build_paralog_reference_from_tsv(
                analyzer=analyzer,
                tsv_path=tsv_path,
                output_dir=root / "paralogs",
                verbose=False,
                ensembl_client=FakeEnsemblClient(),
            )
            payload_a = json.loads(Path(str(summaries[0]["json_path"])).read_text(encoding="utf-8"))
            payload_b = json.loads(Path(str(summaries[1]["json_path"])).read_text(encoding="utf-8"))
        self.assertEqual(payload_a["canonical_family_accession"], "IPR123456")
        self.assertEqual(payload_b["canonical_family_accession"], "IPR123456")
        self.assertEqual(payload_a["identity_matrix_labels"], payload_b["identity_matrix_labels"])
        self.assertEqual(payload_a["identity_matrix"], payload_b["identity_matrix"])
        self.assertEqual(payload_a["coverage_matrix"], payload_b["coverage_matrix"])


if __name__ == "__main__":
    unittest.main()
