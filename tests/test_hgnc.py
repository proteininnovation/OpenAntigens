from __future__ import annotations

import sys
import os
from unittest import mock
import json
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from agdesign2.config import AnalysisConfig
from agdesign2.exceptions import ExternalServiceError
from agdesign2.hgnc import HGNCClient, INTEGRIN_ALPHA_GENES
from agdesign2.models import FamilyMember


def build_entry(entry_name: str, accession: str, sequence: str, gene: str) -> dict:
    return {
        "primaryAccession": accession,
        "uniProtkbId": entry_name,
        "sequence": {"value": sequence},
        "organism": {"scientificName": "Homo sapiens", "taxonId": 9606},
        "genes": [{"geneName": {"value": gene}}],
        "proteinDescription": {"recommendedName": {"fullName": {"value": gene}}},
        "features": [
            {"type": "Signal", "location": {"start": {"value": 1}, "end": {"value": 30}}, "description": "Signal peptide"},
            {"type": "Topological domain", "location": {"start": {"value": 31}, "end": {"value": 300}}, "description": "Extracellular"},
            {"type": "Transmembrane", "location": {"start": {"value": 301}, "end": {"value": 320}}, "description": "Helical"},
            {"type": "Topological domain", "location": {"start": {"value": 321}, "end": {"value": 350}}, "description": "Cytoplasmic"},
        ],
    }


class FakeHttpClient:
    def fetch_json(self, url: str, *, headers=None, cache_namespace: str = "json"):
        if url.endswith("/fetch/symbol/ITGAV"):
            return {
                "response": {
                    "docs": [
                        {
                            "symbol": "ITGAV",
                            "hgnc_id": "HGNC:6141",
                            "gene_group": ["Integrin alpha chains"],
                            "gene_group_id": [1160],
                        }
                    ]
                }
            }
        if url.endswith("/fetch/gene_group_id/1160"):
            raise ExternalServiceError("simulated HGNC family endpoint outage")
        raise AssertionError(f"Unexpected URL: {url}")


class FailingSymbolHttpClient:
    def fetch_json(self, url: str, *, headers=None, cache_namespace: str = "json"):
        if url.endswith("/fetch/symbol/ITGAV"):
            raise ExternalServiceError("simulated HGNC symbol endpoint outage")
        raise AssertionError(f"Unexpected URL: {url}")


class MultiFamilyHttpClient:
    def fetch_json(self, url: str, *, headers=None, cache_namespace: str = "json"):
        if url.endswith("/fetch/symbol/ITGAV"):
            return {
                "response": {
                    "docs": [
                        {
                            "symbol": "ITGAV",
                            "hgnc_id": "HGNC:6141",
                            "gene_group": ["CD molecules", "Integrin alpha subunits"],
                            "gene_group_id": [471, 1160],
                        }
                    ]
                }
            }
        if url.endswith("/fetch/gene_group_id/471"):
            return {
                "response": {
                    "docs": [{"symbol": f"CDX{i}", "name": f"CD molecule {i}"} for i in range(1, 6)]
                }
            }
        if url.endswith("/fetch/gene_group_id/1160"):
            return {
                "response": {
                    "docs": [
                        {"symbol": "ITGAV", "name": "integrin subunit alpha V"},
                        {"symbol": "ITGA5", "name": "integrin subunit alpha 5"},
                    ]
                }
            }
        raise AssertionError(f"Unexpected URL: {url}")


class LargeFamilyHttpClient:
    def fetch_json(self, url: str, *, headers=None, cache_namespace: str = "json"):
        if url.endswith("/fetch/symbol/BIGFAM"):
            return {
                "response": {
                    "docs": [
                        {
                            "symbol": "BIGFAM",
                            "hgnc_id": "HGNC:99999",
                            "gene_group": ["Very large family"],
                            "gene_group_id": [4242],
                        }
                    ]
                }
            }
        if url.endswith("/fetch/gene_group_id/4242"):
            return {
                "response": {
                    "docs": [{"symbol": f"BIG{i}", "name": f"Big family member {i}"} for i in range(1, 23)]
                }
            }
        if "/cgi-bin/orthologs/hcop?" in url and "q=EGFR" in url:
            return {
                "gene_symbol": "EGFR",
                "ncbi_gene_id": "1956",
                "predictions": {
                    "mouse": [
                        {
                            "gene_symbol": "Egfr",
                            "gene_name": "epidermal growth factor receptor",
                            "mod_id": "MGI:95294",
                            "ncbi_gene_id": "13649",
                            "ensembl_id": "ENSMUSG00000020122",
                            "gene_symbol_source": "MGI",
                            "evidence": [{"source": "NCBI"}, {"source": "Ensembl"}],
                        },
                        {
                            "gene_symbol": "Egfr2",
                            "gene_name": "other candidate",
                            "mod_id": "MGI:00000",
                            "ncbi_gene_id": "99999",
                            "ensembl_id": "ENSMUSG99999999999",
                            "gene_symbol_source": "NCBI",
                            "evidence": [{"source": "NCBI"}],
                        },
                    ],
                    "macaque": [
                        {
                            "gene_symbol": "EGFR",
                            "gene_name": "epidermal growth factor receptor",
                            "mod_id": None,
                            "ncbi_gene_id": "613027",
                            "ensembl_id": "ENSMMUG00000022394",
                            "gene_symbol_source": "VGNC",
                            "evidence": [{"source": "NCBI"}, {"source": "OMA"}, {"source": "Ensembl"}],
                        }
                    ],
                },
            }
        raise AssertionError(f"Unexpected URL: {url}")


class FakeUniProtClient:
    def __init__(self) -> None:
        self.entries = {
            "P06756": build_entry("ITAV_HUMAN", "P06756", "M" * 350, "ITGAV"),
            "P08648": build_entry("ITA5_HUMAN", "P08648", "M" * 360, "ITGA5"),
        }
        self.by_symbol = {
            "ITGAV": "P06756",
            "ITGA5": "P08648",
        }

    def search(self, query: str, *, size: int = 10):
        for symbol, accession in self.by_symbol.items():
            if f"gene_exact:{symbol}" in query:
                return {"results": [{"primaryAccession": accession}]}
        return {"results": []}

    def _fetch_entry(self, accession: str):
        return self.entries[accession]

    def fetch_entry_by_name(self, entry_name: str):
        for entry in self.entries.values():
            if entry["uniProtkbId"] == entry_name:
                return entry
        return None

    def get_features(self, entry: dict):
        from agdesign2.uniprot import UniProtClient

        return UniProtClient(http=None).get_features(entry)  # type: ignore[arg-type]


class HGNCTests(unittest.TestCase):
    def test_falls_back_to_integrin_alpha_family_when_hgnc_family_members_fail(self) -> None:
        client = HGNCClient(FakeHttpClient(), FakeUniProtClient(), AnalysisConfig())
        family = client.fetch_family_context("ITGAV")
        self.assertIsNotNone(family)
        assert family is not None
        self.assertEqual(family.source, "HGNC fallback")
        self.assertIn("Integrin alpha chains", family.family_names)
        self.assertEqual(family.metadata.get("fallback"), True)
        self.assertEqual(len(family.members), len(INTEGRIN_ALPHA_GENES))
        member_symbols = {member.gene_symbol for member in family.members}
        self.assertIn("ITGAV", member_symbols)
        self.assertIn("ITGA5", member_symbols)

    def test_falls_back_to_integrin_alpha_family_when_hgnc_symbol_lookup_fails(self) -> None:
        client = HGNCClient(FailingSymbolHttpClient(), FakeUniProtClient(), AnalysisConfig())
        family = client.fetch_family_context("ITGAV")
        self.assertIsNotNone(family)
        assert family is not None
        self.assertEqual(family.source, "HGNC fallback")
        self.assertIn("Integrin alpha chains", family.family_names)
        self.assertEqual(len(family.members), len(INTEGRIN_ALPHA_GENES))
        self.assertTrue(all(member.accession is None for member in family.members))

    def test_prefers_smallest_hgnc_family_when_multiple_groups_are_present(self) -> None:
        client = HGNCClient(MultiFamilyHttpClient(), FakeUniProtClient(), AnalysisConfig())
        family = client.fetch_family_context("ITGAV")
        self.assertIsNotNone(family)
        assert family is not None
        self.assertEqual(family.source, "HGNC")
        self.assertEqual(family.family_names, ["Integrin alpha subunits"])
        self.assertEqual(family.metadata.get("selected_gene_group_id"), 1160)
        self.assertEqual([member.gene_symbol for member in family.members], ["ITGA5", "ITGAV"])

    def test_ignores_family_context_when_family_is_too_large(self) -> None:
        client = HGNCClient(LargeFamilyHttpClient(), FakeUniProtClient(), AnalysisConfig(max_family_context_members=20))
        family = client.fetch_family_context("BIGFAM")
        self.assertIsNone(family)

    def test_fetches_hcop_orthologs_for_mouse_and_macaca(self) -> None:
        client = HGNCClient(LargeFamilyHttpClient(), FakeUniProtClient(), AnalysisConfig())
        lookup = client.fetch_hcop_orthologs("EGFR")
        self.assertIsNotNone(lookup)
        assert lookup is not None
        self.assertEqual(lookup.human_ncbi_gene_id, "1956")
        self.assertEqual(lookup.predictions["mouse"].gene_symbol, "Egfr")
        self.assertEqual(lookup.predictions["mouse"].ncbi_gene_id, "13649")
        self.assertEqual(lookup.predictions["macaca_fascicularis"].gene_symbol, "EGFR")
        self.assertEqual(lookup.predictions["macaca_fascicularis"].ncbi_gene_id, "613027")

    @mock.patch.dict(os.environ, {"AGDESIGN2_ALIGNMENT_BACKEND": "python"})
    def test_long_sequence_identity_uses_global_alignment_in_both_orders(self) -> None:
        fixture = json.loads((Path(__file__).parent / "fixtures/ptprc_mapping.json").read_text())
        query, subject = fixture["human_sequence"], fixture["mouse_sequence"]
        client = HGNCClient(FakeHttpClient(), FakeUniProtClient(), AnalysisConfig())
        forward = client._pairwise_directional_identity(query, subject)
        self.assertEqual(forward, (42.57, 43.44))
        client._identity_cache.clear()
        reverse = client._pairwise_directional_identity(subject, query)
        self.assertEqual(reverse, forward[::-1])

    def test_identity_matrix_is_directional_for_contained_ectodomains(self) -> None:
        client = HGNCClient(FakeHttpClient(), FakeUniProtClient(), AnalysisConfig())
        client.uniprot_client.entries["P11111"] = build_entry("LONG_HUMAN", "P11111", "M" * 30 + "A" * 260 + "M" * 20, "LONG")
        client.uniprot_client.entries["P22222"] = build_entry("SHORT_HUMAN", "P22222", "M" * 30 + "A" * 200 + "M" * 20, "SHORT")
        matrix = client._build_identity_matrix(
            [
                FamilyMember(
                    gene_symbol="LONG",
                    gene_name="long protein",
                    accession="P11111",
                    entry_name="LONG_HUMAN",
                    ectodomain_start=31,
                    ectodomain_end=290,
                    ectodomain_length=260,
                ),
                FamilyMember(
                    gene_symbol="SHORT",
                    gene_name="short protein",
                    accession="P22222",
                    entry_name="SHORT_HUMAN",
                    ectodomain_start=31,
                    ectodomain_end=230,
                    ectodomain_length=200,
                ),
            ]
        )
        self.assertEqual(matrix[0][1], 76.92)
        self.assertEqual(matrix[1][0], 100.0)


if __name__ == "__main__":
    unittest.main()
