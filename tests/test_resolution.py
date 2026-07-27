from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from agdesign2.exceptions import ResolutionError
from agdesign2.uniprot import UniProtClient


class FakeHttp:
    def __init__(self, mapping: dict[str, object]) -> None:
        self.mapping = mapping

    def fetch_json(self, url: str, **_: object) -> object:
        for key, value in self.mapping.items():
            if key in url:
                return value
        raise AssertionError(f"Unexpected URL: {url}")


ENTRY = {
    "primaryAccession": "P00533",
    "uniProtkbId": "EGFR_HUMAN",
    "sequence": {"value": "M" * 120},
    "organism": {"scientificName": "Homo sapiens", "taxonId": 9606},
    "genes": [{"geneName": {"value": "EGFR"}}],
    "proteinDescription": {"recommendedName": {"fullName": {"value": "EGFR"}}},
}


class UniProtResolutionTests(unittest.TestCase):
    def test_resolves_accession(self) -> None:
        client = UniProtClient(FakeHttp({"P00533.json": ENTRY}))
        resolution = client.resolve_target("P00533")
        self.assertEqual(resolution.target.entry_name, "EGFR_HUMAN")

    def test_resolves_entry_name(self) -> None:
        client = UniProtClient(FakeHttp({"EGFR_HUMAN.json": ENTRY}))
        resolution = client.resolve_target("EGFR_HUMAN")
        self.assertEqual(resolution.query_type, "entry_name")

    def test_resolves_gene_symbol(self) -> None:
        client = UniProtClient(
            FakeHttp(
                {
                    "search?query=%28gene_exact%3AEGFR": {"results": [{"primaryAccession": "P00533"}]},
                    "P00533.json": ENTRY,
                }
            )
        )
        resolution = client.resolve_target("EGFR")
        self.assertEqual(resolution.target.accession, "P00533")

    def test_uses_displayed_uniprot_isoform_as_canonical_sequence(self) -> None:
        client = UniProtClient(
            FakeHttp(
                {
                    "P08195.json": {
                        "primaryAccession": "P08195",
                        "uniProtkbId": "4F2_HUMAN",
                        "sequence": {"value": "C" * 529},
                        "organism": {"scientificName": "Homo sapiens", "taxonId": 9606},
                        "genes": [{"geneName": {"value": "SLC3A2"}}],
                        "proteinDescription": {"recommendedName": {"fullName": {"value": "4F2"}}},
                        "comments": [
                            {
                                "commentType": "ALTERNATIVE PRODUCTS",
                                "isoforms": [
                                    {"isoformIds": ["P08195-2"], "isoformSequenceStatus": "Displayed"},
                                    {"isoformIds": ["P08195-1"], "isoformSequenceStatus": "Described"},
                                ],
                            }
                        ],
                    }
                }
            )
        )
        resolution = client.resolve_target("P08195")
        self.assertEqual(resolution.target.canonical_isoform_id, "P08195-2")
        self.assertEqual(len(resolution.target.sequence), 529)

    def test_rejects_ambiguous_gene_symbol(self) -> None:
        client = UniProtClient(
            FakeHttp(
                {
                    "search?query=%28gene_exact%3AAKT1": {
                        "results": [{"primaryAccession": "P31749"}, {"primaryAccession": "Q9Y243"}]
                    }
                }
            )
        )
        with self.assertRaises(ResolutionError):
            client.resolve_target("AKT1")

    def test_gene_symbol_not_misread_as_accession(self) -> None:
        client = UniProtClient(
            FakeHttp(
                {
                    "search?query=%28gene_exact%3ANOTCH1": {"results": [{"primaryAccession": "Q9UKV3"}]},
                    "Q9UKV3.json": {
                        **ENTRY,
                        "primaryAccession": "Q9UKV3",
                        "uniProtkbId": "NOTC1_HUMAN",
                        "genes": [{"geneName": {"value": "NOTCH1"}}],
                    },
                }
            )
        )
        resolution = client.resolve_target("NOTCH1")
        self.assertEqual(resolution.query_type, "gene_symbol")
        self.assertEqual(resolution.target.gene_symbol, "NOTCH1")

    def test_normalizes_uniprot_feature_types(self) -> None:
        client = UniProtClient(FakeHttp({}))
        features = client.get_features(
            {
                "features": [
                    {
                        "type": "Signal",
                        "location": {"start": {"value": 1}, "end": {"value": 24}},
                    },
                    {
                        "type": "Topological domain",
                        "location": {"start": {"value": 25}, "end": {"value": 100}},
                        "description": "Extracellular",
                    },
                    {
                        "type": "Transmembrane",
                        "location": {"start": {"value": 101}, "end": {"value": 123}},
                    },
                ]
            }
        )
        self.assertEqual([feature.type for feature in features], ["SIGNAL", "TOPO_DOM", "TRANSMEM"])


if __name__ == "__main__":
    unittest.main()
