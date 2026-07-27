from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from agdesign2.config import AnalysisConfig
from agdesign2.family_alignments import build_family_alignments_from_tsv
from agdesign2.models import (
    AnalysisNote,
    CanonicalFamily,
    DomainAnnotation,
    Feature,
    Region,
    TargetRecord,
    TopologyAnnotation,
)
from agdesign2.uniprot import TargetResolution


class FakeUniProtClient:
    def __init__(self) -> None:
        self.entries = {
            "A_HUMAN": self._entry("P00001", "A_HUMAN", "GENEA", "MMMMMMMMMMAAAAAATTTTCC"),
            "B_HUMAN": self._entry("P00002", "B_HUMAN", "GENEB", "MMMMMMMMMMAAATTTTCC"),
            "C_HUMAN": self._entry("P00003", "C_HUMAN", "GENEC", "MMMMMMMMMMGGGGTTTTCC"),
            # Secreted (no transmembrane) targets that only get a design region
            # from the mature-chain fallback, like IL1B / FGF2 in production.
            "S1_HUMAN": self._entry("P10001", "S1_HUMAN", "SEC1", "QQQQWWWWEEEERRRRTTTT"),
            "S2_HUMAN": self._entry("P10002", "S2_HUMAN", "SEC2", "QQQQWWWWEEEERRRR"),
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
        if entry["uniProtkbId"].startswith("S"):
            # Secreted: a mature CHAIN but no transmembrane/topology evidence, so
            # derive_ectodomain alone yields no ectodomain.
            return [Feature(type="CHAIN", start=1, end=len(seq), description="Mature chain")]
        ecto_end = len(seq) - 6
        return [
            Feature(type="SIGNAL", start=1, end=10, description="Signal peptide"),
            Feature(type="TOPO_DOM", start=11, end=ecto_end, description="Extracellular"),
            Feature(type="TRANSMEM", start=ecto_end + 1, end=ecto_end + 4, description="Helical"),
            Feature(type="TOPO_DOM", start=ecto_end + 5, end=len(seq), description="Cytoplasmic"),
        ]


class FakeInterProClient:
    def fetch_annotations(self, accession: str):
        if accession == "P00003":
            return [
                DomainAnnotation(
                    accession="PF00001",
                    name="PF-only family",
                    type="family",
                    source_database="PFAM",
                    start=1,
                    end=10,
                )
            ]
        if accession in {"P10001", "P10002"}:
            return [
                DomainAnnotation(
                    accession="IPR777777",
                    name="Secreted test family",
                    type="family",
                    source_database="INTERPRO",
                    start=1,
                    end=10,
                )
            ]
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


class FakeAnalyzer:
    def __init__(self) -> None:
        self.config = AnalysisConfig()
        self.uniprot_client = FakeUniProtClient()
        self.interpro_client = FakeInterProClient()

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

    def _apply_secreted_universe_fallback(self, *, target, topology, features):
        # Faithful to AntigenAnalyzer: no-op when an ectodomain already exists;
        # otherwise derive the design region from the mature CHAIN feature.
        if topology.ectodomain is not None:
            return None
        chains = [f for f in features if f.type.upper() in {"CHAIN", "PEPTIDE"}]
        if not chains:
            return None
        chain = max(chains, key=lambda f: f.end - f.start + 1)
        region = Region(start=chain.start, end=chain.end, label="Secreted chain", source="target-universe", confidence=0.75)
        topology.ectodomain = region
        topology.topology = TopologyAnnotation(
            topology_class="secreted_external_evidence",
            extracellular_regions=[region],
            major_extracellular_region=region,
        )
        return AnalysisNote(severity="warning", message="secreted fallback", source="target-universe")


class FamilyAlignmentTests(unittest.TestCase):
    def test_builds_interpro_only_family_alignments_with_directional_matrix(self) -> None:
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
                        "P00003\tC_HUMAN\tPF-only\tGENEC",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            summary_tsv, summary_json, summaries = build_family_alignments_from_tsv(
                analyzer=analyzer,
                tsv_path=tsv_path,
                output_dir=root / "family_alignments",
                verbose=False,
            )

            self.assertTrue(summary_tsv.exists())
            self.assertTrue(summary_json.exists())
            self.assertEqual(len(summaries), 1)
            self.assertEqual(summaries[0]["family_accession"], "IPR123456")
            self.assertEqual(summaries[0]["member_count"], 2)

            family_json = Path(str(summaries[0]["json_path"]))
            payload = json.loads(family_json.read_text(encoding="utf-8"))
            self.assertEqual(payload["family_accession"], "IPR123456")
            self.assertEqual(payload["identity_matrix_labels"], ["A_HUMAN", "B_HUMAN"])
            self.assertEqual(payload["identity_matrix"][0][1], 50.0)
            self.assertEqual(payload["identity_matrix"][1][0], 100.0)
            self.assertEqual(payload["members"][0]["ectodomain_sequence"], "AAAAAA")
            self.assertEqual(payload["members"][1]["ectodomain_sequence"], "AAA")
            self.assertEqual(len(payload["pairwise_alignments"]), 1)
            self.assertEqual(payload["pairwise_alignments"][0]["query_to_subject_identity"], 50.0)
            self.assertEqual(payload["pairwise_alignments"][0]["subject_to_query_identity"], 100.0)

    def test_secreted_no_tm_targets_included_via_chain_fallback(self) -> None:
        # Regression for the 22 secreted targets (IL1B, FGF2, immunoglobulins...)
        # that derive_ectodomain alone drops, so they vanished from their own
        # family context. _build_family_member must apply the same secreted
        # mature-chain fallback the report pipeline uses.
        analyzer = FakeAnalyzer()
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            tsv_path = root / "secreted.tsv"
            tsv_path.write_text(
                "\n".join(
                    [
                        "uniprot_id\tuniprot_name\tprot_family\tgene",
                        "P10001\tS1_HUMAN\tSecreted family\tSEC1",
                        "P10002\tS2_HUMAN\tSecreted family\tSEC2",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            _, _, summaries = build_family_alignments_from_tsv(
                analyzer=analyzer,
                tsv_path=tsv_path,
                output_dir=root / "family_alignments",
                verbose=False,
            )
            self.assertEqual(len(summaries), 1)
            self.assertEqual(summaries[0]["family_accession"], "IPR777777")
            self.assertEqual(summaries[0]["member_count"], 2)
            payload = json.loads(Path(str(summaries[0]["json_path"])).read_text(encoding="utf-8"))
            # Both secreted targets present, with the full mature chain as design region.
            self.assertEqual(payload["identity_matrix_labels"], ["S1_HUMAN", "S2_HUMAN"])
            self.assertEqual(payload["members"][0]["ectodomain_sequence"], "QQQQWWWWEEEERRRRTTTT")
            self.assertEqual(payload["members"][1]["ectodomain_sequence"], "QQQQWWWWEEEERRRR")


if __name__ == "__main__":
    unittest.main()
