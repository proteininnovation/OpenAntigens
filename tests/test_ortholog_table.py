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

from unittest import mock

import agdesign2.ortholog_table as ortholog_table_module
from agdesign2.hgnc import OrthologLookup, OrthologPrediction
from agdesign2.pipeline import AntigenAnalyzer
from agdesign2.models import AnalysisNote, DomainAnnotation, TargetRecord
from agdesign2.ortholog_proteomes import ResolvedOrthologProtein
from agdesign2.ortholog_table import build_ortholog_table_from_tsv, load_ortholog_rows
from agdesign2.refseq import RefSeqProtein
from agdesign2.uniprot import TargetResolution


class FakeOrthologResolver:
    """Stubs the local proteome resolver: mouse via UniProt, macaque via RefSeq."""

    def __init__(self, **kwargs) -> None:
        pass

    def prewarm(self, species) -> None:
        pass

    def resolve(self, species, *, gene_symbol=None, gene_id=None, fallback_symbol=None):
        if species == "mouse":
            return ResolvedOrthologProtein(accession="Q01279", sequence="MAAACCG", source="uniprot", organism="Mus musculus")
        if species == "macaca_fascicularis":
            return ResolvedOrthologProtein(accession="XP_005549616", sequence="MAAACCC", source="refseq", organism="Macaca fascicularis")
        return None


class FakeUniProtClient:
    def __init__(self) -> None:
        self.calls: list[str] = []

    def resolve_target(self, query: str, target_species_taxon: int = 9606) -> TargetResolution:
        self.calls.append(query)
        if query == "FAIL_HUMAN":
            raise RuntimeError("simulated UniProt resolution failure")
        gene_symbol = query.split("_", 1)[0]
        return TargetResolution(
            target=TargetRecord(
                accession="P00533",
                entry_name=query,
                gene_symbol=gene_symbol,
                protein_name=f"{gene_symbol} protein",
                organism="Homo sapiens",
                taxon_id=9606,
                sequence="MAAACCC",
            ),
            entry={},
            notes=[AnalysisNote(severity="info", message=f"Resolved from query {query}.", source="UniProt")],
            query_type="uniprot_name",
        )


class FakeHGNCClient:
    def fetch_hcop_orthologs(self, gene_symbol: str | None) -> OrthologLookup | None:
        if gene_symbol == "NOHCOP":
            raise RuntimeError("simulated HCOP failure")
        return OrthologLookup(
            human_gene_symbol=gene_symbol or "EGFR",
            human_ncbi_gene_id="1956",
            predictions={
                "mouse": OrthologPrediction(
                    species="mouse",
                    gene_symbol=(gene_symbol or "EGFR").capitalize(),
                    gene_name=f"{gene_symbol or 'EGFR'} mouse ortholog",
                    ncbi_gene_id="13649",
                    ensembl_id="ENSMUSG00000020122",
                    mod_id="MGI:95294",
                    gene_symbol_source="MGI",
                    evidence_count=8,
                ),
                "macaca_fascicularis": OrthologPrediction(
                    species="macaca_fascicularis",
                    gene_symbol=gene_symbol or "EGFR",
                    gene_name=f"{gene_symbol or 'EGFR'} macaque ortholog",
                    ncbi_gene_id="613027",
                    ensembl_id="ENSMMUG00000022394",
                    mod_id=None,
                    gene_symbol_source="VGNC",
                    evidence_count=5,
                ),
            },
        )


class FakeRefSeqClient:
    def fetch_canonical_protein(
        self,
        *,
        gene_symbol: str,
        organism: str,
        gene_id: str | int | None = None,
    ) -> RefSeqProtein | None:
        if gene_symbol == "NOREFSEQ":
            raise RuntimeError("simulated RefSeq failure")
        records = {
            "Homo sapiens": ("NP_005219", "MAAACCC"),
            "Mus musculus": ("NP_031989", "MAAACCG"),
            "Macaca fascicularis": ("XP_005549616", "MAAACCC"),
        }
        accession, sequence = records[organism]
        return RefSeqProtein(
            target=TargetRecord(
                accession=accession,
                entry_name=accession,
                gene_symbol=gene_symbol,
                protein_name=gene_symbol,
                organism=organism,
                taxon_id={"Homo sapiens": 9606, "Mus musculus": 10090, "Macaca fascicularis": 9541}[organism],
                sequence=sequence,
            ),
            notes=[AnalysisNote(severity="info", message=f"Resolved from RefSeq canonical protein candidate {accession}.", source="RefSeq")],
            title=gene_symbol,
        )


class FakeInterProClient:
    def fetch_annotations(self, accession: str):
        return [
            DomainAnnotation(
                accession="IPR050122",
                name="Receptor Tyrosine Kinase",
                type="family",
                source_database="INTERPRO",
                start=58,
                end=974,
            )
        ]


class _FakeConfig:
    ortholog_fasta_dir = Path("/tmp")
    species_taxonomy = {"human": 9606, "mouse": 10090, "macaca_fascicularis": 9541}


class FakeAnalyzer:
    def __init__(self) -> None:
        self.uniprot_client = FakeUniProtClient()
        self.hgnc_client = FakeHGNCClient()
        self.refseq_client = FakeRefSeqClient()
        self.refseq_client.http = None
        self.interpro_client = FakeInterProClient()
        self.config = _FakeConfig()

    def _select_canonical_family(self, annotations, *, sequence_length: int, protein_name: str, gene_symbol: str | None):
        analyzer = AntigenAnalyzer()
        return analyzer._select_canonical_family(
            annotations,
            sequence_length=sequence_length,
            protein_name=protein_name,
            gene_symbol=gene_symbol,
        )


class OrthologTableTests(unittest.TestCase):
    def setUp(self) -> None:
        patcher = mock.patch.object(ortholog_table_module, "OrthologProteomeResolver", FakeOrthologResolver)
        patcher.start()
        self.addCleanup(patcher.stop)

    def test_interrupted_rows_resume_without_rewriting_prefixes(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source = root / "genes.tsv"
            source.write_text("uniprot_name\tgene\nEGFR_HUMAN\tEGFR\nERBB2_HUMAN\tERBB2\n")
            output = root / "orthologs.tsv"
            original = ortholog_table_module._build_table_row

            def interrupted(**kwargs):
                if kwargs["input_index"] == 2:
                    raise KeyboardInterrupt()
                return original(**kwargs)

            with mock.patch.object(ortholog_table_module, "_build_table_row", side_effect=interrupted):
                with self.assertRaises(KeyboardInterrupt):
                    build_ortholog_table_from_tsv(analyzer=FakeAnalyzer(), tsv_path=source, output_path=output)
            self.assertTrue((root / "orthologs.checkpoints/1.json").is_file())
            analyzer = FakeAnalyzer()
            with mock.patch.object(ortholog_table_module, "_write_outputs", wraps=ortholog_table_module._write_outputs) as writer:
                _, _, rows = build_ortholog_table_from_tsv(analyzer=analyzer, tsv_path=source, output_path=output)
            self.assertEqual(analyzer.uniprot_client.calls, ["ERBB2_HUMAN"])
            self.assertEqual(len(rows), 2)
            self.assertEqual(writer.call_count, 1)
            self.assertEqual(list((root / "orthologs.checkpoints").glob("*.json")), [])
            with mock.patch.object(ortholog_table_module, "_write_outputs", wraps=ortholog_table_module._write_outputs) as writer:
                build_ortholog_table_from_tsv(analyzer=FakeAnalyzer(), tsv_path=source, output_path=output, jobs=2)
            self.assertEqual(writer.call_count, 1)

    def test_load_ortholog_rows_accepts_canonical_csv_columns(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "targets.csv"
            path.write_text(
                "accession,uniprot_entry_name,gene_symbol\n"
                "P00533,EGFR_HUMAN,EGFR\n",
                encoding="utf-8",
            )
            rows = load_ortholog_rows(path)

        self.assertEqual([(row.query, row.source_column) for row in rows], [("EGFR_HUMAN", "uniprot_name")])
        self.assertEqual(rows[0].record["uniprot_id"], "P00533")
        self.assertEqual(rows[0].record["gene"], "EGFR")

    def test_builds_reference_table_from_tsv(self) -> None:
        analyzer = FakeAnalyzer()
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            tsv_path = root / "genes.tsv"
            tsv_path.write_text(
                "uniprot_id\tuniprot_name\tprot_family\tgene\nP00533\tEGFR_HUMAN\tRTK\tEGFR\n",
                encoding="utf-8",
            )
            output_path = root / "ortholog_reference_table.tsv"

            tsv_file, json_file, rows = build_ortholog_table_from_tsv(
                analyzer=analyzer,
                tsv_path=tsv_path,
                output_path=output_path,
                verbose=False,
            )

            self.assertEqual(tsv_file, output_path)
            self.assertTrue(tsv_file.exists())
            self.assertTrue(json_file.exists())
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0].input_index, 1)
            self.assertEqual(rows[0].input_uniprot_id, "P00533")
            self.assertEqual(rows[0].input_uniprot_name, "EGFR_HUMAN")
            self.assertEqual(rows[0].input_gene, "EGFR")
            self.assertEqual(rows[0].input_prot_family, "RTK")
            self.assertEqual(rows[0].query, "EGFR_HUMAN")
            self.assertEqual(rows[0].source_column, "uniprot_name")
            self.assertEqual(rows[0].status, "ok")
            self.assertIsNone(rows[0].error)
            self.assertEqual(rows[0].human_gene_symbol, "EGFR")
            self.assertEqual(rows[0].mouse_gene_symbol, "Egfr")
            self.assertEqual(rows[0].human_refseq_accession, "P00533")  # human reference = UniProt target
            self.assertEqual(rows[0].mouse_refseq_accession, "Q01279")  # mouse = UniProt SwissProt
            self.assertEqual(rows[0].macaca_fascicularis_refseq_accession, "XP_005549616")  # macaque = RefSeq
            self.assertEqual(rows[0].canonical_family_accession, "IPR050122")
            self.assertEqual(rows[0].canonical_family_name, "Receptor Tyrosine Kinase")
            self.assertAlmostEqual(rows[0].mouse_identity_to_human or 0.0, 85.71, places=2)
            self.assertAlmostEqual(rows[0].macaca_fascicularis_identity_to_human or 0.0, 100.0, places=2)

            written = json.loads(json_file.read_text(encoding="utf-8"))
            self.assertEqual(written[0]["human_uniprot_accession"], "P00533")
            self.assertEqual(written[0]["input_uniprot_name"], "EGFR_HUMAN")

    def test_emits_error_row_when_target_resolution_fails(self) -> None:
        analyzer = FakeAnalyzer()
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            tsv_path = root / "genes.tsv"
            tsv_path.write_text(
                "uniprot_id\tuniprot_name\tprot_family\tgene\nP00000\tFAIL_HUMAN\tUnknown\tFAIL\n",
                encoding="utf-8",
            )

            _, _, rows = build_ortholog_table_from_tsv(
                analyzer=analyzer,
                tsv_path=tsv_path,
                output_path=root / "ortholog_reference_table.tsv",
            )

            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0].status, "error")
            self.assertIn("simulated UniProt resolution failure", rows[0].error or "")
            self.assertEqual(rows[0].input_gene, "FAIL")
            self.assertEqual(rows[0].human_gene_symbol, "FAIL")
            self.assertIsNone(rows[0].human_uniprot_accession)

    def test_resume_reuses_existing_rows_and_normalizes_partial_output(self) -> None:
        analyzer = FakeAnalyzer()
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            tsv_path = root / "genes.tsv"
            tsv_path.write_text(
                "\n".join(
                    [
                        "uniprot_id\tuniprot_name\tprot_family\tgene",
                        "P00533\tEGFR_HUMAN\tRTK\tEGFR",
                        "P15056\tERBB2_HUMAN\tRTK\tERBB2",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            output_path = root / "ortholog_reference_table.tsv"
            output_path.with_suffix(".json").write_text(
                json.dumps(
                    [
                        {
                            "human_gene_symbol": "EGFR",
                            "mouse_gene_symbol": "Egfr",
                            "macaca_fascicularis_gene_symbol": "EGFR",
                            "canonical_family_accession": "IPR050122",
                            "canonical_family_name": "Receptor Tyrosine Kinase",
                            "human_refseq_accession": "NP_005219",
                            "mouse_refseq_accession": "NP_031989",
                            "macaca_fascicularis_refseq_accession": "XP_005549616",
                            "human_uniprot_accession": "P00533",
                            "human_refseq_sequence": "MAAACCC",
                            "mouse_refseq_sequence": "MAAACCG",
                            "macaca_fascicularis_refseq_sequence": "MAAACCC",
                            "mouse_identity_to_human": 85.71,
                            "macaca_fascicularis_identity_to_human": 100.0,
                            "human_notes": "old row",
                            "mouse_notes": "",
                            "macaca_fascicularis_notes": "",
                        }
                    ],
                    indent=2,
                )
                + "\n",
                encoding="utf-8",
            )

            _, json_file, rows = build_ortholog_table_from_tsv(
                analyzer=analyzer,
                tsv_path=tsv_path,
                output_path=output_path,
                resume=True,
            )

            self.assertEqual(len(rows), 2)
            self.assertEqual(analyzer.uniprot_client.calls, ["ERBB2_HUMAN"])
            self.assertEqual(rows[0].input_index, 1)
            self.assertEqual(rows[0].input_uniprot_name, "EGFR_HUMAN")
            self.assertEqual(rows[0].status, "ok")
            self.assertEqual(rows[1].input_index, 2)
            self.assertEqual(rows[1].input_uniprot_name, "ERBB2_HUMAN")
            written = json.loads(json_file.read_text(encoding="utf-8"))
            self.assertEqual(len(written), 2)
            self.assertEqual(written[0]["input_index"], 1)
            self.assertEqual(written[1]["input_index"], 2)

    def test_preserves_input_order_and_duplicates(self) -> None:
        analyzer = FakeAnalyzer()
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            tsv_path = root / "genes.tsv"
            tsv_path.write_text(
                "\n".join(
                    [
                        "uniprot_id\tuniprot_name\tprot_family\tgene",
                        "P00533\tEGFR_HUMAN\tRTK\tEGFR",
                        "P00533\tEGFR_HUMAN\tRTK\tEGFR",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            _, _, rows = build_ortholog_table_from_tsv(
                analyzer=analyzer,
                tsv_path=tsv_path,
                output_path=root / "ortholog_reference_table.tsv",
            )

            self.assertEqual(len(rows), 2)
            self.assertEqual([row.input_index for row in rows], [1, 2])
            self.assertEqual([row.query for row in rows], ["EGFR_HUMAN", "EGFR_HUMAN"])
            self.assertEqual([row.input_uniprot_id for row in rows], ["P00533", "P00533"])

    def test_resume_skips_legacy_rows_that_do_not_match_current_input(self) -> None:
        analyzer = FakeAnalyzer()
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            tsv_path = root / "genes.tsv"
            tsv_path.write_text(
                "uniprot_id\tuniprot_name\tprot_family\tgene\nP00533\tEGFR_HUMAN\tRTK\tEGFR\n",
                encoding="utf-8",
            )
            output_path = root / "ortholog_reference_table.tsv"
            output_path.with_suffix(".json").write_text(
                json.dumps([{"human_gene_symbol": "ZPLD1"}], indent=2) + "\n",
                encoding="utf-8",
            )

            _, _, rows = build_ortholog_table_from_tsv(
                analyzer=analyzer,
                tsv_path=tsv_path,
                output_path=output_path,
                resume=True,
            )

            self.assertEqual(analyzer.uniprot_client.calls, ["EGFR_HUMAN"])
            self.assertEqual(rows[0].human_gene_symbol, "EGFR")


if __name__ == "__main__":
    unittest.main()
