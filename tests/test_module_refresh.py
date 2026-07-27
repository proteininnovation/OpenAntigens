from __future__ import annotations

import json
import sys
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from agdesign2.models import AnalysisNote, AssemblyRequirement, BlastHit, ComplexPortalComplex, ComplexPortalParticipant
from agdesign2.module_refresh import refresh_report_modules


class FakeBlastClient:
    def __init__(self) -> None:
        self.queries: list[str] = []

    def search(self, query_sequence: str, *, target_accession: str, target_entry_name: str):
        self.queries.append(query_sequence)
        return [
            BlastHit(
                subject_id="sp|Q01279|EGFR_MOUSE",
                description="EGFR mouse",
                species="Mus musculus",
                identity=88.7,
                coverage=100.0,
                alignment_length=len(query_sequence),
                evalue=1e-30,
                bitscore=200.0,
                query_start=1,
                query_end=len(query_sequence),
                subject_start=1,
                subject_end=len(query_sequence),
                query_alignment=query_sequence,
                subject_alignment=query_sequence,
            )
        ]


class FakeAnalyzer:
    def __init__(self) -> None:
        self.blast_client = FakeBlastClient()

    def _canonicalize_macaca_cross_reactivity_hits(self, query_sequence: str, hits):
        return hits

    def _annotate_blast_hits_with_surface_identity(
        self,
        hits,
        *,
        query_region_start: int,
        extracellular_surface_positions: set[int],
    ):
        annotated = []
        for hit in hits:
            query_position = query_region_start - 1
            classes = []
            for residue in hit.query_alignment or "":
                if residue != "-":
                    query_position += 1
                classes.append("S" if residue != "-" and query_position in extracellular_surface_positions else ".")
            annotated.append(
                replace(
                    hit,
                    query_alignment_annotation="".join(classes) if classes else None,
                    extracellular_surface_aligned_positions=sum(1 for item in classes if item == "S"),
                )
            )
        return annotated


class FakeBulkBlastClient(FakeBlastClient):
    def __init__(self) -> None:
        super().__init__()
        self.batch_queries: list[dict[str, str]] = []

    def search_many(self, queries: list[dict[str, str]]):
        self.batch_queries.extend(queries)
        return {
            item["id"]: [
                BlastHit(
                    subject_id=f"sp|Q01279|{item['id'].upper()}_MOUSE",
                    description="bulk mouse hit",
                    species="Mus musculus",
                    identity=88.7,
                    coverage=100.0,
                    alignment_length=len(item["sequence"]),
                    evalue=1e-30,
                    bitscore=200.0,
                    query_start=1,
                    query_end=len(item["sequence"]),
                    subject_start=1,
                    subject_end=len(item["sequence"]),
                    query_alignment=item["sequence"],
                    subject_alignment=item["sequence"],
                )
            ]
            for item in queries
        }


class FakeBulkAnalyzer(FakeAnalyzer):
    def __init__(self) -> None:
        self.blast_client = FakeBulkBlastClient()


class FakeUniProtClient:
    def _fetch_entry(self, accession_or_id: str):
        self.accession_or_id = accession_or_id
        return {
            "primaryAccession": "P00533",
            "uniProtkbId": "EGFR_HUMAN",
            "organism": {"scientificName": "Homo sapiens", "taxonId": 9606},
            "sequence": {"value": "M" * 80},
            "genes": [
                {
                    "geneName": {"value": "EGFR"},
                    "synonyms": [{"value": "ERBB"}],
                }
            ],
            "proteinDescription": {
                "recommendedName": {"fullName": {"value": "Epidermal growth factor receptor"}},
                "alternativeNames": [
                    {
                        "fullName": {"value": "Receptor tyrosine-protein kinase erbB-1"},
                        "shortNames": [{"value": "ERBB1"}],
                    }
                ],
            },
        }

    def _target_from_entry(self, entry):
        from agdesign2.uniprot import UniProtClient

        return UniProtClient.__new__(UniProtClient)._target_from_entry(entry)


class FakeMetadataAnalyzer:
    def __init__(self) -> None:
        self.uniprot_client = FakeUniProtClient()


class FakeComplexPortalAnalyzer:
    def __init__(self) -> None:
        self.uniprot_client = FakeUniProtClient()
        self.lookups: list[str] = []

    @staticmethod
    def _is_canonical_integrin_chain(gene_symbol: str | None) -> bool:
        return str(gene_symbol or "").startswith("ITGA") or str(gene_symbol or "") in {f"ITGB{i}" for i in range(1, 9)}

    def _fetch_complex_portal_context(self, *, target, homolog_context, notes, lookup_status):
        self.lookups.append(target.gene_symbol or "")
        if target.gene_symbol == "ERRO":
            lookup_status["human"] = {"status": "error", "error": "service unavailable"}
            notes.append(AnalysisNote(severity="warning", message="Human Complex Portal lookup failed: service unavailable", source="Complex Portal"))
            return []
        if target.gene_symbol != "ITGA5":
            lookup_status["human"] = {"status": "no_hits"}
            return []
        lookup_status["human"] = {"status": "ok"}
        return [
            ComplexPortalComplex(
                complex_ac="CPX-ITGA5-ITGB1",
                name="Integrin alpha-5/beta-1 complex",
                species="Homo sapiens; 9606",
                predicted_complex=False,
                participants=[
                    ComplexPortalParticipant(identifier="P08648", name="ITGA5", interactor_type="protein"),
                    ComplexPortalParticipant(identifier="P05556", name="ITGB1", interactor_type="protein"),
                ],
            )
        ]

    def _collect_assembly_requirements(self, **kwargs):
        return [
            AssemblyRequirement(
                summary="Integrin alpha chains require an integrin beta partner.",
                obligatory=True,
                confidence="high",
                source="Complex Portal",
                classification="obligatory_partner_requirement",
                partners=["ITGB1"],
                complex_portal_ids=["CPX-ITGA5-ITGB1"],
            )
        ]


class ModuleRefreshTests(unittest.TestCase):
    def test_complex_portal_default_analyzer_enables_the_module(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            batch_dir = Path(tmpdir)
            report_path = batch_dir / "egfr_human_report.json"
            report_path.write_text(
                json.dumps(
                    {
                        "target": {
                            "accession": "P00533",
                            "entry_name": "EGFR_HUMAN",
                            "gene_symbol": "EGFR",
                            "protein_name": "EGFR",
                            "organism": "Homo sapiens",
                            "taxon_id": 9606,
                            "sequence": "M" * 80,
                        }
                    }
                ),
                encoding="utf-8",
            )
            summary_path = batch_dir / "batch_summary.json"
            summary_path.write_text(json.dumps([{"query": "EGFR", "json_report": str(report_path)}]), encoding="utf-8")
            analyzer = FakeComplexPortalAnalyzer()
            with mock.patch("agdesign2.module_refresh.AntigenAnalyzer", return_value=analyzer) as constructor:
                refresh_report_modules(summary_path, modules=["complex-portal"])

        self.assertTrue(constructor.call_args.kwargs["config"].enable_complex_portal)

    def test_refresh_complex_portal_records_every_human_lookup_and_limits_assembly_to_integrins(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            batch_dir = Path(tmpdir)
            reports = []
            for gene, accession in (("ITGA5", "P08648"), ("EGFR", "P00533"), ("ERRO", "P99999")):
                report_path = batch_dir / f"{gene.lower()}_human_report.json"
                report_path.write_text(
                    json.dumps(
                        {
                            "target": {
                                "accession": accession,
                                "entry_name": f"{gene}_HUMAN",
                                "gene_symbol": gene,
                                "protein_name": gene,
                                "organism": "Homo sapiens",
                                "taxon_id": 9606,
                                "sequence": "M" * 80,
                            },
                            "construct_details": [{"name": "preserve-me"}],
                            "assembly_requirements": [{"classification": "legacy"}],
                            "notes": [{"severity": "info", "message": "keep this note", "source": "test"}],
                        }
                    ),
                    encoding="utf-8",
                )
                reports.append(report_path)
            summary_path = batch_dir / "batch_summary.json"
            summary_path.write_text(
                json.dumps(
                    [
                        {"query": path.stem, "resolved_entry_name": path.stem.upper(), "json_report": str(path), "status": "ok"}
                        for path in reports
                    ]
                ),
                encoding="utf-8",
            )
            analyzer = FakeComplexPortalAnalyzer()
            summary_writes: list[Path] = []
            original_write_text = Path.write_text

            def count_summary_writes(path, *args, **kwargs):
                if path == summary_path:
                    summary_writes.append(path)
                return original_write_text(path, *args, **kwargs)

            with mock.patch.object(Path, "write_text", new=count_summary_writes):
                _, _, updated = refresh_report_modules(summary_path, modules=["complex-portal"], analyzer=analyzer)
            itga5, egfr, erro = [json.loads(path.read_text(encoding="utf-8")) for path in reports]

        self.assertEqual(updated, 3)
        self.assertEqual(summary_writes, [summary_path])
        self.assertEqual(analyzer.lookups, ["ITGA5", "EGFR", "ERRO"])
        self.assertEqual(itga5["complex_portal_lookup"]["human"]["status"], "ok")
        self.assertEqual(itga5["complex_portal_complexes"][0]["complex_ac"], "CPX-ITGA5-ITGB1")
        self.assertTrue(itga5["assembly_requirements"][0]["obligatory"])
        self.assertEqual(egfr["complex_portal_lookup"]["human"]["status"], "no_hits")
        self.assertEqual(egfr["complex_portal_complexes"], [])
        self.assertEqual(egfr["assembly_requirements"], [])
        self.assertEqual(erro["complex_portal_lookup"]["human"]["status"], "error")
        self.assertIn("Human Complex Portal lookup failed", erro["notes"][-1]["message"])
        self.assertEqual(egfr["construct_details"], [{"name": "preserve-me"}])
        self.assertEqual(egfr["notes"], [{"severity": "info", "message": "keep this note", "source": "test"}])
    def test_refresh_cross_reactivity_updates_report_json_in_place(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            batch_dir = tmpdir_path / "batch"
            batch_dir.mkdir()
            report_json = batch_dir / "egfr_human_report.json"
            report_json.write_text(
                json.dumps(
                    {
                        "target": {
                            "accession": "P00533",
                            "entry_name": "EGFR_HUMAN",
                            "gene_symbol": "EGFR",
                            "protein_name": "Epidermal growth factor receptor",
                            "organism": "Homo sapiens",
                            "taxon_id": 9606,
                            "sequence": "M" * 80,
                        },
                        "ectodomain": {"start": 25, "end": 60, "label": "ecto", "source": "topology"},
                        "cross_reactivity_hits": [],
                    }
                ),
                encoding="utf-8",
            )
            summary_path = batch_dir / "batch_summary.json"
            summary_path.write_text(
                json.dumps(
                    [
                        {
                            "query": "EGFR_HUMAN",
                            "source_column": "uniprot_name",
                            "batch_index": 1,
                            "batch_total": 1,
                            "status": "ok",
                            "resolved_entry_name": "EGFR_HUMAN",
                            "json_report": str(report_json),
                            "markdown_report": str(batch_dir / "egfr_human_report.md"),
                            "duration_seconds": 1.0,
                            "error": None,
                        }
                    ]
                ),
                encoding="utf-8",
            )

            _, results, updated = refresh_report_modules(
                summary_path,
                modules=["cross-reactivity"],
                analyzer=FakeAnalyzer(),
            )
            payload = json.loads(report_json.read_text(encoding="utf-8"))

        self.assertEqual(updated, 1)
        self.assertEqual(results[0]["status"], "ok")
        self.assertEqual(len(payload["cross_reactivity_hits"]), 1)
        self.assertEqual(payload["cross_reactivity_hits"][0]["subject_id"], "sp|Q01279|EGFR_MOUSE")

    def test_refresh_target_metadata_updates_uniprot_names_in_place(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            batch_dir = tmpdir_path / "batch"
            batch_dir.mkdir()
            report_json = batch_dir / "egfr_human_report.json"
            report_json.write_text(
                json.dumps(
                    {
                        "target": {
                            "accession": "P00533",
                            "entry_name": "EGFR_HUMAN",
                            "gene_symbol": "EGFR",
                            "protein_name": "Old name",
                            "organism": "Homo sapiens",
                            "taxon_id": 9606,
                            "sequence": "M" * 80,
                        }
                    }
                ),
                encoding="utf-8",
            )
            summary_path = batch_dir / "batch_summary.json"
            summary_path.write_text(
                json.dumps(
                    [
                        {
                            "query": "EGFR_HUMAN",
                            "source_column": "uniprot_name",
                            "batch_index": 1,
                            "batch_total": 1,
                            "status": "ok",
                            "resolved_entry_name": "EGFR_HUMAN",
                            "json_report": str(report_json),
                            "markdown_report": str(batch_dir / "egfr_human_report.md"),
                            "duration_seconds": 1.0,
                            "error": None,
                        }
                    ]
                ),
                encoding="utf-8",
            )

            _, results, updated = refresh_report_modules(
                summary_path,
                modules=["target-metadata"],
                analyzer=FakeMetadataAnalyzer(),
            )
            payload = json.loads(report_json.read_text(encoding="utf-8"))

        self.assertEqual(updated, 1)
        self.assertEqual(results[0]["status"], "ok")
        self.assertEqual(payload["target"]["protein_name"], "Epidermal growth factor receptor")
        self.assertIn("Receptor tyrosine-protein kinase erbB-1", payload["target"]["alternative_names"])
        self.assertIn("ERBB1", payload["target"]["alternative_names"])
        self.assertIn("ERBB", payload["target"]["gene_synonyms"])

    def test_refresh_cross_reactivity_updates_full_length_hits_for_compact_multipass(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            batch_dir = tmpdir_path / "batch"
            batch_dir.mkdir()
            report_json = batch_dir / "oprm_human_report.json"
            report_json.write_text(
                json.dumps(
                    {
                        "target": {
                            "accession": "P35372",
                            "entry_name": "OPRM_HUMAN",
                            "gene_symbol": "OPRM1",
                            "protein_name": "Mu-type opioid receptor",
                            "organism": "Homo sapiens",
                            "taxon_id": 9606,
                            "sequence": "M" * 400,
                        },
                        "ectodomain": None,
                        "topology": {"topology_class": "multipass_compact"},
                        "cross_reactivity_hits": [],
                        "full_length_cross_reactivity_hits": [],
                    }
                ),
                encoding="utf-8",
            )
            summary_path = batch_dir / "batch_summary.json"
            summary_path.write_text(
                json.dumps(
                    [
                        {
                            "query": "OPRM_HUMAN",
                            "source_column": "uniprot_name",
                            "batch_index": 1,
                            "batch_total": 1,
                            "status": "ok",
                            "resolved_entry_name": "OPRM_HUMAN",
                            "json_report": str(report_json),
                            "markdown_report": str(batch_dir / "oprm_human_report.md"),
                            "duration_seconds": 1.0,
                            "error": None,
                        }
                    ]
                ),
                encoding="utf-8",
            )
            analyzer = FakeAnalyzer()

            _, _, updated = refresh_report_modules(
                summary_path,
                modules=["cross-reactivity"],
                only_zero_cross_reactivity=True,
                analyzer=analyzer,
            )
            payload = json.loads(report_json.read_text(encoding="utf-8"))

        self.assertEqual(updated, 1)
        self.assertEqual(analyzer.blast_client.queries, ["M" * 400])
        self.assertEqual(payload["cross_reactivity_hits"], [])
        self.assertEqual(len(payload["full_length_cross_reactivity_hits"]), 1)

    def test_only_zero_cross_reactivity_skips_reports_with_relevant_hits(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            batch_dir = tmpdir_path / "batch"
            batch_dir.mkdir()
            report_json = batch_dir / "egfr_human_report.json"
            report_json.write_text(
                json.dumps(
                    {
                        "target": {
                            "accession": "P00533",
                            "entry_name": "EGFR_HUMAN",
                            "gene_symbol": "EGFR",
                            "protein_name": "Epidermal growth factor receptor",
                            "organism": "Homo sapiens",
                            "taxon_id": 9606,
                            "sequence": "M" * 80,
                        },
                        "ectodomain": {"start": 25, "end": 60, "label": "ecto", "source": "topology"},
                        "cross_reactivity_hits": [{"subject_id": "already_present"}],
                    }
                ),
                encoding="utf-8",
            )
            summary_path = batch_dir / "batch_summary.json"
            summary_path.write_text(
                json.dumps(
                    [
                        {
                            "query": "EGFR_HUMAN",
                            "status": "ok",
                            "resolved_entry_name": "EGFR_HUMAN",
                            "json_report": str(report_json),
                        }
                    ]
                ),
                encoding="utf-8",
            )
            analyzer = FakeAnalyzer()

            _, _, updated = refresh_report_modules(
                summary_path,
                modules=["cross-reactivity"],
                only_zero_cross_reactivity=True,
                analyzer=analyzer,
            )

        self.assertEqual(updated, 0)
        self.assertEqual(analyzer.blast_client.queries, [])

    def test_refresh_cross_reactivity_uses_bulk_blast_when_available(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            batch_dir = tmpdir_path / "batch"
            batch_dir.mkdir()
            reports = []
            summary_rows = []
            for index, entry in enumerate(("AAA_HUMAN", "BBB_HUMAN"), start=1):
                report_json = batch_dir / f"{entry.lower()}_report.json"
                report_json.write_text(
                    json.dumps(
                        {
                            "target": {
                                "accession": f"P{index:05d}",
                                "entry_name": entry,
                                "gene_symbol": entry.split("_", 1)[0],
                                "protein_name": entry,
                                "organism": "Homo sapiens",
                                "taxon_id": 9606,
                                "sequence": "M" * 80,
                            },
                            "ectodomain": {"start": 10, "end": 40, "label": "ecto", "source": "topology"},
                            "residue_annotations": [
                                {"position": 10, "topology_location": "extracellular", "surface_exposed": True},
                                {"position": 12, "topology_location": "extracellular", "surface_exposed": False},
                                {"position": 14, "topology_location": "cytoplasmic", "surface_exposed": True},
                            ],
                            "cross_reactivity_hits": [],
                        }
                    ),
                    encoding="utf-8",
                )
                reports.append(report_json)
                summary_rows.append(
                    {
                        "query": entry,
                        "status": "ok",
                        "resolved_entry_name": entry,
                        "json_report": str(report_json),
                    }
                )
            summary_path = batch_dir / "batch_summary.json"
            summary_path.write_text(json.dumps(summary_rows), encoding="utf-8")
            analyzer = FakeBulkAnalyzer()

            _, _, updated = refresh_report_modules(
                summary_path,
                modules=["cross-reactivity"],
                analyzer=analyzer,
            )
            first_payload = json.loads(reports[0].read_text(encoding="utf-8"))
            second_payload = json.loads(reports[1].read_text(encoding="utf-8"))

        self.assertEqual(updated, 2)
        self.assertEqual(analyzer.blast_client.queries, [])
        self.assertEqual(len(analyzer.blast_client.batch_queries), 2)
        self.assertEqual(first_payload["cross_reactivity_hits"][0]["subject_id"], "sp|Q01279|R0_ECTO_MOUSE")
        self.assertEqual(first_payload["cross_reactivity_hits"][0]["query_alignment_annotation"][0], "S")
        self.assertEqual(first_payload["cross_reactivity_hits"][0]["extracellular_surface_aligned_positions"], 1)
        self.assertEqual(second_payload["cross_reactivity_hits"][0]["subject_id"], "sp|Q01279|R1_ECTO_MOUSE")


if __name__ == "__main__":
    unittest.main()
