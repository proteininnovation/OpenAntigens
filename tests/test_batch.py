from __future__ import annotations

import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from agdesign2.batch import BatchRow, load_batch_rows, rebuild_batch_summary_from_tsv, render_batch_summary_markdown, retry_batch_from_tsv, run_batch_from_tsv, write_batch_summary_markdown
from agdesign2.cache import FileCache
from agdesign2.config import AnalysisConfig
from agdesign2.portal import _PORTAL_SEQUENCE_CACHE, _build_sequence_structure_mapping, _construct_group_label, _copy_report_assets, _ensure_vendor_assets, _features_in_region, _format_alignment_block_html, _load_homolog_mappings, _load_viewer_payload, _portal_align_sequences, _public_gpcr_variant, _pubtator_payload_count, _render_construct_card, _render_cross_reactivity_alignment, _render_obligatory_partner_warning, _render_pubtator_link, _render_structure_widget, build_portal, render_detail_page
from agdesign2.models import AnalysisReport, AnalysisNote, AssemblyRequirement, BlastHit, ConstructDetail, ConstructSuggestion, CysteineFinding, DomainAnnotation, ExperimentalConstruct, FamilyContext, Feature, FurinSite, HomologyRecord, Region, TargetRecord


class PortalHomologMappingTests(unittest.TestCase):
    def test_loads_human_homolog_sequence_from_batch_cache_for_mouse_portal(self) -> None:
        _PORTAL_SEQUENCE_CACHE.clear()
        with tempfile.TemporaryDirectory() as tmpdir:
            batch_dir = Path(tmpdir) / "mouse_data"
            cache = FileCache(batch_dir / ".agdesign2" / "cache")
            url = "https://rest.uniprot.org/uniprotkb/P99999.json"
            cache.path_for("uniprot_entry", url, ".json").write_text(
                json.dumps(
                    {
                        "primaryAccession": "P99999",
                        "uniProtkbId": "TNFA_HUMAN",
                        "genes": [{"geneName": {"value": "TNF"}}],
                        "sequence": {"value": "M" * 56 + "HUMANREGION" + "K" * 20},
                    }
                ),
                encoding="utf-8",
            )
            report = {
                "target": {
                    "accession": "P06804",
                    "entry_name": "TNFA_MOUSE",
                    "gene_symbol": "Tnf",
                    "sequence": "M" * 56 + "MOUSEREGION" + "K" * 20,
                },
                "ectodomain": {"start": 57, "end": 67},
                "ectodomain_homology": [
                    {
                        "species": "human",
                        "accession": "P99999",
                        "entry_name": "TNFA_HUMAN",
                        "available": True,
                        "ectodomain_start": 57,
                        "ectodomain_end": 67,
                    }
                ],
            }

            mappings = _load_homolog_mappings(report=report, batch_dir=batch_dir)

        self.assertEqual(len(mappings), 1)
        self.assertEqual(mappings[0]["species"], "human")
        self.assertEqual(mappings[0]["entryName"], "TNFA_HUMAN")
        self.assertEqual(mappings[0]["ectodomainSequence"], "HUMANREGION")


class FakeAnalyzer:
    def analyze_target(self, query: str, output_dir: str | Path | None = None) -> AnalysisReport:
        output_root = Path(output_dir or ".")
        output_root.mkdir(parents=True, exist_ok=True)
        target = TargetRecord(
            accession="P00001",
            entry_name=f"{query}_HUMAN" if "_" not in query else query,
            gene_symbol=query.split("_", 1)[0],
            protein_name=query,
            organism="Homo sapiens",
            taxon_id=9606,
            sequence="M" * 100,
        )
        report = AnalysisReport(
            target=target,
            resolution={"query": query},
            features=[],
            ectodomain=Region(start=1, end=10, label="ecto", source="test"),
            structural_regions=[],
            experimental_constructs=[],
            construct_recommendations=[],
            furin_sites=[],
            cysteine_analysis=[],
            assembly_requirements=[],
            interpro_annotations=[],
            species_name_matches=[],
            ectodomain_homology=[],
            family_context=None,
            cross_reactivity_hits=[],
            complex_portal_complexes=[],
            ligand_interactions=[],
            construct_details=[],
            notes=[],
        )
        stem = f"{target.entry_name.lower()}_report"
        report.write_json(output_root / f"{stem}.json")
        report.write_markdown(output_root / f"{stem}.md")
        return report


class CountingAnalyzer(FakeAnalyzer):
    def __init__(self) -> None:
        self.calls: list[str] = []

    def analyze_target(self, query: str, output_dir: str | Path | None = None) -> AnalysisReport:
        self.calls.append(query)
        return super().analyze_target(query, output_dir=output_dir)


class RefreshingPortalAnalyzer(FakeAnalyzer):
    def analyze_target(self, query: str, output_dir: str | Path | None = None) -> AnalysisReport:
        output_root = Path(output_dir or ".")
        output_root.mkdir(parents=True, exist_ok=True)
        report_json = output_root / "unc5d_human_report.json"
        report_md = output_root / "unc5d_human_report.md"
        report_json.write_text(
            json.dumps(
                {
                    "target": {
                        "accession": "Q6UXZ4",
                        "entry_name": "UNC5D_HUMAN",
                        "gene_symbol": "UNC5D",
                        "protein_name": "Netrin receptor UNC5D",
                        "organism": "Homo sapiens",
                        "taxon_id": 9606,
                        "sequence": "M" * 120,
                    },
                    "ectodomain": {"start": 25, "end": 400, "label": "ecto", "source": "topology"},
                    "construct_details": [],
                    "construct_recommendations": [],
                    "cross_reactivity_hits": [],
                    "furin_sites": [],
                    "cysteine_analysis": [],
                    "ectodomain_homology": [],
                    "family_context": {
                        "source": "Precomputed InterPro family alignments",
                        "family_names": ["Netrin receptor UNC5A-D"],
                        "members": [
                            {"gene_symbol": "UNC5A", "entry_name": "UNC5A_HUMAN", "ectodomain_start": 30, "ectodomain_end": 380},
                            {"gene_symbol": "UNC5D", "entry_name": "UNC5D_HUMAN", "ectodomain_start": 25, "ectodomain_end": 400},
                        ],
                        "identity_matrix_labels": ["UNC5A_HUMAN", "UNC5D_HUMAN"],
                        "identity_matrix": [[100.0, 51.2], [49.8, 100.0]],
                    },
                    "assembly_requirements": [],
                    "experimental_constructs": [],
                    "interpro_annotations": [],
                    "notes": [],
                    "canonical_family": {
                        "accession": "IPR037936",
                        "name": "Netrin receptor UNC5A-D",
                        "source_database": "INTERPRO",
                        "fragment_count": 1,
                    },
                }
            ),
            encoding="utf-8",
        )
        report_md.write_text("# UNC5D report\n", encoding="utf-8")
        return AnalysisReport(
            target=TargetRecord(
                accession="Q6UXZ4",
                entry_name="UNC5D_HUMAN",
                gene_symbol="UNC5D",
                protein_name="Netrin receptor UNC5D",
                organism="Homo sapiens",
                taxon_id=9606,
                sequence="M" * 120,
            ),
            resolution={"query": query},
            features=[],
            ectodomain=Region(start=25, end=400, label="ecto", source="topology"),
            structural_regions=[],
            experimental_constructs=[],
            construct_recommendations=[],
            furin_sites=[],
            cysteine_analysis=[],
            assembly_requirements=[],
            interpro_annotations=[],
            species_name_matches=[],
            ectodomain_homology=[],
            family_context=None,
            cross_reactivity_hits=[],
            complex_portal_complexes=[],
            ligand_interactions=[],
            construct_details=[],
            notes=[],
        )


class PortalAlignmentTests(unittest.TestCase):
    def test_pubtator_count_helpers_accept_api_and_cached_string_counts(self) -> None:
        self.assertEqual(_pubtator_payload_count({"count": "14923", "results": []}), 14923)
        self.assertIsNone(_pubtator_payload_count({"results": []}))
        self.assertIn(">14,923<", _render_pubtator_link("14923", "https://example.test/pubtator"))

    def test_construct_group_label_prioritizes_specific_evidence(self) -> None:
        self.assertEqual(
            _construct_group_label(
                {
                    "name": "pdb_5w7i_a_c",
                    "classification": "membrane_expression",
                }
            ),
            "PDB-backed",
        )
        self.assertEqual(
            _construct_group_label(
                {
                    "name": "strict_domain_1",
                    "classification": "membrane_expression",
                    "evidence": ["AlphaFold strict region"],
                }
            ),
            "AlphaFold strict",
        )

    def test_detail_page_omits_empty_membrane_construct_sections(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            detail_html = render_detail_page(
                {
                    "entry_name": "TEST_HUMAN",
                    "topology_bucket": "Secreted",
                    "status": "ok",
                    "has_alphafold_structure": False,
                },
                {
                    "target": {
                        "accession": "P00001",
                        "entry_name": "TEST_HUMAN",
                        "gene_symbol": "TEST",
                        "protein_name": "Test protein",
                        "organism": "Homo sapiens",
                        "sequence": "M" * 40,
                    },
                    "ectodomain": {
                        "start": 1,
                        "end": 40,
                        "label": "secreted region",
                        "source": "topology",
                    },
                    "construct_details": [
                        {
                            "name": "full_ectodomain",
                            "start": 1,
                            "end": 40,
                            "length": 40,
                            "sequence": "M" * 40,
                            "classification": "complete_domain",
                            "score": 1.0,
                            "rationale": "Topology-derived design region.",
                            "warnings": [],
                            "evidence": ["Topology-derived design region"],
                        }
                    ],
                },
                batch_dir=Path(tmpdir),
            )

        self.assertNotIn("Native Membrane-Expression Constructs", detail_html)
        self.assertNotIn("No membrane-expression constructs available.", detail_html)

    def test_portal_align_sequences_returns_equal_length_alignment(self) -> None:
        alignment = _portal_align_sequences("ABCDEFGH", "ABCXEFGH")
        self.assertEqual(len(alignment["aligned_query"]), len(alignment["aligned_subject"]))
        self.assertIn("ABCDEFGH", alignment["aligned_query"].replace("-", ""))
        self.assertIn("ABCXEFGH", alignment["aligned_subject"].replace("-", ""))

    def test_alignment_block_marks_surface_accessible_query_positions(self) -> None:
        html = _format_alignment_block_html(
            "ACD-EFG",
            "ACDTEYG",
            "S.S.S..",
            width=10,
        )
        self.assertIn("alignment-surface-ecto", html)
        self.assertIn('class="alignment-residue alignment-surface-ecto">A</span>', html)
        self.assertNotIn("surface ", html)

    def test_cross_reactivity_alignment_uses_mouse_query_label(self) -> None:
        html = _render_cross_reactivity_alignment(
            {
                "query_alignment": "ACD",
                "subject_alignment": "ACD",
                "query_alignment_annotation": "S..",
                "query_start": 1,
                "query_end": 3,
                "subject_start": 2,
                "subject_end": 4,
            },
            rank=1,
            query_label="mouse",
        )
        self.assertIn("mouse  <span", html)
        self.assertIn("mouse 1-3 | hit 2-4", html)
        self.assertIn("mouse extracellular accessible residue", html)
        self.assertNotIn("human extracellular accessible residue", html)

    def test_gpcr_public_compatibility_relabels_only_explicit_topology_only_loops(self) -> None:
        topology_only = _public_gpcr_variant(
            {
                "name": "TEST_HUMAN_BRIL_ICL3_fusion",
                "strategy": "ICL3 cassette replacement",
                "evidence": ["topology-derived internal cytoplasmic loop: 100-130."],
            }
        )
        self.assertEqual(topology_only["name"], "TEST_HUMAN_BRIL_ICL_fusion")
        self.assertEqual(
            topology_only["strategy"],
            "Intracellular-loop cassette replacement",
        )

        genuine_or_unknown = _public_gpcr_variant(
            {
                "name": "OPRM_HUMAN_BRIL_ICL3_fusion",
                "strategy": "ICL3 cassette replacement",
                "evidence": [],
            }
        )
        self.assertEqual(
            genuine_or_unknown["name"],
            "OPRM_HUMAN_BRIL_ICL3_fusion",
        )
        self.assertEqual(
            genuine_or_unknown["strategy"],
            "ICL3 cassette replacement",
        )

    def test_no_structure_construct_exports_map_target_edits_to_homolog_rows(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            portal_dir = root / "portal"
            (portal_dir / "reports").mkdir(parents=True)
            report = {
                "target": {
                    "accession": "P00001",
                    "entry_name": "TEST_MOUSE",
                    "gene_symbol": "Test",
                    "organism": "Mus musculus",
                    "taxon_id": 10090,
                    "sequence": "ACD",
                },
                "ectodomain": {
                    "start": 1,
                    "end": 3,
                    "label": "design region",
                    "source": "test",
                },
                "construct_details": [
                    {
                        "name": "test_construct",
                        "start": 1,
                        "end": 3,
                        "length": 3,
                        "sequence": "ACD",
                        "classification": "full_ectodomain",
                        "cysteine_analysis": [
                            {"position": 2, "paired_with": None}
                        ],
                        "homologs": [
                            {
                                "species": "macaca_fascicularis",
                                "available": True,
                                "start": 10,
                                "end": 13,
                                "sequence": "ACED",
                                "identity_to_human": 75.0,
                            }
                        ],
                    }
                ],
                "construct_recommendations": [],
                "cross_reactivity_hits": [],
                "full_length_cross_reactivity_hits": [],
                "experimental_constructs": [],
                "family_context": {},
                "assembly_requirements": [],
                "notes": [],
                "metadata": {
                    "construct_identity_label": "Identity to mouse source"
                },
            }
            entry = {
                "query": "TEST_MOUSE",
                "entry_name": "TEST_MOUSE",
                "detail_page": "test_mouse.html",
                "status": "ok",
            }

            detail = render_detail_page(
                entry,
                report,
                batch_dir=root,
                portal_dir=portal_dir,
                portal_title="OpenAntigens Mouse",
                include_disease_context=False,
            )

        self.assertIn("mapTargetMutationEditsToRow", detail)
        self.assertIn("alignedTarget", detail)
        self.assertIn("alignedHomolog", detail)
        self.assertIn("targetMutationEdits", detail)
        self.assertNotIn("View in 3D", detail)
        self.assertNotIn("Curated Complex Portal records", detail)

    def test_sequence_structure_mapping_handles_n_terminal_offset(self) -> None:
        pdb_text = """\
ATOM      1  CA  GLY A 101       0.000   0.000   0.000  1.00 90.00           C
ATOM      2  CA  MET A 102       0.000   0.000   0.000  1.00 90.00           C
ATOM      3  CA  SER A 103       0.000   0.000   0.000  1.00 90.00           C
ATOM      4  CA  GLN A 104       0.000   0.000   0.000  1.00 90.00           C
ATOM      5  CA  ASP A 105       0.000   0.000   0.000  1.00 90.00           C
"""
        with tempfile.TemporaryDirectory() as tmpdir:
            pdb_path = Path(tmpdir) / "offset.pdb"
            pdb_path.write_text(pdb_text, encoding="utf-8")
            mapping = _build_sequence_structure_mapping("MSQD", pdb_path)
        self.assertEqual(mapping["sequenceToStructureResidues"], [102, 103, 104, 105])
        self.assertEqual(mapping["structureToSequenceResidues"]["105"], 4)

    def test_multipass_structure_widget_uses_full_length_plot_scope(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            structures = root / "portal" / "structures"
            structures.mkdir(parents=True)
            pdb_path = structures / "test_human.pdb"
            pdb_path.write_text(
                "\n".join(
                    [
                        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 95.00           C",
                        "ATOM      2  CA  ALA A   2       1.000   0.000   0.000  1.00 95.00           C",
                        "END",
                    ]
                ),
                encoding="utf-8",
            )
            report = {
                "target": {"entry_name": "TEST_HUMAN", "accession": "P00001", "sequence": "AA"},
                "ectodomain": {"start": 1, "end": 1, "label": "small ecto", "source": "test"},
                "topology": {"topology_class": "multipass_mixed"},
            }
            html = _render_structure_widget(
                {"portal_structure_path": str(pdb_path)},
                report,
                batch_dir=root,
            )
        self.assertIn("plotScopeFullLength = true", html)
        self.assertIn("Viewing the full-length protein pLDDT trace", html)

    def test_mouse_structure_widget_uses_mouse_as_selected_reference_species(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            structures = root / "portal" / "structures"
            structures.mkdir(parents=True)
            pdb_path = structures / "test_mouse.pdb"
            pdb_path.write_text(
                "\n".join(
                    [
                        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 95.00           C",
                        "ATOM      2  CA  ALA A   2       1.000   0.000   0.000  1.00 95.00           C",
                        "END",
                    ]
                ),
                encoding="utf-8",
            )
            report = {
                "target": {
                    "entry_name": "TEST_MOUSE",
                    "accession": "P00001",
                    "organism": "Mus musculus",
                    "taxon_id": 10090,
                    "sequence": "AA",
                },
                "metadata": {"construct_identity_label": "Identity to human source"},
                "ectodomain": {"start": 1, "end": 2, "label": "small ecto", "source": "test"},
                "topology": {"topology_class": "single_pass"},
            }
            html = _render_structure_widget(
                {"portal_structure_path": str(pdb_path)},
                report,
                batch_dir=root,
            )
        self.assertIn('species: "mouse"', html)
        self.assertIn("Identity to human source", html)
        self.assertIn('identityHeader: "identity_to_human_source"', html)
        self.assertIn("nameWithMutationEdits(targetNameSymbol, targetSpecies", html)
        self.assertIn("if (normalized === targetSpeciesKey()) return mutationEdits", html)

    def test_compact_multipass_viewer_payload_loads_full_length_pae_without_ectodomain(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            af_dir = root / ".agdesign2" / "data" / "alphafold"
            af_dir.mkdir(parents=True)
            (af_dir / "P00001.pdb").write_text(
                "\n".join(
                    [
                        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 95.00           C",
                        "ATOM      2  CA  ALA A   2       1.000   0.000   0.000  1.00 80.00           C",
                        "ATOM      3  CA  ALA A   3       2.000   0.000   0.000  1.00 70.00           C",
                        "END",
                    ]
                ),
                encoding="utf-8",
            )
            (af_dir / "P00001.pae.json").write_text(
                json.dumps({"predicted_aligned_error": [[0, 1, 2], [1, 0, 3], [2, 3, 0]]}),
                encoding="utf-8",
            )
            payload = _load_viewer_payload(
                report={
                    "target": {"entry_name": "TEST_HUMAN", "accession": "P00001", "sequence": "AAA"},
                    "ectodomain": None,
                    "topology": {"topology_class": "multipass_compact"},
                },
                batch_dir=root,
            )
        self.assertEqual(payload["plddt"], [95.0, 80.0, 70.0])
        self.assertEqual(payload["paeMatrix"], [[0.0, 1.0, 2.0], [1.0, 0.0, 3.0], [2.0, 3.0, 0.0]])

    def test_viewer_payload_loads_full_length_pae_when_single_pass_ectodomain_unresolved(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            af_dir = root / ".agdesign2" / "data" / "alphafold"
            af_dir.mkdir(parents=True)
            (af_dir / "P00001.pdb").write_text(
                "\n".join(
                    [
                        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 95.00           C",
                        "ATOM      2  CA  ALA A   2       1.000   0.000   0.000  1.00 80.00           C",
                        "END",
                    ]
                ),
                encoding="utf-8",
            )
            (af_dir / "P00001.pae.json").write_text(
                json.dumps({"predicted_aligned_error": [[0, 4], [4, 0]]}),
                encoding="utf-8",
            )
            payload = _load_viewer_payload(
                report={
                    "target": {"entry_name": "TEST_HUMAN", "accession": "P00001", "sequence": "AA"},
                    "ectodomain": None,
                    "topology": {"topology_class": "single_pass"},
                },
                batch_dir=root,
            )
        self.assertEqual(payload["paeMatrix"], [[0.0, 4.0], [4.0, 0.0]])

    def test_viewer_payload_rejects_pae_with_wrong_sequence_dimensions(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            af_dir = root / ".agdesign2" / "data" / "alphafold"
            af_dir.mkdir(parents=True)
            (af_dir / "P00001.pdb").write_text(
                "\n".join(
                    [
                        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 95.00           C",
                        "ATOM      2  CA  ALA A   2       1.000   0.000   0.000  1.00 80.00           C",
                        "END",
                    ]
                ),
                encoding="utf-8",
            )
            (af_dir / "P00001.pae.json").write_text(
                json.dumps({"predicted_aligned_error": [[0.0]]}),
                encoding="utf-8",
            )

            payload = _load_viewer_payload(
                report={
                    "target": {
                        "entry_name": "TEST_HUMAN",
                        "accession": "P00001",
                        "sequence": "AA",
                    },
                    "ectodomain": None,
                    "topology": {"topology_class": "single_pass"},
                },
                batch_dir=root,
            )

        self.assertEqual(payload["plddt"], [95.0, 80.0])
        self.assertNotIn("paeMatrix", payload)

    def test_homolog_mappings_use_precomputed_ranges_when_target_ectodomain_unresolved(self) -> None:
        report = {
            "target": {"entry_name": "TEST_HUMAN", "accession": "P00001", "sequence": "AAAA"},
            "ectodomain": None,
            "topology": {"topology_class": "single_pass"},
            "ectodomain_homology": [
                {
                    "species": "mouse",
                    "available": True,
                    "accession": "NP_000001",
                    "entry_name": "NP_000001",
                    "ectodomain_start": 1,
                    "ectodomain_end": 4,
                }
            ],
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            (root / "ortholog_reference_table.tsv").write_text(
                "\t".join(
                    [
                        "human_gene_symbol",
                        "human_refseq_accession",
                        "human_refseq_sequence",
                        "mouse_gene_symbol",
                        "mouse_refseq_accession",
                        "mouse_refseq_sequence",
                        "macaca_fascicularis_gene_symbol",
                        "macaca_fascicularis_refseq_accession",
                        "macaca_fascicularis_refseq_sequence",
                    ]
                )
                + "\nTEST\tNP_HUMAN\tAAAA\tTest\tNP_000001\tAAAT\t\t\t\n",
                encoding="utf-8",
            )
            mappings = _load_homolog_mappings(report=report, batch_dir=root)
        self.assertEqual(len(mappings), 1)
        self.assertEqual(mappings[0]["species"], "mouse")
        self.assertEqual(mappings[0]["ectodomainSequence"], "AAAT")


class PortalDisplayClassificationTests(unittest.TestCase):
    def test_detail_page_uses_curated_gpi_track_for_display_labels(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            batch_dir = Path(tmpdir) / "batch"
            portal_dir = Path(tmpdir) / "portal"
            batch_dir.mkdir()
            report = {
                "target": {
                    "accession": "P20827",
                    "entry_name": "EFNA1_HUMAN",
                    "gene_symbol": "EFNA1",
                    "protein_name": "Ephrin-A1",
                    "organism": "Homo sapiens",
                    "sequence": "M" * 205,
                },
                "topology": {"topology_class": "secreted"},
                "ectodomain": {"start": 19, "end": 182, "label": "Secreted mature chain", "source": "topology"},
                "construct_details": [],
                "construct_recommendations": [],
                "cross_reactivity_hits": [],
                "experimental_constructs": [],
                "family_context": {},
                "assembly_requirements": [],
                "notes": [],
            }
            entry = {
                "query": "EFNA1_HUMAN",
                "entry_name": "EFNA1_HUMAN",
                "topology_bucket": "GPI",
                "primary_topology_class": "gpi_anchored_surface",
                "status": "ok",
            }

            html = render_detail_page(entry, report, batch_dir=batch_dir, portal_dir=portal_dir)

        self.assertIn("<dt>Track</dt><dd>GPI</dd>", html)
        self.assertIn("<h2>GPI-Anchored Extracellular Region</h2>", html)
        self.assertIn("Soluble GPI-Anchored Extracellular-region Constructs", html)
        self.assertNotIn("<h2>Secreted Region</h2>", html)
        self.assertNotIn("Soluble Secreted-region Constructs", html)

    def test_detail_page_keeps_true_secreted_track_wording(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            batch_dir = Path(tmpdir) / "batch"
            portal_dir = Path(tmpdir) / "portal"
            batch_dir.mkdir()
            report = {
                "target": {
                    "accession": "O94907",
                    "entry_name": "DKK1_HUMAN",
                    "gene_symbol": "DKK1",
                    "protein_name": "Dickkopf-related protein 1",
                    "organism": "Homo sapiens",
                    "sequence": "M" * 266,
                },
                "topology": {"topology_class": "secreted"},
                "ectodomain": {"start": 32, "end": 266, "label": "Secreted mature chain", "source": "topology"},
                "construct_details": [],
                "construct_recommendations": [],
                "cross_reactivity_hits": [],
                "experimental_constructs": [],
                "family_context": {},
                "assembly_requirements": [],
                "notes": [],
            }
            entry = {
                "query": "DKK1_HUMAN",
                "entry_name": "DKK1_HUMAN",
                "topology_bucket": "Secreted",
                "status": "ok",
            }

            html = render_detail_page(entry, report, batch_dir=batch_dir, portal_dir=portal_dir)

        self.assertIn("<dt>Track</dt><dd>Secreted</dd>", html)
        self.assertIn("<h2>Secreted Region</h2>", html)
        self.assertIn("Soluble Secreted-region Constructs", html)


class PortalFeatureOverlapTests(unittest.TestCase):
    def test_obligatory_partner_warning_only_renders_for_obligatory_requirements(self) -> None:
        no_warning = _render_obligatory_partner_warning(
            [
                {
                    "classification": "interaction_context",
                    "obligatory": False,
                    "partners": ["CIB1"],
                    "summary": "Known interaction context.",
                }
            ]
        )
        warning = _render_obligatory_partner_warning(
            [
                {
                    "classification": "obligatory_partner_requirement",
                    "obligatory": True,
                    "partners": ["ITGB1", "ITGB3"],
                    "summary": "Integrin alpha chains require an integrin beta partner.",
                }
            ]
        )

        self.assertEqual(no_warning, "")
        self.assertIn("Complex-aware design note", warning)
        self.assertIn("Obligatory protein-complex context", warning)
        self.assertIn("ITGB1, ITGB3", warning)
        self.assertIn("co-expression", warning)

    def test_disulfide_overlap_uses_bond_endpoints_not_full_span(self) -> None:
        features = [
            {"type": "DISULFID", "start": 264, "end": 329, "description": ""},
            {"type": "CARBOHYD", "start": 306, "end": 306, "description": "N-linked"},
        ]

        between_endpoints = _features_in_region(features, 298, 328)
        overlapping_endpoint = _features_in_region(features, 298, 330)

        self.assertEqual([feature["type"] for feature in between_endpoints], ["CARBOHYD"])
        self.assertEqual([feature["type"] for feature in overlapping_endpoint], ["DISULFID", "CARBOHYD"])

    def test_multipass_homolog_mappings_use_full_length_sequences_without_ectodomain_boundaries(self) -> None:
        report = {
            "target": {
                "entry_name": "5HT1A_HUMAN",
                "accession": "P08908",
                "sequence": "MAAAACCCC",
            },
            "ectodomain": None,
            "topology": {"topology_class": "multipass_compact"},
            "ectodomain_homology": [
                {"available": True, "species": "mouse", "accession": "Q64264", "entry_name": "5HT1A_MOUSE"},
                {"available": True, "species": "macaca_fascicularis", "accession": "XP_045249695", "entry_name": "XP_045249695"},
            ],
        }

        def fake_fetch(record: dict[str, object], *, batch_dir: Path) -> dict[str, object]:
            accession = str(record.get("accession") or "")
            return {
                "accession": accession,
                "entry_name": record.get("entry_name"),
                "gene_symbol": "HTR1A",
                "sequence": "MAAAACCCC" if accession == "Q64264" else "MAAAATCCC",
            }

        with mock.patch("agdesign2.portal._fetch_homolog_target", side_effect=fake_fetch):
            mappings = _load_homolog_mappings(report=report, batch_dir=Path("."))

        self.assertEqual([item["species"] for item in mappings], ["mouse", "macaca_fascicularis"])
        self.assertEqual([item["humanRegionStart"] for item in mappings], [1, 1])
        self.assertEqual([item["humanRegionEnd"] for item in mappings], [9, 9])
        self.assertEqual([item["homologRegionEnd"] for item in mappings], [9, 9])
        self.assertEqual([item["ectodomainSequence"] for item in mappings], ["MAAAACCCC", "MAAAATCCC"])

    def test_homolog_mappings_use_ortholog_table_when_refseq_cache_is_missing(self) -> None:
        report = {
            "target": {
                "entry_name": "FZD1_HUMAN",
                "accession": "Q9UP38",
                "sequence": "MAAABBBBCCCC",
            },
            "ectodomain": {"start": 2, "end": 9},
            "topology": {"topology_class": "single_pass"},
            "ectodomain_homology": [
                {
                    "available": True,
                    "species": "mouse",
                    "accession": "NP_067432",
                    "entry_name": "NP_067432",
                    "ectodomain_start": 3,
                    "ectodomain_end": 10,
                },
            ],
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            batch_dir = Path(tmpdir)
            (batch_dir / "ortholog_reference_table.tsv").write_text(
                "\t".join(
                    [
                        "human_gene_symbol",
                        "mouse_gene_symbol",
                        "macaca_fascicularis_gene_symbol",
                        "human_refseq_accession",
                        "mouse_refseq_accession",
                        "macaca_fascicularis_refseq_accession",
                        "human_refseq_sequence",
                        "mouse_refseq_sequence",
                        "macaca_fascicularis_refseq_sequence",
                    ]
                )
                + "\n"
                + "\t".join(["FZD1", "Fzd1", "FZD1", "NP_003496", "NP_067432", "", "MAAABBBBCCCC", "XXAABBBBCCYY", ""])
                + "\n",
                encoding="utf-8",
            )

            mappings = _load_homolog_mappings(report=report, batch_dir=batch_dir)

        self.assertEqual(len(mappings), 1)
        self.assertEqual(mappings[0]["species"], "mouse")
        self.assertEqual(mappings[0]["geneSymbol"], "Fzd1")
        self.assertEqual(mappings[0]["accession"], "NP_067432")
        self.assertEqual(mappings[0]["ectodomainSequence"], "AABBBBCC")


class PortalRefreshTests(unittest.TestCase):
    def test_build_portal_refreshes_stale_reports_before_rendering(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            batch_dir = tmpdir_path / "batch"
            batch_dir.mkdir()
            report_json = batch_dir / "unc5d_human_report.json"
            report_md = batch_dir / "unc5d_human_report.md"
            report_json.write_text(
                json.dumps(
                    {
                        "target": {
                            "accession": "Q6UXZ4",
                            "entry_name": "UNC5D_HUMAN",
                            "gene_symbol": "UNC5D",
                            "protein_name": "Old UNC5D",
                            "organism": "Homo sapiens",
                            "taxon_id": 9606,
                            "sequence": "M" * 120,
                        },
                        "ectodomain": {"start": 25, "end": 400, "label": "ecto", "source": "topology"},
                        "construct_details": [],
                        "construct_recommendations": [],
                        "cross_reactivity_hits": [
                            {
                                "subject_id": "sp|Q01279|EGFR_MOUSE",
                                "description": "sp|Q01279|EGFR_MOUSE Epidermal growth factor receptor OS=Mus musculus OX=10090 GN=Egfr",
                                "species": "Mus musculus",
                                "identity": 88.7,
                                "coverage": 100.0,
                                "alignment_length": 56,
                                "evalue": 1e-30,
                                "bitscore": 250.0,
                                "query_start": 1,
                                "query_end": 56,
                                "subject_start": 5,
                                "subject_end": 60,
                                "query_alignment": "ACDEFGHIKLMNPQRSTVWY",
                                "subject_alignment": "ACD-FGHIKLMNPQ-STVWY",
                                "alignment_source": "BLAST",
                            }
                        ],
                        "furin_sites": [],
                        "cysteine_analysis": [],
                        "ectodomain_homology": [],
                        "family_context": {
                            "source": "HGNC",
                            "family_names": ["ZU5 domain containing "],
                            "members": [],
                            "identity_matrix_labels": ["ANK1", "UNC5D"],
                            "identity_matrix": [[100.0, 50.0], [60.0, 100.0]],
                        },
                        "assembly_requirements": [],
                        "experimental_constructs": [],
                        "interpro_annotations": [],
                        "notes": [],
                    }
                ),
                encoding="utf-8",
            )
            report_md.write_text("# old\n", encoding="utf-8")
            dependency = batch_dir / "family_alignment_index.tsv"
            dependency.write_text("family_accession\tfamily_name\nIPR037936\tNetrin receptor UNC5A-D\n", encoding="utf-8")
            summary_path = batch_dir / "batch_summary.json"
            summary_path.write_text(
                json.dumps(
                    [
                        {
                            "query": "UNC5D_HUMAN",
                            "source_column": "uniprot_name",
                            "batch_index": 1,
                            "batch_total": 1,
                            "status": "ok",
                            "resolved_entry_name": "UNC5D_HUMAN",
                            "json_report": str(report_json),
                            "markdown_report": str(report_md),
                            "duration_seconds": 1.0,
                            "error": None,
                        }
                    ]
                ),
                encoding="utf-8",
            )
            old_time = 1_700_000_000
            new_time = old_time + 100
            os.utime(report_json, (old_time, old_time))
            os.utime(report_md, (old_time, old_time))
            os.utime(dependency, (new_time, new_time))

            with mock.patch("agdesign2.portal.AntigenAnalyzer", return_value=RefreshingPortalAnalyzer()):
                with mock.patch("agdesign2.portal._portal_refresh_dependencies", return_value=[dependency]):
                    with mock.patch(
                        "agdesign2.portal._get_pubtator_literature_info",
                        return_value={
                            "gene_symbol": "UNC5D",
                            "query_text": "@GENE_UNC5D",
                            "query_url": "https://www.ncbi.nlm.nih.gov/research/pubtator3/docsum?text=%40GENE_UNC5D",
                            "count": 123,
                        },
                    ):
                        index_path = build_portal(summary_path)

            detail_html = (index_path.parent / "reports" / "unc5d_human.html").read_text(encoding="utf-8")
            refreshed_report = json.loads(report_json.read_text(encoding="utf-8"))
            refreshed_summary = json.loads(summary_path.read_text(encoding="utf-8"))

        self.assertIn("Netrin receptor UNC5A-D", detail_html)
        self.assertNotIn("ZU5 domain containing", detail_html)
        self.assertIn("Directional Identity Matrix", detail_html)
        self.assertEqual(refreshed_report["family_context"]["source"], "Precomputed InterPro family alignments")
        self.assertEqual(refreshed_summary[0]["status"], "ok")


class BatchTests(unittest.TestCase):
    def test_public_portal_rejects_af3_derived_report(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            batch_dir = Path(tmpdir)
            report_path = batch_dir / "target_human_report.json"
            report_path.write_text(
                json.dumps(
                    {
                        "target": {
                            "accession": "P00001",
                            "entry_name": "TARGET_HUMAN",
                            "gene_symbol": "TARGET",
                            "sequence": "M",
                        },
                        "notes": [
                            {
                                "severity": "info",
                                "message": "Using a local AlphaFold 3 model.",
                                "source": "AlphaFold Server",
                            }
                        ],
                    }
                ),
                encoding="utf-8",
            )
            summary_path = batch_dir / "batch_summary.json"
            summary_path.write_text(
                json.dumps(
                    [
                        {
                            "query": "TARGET_HUMAN",
                            "status": "ok",
                            "json_report": str(report_path),
                        }
                    ]
                ),
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "AF3-derived report"):
                build_portal(summary_path, fetch_literature=False)

    def test_load_batch_rows_prefers_gene_then_uniprot_name_then_accession(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "targets.tsv"
            path.write_text(
                "gene\tuniprot_name\tuniprot_id\n"
                "EGFR\tEGFR_HUMAN\tP00533\n"
                "\tITGAV_HUMAN\tP06756\n"
                "\t\tQ9NRM6\n",
                encoding="utf-8",
            )
            rows = load_batch_rows(path)
        self.assertEqual([(row.query, row.source_column) for row in rows], [("EGFR", "gene"), ("ITGAV_HUMAN", "uniprot_name"), ("Q9NRM6", "uniprot_id")])

    def test_load_batch_rows_accepts_canonical_csv_columns(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "targets.csv"
            path.write_text(
                "accession,uniprot_entry_name,gene_symbol\n"
                "P00533,EGFR_HUMAN,EGFR\n"
                "P06756,ITAV_HUMAN,\n",
                encoding="utf-8",
            )
            rows = load_batch_rows(path)

        self.assertEqual([(row.query, row.source_column) for row in rows], [("EGFR_HUMAN", "uniprot_name"), ("ITAV_HUMAN", "uniprot_name")])
        self.assertEqual(rows[0].record["uniprot_id"], "P00533")
        self.assertEqual(rows[0].record["uniprot_name"], "EGFR_HUMAN")

    def test_canonical_csv_keeps_distinct_entries_with_the_same_gene(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "targets.csv"
            path.write_text(
                "accession,uniprot_entry_name,gene_symbol\n"
                "P11111,NRX1A_HUMAN,NRXN1\n"
                "P22222,NRX1B_HUMAN,NRXN1\n",
                encoding="utf-8",
            )
            rows = load_batch_rows(path)

        self.assertEqual([row.query for row in rows], ["NRX1A_HUMAN", "NRX1B_HUMAN"])

    def test_run_batch_from_tsv_writes_summary(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            tsv_path = tmpdir_path / "targets.tsv"
            out_dir = tmpdir_path / "out"
            tsv_path.write_text("uniprot_name\nEGFR_HUMAN\n", encoding="utf-8")
            results = run_batch_from_tsv(
                analyzer=FakeAnalyzer(),
                tsv_path=tsv_path,
                output_dir=out_dir,
                limit=None,
                resume=True,
            )
            summary = json.loads((out_dir / "batch_summary.json").read_text(encoding="utf-8"))
            self.assertTrue((out_dir / "egfr_human_report.json").exists())
            self.assertTrue((out_dir / "batch_summary.md").exists())
        self.assertEqual(len(results), 1)
        self.assertEqual(results[0]["status"], "ok")
        self.assertEqual(summary[0]["resolved_entry_name"], "EGFR_HUMAN")

    def test_render_batch_summary_markdown(self) -> None:
        markdown = render_batch_summary_markdown(
            [
                {
                    "query": "EGFR_HUMAN",
                    "batch_index": 1,
                    "batch_total": 3,
                    "status": "ok",
                    "resolved_entry_name": "EGFR_HUMAN",
                    "markdown_report": "outputs/batch/egfr_human_report.md",
                    "duration_seconds": 12.3,
                    "error": None,
                },
                {
                    "query": "ITGAV_HUMAN",
                    "batch_index": 2,
                    "batch_total": 3,
                    "status": "error",
                    "resolved_entry_name": None,
                    "markdown_report": None,
                    "duration_seconds": 1.2,
                    "error": "Connection error",
                },
            ]
        )
        self.assertIn("# Batch Summary", markdown)
        self.assertIn("Progress: 2/3", markdown)
        self.assertIn("`EGFR_HUMAN`", markdown)
        self.assertIn("Recent Errors", markdown)
        self.assertIn("Connection error", markdown)

    def test_retry_batch_from_tsv_reruns_errors_and_missing_alphafold(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            tsv_path = tmpdir_path / "targets.tsv"
            out_dir = tmpdir_path / "out"
            tsv_path.write_text("uniprot_name\nEGFR_HUMAN\nITGAV_HUMAN\n", encoding="utf-8")

            existing_report = out_dir / "egfr_human_report.json"
            existing_markdown = out_dir / "egfr_human_report.md"
            out_dir.mkdir(parents=True, exist_ok=True)
            existing_report.write_text(
                json.dumps(
                    {
                        "target": {
                            "accession": "P00001",
                            "entry_name": "EGFR_HUMAN",
                            "gene_symbol": "EGFR",
                            "protein_name": "EGFR",
                            "organism": "Homo sapiens",
                            "taxon_id": 9606,
                            "sequence": "M" * 100,
                        }
                    }
                ),
                encoding="utf-8",
            )
            existing_markdown.write_text("# Existing\n", encoding="utf-8")
            summary_path = out_dir / "batch_summary.json"
            summary_path.write_text(
                json.dumps(
                    [
                        {
                            "query": "EGFR_HUMAN",
                            "source_column": "uniprot_name",
                            "batch_index": 1,
                            "batch_total": 2,
                            "status": "ok",
                            "resolved_entry_name": "EGFR_HUMAN",
                            "json_report": str(existing_report),
                            "markdown_report": str(existing_markdown),
                            "duration_seconds": 1.0,
                            "error": None,
                        },
                        {
                            "query": "ITGAV_HUMAN",
                            "source_column": "uniprot_name",
                            "batch_index": 2,
                            "batch_total": 2,
                            "status": "error",
                            "resolved_entry_name": None,
                            "json_report": None,
                            "markdown_report": None,
                            "duration_seconds": 2.0,
                            "error": "boom",
                        },
                    ]
                ),
                encoding="utf-8",
            )

            analyzer = CountingAnalyzer()
            results = retry_batch_from_tsv(
                analyzer=analyzer,
                tsv_path=tsv_path,
                output_dir=out_dir,
                summary_path=summary_path,
                require_alphafold=True,
            )

        self.assertEqual(analyzer.calls, ["EGFR_HUMAN", "ITGAV_HUMAN"])
        self.assertEqual(len(results), 2)
        self.assertTrue(all(item["status"] == "ok" for item in results))

    def test_write_batch_summary_markdown(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            summary_path = tmpdir_path / "batch_summary.json"
            summary_path.write_text(
                json.dumps(
                    [
                        {
                            "query": "EGFR_HUMAN",
                            "batch_index": 1,
                            "batch_total": 1,
                            "status": "ok",
                            "resolved_entry_name": "EGFR_HUMAN",
                            "markdown_report": "outputs/batch/egfr_human_report.md",
                            "duration_seconds": 10.0,
                            "error": None,
                        }
                    ]
                ),
                encoding="utf-8",
            )
            output_path = write_batch_summary_markdown(summary_path)
            markdown = output_path.read_text(encoding="utf-8")
        self.assertIn("EGFR_HUMAN", markdown)

    def test_rebuild_batch_summary_from_tsv_restores_full_input_order(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            tsv_path = tmpdir_path / "targets.tsv"
            out_dir = tmpdir_path / "out"
            out_dir.mkdir()
            tsv_path.write_text(
                "uniprot_id\tuniprot_name\n"
                "P00533\tEGFR_HUMAN\n"
                "P06756\tITAV_HUMAN\n"
                "P07766\tCD3E_HUMAN\n",
                encoding="utf-8",
            )
            for query in ("EGFR_HUMAN", "ITAV_HUMAN"):
                FakeAnalyzer().analyze_target(query, output_dir=out_dir)

            results = rebuild_batch_summary_from_tsv(
                tsv_path=tsv_path,
                output_dir=out_dir,
            )
            summary = json.loads((out_dir / "batch_summary.json").read_text(encoding="utf-8"))

        self.assertEqual(len(results), 3)
        self.assertEqual([item["query"] for item in summary], ["EGFR_HUMAN", "ITAV_HUMAN", "CD3E_HUMAN"])
        self.assertEqual(summary[0]["status"], "ok")
        self.assertEqual(summary[1]["status"], "ok")
        self.assertEqual(summary[2]["status"], "error")
        self.assertIn("Missing report files", str(summary[2]["error"]))

    def test_rebuild_batch_summary_uses_refresh_state_for_pending_rows(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            tsv_path = tmpdir_path / "targets.tsv"
            out_dir = tmpdir_path / "out"
            out_dir.mkdir()
            tsv_path.write_text(
                "uniprot_id\tuniprot_name\n"
                "P00533\tEGFR_HUMAN\n"
                "P07766\tCD3E_HUMAN\n",
                encoding="utf-8",
            )
            FakeAnalyzer().analyze_target("EGFR_HUMAN", output_dir=out_dir)
            (out_dir / "accessible_portal_refresh_state.json").write_text(
                json.dumps(
                    {
                        "items": [
                            {
                                "query": "EGFR_HUMAN",
                                "source_column": "uniprot_name",
                                "status": "ok",
                                "resolved_entry_name": "EGFR_HUMAN",
                                "json_report": str(out_dir / "egfr_human_report.json"),
                                "markdown_report": str(out_dir / "egfr_human_report.md"),
                                "error": None,
                            },
                            {
                                "query": "CD3E_HUMAN",
                                "source_column": "uniprot_name",
                                "status": "pending",
                                "resolved_entry_name": None,
                                "json_report": None,
                                "markdown_report": None,
                                "error": None,
                            },
                        ]
                    }
                ),
                encoding="utf-8",
            )

            summary = rebuild_batch_summary_from_tsv(
                tsv_path=tsv_path,
                output_dir=out_dir,
            )

        self.assertEqual(summary[0]["status"], "ok")
        self.assertEqual(summary[1]["status"], "pending")
        self.assertIsNone(summary[1]["error"])

    def test_build_portal_creates_index_and_detail_pages(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            batch_dir = tmpdir_path / "batch"
            batch_dir.mkdir()
            report_json = batch_dir / "egfr_human_report.json"
            report_md = batch_dir / "egfr_human_report.md"
            assets_dir = batch_dir / "egfr_human_report_assets"
            assets_dir.mkdir()
            (assets_dir / "full_ectodomain_structure.png").write_bytes(b"png")
            (assets_dir / "full_ectodomain_quality.png").write_bytes(b"png")
            alphafold_dir = batch_dir / ".agdesign2" / "data" / "alphafold"
            alphafold_dir.mkdir(parents=True)
            (alphafold_dir / "P00533.pdb").write_text(
                "".join(
                    f"ATOM  {index:5d}  CA  MET A{index:4d}    {float(index):8.3f}{13.207:8.3f}{2.100:8.3f}  1.00 91.00           C\n"
                    for index in range(1, 101)
                ) + "END\n",
                encoding="utf-8",
            )
            report_json.write_text(
                json.dumps(
                    {
                        "target": {
                            "accession": "P00533",
                            "entry_name": "EGFR_HUMAN",
                            "gene_symbol": "EGFR",
                            "protein_name": "Epidermal growth factor receptor",
                            "gene_synonyms": ["ERBB1"],
                            "organism": "Homo sapiens",
                            "taxon_id": 9606,
                            "sequence": "M" * 100,
                        },
                        "ectodomain": {"start": 25, "end": 645, "label": "ecto", "source": "topology"},
                        "construct_details": [
                            {
                                "name": "full_ectodomain",
                                "start": 25,
                                "end": 645,
                                "length": 621,
                                "sequence": "A" * 50,
                                "score": 93.0,
                                "rationale": "Complete ectodomain.",
                                "classification": "multi_domain_unit",
                                "warnings": [],
                                "evidence": ["Topology-derived ectodomain"],
                                "human_entry_name": "EGFR_HUMAN",
                                "domain_summary": "Multiple domains",
                                "homologs": [
                                    {
                                        "species": "mouse",
                                        "accession": "Q01279",
                                        "entry_name": "EGFR_MOUSE",
                                        "available": True,
                                        "start": 25,
                                        "end": 647,
                                        "sequence": "MOUSESEQ",
                                        "identity_to_human": 88.73,
                                    },
                                    {
                                        "species": "macaca_fascicularis",
                                        "accession": "XP_005549616",
                                        "entry_name": "XP_005549616",
                                        "available": True,
                                        "start": 25,
                                        "end": 645,
                                        "sequence": "MONKEYSEQ",
                                        "identity_to_human": 98.71,
                                    },
                                ],
                                "cysteine_analysis": [
                                    {"position": 307, "paired_with": None, "surface_exposed": True, "warning": "Unpaired cysteine appears surface exposed."}
                                ],
                                "ptms": [
                                    {
                                        "type": "CARBOHYD",
                                        "start": 32,
                                        "end": 32,
                                        "description": "N-linked glycosylation",
                                        "metadata": {"ptm_category": "glycosylation"},
                                    },
                                    {
                                        "type": "SITE",
                                        "start": 45,
                                        "end": 45,
                                        "description": "Cleavage site",
                                        "metadata": {"ptm_category": "cleavage site"},
                                    },
                                ],
                                "structural_metrics": {"mean_plddt": 88.1, "structured_fraction": 0.91, "mean_intra_pae": 6.2},
                                "structure_image": None,
                                "quality_plot": None,
                            }
                        ],
                        "construct_recommendations": [],
                        "cross_reactivity_hits": [
                            {
                                "subject_id": "sp|Q01279|EGFR_MOUSE",
                                "description": "sp|Q01279|EGFR_MOUSE Epidermal growth factor receptor OS=Mus musculus OX=10090 GN=Egfr",
                                "species": "Mus musculus",
                                "identity": 88.7,
                                "coverage": 100.0,
                                "alignment_length": 56,
                                "evalue": 1e-30,
                                "bitscore": 250.0,
                                "query_start": 1,
                                "query_end": 56,
                                "subject_start": 5,
                                "subject_end": 60,
                                "query_alignment": "ACDEFGHIKLMNPQRSTVWY",
                                "subject_alignment": "ACD-FGHIKLMNPQ-STVWY",
                                "alignment_source": "BLAST",
                            }
                        ],
                        "furin_sites": [
                            {"start": 40, "end": 43, "motif": "RHRR", "suggestions": ["R40A", "R43A"]}
                        ],
                        "ptms": [
                            {
                                "type": "CARBOHYD",
                                "start": 32,
                                "end": 32,
                                "description": "N-linked glycosylation",
                                "metadata": {"ptm_category": "glycosylation"},
                            },
                            {
                                "type": "SITE",
                                "start": 45,
                                "end": 45,
                                "description": "Cleavage site",
                                "metadata": {"ptm_category": "cleavage site"},
                            },
                        ],
                        "cysteine_analysis": [
                            {"position": 31, "paired_with": 58, "surface_exposed": False, "warning": None},
                            {"position": 58, "paired_with": 31, "surface_exposed": False, "warning": None},
                            {"position": 307, "paired_with": None, "surface_exposed": True, "warning": "Unpaired cysteine appears surface exposed."}
                        ],
                        "ectodomain_homology": [],
                        "family_context": {
                            "source": "HGNC",
                            "family_names": ["ERBB family"],
                            "members": [{"gene_symbol": "EGFR", "entry_name": "EGFR_HUMAN", "ectodomain_start": 25, "ectodomain_end": 645}],
                            "identity_matrix_labels": ["EGFR"],
                            "identity_matrix": [[100.0]],
                        },
                        "assembly_requirements": [
                            {
                                "classification": "heteromeric_assembly",
                                "obligatory": False,
                                "confidence": "medium",
                                "summary": "Can heterodimerize with ERBB2.",
                                "partners": ["ERBB2"],
                            }
                        ],
                        "advanced_membrane_suggestions": [
                            {
                                "category": "gpcr_loop_engineering",
                                "title": "ICL-focused loop engineering screen",
                                "summary": "Consider BRIL or T4L insertion into a long internal cytoplasmic loop.",
                                "rationale": ["Longest internal cytoplasmic loop is suitable for engineering."],
                                "suggested_actions": ["Test BRIL insertion.", "Test T4 lysozyme insertion."],
                                "evidence": ["Topology-derived internal cytoplasmic loop."],
                                "warnings": ["Advisory only."],
                                "start": 250,
                                "end": 310,
                            }
                        ],
                        "experimental_constructs": [
                            {
                                "pdb_id": "1ABC",
                                "method": "X-ray",
                                "resolution": "2.0 A",
                                "chains": [{"label": "Chain A", "start": 25, "end": 525}],
                            }
                        ],
                        "interpro_annotations": [
                            {
                                "start": 57,
                                "end": 167,
                                "source_database": "INTERPRO",
                                "accession": "IPR000719",
                                "type": "domain",
                                "name": "EGF receptor L domain",
                            }
                        ],
                        "notes": [],
                    }
                ),
                encoding="utf-8",
            )
            report_md.write_text("# EGFR report\n", encoding="utf-8")
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
                            "markdown_report": str(report_md),
                            "duration_seconds": 10.0,
                            "error": None,
                        }
                    ]
                ),
                encoding="utf-8",
            )
            (batch_dir / "open_targets_disease_associations.json").write_text(
                json.dumps(
                    {
                        "targets": [
                            {
                                "entry_name": "EGFR_HUMAN",
                                "gene_symbol": "EGFR",
                                "uniprot_accession": "P00533",
                                "ensembl_id": "ENSG00000146648",
                                "status": "ok",
                                "association_count_downloaded": 2,
                                "top_disease_id": "EFO_0000616",
                                "top_disease_name": "neoplasm",
                                "top_disease_score": 0.91,
                                "open_targets_url": "https://platform.opentargets.org/target/ENSG00000146648/associations",
                            }
                        ],
                        "associations": [
                            {
                                "entry_name": "EGFR_HUMAN",
                                "gene_symbol": "EGFR",
                                "ensembl_id": "ENSG00000146648",
                                "rank": 1,
                                "disease_id": "EFO_0000616",
                                "disease_name": "neoplasm",
                                "score": 0.91,
                                "datasource_scores": "chembl:0.8",
                                "open_targets_url": "https://platform.opentargets.org/evidence/ENSG00000146648/EFO_0000616",
                            },
                            {
                                "entry_name": "EGFR_HUMAN",
                                "gene_symbol": "EGFR",
                                "ensembl_id": "ENSG00000146648",
                                "rank": 2,
                                "disease_id": "MONDO_0007254",
                                "disease_name": "breast cancer",
                                "score": 0.72,
                                "datasource_scores": "europepmc:0.6",
                                "open_targets_url": "https://platform.opentargets.org/evidence/ENSG00000146648/MONDO_0007254",
                            },
                        ],
                    }
                ),
                encoding="utf-8",
            )
            with mock.patch(
                "agdesign2.portal._get_pubtator_literature_info",
                return_value={
                    "gene_symbol": "EGFR",
                    "query_text": "@GENE_EGFR",
                    "query_url": "https://www.ncbi.nlm.nih.gov/research/pubtator3/docsum?text=%40GENE_EGFR",
                    "count": 204769,
                },
            ):
                index_path = build_portal(summary_path)
                index_html = index_path.read_text(encoding="utf-8")
                index_js = (index_path.parent / "portal-index.js").read_text(encoding="utf-8")
                index_data_js = (index_path.parent / "portal-index-data.js").read_text(encoding="utf-8")
                disease_index_js = (index_path.parent / "portal-disease-index.js").read_text(encoding="utf-8")
                detail_html = (index_path.parent / "reports" / "egfr_human.html").read_text(encoding="utf-8")
                detail_js = (index_path.parent / "report_scripts" / "egfr_human.js").read_text(encoding="utf-8")
                report_viewer_js = (index_path.parent / "report-viewer.js").read_text(encoding="utf-8")
                help_html = (index_path.parent / "help.html").read_text(encoding="utf-8")
                builder_html = (index_path.parent / "builder.html").read_text(encoding="utf-8")
                constructs_html = (index_path.parent / "constructs.html").read_text(encoding="utf-8")
                methods_html = (index_path.parent / "methods.html").read_text(encoding="utf-8")
                downloads_html = (index_path.parent / "downloads.html").read_text(encoding="utf-8")
                calculator_html = (index_path.parent / "calculator.html").read_text(encoding="utf-8")
                terms_html = (index_path.parent / "terms.html").read_text(encoding="utf-8")
                privacy_html = (index_path.parent / "privacy.html").read_text(encoding="utf-8")
                portal_css = (index_path.parent / "portal.css").read_text(encoding="utf-8")
                download_tsv = (index_path.parent / "downloads" / "agdesign2_portal_index.tsv").read_text(encoding="utf-8")
                download_json = json.loads((index_path.parent / "downloads" / "agdesign2_portal_index.json").read_text(encoding="utf-8"))
                download_manifest = json.loads((index_path.parent / "downloads" / "download_manifest.json").read_text(encoding="utf-8"))
                open_targets_download_json_exists = (index_path.parent / "open_targets_disease_associations.json").exists()
                favicon_ico_exists = (index_path.parent / "favicon.ico").exists()
                favicon_png_exists = (index_path.parent / "favicon-32.png").exists()
                apple_touch_icon_exists = (index_path.parent / "apple-touch-icon.png").exists()
                webmanifest_exists = (index_path.parent / "site.webmanifest").exists()
        self.assertIn("OpenAntigens", index_html)
        self.assertIn("Structure-guided antigen construct design", index_html)
        self.assertIn("cell-surface and secreted proteins", index_html)
        self.assertIn("Completed reports", index_html)
        self.assertIn("Secreted", index_html)
        self.assertIn("Single-pass", index_html)
        self.assertIn("GPI-anchored", index_html)
        self.assertIn("site-footer", index_html)
        self.assertIn("Created by Andre A. R. Teixeira", index_html)
        self.assertIn("Search gene, entry, protein, alias, family, or topology", index_html)
        self.assertIn('id="pageSizeSelect"', index_html)
        self.assertIn('<option value="25" selected>25</option>', index_html)
        self.assertIn('src="portal-index.js?v=', index_html)
        self.assertNotIn("let pageSize = Number(pageSizeSelect.value) || 25", index_html)
        # Rows stay as lightweight data objects; <tr> is built only for the
        # visible page (renderCurrentPage), not for all rows up front.
        self.assertIn("const rows = (window.OpenAntigenIndexRows || []).map(enrichPayload)", index_js)
        self.assertIn("filteredRows.slice(start, end).map(buildRowElement)", index_js)
        self.assertIn("function buildRowElement(payload)", index_js)
        self.assertIn("let pageSize = Number(pageSizeSelect.value) || 25", index_js)
        self.assertIn("pageSizeSelect.addEventListener('change'", index_js)
        self.assertIn("${d.aliases || ''} ${d.family || ''} ${d.track || ''}", index_js)
        self.assertNotIn("row.dataset.family} ${row.dataset.diseases} ${row.dataset.track}", index_js)
        self.assertIn("ipi-logo-light-800.png", index_html)
        self.assertIn("builder.html", index_html)
        self.assertIn("calculator.html", index_html)
        self.assertIn("help.html", index_html)
        self.assertIn("methods.html", index_html)
        self.assertIn("downloads.html", index_html)
        self.assertIn("terms.html", index_html)
        self.assertIn("privacy.html", index_html)
        self.assertIn("portal-index-data.js", index_html)
        self.assertIn('rel="icon" href="favicon.ico"', index_html)
        self.assertIn('rel="apple-touch-icon" sizes="180x180" href="apple-touch-icon.png?v=', index_html)
        self.assertTrue(favicon_ico_exists)
        self.assertTrue(favicon_png_exists)
        self.assertTrue(apple_touch_icon_exists)
        self.assertTrue(webmanifest_exists)
        self.assertIn("EGFR_HUMAN", index_data_js)
        self.assertIn('"aliases":"erbb1"', index_data_js)
        self.assertIn('class="sort-header active" data-sort-key="pubtator"', index_html)
        self.assertIn("Sorted by ${labels[sortState.key] || sortState.key}", index_js)
        self.assertNotIn('id="hasAlphaFoldFilter"', index_html)
        self.assertNotIn("hides targets without a matching canonical AlphaFold structure by default", index_html)
        self.assertNotIn("Status <span", index_html)
        self.assertNotIn("Seconds <span", index_html)
        self.assertNotIn("<th>Error</th>", index_html)
        self.assertNotIn('id="statusFilter"', index_html)
        self.assertNotIn("matchesToggle", index_js)
        self.assertIn("PubTator hits <span class=\"sort-arrow\">↓</span>", index_html)
        self.assertIn("AlphaFold <span class=\"sort-arrow\">↕</span>", index_html)
        self.assertIn("PDB <span class=\"sort-arrow\">↕</span>", index_html)
        self.assertIn("Entry <span class=\"sort-arrow\">↕</span>", index_html)
        self.assertIn("Top disease <span class=\"sort-arrow\">↕</span>", index_html)
        self.assertIn("Disease-focused filter", index_html)
        self.assertIn("neoplasm", index_data_js)  # top-disease cell stays embedded
        # The disease-name filter haystack moved out of the row payload into a
        # deferred script, merged client-side; the index no longer fetches the
        # full (multi-hundred-MB) Open Targets association export.
        self.assertNotIn('"diseases":', index_data_js)
        self.assertNotIn("breast cancer", index_data_js)
        self.assertIn("breast cancer", disease_index_js)
        self.assertIn("window.OpenAntigenDiseaseIndex", disease_index_js)
        self.assertIn('defer src="portal-disease-index.js', index_html)
        self.assertIn("OpenAntigenDiseaseIndex", index_js)
        self.assertNotIn("open_targets_disease_associations.json", index_js)
        self.assertIn(">204,769<", index_data_js)
        self.assertIn("pubtator3/docsum?text=%40GENE_EGFR", index_data_js)
        self.assertNotIn("Any structure", index_html)
        self.assertIn("AlphaFold", index_html)
        self.assertIn('"hasAlphafold":"1"', index_data_js)
        self.assertIn('"hasStructure":"1"', index_data_js)
        self.assertNotIn("Family context", index_html)
        self.assertNotIn("visible /", index_js)
        self.assertIn("Constructs", detail_html)
        self.assertIn("<h2>Construct Builder</h2>", detail_html)
        self.assertIn('id="interactive-construct-builder"', detail_html)
        self.assertIn("scrollIntoView({ behavior: 'smooth', block: 'start' })", report_viewer_js)
        self.assertIn("../ipi-logo-light-800.png", detail_html)
        self.assertIn('rel="icon" href="../favicon.ico"', detail_html)
        self.assertIn('rel="apple-touch-icon" sizes="180x180" href="../apple-touch-icon.png?v=', detail_html)
        self.assertIn("3Dmol-min.js", detail_html)
        self.assertNotIn("3Dmol.org/build", detail_html)
        self.assertIn("cdn.jsdelivr.net/npm/3dmol@2.5.5", detail_html)
        self.assertIn("sha256-98x4khrnLnYj6JzdERQ09Ywu/d0v/aHNISZEtAb7gBY=", detail_html)
        self.assertIn("../report-viewer.js?v=", detail_html)
        self.assertIn('../report_scripts/egfr_human.js', detail_html)
        self.assertNotIn("const sequence = ", detail_html)
        self.assertIn("const sequence = ", detail_js)
        self.assertIn("const viewerData = ", detail_js)
        self.assertIn("const targetInfo = ", detail_js)
        self.assertIn('entryName: "EGFR_HUMAN"', detail_js)
        self.assertIn('geneSymbol: "EGFR"', detail_js)
        self.assertIn('nameSymbol: "EGFR"', detail_js)
        self.assertIn("OpenAntigenReportViewer?.init", detail_js)
        self.assertIn("window.OpenAntigenReportViewer.init = function(payload)", report_viewer_js)
        self.assertIn("scheduleConstructWorkbench();", report_viewer_js)
        self.assertIn('document.addEventListener("DOMContentLoaded", installConstructWorkbench', report_viewer_js)
        self.assertIn("function installConstructWorkbench()", report_viewer_js)
        self.assertIn('builder.classList.add("construct-workbench")', report_viewer_js)
        self.assertIn('["Construct", "Variant"]', report_viewer_js)
        self.assertIn('["Boundary", "Focus region"]', report_viewer_js)
        self.assertIn("const targetInfo = payload.targetInfo || {}", report_viewer_js)
        self.assertIn("targetInfo.nameSymbol", report_viewer_js)
        self.assertNotIn("EGFR_HUMAN", report_viewer_js)
        self.assertIn("View in 3D", detail_html)
        self.assertIn("construct-highlight", detail_html)
        self.assertIn('id="viewerDownloadPdb"', detail_html)
        self.assertIn('id="viewerDownloadSelectedPdb"', detail_html)
        self.assertIn('id="viewerScreenshotPng"', detail_html)
        self.assertIn("function downloadSelectedPdb()", report_viewer_js)
        self.assertIn("function downloadViewerScreenshot()", report_viewer_js)
        self.assertIn("preserveDrawingBuffer: true", report_viewer_js)
        self.assertIn("_3d_view_${exportWidth}x${exportHeight}.png", report_viewer_js)
        self.assertIn("Interactive pLDDT", detail_html)
        self.assertIn("Interactive PAE", detail_html)
        self.assertIn("Live Selected Region", detail_html)
        self.assertIn("Copy TSV", detail_html)
        self.assertIn("Copy FASTA", detail_html)
        self.assertIn(".construct-workbench-shell", portal_css)
        self.assertIn("grid-template-columns: 190px minmax(560px, 1fr) minmax(312px, 330px);", portal_css)
        self.assertIn("flex-basis: auto;", portal_css)
        self.assertIn("max-height: 260px;", portal_css)
        self.assertNotIn("body.construct-workbench-active .page.detail", portal_css)
        self.assertIn("padding: 5px 6px;", portal_css)
        for event_name in (
            "Download PDB Full",
            "Download PDB Selection",
            "Download PNG",
            "Copy TSV",
            "Copy FASTA",
        ):
            self.assertIn(f"trackOpenAntigensEvent('{event_name}')", report_viewer_js)
        self.assertIn("selectionFasta", report_viewer_js)
        self.assertIn("selectionWarnings", detail_html)
        self.assertIn("selection-mutation-checkbox", report_viewer_js)
        self.assertIn("Unpaired cysteine warning", report_viewer_js)
        self.assertIn("Candidate furin-like motif warning", report_viewer_js)
        self.assertIn("Optional furin-site edits", report_viewer_js)
        self.assertIn("construct-furin-mutation-checkbox", report_viewer_js)
        self.assertIn("selection-furin-mutation-checkbox", report_viewer_js)
        self.assertIn("R40A", detail_html)
        self.assertIn("R43A", detail_html)
        self.assertIn("HUMAN", detail_html)
        self.assertIn("MOUSE", detail_html)
        self.assertIn("MACFA", detail_html)
        self.assertIn("full_ectodomain", detail_html)
        self.assertIn("full_ectodomain_structure.png", detail_html)
        self.assertIn("EGFR_HUMAN_25-645", detail_html)
        self.assertIn("EGFR_MOUSE_25-647", detail_html)
        self.assertIn("EGFR_MACFA_25-645", detail_html)
        self.assertIn("MOUSESEQ", detail_html)
        self.assertIn("MONKEYSEQ", detail_html)
        self.assertIn("species=HUMAN boundary=25-645 identity_to_human=100.00%", detail_html)
        self.assertIn("construct-mutation-checkbox", report_viewer_js)
        self.assertIn("Optional Cys→Ser edits", report_viewer_js)
        self.assertIn("Cysteines", detail_html)
        self.assertIn("Candidate Furin-like Motifs", detail_html)
        self.assertIn("Post-translational Modifications", detail_html)
        self.assertIn("N-linked glycosylation", detail_html)
        self.assertIn("Included PTMs", detail_html)
        self.assertIn("PTMs included in selected region", report_viewer_js)
        self.assertIn("Family Context", detail_html)
        self.assertIn("ERBB family", detail_html)
        self.assertIn('class="grid two"', detail_html)
        self.assertIn('class="matrix-cell"', detail_html)
        self.assertIn('data-identity="100.00"', detail_html)
        self.assertIn("Sequence similarity and cross-reactivity context", detail_html)
        self.assertIn("ranked by bit score", detail_html)
        self.assertIn("<th>Rank</th>", detail_html)
        self.assertIn("<th>Bit score</th>", detail_html)
        self.assertIn("Alignment source:</strong> BLAST", detail_html)
        self.assertIn("human  ACDEFGHIKLMNPQRSTVWY", detail_html)
        self.assertIn("hit    ACD-FGHIKLMNPQ-STVWY", detail_html)
        self.assertIn("Interactions / Assembly", detail_html)
        self.assertIn("Membrane Engineering Test Variants", detail_html)
        self.assertIn("ICL-focused loop engineering screen", detail_html)
        self.assertIn("Test BRIL insertion.", detail_html)
        self.assertNotIn("Can heterodimerize with ERBB2.", detail_html)
        self.assertIn("Precomputed Construct Summary", detail_html)
        self.assertIn("Detailed construct evidence", detail_html)
        self.assertIn("#construct-details", detail_html)
        self.assertIn("<th>Action</th>", detail_html)
        self.assertIn("Soluble Extracellular-region Constructs", detail_html)
        self.assertNotIn("Native Membrane-Expression Constructs", detail_html)
        self.assertIn("Topology", detail_html)
        self.assertIn("Topology-derived design region", detail_html)
        self.assertIn("InterPro / Pfam", detail_html)
        self.assertIn("Construct Details", detail_html)
        self.assertIn("IPR000719", detail_html)
        self.assertIn("UniProt", detail_html)
        self.assertIn("HGNC", detail_html)
        self.assertIn("PubTator hits", detail_html)
        self.assertIn("CiteAb", detail_html)
        self.assertIn("citeab.com/antibodies/search?q=EGFR", detail_html)
        self.assertIn("Open Targets Disease Associations", detail_html)
        self.assertIn("chembl:0.8", detail_html)
        self.assertIn(">204,769<", detail_html)
        self.assertIn("../builder.html", detail_html)
        self.assertIn("../help.html", detail_html)
        self.assertIn("site-footer", detail_html)
        self.assertIn("How to use OpenAntigens reports", help_html)
        self.assertIn("Quick start", help_html)
        self.assertIn("Report section map", help_html)
        self.assertIn("Recommended construct-review order", help_html)
        self.assertIn("Common workflows", help_html)
        self.assertIn("Troubleshooting", help_html)
        self.assertIn("The Browse page defaults to all rows sorted by PubTator hits from highest to lowest", help_html)
        self.assertIn("structure/pLDDT/PAE panels when compatible local structure assets exist", help_html)
        self.assertIn("Treat scores as disease-association context", help_html)
        self.assertNotIn("structure filters such as", help_html)
        self.assertNotIn("ranked construct sections", help_html)
        self.assertNotIn("Need raw data?", help_html)
        self.assertIn("Example: EGFR", help_html)
        self.assertIn("Example: ITGAV", help_html)
        self.assertIn("cynomolgus monkey", help_html)
        self.assertIn("How to use the Interactive Construct Builder", builder_html)
        self.assertIn("How to interpret pLDDT", builder_html)
        self.assertIn("How to interpret PAE", builder_html)
        self.assertIn("Live Selected Region", builder_html)
        self.assertIn("Reports with a local AlphaFold model also include", builder_html)
        self.assertIn("Start from the construct classes present in the report", builder_html)
        self.assertIn("export-ready sequence records", builder_html)
        self.assertIn("Manual overrides", builder_html)
        self.assertIn("generates construct candidates from resolved design scope", constructs_html)
        self.assertIn("when a compatible AlphaFold model and PAE matrix are available", constructs_html)
        self.assertIn("Every deposited structure is listed individually", constructs_html)
        self.assertIn("export-ready TSV/FASTA", constructs_html)
        self.assertIn("How OpenAntigens generates portal annotations", methods_html)
        self.assertIn("Evidence class", methods_html)
        self.assertIn("Construct generation", methods_html)
        self.assertIn("Sequence canonicality", methods_html)
        self.assertIn("Family and paralog context", methods_html)
        self.assertIn("Update cadence and limitations", methods_html)
        self.assertNotIn("AlphaFold 3", methods_html)
        self.assertIn("Exact duplicate boundaries from different calculated or annotated sources are collapsed", methods_html)
        self.assertIn("Human, mouse, and cynomolgus monkey primary databases use reviewed UniProt protein sets", methods_html)
        self.assertNotIn("canonical RefSeq-derived sequences where available", methods_html)
        self.assertNotIn("The output is intended to support", methods_html)
        self.assertNotIn("Exact duplicate PDB boundaries", methods_html)
        self.assertIn("Download static release data", downloads_html)
        self.assertIn("OpenAntigens release snapshots include static flat files", downloads_html)
        self.assertIn("Downloads are generated files from the current release snapshot", downloads_html)
        self.assertIn("mouse_ortholog_*", downloads_html)
        self.assertIn("human_source_*", downloads_html)
        self.assertIn("Open Targets TSV/JSON files contain downloaded indirect and direct target-disease association scores", downloads_html)
        self.assertNotIn("users who want", downloads_html)
        self.assertNotIn("recommended starting point", downloads_html)
        self.assertIn('href="open_targets_disease_associations.json"', downloads_html)
        self.assertNotIn('href="../open_targets_disease_associations.json"', downloads_html)
        self.assertIn("Protein concentration calculator", calculator_html)
        self.assertIn('id="molecularWeightUnit"', calculator_html)
        self.assertIn('<option value="kDa" selected>kDa</option>', calculator_html)
        self.assertIn("Molecular weight is the only required field", calculator_html)
        self.assertIn("Equivalent value calculated", calculator_html)
        self.assertIn("massConcentrationUnit", calculator_html)
        self.assertIn("molarConcentrationUnit", calculator_html)
        self.assertIn("volumeUnit", calculator_html)
        self.assertIn("site-footer", calculator_html)
        self.assertIn("Terms for OpenAntigens data and software", terms_html)
        self.assertIn("Apache License 2.0", terms_html)
        self.assertIn("CC BY 4.0", terms_html)
        self.assertIn("Legal review and source-database compliance remain", terms_html)
        self.assertNotIn("AlphaFold 3", terms_html)
        self.assertNotIn("AlphaFold Server", terms_html)
        self.assertIn("Indirect and direct target-disease association scores", terms_html)
        self.assertIn("Complex Portal", terms_html)
        self.assertIn("Original source terms remain attached", terms_html)
        self.assertIn("OpenAntigens is provided as-is and without warranties", terms_html)
        self.assertNotIn("These terms are intended", terms_html)
        self.assertNotIn("not legal advice", terms_html)
        self.assertNotIn("not intended for diagnosis", terms_html)
        self.assertNotIn("not a replacement", terms_html)
        self.assertNotIn("does not redistribute CiteAb", terms_html)
        self.assertNotIn("Users are responsible", terms_html)
        self.assertIn("andre.teixeira@proteininnovation.org", terms_html)
        self.assertIn("UniProt license", terms_html)
        self.assertIn("AlphaFold DB", terms_html)
        self.assertIn("Privacy and Analytics", privacy_html)
        self.assertIn("Plausible Analytics", privacy_html)
        self.assertIn("does not use analytics cookies", privacy_html)
        self.assertIn("fewer than five unique visitors", privacy_html)
        for page_html in (index_html, calculator_html, help_html, detail_html, privacy_html):
            self.assertEqual(
                page_html.count('src="https://plausible.io/js/pa-Mda6DF7h4_8b17ZUlD_9E.js"'),
                1,
            )
            self.assertEqual(page_html.count("data-openantigens-analytics"), 1)
            self.assertIn('fileExtensions: ["tsv", "json"]', page_html)
            self.assertIn("formSubmissions: false", page_html)
            self.assertIn("outboundLinks: false", page_html)
        self.assertIn("entry_name", download_tsv)
        self.assertIn("EGFR_HUMAN", download_tsv)
        self.assertEqual(download_json[0]["entry_name"], "EGFR_HUMAN")
        self.assertEqual(download_manifest["row_count"], 1)
        self.assertIn("top_disease_name", download_manifest["files"][0]["columns"])
        manifest_paths = {item["path"] for item in download_manifest["files"]}
        self.assertIn("open_targets_disease_associations.json", manifest_paths)
        self.assertNotIn("../open_targets_disease_associations.json", manifest_paths)
        self.assertTrue(open_targets_download_json_exists)

    def test_build_portal_uses_snapshot_local_alphafold_and_reuses_public_structure(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            snapshot_dir = Path(tmpdir) / "snapshot"
            batch_dir = snapshot_dir / "outputs" / "surfy_batch"
            public_site = snapshot_dir / "public_site"
            alphafold_dir = snapshot_dir / "data" / "alphafold"
            batch_dir.mkdir(parents=True)
            alphafold_dir.mkdir(parents=True)
            report_json = batch_dir / "egfr_human_report.json"
            report_md = batch_dir / "egfr_human_report.md"
            pdb_path = alphafold_dir / "P00533.pdb"
            pdb_path.write_text(
                "".join(
                    f"ATOM  {index:5d}  CA  MET A{index:4d}    {float(index):8.3f}{13.207:8.3f}{2.100:8.3f}  1.00 91.00           C\n"
                    for index in range(1, 11)
                ) + "END\n",
                encoding="utf-8",
            )
            report_json.write_text(
                json.dumps(
                    {
                        "target": {
                            "accession": "P00533",
                            "entry_name": "EGFR_HUMAN",
                            "gene_symbol": "EGFR",
                            "protein_name": "Epidermal growth factor receptor",
                            "sequence": "M" * 10,
                        },
                        "ectodomain": {"start": 1, "end": 10, "label": "design region", "source": "test"},
                        "construct_details": [],
                        "construct_recommendations": [],
                        "cross_reactivity_hits": [],
                        "experimental_constructs": [],
                        "family_context": {},
                        "assembly_requirements": [],
                        "notes": [],
                    }
                ),
                encoding="utf-8",
            )
            report_md.write_text("# EGFR\n", encoding="utf-8")
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
                            "json_report": "/stale/original/snapshot/egfr_human_report.json",
                            "markdown_report": "/stale/original/snapshot/egfr_human_report.md",
                            "duration_seconds": 1.0,
                            "error": None,
                        }
                    ]
                ),
                encoding="utf-8",
            )

            first_index = build_portal(
                summary_path,
                output_dir=public_site,
                refresh_reports=False,
                refresh_stale_reports=False,
            )
            first_detail = (first_index.parent / "reports" / "egfr_human.html").read_text(encoding="utf-8")
            first_data = (first_index.parent / "portal-index-data.js").read_text(encoding="utf-8")

            self.assertTrue((public_site / "structures" / "egfr_human.pdb").exists())
            self.assertIn('id="viewerDownloadPdb"', first_detail)
            self.assertIn('id="viewerDownloadSelectedPdb"', first_detail)
            self.assertIn("<h2>Construct Builder</h2>", first_detail)
            self.assertIn('"hasAlphafold":"1"', first_data)
            self.assertNotIn("No local AlphaFold structure is available", first_detail)

            pdb_path.unlink()
            second_index = build_portal(
                summary_path,
                output_dir=public_site,
                refresh_reports=False,
                refresh_stale_reports=False,
            )
            second_detail = (second_index.parent / "reports" / "egfr_human.html").read_text(encoding="utf-8")
            second_data = (second_index.parent / "portal-index-data.js").read_text(encoding="utf-8")

            self.assertIn('id="viewerDownloadPdb"', second_detail)
            self.assertIn('id="viewerDownloadSelectedPdb"', second_detail)
            self.assertIn('"hasAlphafold":"1"', second_data)
            self.assertNotIn("No local AlphaFold structure is available", second_detail)

    def test_build_portal_replaces_invalid_public_alphafold_copy(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            snapshot_dir = Path(tmpdir) / "snapshot"
            batch_dir = snapshot_dir / "outputs" / "surfy_batch"
            public_site = snapshot_dir / "public_site"
            alphafold_dir = snapshot_dir / "data" / "alphafold"
            batch_dir.mkdir(parents=True)
            alphafold_dir.mkdir(parents=True)
            (public_site / "structures").mkdir(parents=True)
            report_json = batch_dir / "egfr_human_report.json"
            report_md = batch_dir / "egfr_human_report.md"
            pdb_path = alphafold_dir / "P00533.pdb"
            pdb_text = (
                "".join(
                    f"ATOM  {index:5d}  CA  MET A{index:4d}    {float(index):8.3f}{13.207:8.3f}{2.100:8.3f}  1.00 91.00           C\n"
                    for index in range(1, 11)
                )
                + "END\n"
            )
            pdb_path.write_text(pdb_text, encoding="utf-8")
            public_pdb = public_site / "structures" / "egfr_human.pdb"
            public_pdb.write_text("", encoding="utf-8")
            os.utime(public_pdb, (pdb_path.stat().st_mtime + 100, pdb_path.stat().st_mtime + 100))
            report_json.write_text(
                json.dumps(
                    {
                        "target": {
                            "accession": "P00533",
                            "entry_name": "EGFR_HUMAN",
                            "gene_symbol": "EGFR",
                            "protein_name": "Epidermal growth factor receptor",
                            "sequence": "M" * 10,
                        },
                        "ectodomain": {"start": 1, "end": 10, "label": "design region", "source": "test"},
                        "construct_details": [],
                        "construct_recommendations": [],
                        "cross_reactivity_hits": [],
                        "experimental_constructs": [],
                        "family_context": {},
                        "assembly_requirements": [],
                        "notes": [],
                    }
                ),
                encoding="utf-8",
            )
            report_md.write_text("# EGFR\n", encoding="utf-8")
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
                            "markdown_report": str(report_md),
                            "duration_seconds": 1.0,
                            "error": None,
                        }
                    ]
                ),
                encoding="utf-8",
            )

            index = build_portal(summary_path, output_dir=public_site, refresh_reports=False, refresh_stale_reports=False)
            detail = (index.parent / "reports" / "egfr_human.html").read_text(encoding="utf-8")

            self.assertEqual(public_pdb.read_text(encoding="utf-8"), pdb_text)
            self.assertIn('id="viewerDownloadPdb"', detail)
            self.assertNotIn("No matching canonical AlphaFold structure is available", detail)

    def test_single_detail_render_copies_snapshot_local_alphafold_when_entry_lacks_structure_path(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            snapshot_dir = Path(tmpdir) / "snapshot"
            batch_dir = snapshot_dir / "outputs" / "surfy_batch"
            public_site = snapshot_dir / "public_site"
            alphafold_dir = snapshot_dir / "data" / "alphafold"
            page_dir = public_site / "reports"
            batch_dir.mkdir(parents=True)
            alphafold_dir.mkdir(parents=True)
            page_dir.mkdir(parents=True)
            pdb_path = alphafold_dir / "P00533.pdb"
            pdb_path.write_text(
                "".join(
                    f"ATOM  {index:5d}  CA  MET A{index:4d}    {float(index):8.3f}{13.207:8.3f}{2.100:8.3f}  1.00 91.00           C\n"
                    for index in range(1, 11)
                )
                + "END\n",
                encoding="utf-8",
            )
            report = {
                "target": {
                    "accession": "P00533",
                    "entry_name": "EGFR_HUMAN",
                    "gene_symbol": "EGFR",
                    "protein_name": "Epidermal growth factor receptor",
                    "sequence": "M" * 10,
                },
                "ectodomain": {"start": 1, "end": 10, "label": "design region", "source": "test"},
                "construct_details": [],
                "construct_recommendations": [],
                "cross_reactivity_hits": [],
                "experimental_constructs": [],
                "family_context": {},
                "assembly_requirements": [],
                "notes": [],
            }
            entry = {
                "query": "EGFR_HUMAN",
                "entry_name": "EGFR_HUMAN",
                "detail_page": "egfr_human.html",
                "status": "ok",
            }

            detail = render_detail_page(entry, report, batch_dir=batch_dir, portal_dir=public_site)

            self.assertTrue((public_site / "structures" / "egfr_human.pdb").exists())
            self.assertIn('id="viewerDownloadPdb"', detail)
            self.assertIn('id="viewerDownloadSelectedPdb"', detail)
            self.assertIn("<h2>Construct Builder</h2>", detail)
            self.assertNotIn("No local AlphaFold structure is available", detail)

    def test_detail_page_renders_gpcr_annotation(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            batch_dir = Path(tmpdir) / "batch"
            portal_dir = Path(tmpdir) / "portal"
            batch_dir.mkdir()
            structures_dir = portal_dir / "structures"
            structures_dir.mkdir(parents=True)
            pdb_path = structures_dir / "oprm_human.pdb"
            pdb_path.write_text(
                "\n".join(
                    [
                        f"ATOM  {index:5d}  CA  MET A{index:4d}    {index * 1.0:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00{95.0:6.2f}           C"
                        for index in range(1, 401)
                    ]
                    + ["END"]
                ),
                encoding="utf-8",
            )
            report = {
                "target": {
                    "accession": "P35372",
                    "entry_name": "OPRM_HUMAN",
                    "gene_symbol": "OPRM1",
                    "protein_name": "Mu-type opioid receptor",
                    "sequence": "M" * 400,
                },
                "topology": {"topology_class": "multipass_compact"},
                "construct_details": [
                    {
                        "name": "membrane_expression_full_length",
                        "start": 1,
                        "end": 400,
                        "length": 400,
                        "sequence": "M" * 400,
                        "score": 86.0,
                        "rationale": "Native baseline.",
                        "classification": "membrane_expression",
                        "warnings": [],
                        "evidence": [],
                        "gpcr_segments": [{"name": "TM1", "start": 81, "end": 112, "generic_start": "1.43x43", "generic_end": "1.74x74"}],
                        "gpcr_generic_range": "1.43x43 to 7.53x53",
                    }
                ],
                "construct_recommendations": [],
                "cross_reactivity_hits": [],
                "full_length_cross_reactivity_hits": [],
                "experimental_constructs": [],
                "family_context": {},
                "assembly_requirements": [],
                "advanced_membrane_suggestions": [],
                "gpcr_engineering_variants": [
                    {
                        "name": "OPRM_HUMAN_BRIL_ICL3_fusion",
                        "category": "gpcr_loop_fusion",
                        "strategy": "ICL3 cassette replacement",
                        "summary": "Experimental BRIL fusion.",
                        "sequence": "AAAAAGGSGBBBBGGSGAAAAA",
                        "sequence_role": "engineered_receptor",
                        "cassette_id": "bril",
                        "cassette_name": "BRIL / apocytochrome b562RIL",
                        "start": 240,
                        "end": 275,
                        "replaced_start": 245,
                        "replaced_end": 270,
                        "n_linker": "GGSG",
                        "c_linker": "GGSG",
                        "engineered_length": 23,
                        "warnings": ["Experimental GPCR engineering variant."],
                        "evidence": ["GPCRdb ICL3: 240-275."],
                        "suggested_actions": ["Compare against native constructs."],
                        "tag_suggestions": {"c_terminal": ["Optional C-terminal tag metadata."]},
                    }
                ],
                "notes": [],
                "gpcr_annotation": {
                    "entry_name": "oprm_human",
                    "accession": "P35372",
                    "name": "mu receptor",
                    "family_path": ["Class A (Rhodopsin)", "Opioid receptors", "mu receptor"],
                    "residue_numbering_scheme": "GPCRdb(A)",
                    "url": "https://gpcrdb.org/protein/oprm_human",
                    "segments": [{"name": "TM1", "start": 81, "end": 112, "generic_start": "1.43x43", "generic_end": "1.74x74"}],
                    "conserved_motifs": [{"name": "DRY / E/DRY activation motif", "positions": [164, 165, 166], "generic_numbers": ["3.49x49", "3.50x50", "3.51x51"], "sequence": "DRY", "status": "detected"}],
                },
            }
            entry = {
                "entry_name": "OPRM_HUMAN",
                "detail_page": "oprm_human.html",
                "status": "ok",
                "portal_structure_path": str(pdb_path),
            }

            detail = render_detail_page(entry, report, batch_dir=batch_dir, portal_dir=portal_dir)

        self.assertIn("GPCR Annotation", detail)
        self.assertIn("oprm_human", detail)
        self.assertIn("Class A (Rhodopsin)", detail)
        self.assertIn("DRY / E/DRY activation motif", detail)
        self.assertIn("GPCR Segments", detail)
        self.assertIn("1.43x43 to 7.53x53", detail)
        self.assertIn("gpcrSegmentBands", detail)
        self.assertNotIn("gpcr-segment-track", detail)
        self.assertIn("viewerToggleSelectionOnly", detail)
        self.assertIn("showSelectionOnly", detail)
        self.assertIn("extracellularTopologyPositions", detail)
        self.assertIn("focusResidueSet", detail)
        self.assertIn("Extracellular topology residues", detail)
        self.assertIn("GPCR Expression and Stabilization Strategy", detail)
        self.assertIn("Experimental GPCR Engineering Variants", detail)
        self.assertIn("OPRM_HUMAN_BRIL_ICL3_fusion", detail)
        self.assertIn('["Group", "Category"]', detail)
        self.assertIn('["Boundary", "Focus region"]', detail)
        self.assertIn("Copy-ready exports", detail)
        self.assertIn("construct-export-fasta", detail)


class BatchWorkerReuseTests(unittest.TestCase):
    def test_worker_analyzer_is_built_once_and_reused(self) -> None:
        import agdesign2.batch as batch_module
        from unittest import mock

        constructed = {"n": 0}

        class CountingFakeAnalyzer(FakeAnalyzer):
            def __init__(self, config=None) -> None:
                constructed["n"] += 1

        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir)
            tsv = out / "t.tsv"
            tsv.write_text("uniprot_name\nAAA_HUMAN\nBBB_HUMAN\n", encoding="utf-8")
            rows = load_batch_rows(tsv)
            saved = batch_module._WORKER_ANALYZER
            batch_module._WORKER_ANALYZER = None
            try:
                with mock.patch("agdesign2.batch.AntigenAnalyzer", CountingFakeAnalyzer):
                    # Pool initializer builds exactly one analyzer per worker...
                    batch_module._init_batch_worker(AnalysisConfig())
                    # ...and every row that worker handles reuses it.
                    batch_module._run_batch_worker(AnalysisConfig(), rows[0], 1, 2, out)
                    batch_module._run_batch_worker(AnalysisConfig(), rows[1], 2, 2, out)
            finally:
                batch_module._WORKER_ANALYZER = saved

        self.assertEqual(constructed["n"], 1)


class PortalBuildFailureTests(unittest.TestCase):
    def test_duplicate_boundary_card_uses_shared_portal_quality_plot(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            page_dir = root / "portal" / "reports"
            assets_dir = (
                root
                / "portal"
                / "report_assets"
                / "test_human_report_assets"
            )
            page_dir.mkdir(parents=True)
            assets_dir.mkdir(parents=True)
            (assets_dir / "pdb_1abc_quality.png").write_bytes(b"png")

            html = _render_construct_card(
                {
                    "name": "pdb_2def",
                    "start": 2,
                    "end": 5,
                    "length": 4,
                    "sequence": "AAAA",
                    "score": 1.0,
                    "rationale": "test",
                    "quality_plot": "test_human_report_assets/pdb_1abc_quality.png",
                },
                target={"entry_name": "TEST_HUMAN"},
                batch_dir=root,
                page_dir=page_dir,
                furin_sites=[],
                ptms=[],
            )

            self.assertIn(
                "../report_assets/test_human_report_assets/pdb_1abc_quality.png",
                html,
            )

    def test_report_asset_copy_removes_stale_destination_files(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source = root / "test_human_report_assets"
            source.mkdir()
            report_assets_dir = root / "portal" / "report_assets"
            destination = report_assets_dir / source.name
            destination.mkdir(parents=True)
            stale = destination / "old_structure.png"
            stale.write_bytes(b"png")

            _copy_report_assets(
                report_data={"target": {"entry_name": "TEST_HUMAN"}},
                report_json_path=root / "test_human_report.json",
                report_assets_dir=report_assets_dir,
                batch_dir=root,
            )

            self.assertFalse(stale.exists())

    def test_report_asset_copy_rejects_symlinks(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source = root / "test_human_report_assets"
            source.mkdir()
            outside = root / "outside-secret.txt"
            outside.write_text("PRIVATE", encoding="utf-8")
            (source / "leak.txt").symlink_to(outside)

            with self.assertRaisesRegex(ValueError, "symbolic links"):
                _copy_report_assets(
                    report_data={"target": {"entry_name": "TEST_HUMAN"}},
                    report_json_path=root / "test_human_report.json",
                    report_assets_dir=root / "portal" / "report_assets",
                    batch_dir=root,
                )

    def test_vendor_bundle_rejects_unpinned_local_asset(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            portal_dir = root / "portal"
            portal_dir.mkdir()
            (portal_dir / "3Dmol-min.js").write_text("stale", encoding="utf-8")

            with self.assertRaisesRegex(RuntimeError, "pinned SHA-256"):
                _ensure_vendor_assets(portal_dir=portal_dir, batch_dir=root)

    def test_build_fails_on_invalid_report_json(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            report_path = root / "broken_report.json"
            report_path.write_text("{", encoding="utf-8")
            summary_path = root / "batch_summary.json"
            summary_path.write_text(
                json.dumps([{"query": "BROKEN_HUMAN", "json_report": str(report_path)}]),
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "Could not read report JSON"):
                build_portal(
                    summary_path,
                    refresh_reports=False,
                    refresh_stale_reports=False,
                    include_disease_context=False,
                )


if __name__ == "__main__":
    unittest.main()
