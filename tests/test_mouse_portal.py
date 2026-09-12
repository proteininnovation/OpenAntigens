from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from agdesign2.mouse_portal import (
    build_mouse_portal_from_human_orthologs,
    _require_blast_tools_for_mouse_cross_reactivity,
    _mouse_resume_decision,
    _mouse_pipeline_stamp,
    _record_mouse_pipeline_version,
    MOUSE_PIPELINE_VERSION,
)


class MouseResumeVersionTests(unittest.TestCase):
    """Guardrail: the mouse build must not resume onto reports left by a different
    pipeline version (the bug that left 0 cross-reactivity after the rewrite)."""

    def test_no_sentinel_disables_resume(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            reports = Path(tmpdir)
            (reports / "x_mouse.json").write_text("{}", encoding="utf-8")  # stale, unstamped
            eff, recorded = _mouse_resume_decision(reports, requested_resume=True)
            self.assertFalse(eff)
            self.assertIsNone(recorded)

    def test_matching_sentinel_allows_resume(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            reports = Path(tmpdir)
            _record_mouse_pipeline_version(reports)
            eff, recorded = _mouse_resume_decision(reports, requested_resume=True)
            self.assertTrue(eff)
            self.assertEqual(recorded, _mouse_pipeline_stamp(True))

    def test_mismatched_sentinel_disables_resume(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            reports = Path(tmpdir)
            (reports / ".mouse_pipeline_version").write_text("old-version\n", encoding="utf-8")
            eff, recorded = _mouse_resume_decision(reports, requested_resume=True)
            self.assertFalse(eff)
            self.assertEqual(recorded, "old-version")

    def test_no_resume_request_is_honored_even_with_matching_sentinel(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            reports = Path(tmpdir)
            _record_mouse_pipeline_version(reports)
            eff, _ = _mouse_resume_decision(reports, requested_resume=False)
            self.assertFalse(eff)

    def test_alphafold_fetch_mode_is_part_of_resume_compatibility(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            reports = Path(tmpdir)
            _record_mouse_pipeline_version(reports, fetch_mouse_alphafold=False)
            matching, _ = _mouse_resume_decision(
                reports,
                requested_resume=True,
                fetch_mouse_alphafold=False,
            )
            mismatched, _ = _mouse_resume_decision(
                reports,
                requested_resume=True,
                fetch_mouse_alphafold=True,
            )
            self.assertTrue(matching)
            self.assertFalse(mismatched)


class MousePortalTests(unittest.TestCase):
    def test_builds_mouse_reports_from_human_ortholog_constructs(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            human_dir = root / "human"
            human_dir.mkdir()
            report_path = human_dir / "egfr_human_report.json"
            report_path.write_text(json.dumps(_human_report()), encoding="utf-8")
            erbb2_report_path = human_dir / "erbb2_human_report.json"
            erbb2_report_path.write_text(json.dumps(_human_report("ERBB2_HUMAN", "P04626", "ERBB2")), encoding="utf-8")
            duplicate_report_path = human_dir / "egfr2_human_report.json"
            duplicate_report_path.write_text(json.dumps(_human_report("EGFR2_HUMAN", "P00001", "EGFR2")), encoding="utf-8")
            (human_dir / "batch_summary.json").write_text(
                json.dumps(
                    [
                        {
                            "batch_index": 1,
                            "query": "EGFR_HUMAN",
                            "status": "ok",
                            "json_report": str(report_path),
                            "topology_bucket": "Single-pass",
                        },
                        {
                            "batch_index": 2,
                            "query": "ERBB2_HUMAN",
                            "status": "ok",
                            "json_report": str(erbb2_report_path),
                            "topology_bucket": "Single-pass",
                        },
                        {
                            "batch_index": 3,
                            "query": "EGFR2_HUMAN",
                            "status": "ok",
                            "json_report": str(duplicate_report_path),
                            "topology_bucket": "Single-pass",
                        },
                        {
                            "batch_index": 4,
                            "query": "NO_MOUSE_HUMAN",
                            "status": "ok",
                            "json_report": str(human_dir / "missing_mouse_report.json"),
                        },
                    ]
                ),
                encoding="utf-8",
            )
            (human_dir / "ortholog_reference_table.tsv").write_text(
                "\t".join(
                    [
                        "input_uniprot_name",
                        "human_uniprot_accession",
                        "human_gene_symbol",
                        "mouse_gene_symbol",
                        "mouse_refseq_accession",
                        "mouse_refseq_sequence",
                    ]
                )
                + "\nEGFR_HUMAN\tP00533\tEGFR\tEgfr\tNP_034255\t"
                + "MOUSEFULLSEQ"
                + "\nERBB2_HUMAN\tP04626\tERBB2\tErbb2\tNP_001003817\t"
                + "ERBB2MOUSEFULLSEQ"
                + "\nEGFR2_HUMAN\tP00001\tEGFR2\tEgfr2\tNP_999999\t"
                + "DUPLICATEMOUSESEQ"
                + "\n",
                encoding="utf-8",
            )

            def fake_run_batch(*, analyzer, tsv_path, output_dir, resume, jobs):
                del tsv_path, resume, jobs
                # Per-target cross-reactivity is off during analyze; it is computed
                # once in bulk via refresh_report_modules (search_many) instead.
                self.assertFalse(analyzer.config.run_cross_reactivity)
                self.assertFalse(analyzer.config.fetch_alphafold)
                self.assertEqual(analyzer.config.target_species, "mouse")
                output_root = Path(output_dir)
                output_root.mkdir(parents=True, exist_ok=True)
                report_path = output_root / "egfr_mouse_report.json"
                markdown_path = output_root / "egfr_mouse_report.md"
                report = _mouse_analyzed_report()
                report_path.write_text(json.dumps(report), encoding="utf-8")
                markdown_path.write_text("# EGFR_MOUSE\n", encoding="utf-8")
                erbb2_report_path = output_root / "erbb2_mouse_report.json"
                erbb2_markdown_path = output_root / "erbb2_mouse_report.md"
                erbb2_report_path.write_text(json.dumps(_mouse_analyzed_report("ERBB2_MOUSE", "P70424", "Erbb2")), encoding="utf-8")
                erbb2_markdown_path.write_text("# ERBB2_MOUSE\n", encoding="utf-8")
                return [
                    {
                        "batch_index": 1,
                        "query": "EGFR_MOUSE",
                        "source_column": "uniprot_name",
                        "status": "ok",
                        "resolved_entry_name": "EGFR_MOUSE",
                        "json_report": str(report_path),
                        "markdown_report": str(markdown_path),
                        "duration_seconds": 1.0,
                        "mouse_ortholog_accession": "Q01279",
                        "mouse_refseq_accession": "NP_034255",
                        "mouse_sequence_source": "uniprot_reviewed_mouse_gene_symbol",
                        "human_source_entry_name": "EGFR_HUMAN",
                    },
                    {
                        "batch_index": 2,
                        "query": "ERBB2_MOUSE",
                        "source_column": "uniprot_name",
                        "status": "ok",
                        "resolved_entry_name": "ERBB2_MOUSE",
                        "json_report": str(erbb2_report_path),
                        "markdown_report": str(erbb2_markdown_path),
                        "duration_seconds": 1.0,
                        "mouse_ortholog_accession": "P70424",
                        "mouse_refseq_accession": "NP_001003817",
                        "mouse_sequence_source": "uniprot_reviewed_mouse_gene_symbol",
                        "human_source_entry_name": "ERBB2_HUMAN",
                    },
                ]

            with (
                patch(
                    "agdesign2.mouse_portal._resolve_mouse_uniprot_target",
                    side_effect=[
                        {
                            "accession": "Q01279",
                            "entry_name": "EGFR_MOUSE",
                            "sequence": "UNIPROTFULLSEQ",
                            "protein_name": "Epidermal growth factor receptor",
                        },
                        {
                            "accession": "P70424",
                            "entry_name": "ERBB2_MOUSE",
                            "sequence": "ERBB2UNIPROTFULLSEQ",
                            "protein_name": "Receptor tyrosine-protein kinase erbB-2",
                        },
                        {
                            "accession": "Q01279",
                            "entry_name": "Q01279",
                            "sequence": "UNIPROTFULLSEQ",
                            "protein_name": "Duplicate epidermal growth factor receptor",
                        },
                    ],
                ),
            patch("agdesign2.mouse_portal.shutil.which", return_value="/usr/bin/tool"),
            patch(
                "agdesign2.mouse_portal.prefetch_targets_from_tsv",
                return_value={},
            ) as mock_prefetch,
            patch("agdesign2.mouse_portal.run_batch_from_tsv", side_effect=fake_run_batch),
                patch("agdesign2.mouse_portal.refresh_report_modules", return_value=None) as mock_refresh,
            ):
                result = build_mouse_portal_from_human_orthologs(
                    human_dir / "batch_summary.json",
                    output_dir=root / "mouse",
                    build_site=True,
                    fetch_mouse_alphafold=False,
                    af3_catalog=root / "missing-af3.sqlite",
                )

            self.assertEqual(result.included_count, 2)
            self.assertEqual(result.skipped_count, 2)
            self.assertEqual(
                (root / "mouse" / "ortholog_reference_table.tsv").read_text(encoding="utf-8"),
                (human_dir / "ortholog_reference_table.tsv").read_text(encoding="utf-8"),
            )
            mock_prefetch.assert_called_once_with(
                tsv_path=root / "mouse" / "mouse_targets.tsv",
                cache_dir=root / "mouse" / ".agdesign2" / "cache",
            )
            # Cross-reactivity must go through the bulk search_many path, not per target.
            self.assertTrue(mock_refresh.called, "mouse build should run the bulk cross-reactivity refresh")
            self.assertEqual(mock_refresh.call_args.kwargs.get("modules"), ["cross-reactivity"])
            skipped = json.loads(result.skipped_manifest_path.read_text(encoding="utf-8"))
            self.assertEqual(
                {item["reason"] for item in skipped},
                {"duplicate_mouse_target:Q01279", "missing_human_report"},
            )
            summary = json.loads(result.summary_path.read_text(encoding="utf-8"))
            self.assertEqual(summary[0]["query"], "EGFR_MOUSE")
            self.assertEqual(summary[0]["mouse_ortholog_accession"], "Q01279")
            self.assertEqual(summary[0]["mouse_refseq_accession"], "NP_034255")
            self.assertEqual(summary[0]["mouse_sequence_source"], "uniprot_reviewed_mouse_gene_symbol")
            self.assertEqual(summary[0]["human_source_entry_name"], "EGFR_HUMAN")
            mouse_report = json.loads(Path(summary[0]["json_report"]).read_text(encoding="utf-8"))
            self.assertEqual(mouse_report["target"]["entry_name"], "EGFR_MOUSE")
            self.assertEqual(mouse_report["target"]["accession"], "Q01279")
            self.assertEqual(mouse_report["target"]["protein_name"], "Epidermal growth factor receptor")
            self.assertEqual(mouse_report["target"]["taxon_id"], 10090)
            self.assertEqual(mouse_report["target"]["sequence"], "MOUSESEQ")
            self.assertEqual(mouse_report["metadata"]["mouse_sequence_provenance"]["refseq_accession"], "NP_034255")
            self.assertEqual(mouse_report["metadata"]["human_source_target"]["entry_name"], "EGFR_HUMAN")
            self.assertEqual(mouse_report["metadata"]["source_database"], "mouse_openantigens_structure_analysis")
            self.assertEqual(mouse_report["family_context"]["identity_matrix_labels"], ["EGFR_MOUSE", "ERBB2_MOUSE"])
            self.assertIn("mouse_ortholog_matrix_caveat", mouse_report["family_context"]["metadata"])
            self.assertEqual(mouse_report["construct_details"][0]["sequence"], "MOUSESEQ")
            self.assertEqual(mouse_report["construct_details"][0]["structural_metrics"]["mean_plddt"], 91.2)
            self.assertEqual(mouse_report["metadata"]["construct_identity_label"], "Identity to mouse source")
            self.assertIsNotNone(result.portal_index_path)
            index_html = result.portal_index_path.read_text(encoding="utf-8")
            detail_html = (result.portal_index_path.parent / "reports" / "egfr_mouse.html").read_text(encoding="utf-8")
            index_data = (result.portal_index_path.parent / "portal-index-data.js").read_text(encoding="utf-8")
            downloads_tsv = (result.portal_index_path.parent / "downloads" / "agdesign2_portal_index.tsv").read_text(encoding="utf-8")
            download_manifest = json.loads((result.portal_index_path.parent / "downloads" / "download_manifest.json").read_text(encoding="utf-8"))
            self.assertIn("OpenAntigens Mouse", index_html)
            self.assertIn("openantigens mouse", index_html.lower())
            self.assertIn("OpenAntigens Human", index_html)
            self.assertIn('href="../../index.html"', detail_html)
            self.assertIn("--accent: #bd4d61", (result.portal_index_path.parent / "portal.css").read_text(encoding="utf-8"))
            self.assertNotIn("Filter by disease", index_html)
            self.assertNotIn("Top disease", index_html)
            self.assertNotIn("Open Targets Disease Associations", detail_html)
            self.assertIn("EGFR_MOUSE_1-8", detail_html)
            self.assertIn("structure-aware construct pipeline", detail_html)
            self.assertIn("The completed BLAST search returned no sequence similarity hits.", detail_html)
            self.assertNotIn("Showing 0 hit(s)", detail_html)
            self.assertNotIn("View in 3D", detail_html)
            self.assertNotIn("Advanced Membrane Engineering Suggestions", detail_html)
            self.assertNotIn("Experimental GPCR Engineering Variants", detail_html)
            self.assertNotIn("GPCR Annotation", detail_html)
            self.assertNotIn("HUMAN_BRIL", detail_html)
            self.assertNotIn("RefSeq NP_034255 because", detail_html)
            self.assertIn('"hasAlphafold":"0"', index_data)
            self.assertIn("mouse_ortholog_accession", downloads_tsv.splitlines()[0])
            self.assertIn("mouse_refseq_accession", downloads_tsv.splitlines()[0])
            self.assertIn("human_source_entry_name", downloads_tsv.splitlines()[0])
            self.assertIn("Q01279", downloads_tsv)
            self.assertIn("NP_034255", downloads_tsv)
            self.assertIn("EGFR_HUMAN", downloads_tsv)
            self.assertIn("extracellular_surface_identity_to_mouse_source", detail_html)
            self.assertNotIn("top_disease_name", downloads_tsv.splitlines()[0])
            self.assertNotIn("top_disease_name", download_manifest["files"][0]["columns"])
            help_html = (result.portal_index_path.parent / "help.html").read_text(encoding="utf-8")
            builder_html = (result.portal_index_path.parent / "builder.html").read_text(encoding="utf-8")
            constructs_html = (result.portal_index_path.parent / "constructs.html").read_text(encoding="utf-8")
            for support_html in (help_html, builder_html, constructs_html):
                self.assertNotIn("equivalent human and cynomolgus monkey regions", support_html)
            self.assertNotIn(
                "Check human, mouse, and cynomolgus monkey sequences",
                help_html,
            )
            self.assertNotIn(
                "Do mouse and cynomolgus monkey equivalent sequences",
                builder_html,
            )
            self.assertIn("originating human target", help_html.lower())
            self.assertIn("originating human target", constructs_html.lower())
            methods_html = (result.portal_index_path.parent / "methods.html").read_text(encoding="utf-8")
            self.assertIn("Reviewed UniProt homolog lookup", methods_html)
            self.assertNotIn(
                "Precomputed RefSeq/HCOP-derived ortholog table",
                methods_html,
            )

    def test_mouse_cross_reactivity_requires_blast_tools(self) -> None:
        with patch("agdesign2.mouse_portal.shutil.which", return_value=None):
            with self.assertRaisesRegex(RuntimeError, "Mouse portal cross-reactivity requires BLAST"):
                _require_blast_tools_for_mouse_cross_reactivity()

        with patch("agdesign2.mouse_portal.shutil.which", return_value="/usr/bin/tool"):
            _require_blast_tools_for_mouse_cross_reactivity()


def _human_report(entry_name: str = "EGFR_HUMAN", accession: str = "P00533", gene_symbol: str = "EGFR") -> dict:
    mouse_entry = entry_name.replace("_HUMAN", "_MOUSE")
    return {
        "target": {
            "accession": accession,
            "entry_name": entry_name,
            "gene_symbol": gene_symbol,
            "protein_name": "Epidermal growth factor receptor",
            "organism": "Homo sapiens",
            "taxon_id": 9606,
            "sequence": "HUMANFULLSEQ",
        },
        "ectodomain": {"start": 25, "end": 645, "label": "ecto", "source": "topology"},
        "topology": {"topology_class": "single_pass"},
        "construct_details": [
            {
                "name": "EGFR_HUMAN_25-645",
                "start": 25,
                "end": 645,
                "length": 621,
                "sequence": "HUMANSEQ",
                "score": 93.0,
                "rationale": "Complete ectodomain.",
                "classification": "multi_domain_unit",
                "warnings": [],
                "evidence": ["Topology-derived ectodomain"],
                "ptms": [{"type": "CHAIN", "start": 1, "end": 1338, "description": "human chain"}],
                "cysteine_analysis": [{"position": 1338, "paired_with": None}],
                "gpcr_segments": [{"name": "H8", "start": 1320, "end": 1338}],
                "gpcr_generic_range": "8.47x47 to 8.63x63",
                "structural_metrics": {"mean_plddt": 91.0},
                "human_entry_name": "EGFR_HUMAN",
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
                    }
                ],
            }
        ],
        "construct_recommendations": [],
        "ectodomain_homology": [
            {
                "species": "mouse",
                "accession": "Q01279",
                "entry_name": mouse_entry,
                "ectodomain_start": 25,
                "ectodomain_end": 647,
                "identity": 88.73,
                "coverage": 100.0,
                "available": True,
            }
        ],
        "cross_reactivity_hits": [],
        "full_length_cross_reactivity_hits": [],
        "cross_reactivity_hits": [
            {
                "subject_id": "sp|Q01279|EGFR_MOUSE",
                "description": "sp|Q01279|EGFR_MOUSE Epidermal growth factor receptor OS=Mus musculus OX=10090 GN=Egfr PE=1 SV=2",
            }
        ],
        "furin_sites": [{"start": 1337, "end": 1338, "motif": "RR"}],
        "cysteine_analysis": [{"position": 1338, "warning": "human-source coordinate"}],
        "features": [{"type": "CHAIN", "start": 1, "end": 1338, "description": "human chain"}],
        "ptms": [{"type": "CHAIN", "start": 1, "end": 1338, "description": "human chain"}],
        "residue_annotations": [{"position": 1338, "region": "human tail"}],
        "advanced_membrane_suggestions": [{"title": "Human C-tail trim", "start": 1330, "end": 1338}],
        "gpcr_engineering_variants": [{"name": "EGFR_HUMAN_BRIL_ICL3_fusion", "sequence": "HUMANBRIL"}],
        "gpcr_annotation": {"entry_name": "egfr_human", "segments": [{"name": "H8", "start": 1320, "end": 1338}]},
        "assembly_requirements": [],
        "experimental_constructs": [],
        "interpro_annotations": [],
        "notes": [],
        "family_context": {
            "gene_symbol": gene_symbol,
            "source": "Precomputed InterPro family alignments",
            "family_names": ["ERBB family"],
            "members": [
                {"gene_symbol": "EGFR", "accession": "P00533", "entry_name": "EGFR_HUMAN", "ectodomain_start": 25, "ectodomain_end": 645},
                {"gene_symbol": "ERBB2", "accession": "P04626", "entry_name": "ERBB2_HUMAN", "ectodomain_start": 25, "ectodomain_end": 645},
            ],
            "identity_matrix_labels": ["EGFR_HUMAN", "ERBB2_HUMAN"],
            "identity_matrix": [[100.0, 50.0], [60.0, 100.0]],
            "metadata": {},
        },
    }


def _mouse_analyzed_report(entry_name: str = "EGFR_MOUSE", accession: str = "Q01279", gene_symbol: str = "Egfr") -> dict:
    report = _human_report()
    sequence = "MOUSESEQ" if entry_name == "EGFR_MOUSE" else "MOUSESEQUENCE"
    report["target"] = {
        "accession": accession,
        "entry_name": entry_name,
        "gene_symbol": gene_symbol,
        "protein_name": "Epidermal growth factor receptor",
        "organism": "Mus musculus",
        "taxon_id": 10090,
        "sequence": sequence,
    }
    report["ectodomain"] = {"start": 1, "end": len(sequence), "label": "ecto", "source": "topology"}
    report["construct_details"] = [
        {
            "name": f"{entry_name}_1-{len(sequence)}",
            "start": 1,
            "end": len(sequence),
            "length": len(sequence),
            "sequence": sequence,
            "score": 94.0,
            "rationale": "Mouse structure-derived ectodomain.",
            "classification": "multi_domain_unit",
            "warnings": [],
            "evidence": ["AlphaFold-guided structured region"],
            "ptms": [],
            "cysteine_analysis": [],
            "gpcr_segments": [],
            "gpcr_generic_range": None,
            "structural_metrics": {"mean_plddt": 91.2},
            "human_entry_name": None,
            "homologs": [],
        }
    ]
    report["construct_recommendations"] = []
    report["ectodomain_homology"] = []
    report["cross_reactivity_hits"] = []
    report["full_length_cross_reactivity_hits"] = []
    report["features"] = []
    report["ptms"] = []
    report["residue_annotations"] = []
    report["advanced_membrane_suggestions"] = []
    report["gpcr_engineering_variants"] = []
    report["gpcr_annotation"] = None
    report["notes"] = []
    return report
