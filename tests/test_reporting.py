from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from agdesign2.models import (
    AnalysisNote,
    AnalysisReport,
    AssemblyRequirement,
    BlastHit,
    CanonicalFamily,
    ComplexPortalComplex,
    ComplexPortalParticipant,
    CysteineFinding,
    ConstructDetail,
    ConstructSplitDiagnostic,
    ConstructSuggestion,
    ConstructStructuralMetrics,
    DomainAnnotation,
    ExperimentalConstruct,
    FamilyContext,
    FamilyMember,
    Feature,
    GPCRdbAnnotation,
    GPCRdbMotif,
    GPCRdbSegment,
    HomologyRecord,
    MembraneEngineeringSuggestion,
    Region,
    TargetRecord,
)
from agdesign2.reporting import _construct_name_md


class ReportingTests(unittest.TestCase):
    def test_pdb_construct_name_links_to_rcsb(self) -> None:
        self.assertEqual(
            _construct_name_md("pdb_1ira_x", "1IRA"),
            "[`pdb_1ira_x`](https://www.rcsb.org/structure/1IRA)",
        )
        # Non-PDB constructs stay plain code spans.
        self.assertEqual(_construct_name_md("full_ectodomain", None), "`full_ectodomain`")

    def test_renders_markdown_summary(self) -> None:
        report = AnalysisReport(
            target=TargetRecord(
                accession="P00533",
                entry_name="EGFR_HUMAN",
                gene_symbol="EGFR",
                protein_name="Epidermal growth factor receptor",
                organism="Homo sapiens",
                taxon_id=9606,
                sequence="M" * 100,
            ),
            resolution={"query": "EGFR"},
            features=[],
            ectodomain=Region(start=25, end=80, label="N-terminal ectodomain", source="topology", confidence=0.8),
            structural_regions=[
                Region(
                    start=25,
                    end=80,
                    label="Lenient structured region",
                    source="AlphaFold",
                    confidence=0.77,
                    metadata={
                        "stringency": "lenient",
                        "mean_intra_pae": 6.2,
                        "included_interpro_domains": ["EGF receptor L domain (IPR000719)"],
                    },
                ),
                Region(
                    start=30,
                    end=60,
                    label="Predicted domain",
                    source="AlphaFold",
                    confidence=0.82,
                    metadata={
                        "mean_intra_pae": 4.5,
                        "split_history": [
                            {
                                "parent_start": 25,
                                "parent_end": 80,
                                "boundary_after": 60,
                                "score": 11.2,
                                "inter_block_pae": 18.0,
                                "intra_block_pae": 6.3,
                                "separation": 11.7,
                                "boundary_pae": 19.4,
                                "linker_plddt": 63.5,
                            }
                        ],
                    },
                )
            ],
            experimental_constructs=[
                ExperimentalConstruct(
                    pdb_id="1ABC",
                    method="X-ray",
                    resolution="2.0 A",
                    chains=[Region(start=25, end=80, label="Chain A", source="PDB")],
                ),
                ExperimentalConstruct(
                    pdb_id="2DEF",
                    method="X-ray",
                    resolution="2.1 A",
                    chains=[Region(start=25, end=80, label="Chain B", source="PDB")],
                ),
            ],
            construct_recommendations=[
                ConstructSuggestion(
                    name="full_ectodomain",
                    start=25,
                    end=80,
                    score=90.0,
                    rationale="Complete extracellular domain.",
                    classification="multi_domain_unit",
                ),
                ConstructSuggestion(
                    name="membrane_expression_trimmed",
                    start=1,
                    end=95,
                    score=88.0,
                    rationale="Trimmed membrane construct.",
                    classification="membrane_expression",
                ),
            ],
            furin_sites=[],
            cysteine_analysis=[
                CysteineFinding(position=31, paired_with=58, surface_exposed=False),
                CysteineFinding(position=58, paired_with=31, surface_exposed=False),
                CysteineFinding(
                    position=75,
                    paired_with=None,
                    surface_exposed=True,
                    warning="Unpaired cysteine appears surface exposed.",
                ),
            ],
            ligand_interactions=[
                Feature(
                    type="REGION",
                    start=30,
                    end=55,
                    description="Ligand-binding region for EGF-like ligands",
                    metadata={"interaction_summary": "Ligand-binding region for EGF-like ligands"},
                )
            ],
            assembly_requirements=[
                AssemblyRequirement(
                    summary="Integrin alpha chains require an integrin beta partner for the biologically relevant cell-surface heterodimer.",
                    obligatory=True,
                    confidence="high",
                    source="Family rule + UniProt",
                    classification="obligatory_partner_requirement",
                    partners=["ITGB1"],
                    rationale=["Target gene symbol ITGA5 matches the integrin alpha-chain family."],
                ),
                AssemblyRequirement(
                    summary="Curated UniProt subunit annotation lists interaction partners or complexes, but this is not sufficient evidence for an obligatory partner requirement.",
                    obligatory=False,
                    confidence="medium",
                    source="UniProt SUBUNIT",
                    classification="interaction_context",
                    partners=["CIB1"],
                    rationale=["Interacts with CIB1."],
                )
            ],
            interpro_annotations=[
                DomainAnnotation(
                    accession="IPR000719",
                    name="EGF receptor L domain",
                    type="domain",
                    source_database="INTERPRO",
                    start=30,
                    end=55,
                    representative=True,
                ),
                DomainAnnotation(
                    accession="PF01030",
                    name="Recep_L_domain",
                    type="domain",
                    source_database="PFAM",
                    start=30,
                    end=55,
                    representative=True,
                    integrated_accession="IPR000719",
                    integrated_name="EGF receptor L domain",
                ),
            ],
            species_name_matches=[],
            ectodomain_homology=[
                HomologyRecord(
                    species="mouse",
                    accession="Q01279",
                    entry_name="EGFR_MOUSE",
                    ectodomain_start=25,
                    ectodomain_end=82,
                    identity=88.73,
                    coverage=100.0,
                )
            ],
            family_context=FamilyContext(
                gene_symbol="EGFR",
                source="HGNC",
                family_names=["ERBB family"],
                members=[
                    FamilyMember(
                        gene_symbol="EGFR",
                        gene_name="Epidermal growth factor receptor",
                        accession="P00533",
                        entry_name="EGFR_HUMAN",
                        ectodomain_start=25,
                        ectodomain_end=80,
                        ectodomain_length=56,
                    ),
                    FamilyMember(
                        gene_symbol="ERBB2",
                        gene_name="erb-b2 receptor tyrosine kinase 2",
                        accession="P04626",
                        entry_name="ERBB2_HUMAN",
                        ectodomain_start=23,
                        ectodomain_end=652,
                        ectodomain_length=630,
                    ),
                ],
                identity_matrix_labels=["EGFR", "ERBB2"],
                identity_matrix=[[100.0, 43.9], [43.9, 100.0]],
            ),
            cross_reactivity_hits=[
                BlastHit(
                    subject_id="sp|Q01279|EGFR_MOUSE",
                    description="EGFR mouse",
                    species="Mus musculus",
                    identity=88.7,
                    coverage=100.0,
                    alignment_length=56,
                    evalue=1e-30,
                    bitscore=250.0,
                    query_start=1,
                    query_end=56,
                    subject_start=1,
                    subject_end=56,
                )
            ],
            canonical_family=CanonicalFamily(
                accession="IPR050122",
                name="Receptor Tyrosine Kinase",
                source_database="INTERPRO",
                fragment_count=1,
                start=58,
                end=974,
                covered_length=917,
                coverage_fraction=76.16,
            ),
            complex_portal_complexes=[
                ComplexPortalComplex(
                    complex_ac="CPX-INT1",
                    name="Integrin alpha-5/beta-1 complex",
                    species="Homo sapiens; 9606",
                    predicted_complex=False,
                    evidence_code="ECO:0000353",
                    evidence_description="physical interaction evidence",
                    confidence_score=5,
                    complex_assemblies=["Heterodimer"],
                    participants=[
                        ComplexPortalParticipant(identifier="P08648", name="ITGA5", interactor_type="protein"),
                        ComplexPortalParticipant(identifier="P05556", name="ITGB1", interactor_type="protein"),
                    ],
                )
            ],
            advanced_membrane_suggestions=[
                MembraneEngineeringSuggestion(
                    category="gpcr_loop_engineering",
                    title="ICL-focused loop engineering screen",
                    summary="Consider an internal cytoplasmic loop engineering series.",
                    rationale=["Longest internal cytoplasmic loop is suitable for engineering."],
                    suggested_actions=["Test BRIL or T4L insertion in the loop."],
                    evidence=["Topology-derived internal cytoplasmic loop."],
                    warnings=["Advisory only."],
                    start=250,
                    end=310,
                )
            ],
            gpcr_annotation=GPCRdbAnnotation(
                entry_name="egfr_human_gpcr_test",
                accession="P00533",
                name="Test GPCR",
                family_path=["Class A (Rhodopsin)", "Test receptors"],
                residue_numbering_scheme="GPCRdb(A)",
                url="https://gpcrdb.org/protein/egfr_human_gpcr_test",
                segments=[GPCRdbSegment(name="TM1", start=250, end=280, generic_start="1.32x32", generic_end="1.62x62")],
                conserved_motifs=[
                    GPCRdbMotif(
                        name="DRY / E/DRY activation motif",
                        positions=[300, 301, 302],
                        generic_numbers=["3.49x49", "3.50x50", "3.51x51"],
                        sequence="DRY",
                    )
                ],
            ),
            notes=[AnalysisNote(severity="info", message="Resolved from official gene symbol.", source="UniProt")],
        )

        markdown = report.to_markdown()
        self.assertIn("# EGFR_HUMAN Antigen Design Report", markdown)
        self.assertIn("## Assembly / Partner Requirements", markdown)
        self.assertIn("## Advanced Membrane Engineering Suggestions", markdown)
        self.assertIn("ICL-focused loop engineering screen", markdown)
        self.assertIn("### Soluble ECD Constructs", markdown)
        self.assertIn("### Native Membrane-Expression Constructs", markdown)
        self.assertIn("| Classification | Obligatory | Confidence | Summary | Named Partner(s) | Source | Rationale |", markdown)
        self.assertIn("obligatory_partner_requirement", markdown)
        self.assertIn("ITGB1", markdown)
        self.assertNotIn("interaction_context", markdown)
        self.assertIn("Additional interaction-context annotations omitted from summary: 1", markdown)
        self.assertIn("## Complex Portal Context", markdown)
        self.assertIn("Curated Complex Portal records:", markdown)
        self.assertIn("## GPCR Annotation", markdown)
        self.assertIn("egfr_human_gpcr_test", markdown)
        self.assertIn("DRY / E/DRY activation motif", markdown)
        self.assertIn("`CPX-INT1` Integrin alpha-5/beta-1 complex", markdown)
        self.assertIn("## Construct Recommendations", markdown)
        self.assertIn("| Construct | Start | End | Length | Score | Classification | Rationale |", markdown)
        self.assertNotIn("| Construct | Start | End | Length | Score | Classification | Rationale | Evidence | Warnings |", markdown)
        self.assertIn("88.73%", markdown)
        self.assertIn("Cross-Reactivity Hits", markdown)
        self.assertIn("| Rank | Hit | Gene | Species | Bit score | Identity | Coverage | E-value |", markdown)
        self.assertIn("EGFR_MOUSE", markdown)
        self.assertIn("## InterPro / Pfam Domain Architecture", markdown)
        self.assertIn("## Canonical Family", markdown)
        self.assertIn("IPR050122", markdown)
        self.assertIn("| Start | End | Length | Source | Accession | Type | Name | Representative | Integrated |", markdown)
        self.assertIn("IPR000719", markdown)
        self.assertIn("Ligand Interaction Regions", markdown)
        self.assertIn("| Start | End | Type | Summary |", markdown)
        self.assertIn("Ligand-binding region for EGF-like ligands", markdown)
        self.assertIn("Family ectodomain directional identity matrix (%)", markdown)
        self.assertIn("Rows are normalized to the row protein ectodomain length", markdown)
        self.assertIn("| Protein | EGFR | ERBB2 |", markdown)
        self.assertNotIn("## Family Members", markdown)
        self.assertLess(markdown.index("## Family Context"), markdown.index("## Assembly / Partner Requirements"))
        self.assertIn("| Region | Stringency | Start | End | Length | Source | Confidence | Mean Intra-domain PAE | Included InterPro Domains |", markdown)
        self.assertIn("4.50", markdown)
        self.assertIn("6.20", markdown)
        self.assertIn("EGF receptor L domain (IPR000719)", markdown)
        self.assertIn("| PDB | Chains | Method | Resolution |", markdown)
        self.assertIn("1ABC", markdown)
        self.assertNotIn("2DEF", markdown)
        self.assertIn("PAE split diagnostics:", markdown)
        self.assertIn("| Parent Region | Boundary After | Score | Inter-block PAE | Intra-block PAE | Separation | Local Boundary PAE | Linker pLDDT |", markdown)
        self.assertIn("25-80", markdown)
        self.assertIn("63.50", markdown)
        self.assertIn("| Residue(s) | Status | Exposure | Closest unpaired cysteine | Warning |", markdown)
        self.assertIn("C31-C58", markdown)
        self.assertEqual(markdown.count("C31-C58"), 1)
        self.assertIn("C75", markdown)

    def test_renders_construct_metrics(self) -> None:
        report = AnalysisReport(
            target=TargetRecord(
                accession="P00533",
                entry_name="EGFR_HUMAN",
                gene_symbol="EGFR",
                protein_name="Epidermal growth factor receptor",
                organism="Homo sapiens",
                taxon_id=9606,
                sequence="M" * 100,
            ),
            resolution={"query": "EGFR"},
            features=[],
            ectodomain=Region(start=25, end=80, label="N-terminal ectodomain", source="topology", confidence=0.8),
            structural_regions=[],
            experimental_constructs=[],
            construct_recommendations=[
                ConstructSuggestion(
                    name="full_ectodomain",
                    start=25,
                    end=80,
                    score=90.0,
                    rationale="Complete extracellular domain.",
                    classification="complete_domain",
                )
            ],
            furin_sites=[],
            cysteine_analysis=[],
            assembly_requirements=[],
            interpro_annotations=[],
            species_name_matches=[],
            ectodomain_homology=[],
            family_context=None,
            cross_reactivity_hits=[],
            construct_details=[
                ConstructDetail(
                    name="full_ectodomain",
                    start=25,
                    end=80,
                    length=56,
                    sequence="A" * 56,
                    score=90.0,
                    rationale="Complete extracellular domain.",
                    classification="complete_domain",
                    human_entry_name="EGFR_HUMAN",
                    domain_summary="Captures one complete annotated domain: Domain A (IPR000001)",
                    complete_domain_annotations=[
                        DomainAnnotation(
                            accession="IPR000001",
                            name="Domain A",
                            type="domain",
                            source_database="INTERPRO",
                            start=25,
                            end=80,
                            representative=True,
                        )
                    ],
                    structural_metrics=ConstructStructuralMetrics(
                        mean_plddt=88.5,
                        min_plddt=62.0,
                        max_plddt=97.0,
                        structured_fraction=0.91,
                        mean_intra_pae=3.4,
                        max_intra_pae=12.8,
                    ),
                    split_diagnostics=[
                        ConstructSplitDiagnostic(
                            parent_start=25,
                            parent_end=80,
                            boundary_after=60,
                            score=11.2,
                            inter_block_pae=18.0,
                            intra_block_pae=6.3,
                            separation=11.7,
                            boundary_pae=19.4,
                            linker_plddt=63.5,
                        )
                    ],
                )
            ],
        )
        markdown = report.to_markdown()
        self.assertIn("Classification: complete_domain", markdown)
        self.assertIn("Architecture summary: Captures one complete annotated domain: Domain A (IPR000001)", markdown)
        self.assertIn("Complete domain hits: INTERPRO IPR000001 Domain A 25-80", markdown)
        self.assertIn("AlphaFold quality: mean pLDDT 88.50; min 62.00; max 97.00", markdown)
        self.assertIn("Structured coverage: 91.00% residues at or above pLDDT 70", markdown)
        self.assertIn("Intra-construct PAE: mean 3.40; max 12.80", markdown)
        self.assertIn("| Parent Region | Boundary After | Score | Inter-block PAE | Intra-block PAE | Separation | Local Boundary PAE | Linker pLDDT |", markdown)
        self.assertIn("| 25-80 | 60 | 11.20 | 18.00 | 6.30 | 11.70 | 19.40 | 63.50 |", markdown)


if __name__ == "__main__":
    unittest.main()
