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
from agdesign2.models import (
    AnalysisNote,
    BlastHit,
    ComplexPortalComplex,
    ComplexPortalParticipant,
    ConstructSuggestion,
    DomainAnnotation,
    ExperimentalConstruct,
    Feature,
    GPCRdbAnnotation,
    GPCRdbMotif,
    GPCRdbResidue,
    GPCRdbSegment,
    HomologyRecord,
    Region,
    TargetRecord,
)
from agdesign2.pipeline import AntigenAnalyzer
from agdesign2.refseq import RefSeqProtein
from agdesign2.structure_utils import parse_alphafold_pdb
from agdesign2.topology import derive_ectodomain
from agdesign2.uniprot import TargetResolution


def build_entry(entry_name: str, accession: str, sequence: str, organism: str, taxon_id: int, gene: str) -> dict:
    return {
        "primaryAccession": accession,
        "uniProtkbId": entry_name,
        "sequence": {"value": sequence},
        "organism": {"scientificName": organism, "taxonId": taxon_id},
        "genes": [{"geneName": {"value": gene}}],
        "proteinDescription": {"recommendedName": {"fullName": {"value": gene}}},
        "features": [],
        "uniProtKBCrossReferences": [],
    }


class FakeUniProtClient:
    def __init__(self, target_entry: dict, mouse_entry: dict | None, macfa_entry: dict | None, extra_entries: list[dict] | None = None) -> None:
        self.target_entry = target_entry
        self.mouse_entry = mouse_entry
        self.macfa_entry = macfa_entry
        self.extra_entries = extra_entries or []

    def resolve_target(self, query: str, target_species_taxon: int = 9606) -> TargetResolution:
        from agdesign2.models import AnalysisNote, TargetRecord

        gene_symbol = self.target_entry["genes"][0]["geneName"]["value"]

        return TargetResolution(
            target=TargetRecord(
                accession=self.target_entry["primaryAccession"],
                entry_name=self.target_entry["uniProtkbId"],
                gene_symbol=gene_symbol,
                protein_name="Test receptor",
                organism="Homo sapiens",
                taxon_id=9606,
                sequence=self.target_entry["sequence"]["value"],
            ),
            entry=self.target_entry,
            notes=[AnalysisNote(severity="info", message="resolved")],
            query_type="gene_symbol",
        )

    def get_features(self, entry: dict) -> list:
        from agdesign2.models import Feature

        return [
            Feature(type="SIGNAL", start=1, end=2, description="Signal peptide"),
            Feature(type="TOPO_DOM", start=3, end=8, description="Extracellular"),
            Feature(type="REGION", start=3, end=6, description="Ligand-binding region"),
            Feature(type="DOMAIN", start=3, end=5, description="Domain A"),
            Feature(type="DOMAIN", start=6, end=8, description="Domain B"),
            Feature(type="CARBOHYD", start=4, end=4, description="N-linked glycosylation"),
            Feature(type="SITE", start=6, end=6, description="Cleavage site"),
            Feature(type="TRANSMEM", start=9, end=10, description="Helical"),
            Feature(type="TOPO_DOM", start=11, end=12, description="Cytoplasmic"),
        ]

    def get_experimental_constructs(self, entry: dict) -> list:
        from agdesign2.models import ExperimentalConstruct, Region

        return [
            ExperimentalConstruct(
                pdb_id="1ABC",
                method="X-ray",
                resolution="2.0 A",
                chains=[Region(start=3, end=8, label="Chain A", source="PDB", metadata={"chain": "A"})],
            ),
            ExperimentalConstruct(
                pdb_id="9ZZZ",
                method="X-ray",
                resolution="1.9 A",
                chains=[Region(start=11, end=12, label="Chain B", source="PDB", metadata={"chain": "B"})],
            )
        ]

    def get_comment_texts(self, entry: dict, comment_type: str) -> list[str]:
        texts: list[str] = []
        for comment in entry.get("comments", []):
            if str(comment.get("commentType", "")).upper() != comment_type.upper():
                continue
            for text_block in comment.get("texts", []):
                value = text_block.get("value")
                if isinstance(value, str):
                    texts.append(value)
        return texts

    def fetch_same_name_species_match(self, base_name: str, species_suffix: str) -> dict | None:
        if species_suffix == "MOUSE":
            return self.mouse_entry
        if species_suffix == "MACFA":
            return self.macfa_entry
        return None

    def _fetch_entry(self, accession: str) -> dict:
        for entry in (self.target_entry, self.mouse_entry, self.macfa_entry, *self.extra_entries):
            if entry and entry["primaryAccession"] == accession:
                return entry
        raise KeyError(accession)

    def fetch_entry_by_name(self, entry_name: str) -> dict | None:
        for entry in (self.target_entry, self.mouse_entry, self.macfa_entry, *self.extra_entries):
            if entry and entry["uniProtkbId"] == entry_name:
                return entry
        return None

    def search(self, query: str, *, size: int = 10) -> dict:
        results = []
        if "organism_id:9541" in query and "gene_exact:TEST" in query:
            for entry in (self.macfa_entry, *self.extra_entries):
                if entry is not None:
                    results.append({"primaryAccession": entry["primaryAccession"]})
                    break
        return {"results": results[:size]}

    def _target_from_entry(self, entry: dict):
        from agdesign2.models import TargetRecord

        gene_symbol = entry["genes"][0]["geneName"]["value"]

        return TargetRecord(
            accession=entry["primaryAccession"],
            entry_name=entry["uniProtkbId"],
            gene_symbol=gene_symbol,
            protein_name="Test receptor",
            organism=entry["organism"]["scientificName"],
            taxon_id=entry["organism"]["taxonId"],
            sequence=entry["sequence"]["value"],
        )


class FakeAlphaFoldClient:
    def __init__(self, tmpdir: Path) -> None:
        self.tmpdir = tmpdir

    def ensure_artifacts(
        self,
        accession: str,
        canonical_sequence: str | None = None,
        canonical_isoform_id: str | None = None,
    ) -> tuple[Path, Path]:
        pdb_path = self.tmpdir / f"{accession}.pdb"
        pae_path = self.tmpdir / f"{accession}.json"
        pdb_path.write_text(
            "\n".join(
                [
                    "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 95.00           C",
                    "ATOM      2  CA  ALA A   2       1.000   0.000   0.000  1.00 95.00           C",
                    "ATOM      3  CA  CYS A   3       2.000   0.000   0.000  1.00 95.00           C",
                    "ATOM      4  SG  CYS A   3       2.300   0.000   0.000  1.00 95.00           S",
                    "ATOM      5  CA  ALA A   4       3.000   0.000   0.000  1.00 95.00           C",
                    "ATOM      6  CA  CYS A   5       4.000   0.000   0.000  1.00 95.00           C",
                    "ATOM      7  SG  CYS A   5       4.600   0.000   0.000  1.00 95.00           S",
                    "ATOM      8  CA  ALA A   6       5.000   0.000   0.000  1.00 95.00           C",
                    "ATOM      9  CA  ALA A   7       6.000   0.000   0.000  1.00 95.00           C",
                    "ATOM     10  CA  ALA A   8       7.000   0.000   0.000  1.00 95.00           C",
                    "ATOM     11  CA  ALA A   9       8.000   0.000   0.000  1.00 40.00           C",
                    "ATOM     12  CA  ALA A  10       9.000   0.000   0.000  1.00 40.00           C",
                    "END",
                ]
            ),
            encoding="utf-8",
        )
        pae = {"predicted_aligned_error": [[1.0] * len(canonical_sequence or "M" * 10) for _ in range(len(canonical_sequence or "M" * 10))]}
        pae_path.write_text(json.dumps(pae), encoding="utf-8")
        return pdb_path, pae_path


class FakeHGNCClient:
    def fetch_family_context(self, gene_symbol: str | None):
        from agdesign2.models import FamilyContext

        return FamilyContext(gene_symbol=gene_symbol or "TEST", source="HGNC", family_names=["Test family"])


class RaisingHGNCClient:
    def fetch_family_context(self, gene_symbol: str | None):
        raise AssertionError("HGNC fallback should not be used when precomputed family data is available")


class ConstructPtmOverlapTests(unittest.TestCase):
    def test_construct_ptms_treat_disulfides_as_endpoint_pairs(self) -> None:
        analyzer = AntigenAnalyzer.__new__(AntigenAnalyzer)
        analyzer.config = AnalysisConfig()
        constructs = [
            ConstructSuggestion(
                name="inside_bond_span",
                start=298,
                end=328,
                score=1.0,
                rationale="test",
            ),
            ConstructSuggestion(
                name="overlaps_bond_endpoint",
                start=298,
                end=330,
                score=1.0,
                rationale="test",
            ),
        ]
        ptms = [
            Feature(type="DISULFID", start=264, end=329, description=""),
            Feature(type="CARBOHYD", start=306, end=306, description="N-linked glycosylation"),
        ]

        details = analyzer._build_construct_details(
            target_entry_name="TEST_HUMAN",
            target_sequence="A" * 400,
            constructs=constructs,
            homolog_context={},
            ectodomain=Region(start=1, end=400, label="test", source="test"),
            ligand_interactions=[],
            interpro_annotations=[],
            ptms=ptms,
        )

        self.assertEqual([feature.type for feature in details[0].ptms], ["CARBOHYD"])
        self.assertEqual([feature.type for feature in details[1].ptms], ["DISULFID", "CARBOHYD"])


class FakeRefSeqClient:
    def __init__(self, sequence: str | None) -> None:
        self.sequence = sequence
        self.calls: list[tuple[str, str]] = []

    def fetch_canonical_protein(self, *, gene_symbol: str, organism: str):
        self.calls.append((gene_symbol, organism))
        if organism != "Macaca fascicularis" or self.sequence is None:
            return None
        return RefSeqProtein(
            target=TargetRecord(
                accession="XP_TEST_1",
                entry_name="XP_TEST_1",
                gene_symbol=gene_symbol,
                protein_name=f"{gene_symbol} protein",
                organism=organism,
                taxon_id=9541,
                sequence=self.sequence,
            ),
            notes=[AnalysisNote(severity="info", message="Resolved from RefSeq canonical protein candidate XP_TEST_1.")],
            title=f"{gene_symbol} isoform X1",
        )


class FakeBlastClient:
    def search(self, query_sequence: str, *, target_accession: str, target_entry_name: str) -> list[BlastHit]:
        return [
            BlastHit(
                subject_id="sp|Q01279|TEST_MOUSE",
                description="mock hit",
                species="Mus musculus",
                identity=87.5,
                coverage=95.0,
                alignment_length=len(query_sequence),
                evalue=1e-30,
                bitscore=220.0,
                query_start=1,
                query_end=len(query_sequence),
                subject_start=5,
                subject_end=5 + len(query_sequence) - 1,
            )
        ]

    def search_species(
        self,
        query_sequence: str,
        *,
        species: str,
        target_accession: str = "",
        target_entry_name: str = "",
        ortholog_search: bool = False,
    ) -> list[BlastHit]:
        if species == "macaca_fascicularis":
            return [
                BlastHit(
                    subject_id="sp|A00003|ALTX_MACFA",
                    description="mock macaque hit",
                    species="Macaca fascicularis",
                    identity=92.0,
                    coverage=100.0,
                    alignment_length=len(query_sequence),
                    evalue=1e-50,
                    bitscore=300.0,
                    query_start=1,
                    query_end=len(query_sequence),
                    subject_start=1,
                    subject_end=len(query_sequence),
                )
            ]
        return []

    def ensure_ortholog_databases(self):
        return {"mouse": Path("/tmp/mouse"), "macaca_fascicularis": Path("/tmp/macfa")}


class UnresolvedTopologyUniProtClient(FakeUniProtClient):
    def get_features(self, entry: dict) -> list:
        return [
            Feature(type="TRANSMEM", start=2, end=3, description="Helical; Signal-anchor for type III membrane protein"),
            Feature(type="TOPO_DOM", start=4, end=len(entry["sequence"]["value"]), description="Cytoplasmic"),
        ]


class FakeGPCRdbClient:
    def __init__(self, annotation: GPCRdbAnnotation | None) -> None:
        self.annotation = annotation
        self.calls: list[tuple[str, str]] = []

    def fetch_annotation(self, *, entry_name: str, accession: str, sequence: str):
        self.calls.append((entry_name, accession))
        return self.annotation


class FakeInterProClient:
    def fetch_annotations(self, accession: str) -> list[DomainAnnotation]:
        return [
            DomainAnnotation(
                accession="IPR_FAMILY",
                name="Broad family",
                type="family",
                source_database="INTERPRO",
                start=3,
                end=8,
                representative=True,
            ),
            DomainAnnotation(
                accession="IPR_SUPER_A",
                name="Domain A superfamily",
                type="homologous_superfamily",
                source_database="INTERPRO",
                start=3,
                end=5,
                representative=True,
            ),
            DomainAnnotation(
                accession="IPR000001",
                name="Domain A",
                type="domain",
                source_database="INTERPRO",
                start=3,
                end=5,
                representative=True,
            ),
            DomainAnnotation(
                accession="IPR000002",
                name="Domain B",
                type="domain",
                source_database="INTERPRO",
                start=6,
                end=8,
                representative=True,
            ),
            DomainAnnotation(
                accession="PF00001",
                name="Domain A Pfam",
                type="domain",
                source_database="PFAM",
                start=3,
                end=5,
                representative=True,
                integrated_accession="IPR000001",
                integrated_name="Domain A",
            ),
            DomainAnnotation(
                accession="PF_REPEAT_DUP",
                name="Domain-like repeat duplicate",
                type="repeat",
                source_database="PFAM",
                start=3,
                end=5,
                representative=True,
                integrated_accession="IPR000001",
                integrated_name="Domain A",
            ),
            DomainAnnotation(
                accession="IPR999999",
                name="Outside domain",
                type="domain",
                source_database="INTERPRO",
                start=2,
                end=9,
                representative=True,
            ),
        ]


class FakeComplexPortalClient:
    def __init__(self) -> None:
        self.calls: list[dict[str, object]] = []

    def fetch_complexes_for_target(
        self,
        *,
        accession: str,
        gene_symbol: str | None,
        taxon_id: int | None,
        entry_name: str | None = None,
        max_pages: int = 2,
    ) -> list[ComplexPortalComplex]:
        self.calls.append(
            {
                "accession": accession,
                "gene_symbol": gene_symbol,
                "taxon_id": taxon_id,
                "entry_name": entry_name,
                "max_pages": max_pages,
            }
        )
        if gene_symbol == "ITGA5":
            if taxon_id == 10090:
                return [
                    ComplexPortalComplex(
                        complex_ac="CPX-MOUSE-INT1",
                        name="Mouse integrin alpha-5/beta-1 complex",
                        species="Mus musculus; 10090",
                        predicted_complex=False,
                        evidence_code="ECO:0000353",
                        evidence_description="physical interaction evidence",
                        confidence_score=4,
                        complex_assemblies=["Heterodimer"],
                        properties=["Stable integrin alpha/beta heterodimer."],
                        participants=[
                            ComplexPortalParticipant(identifier="Q9Z0N0", name="ITGA5", interactor_type="protein", stoichiometry="minValue: 1, maxValue: 1"),
                            ComplexPortalParticipant(identifier="P09055", name="ITGB1", interactor_type="protein", stoichiometry="minValue: 1, maxValue: 1"),
                        ],
                    )
                ]
            return [
                ComplexPortalComplex(
                    complex_ac="CPX-INT1",
                    name="Integrin alpha-5/beta-1 complex",
                    species="Homo sapiens; 9606",
                    predicted_complex=False,
                    evidence_code="ECO:0000353",
                    evidence_description="physical interaction evidence",
                    confidence_score=5,
                    complex_assemblies=["Heterodimer"],
                    properties=["Stable integrin alpha/beta heterodimer."],
                    participants=[
                        ComplexPortalParticipant(identifier="P08648", name="ITGA5", interactor_type="protein", stoichiometry="minValue: 1, maxValue: 1"),
                        ComplexPortalParticipant(identifier="P05556", name="ITGB1", interactor_type="protein", stoichiometry="minValue: 1, maxValue: 1", binding_regions=["111-222"]),
                    ],
                )
            ]
        return []


class FakeAssetGenerator:
    def generate_assets(
        self,
        *,
        construct,
        pdb_path,
        residue_plddt,
        pae_matrix,
        ectodomain,
        output_dir,
        render_structure_image=True,
        render_quality_plot=True,
    ):
        return "structure.png", "quality.png"


class PipelineTests(unittest.TestCase):
    def test_prefers_precomputed_ortholog_and_family_records(self) -> None:
        target = build_entry("TEST_HUMAN", "P00001", "MMAAACCCRRRR", "Homo sapiens", 9606, "TEST")
        config = AnalysisConfig(
            min_structured_segment=3,
            min_domain_size=2,
            min_construct_length=2,
            blast_species=("human",),
            generate_assets=False,
            prefer_precomputed_references=True,
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            ortholog_table = tmpdir_path / "ortholog_reference_table.tsv"
            ortholog_table.write_text(
                "\n".join(
                    [
                        "\t".join(
                            [
                                "input_index",
                                "input_uniprot_id",
                                "input_uniprot_name",
                                "input_gene",
                                "input_prot_family",
                                "query",
                                "source_column",
                                "status",
                                "error",
                                "human_gene_symbol",
                                "mouse_gene_symbol",
                                "macaca_fascicularis_gene_symbol",
                                "canonical_family_accession",
                                "canonical_family_name",
                                "human_refseq_accession",
                                "mouse_refseq_accession",
                                "macaca_fascicularis_refseq_accession",
                                "human_uniprot_accession",
                                "human_refseq_sequence",
                                "mouse_refseq_sequence",
                                "macaca_fascicularis_refseq_sequence",
                                "mouse_identity_to_human",
                                "macaca_fascicularis_identity_to_human",
                                "human_notes",
                                "mouse_notes",
                                "macaca_fascicularis_notes",
                            ]
                        ),
                        "\t".join(
                            [
                                "1",
                                "P00001",
                                "TEST_HUMAN",
                                "TEST",
                                "Test family",
                                "TEST_HUMAN",
                                "uniprot_name",
                                "ok",
                                "",
                                "TEST",
                                "Test",
                                "TEST",
                                "IPR_FAMILY",
                                "Broad family",
                                "NP_TEST_HUMAN",
                                "NP_TEST_MOUSE",
                                "XP_TEST_MACFA",
                                "P00001",
                                "MMAAACCCRRRR",
                                "MMAAACCCRRRK",
                                "MMAAACCCRRRQ",
                                "100.0",
                                "100.0",
                                "precomputed human",
                                "precomputed mouse",
                                "precomputed macfa",
                            ]
                        ),
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            family_dir = tmpdir_path / "family_alignments"
            family_dir.mkdir(parents=True, exist_ok=True)
            family_json = family_dir / "ipr_family.json"
            family_json.write_text(
                json.dumps(
                    {
                        "family_accession": "IPR_FAMILY",
                        "family_name": "Broad family",
                        "identity_matrix_labels": ["TEST_HUMAN", "TEST2_HUMAN"],
                        "identity_matrix": [[100.0, 50.0], [100.0, 100.0]],
                        "members": [
                            {
                                "input_index": 1,
                                "input_uniprot_id": "P00001",
                                "input_uniprot_name": "TEST_HUMAN",
                                "input_gene": "TEST",
                                "input_prot_family": "Test family",
                                "accession": "P00001",
                                "entry_name": "TEST_HUMAN",
                                "gene_symbol": "TEST",
                                "protein_name": "Test receptor",
                                "family_accession": "IPR_FAMILY",
                                "family_name": "Broad family",
                                "ectodomain_start": 3,
                                "ectodomain_end": 8,
                                "ectodomain_sequence": "AAACCC",
                                "matrix_label": "TEST_HUMAN",
                            },
                            {
                                "input_index": 2,
                                "input_uniprot_id": "P99999",
                                "input_uniprot_name": "TEST2_HUMAN",
                                "input_gene": "TEST2",
                                "input_prot_family": "Test family",
                                "accession": "P99999",
                                "entry_name": "TEST2_HUMAN",
                                "gene_symbol": "TEST2",
                                "protein_name": "Test receptor 2",
                                "family_accession": "IPR_FAMILY",
                                "family_name": "Broad family",
                                "ectodomain_start": 3,
                                "ectodomain_end": 5,
                                "ectodomain_sequence": "AAA",
                                "matrix_label": "TEST2_HUMAN",
                            },
                        ],
                        "pairwise_alignments": [],
                    },
                    indent=2,
                )
                + "\n",
                encoding="utf-8",
            )
            family_index = family_dir / "family_alignment_index.tsv"
            family_index.write_text(
                "family_accession\tfamily_name\tmember_count\tjson_path\tfasta_path\n"
                f"IPR_FAMILY\tBroad family\t2\t{family_json}\t{family_dir / 'ipr_family.fasta'}\n",
                encoding="utf-8",
            )
            config.precomputed_ortholog_table_path = ortholog_table
            config.precomputed_family_alignment_index_path = family_index
            analyzer = AntigenAnalyzer(
                config=config,
                uniprot_client=FakeUniProtClient(target, None, None),
                alphafold_client=FakeAlphaFoldClient(tmpdir_path),
                hgnc_client=RaisingHGNCClient(),
                blast_client=FakeBlastClient(),
                refseq_client=FakeRefSeqClient("MMAAACCCRRRQ"),
                interpro_client=FakeInterProClient(),
                complex_portal_client=FakeComplexPortalClient(),
                asset_generator=FakeAssetGenerator(),
            )
            report = analyzer.analyze_target("TEST", output_dir=tmpdir_path)
        self.assertIsNotNone(report.canonical_family)
        assert report.canonical_family is not None
        self.assertEqual(report.canonical_family.accession, "IPR_FAMILY")
        self.assertIsNotNone(report.family_context)
        assert report.family_context is not None
        self.assertEqual(report.family_context.source, "Precomputed InterPro family alignments")
        self.assertEqual(report.family_context.identity_matrix_labels, ["TEST_HUMAN", "TEST2_HUMAN"])
        self.assertEqual(report.family_context.identity_matrix[0][1], 50.0)
        self.assertEqual(report.family_context.identity_matrix[1][0], 100.0)
        mouse_record = next(item for item in report.ectodomain_homology if item.species == "mouse")
        self.assertEqual(mouse_record.accession, "NP_TEST_MOUSE")
        self.assertTrue(any("precomputed ortholog reference table" in note.message.lower() for note in report.notes))
        self.assertEqual(report.construct_details[0].homologs[0].accession, "NP_TEST_MOUSE")
        self.assertTrue(report.construct_details[0].homologs[0].available)

    def test_prefers_precomputed_paralog_context_over_family_alignment_or_hgnc(self) -> None:
        target = build_entry("TEST_HUMAN", "P00001", "MMAAACCCRRRR", "Homo sapiens", 9606, "TEST")
        config = AnalysisConfig(
            min_structured_segment=3,
            min_domain_size=2,
            min_construct_length=2,
            blast_species=("human",),
            generate_assets=False,
            prefer_precomputed_references=True,
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            ortholog_table = tmpdir_path / "ortholog_reference_table.tsv"
            ortholog_table.write_text(
                "input_index\tinput_uniprot_id\tinput_uniprot_name\tinput_gene\tinput_prot_family\tquery\tsource_column\tstatus\terror\thuman_gene_symbol\tcanonical_family_accession\tcanonical_family_name\thuman_uniprot_accession\thuman_refseq_sequence\n"
                "1\tP00001\tTEST_HUMAN\tTEST\tTest family\tTEST_HUMAN\tuniprot_name\tok\t\tTEST\tIPR_FAMILY\tBroad family\tP00001\tMMAAACCCRRRR\n",
                encoding="utf-8",
            )
            family_dir = tmpdir_path / "family_alignments"
            family_dir.mkdir(parents=True, exist_ok=True)
            family_index = family_dir / "family_alignment_index.tsv"
            family_index.write_text(
                "family_accession\tfamily_name\tmember_count\tjson_path\tfasta_path\n",
                encoding="utf-8",
            )
            paralog_dir = tmpdir_path / "paralogs"
            paralog_dir.mkdir(parents=True, exist_ok=True)
            paralog_json = paralog_dir / "test_human.json"
            paralog_json.write_text(
                json.dumps(
                    {
                        "target_entry_name": "TEST_HUMAN",
                        "target_accession": "P00001",
                        "target_gene_symbol": "TEST",
                        "canonical_family_accession": "IPR_FAMILY",
                        "canonical_family_name": "Broad family",
                        "family_names": ["Broader paralog set"],
                        "identity_matrix_labels": ["TEST_HUMAN", "PARA1_HUMAN"],
                        "identity_matrix": [[100.0, 70.0], [65.0, 100.0]],
                        "coverage_matrix_labels": ["TEST_HUMAN", "PARA1_HUMAN"],
                        "coverage_matrix": [[100.0, 95.0], [90.0, 100.0]],
                        "members": [
                            {
                                "gene_symbol": "TEST",
                                "gene_name": "Test protein",
                                "accession": "P00001",
                                "entry_name": "TEST_HUMAN",
                                "sequence_length": 12,
                                "ectodomain_start": 1,
                                "ectodomain_end": 12,
                                "ectodomain_length": 12,
                                "sources": ["self"],
                                "notes": ["Target protein."],
                            },
                            {
                                "gene_symbol": "PARA1",
                                "gene_name": "Paralog 1",
                                "accession": "Q00002",
                                "entry_name": "PARA1_HUMAN",
                                "sequence_length": 11,
                                "ectodomain_start": 1,
                                "ectodomain_end": 11,
                                "ectodomain_length": 11,
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
            paralog_index = paralog_dir / "paralog_index.tsv"
            paralog_index.write_text(
                "input_index\tquery\ttarget_entry_name\ttarget_accession\ttarget_gene_symbol\tcanonical_family_accession\tcanonical_family_name\tmember_count\tstatus\terror\tjson_path\n"
                f"1\tTEST_HUMAN\tTEST_HUMAN\tP00001\tTEST\tIPR_FAMILY\tBroad family\t2\tok\t\t{paralog_json}\n",
                encoding="utf-8",
            )
            config.precomputed_ortholog_table_path = ortholog_table
            config.precomputed_family_alignment_index_path = family_index
            config.precomputed_paralog_index_path = paralog_index
            analyzer = AntigenAnalyzer(
                config=config,
                uniprot_client=FakeUniProtClient(target, None, None),
                alphafold_client=FakeAlphaFoldClient(tmpdir_path),
                hgnc_client=RaisingHGNCClient(),
                blast_client=FakeBlastClient(),
                refseq_client=FakeRefSeqClient("MMAAACCCRRRQ"),
                interpro_client=FakeInterProClient(),
                complex_portal_client=FakeComplexPortalClient(),
                asset_generator=FakeAssetGenerator(),
            )
            report = analyzer.analyze_target("TEST", output_dir=tmpdir_path)
        self.assertIsNotNone(report.family_context)
        assert report.family_context is not None
        self.assertEqual(report.family_context.source, "Precomputed paralog reference set")
        self.assertEqual(report.family_context.family_names, ["Broad family", "Broader paralog set"])
        self.assertEqual(report.family_context.coverage_matrix[0][1], 95.0)
        self.assertTrue(any("precomputed paralog reference set" in note.message.lower() for note in report.notes))

    def test_uses_hgnc_family_fallback_when_precomputed_target_has_only_pfam_family(self) -> None:
        target = build_entry("TEST_HUMAN", "P00001", "MMAAACCCRRRR", "Homo sapiens", 9606, "TEST")
        config = AnalysisConfig(
            min_structured_segment=3,
            min_domain_size=2,
            min_construct_length=2,
            blast_species=("human",),
            generate_assets=False,
            prefer_precomputed_references=True,
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            ortholog_table = tmpdir_path / "ortholog_reference_table.tsv"
            ortholog_table.write_text(
                "\n".join(
                    [
                        "\t".join(
                            [
                                "input_index",
                                "input_uniprot_id",
                                "input_uniprot_name",
                                "input_gene",
                                "input_prot_family",
                                "query",
                                "source_column",
                                "status",
                                "error",
                                "human_gene_symbol",
                                "mouse_gene_symbol",
                                "macaca_fascicularis_gene_symbol",
                                "canonical_family_accession",
                                "canonical_family_name",
                                "human_refseq_accession",
                                "mouse_refseq_accession",
                                "macaca_fascicularis_refseq_accession",
                                "human_uniprot_accession",
                                "human_refseq_sequence",
                                "mouse_refseq_sequence",
                                "macaca_fascicularis_refseq_sequence",
                            ]
                        ),
                        "\t".join(
                            [
                                "1",
                                "P00001",
                                "TEST_HUMAN",
                                "TEST",
                                "Test family",
                                "TEST_HUMAN",
                                "uniprot_name",
                                "ok",
                                "",
                                "TEST",
                                "Test",
                                "TEST",
                                "PF99999",
                                "Pfam-only family",
                                "NP_TEST_HUMAN",
                                "NP_TEST_MOUSE",
                                "XP_TEST_MACFA",
                                "P00001",
                                "MMAAACCCRRRR",
                                "MMAAACCCRRRK",
                                "MMAAACCCRRRQ",
                            ]
                        ),
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            family_dir = tmpdir_path / "family_alignments"
            family_dir.mkdir(parents=True, exist_ok=True)
            family_index = family_dir / "family_alignment_index.tsv"
            family_index.write_text(
                "family_accession\tfamily_name\tmember_count\tjson_path\tfasta_path\n",
                encoding="utf-8",
            )
            config.precomputed_ortholog_table_path = ortholog_table
            config.precomputed_family_alignment_index_path = family_index
            analyzer = AntigenAnalyzer(
                config=config,
                uniprot_client=FakeUniProtClient(target, None, None),
                alphafold_client=FakeAlphaFoldClient(tmpdir_path),
                hgnc_client=FakeHGNCClient(),
                blast_client=FakeBlastClient(),
                refseq_client=FakeRefSeqClient("MMAAACCCRRRQ"),
                interpro_client=FakeInterProClient(),
                complex_portal_client=FakeComplexPortalClient(),
                asset_generator=FakeAssetGenerator(),
            )
            report = analyzer.analyze_target("TEST", output_dir=tmpdir_path)
        self.assertIsNone(report.canonical_family)
        self.assertIsNotNone(report.family_context)
        assert report.family_context is not None
        self.assertEqual(report.family_context.source, "HGNC")
        self.assertEqual(report.family_context.family_names, ["Test family"])

    def test_skips_hgnc_family_fallback_when_precomputed_interpro_family_is_present(self) -> None:
        target = build_entry("TEST_HUMAN", "P00001", "MMAAACCCRRRR", "Homo sapiens", 9606, "TEST")
        config = AnalysisConfig(
            min_structured_segment=3,
            min_domain_size=2,
            min_construct_length=2,
            blast_species=("human",),
            generate_assets=False,
            prefer_precomputed_references=True,
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            ortholog_table = tmpdir_path / "ortholog_reference_table.tsv"
            ortholog_table.write_text(
                "\n".join(
                    [
                        "\t".join(
                            [
                                "input_index",
                                "input_uniprot_id",
                                "input_uniprot_name",
                                "input_gene",
                                "input_prot_family",
                                "query",
                                "source_column",
                                "status",
                                "error",
                                "human_gene_symbol",
                                "canonical_family_accession",
                                "canonical_family_name",
                                "human_uniprot_accession",
                                "human_refseq_sequence",
                            ]
                        ),
                        "\t".join(
                            [
                                "1",
                                "P00001",
                                "TEST_HUMAN",
                                "TEST",
                                "Test family",
                                "TEST_HUMAN",
                                "uniprot_name",
                                "ok",
                                "",
                                "TEST",
                                "IPR_FAMILY",
                                "Broad family",
                                "P00001",
                                "MMAAACCCRRRR",
                            ]
                        ),
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            family_dir = tmpdir_path / "family_alignments"
            family_dir.mkdir(parents=True, exist_ok=True)
            family_index = family_dir / "family_alignment_index.tsv"
            family_index.write_text(
                "family_accession\tfamily_name\tmember_count\tjson_path\tfasta_path\n",
                encoding="utf-8",
            )
            config.precomputed_ortholog_table_path = ortholog_table
            config.precomputed_family_alignment_index_path = family_index
            analyzer = AntigenAnalyzer(
                config=config,
                uniprot_client=FakeUniProtClient(target, None, None),
                alphafold_client=FakeAlphaFoldClient(tmpdir_path),
                hgnc_client=RaisingHGNCClient(),
                blast_client=FakeBlastClient(),
                refseq_client=FakeRefSeqClient("MMAAACCCRRRQ"),
                interpro_client=FakeInterProClient(),
                complex_portal_client=FakeComplexPortalClient(),
                asset_generator=FakeAssetGenerator(),
            )
            report = analyzer.analyze_target("TEST", output_dir=tmpdir_path)
        self.assertIsNotNone(report.canonical_family)
        self.assertIsNone(report.family_context)
        self.assertTrue(any("skipping live family-context fallback" in note.message.lower() for note in report.notes))

    def test_end_to_end_analysis(self) -> None:
        target = build_entry("TEST_HUMAN", "P00001", "MMAAACCCRRRR", "Homo sapiens", 9606, "TEST")
        mouse = build_entry("TEST_MOUSE", "Q00002", "MMAAACCCRRRK", "Mus musculus", 10090, "TEST")
        macfa = build_entry("TEST_MACFA", "A00003", "MMAAACCCRRRQ", "Macaca fascicularis", 9541, "TEST")
        config = AnalysisConfig(
            min_structured_segment=3,
            min_domain_size=2,
            min_construct_length=2,
            blast_species=("human",),
            generate_assets=False,
            enable_complex_portal=True,
            prefer_precomputed_references=False,
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            analyzer = AntigenAnalyzer(
                config=config,
                uniprot_client=FakeUniProtClient(target, mouse, macfa),
                alphafold_client=FakeAlphaFoldClient(tmpdir_path),
                hgnc_client=FakeHGNCClient(),
                blast_client=FakeBlastClient(),
                refseq_client=FakeRefSeqClient(macfa["sequence"]["value"]),
                interpro_client=FakeInterProClient(),
                complex_portal_client=FakeComplexPortalClient(),
                asset_generator=FakeAssetGenerator(),
            )
            report = analyzer.analyze_target("TEST", output_dir=tmpdir_path)
        self.assertEqual(report.target.entry_name, "TEST_HUMAN")
        self.assertIsNotNone(report.ectodomain)
        self.assertTrue(report.construct_recommendations)
        self.assertTrue(report.construct_details)
        self.assertEqual(report.construct_details[0].sequence, "AAACCC")
        self.assertEqual(report.construct_details[0].homologs[0].species, "mouse")
        self.assertTrue(report.construct_details[0].homologs[0].available)
        self.assertEqual(report.construct_details[0].human_entry_name, "TEST_HUMAN")
        self.assertEqual(report.construct_details[0].classification, "multi_domain_unit")
        self.assertTrue(report.interpro_annotations)
        self.assertTrue(all(annotation.start >= 3 and annotation.end <= 8 for annotation in report.interpro_annotations))
        self.assertFalse(any(annotation.accession == "IPR999999" for annotation in report.interpro_annotations))
        self.assertFalse(any(annotation.accession == "IPR_FAMILY" for annotation in report.interpro_annotations))
        self.assertFalse(any(annotation.accession == "PF00001" for annotation in report.interpro_annotations))
        self.assertFalse(any(annotation.accession == "PF_REPEAT_DUP" for annotation in report.interpro_annotations))
        self.assertFalse(any(annotation.accession == "IPR_SUPER_A" for annotation in report.interpro_annotations))
        self.assertIsNotNone(report.canonical_family)
        assert report.canonical_family is not None
        self.assertEqual(report.canonical_family.accession, "IPR_FAMILY")
        self.assertIsNotNone(report.construct_details[0].structural_metrics)
        self.assertIsNotNone(report.construct_details[0].structural_metrics.mean_plddt)
        self.assertIsNotNone(report.construct_details[0].structural_metrics.mean_intra_pae)
        self.assertIsInstance(report.construct_details[0].split_diagnostics, list)
        self.assertIsNotNone(report.construct_details[0].homologs[0].identity_to_human)
        self.assertTrue(report.ligand_interactions)
        self.assertTrue(report.construct_details[0].ligand_interactions)
        self.assertTrue(report.ptms)
        self.assertTrue(report.construct_details[0].ptms)
        self.assertTrue(any("cleavage" in warning.lower() for warning in report.construct_details[0].warnings))
        self.assertEqual(len(report.cross_reactivity_hits), 1)
        self.assertEqual(report.ectodomain_homology[0].species, "mouse")
        self.assertEqual([construct.pdb_id for construct in report.experimental_constructs], ["1ABC"])
        self.assertTrue(any(construct.classification == "complete_domain" for construct in report.construct_recommendations))
        self.assertTrue(all(construct.start >= 3 and construct.end <= 8 for construct in report.construct_recommendations))
        full_ectodomain = next(construct for construct in report.construct_recommendations if construct.name == "full_ectodomain")
        self.assertNotEqual(full_ectodomain.classification, "partial_domain")
        self.assertFalse(any("IPR999999" in warning for warning in full_ectodomain.warnings))
        # PDB-backed constructs are retained even when their boundaries duplicate
        # another construct (here full_ectodomain), so every PDB reference reaches
        # the report, and each carries its PDB id for RCSB linking.
        pdb_backed = [construct for construct in report.construct_recommendations if construct.pdb_id == "1ABC"]
        self.assertEqual(len(pdb_backed), 1)
        self.assertEqual((pdb_backed[0].start, pdb_backed[0].end), (full_ectodomain.start, full_ectodomain.end))
        self.assertTrue(pdb_backed[0].name.startswith("pdb_"))

    def test_unresolved_single_pass_target_runs_full_length_cross_reactivity(self) -> None:
        target = build_entry("TEST_HUMAN", "P00001", "MAAAACCC", "Homo sapiens", 9606, "TEST")
        mouse = build_entry("TEST_MOUSE", "Q00002", "MAAAACCK", "Mus musculus", 10090, "TEST")
        macfa = build_entry("TEST_MACFA", "A00003", "MAAAACCQ", "Macaca fascicularis", 9541, "TEST")
        config = AnalysisConfig(
            min_structured_segment=3,
            min_domain_size=2,
            min_construct_length=2,
            blast_species=("human",),
            generate_assets=False,
            prefer_precomputed_references=False,
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = AntigenAnalyzer(
                config=config,
                uniprot_client=UnresolvedTopologyUniProtClient(target, mouse, macfa),
                alphafold_client=FakeAlphaFoldClient(Path(tmpdir)),
                hgnc_client=FakeHGNCClient(),
                blast_client=FakeBlastClient(),
                refseq_client=FakeRefSeqClient(macfa["sequence"]["value"]),
                interpro_client=FakeInterProClient(),
                complex_portal_client=FakeComplexPortalClient(),
                asset_generator=FakeAssetGenerator(),
            )
            report = analyzer.analyze_target("TEST", output_dir=Path(tmpdir))
        self.assertIsNone(report.ectodomain)
        self.assertFalse(report.cross_reactivity_hits)
        self.assertEqual(len(report.full_length_cross_reactivity_hits), 1)
        self.assertEqual(report.full_length_cross_reactivity_hits[0].alignment_length, len(target["sequence"]["value"]))

    def test_blast_fallback_for_missing_same_name_ortholog(self) -> None:
        target = build_entry("TEST_HUMAN", "P00001", "MMAAACCCRRRR", "Homo sapiens", 9606, "TEST")
        mouse = build_entry("TEST_MOUSE", "Q00002", "MMAAACCCRRRK", "Mus musculus", 10090, "TEST")
        macfa = build_entry("ALTX_MACFA", "A00003", "MMAAACCCRRRQ", "Macaca fascicularis", 9541, "TEST")
        config = AnalysisConfig(
            min_structured_segment=3,
            min_domain_size=2,
            min_construct_length=2,
            blast_species=("human",),
            generate_assets=False,
            enable_complex_portal=True,
            prefer_precomputed_references=False,
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = AntigenAnalyzer(
                config=config,
                uniprot_client=FakeUniProtClient(target, mouse, None, extra_entries=[macfa]),
                alphafold_client=FakeAlphaFoldClient(Path(tmpdir)),
                hgnc_client=FakeHGNCClient(),
                blast_client=FakeBlastClient(),
                refseq_client=FakeRefSeqClient(macfa["sequence"]["value"]),
                interpro_client=FakeInterProClient(),
                complex_portal_client=FakeComplexPortalClient(),
                asset_generator=FakeAssetGenerator(),
            )
            report = analyzer.analyze_target("TEST", output_dir=tmpdir)
        macfa_record = next(item for item in report.ectodomain_homology if item.species == "macaca_fascicularis")
        self.assertTrue(macfa_record.available)
        self.assertEqual(macfa_record.entry_name, "XP_TEST_1")
        self.assertTrue(any("RefSeq canonical protein candidate" in note for note in macfa_record.notes))

    def test_construct_warns_on_unpaired_cysteine(self) -> None:
        target = build_entry("TEST_HUMAN", "P00001", "MMAAACCCRRRR", "Homo sapiens", 9606, "TEST")
        mouse = build_entry("TEST_MOUSE", "Q00002", "MMAAACCCRRRK", "Mus musculus", 10090, "TEST")
        macfa = build_entry("TEST_MACFA", "A00003", "MMAAACCCRRRQ", "Macaca fascicularis", 9541, "TEST")
        config = AnalysisConfig(generate_assets=False, prefer_precomputed_references=False)
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            analyzer = AntigenAnalyzer(
                config=config,
                uniprot_client=FakeUniProtClient(target, mouse, macfa),
                alphafold_client=FakeAlphaFoldClient(tmpdir_path),
                hgnc_client=FakeHGNCClient(),
                blast_client=FakeBlastClient(),
                refseq_client=FakeRefSeqClient(macfa["sequence"]["value"]),
                interpro_client=FakeInterProClient(),
                complex_portal_client=FakeComplexPortalClient(),
                asset_generator=FakeAssetGenerator(),
            )
            pdb_path, _ = analyzer.alphafold_client.ensure_artifacts("P00001")
            structure = parse_alphafold_pdb(pdb_path)
            details = analyzer._build_construct_details(
                target_entry_name="TEST_HUMAN",
                target_sequence=target["sequence"]["value"],
                constructs=[
                    ConstructSuggestion(
                        name="single_cys_construct",
                        start=3,
                        end=4,
                        score=50.0,
                        rationale="Test construct with one cysteine.",
                    )
                ],
                homolog_context={},
                ectodomain=Region(start=3, end=8, label="ecto", source="test"),
                ligand_interactions=[],
                interpro_annotations=[],
                structure_atoms=structure["atoms"],
                residue_names=structure["residue_names"],
            )
        self.assertEqual(len(details), 1)
        self.assertTrue(any("Contains unpaired cysteine(s): C3." == warning for warning in details[0].warnings))
        self.assertEqual([finding.position for finding in details[0].cysteine_analysis if finding.paired_with is None], [3])

    def test_suggest_constructs_for_multipass_target_includes_membrane_expression_variants(self) -> None:
        config = AnalysisConfig(min_construct_length=50, generate_assets=False, prefer_precomputed_references=False)
        analyzer = AntigenAnalyzer(config=config)
        features = [
            Feature(type="SIGNAL", start=1, end=24, description="Signal peptide"),
            Feature(type="TOPO_DOM", start=25, end=180, description="Extracellular"),
            Feature(type="TRANSMEM", start=181, end=201, description="Helical"),
            Feature(type="TOPO_DOM", start=202, end=220, description="Cytoplasmic"),
            Feature(type="TRANSMEM", start=221, end=241, description="Helical"),
            Feature(type="TOPO_DOM", start=242, end=260, description="Cytoplasmic"),
            Feature(type="TRANSMEM", start=261, end=281, description="Helical"),
            Feature(type="TOPO_DOM", start=282, end=300, description="Cytoplasmic"),
            Feature(type="TRANSMEM", start=301, end=321, description="Helical"),
            Feature(type="TOPO_DOM", start=322, end=340, description="Cytoplasmic"),
            Feature(type="TRANSMEM", start=341, end=361, description="Helical"),
            Feature(type="TOPO_DOM", start=362, end=380, description="Cytoplasmic"),
            Feature(type="TRANSMEM", start=381, end=401, description="Helical"),
            Feature(type="TOPO_DOM", start=402, end=420, description="Cytoplasmic"),
            Feature(type="TRANSMEM", start=421, end=441, description="Helical"),
            Feature(type="TOPO_DOM", start=442, end=520, description="Cytoplasmic"),
        ]
        topology = derive_ectodomain(features, 520, config)
        constructs = analyzer._suggest_constructs(
            features=features,
            topology=topology.topology,
            ectodomain=topology.ectodomain,
            sequence_length=520,
            structural_regions=[],
            experimental_constructs=[],
            cysteine_analysis=[],
            interpro_annotations=[],
            residue_plddt={position: 40.0 for position in range(454, 521)},
        )
        names = {construct.name: construct for construct in constructs}
        self.assertIn("full_ectodomain", names)
        self.assertIn("membrane_expression_c_tail_trimmed", names)
        self.assertIn("membrane_expression_full_length", names)
        trimmed = names["membrane_expression_c_tail_trimmed"]
        full_length = names["membrane_expression_full_length"]
        self.assertEqual(trimmed.classification, "membrane_expression")
        self.assertEqual((trimmed.start, trimmed.end), (1, 453))
        self.assertEqual(full_length.classification, "membrane_expression")
        self.assertEqual((full_length.start, full_length.end), (1, 520))
        self.assertTrue(any("mean pLDDT" in item for item in trimmed.evidence))

    def test_short_untrimmed_gpi_region_keeps_full_construct(self) -> None:
        config = AnalysisConfig(min_construct_length=50, generate_assets=False, prefer_precomputed_references=False)
        analyzer = AntigenAnalyzer(config=config)
        features = [
            Feature(type="SIGNAL", start=1, end=24, description="Signal peptide"),
            Feature(type="LIPIDATION", start=36, end=36, description="GPI-anchor amidated serine"),
        ]
        topology = derive_ectodomain(features, 61, config)
        constructs = analyzer._suggest_constructs(
            features=features,
            topology=topology.topology,
            ectodomain=topology.ectodomain,
            sequence_length=61,
            structural_regions=[],
            experimental_constructs=[],
            cysteine_analysis=[],
            interpro_annotations=[],
        )

        self.assertEqual([(item.name, item.start, item.end) for item in constructs], [("full_ectodomain", 25, 61)])

    def test_compact_multipass_target_includes_membrane_expression_variants(self) -> None:
        config = AnalysisConfig(min_construct_length=50, generate_assets=False, prefer_precomputed_references=False)
        analyzer = AntigenAnalyzer(config=config)
        features = [
            Feature(type="TOPO_DOM", start=1, end=68, description="Extracellular"),
            Feature(type="TRANSMEM", start=69, end=93, description="Helical; Name=1"),
            Feature(type="TOPO_DOM", start=94, end=106, description="Cytoplasmic"),
            Feature(type="TRANSMEM", start=107, end=131, description="Helical; Name=2"),
            Feature(type="TOPO_DOM", start=132, end=142, description="Extracellular"),
            Feature(type="TRANSMEM", start=143, end=165, description="Helical; Name=3"),
            Feature(type="TOPO_DOM", start=166, end=185, description="Cytoplasmic"),
            Feature(type="TRANSMEM", start=186, end=207, description="Helical; Name=4"),
            Feature(type="TOPO_DOM", start=208, end=230, description="Extracellular"),
            Feature(type="TRANSMEM", start=231, end=255, description="Helical; Name=5"),
            Feature(type="TOPO_DOM", start=256, end=279, description="Cytoplasmic"),
            Feature(type="TRANSMEM", start=280, end=306, description="Helical; Name=6"),
            Feature(type="TOPO_DOM", start=307, end=314, description="Extracellular"),
            Feature(type="TRANSMEM", start=315, end=338, description="Helical; Name=7"),
            Feature(type="TOPO_DOM", start=339, end=400, description="Cytoplasmic"),
        ]
        topology = derive_ectodomain(features, 400, config)
        constructs = analyzer._suggest_constructs(
            features=features,
            topology=topology.topology,
            ectodomain=topology.ectodomain,
            sequence_length=400,
            structural_regions=[],
            experimental_constructs=[],
            cysteine_analysis=[],
            interpro_annotations=[],
            residue_plddt={position: 40.0 for position in range(339, 401)},
        )
        names = {construct.name: construct for construct in constructs}
        self.assertEqual(topology.topology.topology_class, "multipass_compact")
        self.assertNotIn("full_ectodomain", names)
        self.assertIn("membrane_expression_c_tail_trimmed", names)
        self.assertIn("membrane_expression_full_length", names)
        self.assertEqual((names["membrane_expression_c_tail_trimmed"].start, names["membrane_expression_c_tail_trimmed"].end), (1, 350))
        self.assertEqual(names["membrane_expression_c_tail_trimmed"].classification, "membrane_expression")

    def test_compact_multipass_target_includes_pdb_backed_membrane_constructs(self) -> None:
        config = AnalysisConfig(min_construct_length=50, generate_assets=False, prefer_precomputed_references=False)
        analyzer = AntigenAnalyzer(config=config)
        features = [
            Feature(type="TOPO_DOM", start=1, end=68, description="Extracellular"),
            Feature(type="TRANSMEM", start=69, end=93, description="Helical; Name=1"),
            Feature(type="TOPO_DOM", start=94, end=106, description="Cytoplasmic"),
            Feature(type="TRANSMEM", start=107, end=131, description="Helical; Name=2"),
            Feature(type="TOPO_DOM", start=132, end=142, description="Extracellular"),
            Feature(type="TRANSMEM", start=143, end=165, description="Helical; Name=3"),
            Feature(type="TOPO_DOM", start=166, end=185, description="Cytoplasmic"),
            Feature(type="TRANSMEM", start=186, end=207, description="Helical; Name=4"),
            Feature(type="TOPO_DOM", start=208, end=230, description="Extracellular"),
            Feature(type="TRANSMEM", start=231, end=255, description="Helical; Name=5"),
            Feature(type="TOPO_DOM", start=256, end=279, description="Cytoplasmic"),
            Feature(type="TRANSMEM", start=280, end=306, description="Helical; Name=6"),
            Feature(type="TOPO_DOM", start=307, end=314, description="Extracellular"),
            Feature(type="TRANSMEM", start=315, end=338, description="Helical; Name=7"),
            Feature(type="TOPO_DOM", start=339, end=400, description="Cytoplasmic"),
        ]
        topology = derive_ectodomain(features, 400, config)
        constructs = analyzer._suggest_constructs(
            features=features,
            topology=topology.topology,
            ectodomain=topology.ectodomain,
            sequence_length=400,
            structural_regions=[],
            experimental_constructs=[
                ExperimentalConstruct(
                    pdb_id="6ABC",
                    method="EM",
                    resolution="3.1 A",
                    chains=[Region(start=45, end=360, label="Chain A/B", source="PDB", metadata={"chain": "A/B"})],
                )
            ],
            cysteine_analysis=[],
            interpro_annotations=[],
        )
        pdb_construct = next(construct for construct in constructs if construct.name == "pdb_6abc_a_b")
        self.assertEqual((pdb_construct.start, pdb_construct.end), (45, 360))
        self.assertEqual(pdb_construct.classification, "membrane_expression")
        self.assertTrue(any("PDB 6ABC" in item for item in pdb_construct.evidence))

    def test_compact_multipass_keeps_full_length_baseline_when_pdb_has_same_bounds(self) -> None:
        config = AnalysisConfig(min_construct_length=50, generate_assets=False, prefer_precomputed_references=False)
        analyzer = AntigenAnalyzer(config=config)
        features = [
            Feature(type="TOPO_DOM", start=1, end=60, description="Extracellular"),
            Feature(type="TRANSMEM", start=61, end=82, description="Helical; Name=1"),
            Feature(type="TOPO_DOM", start=83, end=92, description="Cytoplasmic"),
            Feature(type="TRANSMEM", start=93, end=114, description="Helical; Name=2"),
            Feature(type="TOPO_DOM", start=115, end=130, description="Cytoplasmic"),
            Feature(type="TRANSMEM", start=131, end=152, description="Helical; Name=3"),
            Feature(type="TOPO_DOM", start=153, end=190, description="Cytoplasmic"),
        ]
        topology = derive_ectodomain(features, 190, config)
        constructs = analyzer._suggest_constructs(
            features=features,
            topology=topology.topology,
            ectodomain=topology.ectodomain,
            sequence_length=190,
            structural_regions=[],
            experimental_constructs=[
                ExperimentalConstruct(
                    pdb_id="7FUL",
                    method="EM",
                    resolution="3.0 A",
                    chains=[Region(start=1, end=190, label="Chain A", source="PDB", metadata={"chain": "A"})],
                )
            ],
            cysteine_analysis=[],
            interpro_annotations=[],
        )
        names = {construct.name: construct for construct in constructs}
        self.assertIn("membrane_expression_full_length", names)
        self.assertIn("pdb_7ful_a", names)
        self.assertEqual((names["membrane_expression_full_length"].start, names["membrane_expression_full_length"].end), (1, 190))
        self.assertEqual((names["pdb_7ful_a"].start, names["pdb_7ful_a"].end), (1, 190))

    def test_mixed_multipass_target_keeps_soluble_and_membrane_pdb_constructs(self) -> None:
        config = AnalysisConfig(min_construct_length=50, generate_assets=False, prefer_precomputed_references=False)
        analyzer = AntigenAnalyzer(config=config)
        features = [
            Feature(type="SIGNAL", start=1, end=24, description="Signal peptide"),
            Feature(type="TOPO_DOM", start=25, end=180, description="Extracellular"),
            Feature(type="TRANSMEM", start=181, end=201, description="Helical; Name=1"),
            Feature(type="TOPO_DOM", start=202, end=220, description="Cytoplasmic"),
            Feature(type="TRANSMEM", start=221, end=241, description="Helical; Name=2"),
            Feature(type="TOPO_DOM", start=242, end=280, description="Cytoplasmic"),
        ]
        topology = derive_ectodomain(features, 280, config)
        constructs = analyzer._suggest_constructs(
            features=features,
            topology=topology.topology,
            ectodomain=topology.ectodomain,
            sequence_length=280,
            structural_regions=[],
            experimental_constructs=[
                ExperimentalConstruct(
                    pdb_id="1ECD",
                    method="X-ray",
                    resolution="2.0 A",
                    chains=[Region(start=30, end=175, label="Chain E", source="PDB", metadata={"chain": "E"})],
                ),
                ExperimentalConstruct(
                    pdb_id="7MEM",
                    method="EM",
                    resolution="3.5 A",
                    chains=[Region(start=25, end=260, label="Chain M", source="PDB", metadata={"chain": "M"})],
                ),
            ],
            cysteine_analysis=[],
            interpro_annotations=[],
        )
        names = {construct.name: construct for construct in constructs}
        self.assertEqual(topology.topology.topology_class, "multipass_mixed")
        self.assertNotEqual(names["pdb_1ecd_e"].classification, "membrane_expression")
        self.assertEqual(names["pdb_7mem_m"].classification, "membrane_expression")

    def test_non_multipass_pdb_constructs_still_require_design_region_overlap(self) -> None:
        config = AnalysisConfig(min_construct_length=2, generate_assets=False, prefer_precomputed_references=False)
        analyzer = AntigenAnalyzer(config=config)
        features = [
            Feature(type="TOPO_DOM", start=1, end=80, description="Extracellular"),
            Feature(type="TRANSMEM", start=81, end=101, description="Helical"),
            Feature(type="TOPO_DOM", start=102, end=140, description="Cytoplasmic"),
        ]
        topology = derive_ectodomain(features, 140, config)
        constructs = analyzer._suggest_constructs(
            features=features,
            topology=topology.topology,
            ectodomain=topology.ectodomain,
            sequence_length=140,
            structural_regions=[],
            experimental_constructs=[
                ExperimentalConstruct(
                    pdb_id="1OUT",
                    method="X-ray",
                    resolution="2.0 A",
                    chains=[Region(start=105, end=120, label="Chain B", source="PDB", metadata={"chain": "B"})],
                )
            ],
            cysteine_analysis=[],
            interpro_annotations=[],
        )
        self.assertFalse(any(construct.name.startswith("pdb_") for construct in constructs))

    def test_multipass_cytoplasmic_n_tail_trim_variant(self) -> None:
        config = AnalysisConfig(min_construct_length=50, generate_assets=False, prefer_precomputed_references=False)
        analyzer = AntigenAnalyzer(config=config)
        features = [
            Feature(type="TOPO_DOM", start=1, end=60, description="Cytoplasmic"),
            Feature(type="TRANSMEM", start=61, end=82, description="Helical; Name=1"),
            Feature(type="TOPO_DOM", start=83, end=92, description="Extracellular"),
            Feature(type="TRANSMEM", start=93, end=114, description="Helical; Name=2"),
            Feature(type="TOPO_DOM", start=115, end=130, description="Cytoplasmic"),
            Feature(type="TRANSMEM", start=131, end=152, description="Helical; Name=3"),
            Feature(type="TOPO_DOM", start=153, end=190, description="Extracellular"),
        ]
        topology = derive_ectodomain(features, 190, config)
        constructs = analyzer._suggest_constructs(
            features=features,
            topology=topology.topology,
            ectodomain=topology.ectodomain,
            sequence_length=190,
            structural_regions=[],
            experimental_constructs=[],
            cysteine_analysis=[],
            interpro_annotations=[],
            residue_plddt={position: 35.0 for position in range(1, 49)},
        )
        names = {construct.name: construct for construct in constructs}
        self.assertIn("membrane_expression_n_tail_trimmed", names)
        self.assertEqual((names["membrane_expression_n_tail_trimmed"].start, names["membrane_expression_n_tail_trimmed"].end), (49, 190))
        self.assertTrue(any("n-terminal cytoplasmic region 1-48" in item for item in names["membrane_expression_n_tail_trimmed"].evidence))

    def test_multipass_both_terminal_trim_variant(self) -> None:
        config = AnalysisConfig(min_construct_length=50, generate_assets=False, prefer_precomputed_references=False)
        analyzer = AntigenAnalyzer(config=config)
        features = [
            Feature(type="TOPO_DOM", start=1, end=60, description="Cytoplasmic"),
            Feature(type="TRANSMEM", start=61, end=82, description="Helical; Name=1"),
            Feature(type="TOPO_DOM", start=83, end=92, description="Extracellular"),
            Feature(type="TRANSMEM", start=93, end=114, description="Helical; Name=2"),
            Feature(type="TOPO_DOM", start=115, end=130, description="Cytoplasmic"),
            Feature(type="TRANSMEM", start=131, end=152, description="Helical; Name=3"),
            Feature(type="TOPO_DOM", start=153, end=220, description="Cytoplasmic"),
        ]
        topology = derive_ectodomain(features, 220, config)
        low_plddt = {position: 35.0 for position in list(range(1, 49)) + list(range(165, 221))}
        constructs = analyzer._suggest_constructs(
            features=features,
            topology=topology.topology,
            ectodomain=topology.ectodomain,
            sequence_length=220,
            structural_regions=[],
            experimental_constructs=[],
            cysteine_analysis=[],
            interpro_annotations=[],
            residue_plddt=low_plddt,
        )
        names = {construct.name: construct for construct in constructs}
        self.assertIn("membrane_expression_n_tail_trimmed", names)
        self.assertIn("membrane_expression_c_tail_trimmed", names)
        self.assertIn("membrane_expression_terminal_trimmed", names)
        self.assertEqual((names["membrane_expression_terminal_trimmed"].start, names["membrane_expression_terminal_trimmed"].end), (49, 164))

    def test_multipass_extracellular_n_terminus_is_not_silently_trimmed(self) -> None:
        config = AnalysisConfig(min_construct_length=50, generate_assets=False, prefer_precomputed_references=False)
        analyzer = AntigenAnalyzer(config=config)
        features = [
            Feature(type="TOPO_DOM", start=1, end=78, description="Extracellular"),
            Feature(type="TRANSMEM", start=79, end=100, description="Helical; Name=1"),
            Feature(type="TOPO_DOM", start=101, end=120, description="Cytoplasmic"),
            Feature(type="TRANSMEM", start=121, end=142, description="Helical; Name=2"),
            Feature(type="TOPO_DOM", start=143, end=180, description="Cytoplasmic"),
        ]
        topology = derive_ectodomain(features, 180, config)
        constructs = analyzer._suggest_constructs(
            features=features,
            topology=topology.topology,
            ectodomain=topology.ectodomain,
            sequence_length=180,
            structural_regions=[],
            experimental_constructs=[],
            cysteine_analysis=[],
            interpro_annotations=[],
            residue_plddt={position: 30.0 for position in range(1, 67)},
        )
        names = {construct.name: construct for construct in constructs}
        self.assertNotIn("membrane_expression_n_tail_trimmed", names)
        self.assertIn("membrane_expression_full_length", names)
        self.assertTrue(any("topology classifies it as extracellular" in warning for warning in names["membrane_expression_full_length"].warnings))

    def test_gpcrdb_annotation_preserves_helix8_when_trimming_c_tail(self) -> None:
        config = AnalysisConfig(min_construct_length=50, generate_assets=False, prefer_precomputed_references=False)
        analyzer = AntigenAnalyzer(config=config)
        features = [
            Feature(type="TOPO_DOM", start=1, end=68, description="Extracellular"),
            Feature(type="TRANSMEM", start=69, end=93, description="Helical; Name=1"),
            Feature(type="TOPO_DOM", start=94, end=106, description="Cytoplasmic"),
            Feature(type="TRANSMEM", start=107, end=131, description="Helical; Name=2"),
            Feature(type="TOPO_DOM", start=132, end=142, description="Extracellular"),
            Feature(type="TRANSMEM", start=143, end=165, description="Helical; Name=3"),
            Feature(type="TOPO_DOM", start=166, end=185, description="Cytoplasmic"),
            Feature(type="TRANSMEM", start=186, end=207, description="Helical; Name=4"),
            Feature(type="TOPO_DOM", start=208, end=230, description="Extracellular"),
            Feature(type="TRANSMEM", start=231, end=255, description="Helical; Name=5"),
            Feature(type="TOPO_DOM", start=256, end=279, description="Cytoplasmic"),
            Feature(type="TRANSMEM", start=280, end=306, description="Helical; Name=6"),
            Feature(type="TOPO_DOM", start=307, end=314, description="Extracellular"),
            Feature(type="TRANSMEM", start=315, end=338, description="Helical; Name=7"),
            Feature(type="TOPO_DOM", start=339, end=400, description="Cytoplasmic"),
        ]
        topology = derive_ectodomain(features, 400, config)
        gpcr_annotation = GPCRdbAnnotation(
            entry_name="test_human",
            accession="P00001",
            name="Test GPCR",
            residue_numbering_scheme="GPCRdb(A)",
            segments=[GPCRdbSegment(name="H8", start=340, end=350)],
        )
        constructs = analyzer._suggest_constructs(
            features=features,
            topology=topology.topology,
            ectodomain=topology.ectodomain,
            sequence_length=400,
            structural_regions=[],
            experimental_constructs=[],
            cysteine_analysis=[],
            interpro_annotations=[],
            residue_plddt={position: 35.0 for position in range(359, 401)},
            gpcr_annotation=gpcr_annotation,
        )
        names = {construct.name: construct for construct in constructs}
        self.assertEqual((names["membrane_expression_c_tail_trimmed"].start, names["membrane_expression_c_tail_trimmed"].end), (1, 358))
        self.assertTrue(any("helix 8" in item.lower() for item in names["membrane_expression_c_tail_trimmed"].evidence))
        self.assertIn("membrane_expression_gpcr_c_tail_keep_20", names)
        self.assertEqual((names["membrane_expression_gpcr_c_tail_keep_20"].start, names["membrane_expression_gpcr_c_tail_keep_20"].end), (1, 370))
        self.assertTrue(any("GPCR expression/stability screen variant" in item for item in names["membrane_expression_gpcr_c_tail_keep_20"].evidence))

    def test_construct_details_include_gpcr_segments_and_generic_range(self) -> None:
        analyzer = AntigenAnalyzer(config=AnalysisConfig(generate_assets=False, prefer_precomputed_references=False))
        gpcr_annotation = GPCRdbAnnotation(
            entry_name="test_human",
            accession="P00001",
            name="Test GPCR",
            residues=[
                GPCRdbResidue(sequence_number=70, amino_acid="A", protein_segment="TM1", display_generic_number="1.32x32"),
                GPCRdbResidue(sequence_number=80, amino_acid="A", protein_segment="TM1", display_generic_number="1.42x42"),
                GPCRdbResidue(sequence_number=95, amino_acid="A", protein_segment="ICL1"),
            ],
            segments=[
                GPCRdbSegment(name="TM1", start=70, end=93, generic_start="1.32x32", generic_end="1.55x55"),
                GPCRdbSegment(name="ICL1", start=94, end=106),
            ],
        )
        details = analyzer._build_construct_details(
            target_entry_name="TEST_HUMAN",
            target_sequence="A" * 120,
            constructs=[
                ConstructSuggestion(
                    name="membrane_expression_full_length",
                    start=1,
                    end=120,
                    score=86.0,
                    rationale="test",
                    classification="membrane_expression",
                )
            ],
            homolog_context={},
            ectodomain=None,
            ligand_interactions=[],
            interpro_annotations=[],
            gpcr_annotation=gpcr_annotation,
        )
        self.assertEqual([segment.name for segment in details[0].gpcr_segments], ["TM1", "ICL1"])
        self.assertEqual(details[0].gpcr_generic_range, "1.32x32 to 1.42x42")

    def test_gpcr_engineering_variants_generate_bril_and_t4l_loop_fusions(self) -> None:
        analyzer = AntigenAnalyzer(config=AnalysisConfig(generate_assets=False, prefer_precomputed_references=False))
        target = TargetRecord(
            accession="P00001",
            entry_name="OPRM_HUMAN",
            gene_symbol="OPRM1",
            protein_name="Test GPCR",
            organism="Homo sapiens",
            taxon_id=9606,
            sequence="A" * 400,
        )
        features = [
            Feature(type="TOPO_DOM", start=1, end=68, description="Extracellular"),
            Feature(type="TRANSMEM", start=69, end=93, description="Helical; Name=1"),
            Feature(type="TOPO_DOM", start=94, end=106, description="Cytoplasmic"),
            Feature(type="TRANSMEM", start=107, end=131, description="Helical; Name=2"),
            Feature(type="TOPO_DOM", start=132, end=142, description="Extracellular"),
            Feature(type="TRANSMEM", start=143, end=165, description="Helical; Name=3"),
            Feature(type="TOPO_DOM", start=166, end=185, description="Cytoplasmic"),
            Feature(type="TRANSMEM", start=186, end=207, description="Helical; Name=4"),
            Feature(type="TOPO_DOM", start=208, end=230, description="Extracellular"),
            Feature(type="TRANSMEM", start=231, end=255, description="Helical; Name=5"),
            Feature(type="TOPO_DOM", start=256, end=279, description="Cytoplasmic"),
            Feature(type="TRANSMEM", start=280, end=306, description="Helical; Name=6"),
            Feature(type="TOPO_DOM", start=307, end=314, description="Extracellular"),
            Feature(type="TRANSMEM", start=315, end=338, description="Helical; Name=7"),
            Feature(type="TOPO_DOM", start=339, end=400, description="Cytoplasmic"),
        ]
        topology = derive_ectodomain(features, 400, analyzer.config).topology
        gpcr_annotation = GPCRdbAnnotation(
            entry_name="oprm_human",
            accession="P00001",
            name="Test GPCR",
            residue_numbering_scheme="GPCRdb(A)",
            segments=[GPCRdbSegment(name="ICL3", start=256, end=279)],
            residues=[
                GPCRdbResidue(sequence_number=256, amino_acid="A", protein_segment="ICL3", display_generic_number="5.70x70"),
                GPCRdbResidue(sequence_number=279, amino_acid="A", protein_segment="ICL3", display_generic_number="6.20x20"),
            ],
        )

        variants = analyzer._suggest_gpcr_engineering_variants(
            target=target,
            topology=topology,
            features=features,
            ligand_interactions=[],
            ptms=[],
            gpcr_annotation=gpcr_annotation,
        )
        by_id = {variant.cassette_id: variant for variant in variants}

        self.assertIn("bril", by_id)
        self.assertIn("t4l", by_id)
        self.assertEqual((by_id["bril"].start, by_id["bril"].end), (256, 279))
        self.assertEqual((by_id["bril"].replaced_start, by_id["bril"].replaced_end), (261, 274))
        self.assertIn("GGSG", by_id["bril"].sequence or "")
        self.assertGreater(by_id["t4l"].engineered_length or 0, len(target.sequence))
        self.assertTrue(any("foreign cassette" in warning for warning in by_id["bril"].warnings))
        self.assertNotIn("mini_gs", by_id)
        self.assertEqual(by_id["bril"].strategy, "ICL3 cassette replacement")
        self.assertEqual(by_id["bril"].loop_source, "GPCRdb ICL3")
        self.assertIn("_ICL3_fusion", by_id["bril"].name)

        topology_only_variants = analyzer._suggest_gpcr_engineering_variants(
            target=target,
            topology=topology,
            features=features,
            ligand_interactions=[],
            ptms=[],
            gpcr_annotation=None,
        )
        self.assertTrue(topology_only_variants)
        self.assertTrue(
            all("ICL3" not in variant.name for variant in topology_only_variants)
        )
        self.assertTrue(
            all(
                variant.strategy == "Intracellular-loop cassette replacement"
                for variant in topology_only_variants
            )
        )
        self.assertTrue(
            all(
                variant.loop_source == "topology-derived internal cytoplasmic loop"
                for variant in topology_only_variants
            )
        )

    def test_gpcr_engineering_variants_skip_loop_fusion_when_icl3_overlaps_tm(self) -> None:
        analyzer = AntigenAnalyzer(config=AnalysisConfig(generate_assets=False, prefer_precomputed_references=False))
        target = TargetRecord(
            accession="P00001",
            entry_name="BADGPCR_HUMAN",
            gene_symbol="BADGPCR",
            protein_name="Test GPCR",
            organism="Homo sapiens",
            taxon_id=9606,
            sequence="A" * 160,
        )
        features = [
            Feature(type="TRANSMEM", start=1, end=25, description="Helical; Name=1"),
            Feature(type="TOPO_DOM", start=26, end=45, description="Cytoplasmic"),
            Feature(type="TRANSMEM", start=46, end=70, description="Helical; Name=2"),
            Feature(type="TOPO_DOM", start=71, end=90, description="Extracellular"),
            Feature(type="TRANSMEM", start=91, end=115, description="Helical; Name=3"),
            Feature(type="TOPO_DOM", start=116, end=160, description="Cytoplasmic"),
        ]
        topology = derive_ectodomain(features, 160, analyzer.config).topology
        gpcr_annotation = GPCRdbAnnotation(
            entry_name="badgpcr_human",
            accession="P00001",
            name="Test GPCR",
            segments=[GPCRdbSegment(name="ICL3", start=100, end=125)],
        )

        variants = analyzer._suggest_gpcr_engineering_variants(
            target=target,
            topology=topology,
            features=features,
            ligand_interactions=[],
            ptms=[],
            gpcr_annotation=gpcr_annotation,
        )
        skipped = next(variant for variant in variants if variant.name == "BADGPCR_HUMAN_ICL3_engineering_skipped")

        self.assertIsNone(skipped.sequence)
        self.assertTrue(any("membrane-spanning" in warning for warning in skipped.warnings))

    def test_membrane_expression_homolog_mapping_without_context_is_unavailable(self) -> None:
        analyzer = AntigenAnalyzer(config=AnalysisConfig(generate_assets=False, prefer_precomputed_references=False))
        homologs = analyzer._map_construct_to_homologs(
            construct=ConstructSuggestion(
                name="membrane_expression_trimmed",
                start=1,
                end=281,
                score=88.0,
                rationale="test",
                classification="membrane_expression",
            ),
            target_sequence="A" * 330,
            homolog_context={},
            ectodomain=Region(start=25, end=180, label="Major extracellular region", source="topology"),
        )
        self.assertEqual(len(homologs), 2)
        self.assertTrue(all(not homolog.available for homolog in homologs))
        self.assertTrue(all("no same-name uniprot entry" in homolog.notes[0].lower() for homolog in homologs))

    def test_membrane_expression_homolog_mapping_uses_full_length_context(self) -> None:
        analyzer = AntigenAnalyzer(config=AnalysisConfig(generate_assets=False, prefer_precomputed_references=False))
        mouse = TargetRecord(
            accession="NP_MOUSE",
            entry_name="NP_MOUSE",
            gene_symbol="TEST",
            protein_name="Test receptor",
            organism="Mus musculus",
            taxon_id=10090,
            sequence="MAAACCCDDDEEE",
        )
        macfa = TargetRecord(
            accession="NP_MACFA",
            entry_name="NP_MACFA",
            gene_symbol="TEST",
            protein_name="Test receptor",
            organism="Macaca fascicularis",
            taxon_id=9541,
            sequence="MAAACCCDDDEEF",
        )
        human_sequence = "MAAACCCDDDEEE"
        homolog_context = {
            "mouse": {
                "available": True,
                "target": mouse,
                "ectodomain": Region(start=1, end=len(mouse.sequence), label="full length", source="test"),
                "ectodomain_sequence": mouse.sequence,
                "alignment": __import__("agdesign2.sequence_utils", fromlist=["global_align"]).global_align(human_sequence, mouse.sequence),
                "notes": [],
            },
            "macaca_fascicularis": {
                "available": True,
                "target": macfa,
                "ectodomain": Region(start=1, end=len(macfa.sequence), label="full length", source="test"),
                "ectodomain_sequence": macfa.sequence,
                "alignment": __import__("agdesign2.sequence_utils", fromlist=["global_align"]).global_align(human_sequence, macfa.sequence),
                "notes": [],
            },
        }
        homologs = analyzer._map_construct_to_homologs(
            construct=ConstructSuggestion(
                name="membrane_expression_c_tail_trimmed",
                start=2,
                end=12,
                score=88.0,
                rationale="test",
                classification="membrane_expression",
            ),
            target_sequence=human_sequence,
            homolog_context=homolog_context,
            ectodomain=None,
        )
        self.assertTrue(all(homolog.available for homolog in homologs))
        self.assertEqual((homologs[0].start, homologs[0].end, homologs[0].sequence), (2, 12, "AAACCCDDDEE"))
        self.assertIsNotNone(homologs[0].identity_to_human)

    def test_full_length_homology_surface_identity_does_not_require_ectodomain(self) -> None:
        analyzer = AntigenAnalyzer(config=AnalysisConfig(generate_assets=False, prefer_precomputed_references=False))
        alignment = __import__("agdesign2.sequence_utils", fromlist=["global_align"]).global_align("ABCDE", "ABYXE")
        homology = [
            HomologyRecord(
                species="mouse",
                accession="NP_MOUSE",
                entry_name="NP_MOUSE",
                ectodomain_start=1,
                ectodomain_end=5,
                identity=60.0,
                coverage=100.0,
            )
        ]
        annotated = analyzer._annotate_homology_with_surface_identity(
            homology,
            homolog_context={"mouse": {"alignment": alignment}},
            ectodomain=None,
            extracellular_surface_positions={1, 3, 5},
        )
        self.assertEqual(annotated[0].extracellular_surface_aligned_positions, 3)
        self.assertAlmostEqual(annotated[0].extracellular_surface_identity or 0.0, 66.67, places=2)

    def test_generates_advanced_membrane_engineering_suggestions_for_gpcr_like_topology(self) -> None:
        analyzer = AntigenAnalyzer(config=AnalysisConfig(generate_assets=False, prefer_precomputed_references=False))
        topology = derive_ectodomain(
            [
                Feature(type="SIGNAL", start=1, end=24, description="Signal peptide"),
                Feature(type="TOPO_DOM", start=25, end=180, description="Extracellular"),
                Feature(type="TRANSMEM", start=181, end=201, description="Helical"),
                Feature(type="TOPO_DOM", start=202, end=220, description="Cytoplasmic"),
                Feature(type="TRANSMEM", start=221, end=241, description="Helical"),
                Feature(type="TOPO_DOM", start=242, end=260, description="Cytoplasmic"),
                Feature(type="TRANSMEM", start=261, end=281, description="Helical"),
                Feature(type="TOPO_DOM", start=282, end=300, description="Cytoplasmic"),
                Feature(type="TRANSMEM", start=301, end=321, description="Helical"),
                Feature(type="TOPO_DOM", start=322, end=340, description="Cytoplasmic"),
                Feature(type="TRANSMEM", start=341, end=361, description="Helical"),
                Feature(type="TOPO_DOM", start=362, end=405, description="Cytoplasmic"),
                Feature(type="TRANSMEM", start=406, end=426, description="Helical"),
                Feature(type="TOPO_DOM", start=427, end=445, description="Cytoplasmic"),
                Feature(type="TRANSMEM", start=446, end=466, description="Helical"),
                Feature(type="TOPO_DOM", start=467, end=510, description="Cytoplasmic"),
            ],
            510,
            analyzer.config,
        )
        suggestions = analyzer._suggest_advanced_membrane_engineering(
            target=TargetRecord(
                accession="QTEST1",
                entry_name="FZD1_HUMAN",
                gene_symbol="FZD1",
                protein_name="Frizzled-1",
                organism="Homo sapiens",
                taxon_id=9606,
                sequence="A" * 510,
            ),
            features=[
                Feature(type="DOMAIN", start=111, end=232, description="Frizzled domain"),
            ],
            topology=topology.topology,
            interpro_annotations=[],
        )
        categories = {item.category for item in suggestions}
        self.assertIn("gpcr_conservative_engineering", categories)
        self.assertIn("gpcr_loop_engineering", categories)
        self.assertIn("gpcr_tail_engineering", categories)

    def test_canonicalizes_macaque_cross_reactivity_hits_to_refseq(self) -> None:
        analyzer = AntigenAnalyzer(
            config=AnalysisConfig(generate_assets=False, prefer_precomputed_references=False),
            uniprot_client=FakeUniProtClient(
                build_entry("TEST_HUMAN", "P00001", "MMAAACCCRRRR", "Homo sapiens", 9606, "TEST"),
                None,
                None,
            ),
            alphafold_client=None,
            hgnc_client=FakeHGNCClient(),
            blast_client=FakeBlastClient(),
            refseq_client=FakeRefSeqClient("MMAAACCCRRRR"),
            interpro_client=FakeInterProClient(),
            complex_portal_client=FakeComplexPortalClient(),
            asset_generator=FakeAssetGenerator(),
        )
        hits = [
            BlastHit(
                subject_id="tr|A0A2K5WK39|A0A2K5WK39_MACFA",
                description="tr|A0A2K5WK39|A0A2K5WK39_MACFA Receptor protein-tyrosine kinase OS=Macaca fascicularis OX=9541 GN=EGFR",
                species="Macaca fascicularis",
                identity=98.7,
                coverage=100.0,
                alignment_length=10,
                evalue=1e-50,
                bitscore=300.0,
                query_start=1,
                query_end=10,
                subject_start=1,
                subject_end=10,
            ),
            BlastHit(
                subject_id="sp|Q01279|EGFR_MOUSE",
                description="sp|Q01279|EGFR_MOUSE Epidermal growth factor receptor OS=Mus musculus OX=10090 GN=Egfr",
                species="Mus musculus",
                identity=88.7,
                coverage=100.0,
                alignment_length=10,
                evalue=1e-30,
                bitscore=250.0,
                query_start=1,
                query_end=10,
                subject_start=1,
                subject_end=10,
            ),
        ]
        normalized = analyzer._canonicalize_macaca_cross_reactivity_hits("MMAAACCCRR", hits)
        self.assertEqual(normalized[0].subject_id, "XP_TEST_1")
        self.assertEqual(normalized[0].species, "Macaca fascicularis")
        self.assertIn("GN=EGFR", normalized[0].description)
        self.assertEqual(normalized[1].subject_id, "sp|Q01279|EGFR_MOUSE")

    def test_reuses_macaque_refseq_result_for_same_gene(self) -> None:
        refseq_client = FakeRefSeqClient("MMAAACCCRRRR")
        analyzer = AntigenAnalyzer(
            config=AnalysisConfig(generate_assets=False, prefer_precomputed_references=False),
            refseq_client=refseq_client,
        )
        hits = [
            BlastHit(
                subject_id="tr|A0A2K5WK39|A0A2K5WK39_MACFA",
                description=(
                    "tr|A0A2K5WK39|A0A2K5WK39_MACFA Receptor "
                    "OS=Macaca fascicularis OX=9541 GN=EGFR"
                ),
                species="Macaca fascicularis",
                identity=98.7,
                coverage=100.0,
                alignment_length=10,
                evalue=1e-50,
                bitscore=300.0,
                query_start=1,
                query_end=10,
                subject_start=1,
                subject_end=10,
            ),
        ]

        analyzer._canonicalize_macaca_cross_reactivity_hits("MMAAACCCRR", hits)
        analyzer._canonicalize_macaca_cross_reactivity_hits("MMAAACCCRR", hits)

        self.assertEqual(
            refseq_client.calls,
            [("EGFR", "Macaca fascicularis")],
        )

    def test_reuses_macaque_refseq_miss_for_same_gene(self) -> None:
        refseq_client = FakeRefSeqClient(None)
        analyzer = AntigenAnalyzer(
            config=AnalysisConfig(generate_assets=False, prefer_precomputed_references=False),
            refseq_client=refseq_client,
        )
        hits = [
            BlastHit(
                subject_id="tr|TEST|TEST_MACFA",
                description="tr|TEST|TEST_MACFA OS=Macaca fascicularis GN=MISSING",
                species="Macaca fascicularis",
                identity=80.0,
                coverage=90.0,
                alignment_length=9,
                evalue=1e-20,
                bitscore=200.0,
                query_start=1,
                query_end=9,
                subject_start=1,
                subject_end=9,
            ),
        ]

        analyzer._canonicalize_macaca_cross_reactivity_hits("MMAAACCCRR", hits)
        analyzer._canonicalize_macaca_cross_reactivity_hits("MMAAACCCRR", hits)

        self.assertEqual(
            refseq_client.calls,
            [("MISSING", "Macaca fascicularis")],
        )

    def test_retries_macaque_refseq_after_transient_error(self) -> None:
        class FlakyRefSeqClient(FakeRefSeqClient):
            def fetch_canonical_protein(self, *, gene_symbol: str, organism: str):
                if not self.calls:
                    self.calls.append((gene_symbol, organism))
                    raise RuntimeError("temporary RefSeq error")
                return super().fetch_canonical_protein(
                    gene_symbol=gene_symbol,
                    organism=organism,
                )

        refseq_client = FlakyRefSeqClient("MMAAACCCRRRR")
        analyzer = AntigenAnalyzer(
            config=AnalysisConfig(generate_assets=False, prefer_precomputed_references=False),
            refseq_client=refseq_client,
        )
        hits = [
            BlastHit(
                subject_id="tr|TEST|TEST_MACFA",
                description="tr|TEST|TEST_MACFA OS=Macaca fascicularis GN=EGFR",
                species="Macaca fascicularis",
                identity=80.0,
                coverage=90.0,
                alignment_length=9,
                evalue=1e-20,
                bitscore=200.0,
                query_start=1,
                query_end=9,
                subject_start=1,
                subject_end=9,
            ),
        ]

        first = analyzer._canonicalize_macaca_cross_reactivity_hits("MMAAACCCRR", hits)
        second = analyzer._canonicalize_macaca_cross_reactivity_hits("MMAAACCCRR", hits)

        self.assertEqual(first, [])
        self.assertEqual(second[0].subject_id, "XP_TEST_1")
        self.assertEqual(
            refseq_client.calls,
            [
                ("EGFR", "Macaca fascicularis"),
                ("EGFR", "Macaca fascicularis"),
            ],
        )

    def test_prefers_broad_interpro_family_over_target_specific_label(self) -> None:
        analyzer = AntigenAnalyzer(config=AnalysisConfig(generate_assets=False, prefer_precomputed_references=False))
        family = analyzer._select_canonical_family(
            [
                DomainAnnotation(
                    accession="IPR016335",
                    name="Receptor-type tyrosine-protein phosphatase C",
                    type="family",
                    source_database="INTERPRO",
                    start=1,
                    end=1302,
                ),
                DomainAnnotation(
                    accession="IPR050348",
                    name="Protein-Tyrosine Phosphatase",
                    type="family",
                    source_database="INTERPRO",
                    start=153,
                    end=1002,
                ),
            ],
            sequence_length=1302,
            protein_name="Receptor-type tyrosine-protein phosphatase C",
            gene_symbol="PTPRC",
        )
        self.assertIsNotNone(family)
        assert family is not None
        self.assertEqual(family.accession, "IPR050348")

    def test_reports_integrin_partner_requirement(self) -> None:
        target = build_entry("ITGA5_HUMAN", "P08648", "MMAAACCCRRRR", "Homo sapiens", 9606, "ITGA5")
        target["comments"] = [
            {
                "commentType": "SUBUNIT",
                "texts": [
                    {
                        "value": "Heterodimer of an alpha and a beta subunit. Integrin alpha-5 associates with ITGB1."
                    }
                ],
            }
        ]
        mouse = build_entry("ITGA5_MOUSE", "Q9Z0N0", "MMAAACCCRRRK", "Mus musculus", 10090, "ITGA5")
        macfa = build_entry("ITGA5_MACFA", "A00003", "MMAAACCCRRRQ", "Macaca fascicularis", 9541, "ITGA5")
        config = AnalysisConfig(
            min_structured_segment=3,
            min_domain_size=2,
            min_construct_length=2,
            blast_species=("human",),
            generate_assets=False,
            enable_complex_portal=True,
            prefer_precomputed_references=False,
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = AntigenAnalyzer(
                config=config,
                uniprot_client=FakeUniProtClient(target, mouse, macfa),
                alphafold_client=FakeAlphaFoldClient(Path(tmpdir)),
                hgnc_client=FakeHGNCClient(),
                blast_client=FakeBlastClient(),
                refseq_client=FakeRefSeqClient(macfa["sequence"]["value"]),
                interpro_client=FakeInterProClient(),
                complex_portal_client=FakeComplexPortalClient(),
                asset_generator=FakeAssetGenerator(),
            )
            report = analyzer.analyze_target("ITGA5", output_dir=tmpdir)
        self.assertTrue(report.assembly_requirements)
        requirement = next(
            item
            for item in report.assembly_requirements
            if item.classification == "obligatory_partner_requirement"
        )
        self.assertTrue(requirement.obligatory)
        self.assertEqual(requirement.confidence, "high")
        self.assertIn("ITGB1", requirement.partners)
        self.assertNotIn("ITGA5", requirement.partners)
        self.assertIn("complex portal", requirement.summary.lower())
        self.assertEqual(requirement.classification, "obligatory_partner_requirement")
        self.assertIn("CPX-INT1", requirement.complex_portal_ids)
        self.assertIn("CPX-MOUSE-INT1", requirement.complex_portal_ids)
        self.assertIn("mouse support", requirement.source.lower())
        self.assertTrue(report.complex_portal_complexes)

    def test_reports_integrin_partner_requirement_without_complex_portal_record(self) -> None:
        target = build_entry("ITGA1_HUMAN", "P56199", "MMAAACCCRRRR", "Homo sapiens", 9606, "ITGA1")
        target["comments"] = [
            {
                "commentType": "SUBUNIT",
                "texts": [
                    {
                        "value": "Heterodimer of an alpha and a beta subunit. Alpha-1 associates with beta-1."
                    }
                ],
            }
        ]
        mouse = build_entry("ITGA1_MOUSE", "Q9Z0N1", "MMAAACCCRRRK", "Mus musculus", 10090, "ITGA1")
        macfa = build_entry("ITGA1_MACFA", "A00004", "MMAAACCCRRRQ", "Macaca fascicularis", 9541, "ITGA1")
        config = AnalysisConfig(
            min_structured_segment=3,
            min_domain_size=2,
            min_construct_length=2,
            blast_species=("human",),
            generate_assets=False,
            enable_complex_portal=False,
            prefer_precomputed_references=False,
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = AntigenAnalyzer(
                config=config,
                uniprot_client=FakeUniProtClient(target, mouse, macfa),
                alphafold_client=FakeAlphaFoldClient(Path(tmpdir)),
                hgnc_client=FakeHGNCClient(),
                blast_client=FakeBlastClient(),
                refseq_client=FakeRefSeqClient(macfa["sequence"]["value"]),
                interpro_client=FakeInterProClient(),
                complex_portal_client=FakeComplexPortalClient(),
                asset_generator=FakeAssetGenerator(),
            )
            report = analyzer.analyze_target("ITGA1", output_dir=tmpdir)

        requirement = next(
            item
            for item in report.assembly_requirements
            if item.classification == "obligatory_partner_requirement"
        )
        self.assertTrue(requirement.obligatory)
        self.assertEqual(requirement.source, "Family rule + UniProt")
        self.assertIn("ITGB1", requirement.partners)

    def test_complex_portal_accession_partners_map_to_integrin_symbols(self) -> None:
        analyzer = AntigenAnalyzer.__new__(AntigenAnalyzer)
        requirements = analyzer._requirements_from_complex_portal(
            [
                ComplexPortalComplex(
                    complex_ac="CPX-1794",
                    name="Integrin alpha5-beta1 complex",
                    species="Homo sapiens; 9606",
                    predicted_complex=False,
                    confidence_score=5,
                    complex_assemblies=["Heterodimer"],
                    participants=[
                        ComplexPortalParticipant(
                            identifier="P08648",
                            name="P08648",
                            interactor_type="protein",
                        ),
                        ComplexPortalParticipant(
                            identifier="P05556",
                            name="P05556",
                            interactor_type="protein",
                        ),
                    ],
                )
            ],
            gene_symbol="ITGA5",
            target_accession="P08648",
        )

        self.assertEqual(len(requirements), 1)
        self.assertTrue(requirements[0].obligatory)
        self.assertEqual(requirements[0].classification, "obligatory_partner_requirement")
        self.assertEqual(requirements[0].partners, ["ITGB1"])

    def test_complex_portal_fetches_human_and_mouse_targets(self) -> None:
        target = build_entry("ITGA5_HUMAN", "P08648", "MMAAACCCRRRR", "Homo sapiens", 9606, "ITGA5")
        mouse = build_entry("ITGA5_MOUSE", "Q9Z0N0", "MMAAACCCRRRK", "Mus musculus", 10090, "ITGA5")
        macfa = build_entry("ITGA5_MACFA", "A00003", "MMAAACCCRRRQ", "Macaca fascicularis", 9541, "ITGA5")
        complex_portal_client = FakeComplexPortalClient()
        config = AnalysisConfig(
            min_structured_segment=3,
            min_domain_size=2,
            min_construct_length=2,
            blast_species=("human",),
            generate_assets=False,
            enable_complex_portal=True,
            prefer_precomputed_references=False,
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = AntigenAnalyzer(
                config=config,
                uniprot_client=FakeUniProtClient(target, mouse, macfa),
                alphafold_client=FakeAlphaFoldClient(Path(tmpdir)),
                hgnc_client=FakeHGNCClient(),
                blast_client=FakeBlastClient(),
                refseq_client=FakeRefSeqClient(macfa["sequence"]["value"]),
                interpro_client=FakeInterProClient(),
                complex_portal_client=complex_portal_client,
                asset_generator=FakeAssetGenerator(),
            )
            analyzer.analyze_target("ITGA5", output_dir=tmpdir)
        calls = complex_portal_client.calls
        self.assertTrue(any(call["accession"] == "P08648" and call["taxon_id"] == 9606 for call in calls))
        self.assertTrue(any(call["accession"] == "Q9Z0N0" and call["taxon_id"] == 10090 for call in calls))

    def test_nonintegrin_subunit_comments_do_not_create_assembly_requirements(self) -> None:
        target = build_entry("EGFR_HUMAN", "P00533", "MMAAACCCRRRR", "Homo sapiens", 9606, "EGFR")
        target["comments"] = [
            {
                "commentType": "SUBUNIT",
                "texts": [
                    {
                        "value": "Binding of the ligand triggers homo- and/or heterodimerization of the receptor. Heterodimer with ERBB2. Interacts with GRB2."
                    }
                ],
            }
        ]
        mouse = build_entry("EGFR_MOUSE", "Q01279", "MMAAACCCRRRK", "Mus musculus", 10090, "EGFR")
        macfa = build_entry("EGFR_MACFA", "A00003", "MMAAACCCRRRQ", "Macaca fascicularis", 9541, "EGFR")
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = AntigenAnalyzer(
                config=AnalysisConfig(generate_assets=False, min_structured_segment=3, min_domain_size=2, min_construct_length=2, blast_species=("human",), prefer_precomputed_references=False),
                # Complex Portal is disabled by default; not needed for this test.
                uniprot_client=FakeUniProtClient(target, mouse, macfa),
                alphafold_client=FakeAlphaFoldClient(Path(tmpdir)),
                hgnc_client=FakeHGNCClient(),
                blast_client=FakeBlastClient(),
                refseq_client=FakeRefSeqClient(macfa["sequence"]["value"]),
                interpro_client=FakeInterProClient(),
                complex_portal_client=FakeComplexPortalClient(),
                asset_generator=FakeAssetGenerator(),
            )
            report = analyzer.analyze_target("EGFR", output_dir=tmpdir)
        self.assertEqual(report.assembly_requirements, [])

    def test_nonintegrin_requires_text_does_not_create_assembly_requirements(self) -> None:
        target = build_entry("EGFR_HUMAN", "P00533", "MMAAACCCRRRR", "Homo sapiens", 9606, "EGFR")
        target["comments"] = [
            {
                "commentType": "SUBUNIT",
                "texts": [
                    {
                        "value": "Interacts with PGRMC1; the interaction requires PGRMC1 homodimerization (PubMed:26988023)."
                    }
                ],
            }
        ]
        mouse = build_entry("EGFR_MOUSE", "Q01279", "MMAAACCCRRRK", "Mus musculus", 10090, "EGFR")
        macfa = build_entry("EGFR_MACFA", "A00003", "MMAAACCCRRRQ", "Macaca fascicularis", 9541, "EGFR")
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = AntigenAnalyzer(
                config=AnalysisConfig(generate_assets=False, min_structured_segment=3, min_domain_size=2, min_construct_length=2, blast_species=("human",), prefer_precomputed_references=False),
                uniprot_client=FakeUniProtClient(target, mouse, macfa),
                alphafold_client=FakeAlphaFoldClient(Path(tmpdir)),
                hgnc_client=FakeHGNCClient(),
                blast_client=FakeBlastClient(),
                refseq_client=FakeRefSeqClient(macfa["sequence"]["value"]),
                interpro_client=FakeInterProClient(),
                complex_portal_client=FakeComplexPortalClient(),
                asset_generator=FakeAssetGenerator(),
            )
            report = analyzer.analyze_target("EGFR", output_dir=tmpdir)
        self.assertEqual(report.assembly_requirements, [])

    def test_mouse_only_complex_portal_context_does_not_trigger_obligatory_warning(self) -> None:
        analyzer = AntigenAnalyzer.__new__(AntigenAnalyzer)
        requirements = analyzer._requirements_from_complex_portal(
            [
                ComplexPortalComplex(
                    complex_ac="CPX-MOUSE-ONLY",
                    name="Mouse test receptor partner complex",
                    species="Mus musculus; 10090",
                    predicted_complex=False,
                    confidence_score=4,
                    complex_assemblies=["Heterodimer"],
                    participants=[
                        ComplexPortalParticipant(identifier="Q00002", name="TEST", interactor_type="protein"),
                        ComplexPortalParticipant(identifier="Q00003", name="PARTNER", interactor_type="protein"),
                    ],
                )
            ],
            gene_symbol="TEST",
            target_accession="P00001",
        )
        self.assertEqual(len(requirements), 1)
        self.assertFalse(requirements[0].obligatory)
        self.assertEqual(requirements[0].classification, "mouse_supported_complex_context")

    def test_human_complex_portal_membership_without_rule_is_stable_context(self) -> None:
        analyzer = AntigenAnalyzer.__new__(AntigenAnalyzer)
        requirements = analyzer._requirements_from_complex_portal(
            [
                ComplexPortalComplex(
                    complex_ac="CPX-HUMAN-STABLE",
                    name="Test receptor partner complex",
                    species="Homo sapiens; 9606",
                    predicted_complex=False,
                    confidence_score=5,
                    complex_assemblies=["Heterodimer"],
                    participants=[
                        ComplexPortalParticipant(identifier="P00001", name="TEST", interactor_type="protein"),
                        ComplexPortalParticipant(identifier="P00002", name="PARTNER", interactor_type="protein"),
                    ],
                )
            ],
            gene_symbol="TEST",
            target_accession="P00001",
        )
        self.assertEqual(len(requirements), 1)
        self.assertFalse(requirements[0].obligatory)
        self.assertEqual(requirements[0].classification, "stable_complex_context")


if __name__ == "__main__":
    unittest.main()
