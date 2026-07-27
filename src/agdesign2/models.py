from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any


def _normalize(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, list):
        return [_normalize(item) for item in value]
    if isinstance(value, dict):
        return {key: _normalize(item) for key, item in value.items()}
    return value


@dataclass(slots=True)
class Serializable:
    def to_dict(self) -> dict[str, Any]:
        return _normalize(asdict(self))

    def to_json(self, *, indent: int = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent, sort_keys=True)


@dataclass(slots=True)
class AnalysisNote(Serializable):
    severity: str
    message: str
    source: str | None = None


@dataclass(slots=True)
class TargetRecord(Serializable):
    accession: str
    entry_name: str
    gene_symbol: str | None
    protein_name: str
    organism: str
    taxon_id: int | None
    sequence: str
    canonical_isoform_id: str | None = None
    alternative_names: list[str] = field(default_factory=list)
    gene_synonyms: list[str] = field(default_factory=list)


@dataclass(slots=True)
class Feature(Serializable):
    type: str
    start: int
    end: int
    description: str | None = None
    source: str = "UniProt"
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class Region(Serializable):
    start: int
    end: int
    label: str
    source: str
    confidence: float | None = None
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class TopologyAnnotation(Serializable):
    topology_class: str
    extracellular_regions: list[Region] = field(default_factory=list)
    transmembrane_regions: list[Region] = field(default_factory=list)
    intramembrane_regions: list[Region] = field(default_factory=list)
    cytoplasmic_regions: list[Region] = field(default_factory=list)
    major_extracellular_region: Region | None = None
    membrane_core_region: Region | None = None
    cytoplasmic_tail_region: Region | None = None


@dataclass(slots=True)
class ExperimentalConstruct(Serializable):
    pdb_id: str
    method: str | None
    resolution: str | None
    chains: list[Region]


@dataclass(slots=True)
class ConstructSuggestion(Serializable):
    name: str
    start: int
    end: int
    score: float
    rationale: str
    classification: str | None = None
    warnings: list[str] = field(default_factory=list)
    evidence: list[str] = field(default_factory=list)
    pdb_id: str | None = None


@dataclass(slots=True)
class FurinSite(Serializable):
    start: int
    end: int
    motif: str
    suggestions: list[str]


@dataclass(slots=True)
class CysteineFinding(Serializable):
    position: int
    paired_with: int | None
    surface_exposed: bool | None
    warning: str | None = None
    closest_unpaired_cysteine: int | None = None
    closest_unpaired_cysteine_distance: float | None = None
    closest_unpaired_cysteine_distance_atom: str | None = None


@dataclass(slots=True)
class ResidueAnnotation(Serializable):
    position: int
    amino_acid: str
    topology_location: str
    surface_exposed: bool | None = None
    surface_accessibility: str | None = None
    surface_accessibility_reason: str | None = None
    full_model_sasa: float | None = None
    full_model_relative_sasa: float | None = None
    domain_context_sasa: float | None = None
    domain_context_relative_sasa: float | None = None
    pae_block_id: str | None = None
    plddt: float | None = None
    source: str = "topology+AlphaFold"


@dataclass(slots=True)
class AssemblyRequirement(Serializable):
    summary: str
    obligatory: bool
    confidence: str
    source: str
    classification: str | None = None
    partners: list[str] = field(default_factory=list)
    rationale: list[str] = field(default_factory=list)
    complex_portal_ids: list[str] = field(default_factory=list)


@dataclass(slots=True)
class MembraneEngineeringSuggestion(Serializable):
    category: str
    title: str
    summary: str
    rationale: list[str] = field(default_factory=list)
    suggested_actions: list[str] = field(default_factory=list)
    evidence: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    start: int | None = None
    end: int | None = None


@dataclass(slots=True)
class GPCREngineeringCassette(Serializable):
    cassette_id: str
    name: str
    sequence: str | None
    source_note: str
    version: str
    checksum_sha256: str | None = None
    default_n_linker: str = ""
    default_c_linker: str = ""
    allowed_use_modes: list[str] = field(default_factory=list)


@dataclass(slots=True)
class GPCREngineeringVariant(Serializable):
    name: str
    category: str
    strategy: str
    summary: str
    loop_source: str | None = None
    sequence: str | None = None
    sequence_role: str = "engineered_receptor"
    cassette_id: str | None = None
    cassette_name: str | None = None
    cassette_sequence: str | None = None
    start: int | None = None
    end: int | None = None
    replaced_start: int | None = None
    replaced_end: int | None = None
    insertion_after: int | None = None
    retained_n_terminal_loop_residues: int | None = None
    retained_c_terminal_loop_residues: int | None = None
    n_linker: str | None = None
    c_linker: str | None = None
    native_length: int | None = None
    engineered_length: int | None = None
    gpcr_segments: list[GPCRdbSegment] = field(default_factory=list)
    gpcr_generic_range: str | None = None
    tag_suggestions: dict[str, Any] = field(default_factory=dict)
    suggested_actions: list[str] = field(default_factory=list)
    evidence: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)


@dataclass(slots=True)
class GPCRdbResidue(Serializable):
    sequence_number: int
    amino_acid: str
    protein_segment: str | None = None
    display_generic_number: str | None = None


@dataclass(slots=True)
class GPCRdbSegment(Serializable):
    name: str
    start: int
    end: int
    generic_start: str | None = None
    generic_end: str | None = None


@dataclass(slots=True)
class GPCRdbMotif(Serializable):
    name: str
    positions: list[int] = field(default_factory=list)
    generic_numbers: list[str] = field(default_factory=list)
    sequence: str | None = None
    status: str = "detected"


@dataclass(slots=True)
class GPCRdbAnnotation(Serializable):
    entry_name: str
    accession: str | None
    name: str | None
    family_slug: str | None = None
    family_name: str | None = None
    family_path: list[str] = field(default_factory=list)
    gpcr_class: str | None = None
    species: str | None = None
    source: str | None = None
    residue_numbering_scheme: str | None = None
    url: str | None = None
    confidence: str = "high"
    detection_source: str = "GPCRdb"
    segments: list[GPCRdbSegment] = field(default_factory=list)
    residues: list[GPCRdbResidue] = field(default_factory=list)
    conserved_motifs: list[GPCRdbMotif] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)


@dataclass(slots=True)
class ComplexPortalParticipant(Serializable):
    identifier: str
    name: str
    description: str | None = None
    stoichiometry: str | None = None
    bio_role: str | None = None
    interactor_type: str | None = None
    binding_regions: list[str] = field(default_factory=list)


@dataclass(slots=True)
class ComplexPortalComplex(Serializable):
    complex_ac: str
    name: str
    species: str
    predicted_complex: bool
    evidence_code: str | None = None
    evidence_description: str | None = None
    confidence_score: int | None = None
    complex_assemblies: list[str] = field(default_factory=list)
    properties: list[str] = field(default_factory=list)
    participants: list[ComplexPortalParticipant] = field(default_factory=list)


@dataclass(slots=True)
class HomologyRecord(Serializable):
    species: str
    accession: str | None
    entry_name: str | None
    ectodomain_start: int | None
    ectodomain_end: int | None
    identity: float | None = None
    coverage: float | None = None
    extracellular_surface_identity: float | None = None
    extracellular_surface_aligned_positions: int | None = None
    available: bool = True
    notes: list[str] = field(default_factory=list)


@dataclass(slots=True)
class FamilyMember(Serializable):
    gene_symbol: str
    gene_name: str
    accession: str | None
    entry_name: str | None
    sequence_length: int | None = None
    ectodomain_start: int | None = None
    ectodomain_end: int | None = None
    ectodomain_length: int | None = None
    notes: list[str] = field(default_factory=list)


@dataclass(slots=True)
class FamilyContext(Serializable):
    gene_symbol: str
    source: str
    family_names: list[str] = field(default_factory=list)
    members: list[FamilyMember] = field(default_factory=list)
    identity_matrix_labels: list[str] = field(default_factory=list)
    identity_matrix: list[list[float | None]] = field(default_factory=list)
    coverage_matrix_labels: list[str] = field(default_factory=list)
    coverage_matrix: list[list[float | None]] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class BlastHit(Serializable):
    subject_id: str
    description: str
    species: str | None
    identity: float
    coverage: float
    alignment_length: int
    evalue: float
    bitscore: float
    query_start: int
    query_end: int
    subject_start: int
    subject_end: int
    query_alignment: str | None = None
    subject_alignment: str | None = None
    alignment_source: str | None = None
    query_alignment_annotation: str | None = None
    extracellular_surface_identity: float | None = None
    extracellular_surface_aligned_positions: int | None = None
    extracellular_surface_matches: int | None = None


@dataclass(slots=True)
class DomainAnnotation(Serializable):
    accession: str
    name: str
    type: str
    source_database: str
    start: int
    end: int
    representative: bool = False
    integrated_accession: str | None = None
    integrated_name: str | None = None
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass(slots=True)
class CanonicalFamily(Serializable):
    accession: str
    name: str
    source_database: str
    fragment_count: int
    start: int | None = None
    end: int | None = None
    covered_length: int | None = None
    coverage_fraction: float | None = None
    notes: list[str] = field(default_factory=list)


@dataclass(slots=True)
class ConstructHomolog(Serializable):
    species: str
    accession: str | None
    entry_name: str | None
    available: bool
    start: int | None = None
    end: int | None = None
    sequence: str | None = None
    identity_to_human: float | None = None
    extracellular_surface_identity_to_human: float | None = None
    extracellular_surface_aligned_positions: int | None = None
    notes: list[str] = field(default_factory=list)


@dataclass(slots=True)
class ConstructStructuralMetrics(Serializable):
    mean_plddt: float | None = None
    min_plddt: float | None = None
    max_plddt: float | None = None
    structured_fraction: float | None = None
    mean_intra_pae: float | None = None
    max_intra_pae: float | None = None


@dataclass(slots=True)
class ConstructSplitDiagnostic(Serializable):
    parent_start: int
    parent_end: int
    boundary_after: int
    score: float
    inter_block_pae: float
    intra_block_pae: float
    separation: float
    boundary_pae: float
    linker_plddt: float | None = None


@dataclass(slots=True)
class ConstructDetail(Serializable):
    name: str
    start: int
    end: int
    length: int
    sequence: str
    score: float
    rationale: str
    classification: str | None = None
    warnings: list[str] = field(default_factory=list)
    evidence: list[str] = field(default_factory=list)
    pdb_id: str | None = None
    human_entry_name: str | None = None
    domain_summary: str | None = None
    complete_domain_annotations: list[DomainAnnotation] = field(default_factory=list)
    clipped_domain_annotations: list[DomainAnnotation] = field(default_factory=list)
    homologs: list[ConstructHomolog] = field(default_factory=list)
    cysteine_analysis: list[CysteineFinding] = field(default_factory=list)
    ptms: list[Feature] = field(default_factory=list)
    ligand_interactions: list[Feature] = field(default_factory=list)
    gpcr_segments: list[GPCRdbSegment] = field(default_factory=list)
    gpcr_generic_range: str | None = None
    structural_metrics: ConstructStructuralMetrics | None = None
    split_diagnostics: list[ConstructSplitDiagnostic] = field(default_factory=list)
    structure_image: str | None = None
    quality_plot: str | None = None


@dataclass(slots=True)
class AnalysisReport(Serializable):
    target: TargetRecord
    resolution: dict[str, Any]
    features: list[Feature]
    ectodomain: Region | None
    structural_regions: list[Region]
    experimental_constructs: list[ExperimentalConstruct]
    construct_recommendations: list[ConstructSuggestion]
    furin_sites: list[FurinSite]
    cysteine_analysis: list[CysteineFinding]
    assembly_requirements: list[AssemblyRequirement]
    interpro_annotations: list[DomainAnnotation]
    species_name_matches: list[HomologyRecord]
    ectodomain_homology: list[HomologyRecord]
    family_context: FamilyContext | None
    cross_reactivity_hits: list[BlastHit]
    complex_portal_complexes: list[ComplexPortalComplex] = field(default_factory=list)
    complex_portal_lookup: dict[str, dict[str, str]] = field(default_factory=dict)
    full_length_family_context: FamilyContext | None = None
    full_length_cross_reactivity_hits: list[BlastHit] = field(default_factory=list)
    ligand_interactions: list[Feature] = field(default_factory=list)
    construct_details: list[ConstructDetail] = field(default_factory=list)
    advanced_membrane_suggestions: list[MembraneEngineeringSuggestion] = field(default_factory=list)
    gpcr_engineering_variants: list[GPCREngineeringVariant] = field(default_factory=list)
    notes: list[AnalysisNote] = field(default_factory=list)
    canonical_family: CanonicalFamily | None = None
    topology: TopologyAnnotation | None = None
    ptms: list[Feature] = field(default_factory=list)
    residue_annotations: list[ResidueAnnotation] = field(default_factory=list)
    gpcr_annotation: GPCRdbAnnotation | None = None

    def write_json(self, path: str | Path) -> None:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(self.to_json() + "\n", encoding="utf-8")

    def to_markdown(self) -> str:
        from .reporting import render_markdown_report

        return render_markdown_report(self)

    def write_markdown(self, path: str | Path) -> None:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(self.to_markdown(), encoding="utf-8")
