from __future__ import annotations

import csv
from dataclasses import dataclass, field, replace
from difflib import SequenceMatcher
from pathlib import Path
import re
from typing import Any

from .alphafold import AlphaFoldClient
from .af3 import is_af3_artifact
from .blast import BlastClient, extract_accession, extract_entry_name
from .cache import FileCache
from .complex_portal import ComplexPortalClient
from .config import AnalysisConfig
from .gpcr_engineering import load_gpcr_engineering_registry
from .gpcrdb import GPCRdbClient
from .hgnc import HGNCClient
from .http import HttpClient
from .interpro import InterProClient
from .models import (
    AnalysisNote,
    AnalysisReport,
    AssemblyRequirement,
    CanonicalFamily,
    ComplexPortalComplex,
    ConstructDetail,
    ConstructHomolog,
    ConstructSuggestion,
    CysteineFinding,
    DomainAnnotation,
    ExperimentalConstruct,
    FamilyContext,
    FamilyMember,
    Feature,
    FurinSite,
    GPCRdbAnnotation,
    GPCREngineeringVariant,
    GPCRdbSegment,
    HomologyRecord,
    MembraneEngineeringSuggestion,
    Region,
    ResidueAnnotation,
    TopologyAnnotation,
)
from .precomputed import PrecomputedReferenceStore
from .refseq import RefSeqClient
from .sequence_utils import find_furin_sites, global_align, identity_for_query_positions, map_query_region_to_subject, slice_sequence
from .structure_utils import (
    analyze_cysteines,
    build_pae_stats,
    classify_surface_accessibility,
    collect_region_split_diagnostics,
    find_structured_regions,
    load_pae_matrix,
    parse_alphafold_pdb,
    split_domains_from_pae,
    summarize_construct_quality,
)
from .topology import TopologyResult, apply_surfy_topology_fallback, derive_ectodomain
from .uniprot import UniProtClient
from .visualization import ConstructAssetGenerator


_PTM_FEATURE_TYPES = {
    "CARBOHYD",
    "MOD_RES",
    "LIPID",
    "DISULFID",
    "CROSSLNK",
    "INIT_MET",
    "SIGNAL",
    "PROPEP",
    "PEPTIDE",
    "CHAIN",
}
_PROCESSING_FEATURE_TYPES = {"SIGNAL", "PROPEP", "PEPTIDE"}

_CANONICAL_INTEGRIN_SYMBOL_BY_ACCESSION = {
    "P56199": "ITGA1",
    "O75578": "ITGA10",
    "Q9UKX5": "ITGA11",
    "P17301": "ITGA2",
    "P08514": "ITGA2B",
    "P26006": "ITGA3",
    "P13612": "ITGA4",
    "P08648": "ITGA5",
    "P23229": "ITGA6",
    "Q13683": "ITGA7",
    "P53708": "ITGA8",
    "Q13797": "ITGA9",
    "Q13349": "ITGAD",
    "P38570": "ITGAE",
    "P20701": "ITGAL",
    "P11215": "ITGAM",
    "P06756": "ITGAV",
    "P20702": "ITGAX",
    "P05556": "ITGB1",
    "P05107": "ITGB2",
    "P05106": "ITGB3",
    "P16144": "ITGB4",
    "P18084": "ITGB5",
    "P18564": "ITGB6",
    "P26010": "ITGB7",
    "P26012": "ITGB8",
}
_CANONICAL_INTEGRIN_SYMBOLS = frozenset(_CANONICAL_INTEGRIN_SYMBOL_BY_ACCESSION.values())
_CANONICAL_INTEGRIN_ALPHA_SYMBOLS = frozenset(
    symbol for symbol in _CANONICAL_INTEGRIN_SYMBOLS if symbol.startswith("ITGA")
)
_CANONICAL_INTEGRIN_BETA_SYMBOLS = frozenset(
    symbol for symbol in _CANONICAL_INTEGRIN_SYMBOLS if symbol.startswith("ITGB")
)


@dataclass(slots=True)
class _TerminalTrimDecision:
    enabled: bool = False
    start: int | None = None
    end: int | None = None
    evidence: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    blocked_reason: str | None = None


class AntigenAnalyzer:
    def __init__(
        self,
        config: AnalysisConfig | None = None,
        *,
        uniprot_client: UniProtClient | None = None,
        alphafold_client: AlphaFoldClient | None = None,
        hgnc_client: HGNCClient | None = None,
        complex_portal_client: ComplexPortalClient | None = None,
        blast_client: BlastClient | None = None,
        refseq_client: RefSeqClient | None = None,
        interpro_client: InterProClient | None = None,
        gpcrdb_client: GPCRdbClient | None = None,
        asset_generator: ConstructAssetGenerator | None = None,
        precomputed_store: PrecomputedReferenceStore | None = None,
    ) -> None:
        self.config = config or AnalysisConfig()
        self.config.ensure_directories()
        cache = FileCache(self.config.cache_dir)
        http = HttpClient(cache=cache)
        complex_portal_http = HttpClient(cache=cache, timeout=self.config.complex_portal_timeout)
        self.uniprot_client = uniprot_client or UniProtClient(http)
        self.alphafold_client = alphafold_client or AlphaFoldClient(http, self.config)
        self.hgnc_client = hgnc_client or HGNCClient(http, self.uniprot_client, self.config)
        self.complex_portal_client = complex_portal_client or ComplexPortalClient(complex_portal_http)
        self.blast_client = blast_client or BlastClient(http, self.config)
        self.refseq_client = refseq_client or RefSeqClient(http)
        self._macaca_refseq_by_gene: dict[str, Any] = {}
        self.interpro_client = interpro_client or InterProClient(http)
        gpcrdb_cache = FileCache(self.config.gpcrdb_cache_dir) if self.config.gpcrdb_cache_dir else cache
        self.gpcrdb_client = gpcrdb_client or GPCRdbClient(HttpClient(cache=gpcrdb_cache))
        self.asset_generator = asset_generator or ConstructAssetGenerator()
        self.precomputed_store = precomputed_store or PrecomputedReferenceStore(self.config)
        self._surfy_topology_index: dict[str, str] | None = None
        self._secreted_universe_index: set[str] | None = None

    def analyze_target(self, query: str, output_dir: str | Path | None = None) -> AnalysisReport:
        self._progress(f"Starting analysis for {query}")
        self._progress("Resolving target in UniProt")
        resolution = self.uniprot_client.resolve_target(
            query,
            target_species_taxon=self.config.species_taxonomy[self.config.target_species],
        )
        target = resolution.target
        self._progress(f"Resolved target to {target.entry_name} ({target.accession})")
        self._progress("Fetching UniProt features and experimental constructs")
        features = self.uniprot_client.get_features(resolution.entry)
        experimental_constructs = self.uniprot_client.get_experimental_constructs(resolution.entry)
        self._progress("Deriving ectodomain / extracellular design region")
        topology = derive_ectodomain(features, len(target.sequence), self.config)
        secreted_fallback_note = self._apply_secreted_universe_fallback(
            target=target,
            topology=topology,
            features=features,
        )
        if secreted_fallback_note is not None:
            topology.notes.append(secreted_fallback_note)
        if self._apply_surfy_topology_fallback(target=target, topology=topology.topology):
            topology.notes.append(
                AnalysisNote(
                    severity="info",
                    message=(
                        "UniProt topological-domain annotations did not define extracellular/cytoplasmic loops; "
                        "using SURFY topology segments for residue-level extracellular accessibility annotation."
                    ),
                    source="SURFY",
                )
            )
        notes = list(resolution.notes) + topology.notes
        is_multipass = bool(topology.topology and topology.topology.topology_class.startswith("multipass"))
        is_multipass_mixed = (
            bool(topology.ectodomain)
            and topology.topology is not None
            and topology.topology.topology_class == "multipass_mixed"
            and (topology.ectodomain.end - topology.ectodomain.start + 1)
            >= self.config.multipass_mixed_ectodomain_min_length
        )
        experimental_constructs = self._filter_experimental_constructs(
            experimental_constructs,
            topology.ectodomain,
            topology.topology,
        )
        interpro_annotations: list[DomainAnnotation] = []
        canonical_family: CanonicalFamily | None = None
        precomputed_ortholog = self.precomputed_store.find_ortholog_record(target=target)
        precomputed_canonical_family = self.precomputed_store.canonical_family_from_record(precomputed_ortholog)
        try:
            self._progress("Fetching InterPro / Pfam annotations")
            interpro_annotations = self.interpro_client.fetch_annotations(target.accession)
        except Exception as exc:
            notes.append(
                AnalysisNote(
                    severity="warning",
                    message=f"InterPro/Pfam annotation lookup failed: {exc}",
                    source="InterPro",
                )
            )
        live_canonical_family = self._select_canonical_family(
            interpro_annotations,
            sequence_length=len(target.sequence),
            protein_name=target.protein_name,
            gene_symbol=target.gene_symbol,
        )
        if precomputed_ortholog is not None:
            canonical_family = self._merge_canonical_family_sources(
                precomputed=precomputed_canonical_family,
                live=live_canonical_family if precomputed_canonical_family is not None else None,
            )
        else:
            canonical_family = live_canonical_family
        interpro_annotations = self._filter_design_annotations(interpro_annotations, topology.ectodomain)
        gpcr_annotation = self._fetch_gpcr_annotation(
            target=target,
            topology=topology.topology,
            interpro_annotations=interpro_annotations,
            notes=notes,
        )

        structural_regions: list[Region] = []
        cysteine_analysis: list[CysteineFinding] = []
        residue_annotations: list[ResidueAnnotation] = []
        surface_exposure: dict[int, bool | None] = {}
        surface_accessibility = {}
        structure: dict[str, Any] | None = None
        pae_matrix: list[list[float]] | None = None
        pae_stats = None
        alphafold_pdb_path: Path | None = None
        if topology.ectodomain or is_multipass:
            try:
                self._progress("Fetching AlphaFold structure and PAE")
                pdb_path, pae_path = self.alphafold_client.ensure_artifacts(
                    target.accession,
                    canonical_sequence=target.sequence,
                    canonical_isoform_id=target.canonical_isoform_id,
                )
                alphafold_pdb_path = pdb_path
                if is_af3_artifact(pdb_path):
                    notes.append(
                        AnalysisNote(
                            severity="info",
                            message=(
                                "Using a local AlphaFold 3 model because no matching canonical "
                                "AlphaFold DB model was available."
                            ),
                            source="AlphaFold Server",
                        )
                    )
                self._progress("Parsing AlphaFold structure")
                structure = parse_alphafold_pdb(pdb_path)
                pae_matrix = load_pae_matrix(pae_path)
                pae_stats = build_pae_stats(pae_matrix)
                if topology.ectodomain:
                    self._progress("Calling structured regions and domains")
                    lenient_structured_regions = find_structured_regions(
                        topology.ectodomain.start,
                        topology.ectodomain.end,
                        structure["plddt"],
                        self.config,
                        plddt_threshold=self.config.lenient_plddt_structured_threshold,
                        max_disordered_gap=self.config.lenient_max_disordered_gap,
                        label="Lenient structured region",
                        metadata={"stringency": "lenient"},
                    )
                    strict_seed_segments = find_structured_regions(
                        topology.ectodomain.start,
                        topology.ectodomain.end,
                        structure["plddt"],
                        self.config,
                        label="Strict seed region",
                        metadata={"stringency": "strict_seed"},
                    )
                    structural_domains = []
                    for segment in strict_seed_segments:
                        structural_domains.extend(
                            split_domains_from_pae(
                                segment,
                                pae_stats,
                                self.config,
                                residue_plddt=structure["plddt"],
                            )
                        )
                    strict_domains = [
                        replace(
                            domain,
                            label="Strict predicted domain",
                            metadata={**domain.metadata, "stringency": "strict"},
                        )
                        for domain in structural_domains
                    ]
                    structural_regions = self._annotate_structural_regions(
                        self._dedupe_regions(lenient_structured_regions + strict_domains),
                        interpro_annotations=interpro_annotations,
                        residue_plddt=structure["plddt"],
                        pae_stats=pae_stats,
                    )
                    self._progress("Classifying residue solvent accessibility")
                    surface_target_residues = self._surface_accessibility_target_positions(
                        sequence_length=len(target.sequence),
                        topology=topology.topology,
                        residue_names=structure["residue_names"],
                        scope=topology.ectodomain,
                    )
                    surface_accessibility = classify_surface_accessibility(
                        structure["atoms"],
                        structure["residue_names"],
                        structure["plddt"],
                        self._surface_context_regions(
                            structural_regions=structural_regions,
                            topology=topology.topology,
                            sequence_length=len(target.sequence),
                        ),
                        self.config,
                        target_residues=surface_target_residues
                        if self.config.sasa_target_accessibility_only
                        else None,
                    )
                    surface_exposure = {
                        residue: result.surface_exposed
                        for residue, result in surface_accessibility.items()
                    }
                    residue_annotations = self._build_residue_annotations(
                        sequence=target.sequence,
                        topology=topology.topology,
                        residue_plddt=structure["plddt"],
                        surface_accessibility=surface_accessibility,
                    )
                    cysteine_analysis = analyze_cysteines(
                        structure["atoms"],
                        structure["residue_names"],
                        topology.ectodomain.start,
                        topology.ectodomain.end,
                        surface_exposure=surface_exposure,
                    )
                    self._progress(
                        f"Derived {len(structural_regions)} structural regions and {len(cysteine_analysis)} cysteine findings"
                    )
                elif is_multipass:
                    self._progress("Classifying full-length residue solvent accessibility")
                    surface_target_residues = self._surface_accessibility_target_positions(
                        sequence_length=len(target.sequence),
                        topology=topology.topology,
                        residue_names=structure["residue_names"],
                        scope=Region(start=1, end=len(target.sequence), label="Full length", source="topology"),
                    )
                    surface_accessibility = classify_surface_accessibility(
                        structure["atoms"],
                        structure["residue_names"],
                        structure["plddt"],
                        self._surface_context_regions(
                            structural_regions=structural_regions,
                            topology=topology.topology,
                            sequence_length=len(target.sequence),
                        ),
                        self.config,
                        target_residues=surface_target_residues
                        if self.config.sasa_target_accessibility_only
                        else None,
                    )
                    surface_exposure = {
                        residue: result.surface_exposed
                        for residue, result in surface_accessibility.items()
                    }
                    residue_annotations = self._build_residue_annotations(
                        sequence=target.sequence,
                        topology=topology.topology,
                        residue_plddt=structure["plddt"],
                        surface_accessibility=surface_accessibility,
                    )
                    cysteine_analysis = analyze_cysteines(
                        structure["atoms"],
                        structure["residue_names"],
                        1,
                        len(target.sequence),
                        surface_exposure=surface_exposure,
                    )
                    notes.append(
                        AnalysisNote(
                            severity="info",
                            message=(
                                "AlphaFold artifacts were loaded for the full multipass protein; soluble ECD structural "
                                "region calling was skipped because no extracellular region met the mixed-case length threshold."
                            ),
                            source="AlphaFold",
                        )
                    )
            except Exception as exc:
                notes.append(
                    AnalysisNote(
                        severity="warning",
                        message=f"Structure-driven analyses were skipped: {exc}",
                        source="AlphaFold",
                    )
                )
        else:
            notes.append(
                AnalysisNote(
                    severity="warning",
                    message="Skipping structure-driven analyses because the ectodomain could not be determined.",
                    source="pipeline",
                )
                )

        ligand_interactions = self._collect_ligand_interactions(features, topology.ectodomain)
        ptms = self._collect_ptm_annotations(features)
        furin_sites = self._find_furin_sites(target.sequence, topology.ectodomain)

        self._progress("Building construct recommendations")
        constructs = self._suggest_constructs(
            features=features,
            topology=topology.topology,
            ectodomain=topology.ectodomain,
            sequence_length=len(target.sequence),
            structural_regions=structural_regions,
            experimental_constructs=experimental_constructs,
            cysteine_analysis=cysteine_analysis,
            interpro_annotations=interpro_annotations,
            residue_plddt=structure["plddt"] if structure is not None else None,
            ligand_interactions=ligand_interactions,
            ptms=ptms,
            gpcr_annotation=gpcr_annotation,
        )
        if precomputed_ortholog is not None:
            self._progress("Resolving species homologs from precomputed ortholog table")
            species_matches, homology, homolog_context = self.precomputed_store.build_species_matches(
                target=target,
                ectodomain=topology.ectodomain,
                record=precomputed_ortholog,
            )
            notes.append(
                AnalysisNote(
                    severity="info",
                    message="Species homologs loaded from precomputed ortholog reference table.",
                    source="precomputed",
                )
            )
        else:
            try:
                self._progress("Resolving species homologs")
                species_matches, homology, homolog_context = self._collect_species_matches(
                    target.entry_name, target.gene_symbol or target.entry_name.split("_", 1)[0], topology.ectodomain, target.sequence
                )
            except Exception as exc:
                notes.append(
                    AnalysisNote(
                        severity="warning",
                        message=f"Species homolog lookup failed: {exc}",
                        source="UniProt",
                    )
                )
                species_matches, homology = [], []
                homolog_context = {}
        homology = self._annotate_homology_with_surface_identity(
            homology,
            homolog_context=homolog_context,
            ectodomain=topology.ectodomain,
            extracellular_surface_positions=self._extracellular_surface_positions(residue_annotations),
        )
        family_context = self.precomputed_store.load_paralog_context(
            target=target,
            max_members=self.config.max_paralog_context_members,
        )
        if family_context is not None:
            family_context = self._annotate_family_context_with_surface_identity(
                family_context=family_context,
                target_entry_name=target.entry_name,
                ectodomain=topology.ectodomain,
                residue_annotations=residue_annotations,
            )
            notes.append(
                AnalysisNote(
                    severity="info",
                    message="Family context loaded from the precomputed paralog reference set.",
                    source="precomputed",
                )
            )
        else:
            family_context = self.precomputed_store.load_family_context(
                target=target,
                canonical_family=canonical_family,
                max_members=self.config.max_family_context_members,
            )
            if family_context is not None and family_context.source == "Precomputed InterPro family alignments":
                family_context = self._annotate_family_context_with_surface_identity(
                    family_context=family_context,
                    target_entry_name=target.entry_name,
                    ectodomain=topology.ectodomain,
                    residue_annotations=residue_annotations,
                )
                notes.append(
                    AnalysisNote(
                        severity="info",
                        message="Family context loaded from precomputed family alignments.",
                        source="precomputed",
                    )
                )
            elif self._should_skip_live_family_context_fallback(
                precomputed_ortholog=precomputed_ortholog,
                canonical_family=canonical_family,
            ):
                notes.append(
                    AnalysisNote(
                        severity="info",
                        message="Skipping live family-context fallback because this target is covered by precomputed reference tables.",
                        source="precomputed",
                    )
                )
            else:
                try:
                    self._progress("Fetching family context")
                    family_context = self.hgnc_client.fetch_family_context(target.gene_symbol)
                except Exception as exc:
                    notes.append(
                        AnalysisNote(
                            severity="warning",
                            message=f"Family context lookup failed: {exc}",
                            source="HGNC",
                        )
                    )
                    family_context = None
        complex_portal_lookup: dict[str, dict[str, str]] = {}
        complex_portal_complexes = self._fetch_complex_portal_context(
            target=target,
            homolog_context=homolog_context,
            notes=notes,
            lookup_status=complex_portal_lookup,
        )
        cross_reactivity_hits = []
        full_length_cross_reactivity_hits = []
        assembly_requirements = self._collect_assembly_requirements(
            entry=resolution.entry,
            gene_symbol=target.gene_symbol,
            family_context=family_context,
            complex_portal_complexes=complex_portal_complexes,
            target_accession=target.accession,
        )
        advanced_membrane_suggestions = self._suggest_advanced_membrane_engineering(
            target=target,
            features=features,
            topology=topology.topology,
            interpro_annotations=interpro_annotations,
            gpcr_annotation=gpcr_annotation,
        )
        gpcr_engineering_variants = self._suggest_gpcr_engineering_variants(
            target=target,
            topology=topology.topology,
            features=features,
            ligand_interactions=ligand_interactions,
            ptms=ptms,
            gpcr_annotation=gpcr_annotation,
        )
        if topology.ectodomain and self.config.run_cross_reactivity:
            ectodomain_sequence = slice_sequence(target.sequence, topology.ectodomain.start, topology.ectodomain.end)
            try:
                self._progress("Running cross-reactivity BLAST")
                cross_reactivity_hits = self.blast_client.search(
                    ectodomain_sequence,
                    target_accession=target.accession,
                    target_entry_name=target.entry_name,
                )
                cross_reactivity_hits = self._canonicalize_macaca_cross_reactivity_hits(
                    ectodomain_sequence,
                    cross_reactivity_hits,
                )
                cross_reactivity_hits = self._annotate_blast_hits_with_surface_identity(
                    cross_reactivity_hits,
                    query_region_start=topology.ectodomain.start,
                    extracellular_surface_positions=self._extracellular_surface_positions(residue_annotations),
                )
            except Exception as exc:
                notes.append(
                    AnalysisNote(
                        severity="warning",
                        message=f"Cross-reactivity BLAST search failed: {exc}",
                    source="BLAST",
                )
            )
        elif topology.ectodomain:
            self._progress("Skipping cross-reactivity BLAST (disabled)")
            notes.append(
                AnalysisNote(
                    severity="info",
                    message="Design-region BLAST sequence similarity search was not run.",
                    source="BLAST",
                )
            )

        self._progress("Assembling construct detail records")
        construct_details = self._build_construct_details(
            target_entry_name=target.entry_name,
            target_sequence=target.sequence,
            constructs=constructs,
            homolog_context=homolog_context,
            ectodomain=topology.ectodomain,
            ligand_interactions=ligand_interactions,
            ptms=ptms,
            furin_sites=furin_sites,
            interpro_annotations=interpro_annotations,
            structure_atoms=structure["atoms"] if structure is not None else None,
            residue_names=structure["residue_names"] if structure is not None else None,
            residue_plddt=structure["plddt"] if structure is not None else None,
            surface_exposure=surface_exposure,
            residue_annotations=residue_annotations,
            pae_matrix=pae_matrix,
            pae_stats=pae_stats,
            gpcr_annotation=gpcr_annotation,
        )
        if (
            output_dir
            and construct_details
            and self.config.generate_assets
            and topology.ectodomain
            and structure is not None
            and pae_matrix is not None
            and alphafold_pdb_path is not None
        ):
            output_root = Path(output_dir)
            stem = f"{target.entry_name.lower()}_report"
            assets_dir = output_root / f"{stem}_assets"
            try:
                self._progress(f"Generating construct assets in {assets_dir}")
                self._generate_construct_assets(
                    construct_details=construct_details,
                    assets_dir=assets_dir,
                    alphafold_pdb_path=alphafold_pdb_path,
                    residue_plddt=structure["plddt"],
                    pae_matrix=pae_matrix,
                    ectodomain=topology.ectodomain,
                    residue_annotations=residue_annotations,
                )
            except Exception as exc:
                notes.append(
                    AnalysisNote(
                        severity="warning",
                        message=f"Construct asset generation failed: {exc}",
                        source="rendering",
                    )
                )
        run_full_length_cross_reactivity = is_multipass or topology.ectodomain is None
        if run_full_length_cross_reactivity and self.config.run_cross_reactivity:
            try:
                self._progress("Running full-length cross-reactivity BLAST")
                full_length_cross_reactivity_hits = self.blast_client.search(
                    target.sequence,
                    target_accession=target.accession,
                    target_entry_name=target.entry_name,
                )
                full_length_cross_reactivity_hits = self._canonicalize_macaca_cross_reactivity_hits(
                    target.sequence,
                    full_length_cross_reactivity_hits,
                )
                full_length_cross_reactivity_hits = self._annotate_blast_hits_with_surface_identity(
                    full_length_cross_reactivity_hits,
                    query_region_start=1,
                    extracellular_surface_positions=self._extracellular_surface_positions(residue_annotations),
                )
            except Exception as exc:
                notes.append(
                    AnalysisNote(
                        severity="warning",
                        message=f"Full-length cross-reactivity BLAST search failed: {exc}",
                    source="BLAST",
                )
            )
        elif run_full_length_cross_reactivity:
            self._progress("Skipping full-length cross-reactivity BLAST (disabled)")
            notes.append(
                AnalysisNote(
                    severity="info",
                    message="Full-length BLAST sequence similarity search was not run.",
                    source="BLAST",
                )
            )
        full_length_family_context = None
        if is_multipass_mixed or (is_multipass and topology.ectodomain is None):
            full_length_family_context = self._build_full_length_family_context(
                target=target,
                family_context=family_context,
                canonical_family=canonical_family,
            )

        report = AnalysisReport(
            target=target,
            resolution={
                "query": query,
                "query_type": resolution.query_type,
                "resolved_accession": target.accession,
                "resolved_entry_name": target.entry_name,
            },
            features=features,
            ectodomain=topology.ectodomain,
            structural_regions=structural_regions,
            experimental_constructs=experimental_constructs,
            construct_recommendations=constructs,
            furin_sites=furin_sites,
            cysteine_analysis=cysteine_analysis,
            residue_annotations=residue_annotations,
            ptms=ptms,
            assembly_requirements=assembly_requirements,
            complex_portal_complexes=complex_portal_complexes,
            complex_portal_lookup=complex_portal_lookup,
            interpro_annotations=interpro_annotations,
            ligand_interactions=ligand_interactions,
            species_name_matches=species_matches,
            ectodomain_homology=homology,
            family_context=family_context,
            cross_reactivity_hits=cross_reactivity_hits,
            full_length_family_context=full_length_family_context,
            full_length_cross_reactivity_hits=full_length_cross_reactivity_hits,
            construct_details=construct_details,
            advanced_membrane_suggestions=advanced_membrane_suggestions,
            gpcr_engineering_variants=gpcr_engineering_variants,
            notes=notes,
            canonical_family=canonical_family,
            topology=topology.topology,
            gpcr_annotation=gpcr_annotation,
        )
        if output_dir:
            output_root = Path(output_dir)
            stem = f"{target.entry_name.lower()}_report"
            self._progress("Writing JSON and Markdown reports")
            report.write_json(output_root / f"{stem}.json")
            report.write_markdown(output_root / f"{stem}.md")
        self._progress(f"Finished analysis for {target.entry_name}")
        return report

    def _fetch_complex_portal_context(
        self,
        *,
        target,
        homolog_context: dict[str, dict[str, Any]],
        notes: list[AnalysisNote],
        lookup_status: dict[str, dict[str, str]] | None = None,
    ) -> list[ComplexPortalComplex]:
        if not self.config.enable_complex_portal:
            self._progress("Skipping Complex Portal evidence (disabled)")
            if lookup_status is not None:
                lookup_status["human"] = {"status": "disabled"}
            return []

        complexes: dict[tuple[str, str], ComplexPortalComplex] = {}

        def add_records(records: list[ComplexPortalComplex]) -> None:
            for complex_item in records:
                key = (complex_item.species, complex_item.complex_ac)
                complexes[key] = complex_item

        try:
            self._progress("Fetching human Complex Portal evidence")
            bulk_fetch = getattr(self.complex_portal_client, "fetch_human_complexes_for_target", None)
            human_records = (
                bulk_fetch(accession=target.accession)
                if callable(bulk_fetch)
                else self.complex_portal_client.fetch_complexes_for_target(
                    accession=target.accession,
                    gene_symbol=target.gene_symbol,
                    taxon_id=9606,
                    entry_name=target.entry_name,
                    max_pages=self.config.complex_portal_max_pages,
                )
            )
            add_records(human_records)
            if lookup_status is not None:
                lookup_status["human"] = {"status": "ok" if human_records else "no_hits"}
        except Exception as exc:
            if lookup_status is not None:
                lookup_status["human"] = {"status": "error", "error": str(exc)}
            notes.append(
                AnalysisNote(
                    severity="warning",
                    message=f"Human Complex Portal lookup failed: {exc}",
                    source="Complex Portal",
                )
            )

        if not self._is_canonical_integrin_chain(target.gene_symbol):
            return sorted(complexes.values(), key=self._complex_portal_sort_key)

        mouse_target = (homolog_context.get("mouse") or {}).get("target")
        if mouse_target is None:
            return sorted(complexes.values(), key=self._complex_portal_sort_key)

        try:
            self._progress("Fetching mouse Complex Portal evidence")
            mouse_records = self.complex_portal_client.fetch_complexes_for_target(
                accession=mouse_target.accession,
                gene_symbol=mouse_target.gene_symbol,
                taxon_id=10090,
                entry_name=mouse_target.entry_name,
                max_pages=self.config.complex_portal_max_pages,
            )
            add_records(mouse_records)
            if lookup_status is not None:
                lookup_status["mouse"] = {"status": "ok" if mouse_records else "no_hits"}
        except Exception as exc:
            if lookup_status is not None:
                lookup_status["mouse"] = {"status": "error", "error": str(exc)}
            notes.append(
                AnalysisNote(
                    severity="warning",
                    message=f"Mouse Complex Portal lookup failed: {exc}",
                    source="Complex Portal",
                )
            )
        return sorted(complexes.values(), key=self._complex_portal_sort_key)

    def _complex_portal_sort_key(self, complex_item: ComplexPortalComplex) -> tuple[int, int, int, str]:
        species_rank = 0 if self._complex_portal_taxon_id(complex_item) == 9606 else 1
        return (species_rank, int(complex_item.predicted_complex), -(complex_item.confidence_score or 0), complex_item.name)

    def _progress(self, message: str) -> None:
        if self.config.verbose_progress:
            print(f"[agdesign2] {message}", flush=True)

    def _apply_surfy_topology_fallback(
        self,
        *,
        target: TargetRecord,
        topology: TopologyAnnotation | None,
    ) -> bool:
        surfy_topology = self._surfy_topology_for_target(target)
        return apply_surfy_topology_fallback(
            topology,
            surfy_topology,
            sequence_length=len(target.sequence),
        )

    def _surfy_topology_for_target(self, target: TargetRecord) -> str | None:
        index = self._load_surfy_topology_index()
        keys = [
            str(target.accession or "").upper(),
            str(target.entry_name or "").upper(),
            str(target.gene_symbol or "").upper(),
        ]
        return next((index[key] for key in keys if key and key in index), None)

    def _load_surfy_topology_index(self) -> dict[str, str]:
        if self._surfy_topology_index is not None:
            return self._surfy_topology_index
        self._surfy_topology_index = {}
        for path in (
            Path("data/inputs/accessible_human_topology_refined_accessible_only.csv"),
            Path.cwd() / "data" / "inputs" / "accessible_human_topology_refined_accessible_only.csv",
        ):
            if not path.exists():
                continue
            try:
                with path.open(newline="", encoding="utf-8") as handle:
                    for row in csv.DictReader(handle):
                        topology_text = str(row.get("surfy_topology") or "").strip()
                        if not topology_text:
                            continue
                        for key in (
                            row.get("accession"),
                            row.get("uniprot_entry_name"),
                            row.get("gene_symbol"),
                        ):
                            key_text = str(key or "").strip().upper()
                            if key_text:
                                self._surfy_topology_index[key_text] = topology_text
            except Exception:
                continue
            break
        return self._surfy_topology_index

    def _apply_secreted_universe_fallback(
        self,
        *,
        target: TargetRecord,
        topology: TopologyResult,
        features: list[Feature],
    ) -> AnalysisNote | None:
        if topology.ectodomain is not None:
            return None
        if topology.topology is not None and topology.topology.transmembrane_regions:
            return None
        if not self._target_is_curated_secreted(target):
            return None

        region = self._secreted_design_region_from_features(features, len(target.sequence))
        topology.ectodomain = region
        topology.topology = TopologyAnnotation(
            topology_class="secreted_external_evidence",
            extracellular_regions=[region],
            major_extracellular_region=region,
        )
        if region.metadata.get("fallback_reason") == "full_length":
            message = (
                "Curated target-universe evidence classifies this no-transmembrane target as secreted, but UniProt did "
                "not provide signal-peptide, chain, peptide, propeptide, or extracellular-topology boundaries; using the "
                "full canonical sequence as the secreted design region for construct review."
            )
        else:
            message = (
                "Curated target-universe evidence classifies this no-transmembrane target as secreted; using UniProt "
                f"{region.metadata.get('fallback_reason')} boundaries as the secreted design region."
            )
        return AnalysisNote(severity="warning", message=message, source="target-universe")

    def _target_is_curated_secreted(self, target: TargetRecord) -> bool:
        index = self._load_secreted_universe_index()
        keys = [
            str(target.accession or "").upper(),
            str(target.entry_name or "").upper(),
            str(target.gene_symbol or "").upper(),
        ]
        return any(key and key in index for key in keys)

    def _load_secreted_universe_index(self) -> set[str]:
        if self._secreted_universe_index is not None:
            return self._secreted_universe_index
        self._secreted_universe_index = set()
        for path in (
            Path("data/inputs/accessible_human_topology_refined_accessible_only.csv"),
            Path.cwd() / "data" / "inputs" / "accessible_human_topology_refined_accessible_only.csv",
        ):
            if not path.exists():
                continue
            try:
                with path.open(newline="", encoding="utf-8") as handle:
                    for row in csv.DictReader(handle):
                        labels = " ".join(
                            str(row.get(key) or "")
                            for key in (
                                "primary_topology_class",
                                "simplified_topology_bucket",
                                "final_accessibility_bucket",
                                "topology_class",
                                "topology_bucket",
                                "accessibility_class",
                                "classification",
                            )
                        ).lower()
                        if "secreted" not in labels:
                            continue
                        # 'secreted_isoform_only' classifies targets whose secreted
                        # evidence exists only for a non-canonical isoform. The
                        # canonical sequence we design is not the secreted species,
                        # so these must not enter the secreted-universe fallback
                        # (which would otherwise fabricate a canonical secreted
                        # design region). 'surface_isoform_only' lacks 'secreted'
                        # and is already excluded above.
                        if "isoform_only" in labels:
                            continue
                        for key in (
                            row.get("accession"),
                            row.get("uniprot_entry_name"),
                            row.get("gene_symbol"),
                        ):
                            key_text = str(key or "").strip().upper()
                            if key_text:
                                self._secreted_universe_index.add(key_text)
            except Exception:
                continue
            break
        return self._secreted_universe_index

    def _secreted_design_region_from_features(self, features: list[Feature], sequence_length: int) -> Region:
        chains = [
            feature
            for feature in features
            if feature.type.upper() in {"CHAIN", "PEPTIDE"}
            and 1 <= feature.start <= feature.end <= sequence_length
            and feature.end - feature.start + 1 >= self.config.min_construct_length
        ]
        if chains:
            selected = max(chains, key=lambda feature: (feature.end - feature.start + 1, -feature.start))
            return Region(
                start=selected.start,
                end=selected.end,
                label=selected.description or "Secreted chain",
                source="target-universe",
                confidence=0.75,
                metadata={"fallback_reason": selected.type.upper().lower()},
            )
        propeptide_end = max(
            (
                feature.end
                for feature in features
                if feature.type.upper() == "PROPEP" and 1 <= feature.start <= feature.end < sequence_length
            ),
            default=0,
        )
        if propeptide_end and sequence_length - propeptide_end >= self.config.min_construct_length:
            return Region(
                start=propeptide_end + 1,
                end=sequence_length,
                label="Secreted mature chain",
                source="target-universe",
                confidence=0.65,
                metadata={"fallback_reason": "propeptide"},
            )
        return Region(
            start=1,
            end=sequence_length,
            label="Secreted canonical sequence",
            source="target-universe",
            confidence=0.45,
            metadata={"fallback_reason": "full_length"},
        )

    def _build_residue_annotations(
        self,
        *,
        sequence: str,
        topology: TopologyAnnotation | None,
        residue_plddt: dict[int, float],
        surface_accessibility: dict[int, Any],
    ) -> list[ResidueAnnotation]:
        locations = self._topology_location_by_residue(len(sequence), topology)
        annotations: list[ResidueAnnotation] = []
        for position, amino_acid in enumerate(sequence, start=1):
            plddt = residue_plddt.get(position)
            accessibility = surface_accessibility.get(position)
            annotations.append(
                ResidueAnnotation(
                    position=position,
                    amino_acid=amino_acid,
                    topology_location=locations.get(position, "unknown"),
                    surface_exposed=accessibility.surface_exposed if accessibility is not None else None,
                    surface_accessibility=accessibility.classification if accessibility is not None else None,
                    surface_accessibility_reason=accessibility.reason if accessibility is not None else None,
                    full_model_sasa=accessibility.full_model_sasa if accessibility is not None else None,
                    full_model_relative_sasa=accessibility.full_model_relative_sasa if accessibility is not None else None,
                    domain_context_sasa=accessibility.domain_context_sasa if accessibility is not None else None,
                    domain_context_relative_sasa=accessibility.domain_context_relative_sasa if accessibility is not None else None,
                    pae_block_id=accessibility.pae_block_id if accessibility is not None else None,
                    plddt=round(float(plddt), 2) if plddt is not None else None,
                )
            )
        return annotations

    def _topology_location_by_residue(
        self,
        sequence_length: int,
        topology: TopologyAnnotation | None,
    ) -> dict[int, str]:
        locations = {position: "unknown" for position in range(1, sequence_length + 1)}
        if topology is None:
            return locations

        def mark(regions: list[Region], label: str) -> None:
            for region in regions:
                for position in range(max(1, region.start), min(sequence_length, region.end) + 1):
                    locations[position] = label

        mark(topology.extracellular_regions, "extracellular")
        if topology.major_extracellular_region is not None:
            mark([topology.major_extracellular_region], "extracellular")
        mark(topology.cytoplasmic_regions, "intracellular")
        mark(topology.transmembrane_regions, "membrane")
        mark(topology.intramembrane_regions, "membrane")

        if topology.topology_class.startswith("secreted") and topology.major_extracellular_region is not None:
            mark([topology.major_extracellular_region], "extracellular")
        return locations

    def _surface_accessibility_target_positions(
        self,
        *,
        sequence_length: int,
        topology: TopologyAnnotation | None,
        residue_names: dict[int, str],
        scope: Region | None,
    ) -> set[int]:
        locations = self._topology_location_by_residue(sequence_length, topology)
        if scope is None:
            scope_start = 1
            scope_end = sequence_length
        else:
            scope_start = max(1, scope.start)
            scope_end = min(sequence_length, scope.end)
        positions = {
            position
            for position in range(scope_start, scope_end + 1)
            if locations.get(position) == "extracellular"
        }
        positions.update(
            position
            for position, residue_name in residue_names.items()
            if scope_start <= position <= scope_end and residue_name == "CYS"
        )
        return positions

    def _surface_context_regions(
        self,
        *,
        structural_regions: list[Region],
        topology: TopologyAnnotation | None,
        sequence_length: int,
    ) -> list[Region]:
        strict_blocks = [
            region
            for region in structural_regions
            if str(region.metadata.get("stringency") or "").lower() == "strict"
            or str(region.metadata.get("segmentation_method") or "").lower() == "pae_boundary_scoring"
        ]
        blocks = self._non_overlapping_regions(strict_blocks)
        if blocks:
            return blocks
        fallback: list[Region] = []
        if topology is not None:
            fallback.extend(topology.extracellular_regions)
            fallback.extend(topology.transmembrane_regions)
            fallback.extend(topology.cytoplasmic_regions)
            if topology.major_extracellular_region is not None:
                fallback.append(topology.major_extracellular_region)
        blocks = self._non_overlapping_regions(fallback)
        if blocks:
            return blocks
        return [Region(start=1, end=sequence_length, label="Full AlphaFold model", source="AlphaFold")]

    def _non_overlapping_regions(self, regions: list[Region]) -> list[Region]:
        selected: list[Region] = []
        occupied: set[int] = set()
        for region in sorted(regions, key=lambda item: (item.start, item.end - item.start)):
            residues = set(range(region.start, region.end + 1))
            if residues & occupied:
                continue
            selected.append(region)
            occupied.update(residues)
        return selected

    def _extracellular_surface_positions(
        self,
        residue_annotations: list[ResidueAnnotation],
    ) -> set[int]:
        return {
            annotation.position
            for annotation in residue_annotations
            if annotation.topology_location == "extracellular" and annotation.surface_exposed is True
        }

    def _annotate_blast_hits_with_surface_identity(
        self,
        hits: list,
        *,
        query_region_start: int,
        extracellular_surface_positions: set[int],
    ) -> list:
        if not extracellular_surface_positions:
            return hits
        annotated = []
        for hit in hits:
            if not hit.query_alignment or not hit.subject_alignment:
                annotated.append(hit)
                continue
            alignment = type("BlastAlignment", (), {})()
            alignment.aligned_query = hit.query_alignment
            alignment.aligned_subject = hit.subject_alignment
            surface_result = identity_for_query_positions(
                alignment,
                query_region_start=query_region_start,
                query_positions=extracellular_surface_positions,
            )
            annotation = self._blast_query_alignment_annotation(
                hit.query_alignment,
                query_region_start=query_region_start,
                extracellular_surface_positions=extracellular_surface_positions,
            )
            annotated.append(
                replace(
                    hit,
                    query_alignment_annotation=annotation,
                    extracellular_surface_identity=(
                        round(surface_result.identity, 2)
                        if surface_result.identity is not None
                        else None
                    ),
                    extracellular_surface_aligned_positions=surface_result.aligned_positions,
                    extracellular_surface_matches=surface_result.matches,
                )
            )
        return annotated

    def _blast_query_alignment_annotation(
        self,
        query_alignment: str | None,
        *,
        query_region_start: int,
        extracellular_surface_positions: set[int],
    ) -> str | None:
        if not query_alignment:
            return None
        query_position = query_region_start - 1
        classes: list[str] = []
        for query_residue in query_alignment:
            if query_residue != "-":
                query_position += 1
            if query_residue != "-" and query_position in extracellular_surface_positions:
                classes.append("S")
            else:
                classes.append(".")
        return "".join(classes)

    def _annotate_homology_with_surface_identity(
        self,
        homology: list[HomologyRecord],
        *,
        homolog_context: dict[str, dict[str, Any]],
        ectodomain: Region | None,
        extracellular_surface_positions: set[int],
    ) -> list[HomologyRecord]:
        if not extracellular_surface_positions:
            return homology
        query_region_start = ectodomain.start if ectodomain is not None else 1
        annotated: list[HomologyRecord] = []
        for record in homology:
            info = homolog_context.get(record.species)
            alignment = info.get("alignment") if info else None
            if alignment is None:
                annotated.append(record)
                continue
            surface_result = identity_for_query_positions(
                alignment,
                query_region_start=query_region_start,
                query_positions=extracellular_surface_positions,
            )
            annotated.append(
                replace(
                    record,
                    extracellular_surface_identity=(
                        round(surface_result.identity, 2)
                        if surface_result.identity is not None
                        else None
                    ),
                    extracellular_surface_aligned_positions=surface_result.aligned_positions,
                )
            )
        return annotated

    def _annotate_family_context_with_surface_identity(
        self,
        *,
        family_context: FamilyContext | None,
        target_entry_name: str,
        ectodomain: Region | None,
        residue_annotations: list[ResidueAnnotation],
    ) -> FamilyContext | None:
        if family_context is None or ectodomain is None:
            return family_context
        accessible_positions = self._extracellular_surface_positions(residue_annotations)
        if not accessible_positions:
            return family_context
        pairwise = family_context.metadata.get("pairwise_alignments")
        if not isinstance(pairwise, list):
            return family_context
        rows: list[dict[str, Any]] = []
        for record in pairwise:
            if not isinstance(record, dict):
                continue
            query_label = str(record.get("query_label") or "")
            subject_label = str(record.get("subject_label") or "")
            aligned_query = str(record.get("aligned_query") or "")
            aligned_subject = str(record.get("aligned_subject") or "")
            if not aligned_query or not aligned_subject:
                continue
            target_is_query = query_label == target_entry_name
            target_is_subject = subject_label == target_entry_name
            if not target_is_query and not target_is_subject:
                continue
            alignment = type("FamilyAlignment", (), {})()
            if target_is_query:
                alignment.aligned_query = aligned_query
                alignment.aligned_subject = aligned_subject
                other_label = subject_label
            else:
                alignment.aligned_query = aligned_subject
                alignment.aligned_subject = aligned_query
                other_label = query_label
            surface_result = identity_for_query_positions(
                alignment,
                query_region_start=ectodomain.start,
                query_positions=accessible_positions,
            )
            rows.append(
                {
                    "target_label": target_entry_name,
                    "member_label": other_label,
                    "identity": round(surface_result.identity, 2)
                    if surface_result.identity is not None
                    else None,
                    "aligned_positions": surface_result.aligned_positions,
                    "matches": surface_result.matches,
                    "identity_denominator": "target extracellular accessible residues aligned to member",
                }
            )
        if not rows:
            return family_context
        metadata = dict(family_context.metadata)
        metadata["target_extracellular_accessible_identity"] = rows
        return replace(family_context, metadata=metadata)

    def _merge_canonical_family_sources(
        self,
        *,
        precomputed: CanonicalFamily | None,
        live: CanonicalFamily | None,
    ) -> CanonicalFamily | None:
        if precomputed is None:
            return live
        if live is None:
            return precomputed
        if precomputed.accession != live.accession:
            return precomputed
        merged_notes = list(precomputed.notes)
        merged_notes.extend(note for note in live.notes if note not in merged_notes)
        return CanonicalFamily(
            accession=precomputed.accession,
            name=precomputed.name or live.name,
            source_database=precomputed.source_database or live.source_database,
            fragment_count=live.fragment_count or precomputed.fragment_count,
            start=live.start if live.start is not None else precomputed.start,
            end=live.end if live.end is not None else precomputed.end,
            covered_length=live.covered_length if live.covered_length is not None else precomputed.covered_length,
            coverage_fraction=live.coverage_fraction if live.coverage_fraction is not None else precomputed.coverage_fraction,
            notes=merged_notes,
        )

    def _should_skip_live_family_context_fallback(
        self,
        *,
        precomputed_ortholog: Any,
        canonical_family: CanonicalFamily | None,
    ) -> bool:
        if precomputed_ortholog is None:
            return False
        if canonical_family is None:
            return False
        return str(canonical_family.accession or "").startswith("IPR")

    def _build_full_length_family_context(
        self,
        *,
        target,
        family_context: FamilyContext | None,
        canonical_family: CanonicalFamily | None,
    ) -> FamilyContext | None:
        if family_context is None or not family_context.members:
            return None
        resolved: list[tuple[FamilyMember, str]] = []
        seen: set[str] = set()

        def add_member(member: FamilyMember, sequence: str) -> None:
            label = member.entry_name or member.accession or member.gene_symbol
            if not label or label in seen or not sequence:
                return
            seen.add(label)
            resolved.append((member, sequence))

        for member in family_context.members:
            if member.entry_name == target.entry_name or member.accession == target.accession:
                add_member(
                    FamilyMember(
                        gene_symbol=target.gene_symbol or member.gene_symbol,
                        gene_name=target.protein_name,
                        accession=target.accession,
                        entry_name=target.entry_name,
                        sequence_length=len(target.sequence),
                        ectodomain_start=member.ectodomain_start,
                        ectodomain_end=member.ectodomain_end,
                        ectodomain_length=member.ectodomain_length,
                        notes=list(member.notes),
                    ),
                    target.sequence,
                )
                continue
            identifier = member.entry_name or member.accession
            if not identifier:
                continue
            try:
                entry = self.uniprot_client.fetch_entry_by_name(identifier)
                if entry is None:
                    continue
                member_target = self.uniprot_client._target_from_entry(entry)
            except Exception:
                continue
            add_member(
                FamilyMember(
                    gene_symbol=member_target.gene_symbol or member.gene_symbol,
                    gene_name=member_target.protein_name,
                    accession=member_target.accession,
                    entry_name=member_target.entry_name,
                    sequence_length=len(member_target.sequence),
                    ectodomain_start=member.ectodomain_start,
                    ectodomain_end=member.ectodomain_end,
                    ectodomain_length=member.ectodomain_length,
                    notes=list(member.notes),
                ),
                member_target.sequence,
            )

        if target.entry_name not in seen:
            add_member(
                FamilyMember(
                    gene_symbol=target.gene_symbol or target.entry_name.split("_", 1)[0],
                    gene_name=target.protein_name,
                    accession=target.accession,
                    entry_name=target.entry_name,
                    sequence_length=len(target.sequence),
                ),
                target.sequence,
            )
        if len(resolved) < 2:
            return None

        labels = [member.entry_name or member.accession or member.gene_symbol for member, _ in resolved]
        identity_matrix: list[list[float | None]] = []
        coverage_matrix: list[list[float | None]] = []
        for _, row_sequence in resolved:
            identity_row: list[float | None] = []
            coverage_row: list[float | None] = []
            row_length = len(row_sequence)
            for _, column_sequence in resolved:
                if not row_sequence or not column_sequence:
                    identity_row.append(None)
                    coverage_row.append(None)
                    continue
                alignment = global_align(row_sequence, column_sequence)
                identity_row.append(round(100.0 * alignment.matches / row_length, 2) if row_length else None)
                coverage_row.append(round(100.0 * alignment.aligned_positions / row_length, 2) if row_length else None)
            identity_matrix.append(identity_row)
            coverage_matrix.append(coverage_row)

        family_names = list(family_context.family_names)
        if not family_names and canonical_family is not None:
            family_names = [canonical_family.name]
        return FamilyContext(
            gene_symbol=target.gene_symbol or target.entry_name.split("_", 1)[0],
            source=f"{family_context.source}; full-length matrix recomputed locally",
            family_names=family_names,
            members=[member for member, _ in resolved],
            identity_matrix_labels=labels,
            identity_matrix=identity_matrix,
            coverage_matrix_labels=labels,
            coverage_matrix=coverage_matrix,
            metadata={
                **family_context.metadata,
                "matrix_scope": "full_length",
                "identity_denominator": "row full-length sequence length",
                "coverage_denominator": "row full-length sequence length",
            },
        )

    def _collect_species_matches(
        self,
        entry_name: str,
        gene_symbol: str,
        ectodomain: Region | None,
        sequence: str,
    ) -> tuple[list[HomologyRecord], list[HomologyRecord], dict[str, dict[str, Any]]]:
        matches: list[HomologyRecord] = []
        homology: list[HomologyRecord] = []
        context: dict[str, dict[str, Any]] = {}
        base_name = entry_name.split("_", 1)[0]
        for species, suffix in self.config.species_suffixes.items():
            if suffix == entry_name.split("_", 1)[-1]:
                continue
            if species == "macaca_fascicularis":
                refseq_match = self._find_macaca_refseq_match(
                    gene_symbol=gene_symbol,
                    human_sequence=sequence,
                    human_ectodomain=ectodomain,
                )
                if refseq_match is None:
                    record = HomologyRecord(
                        species=species,
                        accession=None,
                        entry_name=None,
                        ectodomain_start=None,
                        ectodomain_end=None,
                        available=False,
                        notes=["No RefSeq canonical isoform could be resolved for cynomolgus monkey."],
                    )
                    matches.append(record)
                    homology.append(record)
                    context[species] = {"available": False, "notes": record.notes}
                    continue
                target = refseq_match["target"]
                transferred_ectodomain = refseq_match["ectodomain"]
                base_notes = refseq_match["notes"]
                match = HomologyRecord(
                    species=species,
                    accession=target.accession,
                    entry_name=target.entry_name,
                    ectodomain_start=transferred_ectodomain.start if transferred_ectodomain else None,
                    ectodomain_end=transferred_ectodomain.end if transferred_ectodomain else None,
                    available=True,
                    notes=base_notes,
                )
                matches.append(match)
                if ectodomain and transferred_ectodomain and refseq_match["ectodomain_sequence"]:
                    query_ecto = slice_sequence(sequence, ectodomain.start, ectodomain.end)
                    subject_ecto = refseq_match["ectodomain_sequence"]
                    alignment = global_align(query_ecto, subject_ecto)
                    context[species] = {
                        "available": True,
                        "target": target,
                        "ectodomain": transferred_ectodomain,
                        "ectodomain_sequence": subject_ecto,
                        "alignment": alignment,
                        "notes": match.notes,
                    }
                    homology.append(
                        replace(
                            match,
                            identity=round(alignment.identity, 2),
                            coverage=round(alignment.coverage, 2),
                        )
                    )
                else:
                    context[species] = {
                        "available": True,
                        "target": target,
                        "ectodomain": transferred_ectodomain,
                        "ectodomain_sequence": refseq_match["ectodomain_sequence"],
                        "alignment": None,
                        "notes": match.notes,
                    }
                    homology.append(match)
                continue
            entry = self.uniprot_client.fetch_same_name_species_match(base_name, suffix)
            if not entry:
                gene_match = self._find_species_match_by_gene_symbol(species=species, gene_symbol=gene_symbol)
                if gene_match is not None:
                    entry = gene_match
                    base_notes = [
                        "No same-name UniProt entry found for species.",
                        f"Using UniProt gene-symbol search match in taxon {self.config.species_taxonomy[species]}.",
                    ]
                else:
                    fallback = self._find_species_match_by_blast(species=species, query_sequence=sequence)
                    if fallback is None:
                        record = HomologyRecord(
                            species=species,
                            accession=None,
                            entry_name=None,
                            ectodomain_start=None,
                            ectodomain_end=None,
                            available=False,
                            notes=[
                                "No same-name UniProt entry found for species.",
                                "No gene-symbol UniProt match found in target taxon.",
                                "No proteome BLAST fallback hit could be resolved.",
                            ],
                        )
                        matches.append(record)
                        homology.append(record)
                        context[species] = {"available": False, "notes": record.notes}
                        continue
                    entry = fallback["entry"]
                    base_notes = [
                        "No same-name UniProt entry found for species.",
                        "No gene-symbol UniProt match found in target taxon.",
                        f"Using top proteome BLAST hit as homolog candidate: {fallback['hit'].subject_id}.",
                    ]
            else:
                base_notes = []
            target = self.uniprot_client._target_from_entry(entry)
            species_features = self.uniprot_client.get_features(entry)
            topology = derive_ectodomain(species_features, len(target.sequence), self.config)
            match = HomologyRecord(
                species=species,
                accession=target.accession,
                entry_name=target.entry_name,
                ectodomain_start=topology.ectodomain.start if topology.ectodomain else None,
                ectodomain_end=topology.ectodomain.end if topology.ectodomain else None,
                available=True,
                notes=base_notes + [note.message for note in topology.notes],
            )
            matches.append(match)
            if ectodomain and topology.ectodomain:
                query_ecto = slice_sequence(sequence, ectodomain.start, ectodomain.end)
                subject_ecto = slice_sequence(target.sequence, topology.ectodomain.start, topology.ectodomain.end)
                alignment = global_align(query_ecto, subject_ecto)
                context[species] = {
                    "available": True,
                    "target": target,
                    "ectodomain": topology.ectodomain,
                    "ectodomain_sequence": subject_ecto,
                    "alignment": alignment,
                    "notes": match.notes,
                }
                homology.append(
                    replace(
                        match,
                        identity=round(alignment.identity, 2),
                        coverage=round(alignment.coverage, 2),
                    )
                )
            else:
                context[species] = {
                    "available": True,
                    "target": target,
                    "ectodomain": topology.ectodomain,
                    "ectodomain_sequence": None,
                    "alignment": None,
                    "notes": match.notes,
                }
                homology.append(match)
        return matches, homology, context

    def _find_species_match_by_blast(self, *, species: str, query_sequence: str) -> dict[str, Any] | None:
        hits = self.blast_client.search_species(
            query_sequence,
            species=species,
            ortholog_search=(species == "macaca_fascicularis"),
        )
        for hit in hits:
            accession = extract_accession(hit.subject_id)
            entry_name = extract_entry_name(hit.subject_id)
            try:
                entry = None
                if accession:
                    entry = self.uniprot_client._fetch_entry(accession)
                elif entry_name:
                    entry = self.uniprot_client.fetch_entry_by_name(entry_name)
                if entry is not None:
                    return {"hit": hit, "entry": entry}
            except Exception:
                continue
        return None

    def _find_species_match_by_gene_symbol(self, *, species: str, gene_symbol: str) -> dict[str, Any] | None:
        taxon_id = self.config.species_taxonomy[species]
        queries = [f"(gene_exact:{gene_symbol}) AND (organism_id:{taxon_id}) AND (reviewed:true)"]
        for query in queries:
            payload = self.uniprot_client.search(query, size=5)
            results = payload.get("results", [])
            for result in results:
                accession = result.get("primaryAccession")
                if not accession:
                    continue
                try:
                    return self.uniprot_client._fetch_entry(accession)
                except Exception:
                    continue
        return None

    def _find_macaca_refseq_match(
        self,
        *,
        gene_symbol: str,
        human_sequence: str,
        human_ectodomain: Region | None,
    ) -> dict[str, Any] | None:
        result = self.refseq_client.fetch_canonical_protein(
            gene_symbol=gene_symbol,
            organism="Macaca fascicularis",
        )
        if result is None:
            return None
        notes = [note.message for note in result.notes]
        if not human_sequence or human_ectodomain is None:
            notes.append("Unable to transfer ectodomain boundaries from the human target.")
            return {
                "target": result.target,
                "ectodomain": None,
                "ectodomain_sequence": None,
                "notes": notes,
            }
        full_alignment = global_align(human_sequence, result.target.sequence)
        mapped_start, mapped_end, mapped_sequence, mapping_notes = map_query_region_to_subject(
            full_alignment,
            query_start=human_ectodomain.start,
            query_end=human_ectodomain.end,
            subject_sequence=result.target.sequence,
        )
        notes.extend(mapping_notes)
        if mapped_start is None or mapped_end is None or mapped_sequence is None:
            notes.append("Unable to map the human ectodomain onto the RefSeq cynomolgus monkey protein.")
            return {
                "target": result.target,
                "ectodomain": None,
                "ectodomain_sequence": None,
                "notes": notes,
            }
        return {
            "target": result.target,
            "ectodomain": Region(
                start=mapped_start,
                end=mapped_end,
                label="Transferred ectodomain",
                source="RefSeq-alignment",
                confidence=0.7,
            ),
            "ectodomain_sequence": mapped_sequence,
            "notes": notes,
        }

    def _filter_experimental_constructs(
        self,
        constructs: list[ExperimentalConstruct],
        ectodomain: Region | None,
        topology: TopologyAnnotation | None = None,
    ) -> list[ExperimentalConstruct]:
        if topology is not None and topology.topology_class.startswith("multipass"):
            return constructs
        if ectodomain is None:
            return constructs
        filtered: list[ExperimentalConstruct] = []
        for construct in constructs:
            ectodomain_chains = [
                chain
                for chain in construct.chains
                if chain.start >= ectodomain.start and chain.end <= ectodomain.end
            ]
            if not ectodomain_chains:
                continue
            filtered.append(replace(construct, chains=ectodomain_chains))
        return filtered

    def _filter_design_annotations(
        self,
        annotations: list[DomainAnnotation],
        ectodomain: Region | None,
    ) -> list[DomainAnnotation]:
        if ectodomain is None:
            return []
        contained = [
            annotation
            for annotation in annotations
            if annotation.start >= ectodomain.start
            and annotation.end <= ectodomain.end
            and annotation.type in {"domain", "repeat", "homologous_superfamily"}
        ]
        primary: dict[tuple[int, int, str], DomainAnnotation] = {}
        for annotation in contained:
            key = (annotation.start, annotation.end, annotation.type)
            current = primary.get(key)
            if current is None:
                primary[key] = annotation
                continue
            current_priority = self._annotation_priority(current)
            new_priority = self._annotation_priority(annotation)
            if new_priority > current_priority:
                primary[key] = annotation
        ranked = sorted(
            primary.values(),
            key=lambda item: (
                -self._annotation_priority(item)[0],
                -self._annotation_priority(item)[1],
                item.end - item.start + 1,
                item.start,
                item.end,
            ),
        )
        curated: list[DomainAnnotation] = []
        for annotation in ranked:
            if any(self._annotations_are_redundant(annotation, kept) for kept in curated):
                continue
            curated.append(annotation)
        return sorted(curated, key=lambda item: (item.start, item.end, item.source_database, item.accession))

    def _fetch_gpcr_annotation(
        self,
        *,
        target,
        topology: TopologyAnnotation | None,
        interpro_annotations: list[DomainAnnotation],
        notes: list[AnalysisNote],
    ) -> GPCRdbAnnotation | None:
        if not self.config.enable_gpcrdb:
            return None
        if not self._is_gpcr_candidate(target=target, topology=topology, interpro_annotations=interpro_annotations):
            return None
        try:
            self._progress("Fetching GPCRdb annotations")
            annotation = self.gpcrdb_client.fetch_annotation(
                entry_name=target.entry_name,
                accession=target.accession,
                sequence=target.sequence,
            )
        except Exception as exc:
            notes.append(
                AnalysisNote(
                    severity="warning",
                    message=f"GPCRdb annotation lookup failed: {exc}",
                    source="GPCRdb",
                )
            )
            return None
        if annotation is None:
            if topology is not None and len(topology.transmembrane_regions) >= self.config.multipass_gpcr_like_tm_count:
                notes.append(
                    AnalysisNote(
                        severity="warning",
                        message="Target has 7TM-like topology, but no GPCRdb annotation was available; using generic multipass logic.",
                        source="GPCRdb",
                    )
                )
            return None
        return annotation

    def _is_gpcr_candidate(
        self,
        *,
        target,
        topology: TopologyAnnotation | None,
        interpro_annotations: list[DomainAnnotation],
    ) -> bool:
        text = " ".join(
            [
                str(target.entry_name or ""),
                str(target.gene_symbol or ""),
                str(target.protein_name or ""),
                " ".join(annotation.name for annotation in interpro_annotations),
            ]
        ).lower()
        if any(token in text for token in ("g protein-coupled receptor", "gpcr", "frizzled", "adhesion g")):
            return True
        return bool(topology is not None and len(topology.transmembrane_regions) >= self.config.multipass_gpcr_like_tm_count)

    def _select_canonical_family(
        self,
        annotations: list[DomainAnnotation],
        *,
        sequence_length: int,
        protein_name: str,
        gene_symbol: str | None,
    ) -> CanonicalFamily | None:
        family_annotations = [annotation for annotation in annotations if annotation.type == "family"]
        if not family_annotations:
            return None

        ranked_pool = [annotation for annotation in family_annotations if annotation.source_database == "INTERPRO"]
        if not ranked_pool:
            return None
        grouped: dict[tuple[str, str, str], list[DomainAnnotation]] = {}
        for annotation in ranked_pool:
            accession = annotation.accession
            name = annotation.name
            source_database = annotation.source_database
            if source_database != "INTERPRO" and annotation.integrated_accession and annotation.integrated_name:
                accession = annotation.integrated_accession
                name = annotation.integrated_name
                source_database = "INTERPRO"
            grouped.setdefault((source_database, accession, name), []).append(annotation)

        candidates: list[tuple[float, CanonicalFamily]] = []
        for (source_database, accession, name), items in grouped.items():
            merged_intervals = self._merge_annotation_intervals(items)
            if not merged_intervals:
                continue
            covered_length = sum(end - start + 1 for start, end in merged_intervals)
            coverage_fraction = (covered_length / sequence_length) if sequence_length > 0 else 0.0
            name_similarity = self._annotation_name_similarity(
                family_name=name,
                protein_name=protein_name,
                gene_symbol=gene_symbol,
            )
            score = 0.0
            score += 100.0 if source_database == "INTERPRO" else 50.0
            score += min(25.0, coverage_fraction * 25.0)
            score += min(5.0, len(merged_intervals) - 1)
            score += 3.0 if accession.startswith("IPR") else 0.0
            if len(grouped) > 1:
                if coverage_fraction >= 0.85 and name_similarity >= 0.60:
                    score -= 30.0
                elif coverage_fraction >= 0.65 and name_similarity >= 0.75:
                    score -= 18.0
            candidates.append(
                (
                    score,
                    CanonicalFamily(
                        accession=accession,
                        name=name,
                        source_database=source_database,
                        fragment_count=len(merged_intervals),
                        start=min(start for start, _ in merged_intervals),
                        end=max(end for _, end in merged_intervals),
                        covered_length=covered_length,
                        coverage_fraction=round(coverage_fraction * 100.0, 2),
                        notes=[
                            f"Coverage across merged family fragments: {covered_length} aa ({coverage_fraction * 100.0:.1f}% of sequence).",
                            f"Name similarity to target protein label: {name_similarity:.2f}.",
                        ],
                    ),
                )
            )

        if not candidates:
            return None
        candidates.sort(
            key=lambda item: (
                item[0],
                item[1].coverage_fraction or 0.0,
                -(item[1].fragment_count),
                item[1].accession,
            ),
            reverse=True,
        )
        return candidates[0][1]

    def _merge_annotation_intervals(
        self,
        annotations: list[DomainAnnotation],
    ) -> list[tuple[int, int]]:
        intervals = sorted((annotation.start, annotation.end) for annotation in annotations)
        if not intervals:
            return []
        merged = [intervals[0]]
        for start, end in intervals[1:]:
            last_start, last_end = merged[-1]
            if start <= last_end + 1:
                merged[-1] = (last_start, max(last_end, end))
            else:
                merged.append((start, end))
        return merged

    def _annotation_name_similarity(
        self,
        *,
        family_name: str,
        protein_name: str,
        gene_symbol: str | None,
    ) -> float:
        normalized_family = self._normalize_label_for_similarity(family_name)
        normalized_protein = self._normalize_label_for_similarity(protein_name)
        if not normalized_family or not normalized_protein:
            return 0.0
        similarity = SequenceMatcher(None, normalized_family, normalized_protein).ratio()
        if gene_symbol:
            normalized_symbol = self._normalize_label_for_similarity(gene_symbol)
            if normalized_symbol and normalized_symbol in normalized_family:
                similarity += 0.05
        return min(similarity, 1.0)

    def _normalize_label_for_similarity(self, value: str | None) -> str:
        text = re.sub(r"[^a-z0-9]+", " ", str(value or "").lower())
        return " ".join(token for token in text.split() if token)

    def _annotation_priority(self, annotation: DomainAnnotation) -> tuple[int, int]:
        type_priority = {
            "domain": 4,
            "repeat": 3,
            "homologous_superfamily": 2,
        }.get(annotation.type, 1)
        source_priority = 2 if annotation.source_database == "INTERPRO" else 1
        return (type_priority * 10 + source_priority, 1 if annotation.representative else 0)

    def _annotations_are_redundant(
        self,
        candidate: DomainAnnotation,
        selected: DomainAnnotation,
    ) -> bool:
        overlap = self._annotation_overlap_fraction(candidate, selected)
        if overlap <= 0.0:
            return False
        candidate_key = candidate.integrated_accession or candidate.accession
        selected_key = selected.integrated_accession or selected.accession
        if candidate_key == selected_key and overlap >= 0.6:
            return True
        if candidate.source_database == "PFAM" and selected.source_database == "INTERPRO" and overlap >= 0.8:
            return True
        if candidate.type == "homologous_superfamily" and selected.type in {"domain", "repeat"} and overlap >= 0.6:
            return True
        if candidate.type == "repeat" and selected.type == "domain" and overlap >= 0.85:
            return True
        return False

    def _annotation_overlap_fraction(
        self,
        left: DomainAnnotation,
        right: DomainAnnotation,
    ) -> float:
        overlap_start = max(left.start, right.start)
        overlap_end = min(left.end, right.end)
        if overlap_start > overlap_end:
            return 0.0
        overlap_length = overlap_end - overlap_start + 1
        shorter_length = min(left.end - left.start + 1, right.end - right.start + 1)
        if shorter_length <= 0:
            return 0.0
        return overlap_length / shorter_length

    def _annotation_construct_name(self, annotation: DomainAnnotation) -> str:
        base = re.sub(r"[^a-z0-9]+", "_", annotation.name.lower()).strip("_")
        prefix = "repeat" if annotation.type == "repeat" else "domain"
        return f"{prefix}_{base or annotation.accession.lower()}"

    def _annotation_seed_score(self, annotation: DomainAnnotation) -> float:
        base_score = 82.0 if annotation.source_database == "INTERPRO" else 80.0
        if annotation.representative:
            base_score += 2.0
        if annotation.type == "repeat":
            base_score -= 1.0
        return base_score

    def _annotation_seed_classification(self, annotation: DomainAnnotation) -> str:
        if annotation.type == "repeat":
            return "repeat_module"
        return "complete_domain"

    def _annotation_label(self, annotation: DomainAnnotation) -> str:
        suffix = " representative" if annotation.representative else ""
        return f"{annotation.source_database} {annotation.accession} {annotation.name}{suffix}"

    def _classify_construct(
        self,
        start: int,
        end: int,
        annotations: list[DomainAnnotation],
    ) -> tuple[str, str, list[DomainAnnotation], list[DomainAnnotation]]:
        overlapping = [
            annotation
            for annotation in annotations
            if annotation.end >= start and annotation.start <= end
        ]
        complete = [
            annotation
            for annotation in overlapping
            if annotation.start >= start and annotation.end <= end
        ]
        clipped = [
            annotation
            for annotation in overlapping
            if not (annotation.start >= start and annotation.end <= end)
        ]
        complete_domains = [annotation for annotation in complete if annotation.type == "domain"]
        complete_repeats = [annotation for annotation in complete if annotation.type == "repeat"]
        complete_superfamilies = [
            annotation for annotation in complete if annotation.type in {"homologous_superfamily", "family"}
        ]
        if clipped:
            label = "partial_domain"
            summary = (
                "Truncates annotated architecture: "
                + ", ".join(self._short_annotation_label(annotation) for annotation in clipped[:4])
            )
        elif complete_repeats and not complete_domains:
            label = "repeat_module"
            summary = (
                "Captures annotated repeat module(s): "
                + ", ".join(self._short_annotation_label(annotation) for annotation in complete_repeats[:4])
            )
        elif len(complete_domains) >= 2:
            label = "multi_domain_unit"
            summary = (
                "Captures multiple complete annotated domains: "
                + ", ".join(self._short_annotation_label(annotation) for annotation in complete_domains[:4])
            )
        elif len(complete_domains) == 1:
            label = "complete_domain"
            summary = f"Captures one complete annotated domain: {self._short_annotation_label(complete_domains[0])}"
        elif complete_superfamilies:
            label = "family_module"
            summary = (
                "Matches broader annotated family/superfamily region: "
                + ", ".join(self._short_annotation_label(annotation) for annotation in complete_superfamilies[:4])
            )
        else:
            label = "structured_region"
            summary = "No complete InterPro/Pfam domain match; treat as structure-led region."
        return label, summary, complete, clipped

    def _short_annotation_label(self, annotation: DomainAnnotation) -> str:
        return f"{annotation.name} ({annotation.accession})"

    def _apply_domain_annotation_context(
        self,
        construct: ConstructSuggestion,
        annotations: list[DomainAnnotation],
    ) -> ConstructSuggestion:
        classification, summary, complete_annotations, clipped_annotations = self._classify_construct(
            construct.start,
            construct.end,
            annotations,
        )
        evidence = list(construct.evidence)
        warnings = list(construct.warnings)
        score = construct.score
        if complete_annotations:
            evidence.append(summary)
        if clipped_annotations:
            warnings.append(
                "Clips annotated domain architecture: "
                + ", ".join(self._short_annotation_label(annotation) for annotation in clipped_annotations[:4])
            )
            score -= 8.0
        elif classification == "complete_domain":
            score += 5.0
        elif classification in {"multi_domain_unit", "repeat_module"}:
            score += 3.0
        return replace(
            construct,
            score=round(score, 1),
            classification=classification,
            evidence=self._unique_strings(evidence),
            warnings=self._unique_strings(warnings),
        )

    def _unique_strings(self, values: list[str]) -> list[str]:
        seen: set[str] = set()
        unique: list[str] = []
        for value in values:
            if value in seen:
                continue
            seen.add(value)
            unique.append(value)
        return unique

    def _suggest_constructs(
        self,
        *,
        features: list[Feature],
        topology: TopologyAnnotation | None,
        ectodomain: Region | None,
        sequence_length: int,
        structural_regions: list[Region],
        experimental_constructs: list[ExperimentalConstruct],
        cysteine_analysis: list[CysteineFinding],
        interpro_annotations: list[DomainAnnotation],
        residue_plddt: dict[int, float] | None = None,
        ligand_interactions: list[Feature] | None = None,
        ptms: list[Feature] | None = None,
        gpcr_annotation: GPCRdbAnnotation | None = None,
    ) -> list[ConstructSuggestion]:
        if not ectodomain and not topology:
            return []

        warnings = [finding.warning for finding in cysteine_analysis if finding.warning]
        suggestions: list[ConstructSuggestion] = []
        if ectodomain is not None:
            suggestions.append(
                ConstructSuggestion(
                    name="full_ectodomain",
                    start=ectodomain.start,
                    end=ectodomain.end,
                    score=90.0 - (5.0 if warnings else 0.0),
                    rationale="Complete extracellular domain based on topology annotations.",
                    warnings=warnings,
                    evidence=["Topology-derived ectodomain"],
                    classification="soluble_ecd",
                )
            )

        if ectodomain is not None:
            design_annotations = self._filter_design_annotations(interpro_annotations, ectodomain)
            for annotation in design_annotations:
                suggestions.append(
                    ConstructSuggestion(
                        name=self._annotation_construct_name(annotation),
                        start=annotation.start,
                        end=annotation.end,
                        score=self._annotation_seed_score(annotation),
                        rationale=(
                            f"Representative {annotation.source_database} {annotation.type.replace('_', ' ')} "
                            f"annotation for {annotation.name}."
                        ),
                        classification=self._annotation_seed_classification(annotation),
                        evidence=[self._annotation_label(annotation)],
                    )
                )

            domain_features = [
                feature
                for feature in features
                if feature.type.upper() in {"DOMAIN", "REGION"}
                and feature.start >= ectodomain.start
                and feature.end <= ectodomain.end
            ]
            for feature in domain_features:
                suggestions.append(
                    ConstructSuggestion(
                        name=(feature.description or feature.type).lower().replace(" ", "_"),
                        start=feature.start,
                        end=feature.end,
                        score=80.0,
                        rationale="Curated extracellular domain annotation from UniProt.",
                        evidence=["UniProt domain annotation"],
                    )
                )
            lenient_index = 1
            strict_index = 1
            for region in structural_regions:
                if region.start < ectodomain.start or region.end > ectodomain.end:
                    continue
                stringency = region.metadata.get("stringency")
                if stringency == "lenient":
                    name = f"lenient_region_{lenient_index}"
                    lenient_index += 1
                    rationale = (
                        "Less-stringent AlphaFold structured region intended to preserve larger "
                        "cohesive extracellular units, potentially spanning multiple domains."
                    )
                    evidence = ["AlphaFold lenient structured region"]
                    score = 74.0
                elif stringency == "strict":
                    name = f"strict_domain_{strict_index}"
                    strict_index += 1
                    rationale = (
                        "More-stringent AlphaFold segmentation intended to isolate likely "
                        "single extracellular domains."
                    )
                    evidence = ["AlphaFold strict domain segmentation"]
                    score = 76.0
                else:
                    name = f"structural_region_{strict_index}"
                    strict_index += 1
                    rationale = "AlphaFold-defined extracellular structural region."
                    evidence = ["AlphaFold structure segmentation"]
                    score = 75.0
                suggestions.append(
                    ConstructSuggestion(
                        name=name,
                        start=region.start,
                        end=region.end,
                        score=score,
                        rationale=rationale,
                        evidence=evidence,
                    )
                )
        suggestions.extend(
            self._suggest_pdb_backed_constructs(
                experimental_constructs=experimental_constructs,
                ectodomain=ectodomain,
                topology=topology,
                sequence_length=sequence_length,
            )
        )
        if topology is not None:
            suggestions.extend(
                self._suggest_membrane_expression_constructs(
                    topology=topology,
                    sequence_length=sequence_length,
                    residue_plddt=residue_plddt,
                    features=features,
                    interpro_annotations=interpro_annotations,
                    ligand_interactions=ligand_interactions or [],
                    ptms=ptms or [],
                    gpcr_annotation=gpcr_annotation,
                )
            )
        deduped = self._dedupe_constructs(suggestions)
        if ectodomain is not None:
            design_annotations = self._filter_design_annotations(interpro_annotations, ectodomain)
            deduped = [
                construct
                if construct.classification == "membrane_expression"
                else self._apply_domain_annotation_context(construct, design_annotations)
                for construct in deduped
            ]
        filtered: list[ConstructSuggestion] = []
        for construct in deduped:
            if construct.classification == "membrane_expression":
                filtered.append(construct)
                continue
            if ectodomain is None:
                continue
            if construct.start >= ectodomain.start and construct.end <= ectodomain.end:
                filtered.append(construct)
        deduped = [
            construct
            for construct in filtered
            if (construct.end - construct.start + 1) >= self.config.min_construct_length
        ]
        return sorted(deduped, key=lambda item: (-item.score, item.start, item.end))

    def _suggest_pdb_backed_constructs(
        self,
        *,
        experimental_constructs: list[ExperimentalConstruct],
        ectodomain: Region | None,
        topology: TopologyAnnotation | None,
        sequence_length: int,
    ) -> list[ConstructSuggestion]:
        suggestions: list[ConstructSuggestion] = []
        is_multipass = bool(topology and topology.topology_class.startswith("multipass"))
        for construct in experimental_constructs:
            for chain in construct.chains:
                if chain.start < 1 or chain.end > sequence_length or chain.start > chain.end:
                    continue
                classification: str | None = None
                rationale = "Experimental structural construct mapped from PDB chain coverage."
                inside_design_region = (
                    ectodomain is not None
                    and chain.start >= ectodomain.start
                    and chain.end <= ectodomain.end
                )
                if inside_design_region:
                    pass
                elif is_multipass:
                    classification = "membrane_expression"
                    rationale = "Experimental membrane-protein construct mapped from PDB chain coverage."
                else:
                    continue
                evidence = [f"PDB {construct.pdb_id}"]
                if construct.method:
                    evidence.append(f"Method: {construct.method}")
                if construct.resolution:
                    evidence.append(f"Resolution: {construct.resolution}")
                suggestions.append(
                    ConstructSuggestion(
                        name=self._pdb_construct_name(construct, chain),
                        start=chain.start,
                        end=chain.end,
                        score=85.0,
                        rationale=rationale,
                        classification=classification,
                        evidence=evidence,
                        pdb_id=construct.pdb_id,
                    )
                )
        return suggestions

    def _pdb_construct_name(self, construct: ExperimentalConstruct, chain: Region) -> str:
        chain_label = str(chain.metadata.get("chain") or chain.label or "x").lower()
        chain_label = re.sub(r"[^a-z0-9]+", "_", chain_label).strip("_") or "x"
        return f"pdb_{construct.pdb_id.lower()}_{chain_label}"

    def _suggest_membrane_expression_constructs(
        self,
        *,
        topology: TopologyAnnotation,
        sequence_length: int,
        residue_plddt: dict[int, float] | None = None,
        features: list[Feature] | None = None,
        interpro_annotations: list[DomainAnnotation] | None = None,
        ligand_interactions: list[Feature] | None = None,
        ptms: list[Feature] | None = None,
        gpcr_annotation: GPCRdbAnnotation | None = None,
    ) -> list[ConstructSuggestion]:
        if not topology.topology_class.startswith("multipass") or topology.membrane_core_region is None:
            return []
        features = features or []
        interpro_annotations = interpro_annotations or []
        ligand_interactions = ligand_interactions or []
        ptms = ptms or []
        suggestions: list[ConstructSuggestion] = []
        membrane_core = topology.membrane_core_region
        transmembrane_count = len(topology.transmembrane_regions)
        n_trim = self._determine_multipass_n_trim_start(
            topology=topology,
            sequence_length=sequence_length,
            residue_plddt=residue_plddt,
            features=features,
            interpro_annotations=interpro_annotations,
            ligand_interactions=ligand_interactions,
            ptms=ptms,
        )
        c_trim = self._determine_multipass_c_trim_end(
            topology=topology,
            sequence_length=sequence_length,
            residue_plddt=residue_plddt,
            features=features,
            interpro_annotations=interpro_annotations,
            ligand_interactions=ligand_interactions,
            ptms=ptms,
            gpcr_annotation=gpcr_annotation,
        )
        full_length_warnings = []
        if n_trim.blocked_reason:
            full_length_warnings.append(n_trim.blocked_reason)
        if c_trim.blocked_reason:
            full_length_warnings.append(c_trim.blocked_reason)
        suggestions.append(
            ConstructSuggestion(
                name="membrane_expression_full_length",
                start=1,
                end=sequence_length,
                score=82.0 if (n_trim.enabled or c_trim.enabled) else 86.0,
                rationale=(
                    "Baseline native membrane-expression construct retaining the full membrane-spanning protein sequence."
                ),
                classification="membrane_expression",
                evidence=[
                    "Native membrane-expression baseline.",
                    f"Topology-derived membrane core: {membrane_core.start}-{membrane_core.end}.",
                    f"Detected {transmembrane_count} transmembrane segment(s).",
                ],
                warnings=self._unique_strings(full_length_warnings),
            )
        )
        if c_trim.enabled and c_trim.end is not None and c_trim.end < sequence_length:
            suggestions.append(
                ConstructSuggestion(
                    name="membrane_expression_c_tail_trimmed",
                    start=1,
                    end=c_trim.end,
                    score=88.0,
                    rationale=(
                        "Native membrane-expression construct retaining the complete transmembrane core while "
                        "removing a long, likely disordered cytoplasmic C-terminal tail."
                    ),
                    classification="membrane_expression",
                    evidence=c_trim.evidence,
                    warnings=c_trim.warnings,
                )
            )
        if n_trim.enabled and n_trim.start is not None and n_trim.start > 1:
            suggestions.append(
                ConstructSuggestion(
                    name="membrane_expression_n_tail_trimmed",
                    start=n_trim.start,
                    end=sequence_length,
                    score=86.0,
                    rationale=(
                        "Native membrane-expression construct retaining the complete transmembrane core while "
                        "removing a long, likely disordered cytoplasmic N-terminal tail."
                    ),
                    classification="membrane_expression",
                    evidence=n_trim.evidence,
                    warnings=n_trim.warnings,
                )
            )
        if (
            n_trim.enabled
            and c_trim.enabled
            and n_trim.start is not None
            and c_trim.end is not None
            and n_trim.start > 1
            and c_trim.end < sequence_length
            and n_trim.start <= c_trim.end
        ):
            suggestions.append(
                ConstructSuggestion(
                    name="membrane_expression_terminal_trimmed",
                    start=n_trim.start,
                    end=c_trim.end,
                    score=89.0,
                    rationale=(
                        "Native membrane-expression construct retaining the membrane-spanning core while "
                        "removing both long, likely disordered cytoplasmic terminal regions."
                    ),
                    classification="membrane_expression",
                    evidence=self._unique_strings(n_trim.evidence + c_trim.evidence),
                    warnings=self._unique_strings(n_trim.warnings + c_trim.warnings),
                )
            )
        suggestions.extend(
            self._suggest_gpcr_c_tail_series(
                topology=topology,
                sequence_length=sequence_length,
                features=features,
                interpro_annotations=interpro_annotations,
                ligand_interactions=ligand_interactions,
                ptms=ptms,
                gpcr_annotation=gpcr_annotation,
            )
        )
        return suggestions

    def _suggest_gpcr_c_tail_series(
        self,
        *,
        topology: TopologyAnnotation,
        sequence_length: int,
        features: list[Feature],
        interpro_annotations: list[DomainAnnotation],
        ligand_interactions: list[Feature],
        ptms: list[Feature],
        gpcr_annotation: GPCRdbAnnotation | None,
    ) -> list[ConstructSuggestion]:
        if gpcr_annotation is None:
            return []
        terminal_tm = max(topology.transmembrane_regions, key=lambda region: region.end, default=None)
        if terminal_tm is None:
            return []
        anchor = terminal_tm.end
        anchor_label = f"TM{len(topology.transmembrane_regions)}"
        h8 = self._gpcr_segment(gpcr_annotation, "H8")
        if h8 is not None and h8.end > anchor:
            anchor = h8.end
            anchor_label = f"GPCRdb helix 8 ({h8.start}-{h8.end})"
        variants: list[ConstructSuggestion] = []
        for buffer_size in self.config.gpcr_c_tail_series_buffers:
            end = min(sequence_length, anchor + int(buffer_size))
            removed_start = end + 1
            removed_length = sequence_length - end
            if removed_length < 10:
                continue
            if end < terminal_tm.end or (end - 1) < self.config.min_construct_length:
                continue
            location = self._classify_terminal_region(topology, removed_start, sequence_length)
            warnings: list[str] = []
            if location not in {"cytoplasmic", "membrane-adjacent", "unknown"}:
                warnings.append(
                    f"Caution: removed C-terminal region {removed_start}-{sequence_length} is classified as {location}."
                )
            for reason in self._terminal_trim_protection_reasons(
                removed_start,
                sequence_length,
                features=features,
                interpro_annotations=interpro_annotations,
                ligand_interactions=ligand_interactions,
            ):
                warnings.append(f"Caution: removed C-terminal region {reason}.")
            excluded_ptms = [
                feature
                for feature in ptms
                if feature.start <= sequence_length and feature.end >= removed_start
            ]
            if excluded_ptms:
                warnings.append(
                    f"Trimmed region excludes {len(excluded_ptms)} PTM annotation(s); review trafficking/regulatory biology."
                )
            variants.append(
                ConstructSuggestion(
                    name=f"membrane_expression_gpcr_c_tail_keep_{int(buffer_size)}",
                    start=1,
                    end=end,
                    score=84.0,
                    rationale=(
                        "GPCR native membrane-expression tail-series construct retaining the complete "
                        "7TM bundle while testing a more conservative C-terminal tail truncation."
                    ),
                    classification="membrane_expression",
                    evidence=[
                        f"Retains {anchor_label} plus {int(buffer_size)} residue(s).",
                        f"Removes C-terminal region {removed_start}-{sequence_length} ({removed_length} aa).",
                        "GPCR expression/stability screen variant; receptor sequence only, no foreign fusion.",
                    ],
                    warnings=self._unique_strings(warnings),
                )
            )
        return variants

    def _determine_multipass_c_trim_end(
        self,
        *,
        topology: TopologyAnnotation,
        sequence_length: int,
        residue_plddt: dict[int, float] | None,
        features: list[Feature],
        interpro_annotations: list[DomainAnnotation],
        ligand_interactions: list[Feature],
        ptms: list[Feature],
        gpcr_annotation: GPCRdbAnnotation | None = None,
    ) -> _TerminalTrimDecision:
        terminal_tm = max(topology.transmembrane_regions, key=lambda region: region.end)
        start = terminal_tm.end + 1
        end = sequence_length
        if start > end:
            return _TerminalTrimDecision()
        location = self._classify_terminal_region(topology, start, end)
        retained_anchor = terminal_tm.end
        coordinate_label = f"retaining {self.config.multipass_juxtamembrane_tail_keep} juxtamembrane residue(s) after TM{len(topology.transmembrane_regions)}"
        if gpcr_annotation is not None:
            h8 = self._gpcr_segment(gpcr_annotation, "H8")
            if h8 is not None and h8.end > retained_anchor:
                retained_anchor = h8.end
                coordinate_label = (
                    f"retaining GPCRdb helix 8 ({h8.start}-{h8.end}) plus "
                    f"{self.config.gpcr_helix8_retention_buffer} residue(s)."
                )
        buffer_size = (
            self.config.gpcr_helix8_retention_buffer
            if gpcr_annotation is not None and retained_anchor != terminal_tm.end
            else self.config.multipass_juxtamembrane_tail_keep
        )
        trim_end = min(sequence_length, retained_anchor + buffer_size)
        removed_start = trim_end + 1
        if removed_start > end:
            return _TerminalTrimDecision()
        return self._evaluate_terminal_trim(
            side="C-terminal",
            removed_start=removed_start,
            removed_end=end,
            retained_boundary=trim_end,
            location=location,
            residue_plddt=residue_plddt,
            features=features,
            interpro_annotations=interpro_annotations,
            ligand_interactions=ligand_interactions,
            ptms=ptms,
            allow_unknown_without_plddt=len(topology.transmembrane_regions) >= self.config.multipass_gpcr_like_tm_count,
            coordinate_label=coordinate_label,
            start=None,
            end=trim_end,
            plddt_threshold=(
                self.config.gpcr_c_tail_trim_plddt_threshold
                if gpcr_annotation is not None
                else self.config.multipass_terminal_disorder_plddt_threshold
            ),
            min_length=(
                self.config.gpcr_terminal_trim_min_length
                if gpcr_annotation is not None
                else self.config.multipass_terminal_trim_min_length
            ),
        )

    def _determine_multipass_n_trim_start(
        self,
        *,
        topology: TopologyAnnotation,
        sequence_length: int,
        residue_plddt: dict[int, float] | None,
        features: list[Feature],
        interpro_annotations: list[DomainAnnotation],
        ligand_interactions: list[Feature],
        ptms: list[Feature],
    ) -> _TerminalTrimDecision:
        first_tm = min(topology.transmembrane_regions, key=lambda region: region.start)
        start = 1
        end = first_tm.start - 1
        if start > end:
            return _TerminalTrimDecision()
        location = self._classify_terminal_region(topology, start, end)
        trim_start = max(1, first_tm.start - self.config.multipass_juxtamembrane_tail_keep)
        removed_end = trim_start - 1
        if removed_end < start:
            return _TerminalTrimDecision()
        return self._evaluate_terminal_trim(
            side="N-terminal",
            removed_start=start,
            removed_end=removed_end,
            retained_boundary=trim_start,
            location=location,
            residue_plddt=residue_plddt,
            features=features,
            interpro_annotations=interpro_annotations,
            ligand_interactions=ligand_interactions,
            ptms=ptms,
            allow_unknown_without_plddt=False,
            coordinate_label=f"retaining {self.config.multipass_juxtamembrane_tail_keep} juxtamembrane residue(s) before TM1",
            start=trim_start,
            end=None,
        )

    def _evaluate_terminal_trim(
        self,
        *,
        side: str,
        removed_start: int,
        removed_end: int,
        retained_boundary: int,
        location: str,
        residue_plddt: dict[int, float] | None,
        features: list[Feature],
        interpro_annotations: list[DomainAnnotation],
        ligand_interactions: list[Feature],
        ptms: list[Feature],
        allow_unknown_without_plddt: bool,
        coordinate_label: str,
        start: int | None,
        end: int | None,
        plddt_threshold: float | None = None,
        min_length: int | None = None,
    ) -> _TerminalTrimDecision:
        removed_length = removed_end - removed_start + 1
        min_length = min_length or self.config.multipass_terminal_trim_min_length
        plddt_threshold = plddt_threshold or self.config.multipass_terminal_disorder_plddt_threshold
        if removed_length < min_length:
            return _TerminalTrimDecision()
        protected = self._terminal_trim_protection_reasons(
            removed_start,
            removed_end,
            features=features,
            interpro_annotations=interpro_annotations,
            ligand_interactions=ligand_interactions,
        )
        ptm_count = sum(1 for feature in ptms if self._feature_overlaps_construct(feature, removed_start, removed_end))
        if side == "N-terminal" and ptm_count >= 3:
            protected.append(f"PTM-rich N-terminal region ({ptm_count} PTM annotation(s)).")
        if protected:
            return _TerminalTrimDecision(
                blocked_reason=(
                    f"Retained {side.lower()} region {removed_start}-{removed_end}; trimming was blocked because "
                    + "; ".join(protected[:4])
                )
            )
        if location != "cytoplasmic":
            return _TerminalTrimDecision(
                blocked_reason=(
                    f"Retained {side.lower()} region {removed_start}-{removed_end}; topology classifies it as "
                    f"{location}, so OpenAntigens does not trim it automatically."
                )
            )
        mean_plddt = self._mean_plddt(residue_plddt, removed_start, removed_end)
        low_confidence = (
            mean_plddt is not None
            and mean_plddt <= plddt_threshold
        )
        length_only_allowed = mean_plddt is None and (
            allow_unknown_without_plddt or removed_length >= max(40, min_length)
        )
        if not low_confidence and not length_only_allowed:
            if mean_plddt is not None:
                return _TerminalTrimDecision(
                    blocked_reason=(
                        f"Retained {side.lower()} cytoplasmic region {removed_start}-{removed_end}; mean pLDDT "
                        f"{mean_plddt:.1f} is above the terminal-disorder trim threshold."
                    )
                )
            return _TerminalTrimDecision()
        evidence = [
            "Native membrane-expression terminal-disorder optimization.",
            f"Trimmed {side.lower()} cytoplasmic region {removed_start}-{removed_end} ({removed_length} aa), {coordinate_label}.",
        ]
        if mean_plddt is not None:
            evidence.append(f"Trimmed-region mean pLDDT: {mean_plddt:.1f}.")
        else:
            evidence.append("Trimmed by topology/length heuristic because residue-level pLDDT was unavailable.")
        warnings = []
        if ptm_count:
            warnings.append(f"Trimmed region excludes {ptm_count} PTM annotation(s); review regulatory biology.")
        return _TerminalTrimDecision(
            enabled=True,
            start=start,
            end=end,
            evidence=evidence,
            warnings=warnings,
        )

    def _classify_terminal_region(self, topology: TopologyAnnotation, start: int, end: int) -> str:
        if self._regions_overlap_any(start, end, topology.transmembrane_regions + topology.intramembrane_regions):
            return "membrane"
        if self._regions_overlap_any(start, end, topology.extracellular_regions):
            return "extracellular"
        if self._regions_overlap_any(start, end, topology.cytoplasmic_regions):
            return "cytoplasmic"
        for tm in topology.transmembrane_regions:
            if abs(end - tm.start) <= self.config.multipass_juxtamembrane_tail_keep or abs(start - tm.end) <= self.config.multipass_juxtamembrane_tail_keep:
                return "membrane-adjacent"
        return "unknown"

    def _terminal_trim_protection_reasons(
        self,
        start: int,
        end: int,
        *,
        features: list[Feature],
        interpro_annotations: list[DomainAnnotation],
        ligand_interactions: list[Feature],
    ) -> list[str]:
        reasons: list[str] = []
        protected_feature_types = {"DOMAIN", "REPEAT", "MOTIF"}
        for annotation in interpro_annotations:
            if annotation.end >= start and annotation.start <= end:
                reasons.append(f"overlaps annotated domain {annotation.name} ({annotation.accession})")
                break
        for feature in features:
            feature_type = feature.type.upper()
            if feature_type == "CHAIN" and (
                start < feature.start < end or start < feature.end < end
            ):
                reasons.append(f"crosses mature-chain boundary {feature.start}-{feature.end}")
                break
            if feature_type in protected_feature_types and feature.start <= end and feature.end >= start:
                reasons.append(f"overlaps UniProt {feature_type} annotation {feature.start}-{feature.end}")
                break
            if feature_type == "DISULFID" and (start <= feature.start <= end or start <= feature.end <= end):
                reasons.append(f"contains disulfide endpoint {feature.start}-{feature.end}")
                break
        for feature in ligand_interactions:
            if feature.start <= end and feature.end >= start:
                reasons.append(f"overlaps ligand/interaction annotation {feature.start}-{feature.end}")
                break
        return self._unique_strings(reasons)

    def _regions_overlap_any(self, start: int, end: int, regions: list[Region]) -> bool:
        return any(region.start <= end and region.end >= start for region in regions)

    def _mean_plddt(self, residue_plddt: dict[int, float] | None, start: int, end: int) -> float | None:
        if residue_plddt is None:
            return None
        scores = [float(residue_plddt[position]) for position in range(start, end + 1) if position in residue_plddt]
        if not scores:
            return None
        return sum(scores) / len(scores)

    def _suggest_advanced_membrane_engineering(
        self,
        *,
        target,
        features: list[Feature],
        topology: TopologyAnnotation | None,
        interpro_annotations: list[DomainAnnotation],
        gpcr_annotation: GPCRdbAnnotation | None = None,
    ) -> list[MembraneEngineeringSuggestion]:
        if topology is None or not topology.topology_class.startswith("multipass"):
            return []
        tm_count = len(topology.transmembrane_regions)
        suggestions: list[MembraneEngineeringSuggestion] = []
        is_gpcr = gpcr_annotation is not None or tm_count >= self.config.multipass_gpcr_like_tm_count
        if is_gpcr:
            internal_cytoplasmic_loops = self._find_internal_cytoplasmic_loops(topology)
            longest_loop = max(
                internal_cytoplasmic_loops,
                key=lambda region: (region.end - region.start + 1, region.start),
                default=None,
            )
            tail = topology.cytoplasmic_tail_region
            gpcr_evidence = []
            if gpcr_annotation is not None:
                gpcr_evidence.append(
                    f"GPCRdb: {gpcr_annotation.entry_name}"
                    + (f"; {gpcr_annotation.family_name}" if gpcr_annotation.family_name else "")
                    + (f"; {gpcr_annotation.residue_numbering_scheme}" if gpcr_annotation.residue_numbering_scheme else "")
                )
            suggestions.append(
                MembraneEngineeringSuggestion(
                    category="gpcr_conservative_engineering",
                    title="Conservative 7TM expression series",
                    summary=(
                        "Treat this target as a 7TM/GPCR-like membrane protein and evaluate a small expression series "
                        "separate from the native construct recommendations."
                    ),
                    rationale=[
                        f"Topology contains {tm_count} transmembrane helices.",
                        "7TM receptors often benefit from separate native-like and expression-engineered construct screens.",
                    ],
                    suggested_actions=[
                        "Test the native full-length membrane construct alongside the trimmed membrane construct as the baseline pair.",
                        "Preserve the full 7TM core and a short juxtamembrane segment after the terminal helix before attempting more aggressive edits.",
                    ],
                    evidence=gpcr_evidence + [f"Topology-derived membrane core: {topology.membrane_core_region.start}-{topology.membrane_core_region.end}"],
                    warnings=["Advisory only; not part of the main native construct set."],
                    start=topology.membrane_core_region.start,
                    end=topology.membrane_core_region.end,
                )
            )
            if longest_loop is not None and (longest_loop.end - longest_loop.start + 1) >= 15:
                loop_length = longest_loop.end - longest_loop.start + 1
                suggestions.append(
                    MembraneEngineeringSuggestion(
                        category="gpcr_loop_engineering",
                        title="ICL-focused loop engineering screen",
                        summary="Consider an engineering series centered on the longest internal cytoplasmic loop.",
                        rationale=[
                            f"Longest internal cytoplasmic loop spans {longest_loop.start}-{longest_loop.end} ({loop_length} aa).",
                            "This topology-derived loop is a candidate insertion or truncation site for experimental testing.",
                        ],
                        suggested_actions=[
                            f"Create a conservative truncation variant preserving approximately 5-10 residues at each end of loop {longest_loop.start}-{longest_loop.end}.",
                            f"Create a fusion-insertion variant placing BRIL or T4 lysozyme into loop {longest_loop.start}-{longest_loop.end} rather than altering the main membrane core.",
                        ],
                        evidence=["Topology-derived internal cytoplasmic loop."],
                        warnings=[
                            "Advisory only; foreign-fusion constructs are not represented in the main construct tables.",
                            "Loop engineering may alter signaling or G-protein coupling and should be interpreted as an expression/stability optimization path.",
                        ],
                        start=longest_loop.start,
                        end=longest_loop.end,
                    )
                )
            if gpcr_annotation is not None:
                icl3 = self._gpcr_segment(gpcr_annotation, "ICL3")
                if icl3 is not None and (icl3.end - icl3.start + 1) >= 12:
                    suggestions.append(
                        MembraneEngineeringSuggestion(
                            category="gpcr_icl3_engineering",
                            title="GPCRdb ICL3 engineering checkpoint",
                            summary=(
                                "Use the GPCRdb-defined ICL3 boundaries for receptor-specific loop-engineering designs."
                            ),
                            rationale=[
                                f"GPCRdb maps ICL3 to {icl3.start}-{icl3.end}.",
                                "ICL3 is a common site for receptor stabilization screens, but edits should be treated as expression/stability engineering rather than antigen design.",
                            ],
                            suggested_actions=[
                                f"If a fusion screen is needed, center BRIL/T4L insertion design on ICL3 {icl3.start}-{icl3.end}.",
                                "Include native full-length and terminal-trimmed constructs as comparators.",
                            ],
                            evidence=["GPCRdb segment annotation."],
                            warnings=[
                                "Advisory only; foreign-fusion designs are not included in the precomputed construct table.",
                                "ICL3 engineering can alter signaling-state preference and intracellular-effector coupling.",
                            ],
                            start=icl3.start,
                            end=icl3.end,
                        )
                    )
            if tail is not None:
                suggestions.append(
                    MembraneEngineeringSuggestion(
                        category="gpcr_tail_engineering",
                        title="Terminal tail engineering options",
                        summary="Topology-derived terminal cytoplasmic-tail variants are candidates for experimental testing.",
                        rationale=[
                            f"Cytoplasmic tail spans {tail.start}-{tail.end}.",
                            "Tail edits leave the transmembrane bundle unchanged but can affect trafficking and signaling.",
                        ],
                        suggested_actions=[
                            "Compare terminal-tail variants with the native full-length construct.",
                            "If further stabilization is needed, test stepwise C-terminal truncations rather than immediate internal-loop edits.",
                        ],
                        evidence=["Topology-derived cytoplasmic tail."],
                        warnings=["Evaluate trafficking and signaling, because tail edits can alter regulation and internalization."],
                        start=tail.start,
                        end=tail.end,
                    )
                )

        if (
            (target.gene_symbol and str(target.gene_symbol).upper().startswith("ADGR"))
            or "adhesion g protein-coupled receptor" in str(target.protein_name).lower()
        ):
            gps_feature = self._find_feature_by_keyword(features, keyword="gps")
            gain_annotations = [
                annotation
                for annotation in interpro_annotations
                if "gain" in str(annotation.name or "").lower()
            ]
            start = gps_feature.start if gps_feature is not None else None
            end = gps_feature.end if gps_feature is not None else None
            evidence = []
            if gps_feature is not None:
                evidence.append(f"UniProt GPS feature: {gps_feature.start}-{gps_feature.end}")
            if gain_annotations:
                evidence.append(
                    "GAIN annotations: "
                    + ", ".join(f"{item.name} ({item.start}-{item.end})" for item in gain_annotations[:3])
                )
            suggestions.append(
                MembraneEngineeringSuggestion(
                    category="adhesion_gpcr_processing",
                    title="Adhesion GPCR cleavage-state series",
                    summary=(
                        "For adhesion GPCR-like targets, evaluate construct behavior across different GAIN/GPS processing states separately from the baseline membrane constructs."
                    ),
                    rationale=[
                        "Adhesion GPCRs often have large extracellular regions linked to a 7TM module through the GAIN/GPS region.",
                        "Proteolysis or local rearrangement around the GPS region can materially affect expression and biochemical behavior.",
                    ],
                    suggested_actions=[
                        "Track whether the baseline membrane constructs are cleaved near the GPS region.",
                        "If expression is poor, compare constructs that preserve the native GAIN/GPS neighborhood before attempting foreign-fusion engineering.",
                    ],
                    evidence=evidence,
                    warnings=["Treat GPS/GAIN engineering as a secondary optimization path after establishing the baseline membrane construct behavior."],
                    start=start,
                    end=end,
                )
            )
        return suggestions

    def _suggest_gpcr_engineering_variants(
        self,
        *,
        target,
        topology: TopologyAnnotation | None,
        features: list[Feature],
        ligand_interactions: list[Feature],
        ptms: list[Feature],
        gpcr_annotation: GPCRdbAnnotation | None = None,
    ) -> list[GPCREngineeringVariant]:
        if topology is None or not topology.topology_class.startswith("multipass"):
            return []
        tm_count = len(topology.transmembrane_regions)
        if gpcr_annotation is None and tm_count < self.config.multipass_gpcr_like_tm_count:
            return []
        registry = load_gpcr_engineering_registry(self.config.gpcr_engineering_cassette_registry_path)
        variants: list[GPCREngineeringVariant] = []
        loop, loop_source = self._select_gpcr_engineering_loop(topology, gpcr_annotation)
        loop_has_icl3_annotation = "ICL3" in loop_source
        loop_name = "ICL3" if loop_has_icl3_annotation else "ICL"
        loop_strategy = "ICL3 cassette replacement" if loop_has_icl3_annotation else "Intracellular-loop cassette replacement"
        blockers = self._gpcr_loop_engineering_blockers(
            loop=loop,
            topology=topology,
            features=features,
            ligand_interactions=ligand_interactions,
            ptms=ptms,
        )
        loop_length = (loop.end - loop.start + 1) if loop is not None else 0
        if loop is not None and not blockers:
            for cassette_id in ("bril", "t4l"):
                cassette, warning = registry.validated_for_mode(cassette_id, "loop_fusion")
                if warning:
                    variants.append(
                        self._gpcr_engineering_warning_variant(
                            target_entry_name=target.entry_name,
                            cassette_id=cassette_id,
                            warning=warning,
                            loop=loop,
                            loop_source=loop_source,
                        )
                    )
                    continue
                assert cassette is not None
                retain = self.config.gpcr_loop_fusion_retain_residues
                replaced_start = loop.start + retain
                replaced_end = loop.end - retain
                if replaced_start > replaced_end:
                    variants.append(
                        self._gpcr_engineering_warning_variant(
                            target_entry_name=target.entry_name,
                            cassette_id=cassette_id,
                            warning=(
                                f"Selected GPCR loop {loop.start}-{loop.end} is too short for "
                                f"{retain}-residue flank retention."
                            ),
                            loop=loop,
                            loop_source=loop_source,
                        )
                    )
                    continue
                engineered_sequence = (
                    target.sequence[: replaced_start - 1]
                    + cassette.default_n_linker
                    + str(cassette.sequence)
                    + cassette.default_c_linker
                    + target.sequence[replaced_end:]
                )
                cassette_label = cassette.cassette_id.upper().replace("_", "-")
                variants.append(
                    GPCREngineeringVariant(
                        name=f"{target.entry_name}_{cassette_label}_{loop_name}_fusion",
                        category="gpcr_loop_fusion",
                        strategy=loop_strategy,
                        summary=(
                            f"Experimental {cassette.name} fusion replacing the central part of "
                            f"{loop_source} {loop.start}-{loop.end} while retaining {retain} receptor "
                            "residues on each side."
                        ),
                        loop_source=loop_source,
                        sequence=engineered_sequence,
                        sequence_role="engineered_receptor",
                        cassette_id=cassette.cassette_id,
                        cassette_name=cassette.name,
                        cassette_sequence=cassette.sequence,
                        start=loop.start,
                        end=loop.end,
                        replaced_start=replaced_start,
                        replaced_end=replaced_end,
                        insertion_after=replaced_start - 1,
                        retained_n_terminal_loop_residues=retain,
                        retained_c_terminal_loop_residues=retain,
                        n_linker=cassette.default_n_linker,
                        c_linker=cassette.default_c_linker,
                        native_length=len(target.sequence),
                        engineered_length=len(engineered_sequence),
                        gpcr_segments=self._gpcr_segments_in_region(gpcr_annotation, loop.start, loop.end),
                        gpcr_generic_range=self._gpcr_generic_range(gpcr_annotation, loop.start, loop.end),
                        tag_suggestions=self._gpcr_tag_suggestions(),
                        suggested_actions=[
                            "Treat as an expression/stability rescue construct, not as the primary antigen construct.",
                            "Compare against native full-length and terminal-trimmed membrane-expression constructs.",
                            "Validate ligand binding and cell-surface expression before antibody-campaign use.",
                        ],
                        evidence=[
                            f"{loop_source}: {loop.start}-{loop.end} ({loop_length} aa).",
                            f"Cassette registry: {cassette.name}; {cassette.version}; {cassette.source_note}",
                        ],
                        warnings=[
                            "Experimental GPCR engineering variant; may alter receptor conformational ensemble.",
                            f"Receptor residues {replaced_start}-{replaced_end} are replaced by a foreign cassette.",
                            "Intracellular loop fusions can alter G-protein/effector coupling and receptor-state preference.",
                        ],
                    )
                )
        elif loop is not None:
            variants.append(
                GPCREngineeringVariant(
                    name=f"{target.entry_name}_{loop_name}_engineering_skipped",
                    category="gpcr_loop_fusion_skipped",
                    strategy=loop_strategy,
                    summary=f"Experimental loop-fusion variants were not generated for {loop_source} {loop.start}-{loop.end}.",
                    loop_source=loop_source,
                    sequence=None,
                    sequence_role="warning",
                    start=loop.start,
                    end=loop.end,
                    gpcr_segments=self._gpcr_segments_in_region(gpcr_annotation, loop.start, loop.end),
                    gpcr_generic_range=self._gpcr_generic_range(gpcr_annotation, loop.start, loop.end),
                    evidence=[f"{loop_source}: {loop.start}-{loop.end} ({loop_length} aa)."],
                    warnings=blockers,
                )
            )
        else:
            variants.append(
                GPCREngineeringVariant(
                    name=f"{target.entry_name}_ICL_engineering_skipped",
                    category="gpcr_loop_fusion_skipped",
                    strategy="Intracellular-loop cassette replacement",
                    summary="Experimental loop-fusion variants were not generated because no suitable intracellular loop was identified.",
                    loop_source=loop_source,
                    sequence=None,
                    sequence_role="warning",
                    warnings=["No GPCRdb ICL3 or topology-derived internal cytoplasmic loop was available."],
                )
            )

        return variants

    def _select_gpcr_engineering_loop(
        self,
        topology: TopologyAnnotation,
        gpcr_annotation: GPCRdbAnnotation | None,
    ) -> tuple[Region | None, str]:
        internal_loops = self._find_internal_cytoplasmic_loops(topology)
        if gpcr_annotation is not None:
            icl3 = self._gpcr_segment(gpcr_annotation, "ICL3")
            if icl3 is not None and (icl3.end - icl3.start + 1) >= self.config.gpcr_loop_fusion_min_loop_length:
                return Region(start=icl3.start, end=icl3.end, label="ICL3", source="GPCRdb"), "GPCRdb ICL3"
            if icl3 is not None and internal_loops:
                overlapping = [
                    loop
                    for loop in internal_loops
                    if loop.start <= icl3.end and loop.end >= icl3.start
                ]
                if overlapping:
                    loop = max(overlapping, key=lambda region: (region.end - region.start + 1, -region.start))
                    return loop, "topology-derived ICL3 loop anchored by GPCRdb"
        if not internal_loops:
            return None, "topology-derived internal cytoplasmic loop"
        loop = max(internal_loops, key=lambda region: (region.end - region.start + 1, -region.start))
        return loop, "topology-derived internal cytoplasmic loop"

    def _gpcr_loop_engineering_blockers(
        self,
        *,
        loop: Region | None,
        topology: TopologyAnnotation,
        features: list[Feature],
        ligand_interactions: list[Feature],
        ptms: list[Feature],
    ) -> list[str]:
        if loop is None:
            return ["No GPCRdb ICL3 or topology-derived internal cytoplasmic loop was available."]
        blockers: list[str] = []
        loop_length = loop.end - loop.start + 1
        if loop_length < self.config.gpcr_loop_fusion_min_loop_length:
            blockers.append(
                f"Selected loop is only {loop_length} aa; minimum is {self.config.gpcr_loop_fusion_min_loop_length} aa."
            )
        if self._regions_overlap_any(loop.start, loop.end, topology.transmembrane_regions + topology.intramembrane_regions):
            blockers.append("Selected loop overlaps a membrane-spanning region.")
        if self._regions_overlap_any(loop.start, loop.end, topology.extracellular_regions):
            blockers.append("Selected loop overlaps extracellular topology annotation.")
        for feature in ligand_interactions:
            feature_kind = str(feature.metadata.get("raw_type") or feature.type or "").upper()
            feature_text = f"{feature.type} {feature.description or ''}".lower()
            is_ligand_like = feature_kind in {"BINDING", "SITE"} or "ligand" in feature_text
            if is_ligand_like and feature.start <= loop.end and feature.end >= loop.start:
                blockers.append(f"Selected loop overlaps ligand-interaction annotation {feature.start}-{feature.end}.")
                break
        for feature in features + ptms:
            feature_type = str(feature.type or "").upper()
            if feature_type == "DISULFID" and (loop.start <= feature.start <= loop.end or loop.start <= feature.end <= loop.end):
                blockers.append(f"Selected loop contains disulfide endpoint {feature.start}-{feature.end}.")
                break
            if feature_type in {"CARBOHYD", "LIPID", "MOD_RES", "CROSSLNK", "SITE"} and feature.start <= loop.end and feature.end >= loop.start:
                blockers.append(f"Selected loop overlaps {feature_type} annotation {feature.start}-{feature.end}.")
                break
        return self._unique_strings(blockers)

    def _gpcr_engineering_warning_variant(
        self,
        *,
        target_entry_name: str,
        cassette_id: str,
        warning: str,
        loop: Region,
        loop_source: str,
    ) -> GPCREngineeringVariant:
        loop_name = "ICL3" if "ICL3" in loop_source else "ICL"
        loop_strategy = "ICL3 cassette replacement" if "ICL3" in loop_source else "Intracellular-loop cassette replacement"
        return GPCREngineeringVariant(
            name=f"{target_entry_name}_{cassette_id}_{loop_name}_fusion_unavailable",
            category="gpcr_loop_fusion_skipped",
            strategy=loop_strategy,
            summary=f"{cassette_id} loop-fusion variant was not generated for {loop_source} {loop.start}-{loop.end}.",
            loop_source=loop_source,
            sequence=None,
            sequence_role="warning",
            cassette_id=cassette_id,
            start=loop.start,
            end=loop.end,
            evidence=[f"{loop_source}: {loop.start}-{loop.end}."],
            warnings=[warning],
        )

    def _gpcr_tag_suggestions(self) -> dict[str, Any]:
        return {
            "n_terminal": [
                "Optional secretion/trafficking signal or epitope tag should be chosen in the expression vector, not appended by OpenAntigens.",
                "For antibody campaigns, avoid masking extracellular epitopes with large N-terminal tags unless experimentally required.",
            ],
            "c_terminal": [
                "Optional detection/purification tags can be placed after retained C-terminal receptor sequence.",
                "C-terminal tags may affect trafficking or internalization for some GPCRs.",
            ],
            "notes": [
                "OpenAntigens reports tag suggestions as metadata only; generated receptor sequences are untagged."
            ],
        }

    def _find_internal_cytoplasmic_loops(self, topology: TopologyAnnotation) -> list[Region]:
        if topology.membrane_core_region is None:
            return []
        core = topology.membrane_core_region
        tail = topology.cytoplasmic_tail_region
        loops: list[Region] = []
        for region in topology.cytoplasmic_regions:
            if region.end < core.start or region.start > core.end:
                continue
            if tail is not None and region.start == tail.start and region.end == tail.end:
                continue
            loops.append(region)
        return loops

    def _gpcr_segment(self, annotation: GPCRdbAnnotation, name: str) -> GPCRdbSegment | None:
        normalized = name.strip().lower()
        for segment in annotation.segments:
            if str(segment.name or "").strip().lower() == normalized:
                return segment
        return None

    def _find_feature_by_keyword(self, features: list[Feature], *, keyword: str) -> Feature | None:
        lowered = keyword.lower()
        for feature in features:
            text = " ".join(
                [
                    str(feature.type or ""),
                    str(feature.description or ""),
                ]
            ).lower()
            if lowered in text:
                return feature
        return None

    def _annotate_structural_regions(
        self,
        structural_regions: list[Region],
        *,
        interpro_annotations: list[DomainAnnotation],
        residue_plddt: dict[int, float],
        pae_stats,
    ) -> list[Region]:
        annotated: list[Region] = []
        domain_annotations = [
            annotation
            for annotation in interpro_annotations
            if annotation.type == "domain"
        ]
        for region in structural_regions:
            metrics = summarize_construct_quality(
                region.start,
                region.end,
                residue_plddt,
                pae_stats,
                self.config,
            )
            mean_intra_pae = metrics.mean_intra_pae
            confidence = region.confidence
            if confidence is None and mean_intra_pae is not None:
                confidence = max(
                    0.4,
                    min(0.95, 1.0 - (mean_intra_pae / max(self.config.pae_domain_threshold * 2.0, 1.0))),
                )
                confidence = round(confidence, 2)
            included_domain_labels = [
                self._short_annotation_label(annotation)
                for annotation in domain_annotations
                if annotation.start >= region.start and annotation.end <= region.end
            ]
            annotated.append(
                replace(
                    region,
                    confidence=confidence,
                    metadata={
                        **region.metadata,
                        "mean_intra_pae": mean_intra_pae,
                        "included_interpro_domains": included_domain_labels,
                    },
                )
            )
        return annotated

    def _collect_assembly_requirements(
        self,
        *,
        entry: dict[str, Any],
        gene_symbol: str | None,
        family_context,
        complex_portal_complexes: list[ComplexPortalComplex],
        target_accession: str,
    ) -> list[AssemblyRequirement]:
        if not self._is_canonical_integrin_chain(gene_symbol):
            return []
        requirements: list[AssemblyRequirement] = []
        subunit_texts = self.uniprot_client.get_comment_texts(entry, "SUBUNIT")
        requirements.extend(self._requirements_from_subunit_comments(subunit_texts, gene_symbol=gene_symbol))
        requirements.extend(
            self._requirements_from_family_rules(
                gene_symbol,
                family_context,
                subunit_texts,
                complex_portal_complexes=complex_portal_complexes,
                target_accession=target_accession,
            )
        )
        requirements.extend(
            self._requirements_from_complex_portal(
                complex_portal_complexes,
                gene_symbol=gene_symbol,
                target_accession=target_accession,
            )
        )
        return self._dedupe_assembly_requirements(requirements)

    def _requirements_from_subunit_comments(
        self,
        comments: list[str],
        *,
        gene_symbol: str | None,
    ) -> list[AssemblyRequirement]:
        requirements: list[AssemblyRequirement] = []
        for comment in comments:
            for sentence in self._split_subunit_comment(comment):
                lowered = sentence.lower()
                if not any(
                    token in lowered
                    for token in (
                        "heterodimer",
                        "heteromer",
                        "homodimer",
                        "interacts with",
                        "forms a complex",
                        "complex with",
                        "associates with",
                    )
                ):
                    continue
                partners = self._extract_partner_symbols(sentence, target_gene_symbol=gene_symbol)
                classification, obligatory, confidence, summary = self._classify_subunit_comment(sentence)
                if obligatory:
                    obligatory = False
                    classification = "interaction_context"
                    summary = (
                        "Curated UniProt subunit annotation provides interaction or assembly context, "
                        "but Complex Portal evidence is required before calling an obligatory partner requirement."
                    )
                requirements.append(
                    AssemblyRequirement(
                        summary=summary,
                        obligatory=obligatory,
                        confidence=confidence,
                        source="UniProt SUBUNIT",
                        classification=classification,
                        partners=partners,
                        rationale=[sentence],
                    )
                )
        return requirements

    def _split_subunit_comment(self, comment: str) -> list[str]:
        text = re.sub(r"\s+", " ", comment).strip()
        if not text:
            return []
        sentences = [
            sentence.strip()
            for sentence in re.split(r"(?<=[.;])\s+", text)
            if sentence.strip()
        ]
        return sentences or [text]

    def _classify_subunit_comment(self, text: str) -> tuple[str, bool, str, str]:
        lowered = text.lower()
        explicit_obligate_tokens = (
            "obligate",
            "obligatory",
            "requires ",
            "required for cell surface expression",
            "required for surface expression",
            "disulfide-linked heterodimer",
            "heterodimer of an alpha and a beta subunit",
            "stable heterodimer with",
        )
        conditional_tokens = (
            "binding of the ligand triggers",
            "ligand triggers",
            "ligand-induced",
            "upon ligand",
            "after ligand stimulation",
            "homo- and/or heterodimerization",
        )
        generic_interaction_tokens = (
            "interacts with",
            "forms a complex with",
            "part of a complex with",
            "complex with",
        )
        if any(token in lowered for token in explicit_obligate_tokens):
            return (
                "obligatory_partner_requirement",
                True,
                "high",
                "Curated UniProt subunit annotation indicates a likely obligatory assembly partner requirement.",
            )
        if any(token in lowered for token in conditional_tokens):
            return (
                "conditional_assembly",
                False,
                "medium",
                "Curated UniProt subunit annotation indicates conditional or ligand-triggered dimerization, not an obligatory partner requirement by itself.",
            )
        if any(token in lowered for token in generic_interaction_tokens):
            return (
                "interaction_context",
                False,
                "medium",
                "Curated UniProt subunit annotation lists interaction partners or complexes, but this is not sufficient evidence for an obligatory partner requirement.",
            )
        if "heterodimer" in lowered or "heteromer" in lowered:
            return (
                "heteromeric_assembly",
                False,
                "medium",
                "Curated UniProt subunit annotation indicates heteromeric assembly, but the text does not establish that the partner is obligatory.",
            )
        return (
            "interaction_context",
            False,
            "low",
            "Curated UniProt subunit annotation mentions assembly-related context without establishing an obligatory partner requirement.",
        )

    def _requirements_from_family_rules(
        self,
        gene_symbol: str | None,
        family_context,
        subunit_texts: list[str],
        *,
        complex_portal_complexes: list[ComplexPortalComplex],
        target_accession: str,
    ) -> list[AssemblyRequirement]:
        if not gene_symbol:
            return []
        requirements: list[AssemblyRequirement] = []
        combined_text = " ".join(subunit_texts)
        partners = self._extract_partner_symbols(combined_text)
        partners.extend(
            partner for partner in self._infer_integrin_partner_symbols(combined_text) if partner not in partners
        )
        complex_portal_partners = self._collect_complex_portal_partner_symbols(
            complex_portal_complexes,
            target_accession=target_accession,
        )
        if self._is_canonical_integrin_chain(gene_symbol):
            if gene_symbol.startswith("ITGA"):
                beta_partners = [partner for partner in partners if partner.startswith("ITGB")]
                beta_partners.extend(
                    partner for partner in complex_portal_partners if partner.startswith("ITGB") and partner not in beta_partners
                )
                complex_ids = self._complex_portal_ids_for_target(complex_portal_complexes, target_accession)
                requirements.append(
                    AssemblyRequirement(
                        summary="Integrin alpha chains require an integrin beta partner for the biologically relevant cell-surface heterodimer.",
                        obligatory=True,
                        confidence="high",
                        source="Family rule + UniProt" + (" + Complex Portal" if complex_ids else ""),
                        classification="obligatory_partner_requirement",
                        partners=beta_partners or ["integrin beta chain"],
                        complex_portal_ids=complex_ids,
                        rationale=[
                            f"Target gene symbol {gene_symbol} matches the integrin alpha-chain family.",
                            "Integrins are functional as alpha/beta heterodimers at the cell surface.",
                        ]
                        + ([subunit_texts[0]] if subunit_texts else []),
                    )
                )
            elif gene_symbol.startswith("ITGB"):
                alpha_partners = [partner for partner in partners if partner.startswith("ITGA")]
                alpha_partners.extend(
                    partner for partner in complex_portal_partners if partner.startswith("ITGA") and partner not in alpha_partners
                )
                complex_ids = self._complex_portal_ids_for_target(complex_portal_complexes, target_accession)
                requirements.append(
                    AssemblyRequirement(
                        summary="Integrin beta chains require an integrin alpha partner for the biologically relevant cell-surface heterodimer.",
                        obligatory=True,
                        confidence="high",
                        source="Family rule + UniProt" + (" + Complex Portal" if complex_ids else ""),
                        classification="obligatory_partner_requirement",
                        partners=alpha_partners or ["integrin alpha chain"],
                        complex_portal_ids=complex_ids,
                        rationale=[
                            f"Target gene symbol {gene_symbol} matches the integrin beta-chain family.",
                            "Integrins are functional as alpha/beta heterodimers at the cell surface.",
                        ]
                        + ([subunit_texts[0]] if subunit_texts else []),
                    )
                )
        return requirements

    def _infer_integrin_partner_symbols(self, text: str) -> list[str]:
        partners: list[str] = []
        for match in re.finditer(r"\bbeta[-\s]?([0-9]+|[a-z])\b", text, flags=re.IGNORECASE):
            suffix = match.group(1).upper()
            symbol = f"ITGB{suffix}"
            if symbol not in partners:
                partners.append(symbol)
        for match in re.finditer(r"\balpha[-\s]?([0-9]+)\b", text, flags=re.IGNORECASE):
            suffix = match.group(1).upper()
            symbol = f"ITGA{suffix}"
            if symbol not in partners:
                partners.append(symbol)
        return partners

    def _complex_portal_ids_for_target(
        self,
        complexes: list[ComplexPortalComplex],
        target_accession: str,
    ) -> list[str]:
        return [
            complex_item.complex_ac
            for complex_item in complexes
            if any(
                self._complex_participant_matches_target(
                    participant,
                    target_accession=target_accession,
                    gene_symbol=None,
                )
                for participant in complex_item.participants
            )
            and any(
                not self._complex_participant_matches_target(
                    participant,
                    target_accession=target_accession,
                    gene_symbol=None,
                )
                for participant in complex_item.participants
            )
        ]

    def _requirements_from_complex_portal(
        self,
        complexes: list[ComplexPortalComplex],
        *,
        gene_symbol: str | None,
        target_accession: str,
    ) -> list[AssemblyRequirement]:
        requirements: list[AssemblyRequirement] = []
        human_requirements: list[AssemblyRequirement] = []
        mouse_support: list[AssemblyRequirement] = []
        for complex_item in complexes:
            if complex_item.predicted_complex:
                continue
            target_participants, partner_participants = self._partition_complex_participants(
                complex_item,
                target_accession=target_accession,
                gene_symbol=gene_symbol,
            )
            if not target_participants or not partner_participants:
                continue
            partner_symbols = [
                self._complex_participant_symbol(participant)
                for participant in partner_participants
                if self._complex_participant_symbol(participant)
            ]
            score = complex_item.confidence_score or 0
            taxon_id = self._complex_portal_taxon_id(complex_item)
            is_mouse = taxon_id == 10090
            if is_mouse:
                mouse_support.append(
                    AssemblyRequirement(
                        summary="Mouse Complex Portal complex membership provides ortholog support/context but does not trigger a human obligatory-partner warning by itself.",
                        obligatory=False,
                        confidence="high" if score >= 4 else "medium",
                        source="Complex Portal (mouse)",
                        classification="mouse_supported_complex_context",
                        partners=partner_symbols,
                        complex_portal_ids=[complex_item.complex_ac],
                        rationale=[
                            f"{complex_item.name} ({complex_item.complex_ac})",
                            f"Species: {complex_item.species or 'n/a'}",
                            f"Assembly: {', '.join(complex_item.complex_assemblies) or 'n/a'}",
                            f"Evidence: {complex_item.evidence_code or 'n/a'} {complex_item.evidence_description or ''}".strip(),
                        ],
                    )
                )
                continue
            if taxon_id not in (None, 9606):
                continue

            if self._complex_portal_supports_obligatory_call(
                complex_item,
                gene_symbol=gene_symbol,
                partner_symbols=partner_symbols,
            ):
                human_requirements.append(
                    AssemblyRequirement(
                        summary="Human Complex Portal evidence supports an obligatory partner relationship for the biologically relevant complex.",
                        obligatory=True,
                        confidence="high" if score >= 4 else "medium",
                        source="Complex Portal",
                        classification="obligatory_partner_requirement",
                        partners=partner_symbols,
                        complex_portal_ids=[complex_item.complex_ac],
                        rationale=[
                            f"{complex_item.name} ({complex_item.complex_ac})",
                            f"Species: {complex_item.species or 'n/a'}",
                            f"Assembly: {', '.join(complex_item.complex_assemblies) or 'n/a'}",
                            f"Evidence: {complex_item.evidence_code or 'n/a'} {complex_item.evidence_description or ''}".strip(),
                        ],
                    )
                )
            else:
                human_requirements.append(
                    AssemblyRequirement(
                        summary="Human Complex Portal complex membership provides stable complex context but does not by itself prove an obligatory partner requirement.",
                        obligatory=False,
                        confidence="high" if score >= 4 else "medium",
                        source="Complex Portal",
                        classification="stable_complex_context",
                        partners=partner_symbols,
                        complex_portal_ids=[complex_item.complex_ac],
                        rationale=[
                            f"{complex_item.name} ({complex_item.complex_ac})",
                            f"Species: {complex_item.species or 'n/a'}",
                            f"Assembly: {', '.join(complex_item.complex_assemblies) or 'n/a'}",
                            f"Evidence: {complex_item.evidence_code or 'n/a'} {complex_item.evidence_description or ''}".strip(),
                        ],
                    )
                )
        requirements.extend(self._attach_mouse_complex_support(human_requirements, mouse_support))
        human_partner_keys = {self._partner_class_key(requirement.partners) for requirement in human_requirements}
        requirements.extend(
            requirement
            for requirement in mouse_support
            if self._partner_class_key(requirement.partners) not in human_partner_keys
        )
        return requirements

    def _complex_portal_taxon_id(self, complex_item: ComplexPortalComplex) -> int | None:
        species_text = str(complex_item.species or "")
        for taxon_id in (9606, 10090):
            if str(taxon_id) in species_text:
                return taxon_id
        return None

    def _complex_portal_supports_obligatory_call(
        self,
        complex_item: ComplexPortalComplex,
        *,
        gene_symbol: str | None,
        partner_symbols: list[str],
    ) -> bool:
        if complex_item.predicted_complex or not self._is_canonical_integrin_chain(gene_symbol):
            return False
        return bool(gene_symbol and self._is_integrin_alpha_beta_pair(gene_symbol, partner_symbols))

    def _is_integrin_alpha_beta_pair(self, gene_symbol: str, partner_symbols: list[str]) -> bool:
        symbol = gene_symbol.upper()
        if symbol in _CANONICAL_INTEGRIN_ALPHA_SYMBOLS and any(
            partner in _CANONICAL_INTEGRIN_BETA_SYMBOLS for partner in partner_symbols
        ):
            return True
        if symbol in _CANONICAL_INTEGRIN_BETA_SYMBOLS and any(
            partner in _CANONICAL_INTEGRIN_ALPHA_SYMBOLS for partner in partner_symbols
        ):
            return True
        return False

    @staticmethod
    def _canonical_integrin_beta_symbols() -> set[str]:
        return set(_CANONICAL_INTEGRIN_BETA_SYMBOLS)

    def _is_canonical_integrin_chain(self, gene_symbol: str | None) -> bool:
        symbol = str(gene_symbol or "").upper()
        return symbol in _CANONICAL_INTEGRIN_SYMBOLS

    def _attach_mouse_complex_support(
        self,
        human_requirements: list[AssemblyRequirement],
        mouse_support: list[AssemblyRequirement],
    ) -> list[AssemblyRequirement]:
        attached: list[AssemblyRequirement] = []
        for human_requirement in human_requirements:
            matches = [
                requirement
                for requirement in mouse_support
                if self._partner_class_key(requirement.partners) == self._partner_class_key(human_requirement.partners)
            ]
            if not matches:
                attached.append(human_requirement)
                continue
            rationale = list(human_requirement.rationale)
            complex_portal_ids = list(human_requirement.complex_portal_ids)
            for match in matches:
                complex_portal_ids.extend(
                    complex_id for complex_id in match.complex_portal_ids if complex_id not in complex_portal_ids
                )
                rationale.append(f"Mouse Complex Portal support: {match.rationale[0] if match.rationale else match.summary}")
            attached.append(
                replace(
                    human_requirement,
                    source=f"{human_requirement.source}; Complex Portal (mouse support)",
                    rationale=rationale,
                    complex_portal_ids=complex_portal_ids,
                )
            )
        return attached

    def _partner_class_key(self, partners: list[str]) -> tuple[str, ...]:
        classes: list[str] = []
        for partner in partners:
            if partner.startswith("ITGA"):
                item = "ITGA"
            elif partner.startswith("ITGB"):
                item = "ITGB"
            else:
                item = partner
            if item not in classes:
                classes.append(item)
        return tuple(sorted(classes))

    def _partition_complex_participants(
        self,
        complex_item: ComplexPortalComplex,
        *,
        target_accession: str,
        gene_symbol: str | None,
    ) -> tuple[list[Any], list[Any]]:
        target_participants = []
        partner_participants = []
        for participant in complex_item.participants:
            if self._complex_participant_matches_target(participant, target_accession=target_accession, gene_symbol=gene_symbol):
                target_participants.append(participant)
            elif (participant.interactor_type or "").lower() == "protein":
                partner_participants.append(participant)
        return target_participants, partner_participants

    def _complex_participant_matches_target(self, participant, *, target_accession: str, gene_symbol: str | None) -> bool:
        identifier = str(participant.identifier or "").split("-", 1)[0].upper()
        name = str(participant.name or "").upper()
        if identifier == target_accession.upper():
            return True
        if gene_symbol and (name == gene_symbol.upper() or identifier == gene_symbol.upper()):
            return True
        return False

    def _complex_participant_symbol(self, participant) -> str | None:
        name = str(participant.name or "").upper()
        identifier = str(participant.identifier or "").split("-", 1)[0].upper()
        integrin_symbol = _CANONICAL_INTEGRIN_SYMBOL_BY_ACCESSION.get(identifier)
        if integrin_symbol is not None:
            return integrin_symbol
        if re.fullmatch(r"[A-Z0-9-]{2,12}", name):
            return name
        if re.fullmatch(r"(?:ITGA[0-9A-Z]+|ITGB[0-9A-Z]+|[A-Z][A-Z0-9]{2,9})", identifier):
            return identifier
        return name or None

    def _collect_complex_portal_partner_symbols(
        self,
        complexes: list[ComplexPortalComplex],
        *,
        target_accession: str,
    ) -> list[str]:
        symbols: list[str] = []
        for complex_item in complexes:
            for participant in complex_item.participants:
                identifier = str(participant.identifier or "").split("-", 1)[0].upper()
                if identifier == target_accession.upper():
                    continue
                symbol = self._complex_participant_symbol(participant)
                if symbol and symbol not in symbols:
                    symbols.append(symbol)
        return symbols

    def _extract_partner_symbols(self, text: str, target_gene_symbol: str | None = None) -> list[str]:
        symbols = re.findall(r"\b(?:ITGA[0-9A-Z]+|ITGB[0-9A-Z]+|[A-Z][A-Z0-9]{2,9})\b", text or "")
        ignored = {
            "DNA",
            "RNA",
            "EGF",
            "ECM",
            "ATP",
            "GTP",
            "SUBUNIT",
        }
        unique: list[str] = []
        for symbol in symbols:
            if symbol in ignored:
                continue
            if target_gene_symbol and symbol == target_gene_symbol:
                continue
            if symbol not in unique:
                unique.append(symbol)
        return unique

    def _dedupe_assembly_requirements(
        self,
        requirements: list[AssemblyRequirement],
    ) -> list[AssemblyRequirement]:
        deduped: dict[tuple[str, tuple[str, ...], bool], AssemblyRequirement] = {}
        for requirement in requirements:
            key = (requirement.classification or "", tuple(sorted(requirement.partners)), requirement.obligatory)
            current = deduped.get(key)
            if current is None:
                deduped[key] = requirement
                continue
            rationale = current.rationale + [item for item in requirement.rationale if item not in current.rationale]
            partners = current.partners + [item for item in requirement.partners if item not in current.partners]
            confidence = current.confidence
            if requirement.confidence == "high" or current.confidence != "high":
                confidence = requirement.confidence
            summary = current.summary
            if requirement.source == "Complex Portal" or (
                requirement.source.startswith("Complex Portal")
                and not current.source.startswith("Complex Portal")
            ):
                summary = requirement.summary
            elif (
                "Family rule" in requirement.source
                and "Family rule" not in current.source
                and "Complex Portal" not in current.source
            ):
                summary = requirement.summary
            sources = current.source.split("; ")
            for source in requirement.source.split("; "):
                if source == "Complex Portal" and any("Complex Portal" in item for item in sources):
                    continue
                if source not in sources:
                    sources.append(source)
            deduped[key] = replace(
                current,
                summary=summary,
                obligatory=current.obligatory or requirement.obligatory,
                confidence=confidence,
                source="; ".join(sources),
                classification=current.classification or requirement.classification,
                partners=partners,
                rationale=rationale,
                complex_portal_ids=current.complex_portal_ids
                + [item for item in requirement.complex_portal_ids if item not in current.complex_portal_ids],
            )
        return sorted(
            deduped.values(),
            key=lambda item: (
                0 if item.obligatory else 1,
                0 if item.confidence == "high" else 1,
                item.summary,
            ),
        )

    def _find_furin_sites(self, sequence: str, ectodomain: Region | None) -> list[FurinSite]:
        if not ectodomain:
            return []
        ectodomain_sequence = slice_sequence(sequence, ectodomain.start, ectodomain.end)
        sites = []
        for start, end, motif in find_furin_sites(ectodomain_sequence, self.config.furin_pattern):
            absolute_start = ectodomain.start + start - 1
            absolute_end = ectodomain.start + end - 1
            suggestions = [
                f"{motif[0]}{absolute_start}A",
                f"{motif[-1]}{absolute_end}A",
            ]
            sites.append(FurinSite(start=absolute_start, end=absolute_end, motif=motif, suggestions=suggestions))
        return sites

    def _collect_ptm_annotations(self, features: list[Feature]) -> list[Feature]:
        ptms: list[Feature] = []
        for feature in features:
            if self._is_ptm_feature(feature):
                metadata = dict(feature.metadata)
                metadata["ptm_category"] = self._ptm_category(feature)
                ptms.append(
                    Feature(
                        type=feature.type,
                        start=feature.start,
                        end=feature.end,
                        description=feature.description,
                        source=feature.source,
                        metadata=metadata,
                    )
                )
        return sorted(ptms, key=lambda item: (item.start, item.end, item.type))

    def _is_ptm_feature(self, feature: Feature) -> bool:
        feature_type = feature.type.upper()
        description = (feature.description or "").lower()
        if feature_type in _PTM_FEATURE_TYPES:
            return True
        if feature_type == "SITE" and any(
            keyword in description
            for keyword in ("cleavage", "glycosyl", "phospho", "acetyl", "methyl", "ubiquitin", "sumoyl", "lipid")
        ):
            return True
        return False

    def _ptm_category(self, feature: Feature) -> str:
        feature_type = feature.type.upper()
        description = (feature.description or "").lower()
        if feature_type == "CARBOHYD":
            return "glycosylation"
        if feature_type == "MOD_RES":
            return "modified residue"
        if feature_type == "LIPID":
            return "lipidation"
        if feature_type == "DISULFID":
            return "disulfide"
        if feature_type == "CROSSLNK":
            return "cross-link"
        if feature_type == "INIT_MET":
            return "initiator methionine processing"
        if feature_type in {"SIGNAL", "PROPEP", "PEPTIDE", "CHAIN"}:
            return "molecule processing"
        if "cleavage" in description:
            return "cleavage site"
        return "post-translational modification"

    def _ptm_warning(self, feature: Feature) -> str | None:
        feature_type = feature.type.upper()
        category = str(feature.metadata.get("ptm_category") or self._ptm_category(feature)).lower()
        description = (feature.description or "").lower()
        if feature_type in _PROCESSING_FEATURE_TYPES or "cleavage" in category or "cleavage" in description:
            return f"Contains processing/cleavage annotation {feature.start}-{feature.end}: {feature.description or feature.type}."
        return None

    def _build_construct_details(
        self,
        *,
        target_entry_name: str,
        target_sequence: str,
        constructs: list[ConstructSuggestion],
        homolog_context: dict[str, dict[str, Any]],
        ectodomain: Region | None,
        ligand_interactions: list[Feature],
        interpro_annotations: list[DomainAnnotation],
        ptms: list[Feature] | None = None,
        furin_sites: list[FurinSite] | None = None,
        structure_atoms: dict[int, dict[str, tuple[float, float, float]]] | None = None,
        residue_names: dict[int, str] | None = None,
        residue_plddt: dict[int, float] | None = None,
        surface_exposure: dict[int, bool | None] | None = None,
        residue_annotations: list[ResidueAnnotation] | None = None,
        pae_matrix: list[list[float]] | None = None,
        pae_stats=None,
        gpcr_annotation: GPCRdbAnnotation | None = None,
    ) -> list[ConstructDetail]:
        details: list[ConstructDetail] = []
        ptms = ptms or []
        furin_sites = furin_sites or []
        design_annotations = self._filter_design_annotations(interpro_annotations, ectodomain)
        structural_metrics_cache: dict[tuple[int, int], Any] = {}
        split_diagnostics_cache: dict[tuple[int, int], list[Any]] = {}
        for construct in constructs:
            homologs = self._map_construct_to_homologs(
                construct=construct,
                target_sequence=target_sequence,
                homolog_context=homolog_context,
                ectodomain=ectodomain,
                extracellular_surface_positions=self._extracellular_surface_positions(residue_annotations or []),
            )
            derived_classification, domain_summary, complete_annotations, clipped_annotations = self._classify_construct(
                construct.start,
                construct.end,
                design_annotations,
            )
            classification = construct.classification or derived_classification
            construct_warnings = list(construct.warnings)
            construct_cysteines: list[CysteineFinding] = []
            structural_metrics = None
            split_diagnostics = []
            construct_ligand_interactions = [
                feature
                for feature in ligand_interactions
                if feature.start <= construct.end and feature.end >= construct.start
            ]
            construct_ptms = [
                feature
                for feature in ptms
                if self._feature_overlaps_construct(feature, construct.start, construct.end)
            ]
            gpcr_segments = self._gpcr_segments_in_region(gpcr_annotation, construct.start, construct.end)
            gpcr_generic_range = self._gpcr_generic_range(gpcr_annotation, construct.start, construct.end)
            for feature in construct_ptms:
                warning = self._ptm_warning(feature)
                if warning and warning not in construct_warnings:
                    construct_warnings.append(warning)
            construct_furin_sites = [
                site
                for site in furin_sites
                if site.start >= construct.start and site.end <= construct.end
            ]
            for site in construct_furin_sites:
                construct_warnings.append(
                    f"Contains candidate furin cleavage motif {site.motif} at {site.start}-{site.end}."
                )
            if residue_plddt is not None and pae_stats is not None:
                cache_key = (construct.start, construct.end)
                structural_metrics = structural_metrics_cache.get(cache_key)
                if structural_metrics is None:
                    structural_metrics = summarize_construct_quality(
                        construct.start,
                        construct.end,
                        residue_plddt,
                        pae_stats,
                        self.config,
                    )
                    structural_metrics_cache[cache_key] = structural_metrics
                if cache_key not in split_diagnostics_cache:
                    split_diagnostics = collect_region_split_diagnostics(
                        Region(
                            start=construct.start,
                            end=construct.end,
                            label=construct.name,
                            source="Construct",
                        ),
                        pae_stats,
                        self.config,
                        residue_plddt=residue_plddt,
                    )
                    split_diagnostics_cache[cache_key] = split_diagnostics
                split_diagnostics = split_diagnostics_cache[cache_key]
            if structure_atoms is not None and residue_names is not None:
                construct_cysteines = analyze_cysteines(
                    structure_atoms,
                    residue_names,
                    construct.start,
                    construct.end,
                    surface_exposure=surface_exposure,
                )
                unpaired = [finding.position for finding in construct_cysteines if finding.paired_with is None]
                if unpaired:
                    positions = ", ".join(
                        self._unpaired_cysteine_warning_label(finding)
                        for finding in construct_cysteines
                        if finding.paired_with is None
                    )
                    construct_warnings.append(f"Contains unpaired cysteine(s): {positions}.")
            details.append(
                ConstructDetail(
                    name=construct.name,
                    start=construct.start,
                    end=construct.end,
                    length=construct.end - construct.start + 1,
                    sequence=slice_sequence(target_sequence, construct.start, construct.end),
                    score=construct.score,
                    rationale=construct.rationale,
                    classification=classification,
                    warnings=construct_warnings,
                    evidence=list(construct.evidence),
                    pdb_id=construct.pdb_id,
                    human_entry_name=target_entry_name,
                    domain_summary=domain_summary,
                    complete_domain_annotations=complete_annotations,
                    clipped_domain_annotations=clipped_annotations,
                    homologs=homologs,
                    cysteine_analysis=construct_cysteines,
                    ptms=construct_ptms,
                    ligand_interactions=construct_ligand_interactions,
                    gpcr_segments=gpcr_segments,
                    gpcr_generic_range=gpcr_generic_range,
                    structural_metrics=structural_metrics,
                    split_diagnostics=split_diagnostics,
                )
            )
        return details

    def _unpaired_cysteine_warning_label(self, finding: CysteineFinding) -> str:
        label = f"C{finding.position}"
        if finding.closest_unpaired_cysteine is None or finding.closest_unpaired_cysteine_distance is None:
            return label
        atom = finding.closest_unpaired_cysteine_distance_atom or "distance"
        return (
            f"{label} (nearest unpaired C{finding.closest_unpaired_cysteine}, "
            f"{finding.closest_unpaired_cysteine_distance:.2f} A {atom})"
        )

    def _gpcr_segments_in_region(
        self,
        annotation: GPCRdbAnnotation | None,
        start: int,
        end: int,
    ) -> list[GPCRdbSegment]:
        if annotation is None:
            return []
        return [
            segment
            for segment in annotation.segments
            if segment.start <= end and segment.end >= start
        ]

    def _gpcr_generic_range(
        self,
        annotation: GPCRdbAnnotation | None,
        start: int,
        end: int,
    ) -> str | None:
        if annotation is None:
            return None
        generics = [
            residue.display_generic_number
            for residue in annotation.residues
            if start <= residue.sequence_number <= end and residue.display_generic_number
        ]
        if not generics:
            return None
        return f"{generics[0]} to {generics[-1]}"

    def _feature_overlaps_construct(self, feature: Feature, start: int, end: int) -> bool:
        if feature.type.upper() == "DISULFID":
            # UniProt disulfide features encode two bonded cysteine positions,
            # not a continuous residue interval between them.
            return feature.start in range(start, end + 1) or feature.end in range(start, end + 1)
        return feature.start <= end and feature.end >= start

    def _canonicalize_macaca_cross_reactivity_hits(
        self,
        query_sequence: str,
        hits: list,
    ) -> list:
        macaca_hits = [hit for hit in hits if (hit.species or "").lower() == "macaca fascicularis"]
        if not macaca_hits:
            return hits

        best_by_gene: dict[str, Any] = {}
        for hit in macaca_hits:
            gene_symbol = self._extract_gene_symbol_from_hit(hit)
            if not gene_symbol:
                continue
            current = best_by_gene.get(gene_symbol)
            if current is None or (hit.bitscore, hit.coverage, hit.identity) > (
                current.bitscore,
                current.coverage,
                current.identity,
            ):
                best_by_gene[gene_symbol] = hit

        canonical_hits = []
        for gene_symbol, seed_hit in best_by_gene.items():
            try:
                if gene_symbol not in self._macaca_refseq_by_gene:
                    self._macaca_refseq_by_gene[gene_symbol] = (
                        self.refseq_client.fetch_canonical_protein(
                            gene_symbol=gene_symbol,
                            organism="Macaca fascicularis",
                        )
                    )
                result = self._macaca_refseq_by_gene[gene_symbol]
            except Exception:
                result = None
            if result is None:
                continue
            alignment = global_align(query_sequence, result.target.sequence)
            query_start, query_end, subject_start, subject_end = self._alignment_bounds(alignment)
            canonical_hits.append(
                replace(
                    seed_hit,
                    subject_id=result.target.accession,
                    description=f"GN={gene_symbol} {result.title or result.target.protein_name}",
                    species=result.target.organism,
                    identity=round(alignment.identity, 3),
                    coverage=round(alignment.coverage, 3),
                    alignment_length=alignment.aligned_positions,
                    query_start=query_start,
                    query_end=query_end,
                    subject_start=subject_start,
                    subject_end=subject_end,
                    query_alignment=alignment.aligned_query,
                    subject_alignment=alignment.aligned_subject,
                    alignment_source="Canonical remap",
                )
            )

        non_macaca_hits = [hit for hit in hits if (hit.species or "").lower() != "macaca fascicularis"]
        combined = non_macaca_hits + canonical_hits
        return sorted(combined, key=lambda hit: (-hit.bitscore, hit.evalue))

    def _extract_gene_symbol_from_hit(self, hit) -> str | None:
        for text in (hit.description, hit.subject_id):
            if not text:
                continue
            match = re.search(r"\bGN=([A-Za-z0-9_-]+)\b", text)
            if match:
                return match.group(1)
        return None

    def _alignment_bounds(self, alignment) -> tuple[int, int, int, int]:
        query_position = 0
        subject_position = 0
        query_positions: list[int] = []
        subject_positions: list[int] = []
        for query_residue, subject_residue in zip(
            alignment.aligned_query,
            alignment.aligned_subject,
            strict=True,
        ):
            if query_residue != "-":
                query_position += 1
            if subject_residue != "-":
                subject_position += 1
            if query_residue != "-" and subject_residue != "-":
                query_positions.append(query_position)
                subject_positions.append(subject_position)
        if not query_positions or not subject_positions:
            return 1, len(alignment.aligned_query.replace("-", "")), 1, len(alignment.aligned_subject.replace("-", ""))
        return min(query_positions), max(query_positions), min(subject_positions), max(subject_positions)

    def _collect_ligand_interactions(
        self,
        features: list[Feature],
        ectodomain: Region | None,
    ) -> list[Feature]:
        interactions: list[Feature] = []
        for feature in features:
            if ectodomain and (feature.end < ectodomain.start or feature.start > ectodomain.end):
                continue
            if not self._is_ligand_interaction_feature(feature):
                continue
            interactions.append(
                replace(
                    feature,
                    metadata={
                        **feature.metadata,
                        "interaction_summary": self._describe_ligand_interaction(feature),
                    },
                )
            )
        return interactions

    def _is_ligand_interaction_feature(self, feature: Feature) -> bool:
        feature_type = feature.type.upper()
        if feature_type in {"BINDING", "BINDING_SITE"}:
            return True
        if "ligand" in feature.metadata:
            return True
        if feature_type not in {"SITE", "REGION", "DOMAIN", "MUTAGENESIS"}:
            return False
        text = " ".join(
            part
            for part in [
                feature.description or "",
                str(feature.metadata.get("raw_type", "")),
            ]
            if part
        ).lower()
        keywords = (
            "ligand",
            "binding",
            "growth factor",
            "binds",
            "interaction with",
        )
        return any(keyword in text for keyword in keywords)

    def _describe_ligand_interaction(self, feature: Feature) -> str:
        ligand = feature.metadata.get("ligand")
        if isinstance(ligand, dict):
            name = ligand.get("name") or ligand.get("label")
            if name:
                return str(name)
        if isinstance(ligand, str) and ligand.strip():
            return ligand.strip()
        if feature.description:
            return feature.description
        return feature.type.replace("_", " ").title()

    def _map_construct_to_homologs(
        self,
        *,
        construct: ConstructSuggestion,
        target_sequence: str,
        homolog_context: dict[str, dict[str, Any]],
        ectodomain: Region | None,
        extracellular_surface_positions: set[int] | None = None,
    ) -> list[ConstructHomolog]:
        homologs: list[ConstructHomolog] = []
        human_construct_sequence = slice_sequence(target_sequence, construct.start, construct.end)
        for species in ("mouse", "macaca_fascicularis"):
            info = homolog_context.get(species)
            if not info or not info.get("available"):
                notes = info.get("notes", ["No same-name UniProt entry found for species."]) if info else [
                    "No same-name UniProt entry found for species."
                ]
                homologs.append(
                    ConstructHomolog(
                        species=species,
                        accession=None,
                        entry_name=None,
                        available=False,
                        notes=list(notes),
                    )
                )
                continue

            target = info["target"]
            homolog_ectodomain = info.get("ectodomain")
            alignment = info.get("alignment")
            homolog_ecto_sequence = info.get("ectodomain_sequence")
            notes = list(info.get("notes", []))
            if not homolog_ectodomain or not alignment or not homolog_ecto_sequence:
                homologs.append(
                    ConstructHomolog(
                        species=species,
                        accession=target.accession,
                        entry_name=target.entry_name,
                        available=False,
                        notes=notes or ["Homolog ectodomain mapping unavailable."],
                    )
                )
                continue

            query_region_start = ectodomain.start if ectodomain is not None else 1
            relative_start = construct.start - query_region_start + 1
            relative_end = construct.end - query_region_start + 1
            mapped_start, mapped_end, mapped_sequence, mapping_notes = map_query_region_to_subject(
                alignment,
                query_start=relative_start,
                query_end=relative_end,
                subject_sequence=homolog_ecto_sequence,
            )
            if mapped_start is None or mapped_end is None or mapped_sequence is None:
                homologs.append(
                    ConstructHomolog(
                        species=species,
                        accession=target.accession,
                        entry_name=target.entry_name,
                        available=False,
                        notes=notes + mapping_notes,
                    )
                )
                continue
            surface_identity = None
            surface_aligned = None
            if extracellular_surface_positions:
                subset_positions = {
                    position
                    for position in extracellular_surface_positions
                    if construct.start <= position <= construct.end
                }
                if subset_positions:
                    surface_result = identity_for_query_positions(
                        alignment,
                        query_region_start=query_region_start,
                        query_positions=subset_positions,
                    )
                    surface_identity = (
                        round(surface_result.identity, 2)
                        if surface_result.identity is not None
                        else None
                    )
                    surface_aligned = surface_result.aligned_positions
            homologs.append(
                ConstructHomolog(
                    species=species,
                    accession=target.accession,
                    entry_name=target.entry_name,
                    available=True,
                    start=homolog_ectodomain.start + mapped_start - 1,
                    end=homolog_ectodomain.start + mapped_end - 1,
                    sequence=mapped_sequence,
                    identity_to_human=round(
                        global_align(human_construct_sequence, mapped_sequence).identity,
                        2,
                    ),
                    extracellular_surface_identity_to_human=surface_identity,
                    extracellular_surface_aligned_positions=surface_aligned,
                    notes=notes + mapping_notes,
                )
            )
        return homologs

    def _generate_construct_assets(
        self,
        *,
        construct_details: list[ConstructDetail],
        assets_dir: Path,
        alphafold_pdb_path: Path,
        residue_plddt: dict[int, float],
        pae_matrix: list[list[float]],
        ectodomain: Region,
        residue_annotations: list[ResidueAnnotation] | None = None,
    ) -> None:
        total = len(construct_details)
        eligible_constructs = [
            construct
            for construct in construct_details
            if construct.start >= ectodomain.start and construct.end <= ectodomain.end
        ]
        if not eligible_constructs:
            return
        for construct in eligible_constructs:
            if not self.config.render_structure_images:
                construct.structure_image = None
            if not self.config.render_quality_plots:
                construct.quality_plot = None
        self._progress(f"Rendering assets for {len(eligible_constructs)}/{total} construct(s)")
        asset_results = self.asset_generator.generate_assets_for_constructs(
            constructs=eligible_constructs,
            pdb_path=alphafold_pdb_path,
            residue_plddt=residue_plddt,
            pae_matrix=pae_matrix,
            ectodomain=ectodomain,
            residue_annotations=residue_annotations,
            output_dir=assets_dir,
            render_structure_image=self.config.render_structure_images,
            render_quality_plot=self.config.render_quality_plots,
        )
        for construct in eligible_constructs:
            structure_image, quality_plot = asset_results.get(construct.name, (None, None))
            if structure_image is not None:
                construct.structure_image = f"{assets_dir.name}/{structure_image}"
            if quality_plot is not None:
                construct.quality_plot = f"{assets_dir.name}/{quality_plot}"

    def _dedupe_regions(self, regions: list[Region]) -> list[Region]:
        seen: set[tuple[int, int, str]] = set()
        deduped: list[Region] = []
        for region in sorted(regions, key=lambda item: (item.start, item.end, item.label)):
            key = (region.start, region.end, region.label)
            if key in seen:
                continue
            seen.add(key)
            deduped.append(region)
        return deduped

    def _dedupe_constructs(self, constructs: list[ConstructSuggestion]) -> list[ConstructSuggestion]:
        merged: dict[tuple[int, int] | tuple[int, int, str], ConstructSuggestion] = {}
        for construct in constructs:
            key: tuple[int, int] | tuple[int, int, str]
            # PDB-backed constructs are retained individually, even when their
            # boundaries duplicate another construct, because each deposited
            # structure is independent experimental evidence. Keying on the
            # unique construct name keeps every PDB reference on the report.
            if construct.pdb_id or construct.name == "membrane_expression_full_length":
                key = (construct.start, construct.end, construct.name)
            else:
                key = (construct.start, construct.end)
            current = merged.get(key)
            if current is None or construct.score > current.score:
                merged[key] = construct
        return list(merged.values())


def analyze_target(
    query: str,
    config: AnalysisConfig | None = None,
    output_dir: str | Path | None = None,
) -> AnalysisReport:
    analyzer = AntigenAnalyzer(config=config)
    return analyzer.analyze_target(query, output_dir=output_dir)
