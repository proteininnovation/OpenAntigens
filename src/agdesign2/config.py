from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path


@dataclass(slots=True)
class AnalysisConfig:
    target_species: str = "human"
    cache_dir: Path = Path(".agdesign2/cache")
    data_dir: Path = Path(".agdesign2/data")
    fetch_alphafold: bool = True
    enable_local_af3_fallback: bool = False
    blast_db_dir: Path = Path(".agdesign2/data/blastdb")
    ortholog_fasta_dir: Path = Path(".agdesign2/data/proteomes")
    ortholog_blast_db_dir: Path = Path(".agdesign2/data/ortholog_blastdb")
    precomputed_ortholog_table_path: Path = Path("outputs/surfy_batch/ortholog_reference_table.tsv")
    precomputed_family_alignment_index_path: Path = Path("outputs/surfy_batch/family_alignments/family_alignment_index.tsv")
    precomputed_paralog_index_path: Path = Path("outputs/surfy_batch/paralogs/paralog_index.tsv")
    prefer_precomputed_references: bool = True
    plddt_structured_threshold: float = 70.0
    pae_domain_threshold: float = 12.0
    pae_separation_threshold: float = 4.0
    max_disordered_gap: int = 8
    lenient_plddt_structured_threshold: float = 60.0
    lenient_max_disordered_gap: int = 12
    min_structured_segment: int = 25
    min_domain_size: int = 40
    min_construct_length: int = 50
    domain_linker_window: int = 5
    max_domain_recursion_depth: int = 4
    furin_pattern: str = r"R.[KR]R"
    blast_evalue: float = 1e-5
    blast_max_hits: int = 25
    blast_threads: int = field(default_factory=lambda: max(1, min(8, os.cpu_count() or 1)))
    blast_species: tuple[str, ...] = ("human", "mouse", "macaca_fascicularis")
    ortholog_blast_species: tuple[str, ...] = ("macaca_fascicularis",)
    cross_reactivity_fallback_species: tuple[str, ...] = ("macaca_fascicularis",)
    run_cross_reactivity: bool = True
    generate_assets: bool = True
    render_structure_images: bool = True
    render_quality_plots: bool = True
    verbose_progress: bool = False
    enable_complex_portal: bool = False
    complex_portal_timeout: int = 8
    complex_portal_max_pages: int = 1
    enable_gpcrdb: bool = True
    gpcrdb_cache_dir: Path | None = None
    gpcr_engineering_cassette_registry_path: Path | None = None
    gpcr_c_tail_trim_plddt_threshold: float = 55.0
    gpcr_terminal_trim_min_length: int = 30
    gpcr_helix8_retention_buffer: int = 8
    gpcr_c_tail_series_buffers: tuple[int, ...] = (20, 40)
    gpcr_loop_fusion_retain_residues: int = 5
    gpcr_loop_fusion_min_loop_length: int = 12
    max_family_context_members: int = 20
    max_paralog_context_members: int = 20
    multipass_tail_trim_min_length: int = 20
    multipass_terminal_trim_min_length: int = 30
    multipass_terminal_disorder_plddt_threshold: float = 55.0
    multipass_gpcr_like_tm_count: int = 7
    multipass_juxtamembrane_tail_keep: int = 12
    multipass_mixed_ectodomain_min_length: int = 80
    sasa_probe_radius: float = 1.4
    sasa_slices: int = 20
    sasa_relative_exposed_threshold: float = 0.20
    sasa_low_confidence_plddt_threshold: float = 50.0
    sasa_target_accessibility_only: bool = True
    sasa_skip_low_confidence_residues: bool = True
    topology_keywords_extracellular: tuple[str, ...] = (
        "extracellular",
        "lumenal",
        "luminal",
        "exoplasmic",
    )
    topology_keywords_cytoplasmic: tuple[str, ...] = ("cytoplasmic", "cytosolic")
    species_suffixes: dict[str, str] = field(
        default_factory=lambda: {
            "human": "HUMAN",
            "mouse": "MOUSE",
            "macaca_fascicularis": "MACFA",
        }
    )
    species_taxonomy: dict[str, int] = field(
        default_factory=lambda: {
            "human": 9606,
            "mouse": 10090,
            "macaca_fascicularis": 9541,
        }
    )

    def ensure_directories(self) -> None:
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.data_dir.mkdir(parents=True, exist_ok=True)
        self.blast_db_dir.mkdir(parents=True, exist_ok=True)
        self.ortholog_fasta_dir.mkdir(parents=True, exist_ok=True)
        self.ortholog_blast_db_dir.mkdir(parents=True, exist_ok=True)
