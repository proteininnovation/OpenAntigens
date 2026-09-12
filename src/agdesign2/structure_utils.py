from __future__ import annotations

from dataclasses import dataclass
import json
import gzip
import math
from pathlib import Path
from typing import Any

import freesasa

from .config import AnalysisConfig
from .models import (
    AnalysisNote,
    ConstructSplitDiagnostic,
    ConstructStructuralMetrics,
    CysteineFinding,
    Region,
)

freesasa.setVerbosity(freesasa.nowarnings)

_THREE_TO_ONE = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "PYL": "O",
    "SER": "S",
    "SEC": "U",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "ASX": "B",
    "GLX": "Z",
    "UNK": "X",
}


@dataclass(slots=True)
class PaeStats:
    matrix: list[list[float]]
    prefix_sum: list[list[float]]
    diagonal_prefix: list[float]


@dataclass(slots=True)
class SurfaceAccessibility:
    position: int
    full_model_sasa: float | None
    full_model_relative_sasa: float | None
    domain_context_sasa: float | None
    domain_context_relative_sasa: float | None
    classification: str
    reason: str
    pae_block_id: str | None = None

    @property
    def surface_exposed(self) -> bool | None:
        if self.classification in {"exposed", "conditional"}:
            return True
        if self.classification == "buried":
            return False
        return None


_RESIDUE_MAX_ASA = {
    "ALA": 129.0,
    "ARG": 274.0,
    "ASN": 195.0,
    "ASP": 193.0,
    "CYS": 167.0,
    "GLN": 225.0,
    "GLU": 223.0,
    "GLY": 104.0,
    "HIS": 224.0,
    "ILE": 197.0,
    "LEU": 201.0,
    "LYS": 236.0,
    "MET": 224.0,
    "PHE": 240.0,
    "PRO": 159.0,
    "SER": 155.0,
    "THR": 172.0,
    "TRP": 285.0,
    "TYR": 263.0,
    "VAL": 174.0,
}


def parse_alphafold_pdb(path: Path) -> dict[str, Any]:
    residue_plddt: dict[int, list[float]] = {}
    residue_atoms: dict[int, dict[str, tuple[float, float, float]]] = {}
    residue_names: dict[int, str] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            atom_name = line[12:16].strip()
            residue_name = line[17:20].strip()
            residue_id = int(line[22:26].strip())
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            bfactor = float(line[60:66].strip())
            residue_plddt.setdefault(residue_id, []).append(bfactor)
            residue_atoms.setdefault(residue_id, {})[atom_name] = (x, y, z)
            residue_names[residue_id] = residue_name
    mean_plddt = {
        residue_id: sum(values) / len(values) for residue_id, values in residue_plddt.items()
    }
    return {"plddt": mean_plddt, "atoms": residue_atoms, "residue_names": residue_names}


def residue_names_to_sequence(residue_names: dict[int, str]) -> str:
    return "".join(
        _THREE_TO_ONE.get(str(residue_names[residue_id]).upper(), "X")
        for residue_id in sorted(residue_names)
    )


def load_pae_matrix(path: Path) -> list[list[float]]:
    if path.name.endswith(".gz"):
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            payload = json.load(handle)
    else:
        payload = json.loads(path.read_text(encoding="utf-8"))
    if isinstance(payload, list) and payload and isinstance(payload[0], dict):
        payload = payload[0]
    matrix = payload.get("predicted_aligned_error") or payload.get("pae")
    if matrix is None:
        raise ValueError("PAE JSON does not contain a predicted aligned error matrix.")
    return matrix


def pae_matches_sequence_length(pae_matrix: Any, sequence_length: int) -> bool:
    return (
        sequence_length > 0
        and isinstance(pae_matrix, list)
        and len(pae_matrix) == sequence_length
        and all(isinstance(row, list) and len(row) == sequence_length for row in pae_matrix)
    )


def build_pae_stats(pae_matrix: list[list[float]]) -> PaeStats:
    size = len(pae_matrix)
    if not pae_matches_sequence_length(pae_matrix, size):
        raise ValueError("PAE must be a nonempty square matrix")
    prefix_sum = [[0.0] * (size + 1) for _ in range(size + 1)]
    diagonal_prefix = [0.0] * (size + 1)
    for row in range(size):
        row_running_sum = 0.0
        for col in range(size):
            value = float(pae_matrix[row][col])
            if not math.isfinite(value) or value < 0:
                raise ValueError("PAE values must be finite and nonnegative")
            row_running_sum += value
            prefix_sum[row + 1][col + 1] = prefix_sum[row][col + 1] + row_running_sum
        diagonal_prefix[row + 1] = diagonal_prefix[row] + float(pae_matrix[row][row])
    return PaeStats(matrix=pae_matrix, prefix_sum=prefix_sum, diagonal_prefix=diagonal_prefix)


def find_structured_regions(
    ectodomain_start: int,
    ectodomain_end: int,
    residue_plddt: dict[int, float],
    config: AnalysisConfig,
    *,
    plddt_threshold: float | None = None,
    max_disordered_gap: int | None = None,
    label: str = "Structured segment",
    metadata: dict[str, Any] | None = None,
) -> list[Region]:
    threshold = config.plddt_structured_threshold if plddt_threshold is None else plddt_threshold
    gap_limit = config.max_disordered_gap if max_disordered_gap is None else max_disordered_gap
    region_metadata = dict(metadata or {})
    residues = []
    for residue in range(ectodomain_start, ectodomain_end + 1):
        if residue_plddt.get(residue, 0.0) >= threshold:
            residues.append(residue)
    if not residues:
        return []

    regions: list[Region] = []
    start = residues[0]
    previous = residues[0]
    for residue in residues[1:]:
        if residue - previous <= gap_limit + 1:
            previous = residue
            continue
        if previous - start + 1 >= config.min_structured_segment:
            regions.append(
                Region(
                    start=start,
                    end=previous,
                    label=label,
                    source="AlphaFold",
                    metadata=dict(region_metadata),
                )
            )
        start = residue
        previous = residue

    if previous - start + 1 >= config.min_structured_segment:
        regions.append(
            Region(
                start=start,
                end=previous,
                label=label,
                source="AlphaFold",
                metadata=dict(region_metadata),
            )
        )
    return regions


def split_domains_from_pae(
    region: Region,
    pae_matrix: list[list[float]] | PaeStats,
    config: AnalysisConfig,
    residue_plddt: dict[int, float] | None = None,
) -> list[Region]:
    pae_stats = _coerce_pae_stats(pae_matrix)
    return _split_region(region.start, region.end, pae_stats, config, residue_plddt, depth=0)


def analyze_cysteines(
    atoms: dict[int, dict[str, tuple[float, float, float]]],
    residue_names: dict[int, str],
    construct_start: int,
    construct_end: int,
    surface_exposure: dict[int, bool | None] | None = None,
) -> list[CysteineFinding]:
    cysteines = [
        residue
        for residue, name in residue_names.items()
        if construct_start <= residue <= construct_end and name == "CYS"
    ]
    findings: list[CysteineFinding] = []
    paired: dict[int, int] = {}
    candidates: dict[int, list[int]] = {residue: [] for residue in cysteines}
    for i, residue_a in enumerate(cysteines):
        for residue_b in cysteines[i + 1 :]:
            sg_a = atoms.get(residue_a, {}).get("SG")
            sg_b = atoms.get(residue_b, {}).get("SG")
            if not sg_a or not sg_b:
                continue
            if _distance(sg_a, sg_b) <= 2.4:
                candidates[residue_a].append(residue_b)
                candidates[residue_b].append(residue_a)
    for residue, partners in candidates.items():
        if len(partners) == 1 and len(candidates[partners[0]]) == 1:
            paired[residue] = partners[0]

    unpaired_cysteines = [residue for residue in cysteines if residue not in paired]
    closest_unpaired = {
        residue: _closest_unpaired_cysteine(residue, unpaired_cysteines, atoms)
        for residue in unpaired_cysteines
    }

    for residue in cysteines:
        pair = paired.get(residue)
        exposed = (
            surface_exposure.get(residue)
            if surface_exposure is not None and residue in surface_exposure
            else _is_surface_exposed(residue, atoms)
        )
        warning = None
        if pair is None and candidates[residue]:
            warning = "Ambiguous cysteine pairing: competing sulfur contacts."
        elif pair is None and exposed:
            warning = "Unpaired cysteine appears surface exposed."
        closest = closest_unpaired.get(residue)
        findings.append(
            CysteineFinding(
                position=residue,
                paired_with=pair,
                surface_exposed=exposed,
                warning=warning,
                closest_unpaired_cysteine=closest[0] if closest else None,
                closest_unpaired_cysteine_distance=round(closest[1], 2) if closest else None,
                closest_unpaired_cysteine_distance_atom=closest[2] if closest else None,
            )
        )
    return findings


def _closest_unpaired_cysteine(
    residue: int,
    unpaired_cysteines: list[int],
    atoms: dict[int, dict[str, tuple[float, float, float]]],
) -> tuple[int, float, str] | None:
    best: tuple[int, float, str] | None = None
    for other in unpaired_cysteines:
        if other == residue:
            continue
        distance_info = _cysteine_distance(residue, other, atoms)
        if distance_info is None:
            continue
        distance, atom_name = distance_info
        if best is None or distance < best[1]:
            best = (other, distance, atom_name)
    return best


def _cysteine_distance(
    residue_a: int,
    residue_b: int,
    atoms: dict[int, dict[str, tuple[float, float, float]]],
) -> tuple[float, str] | None:
    atom_map_a = atoms.get(residue_a, {})
    atom_map_b = atoms.get(residue_b, {})
    sg_a = atom_map_a.get("SG")
    sg_b = atom_map_b.get("SG")
    if sg_a and sg_b:
        return _distance(sg_a, sg_b), "SG-SG"
    ca_a = atom_map_a.get("CA")
    ca_b = atom_map_b.get("CA")
    if ca_a and ca_b:
        return _distance(ca_a, ca_b), "CA-CA"
    return None


def estimate_surface_exposure(
    atoms: dict[int, dict[str, tuple[float, float, float]]],
    *,
    residues: set[int] | None = None,
    radius: float = 10.0,
    neighbor_threshold: int = 14,
) -> dict[int, bool | None]:
    """Approximate residue exposure from local C-alpha/C-beta packing density.

    This deliberately avoids heavyweight DSSP/SASA dependencies. It is a
    deterministic proxy: low neighbor count within `radius` is treated as
    surface-exposed, while high neighbor count is treated as buried.
    """

    selected = residues or set(atoms)
    centers: dict[int, tuple[float, float, float]] = {}
    for residue in selected:
        atom_map = atoms.get(residue, {})
        center = atom_map.get("CA") or atom_map.get("CB")
        if center:
            centers[residue] = center

    exposure: dict[int, bool | None] = {residue: None for residue in selected}
    if not centers:
        return exposure

    cell_size = radius
    grid: dict[tuple[int, int, int], list[int]] = {}
    for residue, center in centers.items():
        cell = _spatial_cell(center, cell_size)
        grid.setdefault(cell, []).append(residue)

    radius_squared = radius * radius
    for residue, center in centers.items():
        neighbors = 0
        cell = _spatial_cell(center, cell_size)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for other in grid.get((cell[0] + dx, cell[1] + dy, cell[2] + dz), []):
                        if other == residue:
                            continue
                        if _distance_squared(center, centers[other]) <= radius_squared:
                            neighbors += 1
        exposure[residue] = neighbors < neighbor_threshold
    return exposure


def classify_surface_accessibility(
    atoms: dict[int, dict[str, tuple[float, float, float]]],
    residue_names: dict[int, str],
    residue_plddt: dict[int, float],
    blocks: list[Region],
    config: AnalysisConfig,
    *,
    target_residues: set[int] | None = None,
) -> dict[int, SurfaceAccessibility]:
    residues = set(residue_names)
    if target_residues is not None:
        residues &= set(target_residues)
    if not residues:
        return {}
    low_confidence_residues = {
        residue
        for residue in residues
        if residue_plddt.get(residue) is not None
        and residue_plddt[residue] < config.sasa_low_confidence_plddt_threshold
    }
    full_sasa_residues = residues
    if config.sasa_skip_low_confidence_residues:
        full_sasa_residues = residues - low_confidence_residues
    full_sasa = calculate_residue_sasa(
        atoms,
        residue_names,
        residues=full_sasa_residues,
        context_residues=set(residue_names),
        probe_radius=config.sasa_probe_radius,
        slices=config.sasa_slices,
    ) if full_sasa_residues else {}
    needs_block_context: set[int] = set()
    for residue in residues - low_confidence_residues:
        full_relative = _relative_sasa(residue_names.get(residue), full_sasa.get(residue))
        if full_relative is None or full_relative < config.sasa_relative_exposed_threshold:
            needs_block_context.add(residue)
    block_sasa: dict[int, float] = {}
    block_ids: dict[int, str] = {}
    for index, block in enumerate(blocks, start=1):
        block_context_residues = set(range(block.start, block.end + 1)) & set(residue_names)
        block_residues = block_context_residues & needs_block_context
        if not block_residues:
            continue
        current = calculate_residue_sasa(
            atoms,
            residue_names,
            residues=block_residues,
            context_residues=block_context_residues,
            probe_radius=config.sasa_probe_radius,
            slices=config.sasa_slices,
        )
        block_id = f"pae_block_{index}:{block.start}-{block.end}"
        for residue, sasa in current.items():
            block_sasa[residue] = sasa
            block_ids[residue] = block_id

    accessibility: dict[int, SurfaceAccessibility] = {}
    for residue in sorted(residues):
        full_value = full_sasa.get(residue)
        block_value = block_sasa.get(residue)
        full_relative = _relative_sasa(residue_names.get(residue), full_value)
        block_relative = _relative_sasa(residue_names.get(residue), block_value)
        plddt = residue_plddt.get(residue)
        low_confidence = plddt is not None and plddt < config.sasa_low_confidence_plddt_threshold
        if full_relative is not None and full_relative >= config.sasa_relative_exposed_threshold:
            classification = "exposed"
            if low_confidence:
                reason = f"Solvent-accessible in the full AlphaFold model, but low-confidence/disordered by pLDDT ({plddt:.1f})."
            else:
                reason = "Solvent-accessible in the full AlphaFold model."
        elif block_relative is not None and block_relative >= config.sasa_relative_exposed_threshold:
            classification = "conditional"
            if low_confidence:
                reason = (
                    "Buried in the full model but exposed within the PAE-coherent block; "
                    f"low-confidence/disordered by pLDDT ({plddt:.1f})."
                )
            else:
                reason = "Buried in the full model but exposed when SASA is computed within the PAE-coherent block."
        elif low_confidence:
            classification = "conditional"
            reason = f"Low pLDDT ({plddt:.1f}); treating as accessible but low-confidence/disordered rather than buried."
        elif full_relative is not None:
            classification = "buried"
            reason = "Low solvent accessibility in both full-model and available PAE-block contexts."
        else:
            classification = "uncertain"
            reason = "No usable SASA context was available."
        accessibility[residue] = SurfaceAccessibility(
            position=residue,
            full_model_sasa=round(full_value, 2) if full_value is not None else None,
            full_model_relative_sasa=round(full_relative, 4) if full_relative is not None else None,
            domain_context_sasa=round(block_value, 2) if block_value is not None else None,
            domain_context_relative_sasa=round(block_relative, 4) if block_relative is not None else None,
            classification=classification,
            reason=reason,
            pae_block_id=block_ids.get(residue),
        )
    return accessibility


def calculate_residue_sasa(
    atoms: dict[int, dict[str, tuple[float, float, float]]],
    residue_names: dict[int, str],
    *,
    residues: set[int],
    context_residues: set[int],
    probe_radius: float,
    slices: int = 20,
) -> dict[int, float]:
    structure = freesasa.Structure()
    added_atoms = 0
    for residue in sorted(context_residues):
        residue_name = residue_names.get(residue)
        atom_map = atoms.get(residue, {})
        if not residue_name or not atom_map:
            continue
        for atom_name, coord in atom_map.items():
            try:
                structure.addAtom(
                    atom_name,
                    residue_name,
                    str(residue),
                    "A",
                    float(coord[0]),
                    float(coord[1]),
                    float(coord[2]),
                )
                added_atoms += 1
            except Exception:
                continue
    if added_atoms == 0:
        return {residue: 0.0 for residue in residues}
    parameters = freesasa.Parameters(
        {
            "algorithm": "LeeRichards",
            "probe-radius": float(probe_radius),
            "n-slices": max(1, int(slices)),
        }
    )
    result = freesasa.calc(structure, parameters)
    residue_areas = result.residueAreas().get("A", {})
    residue_sasa = {residue: 0.0 for residue in residues}
    for residue in residues:
        area = residue_areas.get(str(residue))
        if area is not None:
            residue_sasa[residue] = float(area.total)
    return residue_sasa


def summarize_construct_quality(
    construct_start: int,
    construct_end: int,
    residue_plddt: dict[int, float],
    pae_matrix: list[list[float]] | PaeStats,
    config: AnalysisConfig,
) -> ConstructStructuralMetrics:
    pae_stats = _coerce_pae_stats(pae_matrix)
    plddt_values = [
        residue_plddt[residue]
        for residue in range(construct_start, construct_end + 1)
        if residue in residue_plddt
    ]
    mean_plddt = sum(plddt_values) / len(plddt_values) if plddt_values else None
    structured_fraction = None
    if plddt_values:
        structured_fraction = sum(
            1 for value in plddt_values if value >= config.plddt_structured_threshold
        ) / len(plddt_values)

    if construct_start < 1 or construct_end > len(pae_stats.matrix):
        raise ValueError("Construct boundaries exceed the PAE matrix")
    pae_start, pae_end = _clip_region_to_pae_bounds(construct_start, construct_end, pae_stats)
    mean_intra_pae = None
    max_intra_pae = None
    if pae_start is not None and pae_end is not None:
        mean_intra_pae = _mean_intra_block_pae(pae_start, pae_end, pae_stats)
        max_intra_pae = _max_block_pae(pae_start, pae_end, pae_stats)
    return ConstructStructuralMetrics(
        mean_plddt=round(mean_plddt, 2) if mean_plddt is not None else None,
        min_plddt=round(min(plddt_values), 2) if plddt_values else None,
        max_plddt=round(max(plddt_values), 2) if plddt_values else None,
        structured_fraction=round(structured_fraction, 4) if structured_fraction is not None else None,
        mean_intra_pae=round(mean_intra_pae, 2) if mean_intra_pae is not None else None,
        max_intra_pae=round(max_intra_pae, 2) if max_intra_pae is not None else None,
    )


def collect_region_split_diagnostics(
    region: Region,
    pae_matrix: list[list[float]] | PaeStats,
    config: AnalysisConfig,
    residue_plddt: dict[int, float] | None = None,
) -> list[ConstructSplitDiagnostic]:
    pae_stats = _coerce_pae_stats(pae_matrix)
    clipped_start, clipped_end = _clip_region_to_pae_bounds(region.start, region.end, pae_stats)
    if clipped_start is None or clipped_end is None:
        return []
    clipped_region = Region(
        start=clipped_start,
        end=clipped_end,
        label=region.label,
        source=region.source,
        confidence=region.confidence,
        metadata=dict(region.metadata),
    )
    domains = split_domains_from_pae(clipped_region, pae_stats, config, residue_plddt=residue_plddt)
    unique: dict[tuple[int, int, int], ConstructSplitDiagnostic] = {}
    for domain in domains:
        for split in domain.metadata.get("split_history", []):
            boundary_after = split.get("boundary_after")
            parent_start = split.get("parent_start")
            parent_end = split.get("parent_end")
            if boundary_after is None or parent_start is None or parent_end is None:
                continue
            key = (int(parent_start), int(boundary_after), int(parent_end))
            unique[key] = ConstructSplitDiagnostic(
                parent_start=int(parent_start),
                parent_end=int(parent_end),
                boundary_after=int(boundary_after),
                score=float(split.get("score", 0.0)),
                inter_block_pae=float(split.get("inter_block_pae", 0.0)),
                intra_block_pae=float(split.get("intra_block_pae", 0.0)),
                separation=float(split.get("separation", 0.0)),
                boundary_pae=float(split.get("boundary_pae", 0.0)),
                linker_plddt=(
                    float(split["linker_plddt"])
                    if split.get("linker_plddt") is not None
                    else None
                ),
            )
    return [unique[key] for key in sorted(unique)]


def _split_region(
    start: int,
    end: int,
    pae_stats: PaeStats,
    config: AnalysisConfig,
    residue_plddt: dict[int, float] | None,
    split_history: list[dict[str, float | int | None]] | None = None,
    *,
    depth: int,
) -> list[Region]:
    history = list(split_history or [])
    if end - start + 1 < 2 * config.min_domain_size:
        return [_build_domain_region(start, end, pae_stats, config, depth, history)]

    best_candidate: dict[str, float] | None = None
    for cut in range(start + config.min_domain_size - 1, end - config.min_domain_size + 1):
        inter_error = _mean_inter_block_pae(start, cut, cut + 1, end, pae_stats)
        left_intra = _mean_intra_block_pae(start, cut, pae_stats)
        right_intra = _mean_intra_block_pae(cut + 1, end, pae_stats)
        intra_error = (left_intra + right_intra) / 2.0
        separation = inter_error - intra_error
        if inter_error < config.pae_domain_threshold or separation < config.pae_separation_threshold:
            continue

        boundary_error = _mean_boundary_pae(
            start,
            cut,
            end,
            pae_stats,
            window=config.domain_linker_window,
        )
        linker_plddt = _mean_linker_plddt(
            cut,
            residue_plddt,
            start=start,
            end=end,
            window=config.domain_linker_window,
        )
        linker_bonus = 0.0
        if linker_plddt is not None:
            linker_bonus = max(0.0, (config.plddt_structured_threshold - linker_plddt) / 10.0)
        score = separation + max(0.0, boundary_error - intra_error) * 0.35 + linker_bonus
        if best_candidate is None or score > best_candidate["score"]:
            best_candidate = {
                "cut": float(cut),
                "score": score,
                "inter_error": inter_error,
                "intra_error": intra_error,
                "separation": separation,
                "boundary_error": boundary_error,
                "linker_plddt": linker_plddt if linker_plddt is not None else float("nan"),
            }

    if best_candidate is None or depth >= config.max_domain_recursion_depth:
        return [_build_domain_region(start, end, pae_stats, config, depth, history)]

    best_cut = int(best_candidate["cut"])
    boundary_metadata: dict[str, float | int | None] = {
        "parent_start": start,
        "parent_end": end,
        "boundary_after": best_cut,
        "score": round(best_candidate["score"], 2),
        "inter_block_pae": round(best_candidate["inter_error"], 2),
        "intra_block_pae": round(best_candidate["intra_error"], 2),
        "separation": round(best_candidate["separation"], 2),
        "boundary_pae": round(best_candidate["boundary_error"], 2),
        "linker_plddt": None
        if math.isnan(best_candidate["linker_plddt"])
        else round(best_candidate["linker_plddt"], 2),
    }
    left = _split_region(
        start,
        best_cut,
        pae_stats,
        config,
        residue_plddt,
        history + [boundary_metadata],
        depth=depth + 1,
    )
    right = _split_region(
        best_cut + 1,
        end,
        pae_stats,
        config,
        residue_plddt,
        history + [boundary_metadata],
        depth=depth + 1,
    )
    return left + right


def _mean_inter_block_pae(
    left_start: int,
    left_end: int,
    right_start: int,
    right_end: int,
    pae_stats: PaeStats,
) -> float:
    left_count = left_end - left_start + 1
    right_count = right_end - right_start + 1
    if left_count <= 0 or right_count <= 0:
        return 0.0
    forward = _block_sum(pae_stats, left_start, left_end, right_start, right_end)
    reverse = _block_sum(pae_stats, right_start, right_end, left_start, left_end)
    return (forward + reverse) / (2 * left_count * right_count)


def _mean_intra_block_pae(
    start: int,
    end: int,
    pae_stats: PaeStats,
) -> float:
    size = end - start + 1
    if size <= 1:
        return 0.0
    block_sum = _block_sum(pae_stats, start, end, start, end)
    diagonal_sum = _diagonal_sum(pae_stats, start, end)
    return (block_sum - diagonal_sum) / (size * (size - 1))


def _mean_boundary_pae(
    start: int,
    cut: int,
    end: int,
    pae_stats: PaeStats,
    *,
    window: int,
) -> float:
    left_start = max(start, cut - window + 1)
    right_end = min(end, cut + window)
    return _mean_inter_block_pae(left_start, cut, cut + 1, right_end, pae_stats)


def _mean_linker_plddt(
    cut: int,
    residue_plddt: dict[int, float] | None,
    *,
    start: int,
    end: int,
    window: int,
) -> float | None:
    if not residue_plddt:
        return None
    values = [
        residue_plddt[residue]
        for residue in range(max(start, cut - window + 1), min(end, cut + window) + 1)
        if residue in residue_plddt
    ]
    if not values:
        return None
    return sum(values) / len(values)


def _build_domain_region(
    start: int,
    end: int,
    pae_stats: PaeStats,
    config: AnalysisConfig,
    depth: int,
    split_history: list[dict[str, float | int | None]],
) -> Region:
    mean_intra = _mean_intra_block_pae(start, end, pae_stats)
    confidence = max(0.4, min(0.95, 1.0 - (mean_intra / max(config.pae_domain_threshold * 2.0, 1.0))))
    return Region(
        start=start,
        end=end,
        label="Predicted domain",
        source="AlphaFold",
        confidence=round(confidence, 2),
        metadata={
            "mean_intra_pae": round(mean_intra, 2),
            "segmentation_method": "pae_boundary_scoring",
            "recursion_depth": depth,
            "split_history": split_history,
        },
    )


def _coerce_pae_stats(pae_matrix: list[list[float]] | PaeStats) -> PaeStats:
    if isinstance(pae_matrix, PaeStats):
        return pae_matrix
    return build_pae_stats(pae_matrix)


def _clip_region_to_pae_bounds(
    start: int,
    end: int,
    pae_stats: PaeStats,
) -> tuple[int | None, int | None]:
    if start > end:
        return None, None
    max_residue = len(pae_stats.matrix)
    clipped_start = max(1, start)
    clipped_end = min(end, max_residue)
    if clipped_start > clipped_end:
        return None, None
    return clipped_start, clipped_end


def _block_sum(
    pae_stats: PaeStats,
    row_start: int,
    row_end: int,
    col_start: int,
    col_end: int,
) -> float:
    prefix = pae_stats.prefix_sum
    r1 = row_start - 1
    r2 = row_end
    c1 = col_start - 1
    c2 = col_end
    return prefix[r2][c2] - prefix[r1][c2] - prefix[r2][c1] + prefix[r1][c1]


def _diagonal_sum(pae_stats: PaeStats, start: int, end: int) -> float:
    return pae_stats.diagonal_prefix[end] - pae_stats.diagonal_prefix[start - 1]


def _max_block_pae(start: int, end: int, pae_stats: PaeStats) -> float | None:
    max_value: float | None = None
    matrix = pae_stats.matrix
    for row in range(start - 1, end):
        for col in range(start - 1, end):
            value = float(matrix[row][col])
            if max_value is None or value > max_value:
                max_value = value
    return max_value


def _distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b, strict=True)))


def _distance_squared(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return sum((x - y) ** 2 for x, y in zip(a, b, strict=True))


def _relative_sasa(residue_name: str | None, sasa: float | None) -> float | None:
    if sasa is None:
        return None
    max_asa = _RESIDUE_MAX_ASA.get(str(residue_name or "").upper())
    if not max_asa:
        return None
    return max(0.0, min(1.0, sasa / max_asa))


def _spatial_cell(point: tuple[float, float, float], cell_size: float) -> tuple[int, int, int]:
    return (
        math.floor(point[0] / cell_size),
        math.floor(point[1] / cell_size),
        math.floor(point[2] / cell_size),
    )


def _is_surface_exposed(residue: int, atoms: dict[int, dict[str, tuple[float, float, float]]]) -> bool | None:
    center = atoms.get(residue, {}).get("CA")
    if not center:
        center = atoms.get(residue, {}).get("CB")
    if not center:
        return None
    neighbors = 0
    for other_residue, atom_map in atoms.items():
        if other_residue == residue:
            continue
        other_center = atom_map.get("CA") or atom_map.get("CB")
        if not other_center:
            continue
        if _distance(center, other_center) <= 10.0:
            neighbors += 1
    return neighbors < 14
