from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
MPL_DIR = ROOT / ".agdesign2" / "mplconfig"
CACHE_DIR = ROOT / ".agdesign2" / "cache_home"
MPL_DIR.mkdir(parents=True, exist_ok=True)
CACHE_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_DIR))
os.environ.setdefault("XDG_CACHE_HOME", str(CACHE_DIR))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.collections import LineCollection
import numpy as np

from .models import ConstructDetail, Region


@dataclass(slots=True)
class _QualityPlotContext:
    residues: list[int]
    plddt_segments: list[list[tuple[int, float]]]
    plddt_colors: list[str]
    pae_submatrix: Any
    pae_extent: tuple[int, int, int, int]
    pae_max: float
    topology_intervals: list[tuple[int, int, str]]
    ectodomain: Region


class ConstructAssetGenerator:
    def __init__(self, pymol_executable: str = "pymol") -> None:
        self.pymol_executable = self._resolve_pymol_executable(pymol_executable)

    def generate_assets(
        self,
        *,
        construct: ConstructDetail,
        pdb_path: Path,
        residue_plddt: dict[int, float],
        pae_matrix: list[list[float]],
        ectodomain: Region,
        output_dir: Path,
        residue_annotations: list[Any] | None = None,
        render_structure_image: bool = True,
        render_quality_plot: bool = True,
    ) -> tuple[str | None, str | None]:
        output_dir.mkdir(parents=True, exist_ok=True)
        safe_name = _safe_filename(construct.name)
        structure_path = output_dir / f"{safe_name}_structure.png"
        quality_path = output_dir / f"{safe_name}_quality.png"

        structure_relative = None
        quality_relative = None

        if render_structure_image and self.is_pymol_available():
            self.render_structure_image(
                pdb_path=pdb_path,
                start=construct.start,
                end=construct.end,
                output_path=structure_path,
            )
            if structure_path.exists():
                structure_relative = structure_path.name

        if render_quality_plot:
            context = self._build_quality_plot_context(
                residue_plddt=residue_plddt,
                pae_matrix=pae_matrix,
                ectodomain=ectodomain,
                residue_annotations=residue_annotations,
            )
            self._render_quality_plot_from_context(
                context=context,
                construct=construct,
                output_path=quality_path,
            )
            quality_relative = quality_path.name
        return structure_relative, quality_relative

    def generate_assets_for_constructs(
        self,
        *,
        constructs: list[ConstructDetail],
        pdb_path: Path,
        residue_plddt: dict[int, float],
        pae_matrix: list[list[float]],
        ectodomain: Region,
        output_dir: Path,
        residue_annotations: list[Any] | None = None,
        render_structure_image: bool = True,
        render_quality_plot: bool = True,
    ) -> dict[str, tuple[str | None, str | None]]:
        output_dir.mkdir(parents=True, exist_ok=True)
        results: dict[str, tuple[str | None, str | None]] = {}
        primary_by_range: dict[tuple[int, int], ConstructDetail] = {}
        for construct in constructs:
            primary_by_range.setdefault((construct.start, construct.end), construct)

        structure_outputs: dict[str, Path] = {}
        if render_structure_image and self.is_pymol_available():
            missing_structure_ranges: list[tuple[str, int, int]] = []
            for construct in primary_by_range.values():
                safe_name = _safe_filename(construct.name)
                structure_path = output_dir / f"{safe_name}_structure.png"
                structure_outputs[construct.name] = structure_path
                if not structure_path.exists():
                    missing_structure_ranges.append((construct.name, construct.start, construct.end))
            if missing_structure_ranges:
                self.render_structure_images(
                    pdb_path=pdb_path,
                    construct_ranges=missing_structure_ranges,
                    output_paths=structure_outputs,
                )
        quality_context = None
        if render_quality_plot:
            quality_context = self._build_quality_plot_context(
                residue_plddt=residue_plddt,
                pae_matrix=pae_matrix,
                ectodomain=ectodomain,
                residue_annotations=residue_annotations,
            )
        assets_by_range: dict[tuple[int, int], tuple[str | None, str | None]] = {}
        for range_key, construct in primary_by_range.items():
            safe_name = _safe_filename(construct.name)
            structure_relative = None
            quality_relative = None
            structure_path = structure_outputs.get(construct.name)
            if structure_path is not None and structure_path.exists():
                structure_relative = structure_path.name
            if render_quality_plot:
                quality_path = output_dir / f"{safe_name}_quality.png"
                if not quality_path.exists():
                    self._render_quality_plot_from_context(
                        context=quality_context,
                        construct=construct,
                        output_path=quality_path,
                    )
                quality_relative = quality_path.name
            assets_by_range[range_key] = (structure_relative, quality_relative)
        for construct in constructs:
            results[construct.name] = assets_by_range[(construct.start, construct.end)]
        return results

    def is_pymol_available(self) -> bool:
        return shutil.which(self.pymol_executable) is not None

    def _resolve_pymol_executable(self, pymol_executable: str) -> str:
        candidates = [pymol_executable]
        env_override = os.environ.get("PYMOL_EXECUTABLE")
        if env_override:
            candidates.insert(0, env_override)
        candidates.extend(
            [
                "/opt/homebrew/bin/pymol",
                "/usr/local/bin/pymol",
                "/Applications/PyMOL.app/Contents/bin/pymol",
            ]
        )
        seen: set[str] = set()
        for candidate in candidates:
            if candidate in seen:
                continue
            seen.add(candidate)
            if shutil.which(candidate):
                return candidate
        return pymol_executable

    def render_structure_image(
        self,
        *,
        pdb_path: Path,
        start: int,
        end: int,
        output_path: Path,
    ) -> None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as png_handle:
            temp_output = Path(png_handle.name)
        commands = [
            "reinitialize",
            f'load "{pdb_path.resolve()}", af',
            "remove not polymer.protein",
            f"create construct, af and resi {start}-{end}",
            "hide everything, all",
            "show cartoon, construct",
            "spectrum b, red_yellow_green_cyan_blue, construct, minimum=0, maximum=100",
            "select construct_cysteines, construct and resn CYS",
            "show sticks, construct_cysteines",
            "color yellow, construct_cysteines",
            "set stick_radius, 0.3, construct_cysteines",
            "set stick_quality, 16",
            "set sphere_scale, 0.22, construct_cysteines and name SG",
            "bg_color white",
            "set ray_opaque_background, off",
            "set antialias, 2",
            "set cartoon_side_chain_helper, on",
            "orient construct",
            "zoom construct, 6",
            f"png {temp_output}, 750, 750, 150, ray=0",
            "quit",
        ]
        with tempfile.NamedTemporaryFile("w", suffix=".pml", delete=False, encoding="utf-8") as handle:
            handle.write("\n".join(commands))
            script_path = Path(handle.name)
        try:
            subprocess.run(
                [self.pymol_executable, "-cq", str(script_path)],
                check=True,
                capture_output=True,
                text=True,
            )
            if not temp_output.exists():
                raise FileNotFoundError("PyMOL completed without writing a structure image.")
            shutil.move(str(temp_output), str(output_path))
        finally:
            script_path.unlink(missing_ok=True)
            temp_output.unlink(missing_ok=True)

    def render_structure_images(
        self,
        *,
        pdb_path: Path,
        construct_ranges: list[tuple[str, int, int]],
        output_paths: dict[str, Path],
    ) -> None:
        if not construct_ranges:
            return
        commands = [
            "reinitialize",
            f'load "{pdb_path.resolve()}", af',
            "remove not polymer.protein",
            "hide everything, all",
            "bg_color white",
            "set ray_opaque_background, off",
            "set antialias, 2",
            "set cartoon_side_chain_helper, on",
            "set stick_quality, 16",
        ]
        for index, (name, start, end) in enumerate(construct_ranges):
            output_path = output_paths[name]
            output_path.parent.mkdir(parents=True, exist_ok=True)
            object_name = f"construct_{index}"
            cysteine_name = f"{object_name}_cysteines"
            commands.extend(
                [
                    f"delete {object_name}",
                    f"delete {cysteine_name}",
                    f"create {object_name}, af and resi {start}-{end}",
                    f"hide everything, {object_name}",
                    f"show cartoon, {object_name}",
                    f"spectrum b, red_yellow_green_cyan_blue, {object_name}, minimum=0, maximum=100",
                    f"select {cysteine_name}, {object_name} and resn CYS",
                    f"show sticks, {cysteine_name}",
                    f"color yellow, {cysteine_name}",
                    f"set stick_radius, 0.3, {cysteine_name}",
                    f"set sphere_scale, 0.22, {cysteine_name} and name SG",
                    f"orient {object_name}",
                    f"zoom {object_name}, 6",
                    f'png "{output_path.resolve()}", 750, 750, 150, ray=0',
                ]
            )
        commands.append("quit")
        with tempfile.NamedTemporaryFile("w", suffix=".pml", delete=False, encoding="utf-8") as handle:
            handle.write("\n".join(commands))
            script_path = Path(handle.name)
        try:
            subprocess.run(
                [self.pymol_executable, "-cq", str(script_path)],
                check=True,
                capture_output=True,
                text=True,
            )
        finally:
            script_path.unlink(missing_ok=True)

    def render_quality_plot(
        self,
        *,
        residue_plddt: dict[int, float],
        pae_matrix: list[list[float]],
        ectodomain: Region,
        construct: ConstructDetail,
        residue_annotations: list[Any] | None = None,
        output_path: Path,
    ) -> None:
        context = self._build_quality_plot_context(
            residue_plddt=residue_plddt,
            pae_matrix=pae_matrix,
            ectodomain=ectodomain,
            residue_annotations=residue_annotations,
        )
        self._render_quality_plot_from_context(
            context=context,
            construct=construct,
            output_path=output_path,
        )

    def _build_quality_plot_context(
        self,
        *,
        residue_plddt: dict[int, float],
        pae_matrix: list[list[float]],
        ectodomain: Region,
        residue_annotations: list[Any] | None = None,
    ) -> _QualityPlotContext:
        residues = list(range(ectodomain.start, ectodomain.end + 1))
        plddt_values = [residue_plddt.get(residue, float("nan")) for residue in residues]
        pae_submatrix = [
            row[ectodomain.start - 1 : ectodomain.end]
            for row in pae_matrix[ectodomain.start - 1 : ectodomain.end]
        ]
        pae_submatrix, pae_extent = _downsample_pae_for_static_plot(
            pae_submatrix,
            start=ectodomain.start,
            end=ectodomain.end,
            max_size=420,
        )
        pae_max = max(30.0, max(max(row) for row in pae_submatrix) if pae_submatrix else 30.0)
        topology_intervals = _topology_intervals(residue_annotations, ectodomain.start, ectodomain.end)
        segments = []
        colors = []
        for left, right, value_left, value_right in zip(residues[:-1], residues[1:], plddt_values[:-1], plddt_values[1:], strict=False):
            segments.append([(left, value_left), (right, value_right)])
            colors.append(_plddt_color(_mean_numeric(value_left, value_right)))
        return _QualityPlotContext(
            residues=residues,
            plddt_segments=segments,
            plddt_colors=colors,
            pae_submatrix=pae_submatrix,
            pae_extent=pae_extent,
            pae_max=pae_max,
            topology_intervals=topology_intervals,
            ectodomain=ectodomain,
        )

    def _render_quality_plot_from_context(
        self,
        *,
        context: _QualityPlotContext,
        construct: ConstructDetail,
        output_path: Path,
    ) -> None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        ectodomain = context.ectodomain
        topology_intervals = context.topology_intervals
        fig = plt.figure(figsize=(10, 4.6), facecolor="white")
        grid = fig.add_gridspec(1, 2, width_ratios=[1.25, 1], wspace=0.34)

        ax_plddt = fig.add_subplot(grid[0])
        confidence_bands = (
            (0, 50, "#ff7d45", 0.14),
            (50, 70, "#ffdb13", 0.16),
            (70, 90, "#65cbf3", 0.16),
            (90, 100, "#0053d6", 0.12),
        )
        for ymin, ymax, color, alpha in confidence_bands:
            ax_plddt.axhspan(ymin, ymax, color=color, alpha=alpha, linewidth=0)
        for start, end, location in topology_intervals:
            ax_plddt.axvspan(start, end, color=_topology_color(location), alpha=0.08, linewidth=0)
            ax_plddt.axvspan(start, end, ymin=0.0, ymax=0.045, color=_topology_color(location), alpha=0.85, linewidth=0)
        if context.plddt_segments:
            ax_plddt.add_collection(LineCollection(context.plddt_segments, colors=context.plddt_colors, linewidths=2.4, capstyle="round"))
        ax_plddt.axvspan(construct.start, construct.end, color="#4c956c", alpha=0.18)
        ax_plddt.axvline(construct.start, color="#4c956c", linewidth=1.6)
        ax_plddt.axvline(construct.end, color="#4c956c", linewidth=1.6)
        ax_plddt.set_ylabel("pLDDT")
        ax_plddt.set_xlabel("Residue")
        ax_plddt.set_ylim(0, 100)
        ax_plddt.set_xlim(ectodomain.start, ectodomain.end)
        boundary_label = f"Residues {construct.start}-{construct.end}"
        ax_plddt.set_title(f"{boundary_label}\npLDDT with construct highlighted", fontsize=11)
        ax_plddt.grid(axis="y", alpha=0.28, linewidth=0.6, linestyle="--")
        for spine in ("top", "right"):
            ax_plddt.spines[spine].set_visible(False)

        ax_pae = fig.add_subplot(grid[1])
        image = ax_pae.imshow(
            context.pae_submatrix,
            origin="lower",
            cmap="plasma",
            vmin=0,
            vmax=context.pae_max,
            extent=context.pae_extent,
            aspect="equal",
            interpolation="nearest",
        )
        width = construct.end - construct.start + 1
        rect = patches.Rectangle(
            (construct.start, construct.start),
            width,
            width,
            linewidth=2.2,
            edgecolor="#4c956c",
            facecolor="#4c956c",
            alpha=0.18,
        )
        ax_pae.add_patch(rect)
        ax_pae.add_patch(
            patches.Rectangle(
                (construct.start, construct.start),
                width,
                width,
                linewidth=2.2,
                edgecolor="#4c956c",
                facecolor="none",
            )
        )
        if topology_intervals:
            axis_span = ectodomain.end - ectodomain.start + 1
            strip = max(2.0, axis_span * 0.025)
            top_y = ectodomain.end - strip
            for start, end, location in topology_intervals:
                span = end - start + 1
                color = _topology_color(location)
                ax_pae.add_patch(
                    patches.Rectangle(
                        (start, top_y),
                        span,
                        strip,
                        facecolor=color,
                        edgecolor="none",
                        alpha=0.88,
                        zorder=4,
                    )
                )
                ax_pae.add_patch(
                    patches.Rectangle(
                        (ectodomain.start, start),
                        strip,
                        span,
                        facecolor=color,
                        edgecolor="none",
                        alpha=0.88,
                        zorder=4,
                    )
                )
                ax_pae.add_patch(
                    patches.Rectangle(
                        (start, start),
                        span,
                        span,
                        facecolor="none",
                        edgecolor=color,
                        alpha=0.22,
                        linewidth=0.7,
                        zorder=4,
                    )
                )
        ax_pae.set_xlabel("Residue")
        ax_pae.set_ylabel("Residue")
        ax_pae.set_title(f"{boundary_label}\nPAE with construct block highlighted", fontsize=11)
        colorbar = fig.colorbar(image, ax=ax_pae, fraction=0.046, pad=0.04)
        colorbar.set_label("Predicted aligned error", fontsize=9)
        colorbar.ax.tick_params(labelsize=8)

        for axis in (ax_plddt, ax_pae):
            axis.tick_params(labelsize=9)
            axis.set_facecolor("white")
        if topology_intervals:
            legend_locations = []
            for _, _, location in topology_intervals:
                if location not in legend_locations:
                    legend_locations.append(location)
            handles = [
                patches.Patch(facecolor=_topology_color(location), label=_topology_label(location), alpha=0.85)
                for location in legend_locations
            ]
            ax_plddt.legend(
                handles=handles,
                loc="upper center",
                bbox_to_anchor=(0.5, -0.2),
                ncol=min(4, len(handles)),
                frameon=False,
                fontsize=8,
            )

        fig.savefig(output_path, dpi=180, bbox_inches="tight", facecolor="white")
        plt.close(fig)


def _safe_filename(text: str) -> str:
    return "".join(ch if ch.isalnum() or ch in {"-", "_"} else "_" for ch in text)


def _plddt_color(score: float) -> str:
    if score >= 90:
        return "#0053d6"
    if score >= 70:
        return "#65cbf3"
    if score >= 50:
        return "#ffdb13"
    return "#ff7d45"


def _topology_intervals(
    residue_annotations: list[Any] | None,
    start: int,
    end: int,
) -> list[tuple[int, int, str]]:
    if not residue_annotations:
        return []
    locations: dict[int, str] = {}
    for annotation in residue_annotations:
        position = _annotation_value(annotation, "position")
        if position is None:
            continue
        try:
            residue = int(position)
        except (TypeError, ValueError):
            continue
        if start <= residue <= end:
            locations[residue] = _normalize_topology_location(_annotation_value(annotation, "topology_location"))
    if not locations:
        return []
    intervals: list[tuple[int, int, str]] = []
    current_start = start
    current_location = locations.get(start, "unknown")
    for residue in range(start + 1, end + 1):
        location = locations.get(residue, "unknown")
        if location == current_location:
            continue
        intervals.append((current_start, residue - 1, current_location))
        current_start = residue
        current_location = location
    intervals.append((current_start, end, current_location))
    return intervals


def _annotation_value(annotation: Any, key: str) -> Any:
    if isinstance(annotation, dict):
        return annotation.get(key)
    return getattr(annotation, key, None)


def _normalize_topology_location(value: Any) -> str:
    text = str(value or "").lower()
    if "membrane" in text or "transmembrane" in text:
        return "membrane"
    if "intra" in text or "cyto" in text:
        return "intracellular"
    if "extra" in text or "lumen" in text:
        return "extracellular"
    return "unknown"


def _topology_color(location: str) -> str:
    return {
        "extracellular": "#12917e",
        "membrane": "#f09e26",
        "intracellular": "#6952a8",
        "unknown": "#969696",
    }.get(location, "#969696")


def _topology_label(location: str) -> str:
    return {
        "extracellular": "Extracellular",
        "membrane": "Membrane",
        "intracellular": "Intracellular",
        "unknown": "Unknown",
    }.get(location, "Unknown")


def _mean_numeric(*values: float) -> float:
    numeric = [float(value) for value in values if value == value]
    return sum(numeric) / len(numeric) if numeric else 0.0


def _downsample_pae_for_static_plot(
    matrix: list[list[float]],
    *,
    start: int,
    end: int,
    max_size: int,
) -> tuple[list[list[float]], tuple[int, int, int, int]]:
    size = len(matrix)
    extent = (start, end + 1, start, end + 1)
    if size == 0 or size <= max_size:
        return matrix, extent
    sampled_indices = np.linspace(0, size - 1, max_size).astype(int)
    array = np.asarray(matrix, dtype=float)
    return array[np.ix_(sampled_indices, sampled_indices)].tolist(), extent
