from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.lines import Line2D
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "assets" / "figures"

COLORS = {
    "black": "#111111",
    "text": "#222222",
    "muted": "#6B7280",
    "line": "#4B5563",
    "blue": "#2166AC",
    "blue_light": "#DCEBF5",
    "green": "#2E7D32",
    "green_light": "#E7F3E7",
    "orange": "#D9791F",
    "orange_light": "#FFF0DD",
    "purple": "#6A00A8",
    "gray": "#D1D5DB",
    "gray_light": "#F3F4F6",
}


def set_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 7.5,
            "axes.titlesize": 8.5,
            "axes.labelsize": 7.5,
            "xtick.labelsize": 6.8,
            "ytick.labelsize": 6.8,
            "legend.fontsize": 6.8,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.linewidth": 0.8,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def panel_label(ax, label: str, title: str) -> None:
    ax.text(-0.10, 1.08, label, transform=ax.transAxes, fontsize=13, fontweight="bold", va="bottom", ha="left", color=COLORS["black"])
    ax.text(0.00, 1.095, title, transform=ax.transAxes, fontsize=8.8, fontweight="bold", va="bottom", ha="left", color=COLORS["black"])


def box(ax, xy, w, h, text, *, fc, ec, fontsize=7.1, weight="normal") -> None:
    patch = patches.FancyBboxPatch(
        xy,
        w,
        h,
        boxstyle="round,pad=0.012,rounding_size=0.018",
        transform=ax.transAxes,
        facecolor=fc,
        edgecolor=ec,
        linewidth=1.0,
        clip_on=False,
    )
    ax.add_patch(patch)
    ax.text(
        xy[0] + w / 2,
        xy[1] + h / 2,
        text,
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=fontsize,
        fontweight=weight,
        color=COLORS["text"],
        linespacing=1.05,
        clip_on=False,
    )


def arrow(ax, start, end, *, color=COLORS["line"], lw=1.0) -> None:
    ax.annotate(
        "",
        xy=end,
        xycoords=ax.transAxes,
        xytext=start,
        textcoords=ax.transAxes,
        arrowprops=dict(arrowstyle="-|>", lw=lw, color=color, shrinkA=3, shrinkB=3),
    )


def axis_off(ax) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")


def draw_panel_a(ax) -> None:
    axis_off(ax)
    panel_label(ax, "a", "Scope-first workflow")

    box(ax, (0.04, 0.72), 0.18, 0.13, "UniProt\ntopology", fc=COLORS["blue_light"], ec=COLORS["blue"], weight="bold")
    box(ax, (0.29, 0.72), 0.18, 0.13, "AlphaFold\npLDDT/PAE", fc=COLORS["gray_light"], ec=COLORS["line"], weight="bold")
    box(ax, (0.54, 0.72), 0.18, 0.13, "InterPro\nPDB", fc=COLORS["orange_light"], ec=COLORS["orange"], weight="bold")
    box(ax, (0.79, 0.72), 0.17, 0.13, "PTM/Cys\nwarnings", fc=COLORS["gray_light"], ec=COLORS["gray"], weight="bold")

    for x in (0.13, 0.38, 0.63, 0.875):
        arrow(ax, (x, 0.72), (0.50, 0.52))

    box(ax, (0.22, 0.41), 0.56, 0.13, "Design region", fc="white", ec=COLORS["black"], fontsize=8.0, weight="bold")
    ax.text(0.50, 0.36, "mature secreted, extracellular, or membrane-expression scope", transform=ax.transAxes, ha="center", va="top", fontsize=6.7, color=COLORS["muted"])

    box(ax, (0.08, 0.13), 0.35, 0.12, "Strict track", fc=COLORS["green_light"], ec=COLORS["green"], fontsize=8.0, weight="bold")
    box(ax, (0.57, 0.13), 0.35, 0.12, "Lenient track", fc=COLORS["orange_light"], ec=COLORS["orange"], fontsize=8.0, weight="bold")
    arrow(ax, (0.42, 0.41), (0.25, 0.25), color=COLORS["green"])
    arrow(ax, (0.58, 0.41), (0.75, 0.25), color=COLORS["orange"])


def synthetic_plddt() -> tuple[np.ndarray, np.ndarray]:
    x = np.arange(1, 221)
    y = np.empty_like(x, dtype=float)
    y[:16] = np.linspace(43, 68, 16)
    y[16:70] = 84 + 5 * np.sin(np.linspace(0, 2 * np.pi, 54))
    y[70:86] = 55 + 4 * np.sin(np.linspace(0, np.pi, 16))
    y[86:142] = 88 + 4 * np.sin(np.linspace(0, 2 * np.pi, 56))
    y[142:156] = 62 + 5 * np.sin(np.linspace(0, np.pi, 14))
    y[156:] = 84 + 6 * np.sin(np.linspace(0, 2 * np.pi, 64))
    return x, y


def draw_panel_b(ax) -> None:
    panel_label(ax, "b", "pLDDT defines candidate intervals")
    x, y = synthetic_plddt()
    ax.fill_between(x, 0, y, color=COLORS["blue_light"], lw=0)
    ax.plot(x, y, color=COLORS["blue"], lw=1.8)
    ax.axhline(70, color=COLORS["green"], lw=1.1, ls=(0, (4, 3)))
    ax.axhline(60, color=COLORS["orange"], lw=1.1, ls=(0, (4, 3)))

    ax.broken_barh([(17, 53), (87, 55), (157, 63)], (18, 5), facecolors=COLORS["green"], alpha=0.90)
    ax.broken_barh([(17, 203)], (8, 5), facecolors=COLORS["orange"], alpha=0.88)
    ax.text(2, 20.5, "strict", color=COLORS["green"], fontsize=6.5, fontweight="bold", va="center")
    ax.text(2, 10.5, "lenient", color=COLORS["orange"], fontsize=6.5, fontweight="bold", va="center")
    ax.set_xlim(1, 235)
    ax.set_ylim(0, 100)
    ax.set_xlabel("Residue position")
    ax.set_ylabel("pLDDT")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(
        handles=[
            Line2D([0], [0], color=COLORS["green"], lw=1.1, ls=(0, (4, 3)), label="strict threshold"),
            Line2D([0], [0], color=COLORS["orange"], lw=1.1, ls=(0, (4, 3)), label="lenient threshold"),
        ],
        loc="upper left",
        frameon=False,
        handlelength=2.3,
        borderpad=0,
        labelspacing=0.25,
    )


def synthetic_pae() -> np.ndarray:
    rng = np.random.default_rng(11)
    n = 180
    matrix = np.full((n, n), 23.0)
    for start, end in ((0, 58), (58, 118), (118, 180)):
        size = end - start
        matrix[start:end, start:end] = 3.7 + rng.normal(0, 0.45, (size, size))
    matrix = np.clip(matrix + rng.normal(0, 0.5, matrix.shape), 0, 30)
    np.fill_diagonal(matrix, 0)
    return matrix


def draw_panel_c(ax) -> None:
    panel_label(ax, "c", "PAE splits strict intervals")
    matrix = synthetic_pae()
    im = ax.imshow(matrix, cmap="plasma", vmin=0, vmax=30, origin="lower", interpolation="nearest")
    for boundary in (58, 118):
        ax.axvline(boundary - 0.5, color="white", lw=1.5)
        ax.axhline(boundary - 0.5, color="white", lw=1.5)
    for start, end in ((0, 58), (58, 118), (118, 180)):
        ax.add_patch(patches.Rectangle((start - 0.5, start - 0.5), end - start, end - start, fill=False, ec="white", lw=1.6))
    ax.set_xlabel("Aligned residue")
    ax.set_ylabel("Aligned residue")
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.025)
    cbar.set_label("PAE (A)", fontsize=7)
    cbar.ax.tick_params(labelsize=6.5)


def draw_panel_d(ax) -> None:
    axis_off(ax)
    panel_label(ax, "d", "Boundary scoring and outputs")

    ax.text(0.04, 0.83, "Eligible strict cut", transform=ax.transAxes, fontsize=7.7, fontweight="bold", color=COLORS["black"])
    for i, label in enumerate(("size", "inter PAE", "delta PAE", "linker")):
        x = 0.04 + i * 0.23
        fc = [COLORS["gray_light"], COLORS["blue_light"], COLORS["green_light"], COLORS["orange_light"]][i]
        ec = [COLORS["gray"], COLORS["blue"], COLORS["green"], COLORS["orange"]][i]
        box(ax, (x, 0.66), 0.18, 0.11, label, fc=fc, ec=ec, fontsize=7.0, weight="bold")
    ax.text(0.04, 0.57, "score = delta PAE + boundary term + linker bonus", transform=ax.transAxes, fontsize=7.4, color=COLORS["text"])

    ax.text(0.04, 0.43, "strict calculated constructs", transform=ax.transAxes, fontsize=7.3, color=COLORS["green"], fontweight="bold")
    ax.plot([0.05, 0.47], [0.34, 0.34], color=COLORS["gray"], lw=5, solid_capstyle="round", transform=ax.transAxes)
    for x0, x1 in ((0.08, 0.15), (0.23, 0.30), (0.38, 0.45)):
        ax.plot([x0, x1], [0.34, 0.34], color=COLORS["green"], lw=11, solid_capstyle="butt", transform=ax.transAxes)

    ax.text(0.57, 0.43, "lenient calculated construct", transform=ax.transAxes, fontsize=7.3, color=COLORS["orange"], fontweight="bold")
    ax.plot([0.58, 0.94], [0.34, 0.34], color=COLORS["orange"], lw=11, solid_capstyle="butt", transform=ax.transAxes)

    box(ax, (0.11, 0.09), 0.25, 0.12, "compact\nsingle-domain-like", fc=COLORS["green_light"], ec=COLORS["green"], fontsize=7.0)
    box(ax, (0.38, 0.09), 0.22, 0.12, "filter\n>=50 aa", fc=COLORS["gray_light"], ec=COLORS["gray"], fontsize=7.0)
    box(ax, (0.66, 0.09), 0.25, 0.12, "larger\nmulti-domain", fc=COLORS["orange_light"], ec=COLORS["orange"], fontsize=7.0)


def draw_visual_key(ax) -> None:
    axis_off(ax)
    ax.text(0.00, 0.62, "Visual key", transform=ax.transAxes, fontsize=7.2, fontweight="bold", color=COLORS["black"], va="center")
    items = [
        (0.12, COLORS["blue"], "pLDDT"),
        (0.24, COLORS["purple"], "PAE heat map"),
        (0.42, COLORS["green"], "strict track"),
        (0.57, COLORS["orange"], "lenient track"),
        (0.73, COLORS["gray"], "filters/annotations"),
    ]
    for x, color, label in items:
        ax.add_patch(patches.Rectangle((x, 0.43), 0.035, 0.30, transform=ax.transAxes, facecolor=color, edgecolor="none"))
        ax.text(x + 0.045, 0.58, label, transform=ax.transAxes, fontsize=6.9, va="center", color=COLORS["text"])


def write_legend() -> None:
    legend = """# Figure legend

**OpenAntigens strict and lenient calculated construct generation.** **a,** Input annotations are first constrained to a topology-defined design region, such as a mature secreted chain, extracellular region, or membrane-expression scope. **b,** AlphaFold pLDDT is used as a one-dimensional local-confidence signal to call candidate structured intervals. The strict track uses a pLDDT threshold of 70, permits low-confidence gaps up to 8 residues, and requires seed intervals of at least 25 residues. The lenient track uses a pLDDT threshold of 60, permits gaps up to 12 residues, and preserves larger intervals of at least 25 residues. **c,** Strict pLDDT seed intervals are further evaluated with the AlphaFold predicted aligned error (PAE) matrix. Candidate split boundaries are considered only when both resulting fragments are at least 40 residues, interblock PAE is at least 12 A, and the difference between interblock and mean intrablock PAE is at least 4 A. **d,** Eligible strict split boundaries are ranked using a score that combines PAE separation, local boundary PAE enrichment, and a low-pLDDT linker bonus. The highest-scoring split is applied recursively until no eligible split remains or the recursion limit is reached. Final calculated constructs are filtered to remain within the design region, must be at least 50 residues, and are annotated with structural diagnostics, PTMs, cysteine warnings, furin motifs, homolog equivalents, and other evidence.
"""
    (OUT_DIR / "openantigen_construct_algorithm_legend.md").write_text(legend, encoding="utf-8")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    set_style()

    fig = plt.figure(figsize=(7.2, 5.8))
    gs = fig.add_gridspec(
        3,
        2,
        height_ratios=(1.0, 1.0, 0.12),
        left=0.07,
        right=0.985,
        bottom=0.07,
        top=0.94,
        wspace=0.26,
        hspace=0.55,
    )
    axes = [
        fig.add_subplot(gs[0, 0]),
        fig.add_subplot(gs[0, 1]),
        fig.add_subplot(gs[1, 0]),
        fig.add_subplot(gs[1, 1]),
    ]
    key_ax = fig.add_subplot(gs[2, :])

    draw_panel_a(axes[0])
    draw_panel_b(axes[1])
    draw_panel_c(axes[2])
    draw_panel_d(axes[3])
    draw_visual_key(key_ax)

    for suffix, kwargs in {
        "svg": {},
        "pdf": {},
        "png": {"dpi": 450},
    }.items():
        fig.savefig(OUT_DIR / f"openantigen_construct_algorithm.{suffix}", bbox_inches="tight", **kwargs)
    plt.close(fig)
    write_legend()


if __name__ == "__main__":
    main()
