#!/usr/bin/env python3
"""Render an overview figure of the LD-panel study design.

Outputs:
  results/summary/plots/00_study_design.pdf   (vector, for poster)
  results/summary/plots/00_study_design.png   (preview)

Usage:
  python3 scripts/plot_study_design.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OUTDIR = Path("results/summary/plots")
OUTDIR.mkdir(parents=True, exist_ok=True)

# ── Style ────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.linewidth": 0,
})

PALETTE = {
    "source":   "#E8EEF7",   # very light blue
    "source_e": "#3C5A8A",
    "gwas":     "#FFF0DD",   # pale orange
    "gwas_e":   "#B36B00",
    "panel":    "#E6F3E9",   # pale green
    "panel_e":  "#2A7F3A",
    "pheno":    "#F1E6F5",   # pale violet
    "pheno_e":  "#6A3C8C",
    "pgs":      "#FCE6E6",   # pale red
    "pgs_e":    "#A33636",
    "eval":     "#E6EEEE",   # pale teal
    "eval_e":   "#246B6B",
    "arrow":    "#444444",
}

# Canvas: x in [0, 100], y in [0, 100]
fig, ax = plt.subplots(figsize=(11, 7.5))
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.set_aspect("equal")
ax.axis("off")


def box(x, y, w, h, title, body, face, edge, fontsize=9, title_fs=10):
    """Rounded box centred at (x+w/2, y+h/2)."""
    p = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.4,rounding_size=1.2",
        linewidth=1.4, edgecolor=edge, facecolor=face,
    )
    ax.add_patch(p)
    ax.text(x + w / 2, y + h - 2.6, title,
            ha="center", va="top", fontsize=title_fs, weight="bold",
            color=edge)
    ax.text(x + w / 2, y + h - 6.0, body,
            ha="center", va="top", fontsize=fontsize, color="#222")


def arrow(x1, y1, x2, y2, label=None, label_offset=(0, 1.6),
          shrinkA=2, shrinkB=2):
    a = FancyArrowPatch(
        (x1, y1), (x2, y2),
        arrowstyle="-|>", mutation_scale=14,
        linewidth=1.4, color=PALETTE["arrow"],
        shrinkA=shrinkA, shrinkB=shrinkB,
    )
    ax.add_patch(a)
    if label:
        ax.text((x1 + x2) / 2 + label_offset[0],
                (y1 + y2) / 2 + label_offset[1],
                label, ha="center", va="center",
                fontsize=8, color="#333",
                bbox=dict(facecolor="white", edgecolor="none", pad=1.0))


# ── Title ────────────────────────────────────────────────────────────────────
ax.text(50, 96, "LD-panel study design",
        ha="center", va="center", fontsize=14, weight="bold")
ax.text(50, 92,
        "Simulated genotypes  →  GWAS  →  PGS with varying LD panels  →  accuracy",
        ha="center", va="center", fontsize=10, color="#444", style="italic")

# ── Row 1: HAPNEST source ────────────────────────────────────────────────────
box(28, 74, 44, 12,
    "HAPNEST simulated genomes",
    "6 superpopulations  (AFR, AMR, CSA, EAS, EUR, MID)\n"
    "22 autosomes  •  HapMap3 sites  •  MAF ≥ 0.01",
    PALETTE["source"], PALETTE["source_e"], fontsize=9)

# ── Row 2: GWAS cohort  &  LD panels ─────────────────────────────────────────
box(2, 52, 42, 16,
    "GWAS cohort  (AMR)",
    "N_gwas = 100,000  (independent AMR subjects)\n"
    "80 / 20 split:\n"
    "  • train  N = 80,000   → GWAS + 20 PCs\n"
    "  • test   N = 20,000   → held-out evaluation",
    PALETTE["gwas"], PALETTE["gwas_e"], fontsize=9)

box(56, 52, 42, 16,
    "LD reference panels",
    "Panel ancestry × panel N  (6 × 3 = 18)\n"
    "AMR: independent panel cohort\n"
    "AFR, CSA, EAS, EUR, MID: HAPNEST superpop\n"
    "N ∈ {500;  5,000;  50,000}",
    PALETTE["panel"], PALETTE["panel_e"], fontsize=9)

arrow(50, 73.5, 23, 68.5, "subset",  label_offset=(-3, 0),
      shrinkA=0, shrinkB=0)
arrow(50, 73.5, 77, 68.5, "subsets", label_offset=(3, 0),
      shrinkA=0, shrinkB=0)

# ── Row 3: phenotypes ────────────────────────────────────────────────────────
box(2, 28, 42, 18,
    "Phenotype simulation",
    "Trait:  quantitative  •  binary (prev 10%, 30%)\n"
    "h²  ∈ {0.20, 0.40, 0.60, 0.80}\n"
    "p_causal  ∈ {0.001, 0.01, 0.05, 0.20}\n"
    "Effects:  gaussian  •  spike-and-slab",
    PALETTE["pheno"], PALETTE["pheno_e"], fontsize=9)

arrow(23, 52, 23, 46)

# ── Row 3 right: PGS methods ─────────────────────────────────────────────────
box(56, 28, 42, 18,
    "PGS methods",
    "LDpred2-auto  (bigsnpr)\n"
    "PRS-CS\n"
    "\nTrained on GWAS sumstats + LD panel",
    PALETTE["pgs"], PALETTE["pgs_e"], fontsize=9)

arrow(77, 52, 77, 46)

# Cross arrow: sumstats from GWAS feed into PGS methods
arrow(44, 37, 56, 37, "sumstats", label_offset=(0, 2.0))

# ── Row 4: evaluation ────────────────────────────────────────────────────────
box(20, 2, 60, 18,
    "Evaluation on held-out test set  (N = 20,000)",
    "Quantitative:  incremental R²\n"
    "Binary:  incremental AUC\n"
    "Compare across panel ancestry × panel N\n"
    "× h² × p_causal × effect distribution",
    PALETTE["eval"], PALETTE["eval_e"], fontsize=9)

arrow(23, 28, 35, 20)
arrow(77, 28, 65, 20)

# ── Save ─────────────────────────────────────────────────────────────────────
plt.tight_layout(pad=0.2)
pdf_path = OUTDIR / "00_study_design.pdf"
png_path = OUTDIR / "00_study_design.png"
fig.savefig(pdf_path, bbox_inches="tight")
fig.savefig(png_path, dpi=200, bbox_inches="tight")
print(f"[plot] wrote {pdf_path}")
print(f"[plot] wrote {png_path}")
