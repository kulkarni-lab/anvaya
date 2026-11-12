#!/usr/bin/env python3
"""Assembles mock Figures 1-6 as layout templates from the mock pipeline output."""

import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
from pdf2image import convert_from_path

OUT = "mockfigs/figures"
os.makedirs(OUT, exist_ok=True)
DPI = 200


def load(pdf, dpi=DPI):
    return np.array(convert_from_path(pdf, dpi=dpi)[0])


def panel(ax, img, label=None, lx=-0.02, ly=1.04):
    ax.imshow(img)
    ax.axis("off")
    if label:
        ax.text(lx, ly, label, transform=ax.transAxes, fontsize=15,
                fontweight="bold", va="bottom", ha="left")


def save(fig, name):
    path = f"{OUT}/{name}"
    fig.savefig(path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("wrote", path)


# ---------------------------------------------------------------- Figure 1
# Already produced whole by sweep_consensus_params.py --plot
src = "varlg/sweepfig/consensus_thresholds.png"
if os.path.exists(src):
    img = mpimg.imread(src)
    fig, ax = plt.subplots(figsize=(14, 4.4))
    panel(ax, img)
    save(fig, "Figure1_thresholds.png")

# ---------------------------------------------------------------- Figure 2
img = load("mockfigs/out_main/ribbon.pdf")
fig, ax = plt.subplots(figsize=(15, 15 * img.shape[0] / img.shape[1]))
panel(ax, img)
save(fig, "Figure2_ribbon.png")

# ---------------------------------------------------------------- Figure 3
hm = load("mockfigs/out_main/character_matrix_heatmap.pdf")
tr = load("mockfigs/out_main/tree_reversible.pdf")
ti = load("mockfigs/out_main/tree_irreversible.pdf")
fig = plt.figure(figsize=(14, 8.5))
gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1], hspace=0.16, wspace=0.10)
panel(fig.add_subplot(gs[0, :]), hm, "A")
panel(fig.add_subplot(gs[1, 0]), tr, "B")
panel(fig.add_subplot(gs[1, 1]), ti, "C")
save(fig, "Figure3_matrix_trees.png")

# ---------------------------------------------------------------- Figure 4
taxa = ["Sarc", "Ixod", "Neoseiu"]
fig = plt.figure(figsize=(14, 8.2))
gs = fig.add_gridspec(2, 3, height_ratios=[1.5, 1], hspace=0.18, wspace=0.08)
for j, sp in enumerate(taxa):
    panel(fig.add_subplot(gs[0, j]), load(f"mockfigs/out_main/dotplots/{sp}.pdf"),
          "ABC"[j])
    panel(fig.add_subplot(gs[1, j]), load(f"mockfigs/out_main/ideograms/{sp}.pdf"),
          "DEF"[j])
save(fig, "Figure4_synteny_detail.png")

# ---------------------------------------------------------------- Figure 5
before = load("mockfigs/out_norescue/dotplots/Tvan.pdf")
after = load("mockfigs/out_rescue/dotplots/Tvan.wgd.pdf")

# Panel C: direct vs rescued per chromosome, from ohnolog_summary.tsv
rows = [l.rstrip("\n").split("\t")
        for l in open("mockfigs/out_rescue/ohnolog_summary.tsv")][1:]
rows = [r for r in rows if r[0] == "Tvan"][:12]
labels = [r[1].replace("Tvan_", "") for r in rows]
direct = [int(r[3]) for r in rows]
rescued = [int(r[4]) for r in rows]

# Panel D: mock paralog-pair support distribution (stands in for duplication_blocks.tsv)
rng = np.random.default_rng(4)
real = rng.integers(28, 95, 12)                 # confirmed ohnologous pairs
noise = rng.integers(1, 11, 140)                # incidental paralogy
support = np.concatenate([real, noise])

fig = plt.figure(figsize=(14, 9))
gs = fig.add_gridspec(2, 2, height_ratios=[1.35, 1], hspace=0.28, wspace=0.16)
panel(fig.add_subplot(gs[0, 0]), before, "A")
panel(fig.add_subplot(gs[0, 1]), after, "B")

axC = fig.add_subplot(gs[1, 0])
x = np.arange(len(labels))
axC.bar(x, direct, 0.72, label="direct RBH", color="#6E9BC5")
axC.bar(x, rescued, 0.72, bottom=direct, label="rescued ohnolog", color="#C0504D")
axC.set_xticks(x)
axC.set_xticklabels(labels, rotation=90, fontsize=6.5)
axC.set_ylabel("Orthologs")
axC.legend(fontsize=8, frameon=False)
axC.spines[["top", "right"]].set_visible(False)
axC.text(-0.08, 1.04, "C", transform=axC.transAxes, fontsize=15, fontweight="bold")

axD = fig.add_subplot(gs[1, 1])
axD.hist(support, bins=np.arange(0, 100, 4), color="#7F7F7F", edgecolor="white")
axD.axvline(12, color="#C00000", linestyle="--", linewidth=1.6)
axD.text(14, axD.get_ylim()[1] * 0.85, "--wgd-min-paralog-pairs = 12",
         color="#C00000", fontsize=8)
axD.set_xlabel("Supporting paralog pairs per chromosome pair")
axD.set_ylabel("Candidate pairs")
axD.spines[["top", "right"]].set_visible(False)
axD.text(-0.08, 1.04, "D", transform=axD.transAxes, fontsize=15, fontweight="bold")

save(fig, "Figure5_wgd_rescue.png")

# ---------------------------------------------------------------- Figure 6
img = load("mockfigs/out_rescue/ribbon.pdf")
fig, ax = plt.subplots(figsize=(15, 15 * img.shape[0] / img.shape[1]))
panel(ax, img)
save(fig, "Figure6_wgd_context.png")
