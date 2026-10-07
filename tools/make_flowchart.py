#!/usr/bin/env python3
"""Figure 1: Anvaya workflow. Vector PDF with TrueType (editable) text."""
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["font.family"] = "DejaVu Sans"
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

C_IN, C_G2R, C_ALG, C_OUT, C_CAL = "#FCE9D6", "#DCEBDA", "#DCE6F2", "#FFF2CC", "#EDE3F1"
EDGE = "#3A3A3A"

fig, ax = plt.subplots(figsize=(11.2, 12.6))
ax.set_xlim(0, 11.1); ax.set_ylim(0, 12.6); ax.axis("off")


def box(x, y, w, h, title, lines=(), fc="#FFFFFF", tsize=10.2, lsize=8.3, bold=True):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.12",
                                linewidth=1.1, edgecolor=EDGE, facecolor=fc))
    ax.text(x + 0.16, y + h - 0.2, title, ha="left", va="top", fontsize=tsize,
            fontweight="bold" if bold else "normal", color="#1A1A1A")
    for k, ln in enumerate(lines):
        ax.text(x + 0.16, y + h - 0.52 - k * 0.27, ln, ha="left", va="top",
                fontsize=lsize, color="#2B2B2B")


def arrow(x1, y1, x2, y2, dashed=False, color=EDGE):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=13,
                                 linewidth=1.2, color=color,
                                 linestyle=(0, (4, 3)) if dashed else "-"))


def band(y, h, label, fc):
    ax.add_patch(FancyBboxPatch((0.15, y), 7.55, h, boxstyle="round,pad=0.02,rounding_size=0.18",
                                linewidth=0, facecolor=fc, alpha=0.55, zorder=0))
    ax.text(0.32, y + h - 0.16, label, ha="left", va="top", fontsize=11,
            fontweight="bold", color="#333333", style="italic")


# ---------------- inputs ----------------
box(0.3, 11.35, 3.55, 1.05, "Genome assemblies",
    ["CDS + GFF3/GTF (+ protein) per taxon"], fc=C_IN)
box(4.05, 11.35, 3.55, 1.05, "Pre-built .rbh tables",
    ["from an existing ALG reference", "(e.g. BCnS ALGs, odp output)"], fc=C_IN)

# ---------------- genomes_to_rbh.py ----------------
band(6.05, 5.05, "genomes_to_rbh.py", C_G2R)
box(0.35, 9.75, 7.2, 0.95, "1  Assembly quality",
    ["scaffold count, span, N50 reported per genome; warning below N50 1 Mb;",
     "--min-scaffold-length removes short-scaffold genes before the search"], fc="white")
box(0.35, 8.65, 7.2, 0.95, "2  Orthology",
    ["DIAMOND bidirectional reciprocal best hits; alignments cached per",
     "--evalue / --max-target-seqs, so re-runs and sweeps skip alignment"], fc="white")
box(0.35, 7.3, 7.2, 1.2, "3  Linkage-group reference",
    ["--reference-species: one genome's chromosomes define the groups, or",
     "--consensus-species: groups reconstructed de novo from core genomes",
     "(--consensus-min-species, -min-frac, -min-alg-orthogroups)"], fc="white")
box(0.35, 6.2, 7.2, 0.95, "4  Duplication detection and ohnolog rescue",
    ["--detect-wgd: self-comparison finds ohnologous chromosome pairs;",
     "per-species --wgd-min-paralog-pairs; rescued hits flagged wgd_inferred"], fc="white")
for y1, y2 in ((9.75, 9.6), (8.65, 8.5), (7.3, 7.15)):
    arrow(3.95, y1, 3.95, y2)
arrow(3.2, 11.35, 3.2, 10.7)

# ---------------- intermediate ----------------
box(0.35, 5.0, 7.2, 0.85, ".rbh ortholog tables",
    ["+ consensus_linkage_group_sizes.tsv, duplication_blocks.tsv"], fc="#F4F4F4")
arrow(3.95, 6.2, 3.95, 5.85)
# pre-built .rbh bypass: routed down the gutter so it crosses no step
ax.plot([7.3, 7.3, 7.85, 7.85], [11.35, 11.2, 11.2, 5.42], color=EDGE, linewidth=1.2)
arrow(7.85, 5.42, 7.58, 5.42)
ax.text(7.93, 10.1, "bypass", rotation=90, ha="left", va="center", fontsize=7.5,
        color="#555555", style="italic")

# ---------------- alg_parsimony.py ----------------
band(0.95, 3.85, "alg_parsimony.py", C_ALG)
box(0.35, 3.45, 7.2, 1.0, "5  Character calling",
    ["per ALG pair: 0 unfused, 1 fused-clean, 2 fused-mixed, ? missing",
     "(--min-genes, --alpha, --mixing-threshold); compound-fusion report"], fc="white")
box(0.35, 2.3, 7.2, 1.0, "6  Parsimony search and ancestral states",
    ["reversible (Fitch) and irreversible (Camin-Sokal) criteria; exhaustive",
     "or NNI search; backbone --constraint; synapomorphies mapped to branches"], fc="white")
box(0.35, 1.15, 7.2, 1.0, "7  Diagnostics",
    ["CI / RI per character; bootstrap (support written into .nwk);",
     "permutation signal test; branch rate models (--rate-model)"], fc="white")
arrow(3.95, 5.0, 3.95, 4.45)
arrow(3.95, 3.45, 3.95, 3.3)
arrow(3.95, 2.3, 3.95, 2.15)

# ---------------- outputs ----------------
box(0.35, 0.05, 3.5, 0.78, "Trees and matrices",
    ["NEXUS / TNT export"], fc=C_OUT, tsize=9.5, lsize=8)
box(4.05, 0.05, 3.5, 0.78, "Figures (vector, editable text)",
    ["dot plots, ideograms, ribbon, heatmap"], fc=C_OUT, tsize=9.5, lsize=8)
arrow(2.1, 1.15, 2.1, 0.83)
arrow(5.8, 1.15, 5.8, 0.83)

# ---------------- calibration side panel ----------------
box(8.2, 7.3, 2.75, 1.95, "sweep_consensus_params.py",
    ["sweeps the three consensus", "thresholds; one TSV plus a", "threshold-selection figure",
     "(Fig. 6)"], fc=C_CAL, tsize=8.3, lsize=8)
box(8.2, 4.95, 2.75, 1.95, "calibrate_wgd.py",
    ["per-taxon support vs.", "karyotype expectation", "(~n/2 pairs after one WGD);",
     "calibration figure (Fig. S1)"], fc=C_CAL, tsize=8.3, lsize=8)
box(8.2, 2.3, 2.75, 2.2, "Threshold sensitivity",
    ["re-run alg_parsimony.py", "across --min-genes and", "--alpha; compare CI/RI,",
     "character count and", "topology"], fc=C_CAL, tsize=8.3, lsize=8)
arrow(8.2, 8.2, 7.6, 7.95, dashed=True, color="#6A4C93")
arrow(8.2, 5.9, 7.6, 6.65, dashed=True, color="#6A4C93")
arrow(8.2, 3.4, 7.6, 3.95, dashed=True, color="#6A4C93")
ax.text(9.57, 9.55, "Calibration", ha="center", fontsize=11, fontweight="bold",
        style="italic", color="#4B3673")
ax.text(9.57, 1.95, "dashed arrows: thresholds\nchosen from these outputs",
        ha="center", va="top", fontsize=7.6, color="#4B3673")

fig.savefig("/home/claude/Figure1_workflow.pdf", bbox_inches="tight")
fig.savefig("/home/claude/Figure1_workflow.png", dpi=220, bbox_inches="tight")
print("written")
