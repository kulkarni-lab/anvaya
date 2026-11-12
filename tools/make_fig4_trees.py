#!/usr/bin/env python3
"""Figure 4C-D: most-parsimonious trees for the Acari dataset with bootstrap support.

Values are taken from v101b/Fig4/report.txt (Anvaya v1.0.1b,
1,000 replicates per criterion) and, for panel C, from the TNT v1.5 bootstrap
of the same matrix under the same constraint (examples/tnt/script_anvaya.run).
Writes Figure4CD.pdf and Figure4CD.svg (vector, editable text).

    python3 tools/make_fig4_trees.py [outdir]
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["font.family"] = "DejaVu Sans"
import matplotlib.pyplot as plt

NAMES = {
    "Dsil9745": "D. silvarum", "Dand5885": "D. andersoni", "Dret6085": "D. reticulatus",
    "Hyas": "H. asiaticum", "Ixod": "I. scapularis", "Avul": "A. vulgaris",
    "Otur": "O. turicata", "Neoseiu": "N. longispinosus", "Phytos": "P. persimilis",
    "Varr": "V. destructor", "Blom": "B. tropicalis", "Sarc": "S. scabiei",
    "Espinosa": "E. spinosa",
}

# Topology shared by both criteria (first child is drawn above the second).
TREE = (
    (
        (
            ((((("Dsil9745", "Dand5885"), "Dret6085"), "Hyas"), "Ixod"), ("Avul", "Otur")),
            (("Neoseiu", "Phytos"), "Varr"),
        ),
        ("Blom", "Sarc"),
    ),
    "Espinosa",
)


def leaves(t):
    return [t] if isinstance(t, str) else [x for c in t for x in leaves(c)]


def key(*taxa):
    return frozenset(taxa)


IXODIDAE = ("Dand5885", "Dret6085", "Dsil9745", "Hyas", "Ixod")
METASTRIATA = ("Dand5885", "Dret6085", "Dsil9745", "Hyas")
DERMACENTOR = ("Dand5885", "Dret6085", "Dsil9745")
COMMON = {
    key(*IXODIDAE, "Avul", "Otur", "Neoseiu", "Phytos", "Varr"): 100,   # Parasitiformes
    key("Neoseiu", "Phytos", "Varr"): 100,                              # Mesostigmata
    key("Neoseiu", "Phytos"): 100,
    key(*IXODIDAE, "Avul", "Otur"): 100,                                # Ixodida
    key("Avul", "Otur"): 100,                                           # Argasidae
    key(*IXODIDAE): 100,                                                # Ixodidae
    key("Blom", "Sarc"): 100,                                           # Acariformes
}
# Anvaya v1.0.1b (report.txt)
REV = {**COMMON, **{key(*METASTRIATA): 100, key(*DERMACENTOR): 97, key("Dand5885", "Dsil9745"): 98}}
IRR = {**COMMON, **{key(*METASTRIATA): 100, key(*DERMACENTOR): 98, key("Dand5885", "Dsil9745"): 97}}
# TNT v1.5, reversible criterion (values shown in the previous Figure 4C)
TNT = {**COMMON, **{key(*METASTRIATA): 100, key(*DERMACENTOR): 94, key("Dand5885", "Dsil9745"): 97}}

LINE = "#3a3a3a"


def draw(ax, title, labels, second=None):
    order = leaves(TREE)
    ypos = {t: len(order) - 1 - i for i, t in enumerate(order)}

    def depth(t):
        return 0 if isinstance(t, str) else 1 + max(depth(c) for c in t)

    maxd = depth(TREE)
    tip_x = maxd + 0.6

    def walk(t, d):
        if isinstance(t, str):
            return tip_x, ypos[t]
        x = d
        kids = [walk(c, d + 1) for c in t]
        ys = [ky for _, ky in kids]
        y = (min(ys) + max(ys)) / 2
        ax.plot([x, x], [min(ys), max(ys)], color=LINE, lw=2.6, solid_capstyle="butt")
        for kx, ky in kids:
            ax.plot([x, kx], [ky, ky], color=LINE, lw=2.6, solid_capstyle="butt")
        k = frozenset(leaves(t))
        if k in labels:
            text = str(labels[k]) if second is None else f"{labels[k]} | {second[k]}"
            ax.text(x + 0.1, y, text, fontsize=8.5, ha="left", va="center", color="#111111",
                    bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="none"), zorder=3)
        return x, y

    rx, ry = walk(TREE, 0)
    ax.plot([rx - 0.35, rx], [ry, ry], color=LINE, lw=2.6, solid_capstyle="butt")
    for t, y in ypos.items():
        ax.text(tip_x + 0.12, y, NAMES[t], fontsize=12, style="italic", va="center", ha="left")
    ax.set_title(title, loc="left", fontsize=12.5, pad=6)
    ax.set_xlim(-0.6, tip_x + 3.2)
    ax.set_ylim(-0.7, len(order) - 0.3)
    ax.axis("off")


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else "."
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.4))
    draw(axes[0], "Reversible (fusion and mixing can revert)\n(score = 459)", REV, TNT)
    draw(axes[1], "Irreversible (fusion and mixing cannot revert)\n(score = 468)", IRR)
    for ax, letter in zip(axes, "CD"):
        ax.text(-0.02, 1.12, letter, transform=ax.transAxes, fontsize=22, ha="left", va="bottom")
    axes[0].text(0.0, -0.04, "Node values: Anvaya | TNT bootstrap (%)", transform=axes[0].transAxes,
                 fontsize=8.5, color="#555555")
    axes[1].text(0.0, -0.04, "Node values: Anvaya bootstrap (%)", transform=axes[1].transAxes,
                 fontsize=8.5, color="#555555")
    fig.tight_layout()
    for ext in ("pdf", "svg"):
        fig.savefig(os.path.join(outdir, f"Figure4CD.{ext}"), bbox_inches="tight")
    print("wrote", os.path.join(outdir, "Figure4CD.{pdf,svg}"))


if __name__ == "__main__":
    main()
