#!/usr/bin/env python3
"""Anvaya logo. The mark is the character the toolkit codes: two ancestral
linkage groups (teal, amber) that fuse and then interleave, read left to
right as unfused -> fused-clean -> fused-and-mixed. Vector, editable text."""
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib.path import Path
from matplotlib.patches import PathPatch

TEAL, AMBER, INK = "#1F7A6F", "#D08C1E", "#1A2B33"


def chrom(ax, x, y, w, h, segs):
    """segs: list of (fraction, colour) drawn left to right inside a rounded bar."""
    cx = x
    for frac, col in segs:
        ax.add_patch(plt.Rectangle((cx, y), w * frac, h, facecolor=col,
                                   edgecolor="none", zorder=2))
        cx += w * frac
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                                boxstyle="round,pad=0,rounding_size=%.3f" % (h / 2),
                                facecolor="none", edgecolor=INK, linewidth=1.6,
                                zorder=3, mutation_aspect=1))


def ribbon(ax, x0, y0, x1, y1, w, col, alpha=0.30):
    dy = y0 - y1
    c0, c1 = y0 - dy * 0.5, y1 + dy * 0.5
    verts = [(x0 - w, y0), (x0 - w, c0), (x1 - w, c1), (x1 - w, y1),
             (x1 + w, y1), (x1 + w, c1), (x0 + w, c0), (x0 + w, y0), (x0 - w, y0)]
    codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.LINETO,
             Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.CLOSEPOLY]
    ax.add_patch(PathPatch(Path(verts, codes), facecolor=col, edgecolor="none",
                           alpha=alpha, zorder=1))


def build(with_wordmark=True, path="anvaya_logo"):
    if with_wordmark:
        fig, ax = plt.subplots(figsize=(4.6, 2.05))
        ax.set_xlim(0, 4.6); ax.set_ylim(0, 2.05)
    else:
        fig, ax = plt.subplots(figsize=(2.4, 2.4))
        ax.set_xlim(-0.05, 2.45); ax.set_ylim(-0.05, 2.45)
    ax.axis("off")
    ax.set_aspect("equal")

    ox, oy = (0.18, 0.42) if with_wordmark else (0.22, 0.55)
    H, W = 0.20, 1.86

    # top: two separate ancestral linkage groups
    chrom(ax, ox, oy + 1.12, W * 0.45, H, [(1.0, TEAL)])
    chrom(ax, ox + W * 0.55, oy + 1.12, W * 0.45, H, [(1.0, AMBER)])
    # middle: fused, still two contiguous blocks
    chrom(ax, ox, oy + 0.56, W, H, [(0.5, TEAL), (0.5, AMBER)])
    # bottom: fused and mixed
    chrom(ax, ox, oy, W, H,
          [(0.13, TEAL), (0.09, AMBER), (0.10, TEAL), (0.07, AMBER), (0.08, TEAL),
           (0.12, AMBER), (0.07, TEAL), (0.11, AMBER), (0.08, TEAL), (0.15, AMBER)])

    for frac, col in ((0.22, TEAL), (0.78, AMBER)):
        ribbon(ax, ox + W * frac, oy + 1.12, ox + W * (frac * 0.5 + 0.25),
               oy + 0.76, W * 0.075, col)
    for a, b, col in ((0.25, 0.30, TEAL), (0.75, 0.66, AMBER)):
        ribbon(ax, ox + W * a, oy + 0.56, ox + W * b, oy + 0.20, W * 0.06, col)

    if with_wordmark:
        ax.text(2.24, 1.13, "Anvaya", fontsize=31, fontweight="bold",
                color=INK, va="center", ha="left", family="DejaVu Sans")
        ax.text(2.28, 0.66, "macrosynteny as phylogenetic characters",
                fontsize=8.2, color="#5A6B73", va="center", ha="left",
                family="DejaVu Sans")

    for ext, kw in ((".pdf", {}), (".svg", {}), (".png", {"dpi": 400, "transparent": True})):
        fig.savefig(path + ext, bbox_inches="tight", **kw)
    plt.close(fig)
    print("wrote", path + ".{pdf,svg,png}")


build(True, "/home/claude/anvaya_logo")
build(False, "/home/claude/anvaya_mark")
