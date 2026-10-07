#!/usr/bin/env python3
"""
calibrate_wgd.py
================

Per-taxon calibration of --wgd-min-paralog-pairs (and --wgd-alpha) from a
duplication_blocks.tsv written by genomes_to_rbh.py run with
--detect-wgd --wgd-min-paralog-pairs 1, i.e. with the count filter switched off
so that every candidate chromosome pair is reported.

Paralog-pair counts scale with genome size and gene content, so a single
threshold read off the pooled distribution is dominated by whichever genome is
largest. This tool instead applies a statistical criterion to each taxon
separately: each taxon's Fisher's exact p-values are Bonferroni-corrected for
the number of chromosome pairs it could have been tested on, and the suggested
threshold is the lowest paralog-pair count among the pairs that survive. It also
reports how many chromosomes have more than one significant partner - for a
single clean duplication each chromosome should have one, so a connected
network indicates post-duplication rearrangement or more than one round.

USAGE
    python3 calibrate_wgd.py rbh_wgd_diag/duplication_blocks.tsv
    python3 calibrate_wgd.py rbh_wgd_diag/duplication_blocks.tsv \\
        --plot --plot-out FigureS_wgd_thresholds.pdf \\
        --chosen Limulus=40 Luthela=30 Mgig=30

The printed suggestion is ready to paste into genomes_to_rbh.py as
    --wgd-min-paralog-pairs Limulus=.. Luthela=.. Mgig=..
"""

import argparse
import math
import sys
from collections import defaultdict, Counter


def read_blocks(path):
    rows = defaultdict(list)
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 5 or not f[3].strip().isdigit():
            continue                              # header or malformed line
        try:
            p = float(f[4])
        except ValueError:
            p = float("nan")                      # 'None'/'NA': no background to test against
        rows[f[0]].append((f[1], f[2], int(f[3]), p))
    for sp in rows:
        rows[sp].sort(key=lambda r: -r[2])
    return rows


def calibrate(rs, alpha=0.05):
    chroms = {c for r in rs for c in r[:2]}
    n_tests = max(len(chroms) * (len(chroms) - 1) // 2, 1)
    bonf = alpha / n_tests
    sig = [r for r in rs if r[3] == r[3] and r[3] <= bonf]
    # The Fisher test asks whether a chromosome PAIR is enriched for paralogy
    # relative to the genome-wide background. On real genomes that is not the same
    # question as "are these two chromosomes ohnologous": two small chromosomes
    # sharing a handful of paralogs can be significantly enriched, so significance
    # alone will happily accept pairs with 1-3 supporting paralogs. Significance is
    # therefore necessary but nowhere near sufficient, and a threshold read off it
    # is meaningless on its own.
    #
    # Karyotype expectation is the usable constraint. After ONE whole-genome
    # duplication every chromosome has exactly one ohnologous partner, so roughly
    # n/2 pairs are expected; post-duplication rearrangement or further rounds push
    # that up, but a count several times larger than n is not ohnology, it is noise.
    # The suggestion is therefore the lowest threshold whose surviving pair count
    # falls within that expectation.
    n_chr = len(chroms)
    targets = {"single_wgd": max(n_chr // 2, 1), "network": n_chr}
    counts = sorted({r[2] for r in rs}, reverse=True)
    sugg = {}
    for label, target in targets.items():
        pick = None
        for t in counts:
            if sum(1 for r in rs if r[2] >= t) > target:
                break
            pick = t
        sugg[label] = pick
    bonf_floor = min((r[2] for r in sig), default=None)
    return {"chroms": chroms, "n_tests": n_tests, "bonf": bonf, "sig": sig,
            "n_chr": n_chr, "suggested": sugg["single_wgd"],
            "suggested_network": sugg["network"], "bonf_floor": bonf_floor}


def partner_stats(rs, threshold):
    kept = [r for r in rs if r[2] >= threshold]
    partners = Counter(c for r in kept for c in r[:2])
    return len(kept), sum(1 for n in partners.values() if n > 1)


def plot(rows, cal, chosen, out):
    import matplotlib
    matplotlib.use("Agg")
    matplotlib.rcParams["pdf.fonttype"] = 42      # editable text, not outlines
    matplotlib.rcParams["ps.fonttype"] = 42
    matplotlib.rcParams["svg.fonttype"] = "none"
    import matplotlib.pyplot as plt

    taxa = list(rows)
    fig, axes = plt.subplots(len(taxa), 3, figsize=(15.5, 3.9 * len(taxa)),
                             gridspec_kw={"width_ratios": [1.1, 1, 1], "wspace": 0.42,
                                          "hspace": 0.55}, squeeze=False)
    for i, sp in enumerate(taxa):
        rs, c = rows[sp], cal[sp]
        thr = chosen.get(sp, c["suggested"])
        counts = [r[2] for r in rs]
        is_sig = [r[3] == r[3] and r[3] <= c["bonf"] for r in rs]

        # A: rank-ordered counts
        ax = axes[i][0]
        ranks = range(1, len(counts) + 1)
        ax.plot(ranks, counts, color="#7F7F7F", linewidth=1.0, zorder=1)
        ax.scatter([k + 1 for k, s in enumerate(is_sig) if s],
                   [counts[k] for k, s in enumerate(is_sig) if s],
                   s=14, color="#C0504D", zorder=2, label="Bonferroni-significant")
        ax.scatter([k + 1 for k, s in enumerate(is_sig) if not s],
                   [counts[k] for k, s in enumerate(is_sig) if not s],
                   s=8, color="#9DB7D5", zorder=2, label="not significant")
        ax.set_yscale("log")
        if thr:
            ax.axhline(thr, color="#C00000", linestyle="--", linewidth=1.2)
            n_kept, _ = partner_stats(rs, thr)
            ax.text(0.97, 0.93, f"threshold = {thr}\n\u2192 {n_kept} pairs confirmed",
                    transform=ax.transAxes, ha="right", va="top", fontsize=8, color="#C00000",
                    bbox=dict(facecolor="white", edgecolor="none", pad=1.5, alpha=0.85))
        ax.set_xlabel("Candidate chromosome pair (rank)")
        ax.set_ylabel("Supporting paralog pairs")
        ax.legend(fontsize=7, frameon=False, loc="lower left")
        ax.set_title(f"{'ABC'[0]}{i + 1}  {sp}: support distribution", loc="left",
                     fontsize=10, fontweight="bold")

        # B: count vs significance
        ax = axes[i][1]
        xs = [r[2] for r in rs if r[3] == r[3] and r[3] > 0]
        ys = [-math.log10(r[3]) for r in rs if r[3] == r[3] and r[3] > 0]
        ax.scatter(xs, ys, s=10, color="#1F4E79", alpha=0.7)
        ax.axhline(-math.log10(c["bonf"]), color="#7F7F7F", linestyle=":", linewidth=1.1)
        ax.text(0.03, -math.log10(c["bonf"]), f"  Bonferroni \u03b1 = {c['bonf']:.1e}",
                transform=ax.get_yaxis_transform(), va="bottom", fontsize=7, color="#555555")
        if thr:
            ax.axvline(thr, color="#C00000", linestyle="--", linewidth=1.2)
        ax.set_xscale("log")
        ax.set_xlabel("Supporting paralog pairs")
        ax.set_ylabel("\u2212log\u2081\u2080 Fisher's p")
        ax.set_title(f"{'ABC'[1]}{i + 1}  Support vs. significance", loc="left",
                     fontsize=10, fontweight="bold")

        # C: consequences across thresholds
        ax = axes[i][2]
        grid = sorted(set([1, 2, 3, 5, 8, 10, 15, 20, 25, 30, 40, 50, 60, 80, 100]
                          + ([thr] if thr else [])))
        grid = [g for g in grid if g <= max(counts)]
        kept = [partner_stats(rs, g)[0] for g in grid]
        multi = [partner_stats(rs, g)[1] for g in grid]
        ax.plot(grid, kept, "-o", color="#1F4E79", markersize=3.5, label="confirmed pairs")
        ax.plot(grid, multi, "-s", color="#C0504D", markersize=3.5,
                label="chromosomes with >1 partner")
        if thr:
            ax.axvline(thr, color="#C00000", linestyle="--", linewidth=1.2)
        ax.set_xscale("log")
        ax.set_xlabel("--wgd-min-paralog-pairs")
        ax.set_ylabel("Count")
        ax.legend(fontsize=7, frameon=False)
        ax.set_title(f"{'ABC'[2]}{i + 1}  Consequence of the threshold", loc="left",
                     fontsize=10, fontweight="bold")

    fig.suptitle("Calibration of --wgd-min-paralog-pairs per taxon", fontsize=12, y=1.0)
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[1],
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog=__doc__)
    ap.add_argument("blocks", help="duplication_blocks.tsv from a --wgd-min-paralog-pairs 1 run")
    ap.add_argument("--alpha", type=float, default=0.05,
                    help="family-wise alpha before Bonferroni correction [0.05]")
    ap.add_argument("--plot", action="store_true", help="write the calibration figure")
    ap.add_argument("--plot-out", default="wgd_thresholds.pdf",
                    help="figure path; extension sets format [wgd_thresholds.pdf]")
    ap.add_argument("--chosen", nargs="+", default=[],
                    help="thresholds to mark as SPECIES=N [default: the suggested value]")
    args = ap.parse_args()

    rows = read_blocks(args.blocks)
    if not rows:
        sys.exit(f"no candidate pairs read from {args.blocks}")
    chosen = {}
    for v in args.chosen:
        sp, n = v.split("=", 1)
        chosen[sp] = int(n)
    unknown = set(chosen) - set(rows)
    if unknown:
        sys.exit(f"--chosen names {sorted(unknown)}, not present in {args.blocks}: {sorted(rows)}")

    cal = {}
    for sp, rs in rows.items():
        c = cal[sp] = calibrate(rs, args.alpha)
        _, multi = partner_stats(rs, c["suggested"] or 1)
        print(f"=== {sp} ===")
        print(f"  chromosomes: {len(c['chroms'])}   candidate pairs tested: {c['n_tests']}   "
              f"Bonferroni alpha: {c['bonf']:.2e}")
        print(f"  Bonferroni-significant pairs: {len(c['sig'])} "
              f"(lowest count among them: {c['bonf_floor']} \u2014 significance alone is NOT a "
              f"usable threshold; see below)")
        n_at = sum(1 for r in rs if c["suggested"] and r[2] >= c["suggested"])
        n_at_net = sum(1 for r in rs if c["suggested_network"] and r[2] >= c["suggested_network"])
        print(f"  karyotype expectation: {len(c['chroms'])} chromosomes \u2192 ~{len(c['chroms'])//2} "
              f"ohnologous pairs after one WGD, up to ~{len(c['chroms'])} if rearranged/multi-round")
        print(f"  suggested --wgd-min-paralog-pairs: {c['suggested']} "
              f"({n_at} pairs, single-WGD expectation) "
              f"or {c['suggested_network']} ({n_at_net} pairs, permissive)")
        print(f"  chromosomes with >1 significant partner: {multi} "
              f"({'network - post-WGD rearrangement or >1 round' if multi else 'clean pairwise'})")
        print()

    sugg = " ".join(f"{sp}={cal[sp]['suggested']}" for sp in rows if cal[sp]["suggested"])
    if sugg:
        print(f"Suggested: --wgd-min-paralog-pairs {sugg}")
    if args.plot:
        plot(rows, cal, chosen, args.plot_out)
        print(f"wrote {args.plot_out}")


if __name__ == "__main__":
    main()
