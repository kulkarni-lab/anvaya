#!/usr/bin/env python3
"""
sweep_consensus_params.py
=========================

Sweeps the three consensus-mode thresholds of genomes_to_rbh.py and writes a
single TSV summarising how many linkage groups each combination yields, so the
choice can be made by looking at data rather than by guessing.

    --consensus-min-species          how many core genomes must contain an orthogroup
    --consensus-min-frac             how much co-location disagreement to tolerate
    --consensus-min-alg-orthogroups  discard reconstructed groups smaller than this

--------------------------------------------------------------------------
WHY THIS IS FASTER THAN A TRIPLE-NESTED LOOP
--------------------------------------------------------------------------
--consensus-min-alg-orthogroups is a pure POST-reconstruction size filter: it
does not change how groups are formed, only which of the already-formed groups
survive. Every run writes consensus_linkage_group_sizes.tsv listing every
candidate group's member-orthogroup count BEFORE that filter is applied, so the
result for every possible threshold can be derived from one reconstruction by
counting. This script therefore runs genomes_to_rbh.py only once per
(min-species x min-frac) combination and derives the whole min-alg-orthogroups
axis arithmetically - turning an O(a*b*c) sweep into O(a*b) invocations.

On top of that, genomes_to_rbh.py caches its diamond alignments per work
directory (keyed by --evalue/--max-target-seqs), and this script points every
invocation at one shared work directory. The expensive alignment step therefore
runs ONCE for the whole sweep; every subsequent combination reuses it. Expect
the first combination to dominate the total runtime.

--------------------------------------------------------------------------
OUTPUT
--------------------------------------------------------------------------
A TSV with one row per (min_species, min_frac, min_alg_orthogroups):

    min_species  min_frac  min_alg_orthogroups  n_linkage_groups
    n_orthogroups_kept  n_orthogroups_in_kept_groups  largest_group  median_group
    pct_orthogroups_retained

Reading it:
  * n_linkage_groups         - the number you are trying to make biologically
                               plausible (compare against known haploid
                               chromosome counts for your clade)
  * pct_orthogroups_retained - what fraction of the orthogroups that survived
                               the presence filter ended up inside a KEPT group.
                               A combination that yields a nice group count but
                               retains very little of the data is a bad trade.
  * largest_group            - a single enormous group swallowing most
                               orthogroups is a warning that min_frac is too
                               permissive and unrelated groups are merging.

--------------------------------------------------------------------------
USAGE
--------------------------------------------------------------------------
    python3 sweep_consensus_params.py \\
        --input-dir ./Data \\
        --consensus-species Ixod Varr Sarc Neoseiu \\
        --outdir ./sweep \\
        --min-species 3 4 \\
        --min-frac 1.0 0.75 0.5 \\
        --min-alg-orthogroups 1 5 10 15 20 30 \\
        --threads 30

Pass any other genomes_to_rbh.py option through unchanged, e.g.:

        --extra --min-scaffold-length 100000
"""

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path
from statistics import median


def log(msg):
    print(f"[sweep] {msg}", flush=True)


def find_script(explicit):
    if explicit:
        p = Path(explicit)
        if not p.exists():
            sys.exit(f"--genomes-to-rbh {explicit} does not exist")
        return str(p)
    here = Path(__file__).resolve().parent / "genomes_to_rbh.py"
    if here.exists():
        return str(here)
    cwd = Path.cwd() / "genomes_to_rbh.py"
    if cwd.exists():
        return str(cwd)
    found = shutil.which("genomes_to_rbh.py")
    if found:
        return found
    sys.exit("Could not find genomes_to_rbh.py (looked next to this script, in the current "
             "directory, and on PATH). Pass --genomes-to-rbh /path/to/genomes_to_rbh.py.")


def read_sizes(path):
    """consensus_linkage_group_sizes.tsv -> list of member-orthogroup counts,
    covering every candidate group BEFORE any size filtering."""
    sizes = []
    with open(path) as f:
        header = f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                try:
                    sizes.append(int(parts[1]))
                except ValueError:
                    pass
    return sizes


def parse_presence(stdout):
    """Pull 'N/M orthogroups present in >= K core species' out of the log."""
    m = re.search(r"(\d+)/(\d+) orthogroups present in >= (\d+) core species", stdout)
    return (int(m.group(1)), int(m.group(2))) if m else (None, None)


def plot_sweep(rows, outdir, out_png, panel_a_sizes, panel_a_label,
                chosen_frac=None, chosen_minalg=None, chosen_minsp=None):
    """Three-panel threshold-selection figure:
      A  rank-ordered linkage-group size distribution at the chosen setting,
         with the size threshold marked
      B  heatmap of n_linkage_groups over min_frac x min_alg_orthogroups
      C  heatmap of pct_orthogroups_retained over the same grid
    Returns True if written, False if matplotlib is unavailable."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        matplotlib.rcParams["pdf.fonttype"] = 42   # editable text, not outlines
        matplotlib.rcParams["ps.fonttype"] = 42
        matplotlib.rcParams["svg.fonttype"] = "none"
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        log("matplotlib not available; skipping --plot")
        return False

    # Restrict the heatmaps to one min_species value so each cell is unambiguous.
    if chosen_minsp is None:
        chosen_minsp = sorted({r["min_species"] for r in rows})[0]
    sub = [r for r in rows if r["min_species"] == chosen_minsp
           and r["n_linkage_groups"] != "NA"]
    if not sub:
        log("no successful combinations to plot")
        return False

    fracs = sorted({r["min_frac"] for r in sub}, reverse=True)
    maos = sorted({r["min_alg_orthogroups"] for r in sub})
    lut = {(r["min_frac"], r["min_alg_orthogroups"]): r for r in sub}

    ngroups = np.full((len(fracs), len(maos)), np.nan)
    retained = np.full((len(fracs), len(maos)), np.nan)
    for i, f in enumerate(fracs):
        for j, m in enumerate(maos):
            r = lut.get((f, m))
            if r is None:
                continue
            ngroups[i, j] = float(r["n_linkage_groups"])
            if r["pct_orthogroups_retained"] != "NA":
                retained[i, j] = float(r["pct_orthogroups_retained"])

    fig = plt.figure(figsize=(15.5, 4.3))
    # Generous wspace: each heatmap carries a colorbar on its right, whose label
    # sits outside the axes and would otherwise collide with the next panel's
    # y-axis label.
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1, 1], wspace=0.55)

    # ---- Panel A: size distribution ----
    axA = fig.add_subplot(gs[0, 0])
    if panel_a_sizes:
        ranks = range(1, len(panel_a_sizes) + 1)
        axA.plot(ranks, panel_a_sizes, color="#1F4E79", linewidth=1.4)
        axA.fill_between(ranks, panel_a_sizes, 0.7, color="#1F4E79", alpha=0.15)
        axA.set_yscale("log")
        axA.set_xlabel("Linkage group (rank)")
        axA.set_ylabel("Member orthogroups")
        if chosen_minalg is not None:
            axA.axhline(chosen_minalg, color="#C00000", linestyle="--", linewidth=1.2)
            n_kept = sum(1 for s in panel_a_sizes if s >= chosen_minalg)
            axA.text(0.97, 0.93,
                     f"threshold = {chosen_minalg}\n\u2192 {n_kept} groups kept",
                     transform=axA.transAxes, ha="right", va="top", fontsize=8,
                     color="#C00000")
        axA.set_title(f"A  Size distribution ({panel_a_label})", loc="left",
                      fontsize=10, fontweight="bold")
    else:
        axA.axis("off")

    def heat(ax, mat, title, cmap, fmt, cbar_label):
        im = ax.imshow(mat, cmap=cmap, aspect="auto")
        ax.set_xticks(range(len(maos)))
        ax.set_xticklabels(maos, fontsize=8)
        ax.set_yticks(range(len(fracs)))
        ax.set_yticklabels([f"{f:g}" for f in fracs], fontsize=8)
        ax.set_xlabel("--consensus-min-alg-orthogroups", fontsize=8.5)
        ax.set_ylabel("--consensus-min-frac", fontsize=8.5)
        ax.set_title(title, loc="left", fontsize=10, fontweight="bold")
        finite = mat[np.isfinite(mat)]
        mid = (finite.max() + finite.min()) / 2 if finite.size else 0
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                if not np.isfinite(mat[i, j]):
                    continue
                ax.text(j, i, format(mat[i, j], fmt), ha="center", va="center",
                        fontsize=7,
                        color="white" if mat[i, j] > mid else "#222222")
        if chosen_frac is not None and chosen_minalg is not None \
                and chosen_frac in fracs and chosen_minalg in maos:
            ax.add_patch(plt.Rectangle(
                (maos.index(chosen_minalg) - 0.5, fracs.index(chosen_frac) - 0.5), 1, 1,
                fill=False, edgecolor="#C00000", linewidth=2.2))
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cb.set_label(cbar_label, fontsize=7.5)
        cb.ax.tick_params(labelsize=7)

    heat(fig.add_subplot(gs[0, 1]), ngroups, "B  Linkage groups recovered",
         "viridis_r", ".0f", "n linkage groups")
    heat(fig.add_subplot(gs[0, 2]), retained, "C  Orthogroups retained",
         "YlGnBu", ".0f", "% retained")

    fig.suptitle(f"Consensus threshold selection (--consensus-min-species = {chosen_minsp})",
                 fontsize=11, y=1.02)
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return True


def main():
    ap = argparse.ArgumentParser(
        description="Sweep genomes_to_rbh.py consensus thresholds and tabulate the resulting "
                    "linkage-group counts, to guide threshold choice.",
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("--input-dir", required=True)
    ap.add_argument("--consensus-species", nargs="+", required=True)
    ap.add_argument("--outdir", default="consensus_sweep",
                    help="directory for the sweep's intermediate runs and the summary TSV "
                         "[consensus_sweep]")
    ap.add_argument("--out-tsv", default=None,
                    help="path for the summary TSV [<outdir>/consensus_sweep_summary.tsv]")
    ap.add_argument("--min-species", nargs="+", type=int, default=None,
                    help="values to sweep [just the default: all core species]")
    ap.add_argument("--min-frac", nargs="+", type=float, default=[1.0, 0.75, 0.5],
                    help="values to sweep [1.0 0.75 0.5]")
    ap.add_argument("--min-alg-orthogroups", nargs="+", type=int,
                    default=[1, 5, 10, 15, 20, 30],
                    help="values to report; derived arithmetically, so extra values here are "
                         "essentially free [1 5 10 15 20 30]")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--genomes-to-rbh", default=None,
                    help="path to genomes_to_rbh.py [auto-detected]")
    ap.add_argument("--keep-runs", action="store_true",
                    help="keep each combination's .rbh output directory (uses more disk; the "
                         "shared work directory with the cached alignments is kept either way)")
    ap.add_argument("--plot", action="store_true",
                    help="also write a three-panel threshold-selection figure: the linkage-group "
                         "size distribution at the chosen setting, plus heatmaps of the recovered "
                         "group count and the orthogroup retention across the sweep grid "
                         "(requires matplotlib)")
    ap.add_argument("--plot-out", default=None,
                    help="path for the figure [<outdir>/consensus_thresholds.png]; the extension "
                         "sets the format, so use .pdf for vector output")
    ap.add_argument("--chosen-min-species", type=int, default=None,
                    help="which --min-species value the figure's heatmaps should show, since each "
                         "cell must be unambiguous [the lowest swept value]")
    ap.add_argument("--chosen-min-frac", type=float, default=None,
                    help="highlight this --min-frac in the figure, and draw panel A from it")
    ap.add_argument("--chosen-min-alg-orthogroups", type=int, default=None,
                    help="highlight this size threshold in the figure and mark it on panel A")
    ap.add_argument("--extra", nargs=argparse.REMAINDER, default=[],
                    help="everything after this is passed through to genomes_to_rbh.py unchanged")
    args = ap.parse_args()

    script = find_script(args.genomes_to_rbh)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    work_dir = outdir / "_shared_work"     # one work dir -> alignments computed once
    out_tsv = Path(args.out_tsv) if args.out_tsv else outdir / "consensus_sweep_summary.tsv"

    n_core = len(args.consensus_species)
    min_species_values = args.min_species if args.min_species is not None else [n_core]
    for v in min_species_values:
        if v > n_core:
            sys.exit(f"--min-species {v} exceeds the number of core species ({n_core})")
        if v < 2:
            sys.exit(f"--min-species {v} is too low; at least 2 core species must share an "
                     f"orthogroup for co-location to be assessable at all")

    combos = [(ms, mf) for ms in min_species_values for mf in args.min_frac]
    log(f"{len(combos)} reconstruction(s) to run "
        f"({len(min_species_values)} min-species x {len(args.min_frac)} min-frac); "
        f"{len(args.min_alg_orthogroups)} size thresholds derived from each without rerunning")
    log(f"alignments are cached in {work_dir} - the first run does the expensive work, "
        f"the rest reuse it")

    rows = []
    for i, (min_sp, min_frac) in enumerate(combos, start=1):
        tag = f"ms{min_sp}_mf{min_frac}"
        run_dir = outdir / f"run_{tag}"
        log(f"[{i}/{len(combos)}] min-species={min_sp} min-frac={min_frac}")

        cmd = [sys.executable, script,
               "--input-dir", args.input_dir,
               "--consensus-species", *args.consensus_species,
               "--consensus-min-species", str(min_sp),
               "--consensus-min-frac", str(min_frac),
               "--rbh-outdir", str(run_dir),
               "--work-dir", str(work_dir),
               "--threads", str(args.threads),
               "--rbh-only"] + list(args.extra)

        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            log(f"  FAILED (exit {proc.returncode}) - recording as failed and continuing")
            tail = (proc.stdout + proc.stderr).strip().splitlines()
            for line in tail[-3:]:
                log(f"    {line}")
            for mao in args.min_alg_orthogroups:
                rows.append(dict(min_species=min_sp, min_frac=min_frac, min_alg_orthogroups=mao,
                                 n_linkage_groups="NA", n_orthogroups_kept="NA",
                                 n_orthogroups_in_kept_groups="NA", largest_group="NA",
                                 median_group="NA", pct_orthogroups_retained="NA"))
            continue

        kept_ogs, total_ogs = parse_presence(proc.stdout)
        sizes_path = run_dir / "consensus_linkage_group_sizes.tsv"
        if not sizes_path.exists():
            log("  no consensus_linkage_group_sizes.tsv written; skipping")
            continue
        sizes = read_sizes(sizes_path)

        for mao in args.min_alg_orthogroups:
            surviving = [s for s in sizes if s >= mao]
            n_in_groups = sum(surviving)
            pct = (100.0 * n_in_groups / kept_ogs) if kept_ogs else None
            rows.append(dict(
                min_species=min_sp,
                min_frac=min_frac,
                min_alg_orthogroups=mao,
                n_linkage_groups=len(surviving),
                n_orthogroups_kept=kept_ogs if kept_ogs is not None else "NA",
                n_orthogroups_in_kept_groups=n_in_groups,
                largest_group=max(surviving) if surviving else 0,
                median_group=int(median(surviving)) if surviving else 0,
                pct_orthogroups_retained=f"{pct:.1f}" if pct is not None else "NA",
            ))
        log(f"  -> {len(sizes)} candidate groups; "
            f"{kept_ogs}/{total_ogs} orthogroups passed the presence filter")

        # Always preserve the sizes file - it is the actual per-combination result and
        # what panel A of the threshold figure is drawn from - regardless of whether the
        # bulky .rbh output directory is being kept.
        shutil.copy(sizes_path, outdir / f"sizes_{tag}.tsv")
        if not args.keep_runs:
            shutil.rmtree(run_dir, ignore_errors=True)

    if not rows:
        sys.exit("No results produced - check the errors above.")

    cols = ["min_species", "min_frac", "min_alg_orthogroups", "n_linkage_groups",
            "n_orthogroups_kept", "n_orthogroups_in_kept_groups", "largest_group",
            "median_group", "pct_orthogroups_retained"]
    with open(out_tsv, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")

    log(f"wrote {out_tsv} ({len(rows)} rows)")

    if args.plot:
        chosen_sp = args.chosen_min_species
        if chosen_sp is None:
            chosen_sp = sorted({r["min_species"] for r in rows})[0]
        chosen_frac = args.chosen_min_frac
        if chosen_frac is None:
            chosen_frac = sorted({r["min_frac"] for r in rows})[0]
        panel_a_sizes, panel_a_label = [], ""
        sizes_file = outdir / f"sizes_ms{chosen_sp}_mf{chosen_frac}.tsv"
        if sizes_file.exists():
            panel_a_sizes = sorted(read_sizes(sizes_file), reverse=True)
            panel_a_label = f"min-species {chosen_sp}, min-frac {chosen_frac:g}"
        else:
            log(f"  (no sizes file for min-species={chosen_sp} min-frac={chosen_frac}; "
                f"panel A will be blank)")
        plot_path = Path(args.plot_out) if args.plot_out else outdir / "consensus_thresholds.png"
        if plot_sweep(rows, outdir, plot_path, panel_a_sizes, panel_a_label,
                       chosen_frac=chosen_frac,
                       chosen_minalg=args.chosen_min_alg_orthogroups,
                       chosen_minsp=chosen_sp):
            log(f"wrote {plot_path}")
    log("")
    log("How to read it: pick a combination whose n_linkage_groups is biologically")
    log("plausible for your clade AND whose pct_orthogroups_retained is not tiny. Then")
    log("re-run the DOWNSTREAM analysis at 2-3 nearby combinations and check whether your")
    log("conclusions (topology, supported clades) are stable. If the size distribution")
    log("decays smoothly with no natural break, no threshold is objectively correct -")
    log("report the one you chose as an explicit analytical decision, plus that")
    log("stability check.")
    log("")
    log("Warning signs: n_linkage_groups collapsing toward 1 as min_frac drops, or")
    log("largest_group swallowing most of n_orthogroups_in_kept_groups, both mean")
    log("min_frac is too permissive and unrelated groups are being merged.")


if __name__ == "__main__":
    main()
