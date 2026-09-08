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

        if not args.keep_runs:
            # keep the sizes file (it is the actual result), drop the bulky .rbh output
            keep = outdir / f"sizes_{tag}.tsv"
            shutil.copy(sizes_path, keep)
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
