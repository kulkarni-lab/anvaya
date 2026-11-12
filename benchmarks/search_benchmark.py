#!/usr/bin/env python3
"""
search_benchmark.py - does the heuristic search find the optimum?

Reproduces the test that exposed the weakness of the NNI hill-climber used
before version 1.0: random matrices of 16 terminals x 20 three-state fusion
characters with no phylogenetic structure, for which the optimum can be
found exactly. For each dataset it
reports

  * the exact optimum, by branch and bound (best_tree_bandb), seeded with
    the heuristic score so that it only has to prove optimality;
  * the TBR + ratchet search (the default heuristic from version 1.0);
  * the legacy NNI hill-climber (--search nni), at the old default of 10
    restarts and at 500 restarts;
  * optionally TNT, if --tnt-bin points to a TNT executable: the matrix is
    exported and solved with implicit enumeration (ienum), as in the
    editorial check.

    # editorial design, with TNT as the referee (B&B too slow at 16 taxa in Python)
    python3 benchmarks/search_benchmark.py --datasets 4 --taxa 16 --chars 20 \
        --no-exact --long 100 --tnt-bin tnt --outdir bench16

    # exact check at a size branch and bound can prove in Python
    python3 benchmarks/search_benchmark.py --datasets 50 --taxa 12 --chars 20 \
        --legacy-restarts 10 --outdir bench12

Writes bench_out/search_benchmark.tsv and one TNT/NEXUS matrix per dataset.
Datasets are generated from --seed, so the table is reproducible.
"""
import argparse
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import alg_parsimony as ap  # noqa: E402


def make_dataset(rng, n_taxa, n_chars, p_missing, kind="fusion"):
    """kind='fusion' (default): like real fusion matrices, every character is
    derived (state 1 or 2) in a random 2-6 taxa and unfused elsewhere, with
    no phylogenetic structure; the outgroup (first taxon) is always unfused.
    kind='uniform': every cell drawn uniformly from 0/1/2 - far harder and far
    longer trees, useful as a stress test."""
    taxa = [f"T{i:02d}" for i in range(1, n_taxa + 1)]
    chars = [(f"A{j}", f"B{j}") for j in range(n_chars)]
    matrix = {t: {} for t in taxa}
    for c in chars:
        if kind == "uniform":
            for t in taxa:
                matrix[t][c] = rng.choice((0, 1, 2))
        else:
            derived = set(rng.sample(taxa[1:], rng.randint(2, min(6, n_taxa - 1))))
            for t in taxa:
                matrix[t][c] = rng.choice((1, 2)) if t in derived else 0
        for t in taxa:
            if rng.random() < p_missing:
                matrix[t][c] = "?"
    return taxa, chars, matrix


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--datasets", type=int, default=4)
    p.add_argument("--taxa", type=int, default=16)
    p.add_argument("--chars", type=int, default=20)
    p.add_argument("--missing", type=float, default=0.0)
    p.add_argument("--kind", choices=["fusion", "uniform"], default="fusion")
    p.add_argument("--seed", type=int, default=2026)
    p.add_argument("--criterion", choices=["reversible", "irreversible", "both"], default="both")
    p.add_argument("--replicates", type=int, default=20)
    p.add_argument("--ratchet", type=int, default=10)
    p.add_argument("--legacy-restarts", type=int, nargs="+", default=[10, 500])
    p.add_argument("--long", type=int, default=0,
                   help="also run a second, independent TBR search with this many replicates "
                        "and twice the ratchet iterations, as a convergence check where branch "
                        "and bound is too slow (pure Python B&B is practical to ~14 taxa)")
    p.add_argument("--no-exact", action="store_true", help="skip branch and bound")
    p.add_argument("--tnt-bin", default=None)
    p.add_argument("--outdir", default="bench_out")
    p.add_argument("--threads", type=int, default=1,
                   help="worker processes for the TBR replicates (0 = all CPUs)")
    a = p.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    workers = ap.resolve_workers(a.threads)
    modes = ["reversible", "irreversible"] if a.criterion == "both" else [a.criterion]
    rows = []
    data_rng = random.Random(a.seed)
    for d in range(1, a.datasets + 1):
        taxa, chars, matrix = make_dataset(data_rng, a.taxa, a.chars, a.missing, a.kind)
        og = taxa[0]
        stem = os.path.join(a.outdir, f"dataset{d}")
        ap.write_tnt(stem + ".tnt", taxa, chars, matrix, og, search="ienum")
        ap.write_nexus(stem + ".nex", taxa, chars, matrix, og)
        for mode in modes:
            row = {"dataset": d, "criterion": mode}
            t0 = time.time()
            s, trees, st = ap.tbr_search(taxa, chars, matrix, mode, og, random.Random(a.seed + d),
                                         replicates=a.replicates, ratchet_iter=a.ratchet,
                                         workers=workers)
            row.update(tbr_score=s, tbr_hits=f"{st['hits']}/{st['replicates']}",
                       tbr_seconds=round(time.time() - t0, 1))
            if a.long:
                t0 = time.time()
                s2, _, st2 = ap.tbr_search(taxa, chars, matrix, mode, og, random.Random(a.seed + 1000 + d),
                                           replicates=a.long, ratchet_iter=2 * a.ratchet,
                                           workers=workers)
                row.update(long_tbr_score=s2, long_tbr_hits=f"{st2['hits']}/{st2['replicates']}",
                           long_tbr_seconds=round(time.time() - t0, 1))
            for r in a.legacy_restarts:
                t0 = time.time()
                ls, _ = ap.best_tree_nni(taxa, chars, matrix, mode, og, r, random.Random(a.seed + d))
                row[f"nni{r}_score"] = ls
                row[f"nni{r}_seconds"] = round(time.time() - t0, 1)
            if not a.no_exact:
                t0 = time.time()
                es, etrees = ap.best_tree_bandb(taxa, chars, matrix, mode, og, upper_bound=s,
                                                max_trees=1)
                row.update(exact_score=es,
                           exact_seconds=round(time.time() - t0, 1))
            if a.tnt_bin and mode == "reversible":
                row["tnt_score"] = ap.run_tnt_and_parse(a.tnt_bin, stem + ".tnt", a.outdir)
            rows.append(row)
            print(row, flush=True)

    cols = []
    for r in rows:
        for k in r:
            if k not in cols:
                cols.append(k)
    out = os.path.join(a.outdir, "search_benchmark.tsv")
    with open(out, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
