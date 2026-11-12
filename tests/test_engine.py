#!/usr/bin/env python3
"""
Unit checks for the fast scorer and the heuristic search.

    python3 tests/test_engine.py            # quick (~1 min)
    python3 tests/test_engine.py --full     # adds the exhaustive-optimum check at 9 taxa

1. The bit-parallel scorer (fast_score) equals the reference Sankoff scorer
   (score_tree) on random trees and random matrices with missing data, under
   both criteria and for every possible outgroup.
2. On random matrices small enough for exact enumeration, the TBR + ratchet
   search finds the exact optimum, with and without a backbone constraint.
3. Searches are reproducible across Python processes with different
   PYTHONHASHSEED values (see also tests/test_determinism.sh).
"""
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import alg_parsimony as ap  # noqa: E402


def random_matrix(rng, n, m, p_missing=0.1):
    taxa = [f"t{i:02d}" for i in range(n)]
    chars = [(f"A{j}", f"B{j}") for j in range(m)]
    matrix = {}
    for t in taxa:
        matrix[t] = {}
        for c in chars:
            if rng.random() < p_missing:
                matrix[t][c] = "?"
            else:
                matrix[t][c] = rng.choice((0, 0, 1, 2))
    return taxa, chars, matrix


def check_scorer(n_trials=300):
    rng = random.Random(7)
    for trial in range(n_trials):
        n = rng.randint(3, 12)
        taxa, chars, matrix = random_matrix(rng, n, rng.randint(1, 25))
        g = ap.random_tree(taxa, rng)
        fm = ap.FastMatrix(taxa, chars, matrix)
        for mode in ("reversible", "irreversible"):
            for og in taxa:
                ref = ap.score_tree(ap.root_at_leaf(g, og), chars, matrix, mode)
                fast = ap.fast_score_graph(g, taxa, fm, mode, og)
                assert ref == fast, (trial, mode, og, ref, fast)
    print(f"ok   fast scorer == Sankoff reference on {n_trials} random trees x every outgroup")


def check_search_exact(n_sets=12, n_taxa=8, full=False):
    rng = random.Random(11)
    sizes = [n_taxa] * n_sets + ([9] * 4 if full else [])
    for k, n in enumerate(sizes):
        taxa, chars, matrix = random_matrix(rng, n, 20)
        og = taxa[0]
        constraint = None
        if k % 3 == 2:
            root = ap.parse_newick(f"(({taxa[1]},{taxa[2]}),{taxa[3]},{taxa[4]});")
            clades = ap.required_clades_from_constraint(root)
            constraint = {"root_node": root, "taxa": set(root.leaves()),
                          "required_clades": {c for c in clades if len(c) < n}}
        for mode in ("reversible", "irreversible"):
            exact, exact_trees = ap.best_tree_exhaustive(taxa, chars, matrix, mode, og,
                                                          constraint=constraint)
            heur, trees, stats = ap.tbr_search(taxa, chars, matrix, mode, og, random.Random(k),
                                               constraint=constraint, replicates=5,
                                               ratchet_iter=5, hold=200)
            assert heur == exact, (k, mode, heur, exact)
            if constraint:
                for g in trees:
                    assert ap.satisfies_constraints(g, constraint["required_clades"],
                                                    constraint["taxa"])
            exact_keys = {frozenset(ap.bipartitions(t)) for t in exact_trees}
            for g in trees:
                assert frozenset(ap.bipartitions(g)) in exact_keys
    print(f"ok   TBR+ratchet search reaches the exhaustive optimum on {len(sizes)} random "
          f"matrices x 2 criteria (every returned tree is a true MPT)")


def check_irreversible_offset(n_trials=200):
    """The TNT/PAUP* emulation of the irreversible criterion (all-zero ancestor
    attached on the outgroup branch, root state free) must differ from Anvaya's
    score by irreversible_offset(), independently of topology."""
    rng = random.Random(3)
    for trial in range(n_trials):
        n = rng.randint(4, 9)
        taxa, chars, matrix = random_matrix(rng, n, rng.randint(1, 15))
        og = taxa[0]
        off = ap.irreversible_offset(taxa, chars, matrix, og)
        anc = dict((c, 0) for c in chars)
        mat2 = dict(matrix)
        mat2[ap.ANCESTOR] = anc
        for _ in range(3):
            g = ap.random_tree(taxa, rng)
            own = ap.score_tree(ap.root_at_leaf(g, og), chars, matrix, "irreversible")
            # attach the ancestor on the outgroup's branch, root on it, leave root state free
            sub = ap.root_at_leaf(g, og).children[1]
            y = ap.Node(children=[ap.Node(name=og), sub])
            rho = ap.Node(children=[ap.Node(name=ap.ANCESTOR), y])
            ext = 0
            for c in chars:
                col = {t: mat2[t][c] for t in mat2}
                costs = ap.sankoff_costs(rho, col, ap.cost_irreversible, {})
                ext += min(costs.values())
            assert own == ext + off, (trial, own, ext, off)
    print(f"ok   irreversible TNT/PAUP* emulation = Anvaya score - offset on {n_trials} matrices")


def check_external_log_parsing():
    """Tree lengths are read correctly from a real TNT log (TNT 1.6, Linux):
    the column header '0' of the `length` table must not be taken as a length."""
    path = os.path.join(HERE, "data", "tnt_ienum_log.txt")
    assert ap.parse_best_score(open(path).read()) == 459
    multi = "Tree lengths\n\n             0     1\n\n      0    461   459\n"
    assert ap.parse_best_score(multi) == 459
    assert ap.parse_best_score("Score of best tree(s) found = 468\n") == 468
    assert ap.parse_best_score("no search was run\n") is None
    print("ok   TNT / PAUP* tree lengths parsed correctly from real log output")


if __name__ == "__main__":
    full = "--full" in sys.argv
    check_external_log_parsing()
    check_scorer()
    check_irreversible_offset()
    check_search_exact(full=full)
    print("all engine checks passed")
