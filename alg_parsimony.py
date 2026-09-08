#!/usr/bin/env python3
"""
alg_parsimony.py
=================

Treat shared ancestral-linkage-group (ALG) fusions as phylogenetic characters
and find the most parsimonious tree(s) relating a set of genomes, given a
directory of .rbh files (the format produced by odp-style reciprocal-best-hit
/ Oxford-plot pipelines that compare each genome against the Simakov et al.
2022 BCnS ancestral linkage groups).

--------------------------------------------------------------------------
WHAT THIS SCRIPT DOES
--------------------------------------------------------------------------
1. Parses every *.rbh file in an input directory. Each file is expected to
   contain (at least) these columns:
       gene_group                 -> the ancestral linkage group (ALG) ID
                                      for this ortholog (e.g. "Ea", "D")
       <species>_scaf             -> the scaffold/chromosome this ortholog
                                      sits on in the query genome
       <species>_pos              -> its coordinate on that scaffold
       whole_FET                  -> (optional) Fisher's-exact-test p-value
                                      for enrichment of that ALG on that
                                      scaffold (already present in odp-style
                                      output; used as an extra significance
                                      filter if present)
   The species name is auto-detected from the "<species>_gene" column in
   each file's header (falls back to the filename stem if that fails), so
   files do not need to be named to match column prefixes exactly.

2. For each species, calls which ALGs are significantly associated with
   which scaffold (>= --min-genes shared orthologs, and whole_FET <= --alpha
   if that column is present). A scaffold carrying >=2 significantly-
   associated ALGs is a candidate fusion.

3. For every such ALG pair sharing a scaffold, classifies the fusion as
   "clean" (the two ALGs occupy distinct blocks along the scaffold) or
   "mixed" (their genes are interleaved along the scaffold), using the
   fraction of adjacent-gene identity switches when genes are ordered by
   position. This gives each (ALG_i, ALG_j, species) triple a state:
       0 = unfused / independent
       1 = fused, clean (still block-organized)
       2 = fused, mixed (interleaved -> extensively rearranged post-fusion)
   Missing data ("?") is used when a species doesn't have enough orthologs
   assigned to one of the two ALGs to make a call at all.

4. Assembles a taxa x characters matrix (one column per ALG pair that is
   fused in at least one species) and writes it out as a TSV.

5. Searches tree space for the most parsimonious tree twice:
     (a) "reversible"     - ordinary (Fitch/Sankoff) parsimony, state changes
                             are free to go either direction (fusions could,
                             in principle, revert).
     (b) "irreversible"   - Dollo/Camin-Sokal-style parsimony: characters are
                             ordered 0 < 1 < 2 (unfused < clean fusion <
                             fused+mixed) and a lineage is NOT allowed to
                             decrease state (no un-fusing, no un-mixing).
                             This requires designating a root/outgroup (see
                             --outgroup); the ancestral state for every
                             character is fixed at 0.
   For small numbers of taxa (<= --max-exhaustive, default 8) every distinct
   unrooted binary tree topology is scored exactly. Above that, a heuristic
   search (random-addition starting tree + NNI hill-climbing, multiple
   restarts) is used instead - for serious/large analyses consider a
   dedicated phylogenetics package (e.g. TNT, PAUP*, phangorn) instead.

6. For the best tree in each run, reconstructs ancestral states at internal
   nodes and reports, for every clade, which fusion/mixing characters
   originate on the branch leading to that clade - i.e. its macrosyntenic
   synapomorphies.

--------------------------------------------------------------------------
OPTIONAL: BACKBONE / CONSTRAINT TREES
--------------------------------------------------------------------------
Pass --constraint with a Newick string (or a path to a file containing one)
to restrict the search to trees compatible with a partially- or fully-
resolved backbone. This is useful as a positive control (e.g. "force all
centipedes to be monophyletic") while leaving relationships among other,
less certain, taxa as an unresolved polytomy for the ALG characters to
resolve. Rules:
  * The constraint tree may contain some, all, or none of the taxa found
    in --indir. Taxa it omits are left completely free to attach anywhere.
  * Any node with 2+ children asserts that its descendant leaf set is
    monophyletic. A polytomy (3+ children at a node) asserts only that the
    combined descendant set is monophyletic - relationships *within* it are
    left for the character data to resolve.
  * Example: "((Cent1,Cent2,Cent3),Outgrp,Other1,Other2);" forces
    {Cent1,Cent2,Cent3} to be monophyletic (a positive control) while
    leaving Outgrp/Other1/Other2/the centipede clade in an unresolved
    polytomy for the ALG characters to sort out.
  * The same constraint is applied to both the reversible and the
    irreversible parsimony runs.

--------------------------------------------------------------------------
OPTIONAL: TNT / PAUP* EXPORT (AND OPTIONAL EXECUTION)
--------------------------------------------------------------------------
This script's own search is a lightweight, dependency-free implementation
meant for small-to-moderate numbers of taxa. For larger or more demanding
analyses you may prefer a dedicated, validated parsimony program:
  --export-nexus   writes character_matrix.nex, a NEXUS file readable by
                    PAUP* (and many other phylogenetics programs), including
                    a PAUP command block that runs heuristic search under
                    both an unordered (Fitch) and an irrev.up (Camin-Sokal,
                    i.e. irreversible) character-type model, honoring
                    --constraint if given.
  --export-tnt     writes character_matrix.tnt, a Hennig86/TNT-format data
                    file plus a minimal TNT command script for a standard
                    (reversible, Fitch) heuristic search. TNT has no simple
                    built-in irreversible/Dollo character type and the
                    generated script does not attempt to enforce
                    --constraint, so treat the TNT export as a cross-check
                    of the *reversible* search only.
  --run-paup / --run-tnt   after exporting, also try to invoke the `paup`
                    / `tnt` executable (override the binary name/path with
                    --paup-bin / --tnt-bin) on the generated script. This is
                    best-effort: if the binary isn't found on PATH nothing
                    breaks, a warning is printed, and the exported files are
                    left for you to run manually. PAUP* is commercial
                    software you must install/license separately; TNT is
                    free for academic use from its own distributor. Neither
                    is required for the rest of this script to work.

--------------------------------------------------------------------------
IMPORTANT ASSUMPTIONS / SIMPLIFICATIONS (read before trusting the output)
--------------------------------------------------------------------------
* Fusion/mixing calls are made per scaffold using whatever ALG-to-scaffold
  association the whole-scaffold FET already encodes; sub-scaffold
  "breakchrom" region calls (fission/translocation) are NOT modeled here.
* A scaffold carrying 3+ significant ALGs is decomposed into all pairwise
  fusion characters among those ALGs. Biologically this is one compound
  event; treating it pairwise is the common convention in the macrosynteny
  literature but inflates the apparent number of independent characters
  somewhat - keep this in mind when reading the character count.
* "Irreversible" parsimony needs a real root. If you don't pass --outgroup,
  the taxon with the fewest derived (state>0) character calls is picked
  automatically and reported - sanity-check that choice.
* This is a from-scratch, dependency-light parsimony implementation
  (pandas is the only real requirement) for tens of taxa at most. It is not
  a substitute for validated phylogenetics software on large datasets.

--------------------------------------------------------------------------
USAGE
--------------------------------------------------------------------------
    python3 alg_parsimony.py --indir /path/to/rbh_files --outdir results/

    Useful options:
      --alpha 0.05            whole_FET significance cutoff
      --min-genes 3           min shared orthologs to call an ALG present
      --mixing-threshold 0.3  fraction of order-switches -> "mixed" call
      --outgroup SpeciesX      root taxon for the irreversible analysis
      --max-exhaustive 8       taxon-count ceiling for exact tree search
      --nni-restarts 10        random restarts for the heuristic search
      --seed 1                 RNG seed

Output written to --outdir:
      character_matrix.tsv
      alg_scaffold_calls.tsv
      tree_reversible.nwk / .txt
      tree_irreversible.nwk / .txt
      synapomorphies_reversible.txt
      synapomorphies_irreversible.txt
      report.txt   (human-readable summary of everything above)
"""

import argparse
import glob
import itertools
import os
import random
import re
import sys
import zlib
from collections import defaultdict, Counter

import pandas as pd

BIG = 10 ** 6  # effectively "infinite" cost for a forbidden transition
STATES = (0, 1, 2)
STATE_NAME = {0: "unfused", 1: "fused-clean", 2: "fused-mixed", "?": "?"}


# ==========================================================================
# 1. Parsing
# ==========================================================================

def find_gene_col_candidates(header_cols):
    """All '<X>_gene' column prefixes present in a header."""
    out = []
    for c in header_cols:
        m = re.match(r"^(.+)_gene$", c)
        if m:
            out.append(m.group(1))
    return out


def detect_reference_prefix(all_file_candidates, override=None):
    """Across every .rbh file in the directory, each header has exactly two
    '<X>_gene'-style columns: one is the shared ancestral/reference linkage
    group set (same name in every file, e.g. 'BCnSSimakov2022' or
    'odp_LGs_from_Avul_Ixod_Nlon_Varr'), the other is that file's own query
    species (a different name per file). The reference is whichever prefix
    recurs, unchanged, in (almost) every file; returns that name, or
    `override` if explicitly given."""
    if override:
        return override
    counts = defaultdict(int)
    for cands in all_file_candidates:
        for c in set(cands):
            counts[c] += 1
    if not counts:
        return None
    n_files = len(all_file_candidates)
    # the reference prefix should appear in (almost) every file; a species
    # prefix should appear in only one (or a few, if a species is split
    # across multiple files)
    ref_candidates = [c for c, n in counts.items() if n >= max(2, n_files - 1)]
    if len(ref_candidates) == 1:
        return ref_candidates[0]
    if len(ref_candidates) > 1:
        # ambiguous - most-common one wins, but warn
        best = max(ref_candidates, key=lambda c: counts[c])
        sys.stderr.write(
            f"[warn] multiple plausible shared reference-column prefixes found "
            f"{ref_candidates}; guessing '{best}' (use --reference-prefix to override)\n"
        )
        return best
    return None  # couldn't tell - each file will fall back individually


def detect_species_prefix(header_cols, filepath, reference_prefix=None):
    """Find the species-specific column prefix: the '<X>_gene' column that
    is NOT the shared reference-linkage-group column."""
    candidates = find_gene_col_candidates(header_cols)
    if reference_prefix is not None:
        remaining = [c for c in candidates if c != reference_prefix]
        if len(remaining) == 1:
            return remaining[0]
    if len(candidates) == 1:
        return candidates[0]
    # fall back to filename stem
    stem = os.path.splitext(os.path.basename(filepath))[0]
    sys.stderr.write(
        f"[warn] could not uniquely detect species prefix in {filepath} "
        f"(candidates={candidates}); falling back to filename stem '{stem}'\n"
    )
    return stem


def parse_rbh_directory(indir, reference_prefix_override=None):
    """Returns (species_data, reference_prefix):
      species_data: dict species_name -> DataFrame[alg, scaf, pos, fet, ref_pos, color]
      reference_prefix: the detected/overridden shared ALG-reference column prefix (or None)
    ref_pos (the reference/ALG's own internal position, when available) and color (the
    .rbh file's own per-ALG color, when available) are only used for dot plots - all
    fusion/mixing/parsimony logic ignores them."""
    files = sorted(glob.glob(os.path.join(indir, "*.rbh")))
    if not files:
        sys.exit(f"No .rbh files found in {indir}")

    # pass 1: figure out the shared reference-column prefix across the whole
    # directory (e.g. 'BCnSSimakov2022' or 'odp_LGs_from_<taxa>...') so it's
    # never mistaken for a species
    headers = {}
    for fp in files:
        with open(fp) as f:
            headers[fp] = f.readline().rstrip("\n").split("\t")
    all_candidates = [find_gene_col_candidates(headers[fp]) for fp in files]
    reference_prefix = detect_reference_prefix(all_candidates, override=reference_prefix_override)
    if reference_prefix:
        print(f"[info] detected shared reference/ALG column prefix: '{reference_prefix}'")

    species_data = {}
    for fp in files:
        df = pd.read_csv(fp, sep="\t", dtype=str)
        sp = detect_species_prefix(df.columns, fp, reference_prefix=reference_prefix)
        needed = ["gene_group", f"{sp}_scaf", f"{sp}_pos"]
        missing = [c for c in needed if c not in df.columns]
        if missing:
            sys.exit(f"{fp}: missing expected column(s) {missing}")

        out = pd.DataFrame({
            "alg": df["gene_group"],
            "scaf": df[f"{sp}_scaf"],
            "pos": pd.to_numeric(df[f"{sp}_pos"], errors="coerce"),
        })
        if "whole_FET" in df.columns:
            out["fet"] = pd.to_numeric(df["whole_FET"], errors="coerce")
        else:
            out["fet"] = float("nan")
        ref_pos_col = f"{reference_prefix}_pos" if reference_prefix else None
        if ref_pos_col and ref_pos_col in df.columns:
            out["ref_pos"] = pd.to_numeric(df[ref_pos_col], errors="coerce")
        else:
            out["ref_pos"] = float("nan")
        out["color"] = df["color"] if "color" in df.columns else None
        # orthology_type ('direct' vs 'wgd_inferred') is written by genomes_to_rbh.py
        # --detect-wgd; carried through so plots can highlight rescued ohnologs and so
        # the ohnolog summary table can be built. Absent for .rbh files from other
        # sources, in which case everything is treated as 'direct'.
        if "orthology_type" in df.columns:
            out["orthology_type"] = df["orthology_type"].fillna("direct")
        else:
            out["orthology_type"] = "direct"
        out = out.dropna(subset=["alg", "scaf"])
        out = out[out["alg"].str.lower() != "none"]  # unassigned genes (e.g. taxon-vs-taxon rbh files)

        n_inferred = int((out["orthology_type"] == "wgd_inferred").sum())
        if n_inferred:
            print(f"[info] {sp}: {n_inferred}/{len(out)} orthologs are WGD-inferred "
                  f"(rescued ohnolog hits from genomes_to_rbh.py --detect-wgd, not direct RBH)")

        if sp in species_data:
            sys.stderr.write(
                f"[warn] species '{sp}' already seen (from another file); "
                f"appending rows from {fp}\n"
            )
            species_data[sp] = pd.concat([species_data[sp], out], ignore_index=True)
        else:
            species_data[sp] = out

    return species_data, reference_prefix


# ==========================================================================
# 2+3. ALG-scaffold calls, fusion + mixing classification
# ==========================================================================

def call_alg_scaffold_associations(df, alpha, min_genes):
    """Returns dict scaf -> list of algs significantly present on it."""
    calls = defaultdict(list)
    grouped = df.groupby(["scaf", "alg"])
    for (scaf, alg), sub in grouped:
        n = len(sub)
        if n < min_genes:
            continue
        fet = sub["fet"].min()  # if it varies within a group, be conservative
        if pd.notna(fet) and fet > alpha:
            continue
        calls[scaf].append(alg)
    return calls


def mixing_state(df, scaf, alg_x, alg_y, mixing_threshold):
    """0 never used here (caller already knows they're fused); returns 1 or 2."""
    sub = df[(df["scaf"] == scaf) & (df["alg"].isin([alg_x, alg_y]))].dropna(subset=["pos"])
    sub = sub.sort_values("pos")
    labels = list(sub["alg"])
    if len(labels) < 2:
        return 1  # not enough info to call mixing; default to "clean"
    transitions = sum(1 for i in range(len(labels) - 1) if labels[i] != labels[i + 1])
    frac = transitions / (len(labels) - 1)
    return 2 if frac > mixing_threshold else 1


def build_fusion_calls(species_data, alpha, min_genes, mixing_threshold):
    """
    Returns:
      fusion_state[species][(algA,algB)] -> 1 or 2   (only entries that ARE fused)
      alg_total_counts[species][alg]     -> int       (for missing-data calls)
      scaffold_report: list of (species, scaf, algs, state-per-pair) for the TSV
    """
    fusion_state = defaultdict(dict)
    alg_total_counts = {}
    scaffold_rows = []

    for sp, df in species_data.items():
        alg_total_counts[sp] = df.groupby("alg").size().to_dict()
        calls = call_alg_scaffold_associations(df, alpha, min_genes)
        for scaf, algs in calls.items():
            algs = sorted(set(algs))
            if len(algs) < 2:
                continue
            for a, b in itertools.combinations(algs, 2):
                state = mixing_state(df, scaf, a, b, mixing_threshold)
                key = tuple(sorted((a, b)))
                if key in fusion_state[sp]:
                    state = max(state, fusion_state[sp][key])
                fusion_state[sp][key] = state
                scaffold_rows.append((sp, scaf, a, b, state))

    return fusion_state, alg_total_counts, scaffold_rows


def build_character_matrix(species_list, fusion_state, alg_total_counts, min_genes):
    """Returns (characters: list of (algA,algB), matrix: dict[sp][char] -> 0/1/2/'?')"""
    all_pairs = set()
    for sp in species_list:
        all_pairs |= set(fusion_state[sp].keys())
    characters = sorted(all_pairs)

    matrix = {sp: {} for sp in species_list}
    for sp in species_list:
        counts = alg_total_counts.get(sp, {})
        for pair in characters:
            a, b = pair
            if counts.get(a, 0) < min_genes or counts.get(b, 0) < min_genes:
                matrix[sp][pair] = "?"
            else:
                matrix[sp][pair] = fusion_state[sp].get(pair, 0)
    return characters, matrix


# ==========================================================================
# 4. Unrooted binary tree enumeration (sequential taxon addition)
# ==========================================================================

class Graph:
    """Undirected graph representation of an unrooted binary tree.
    Leaf nodes are strings (taxon names); internal nodes are negative ints."""

    def __init__(self):
        self.adj = defaultdict(set)
        self._next_internal = -1

    def copy(self):
        g = Graph()
        g.adj = defaultdict(set, {k: set(v) for k, v in self.adj.items()})
        g._next_internal = self._next_internal
        return g

    def fresh_internal(self):
        nid = self._next_internal
        self._next_internal -= 1
        return nid

    def add_edge(self, a, b):
        self.adj[a].add(b)
        self.adj[b].add(a)

    def remove_edge(self, a, b):
        self.adj[a].discard(b)
        self.adj[b].discard(a)

    def edges(self):
        seen = set()
        out = []
        for a, nbrs in self.adj.items():
            for b in nbrs:
                if (b, a) not in seen:
                    seen.add((a, b))
                    out.append((a, b))
        return out


def initial_tree(t0, t1, t2):
    g = Graph()
    center = g.fresh_internal()
    g.add_edge(center, t0)
    g.add_edge(center, t1)
    g.add_edge(center, t2)
    return g


def insert_taxon_on_edge(g, taxon, edge):
    a, b = edge
    g2 = g.copy()
    g2.remove_edge(a, b)
    new_internal = g2.fresh_internal()
    g2.add_edge(a, new_internal)
    g2.add_edge(new_internal, b)
    g2.add_edge(new_internal, taxon)
    return g2


def enumerate_unrooted_trees(taxa):
    """Exact enumeration of all distinct unrooted binary tree topologies.
    Number of trees = (2n-5)!! for n taxa. Use only for small n."""
    if len(taxa) < 3:
        raise ValueError("need >=3 taxa")
    trees = [initial_tree(taxa[0], taxa[1], taxa[2])]
    for t in taxa[3:]:
        new_trees = []
        for g in trees:
            for e in g.edges():
                new_trees.append(insert_taxon_on_edge(g, t, e))
        trees = new_trees
    return trees


def random_tree(taxa, rng):
    order = taxa[:]
    rng.shuffle(order)
    g = initial_tree(order[0], order[1], order[2])
    for t in order[3:]:
        e = rng.choice(g.edges())
        g = insert_taxon_on_edge(g, t, e)
    return g


def suppress_degree2(g):
    """Remove any internal node of degree 2, splicing its two edges into one
    direct edge. Needed after building a graph from a rooted Newick tree
    (whose root, once unrooted, is a degree-2 artifact)."""
    changed = True
    while changed:
        changed = False
        for node in list(g.adj.keys()):
            if isinstance(node, int) and len(g.adj[node]) == 2:
                a, b = tuple(g.adj[node])
                g.remove_edge(node, a)
                g.remove_edge(node, b)
                g.add_edge(a, b)
                del g.adj[node]
                changed = True
                break
    return g


def internal_edges(g):
    """Edges connecting two internal (non-leaf) nodes -> candidates for NNI."""
    return [(a, b) for (a, b) in g.edges() if isinstance(a, int) and isinstance(b, int)]


# ---- backbone / constraint trees -----------------------------------------

def parse_newick(text):
    """Minimal Newick parser -> a Node tree (see Node class further below).
    Supports (nested) polytomies, ignores branch lengths and node labels."""
    s = text.strip()
    if s.endswith(";"):
        s = s[:-1]
    pos = [0]

    def parse_label():
        start = pos[0]
        while pos[0] < len(s) and s[pos[0]] not in ",():;":
            pos[0] += 1
        label = s[start:pos[0]].strip()
        return label or None

    def parse_node():
        if s[pos[0]] == "(":
            pos[0] += 1
            children = [parse_node()]
            while s[pos[0]] == ",":
                pos[0] += 1
                children.append(parse_node())
            if s[pos[0]] != ")":
                raise ValueError(f"malformed newick near position {pos[0]}: {s[max(0,pos[0]-20):pos[0]+20]}")
            pos[0] += 1
            parse_label()  # discard internal label / branch length, if any
            return Node(children=children)
        else:
            label = parse_label()
            # a leaf label may be followed by ':<branch length>' already
            # consumed by parse_label's stop-set; nothing more to do here
            if pos[0] < len(s) and s[pos[0]] == ":":
                # consume branch length explicitly if parse_label stopped at ':'
                pass
            return Node(name=label)

    root = parse_node()
    return root


def required_clades_from_constraint(root):
    """Every internal node (root included) asserts that its full set of
    descendant leaves is monophyletic. Returns a set of frozensets."""
    clades = set()

    def walk(node):
        if not node.is_leaf():
            leafset = frozenset(node.leaves())
            if len(leafset) >= 2:
                clades.add(leafset)
            for c in node.children:
                walk(c)

    walk(root)
    return clades


def satisfies_constraints(g, required_clades, constraint_taxa):
    """True if every required clade is present as a bipartition of g, once
    restricted to the taxa actually named in the constraint tree (taxa the
    constraint didn't mention are allowed to sit anywhere)."""
    if not required_clades:
        return True
    all_bip = bipartitions(g)
    for clade in required_clades:
        found = False
        for bip in all_bip:
            s1, s2 = tuple(bip)
            if (s1 & constraint_taxa) == clade or (s2 & constraint_taxa) == clade:
                found = True
                break
        if not found:
            return False
    return True


def resolve_polytomies_randomly(node, rng):
    """Return a new, strictly-binary Node tree consistent with `node`
    (every polytomy is broken up arbitrarily/randomly)."""
    if node.is_leaf():
        return Node(name=node.name)
    pool = [resolve_polytomies_randomly(c, rng) for c in node.children]
    while len(pool) > 2:
        rng.shuffle(pool)
        a = pool.pop()
        b = pool.pop()
        pool.append(Node(children=[a, b]))
    return Node(children=pool)


def node_tree_to_graph(node):
    """Convert a strictly-binary rooted Node tree into our unrooted Graph
    representation (root is unrooted away via suppress_degree2)."""
    g = Graph()

    def build(n):
        if n.is_leaf():
            return n.name
        cids = [build(c) for c in n.children]
        if len(cids) == 1:
            return cids[0]
        nid = g.fresh_internal()
        for cid in cids:
            g.add_edge(nid, cid)
        return nid

    build(node)
    suppress_degree2(g)
    return g


def constrained_random_tree(taxa, constraint_root, rng):
    """A random fully-resolved unrooted binary tree over `taxa` that is
    compatible with `constraint_root` (a possibly-partial, possibly-
    polytomous Newick tree parsed with parse_newick)."""
    resolved = resolve_polytomies_randomly(constraint_root, rng)
    g = node_tree_to_graph(resolved)
    constrained_taxa = set(resolved.leaves())
    remaining = [t for t in taxa if t not in constrained_taxa]
    rng.shuffle(remaining)

    if not g.edges():
        # 0 or 1 constrained taxa - bootstrap an ordinary 3-taxon start
        pool = list(constrained_taxa)
        while len(pool) < 3 and remaining:
            pool.append(remaining.pop())
        if len(pool) < 3:
            raise ValueError("Not enough taxa to build a tree")
        g = initial_tree(pool[0], pool[1], pool[2])

    for t in remaining:
        e = rng.choice(g.edges())
        g = insert_taxon_on_edge(g, t, e)
    return g


def nni_neighbors(g):
    """All trees one NNI move away from g."""
    out = []
    for (u, v) in internal_edges(g):
        u_others = [n for n in g.adj[u] if n != v]
        v_others = [n for n in g.adj[v] if n != u]
        if len(u_others) != 2 or len(v_others) != 2:
            continue
        a, b = u_others
        c, d = v_others
        # move 1: swap b and c
        g1 = g.copy()
        g1.remove_edge(u, b); g1.remove_edge(v, c)
        g1.add_edge(u, c); g1.add_edge(v, b)
        out.append(g1)
        # move 2: swap b and d
        g2 = g.copy()
        g2.remove_edge(u, b); g2.remove_edge(v, d)
        g2.add_edge(u, d); g2.add_edge(v, b)
        out.append(g2)
    return out


# ---- root an unrooted graph at a given leaf (outgroup rooting) ----

class Node:
    __slots__ = ("name", "children")

    def __init__(self, name=None, children=None):
        self.name = name  # set only for leaves
        self.children = children or []

    def is_leaf(self):
        return not self.children

    def leaves(self):
        if self.is_leaf():
            return [self.name]
        out = []
        for c in self.children:
            out.extend(c.leaves())
        return out


def root_at_leaf(g, root_leaf):
    """Root the unrooted graph so that `root_leaf` is an outgroup:
    Root -> [root_leaf, rest-of-tree]."""

    def build(node_id, came_from):
        children_ids = [n for n in g.adj[node_id] if n != came_from]
        if isinstance(node_id, str):
            return Node(name=node_id)
        node = Node()
        for cid in children_ids:
            node.children.append(build(cid, node_id))
        return node

    root_neighbor = next(iter(g.adj[root_leaf]))
    subtree = build(root_neighbor, root_leaf)
    root = Node(children=[Node(name=root_leaf), subtree])
    return root


def root_arbitrarily(g, at_leaf):
    """For the reversible analysis the score doesn't depend on
    rooting, but Sankoff still needs *some* rooted tree to recurse over."""
    return root_at_leaf(g, at_leaf)


# ==========================================================================
# 5. Sankoff parsimony
# ==========================================================================

def cost_reversible(a, b):
    return 0 if a == b else 1


def cost_irreversible(parent_state, child_state):
    if child_state == parent_state:
        return 0
    if child_state > parent_state:
        return 1
    return BIG


def sankoff_costs(node, char_matrix_col, cost_fn, memo):
    """Bottom-up Sankoff DP. Returns dict state->min_cost for `node`,
    given the (species -> state) mapping for one character. Caches in memo
    keyed by id(node)."""
    key = id(node)
    if key in memo:
        return memo[key]
    if node.is_leaf():
        obs = char_matrix_col.get(node.name, "?")
        if obs == "?":
            res = {s: 0 for s in STATES}
        else:
            res = {s: (0 if s == obs else BIG) for s in STATES}
    else:
        child_costs = [sankoff_costs(c, char_matrix_col, cost_fn, memo) for c in node.children]
        res = {}
        for s in STATES:
            total = 0
            for cc in child_costs:
                total += min(cost_fn(s, s2) + cc[s2] for s2 in STATES)
            res[s] = total
    memo[key] = res
    return res


def score_tree(root, characters, matrix, mode):
    """mode: 'reversible' or 'irreversible'. Returns total parsimony score."""
    cost_fn = cost_reversible if mode == "reversible" else cost_irreversible
    total = 0
    for char in characters:
        col = {sp: matrix[sp][char] for sp in matrix}
        memo = {}
        costs = sankoff_costs(root, col, cost_fn, memo)
        if mode == "reversible":
            total += min(costs.values())
        else:
            total += costs[0]  # ancestral state fixed at 0 (unfused)
        if total >= BIG:
            return BIG
    return total


def traceback_states(root, char_matrix_col, cost_fn, mode):
    """Full Sankoff traceback: returns dict id(node) -> assigned state,
    plus the node objects themselves (via a side list) so callers can map
    node identity -> chosen state."""
    memo = {}
    sankoff_costs(root, char_matrix_col, cost_fn, memo)
    assigned = {}

    def assign(node, parent_state):
        costs = memo[id(node)]
        if parent_state is None:
            if mode == "irreversible":
                s = 0
            else:
                s = min(costs, key=lambda k: costs[k])
        else:
            s = min(STATES, key=lambda s2: cost_fn(parent_state, s2) + costs[s2])
        assigned[id(node)] = s
        for c in node.children:
            assign(c, s)

    assign(root, None)
    return assigned


# ==========================================================================
# Homoplasy diagnostics: Consistency Index (CI) / Retention Index (RI)
# ==========================================================================

def compute_ci_ri(root, characters, matrix, mode):
    """Per-character and ensemble CI/RI for the given tree under the given
    criterion. CI = m/s (m = minimum possible steps for that character
    alone, s = its actual steps on this tree - a low CI flags homoplasy,
    i.e. the same state arising more than once). RI = (g-s)/(g-m) (g = the
    steps that character would need on the least-resolved/star tree,
    i.e. taxon count minus the largest single-state count - the standard
    way to bound the "worst case" for RI). Missing ('?') taxa are excluded
    from a character's own m/s/g. Returns (per_char: dict char -> dict with
    keys m,s,g,ci,ri; ensemble: dict with keys ci,ri)."""
    cost_fn = cost_reversible if mode == "reversible" else cost_irreversible
    per_char = {}
    sum_m = sum_s = 0
    sum_g_minus_s = sum_g_minus_m = 0
    for char in characters:
        col = {sp: matrix[sp][char] for sp in matrix}
        non_missing = [v for v in col.values() if v != "?"]
        n_taxa = len(non_missing)
        m = max(len(set(non_missing)) - 1, 0)
        memo = {}
        costs = sankoff_costs(root, col, cost_fn, memo)
        s = min(costs.values()) if mode == "reversible" else costs[0]
        g = (n_taxa - max(Counter(non_missing).values())) if n_taxa else 0
        ci = (m / s) if s > 0 else 1.0
        ri = ((g - s) / (g - m)) if g != m else 1.0
        per_char[char] = {"m": m, "s": s, "g": g, "ci": ci, "ri": ri}
        sum_m += m
        sum_s += s
        sum_g_minus_s += (g - s)
        sum_g_minus_m += (g - m)
    ensemble_ci = (sum_m / sum_s) if sum_s > 0 else 1.0
    ensemble_ri = (sum_g_minus_s / sum_g_minus_m) if sum_g_minus_m > 0 else 1.0
    return per_char, {"ci": ensemble_ci, "ri": ensemble_ri}


# ==========================================================================
# Bootstrap support and permutation/signal test
# ==========================================================================

def _search(species_list, characters, matrix, mode, root_leaf, exhaustive, restarts, rng, constraint):
    """Shared entry point used by both the main search and every bootstrap/
    permutation replicate, so replicates use exactly the same search
    strategy (exhaustive vs. heuristic) as the reported analysis."""
    if exhaustive:
        return best_tree_exhaustive(species_list, characters, matrix, mode, root_leaf, constraint=constraint)
    return best_tree_heuristic(species_list, characters, matrix, mode, root_leaf, restarts, rng,
                                constraint=constraint)


def bootstrap_support(species_list, characters, matrix, mode, root_leaf, target_clades, n_reps,
                       exhaustive, restarts, rng, constraint):
    """Nonparametric bootstrap: resample characters with replacement (same
    total count) n_reps times, rerun the search each time, and for each
    clade in `target_clades` (frozensets of taxon names, from the reported
    tree) tally the fraction of replicates whose best tree contains that
    exact clade as one side of a bipartition. Returns dict frozenset ->
    fraction in [0, 1]."""
    counts = {c: 0 for c in target_clades}
    n_chars = len(characters)
    for _ in range(n_reps):
        resampled = [characters[rng.randrange(n_chars)] for _ in range(n_chars)]
        _, trees = _search(species_list, resampled, matrix, mode, root_leaf, exhaustive, restarts,
                            rng, constraint)
        rep_bips = bipartitions(trees[0])
        for clade in target_clades:
            if any(clade in bip for bip in rep_bips):
                counts[clade] += 1
    return {c: counts[c] / n_reps for c in target_clades}


def permutation_test(species_list, characters, matrix, mode, root_leaf, n_reps, exhaustive, restarts,
                      rng, constraint, observed_score):
    """Signal/PTP-style test: for each of n_reps replicates, independently
    shuffle which taxon has which state within every character (preserving
    each character's own state multiset and missing-data pattern, but
    destroying any real taxon-to-taxon correlation among characters), then
    find the best achievable score on the shuffled matrix. Returns
    (permuted_scores: list[int], p_value), where p_value is the fraction of
    replicates at least as parsimonious as the real, observed score (with a
    +1/+1 continuity correction) - a small p-value means the real data are
    more structured than shuffled data of the same composition, i.e. there
    is real phylogenetic signal."""
    permuted_scores = []
    for _ in range(n_reps):
        perm_matrix = {sp: {} for sp in species_list}
        for char in characters:
            values = [matrix[sp][char] for sp in species_list]
            rng.shuffle(values)
            for sp, v in zip(species_list, values):
                perm_matrix[sp][char] = v
        score, _ = _search(species_list, characters, perm_matrix, mode, root_leaf, exhaustive, restarts,
                            rng, constraint)
        permuted_scores.append(score)
    n_at_least_as_good = sum(1 for s in permuted_scores if s <= observed_score)
    p_value = (n_at_least_as_good + 1) / (n_reps + 1)
    return permuted_scores, p_value


# ==========================================================================
# Search drivers
# ==========================================================================

def best_tree_exhaustive(taxa, characters, matrix, mode, root_leaf, constraint=None):
    trees = enumerate_unrooted_trees(taxa)
    best_score = None
    best_trees = []
    skipped = 0
    for g in trees:
        if constraint and not satisfies_constraints(g, constraint["required_clades"], constraint["taxa"]):
            skipped += 1
            continue
        root = root_at_leaf(g, root_leaf)
        s = score_tree(root, characters, matrix, mode)
        if best_score is None or s < best_score:
            best_score = s
            best_trees = [g]
        elif s == best_score:
            best_trees.append(g)
    if constraint and best_score is None:
        raise RuntimeError(
            "No tree topology satisfies the given --constraint (all "
            f"{len(trees)} candidate trees were incompatible). Check the "
            "constraint Newick for a contradiction with itself or with the "
            "taxon names in your .rbh files."
        )
    return best_score, best_trees


def best_tree_heuristic(taxa, characters, matrix, mode, root_leaf, restarts, rng, constraint=None):
    best_score = None
    best_trees = []
    for _ in range(restarts):
        if constraint:
            g = constrained_random_tree(taxa, constraint["root_node"], rng)
        else:
            g = random_tree(taxa, rng)
        cur_score = score_tree(root_at_leaf(g, root_leaf), characters, matrix, mode)
        improved = True
        while improved:
            improved = False
            candidates = nni_neighbors(g)
            if constraint:
                candidates = [ng for ng in candidates
                              if satisfies_constraints(ng, constraint["required_clades"], constraint["taxa"])]
            for ng in candidates:
                s = score_tree(root_at_leaf(ng, root_leaf), characters, matrix, mode)
                if s < cur_score:
                    cur_score, g, improved = s, ng, True
        if best_score is None or cur_score < best_score:
            best_score, best_trees = cur_score, [g]
        elif cur_score == best_score:
            if not any(same_topology(g, bt) for bt in best_trees):
                best_trees.append(g)
    return best_score, best_trees


def same_topology(g1, g2):
    # cheap check: same set of bipartitions (leaf-set splits induced by internal edges)
    return bipartitions(g1) == bipartitions(g2)


def bipartitions(g):
    leaves = set(n for n in g.adj if isinstance(n, str))
    out = set()
    for (a, b) in internal_edges(g):
        side = set()
        stack = [a]
        seen = {a, b}
        # actually need component of `a` when edge (a,b) removed; do a fresh BFS
        seen = {b}
        stack = [a]
        comp = set()
        while stack:
            n = stack.pop()
            if n in seen:
                continue
            seen.add(n)
            if isinstance(n, str):
                comp.add(n)
            for nb in g.adj[n]:
                if nb not in seen:
                    stack.append(nb)
        s1 = frozenset(comp)
        s2 = frozenset(leaves - comp)
        out.add(frozenset([s1, s2]))
    return out


# ==========================================================================
# Rendering
# ==========================================================================

def to_newick(node, support=None):
    """support: optional dict frozenset(leaf names) -> fraction in [0,1].
    When given, every internal node whose leaf set is a key in `support`
    gets a standard Newick internal-node label (e.g. '95') immediately
    after its closing parenthesis, interpreted by tree viewers as that
    clade's support value. No branch lengths are written anywhere - this
    is a discrete-character parsimony tree with no molecular clock, and
    fabricating lengths would misrepresent that."""
    if node.is_leaf():
        return node.name
    inner = "(" + ",".join(to_newick(c, support) for c in node.children) + ")"
    if support is not None:
        leafset = frozenset(node.leaves())
        if leafset in support:
            inner += f"{support[leafset] * 100:.0f}"
    return inner


def to_ascii(node, prefix="", is_last=True, out=None):
    if out is None:
        out = []
    label = node.name if node.is_leaf() else "+"
    connector = "`- " if is_last else "|- "
    out.append(prefix + (connector if prefix else "") + label)
    new_prefix = prefix + ("   " if is_last else "|  ")
    for i, c in enumerate(node.children):
        to_ascii(c, new_prefix, i == len(node.children) - 1, out)
    return out


# ==========================================================================
# Plotting (optional; requires matplotlib)
# ==========================================================================

def _get_matplotlib():
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        return plt
    except ImportError:
        print("[warn] matplotlib is not installed - skipping plot (pip install matplotlib)")
        return None


def get_global_alg_colors(species_data):
    """One consistent hex color per ALG, shared across every plot AND
    stable across separate analyses/invocations. Prefers a color actually
    given in the .rbh files' own 'color' column (which, when written by
    genomes_to_rbh.py, already uses the same deterministic scheme below);
    falls back to a deterministic hash-based palette pick for ALGs that
    never had one. Uses zlib.crc32 rather than the built-in hash() for
    strings, since the latter is randomized per-process by default and so
    would reassign colors on every run - the exact opposite of what's
    needed to compare figures across separate analyses. Palette and hash
    are identical to genomes_to_rbh.py's color_for(), so a color assigned
    by either script for a given linkage-group name always agrees."""
    palette = [
        "#AB7E26", "#47957F", "#7B7B7B", "#3B69B2", "#764A9B", "#C76E00",
        "#8B0000", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E",
        "#E6AB02", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
    ]
    colors = {}
    for df in species_data.values():
        for alg, color in zip(df["alg"], df["color"]):
            if alg not in colors and isinstance(color, str) and color.strip():
                colors[alg] = color.strip()
    all_algs = sorted(set(a for df in species_data.values() for a in df["alg"]))
    for alg in all_algs:
        if alg not in colors:
            colors[alg] = palette[zlib.crc32(alg.encode("utf-8")) % len(palette)]
    return colors


def build_plot_layout(df, group_col, rank_col, group_order, gap=1):
    """Assign each row of df an integer plot coordinate along one axis,
    grouped by `group_col`. `group_order` is an explicit ordered list of
    group names (e.g. alphabetical ALG order, shared across every species'
    plot so the axis lines up between figures) - groups in the data but not
    in `group_order` are appended, sorted, at the end. `gap` is the number
    of blank units left between consecutive groups.
    Returns (coord: pd.Series aligned to df.index, boundaries: list of
    (group_name, start, end))."""
    present = set(df[group_col])
    order = [g for g in group_order if g in present]
    order += sorted(present - set(order))

    coord = pd.Series(index=df.index, dtype=float)
    boundaries = []
    cursor = 0
    for g in order:
        sub = df[df[group_col] == g]
        if rank_col in sub.columns and sub[rank_col].notna().any():
            sub_sorted = sub.sort_values(rank_col, na_position="last")
        else:
            sub_sorted = sub
        n = len(sub_sorted)
        for offset, idx in enumerate(sub_sorted.index):
            coord[idx] = cursor + offset
        boundaries.append((g, cursor, cursor + n - 1))
        cursor += n + gap
    return coord, boundaries


def plot_dotplot(sp_name, ref_label, df, alg_colors, alg_order, outpath):
    """One Oxford-style synteny dot plot for a single species: x = position
    within its ALG (grouped, ALGs in a fixed order shared across every
    species so plots are comparable), y = position within its own scaffold
    (grouped, scaffolds ordered by how many orthologs they carry)."""
    plt = _get_matplotlib()
    if plt is None:
        return False
    if df.empty:
        return False

    x_coord, x_bounds = build_plot_layout(df, "alg", "ref_pos", alg_order)
    scaf_order = list(df["scaf"].value_counts().index)  # species-specific; count order is fine
    y_coord, y_bounds = build_plot_layout(df, "scaf", "pos", scaf_order)
    colors = [alg_colors.get(a, "#333333") for a in df["alg"]]

    fig_w = max(6.0, min(24.0, len(x_bounds) * 0.35))
    fig_h = max(5.0, min(24.0, len(y_bounds) * 0.30))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.scatter(x_coord, y_coord, c=colors, s=10, linewidths=0)

    for _, s, e in x_bounds:
        ax.axvline(s - 0.5, color="0.85", linewidth=0.5, zorder=0)
    ax.axvline(x_bounds[-1][2] + 0.5, color="0.85", linewidth=0.5, zorder=0)
    ax.set_xticks([(s + e) / 2 for _, s, e in x_bounds])
    ax.set_xticklabels([g for g, _, _ in x_bounds], rotation=90, fontsize=6)

    for _, s, e in y_bounds:
        ax.axhline(s - 0.5, color="0.85", linewidth=0.5, zorder=0)
    ax.axhline(y_bounds[-1][2] + 0.5, color="0.85", linewidth=0.5, zorder=0)
    ax.set_yticks([(s + e) / 2 for _, s, e in y_bounds])
    ax.set_yticklabels([g for g, _, _ in y_bounds], fontsize=6)

    ax.set_xlabel(ref_label)
    ax.set_ylabel(f"{sp_name} scaffold")
    ax.set_title(f"{sp_name} vs {ref_label}")
    present_algs = [a for a in alg_order if a in set(df["alg"])]
    legend_handles = _alg_legend_handles(plt, alg_colors, present_algs)
    ax.legend(handles=legend_handles, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7,
              frameon=False, title="Linkage group", title_fontsize=8)
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return True


def plot_dotplot_wgd(sp_name, ref_label, df, alg_colors, alg_order, outpath):
    """Companion dot plot highlighting WGD-rescued orthologs: points whose
    orthology_type is 'wgd_inferred' (recovered by genomes_to_rbh.py
    --detect-wgd, which plain reciprocal-best-hit would have dropped) are
    drawn at double diameter and ringed in black, so the duplicated
    regions stand out against the directly-recovered background. Returns
    False (writing nothing) if this taxon has no rescued orthologs."""
    plt = _get_matplotlib()
    if plt is None:
        return False
    if df.empty or "orthology_type" not in df.columns:
        return False
    is_wgd = (df["orthology_type"] == "wgd_inferred")
    if not is_wgd.any():
        return False

    x_coord, x_bounds = build_plot_layout(df, "alg", "ref_pos", alg_order)
    scaf_order = list(df["scaf"].value_counts().index)
    y_coord, y_bounds = build_plot_layout(df, "scaf", "pos", scaf_order)
    colors = [alg_colors.get(a, "#333333") for a in df["alg"]]

    fig_w = max(6.0, min(24.0, len(x_bounds) * 0.35))
    fig_h = max(5.0, min(24.0, len(y_bounds) * 0.30))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    direct_mask = ~is_wgd
    ax.scatter(x_coord[direct_mask], y_coord[direct_mask],
               c=[c for c, keep in zip(colors, direct_mask) if keep],
               s=10, linewidths=0, zorder=2)
    ax.scatter(x_coord[is_wgd], y_coord[is_wgd],
               c=[c for c, keep in zip(colors, is_wgd) if keep],
               s=40, linewidths=0.9, edgecolors="black", zorder=3)

    for _, s, e in x_bounds:
        ax.axvline(s - 0.5, color="0.85", linewidth=0.5, zorder=0)
    ax.axvline(x_bounds[-1][2] + 0.5, color="0.85", linewidth=0.5, zorder=0)
    ax.set_xticks([(s + e) / 2 for _, s, e in x_bounds])
    ax.set_xticklabels([g for g, _, _ in x_bounds], rotation=90, fontsize=6)
    for _, s, e in y_bounds:
        ax.axhline(s - 0.5, color="0.85", linewidth=0.5, zorder=0)
    ax.axhline(y_bounds[-1][2] + 0.5, color="0.85", linewidth=0.5, zorder=0)
    ax.set_yticks([(s + e) / 2 for _, s, e in y_bounds])
    ax.set_yticklabels([g for g, _, _ in y_bounds], fontsize=6)

    ax.set_xlabel(ref_label)
    ax.set_ylabel(f"{sp_name} scaffold")
    n_wgd = int(is_wgd.sum())
    ax.set_title(f"{sp_name} vs {ref_label} - WGD-rescued orthologs highlighted "
                 f"({n_wgd}/{len(df)})", fontsize=10)
    marker_handles = [
        plt.matplotlib.lines.Line2D([], [], marker="o", linestyle="none", markersize=4,
                                     markerfacecolor="0.5", markeredgecolor="none",
                                     label="direct RBH"),
        plt.matplotlib.lines.Line2D([], [], marker="o", linestyle="none", markersize=7,
                                     markerfacecolor="0.5", markeredgecolor="black",
                                     markeredgewidth=0.9, label="WGD-rescued ohnolog"),
    ]
    present_algs = [a for a in alg_order if a in set(df["alg"])]
    handles = marker_handles + _alg_legend_handles(plt, alg_colors, present_algs)
    ax.legend(handles=handles, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7,
              frameon=False, title="Linkage group", title_fontsize=8)
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return True


def write_ohnolog_summary(species_data, outpath):
    """Per-taxon / per-scaffold / per-linkage-group counts of WGD-rescued
    versus directly-recovered orthologs. Returns the number of taxa with
    any rescued orthologs (0 = nothing written)."""
    rows = []
    for sp in sorted(species_data):
        df = species_data[sp]
        if "orthology_type" not in df.columns:
            continue
        if not (df["orthology_type"] == "wgd_inferred").any():
            continue
        grouped = df.groupby(["scaf", "alg"])["orthology_type"]
        for (scaf, alg), series in grouped:
            n_wgd = int((series == "wgd_inferred").sum())
            n_direct = int((series == "direct").sum())
            if n_wgd == 0:
                continue
            rows.append((sp, scaf, alg, n_direct, n_wgd, n_direct + n_wgd))
    if not rows:
        return 0
    rows.sort(key=lambda r: (r[0], -r[4], r[1], r[2]))
    with open(outpath, "w") as f:
        f.write("taxon\tscaffold\tlinkage_group\tn_direct\tn_wgd_rescued\tn_total\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")
    return len({r[0] for r in rows})


def ladderize(node):
    """Sort children by descendant leaf count (ascending) at every node, for
    a cleaner staircase layout. Mutates node.children in place. Returns the
    node's own leaf count."""
    if node.is_leaf():
        return 1
    sizes = [ladderize(c) for c in node.children]
    order = sorted(range(len(node.children)), key=lambda i: sizes[i])
    node.children = [node.children[i] for i in order]
    return sum(sizes)


def _tip_depth(node, depth):
    if node.is_leaf():
        depth[id(node)] = 0
        return 0
    d = 1 + max(_tip_depth(c, depth) for c in node.children)
    depth[id(node)] = d
    return d


def assign_tree_xy(root):
    """Rectangular-cladogram layout: leaves aligned at the same rightmost x
    (root has x=0), y = leaf order (ladderized) for leaves / mean of
    children for internal nodes. Returns (xy: dict id(node)->(x,y), max_x)."""
    depth = {}
    _tip_depth(root, depth)
    max_depth = depth[id(root)]
    xy = {}
    counter = [0]

    def walk(node):
        x = max_depth - depth[id(node)]
        if node.is_leaf():
            y = counter[0]
            counter[0] += 1
            xy[id(node)] = (x, y)
            return y
        ys = [walk(c) for c in node.children]
        y = sum(ys) / len(ys)
        xy[id(node)] = (x, y)
        return y

    walk(root)
    return xy, max_depth


def collect_leaf_nodes(node):
    """Like Node.leaves(), but returns the actual Node objects (needed to
    key into the xy layout dict, which is keyed by node identity)."""
    if node.is_leaf():
        return [node]
    out = []
    for c in node.children:
        out.extend(collect_leaf_nodes(c))
    return out


def render_tree_figure(root, characters, matrix, mode, outpath, title):
    """Rectangular cladogram with fusion/mixing synapomorphies labeled on
    the branch where each one originates."""
    plt = _get_matplotlib()
    if plt is None:
        return False

    ladderize(root)
    xy, max_depth = assign_tree_xy(root)
    syn = find_synapomorphies(root, characters, matrix, mode)
    events_by_id = {id(node): events for node, _, events in syn if events}
    leaf_nodes = collect_leaf_nodes(root)
    n_leaves = len(leaf_nodes)

    fig_w = max(6.0, (max_depth + 1) * 1.8)
    fig_h = max(3.0, n_leaves * 0.45)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    def draw(node):
        x, y = xy[id(node)]
        if node.is_leaf():
            return
        child_info = [(c, xy[id(c)]) for c in node.children]
        ys = [cy for _, (cx, cy) in child_info]
        ax.plot([x, x], [min(ys), max(ys)], color="black", linewidth=1)
        for c, (cx, cy) in child_info:
            ax.plot([x, cx], [cy, cy], color="black", linewidth=1)
            events = events_by_id.get(id(c))
            if events:
                label = "\n".join(
                    f"{a}+{b}: {STATE_NAME[ps]}\u2192{STATE_NAME[cs]}" for (a, b), ps, cs in events
                )
                mid_x = (x + cx) / 2
                ax.text(mid_x, cy + 0.10, label, fontsize=6, va="bottom", ha="center",
                        color="darkred")
            draw(c)

    draw(root)
    for leaf in leaf_nodes:
        x, y = xy[id(leaf)]
        ax.text(x + 0.05, y, leaf.name, fontsize=9, va="center", ha="left")

    ax.set_xlim(-0.3, max_depth + 2.2)
    ax.set_ylim(-1, n_leaves)
    ax.axis("off")
    ax.set_title(title, fontsize=10)
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return True


def plot_character_matrix_heatmap(species_list, characters, matrix, outpath):
    """Taxa x characters heatmap: 0=unfused, 1=fused-clean, 2=fused-mixed,
    missing ('?') shown as a distinct hatched grey cell. A quick data-QC
    view of the whole matrix before committing to a tree search."""
    plt = _get_matplotlib()
    if plt is None:
        return False
    if not characters:
        return False
    import numpy as np

    state_to_num = {0: 0, 1: 1, 2: 2, "?": -1}
    grid = np.array([[state_to_num[matrix[sp][c]] for c in characters] for sp in species_list], dtype=float)

    cmap = plt.matplotlib.colors.ListedColormap(["#BEBEBE", "#F2F2F2", "#F4A582", "#B2182B"])
    bounds = [-1.5, -0.5, 0.5, 1.5, 2.5]
    norm = plt.matplotlib.colors.BoundaryNorm(bounds, cmap.N)

    fig_w = max(6.0, min(28.0, len(characters) * 0.45 + 2))
    fig_h = max(4.0, min(24.0, len(species_list) * 0.35 + 1.5))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.imshow(grid, cmap=cmap, norm=norm, aspect="auto")

    ax.set_xticks(range(len(characters)))
    ax.set_xticklabels([f"{a}+{b}" for a, b in characters], rotation=90, fontsize=6)
    ax.set_yticks(range(len(species_list)))
    ax.set_yticklabels(species_list, fontsize=7)
    ax.set_xticks([x - 0.5 for x in range(1, len(characters))], minor=True)
    ax.set_yticks([y - 0.5 for y in range(1, len(species_list))], minor=True)
    ax.grid(which="minor", color="white", linewidth=1)
    ax.tick_params(which="minor", bottom=False, left=False)

    legend_handles = [
        plt.matplotlib.patches.Patch(facecolor="#F2F2F2", edgecolor="0.5", label="0 unfused"),
        plt.matplotlib.patches.Patch(facecolor="#F4A582", edgecolor="0.5", label="1 fused-clean"),
        plt.matplotlib.patches.Patch(facecolor="#B2182B", edgecolor="0.5", label="2 fused-mixed"),
        plt.matplotlib.patches.Patch(facecolor="#BEBEBE", edgecolor="0.5", label="? missing"),
    ]
    ax.legend(handles=legend_handles, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7,
              frameon=False)
    ax.set_title("Character matrix", fontsize=10)
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return True


def _alg_legend_handles(plt, alg_colors, algs_present):
    """Legend patches for a set of ALG names, in a stable, useful order:
    the order they're given in (callers pass an already-sorted or
    priority-ordered sequence), deduplicated."""
    seen = []
    for a in algs_present:
        if a not in seen:
            seen.append(a)
    return [plt.matplotlib.patches.Patch(facecolor=alg_colors.get(a, "#333333"), edgecolor="none", label=a)
            for a in seen]


def plot_ideogram(sp_name, df, alg_colors, alg_order, outpath):
    """Chromosome-painting / ideogram view for one species: each of its
    scaffolds drawn as a horizontal bar, painted at each ortholog's
    position by that ortholog's ALG color - the collapsed, publication-
    style alternative to the full Oxford dot plot."""
    plt = _get_matplotlib()
    if plt is None:
        return False
    if df.empty:
        return False

    scaf_order = list(df["scaf"].value_counts().index)
    x_coord, x_bounds = build_plot_layout(df, "scaf", "pos", scaf_order)
    colors = [alg_colors.get(a, "#333333") for a in df["alg"]]
    present_algs = [a for a in alg_order if a in set(df["alg"])]

    n_scaf = len(x_bounds)
    fig_w = max(6.0, min(24.0, max(e for _, _, e in x_bounds) * 0.06 + 2)) + 1.6
    fig_h = max(2.0, n_scaf * 0.5 + 1)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    bar_height = 0.6
    for i, (scaf, s, e) in enumerate(x_bounds):
        y = n_scaf - 1 - i
        ax.add_patch(plt.matplotlib.patches.Rectangle((s - 0.5, y - bar_height / 2), e - s + 1, bar_height,
                                                        facecolor="white", edgecolor="0.15", linewidth=1.6,
                                                        zorder=1))
    sub_by_scaf = {scaf: (s, e) for scaf, s, e in x_bounds}
    for idx in df.index:
        scaf = df.loc[idx, "scaf"]
        i = scaf_order.index(scaf)
        y = n_scaf - 1 - i
        x = x_coord[idx]
        ax.add_patch(plt.matplotlib.patches.Rectangle((x - 0.5, y - bar_height / 2), 1, bar_height,
                                                        facecolor=colors[df.index.get_loc(idx)],
                                                        edgecolor="none", zorder=2))

    ax.set_yticks(range(n_scaf))
    ax.set_yticklabels([scaf for scaf, _, _ in reversed(x_bounds)], fontsize=8)
    ax.set_xlabel(f"{sp_name} position (ortholog rank)")
    ax.set_xlim(-1, max(e for _, _, e in x_bounds) + 1)
    ax.set_ylim(-1, n_scaf)
    ax.set_title(f"{sp_name} chromosome painting", fontsize=10)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    legend_handles = _alg_legend_handles(plt, alg_colors, present_algs)
    ax.legend(handles=legend_handles, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7,
              frameon=False, title="Linkage group", title_fontsize=8)
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return True


def plot_ribbon(taxon_order, species_data, alg_colors, outpath, highlight_wgd=False, title="Ribbon plot"):
    """Multi-track riparian/ribbon plot: one painted horizontal bar per
    taxon (in the given, user-controlled order), stacked top to bottom,
    with curved ribbons connecting the same ortholog (matched on
    (alg, ref_pos)) between each pair of ADJACENT tracks - the standard way
    to visualize where synteny is conserved (parallel ribbons) versus
    rearranged (crossing ribbons) across more than two genomes at once.

    Every taxon's track is independently rescaled to span the same total
    canvas width, so a taxon with fewer total orthologs doesn't leave a
    trailing blank gap - only relative proportions within a track are
    preserved, not absolute gene counts across tracks."""
    plt = _get_matplotlib()
    if plt is None:
        return False
    if len(taxon_order) < 2:
        print("[warn] --ribbon-order needs at least 2 taxa; skipping ribbon plot")
        return False

    SCAF_GAP = 3          # raw-unit gap between scaffold blocks, before rescaling
    CANVAS_WIDTH = 1000.0  # every track is independently rescaled to span this width
    TRACK_SPACING = 1.8    # vertical distance between adjacent track centers
    BAR_HEIGHT = 0.7 * 0.25  # thin chromosome bars; the vertical room this frees up between
                             # tracks is automatically used by the ribbons below, since their
                             # endpoints sit right at the (now much closer-together) bar edges

    raw = {}
    max_n_scaf = 0
    for sp in taxon_order:
        df = species_data[sp]
        if df.empty:
            print(f"[warn] {sp} has no data; skipping ribbon plot")
            return False
        scaf_order = list(df["scaf"].value_counts().index)
        x_coord, x_bounds = build_plot_layout(df, "scaf", "pos", scaf_order, gap=SCAF_GAP)
        raw[sp] = (x_coord, x_bounds, df)
        max_n_scaf = max(max_n_scaf, len(x_bounds))

    layouts = {}
    for sp in taxon_order:
        x_coord, x_bounds, df = raw[sp]
        raw_span = max(e for _, _, e in x_bounds) if x_bounds else 1
        scale = CANVAS_WIDTH / raw_span if raw_span > 0 else 1.0
        layouts[sp] = (x_coord * scale, [(g, s * scale, e * scale) for g, s, e in x_bounds], df, scale)

    n = len(taxon_order)
    fig_w = max(12.0, min(30.0, max_n_scaf * 0.35 + 6))
    fig_h = max(4.0, n * TRACK_SPACING * 0.9 + 1.5)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    track_y = {sp: (n - 1 - i) * TRACK_SPACING for i, sp in enumerate(taxon_order)}
    any_wgd_marked = False
    for sp in taxon_order:
        x_coord, x_bounds, df, scale = layouts[sp]
        y = track_y[sp]
        half_w = scale / 2
        for idx in df.index:
            x = x_coord[idx]
            color = alg_colors.get(df.loc[idx, "alg"], "#333333")
            ax.add_patch(plt.matplotlib.patches.Rectangle((x - half_w, y - BAR_HEIGHT / 2),
                                                            2 * half_w * 0.92, BAR_HEIGHT,
                                                            facecolor=color, edgecolor="none", zorder=3))
        # one outline per scaffold, so individual chromosomes are demarcated rather
        # than just the track as a whole
        for scaf_name, s, e in x_bounds:
            ax.add_patch(plt.matplotlib.patches.Rectangle((s - half_w, y - BAR_HEIGHT / 2),
                                                            (e - s) + 2 * half_w, BAR_HEIGHT,
                                                            facecolor="none", edgecolor="0.15",
                                                            linewidth=1.2, zorder=4))
            ax.text(s - half_w, y + BAR_HEIGHT / 2 + TRACK_SPACING * 0.02, scaf_name,
                    ha="left", va="bottom", fontsize=5, rotation=45, rotation_mode="anchor", zorder=5)
        if highlight_wgd and "orthology_type" in df.columns:
            wgd_idx = df.index[df["orthology_type"] == "wgd_inferred"]
            if len(wgd_idx):
                any_wgd_marked = True
                ax.scatter(x_coord[wgd_idx],
                           [y + BAR_HEIGHT / 2 + TRACK_SPACING * 0.012] * len(wgd_idx),
                           marker="v", s=9, c="black", linewidths=0, zorder=6)
        ax.text(-CANVAS_WIDTH * 0.02, y, sp, ha="right", va="center", fontsize=9)

    def curved_connector(x0, y0, x1, y1, color, alpha, linewidth):
        Path = plt.matplotlib.path.Path
        dy = y0 - y1
        verts = [(x0, y0), (x0, y0 - dy * 0.5), (x1, y1 + dy * 0.5), (x1, y1)]
        codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
        ax.add_patch(plt.matplotlib.patches.PathPatch(
            Path(verts, codes), facecolor="none", edgecolor=color, alpha=alpha,
            linewidth=linewidth, zorder=1))

    for sp_a, sp_b in zip(taxon_order[:-1], taxon_order[1:]):
        xa, _, df_a, _ = layouts[sp_a]
        xb, _, df_b, _ = layouts[sp_b]
        ya, yb = track_y[sp_a], track_y[sp_b]
        key_a = {(row.alg, row.ref_pos): idx for idx, row in df_a.iterrows() if pd.notna(row.ref_pos)}
        key_b = {(row.alg, row.ref_pos): idx for idx, row in df_b.iterrows() if pd.notna(row.ref_pos)}
        shared = set(key_a) & set(key_b)
        for key in shared:
            xi, xj = xa[key_a[key]], xb[key_b[key]]
            color = alg_colors.get(key[0], "#333333")
            curved_connector(xi, ya - BAR_HEIGHT / 2, xj, yb + BAR_HEIGHT / 2,
                              color, alpha=0.35, linewidth=1.0)

    ax.set_xlim(-CANVAS_WIDTH * 0.06, CANVAS_WIDTH * 1.02)
    ax.set_ylim(-TRACK_SPACING * 0.6, (n - 1) * TRACK_SPACING + TRACK_SPACING * 0.6)
    ax.axis("off")
    ax.set_title(title, fontsize=10)
    all_algs_in_plot = []
    for sp in taxon_order:
        for a in layouts[sp][2]["alg"]:
            if a not in all_algs_in_plot:
                all_algs_in_plot.append(a)
    legend_handles = _alg_legend_handles(plt, alg_colors, sorted(all_algs_in_plot))
    if any_wgd_marked:
        legend_handles = [plt.matplotlib.lines.Line2D(
            [], [], marker="v", linestyle="none", markersize=5, markerfacecolor="black",
            markeredgecolor="none", label="WGD-rescued ohnolog")] + legend_handles
    ax.legend(handles=legend_handles, bbox_to_anchor=(1.0, 1), loc="upper left", fontsize=7,
              frameon=False, title="Linkage group", title_fontsize=8)
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return True


# ==========================================================================
# Synapomorphies
# ==========================================================================

def find_synapomorphies(root, characters, matrix, mode):
    """For every internal node (clade), list characters that change state
    (increase) on the branch leading into that node. Returns:
    list of (node, leaf_set, [(char, from_state, to_state), ...])
    `node` is the actual Node object, so callers (e.g. the tree-figure
    renderer) can key off its identity as well as its leaf set."""
    cost_fn = cost_reversible if mode == "reversible" else cost_irreversible
    per_char_assignment = {}
    for char in characters:
        col = {sp: matrix[sp][char] for sp in matrix}
        per_char_assignment[char] = traceback_states(root, col, cost_fn, mode)

    results = []

    def walk(node, parent):
        if not node.is_leaf():
            events = []
            for char in characters:
                assigned = per_char_assignment[char]
                if parent is None:
                    continue
                ps = assigned[id(parent)]
                cs = assigned[id(node)]
                if cs != ps and cs != "?" and ps != "?":
                    events.append((char, ps, cs))
            results.append((node, frozenset(node.leaves()), events))
        for c in node.children:
            walk(c, node)

    walk(root, None)
    return results


# ==========================================================================
# Optional export to PAUP* (NEXUS) / TNT
# ==========================================================================

def _nexus_taxon_label(name):
    # NEXUS is happier without spaces/odd punctuation in bare labels
    return re.sub(r"[^A-Za-z0-9_.]", "_", name)


def write_nexus(path, species_list, characters, matrix, outgroup, constraint_newick=None):
    """Write a NEXUS file with the character matrix and a PAUP* command
    block that runs both an unordered (Fitch) and an irrev.up (Camin-Sokal
    / irreversible) heuristic search, honoring a constraint tree if given.
    Requires PAUP* (not installed by this script) to actually run."""
    labels = {sp: _nexus_taxon_label(sp) for sp in species_list}
    nchar = len(characters)
    state_of = {0: "0", 1: "1", 2: "2", "?": "?"}

    lines = []
    lines.append("#NEXUS")
    lines.append("")
    lines.append("BEGIN TAXA;")
    lines.append(f"    DIMENSIONS NTAX={len(species_list)};")
    lines.append("    TAXLABELS " + " ".join(labels[sp] for sp in species_list) + ";")
    lines.append("END;")
    lines.append("")
    lines.append("BEGIN CHARACTERS;")
    lines.append(f"    DIMENSIONS NCHAR={nchar};")
    lines.append('    FORMAT DATATYPE=STANDARD SYMBOLS="012" MISSING=? GAP=-;')
    lines.append("    CHARSTATELABELS")
    lab_lines = []
    for i, (a, b) in enumerate(characters, start=1):
        cname = _nexus_taxon_label(f"{a}_{b}_fusion")
        lab_lines.append(f"        {i} {cname}")
    lines.append(",\n".join(lab_lines) + ";")
    lines.append("    MATRIX")
    width = max(len(l) for l in labels.values()) + 2
    for sp in species_list:
        row = "".join(state_of[matrix[sp][c]] for c in characters)
        lines.append(f"    {labels[sp].ljust(width)}{row}")
    lines.append("    ;")
    lines.append("END;")
    lines.append("")

    constraint_name = None
    if constraint_newick:
        # translate taxon names in the constraint newick to nexus-safe labels
        translated = constraint_newick.strip()
        if translated.endswith(";"):
            translated = translated[:-1]
        for sp, lab in labels.items():
            translated = re.sub(rf"\b{re.escape(sp)}\b", lab, translated)
        constraint_name = "backbone"
        lines.append("[ backbone constraint, defined below via the PAUP CONSTRAINTS command ]")
        lines.append("")

    og = labels[outgroup]
    lines.append("BEGIN PAUP;")
    lines.append(f"    outgroup {og};")
    lines.append("    set criterion=parsimony;")
    if constraint_name:
        lines.append(f"    constraints {constraint_name} (MONOPHYLY) = {translated};")
    lines.append("")
    lines.append("    [--- unordered / Fitch (\"reversible\") search ---]")
    lines.append("    ctype unord: all;")
    if constraint_name:
        lines.append(f"    hsearch addseq=random nreps=100 enforce=yes converse=no "
                      f"constraints={constraint_name};")
    else:
        lines.append("    hsearch addseq=random nreps=100;")
    lines.append('    savetrees file=paup_reversible.tre replace=yes;')
    lines.append("")
    lines.append("    [--- irrev.up / Camin-Sokal (\"irreversible\") search ---]")
    lines.append("    [state 0 = unfused, 1 = fused-clean, 2 = fused-mixed;]")
    lines.append("    [irrev.up permits only increasing state changes, matching]")
    lines.append("    [this script's --outgroup-rooted irreversible model.]")
    lines.append("    ctype irrev.up: all;")
    if constraint_name:
        lines.append(f"    hsearch addseq=random nreps=100 enforce=yes converse=no "
                      f"constraints={constraint_name};")
    else:
        lines.append("    hsearch addseq=random nreps=100;")
    lines.append('    savetrees file=paup_irreversible.tre replace=yes;')
    lines.append("END;")
    lines.append("")

    with open(path, "w") as f:
        f.write("\n".join(lines))


def write_tnt(path, species_list, characters, matrix):
    """Write a Hennig86/TNT-format data file plus a minimal search script.
    TNT has no simple built-in irreversible/Dollo character type, so this
    covers the standard (reversible, Fitch/Wagner) search only - use it as
    a cross-check of the reversible analysis, not the irreversible one.
    Requires TNT (not installed by this script) to actually run."""
    labels = {sp: _nexus_taxon_label(sp) for sp in species_list}
    nchar = len(characters)
    state_of = {0: "0", 1: "1", 2: "2", "?": "?"}
    width = max(len(l) for l in labels.values()) + 2

    lines = []
    lines.append("xread")
    lines.append("'ALG fusion/mixing characters (0=unfused 1=fused-clean 2=fused-mixed)'")
    lines.append(f"{nchar} {len(species_list)}")
    for sp in species_list:
        row = "".join(state_of[matrix[sp][c]] for c in characters)
        lines.append(f"{labels[sp].ljust(width)}{row}")
    lines.append(";")
    lines.append("")
    lines.append("cc - .;")           # all characters non-additive (Fitch), TNT default
    lines.append("hold 10000;")
    lines.append("mult=tbr 100;")     # 100 random-addition + TBR replicates
    lines.append("bbreak=tbr;")
    lines.append("export - tnt_reversible.tre;")
    lines.append("zzz;")              # end script

    with open(path, "w") as f:
        f.write("\n".join(lines))


def try_run_external(binary, script_path, outdir, label, log_name):
    """Best-effort optional invocation of an external phylogenetics program.
    Never raises - prints a warning and returns False if the binary is
    missing or fails, so the rest of the script is unaffected."""
    import subprocess
    if shutil_which(binary) is None:
        print(f"[warn] '{binary}' not found on PATH - skipping {label} run. "
              f"The exported script at {script_path} can still be run manually "
              f"(or install {label} / pass --{'paup' if 'paup' in label.lower() else 'tnt'}-bin).")
        return False
    log_path = os.path.join(outdir, log_name)
    try:
        with open(log_path, "w") as logf:
            result = subprocess.run(
                [binary, script_path], cwd=outdir,
                stdout=logf, stderr=subprocess.STDOUT, timeout=1800,
            )
        print(f"[info] ran {label} (exit code {result.returncode}); log written to {log_path}")
        print(f"[info] check {outdir}/ for any tree file(s) {label} was told to save "
              f"- this script does not parse them back in automatically.")
        return True
    except Exception as e:
        print(f"[warn] {label} run failed: {e}")
        return False


def shutil_which(name):
    import shutil
    return shutil.which(name)


# ==========================================================================
# Main
# ==========================================================================

def main():
    ap = argparse.ArgumentParser(
        description="ALG fusion/mixing macrosynteny characters -> parsimony tree search",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument("--indir", required=True, help="directory containing .rbh files")
    ap.add_argument("--reference-prefix", default=None,
                     help="the shared ancestral/ALG reference column prefix used across your "
                          ".rbh files (e.g. 'BCnSSimakov2022' or 'odp_LGs_from_Avul_Ixod_Nlon_Varr'). "
                          "Auto-detected from the headers if omitted - only needed if autodetection "
                          "warns or picks the wrong one.")
    ap.add_argument("--outdir", default="alg_parsimony_out", help="output directory")
    ap.add_argument("--alpha", type=float, default=0.05, help="whole_FET significance cutoff")
    ap.add_argument("--min-genes", type=int, default=3,
                     help="min shared orthologs to call an ALG present on a scaffold / in a genome")
    ap.add_argument("--mixing-threshold", type=float, default=0.3,
                     help="fraction of order-switches along a fused scaffold above which the "
                          "fusion is called 'mixed' rather than 'clean'")
    ap.add_argument("--outgroup", default=None,
                     help="taxon to root the irreversible analysis on (auto-detected if omitted)")
    ap.add_argument("--max-exhaustive", type=int, default=8,
                     help="max number of taxa for exact tree enumeration; above this, "
                          "a heuristic NNI search is used")
    ap.add_argument("--nni-restarts", type=int, default=10,
                     help="random restarts for the heuristic search")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--constraint", default=None,
                     help="optional backbone constraint: a Newick string, or a path to a file "
                          "containing one. May cover some, all, or none of the taxa; polytomies "
                          "are left for the ALG characters to resolve. Applied to both runs.")
    ap.add_argument("--export-nexus", action="store_true",
                     help="write character_matrix.nex (PAUP* NEXUS format) to --outdir")
    ap.add_argument("--export-tnt", action="store_true",
                     help="write character_matrix.tnt (TNT/Hennig86 format) to --outdir")
    ap.add_argument("--run-paup", action="store_true",
                     help="also try to invoke PAUP* on the exported NEXUS file (implies --export-nexus)")
    ap.add_argument("--run-tnt", action="store_true",
                     help="also try to invoke TNT on the exported data file (implies --export-tnt)")
    ap.add_argument("--paup-bin", default="paup", help="PAUP* executable name/path (default: paup)")
    ap.add_argument("--tnt-bin", default="tnt", help="TNT executable name/path (default: tnt)")
    ap.add_argument("--plot-dotplots", action="store_true",
                     help="write an Oxford-style synteny dot plot per species to --outdir/dotplots/ "
                          "(requires matplotlib)")
    ap.add_argument("--plot-trees", action="store_true",
                     help="write a rendered tree figure per run, with fusion/mixing "
                          "synapomorphies labeled on the branch where they originate "
                          "(requires matplotlib)")
    ap.add_argument("--plot-matrix", action="store_true",
                     help="write a taxa x characters heatmap of the character matrix to --outdir/ "
                          "(requires matplotlib)")
    ap.add_argument("--plot-ideograms", action="store_true",
                     help="write a chromosome-painting/ideogram plot per species to --outdir/ideograms/ "
                          "- the collapsed, publication-style alternative to the full dot plot "
                          "(requires matplotlib)")
    ap.add_argument("--plot-ribbon", action="store_true",
                     help="write a multi-track ribbon/riparian plot connecting orthologs between "
                          "adjacent taxa to --outdir/ribbon.<ext> (requires matplotlib). Taxon order "
                          "is alphabetical unless --ribbon-order or --ribbon-order-from-tree is given.")
    ap.add_argument("--ribbon-order", nargs="+", default=None,
                     help="explicit top-to-bottom taxon order for the ribbon plot (implies "
                          "--plot-ribbon). Needs at least 2 taxa and may name any subset of them - "
                          "restricting the plot to a handful of taxa of interest is often much more "
                          "readable than showing all of them. Takes priority over "
                          "--ribbon-order-from-tree if both are given.")
    ap.add_argument("--ribbon-order-from-tree", choices=["reversible", "irreversible"], default=None,
                     help="derive the ribbon plot's taxon order from the computed tree instead of "
                          "specifying it by hand (implies --plot-ribbon): a ladderized leaf order, so "
                          "every clade forms a contiguous block and phylogenetically close taxa end up "
                          "on adjacent tracks. Doesn't minimize ribbon crossings from genuine genome "
                          "rearrangement (that reflects real biology), but does eliminate crossings that "
                          "are only an artefact of an arbitrary taxon order - the effect is usually "
                          "biggest on datasets with substantial genome reorganization, where an "
                          "unrelated taxon order can make unrelated genomes look artificially tangled "
                          "next to each other. Choose which criterion's tree to use ('reversible' or "
                          "'irreversible'); ignored if --ribbon-order is also given.")
    ap.add_argument("--plot-format", choices=["pdf", "png", "svg"], default="pdf",
                     help="file format for every plot this script writes (default: pdf)")
    ap.add_argument("--highlight-wgd", action="store_true",
                     help="highlight WGD-rescued orthologs (rows whose orthology_type column, written "
                          "by genomes_to_rbh.py --detect-wgd, is 'wgd_inferred'). Adds a companion "
                          "'<taxon>.wgd.<ext>' dot plot per affected taxon with those points drawn at "
                          "double diameter and ringed, marks them on the main --plot-ribbon plot's "
                          "chromosome bars, writes a dedicated 'ribbon_wgd/<taxon>.wgd.<ext>' ribbon "
                          "plot per affected taxon (needs --plot-ribbon) using the same taxon set as "
                          "the main ribbon plot - i.e. the consensus/core species plus anything named "
                          "in --ribbon-order, with the affected taxon itself guaranteed to be included "
                          "- and writes ohnolog_summary.tsv counting rescued vs. direct orthologs per "
                          "taxon, scaffold and linkage group. Has no effect on .rbh files lacking that "
                          "column (i.e. not produced with --detect-wgd).")
    ap.add_argument("--bootstrap", type=int, default=0,
                     help="nonparametric bootstrap replicates (default: 0 = off). Resamples "
                          "characters with replacement, reruns the search each time, and reports "
                          "the %% of replicates recovering each clade of the reported tree - both as "
                          "internal-node support labels embedded directly in the .nwk tree file and "
                          "as a separate bootstrap_<run>.tsv. Cost scales with the number of "
                          "replicates times the cost of one full search - start small (e.g. 100) and "
                          "increase if runtime allows.")
    ap.add_argument("--bootstrap-restarts", type=int, default=3,
                     help="NNI restarts per bootstrap replicate when the heuristic search applies "
                          "(default: 3 - deliberately lighter than --nni-restarts, since this runs "
                          "once per replicate)")
    ap.add_argument("--permutation-test", type=int, default=0,
                     help="permutation/signal-test replicates (default: 0 = off). Independently "
                          "shuffles which taxon has which state within every character, reruns the "
                          "search each time, and reports what fraction of shuffled replicates score "
                          "at least as well as the real data - a low fraction means the real "
                          "fusion/mixing data carry more structure than chance.")
    ap.add_argument("--permutation-restarts", type=int, default=3,
                     help="NNI restarts per permutation replicate when the heuristic search applies "
                          "(default: 3)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    rng = random.Random(args.seed)

    # ---- 1. parse ----
    species_data, reference_prefix = parse_rbh_directory(args.indir, reference_prefix_override=args.reference_prefix)
    species_list = sorted(species_data.keys())
    if len(species_list) < 3:
        sys.exit("Need at least 3 species/taxa to build a tree.")
    print(f"[info] parsed {len(species_list)} taxa: {', '.join(species_list)}")

    ext = args.plot_format

    if args.ribbon_order:
        args.plot_ribbon = True
        missing = set(args.ribbon_order) - set(species_list)
        if missing:
            sys.exit(f"--ribbon-order names taxa not found among the parsed .rbh species: "
                      f"{sorted(missing)}. Known taxa: {species_list}")
        dupes = [t for t in set(args.ribbon_order) if args.ribbon_order.count(t) > 1]
        if dupes:
            sys.exit(f"--ribbon-order names {sorted(dupes)} more than once; each taxon may appear "
                      f"at most once.")
        if len(args.ribbon_order) < 2:
            sys.exit("--ribbon-order needs at least 2 taxa to draw ribbons between.")
        omitted = [t for t in species_list if t not in set(args.ribbon_order)]
        if omitted:
            print(f"[info] ribbon plot restricted to {len(args.ribbon_order)} of {len(species_list)} "
                  f"taxa; not shown: {', '.join(omitted)}")
    if args.ribbon_order_from_tree:
        args.plot_ribbon = True
        if args.ribbon_order:
            print(f"[info] both --ribbon-order and --ribbon-order-from-tree were given; "
                  f"using the explicit --ribbon-order and ignoring --ribbon-order-from-tree")

    if args.plot_dotplots or args.plot_ideograms or args.plot_ribbon:
        alg_colors = get_global_alg_colors(species_data)
        alg_order = sorted(alg_colors.keys())
        ref_label = reference_prefix or "reference"

    if args.plot_dotplots:
        dotplot_dir = os.path.join(args.outdir, "dotplots")
        os.makedirs(dotplot_dir, exist_ok=True)
        for sp in species_list:
            out_path = os.path.join(dotplot_dir, f"{sp}.{ext}")
            if plot_dotplot(sp, ref_label, species_data[sp], alg_colors, alg_order, out_path):
                print(f"[info] wrote {out_path}")
            if args.highlight_wgd:
                wgd_path = os.path.join(dotplot_dir, f"{sp}.wgd.{ext}")
                if plot_dotplot_wgd(sp, ref_label, species_data[sp], alg_colors, alg_order, wgd_path):
                    print(f"[info] wrote {wgd_path}")

    if args.highlight_wgd:
        summary_path = os.path.join(args.outdir, "ohnolog_summary.tsv")
        n_taxa_with_wgd = write_ohnolog_summary(species_data, summary_path)
        if n_taxa_with_wgd:
            print(f"[info] wrote {summary_path} ({n_taxa_with_wgd} taxon/taxa with WGD-rescued orthologs)")
        else:
            print("[info] --highlight-wgd: no WGD-rescued orthologs found in any .rbh file "
                  "(was genomes_to_rbh.py run with --detect-wgd?); nothing to highlight or summarize")

    if args.plot_ideograms:
        ideogram_dir = os.path.join(args.outdir, "ideograms")
        os.makedirs(ideogram_dir, exist_ok=True)
        for sp in species_list:
            out_path = os.path.join(ideogram_dir, f"{sp}.{ext}")
            if plot_ideogram(sp, species_data[sp], alg_colors, alg_order, out_path):
                print(f"[info] wrote {out_path}")

    # --plot-ribbon itself is deferred until after the tree search (below) when
    # --ribbon-order-from-tree is requested, since it needs the computed tree's
    # leaf order; --ribbon-order and the default (alphabetical) order don't need
    # the tree and could run here, but are handled in the same place for simplicity.

    # ---- 2/3. fusion + mixing calls ----
    fusion_state, alg_total_counts, scaffold_rows = build_fusion_calls(
        species_data, args.alpha, args.min_genes, args.mixing_threshold
    )

    with open(os.path.join(args.outdir, "alg_scaffold_calls.tsv"), "w") as f:
        f.write("species\tscaffold\talgA\talgB\tstate\n")
        for sp, scaf, a, b, state in scaffold_rows:
            f.write(f"{sp}\t{scaf}\t{a}\t{b}\t{STATE_NAME[state]}\n")

    # ---- 4. character matrix ----
    characters, matrix = build_character_matrix(species_list, fusion_state, alg_total_counts, args.min_genes)
    if not characters:
        sys.exit("No fusion characters were called at all - check --alpha / --min-genes, "
                  "or your genomes may simply retain the ancestral (unfused) linkage groups.")
    print(f"[info] {len(characters)} fusion/mixing characters assembled")

    with open(os.path.join(args.outdir, "character_matrix.tsv"), "w") as f:
        f.write("taxon\t" + "\t".join(f"{a}+{b}" for a, b in characters) + "\n")
        for sp in species_list:
            row = [str(matrix[sp][c]) for c in characters]
            f.write(sp + "\t" + "\t".join(row) + "\n")

    if args.plot_matrix:
        out_path = os.path.join(args.outdir, f"character_matrix_heatmap.{ext}")
        if plot_character_matrix_heatmap(species_list, characters, matrix, out_path):
            print(f"[info] wrote {out_path}")

    # ---- choose outgroup for irreversible run ----
    if args.outgroup:
        if args.outgroup not in species_list:
            sys.exit(f"--outgroup '{args.outgroup}' not found among taxa {species_list}")
        outgroup = args.outgroup
        print(f"[info] using user-specified outgroup: {outgroup}")
    else:
        def derived_count(sp):
            return sum(1 for c in characters if matrix[sp][c] not in (0, "?"))
        outgroup = min(species_list, key=derived_count)
        print(f"[info] auto-selected outgroup (fewest derived states): {outgroup} "
              f"- override with --outgroup if this is wrong")

    # ---- optional backbone constraint ----
    n = len(species_list)
    constraint = None
    constraint_newick_text = None
    if args.constraint:
        if os.path.isfile(args.constraint):
            with open(args.constraint) as f:
                constraint_newick_text = f.read()
        else:
            constraint_newick_text = args.constraint
        constraint_root = parse_newick(constraint_newick_text)
        constraint_taxa = set(constraint_root.leaves())
        unknown = constraint_taxa - set(species_list)
        if unknown:
            hint = ""
            if reference_prefix in unknown:
                hint = (
                    f"\n\nNote: '{reference_prefix}' is the detected ALG/reference column prefix for this "
                    f"dataset, not a taxon with its own .rbh file - in single-reference mode "
                    f"(genomes_to_rbh.py --reference-species {reference_prefix}) the reference genome is "
                    f"deliberately excluded from the taxon set (there's nothing to compare it to), so it "
                    f"can never appear in a --constraint tree. Remove it from the constraint, or rebuild "
                    f"the .rbh files with --consensus-species instead if you want every genome (including "
                    f"this one) to be an independent taxon."
                )
            sys.exit(f"--constraint mentions taxa not found among the parsed .rbh species: "
                      f"{sorted(unknown)}. Known taxa: {species_list}{hint}")
        required_clades = required_clades_from_constraint(constraint_root)
        # a clade spanning the *entire* known taxon set is trivially always
        # true (no bipartition can ever put "everything" on one side) - drop it
        required_clades = {c for c in required_clades if len(c) < len(species_list)}
        constraint = {
            "root_node": constraint_root,
            "taxa": constraint_taxa,
            "required_clades": required_clades,
        }
        print(f"[info] backbone constraint loaded: {len(required_clades)} required clade(s) "
              f"covering {len(constraint_taxa)}/{n} taxa")
        for clade in sorted(required_clades, key=len):
            print(f"         - {{{', '.join(sorted(clade))}}}")

    # ---- 5. tree search: reversible / irreversible parsimony ----
    exhaustive = n <= args.max_exhaustive

    def run(mode):
        if exhaustive:
            print(f"[info] exhaustive search over all unrooted topologies for '{mode}' run "
                  f"({n} taxa) ...")
            score, trees = best_tree_exhaustive(species_list, characters, matrix, mode, outgroup,
                                                 constraint=constraint)
        else:
            print(f"[info] heuristic NNI search for '{mode}' run ({n} taxa, "
                  f"{args.nni_restarts} restarts) ...")
            score, trees = best_tree_heuristic(species_list, characters, matrix, mode, outgroup,
                                                args.nni_restarts, rng, constraint=constraint)
        return score, trees

    rev_score, rev_trees = run("reversible")
    irr_score, irr_trees = run("irreversible")

    report_lines = []
    report_lines.append("ALG fusion/mixing parsimony report")
    report_lines.append("=" * 60)
    report_lines.append(f"Taxa ({n}): {', '.join(species_list)}")
    report_lines.append(f"Characters: {len(characters)}")
    report_lines.append(f"Outgroup used for irreversible rooting: {outgroup}")
    if constraint:
        report_lines.append(f"Backbone constraint: {len(constraint['required_clades'])} required "
                             f"clade(s) enforced (from {args.constraint})")
        for clade in sorted(constraint["required_clades"], key=len):
            report_lines.append(f"    - {{{', '.join(sorted(clade))}}}")
    report_lines.append("")

    roots_by_mode = {}

    def emit_run(name, mode, score, trees, fname_prefix):
        report_lines.append("-" * 60)
        report_lines.append(f"{name} parsimony")
        report_lines.append(f"  best score: {score}")
        report_lines.append(f"  number of equally parsimonious trees found: {len(trees)}")
        best = trees[0]
        root = root_at_leaf(best, outgroup)
        roots_by_mode[fname_prefix] = root

        syn = find_synapomorphies(root, characters, matrix, mode)
        # Clades of size >= n_taxa-1 (the whole taxon set, and "everyone except the
        # outgroup") are automatically implied by the *choice* of outgroup rooting
        # itself, for any topology - they can never match a real bipartition (see
        # bipartitions()), so bootstrap "support" for them is always spuriously 0%
        # rather than meaningful. Exclude them; only clades the data could actually
        # confirm or refute are worth a support value.
        target_clades = [leafset for _, leafset, _ in syn if len(leafset) < len(species_list) - 1]

        # Bootstrap is computed here, before the Newick string is built, so support
        # values (when requested) can be embedded as internal-node labels in the
        # .nwk file itself, not just reported in the text summary and .tsv.
        support = None
        if args.bootstrap > 0:
            print(f"[info] running {args.bootstrap} bootstrap replicates for '{name}' "
                  f"({'exhaustive' if exhaustive else 'heuristic'} search per replicate) ...")
            support = bootstrap_support(species_list, characters, matrix, mode, outgroup, target_clades,
                                          args.bootstrap, exhaustive, args.bootstrap_restarts, rng, constraint)

        nwk = to_newick(root, support=support) + ";"
        ascii_tree = "\n".join(to_ascii(root))
        report_lines.append("  newick (first best tree, rooted on outgroup for display" +
                             (", internal-node labels = bootstrap % support" if support else "") + "):")
        report_lines.append("    " + nwk)
        report_lines.append("  ascii:")
        for line in ascii_tree.splitlines():
            report_lines.append("    " + line)

        with open(os.path.join(args.outdir, f"tree_{fname_prefix}.nwk"), "w") as f:
            f.write(nwk + "\n")
        with open(os.path.join(args.outdir, f"tree_{fname_prefix}.txt"), "w") as f:
            f.write(ascii_tree + "\n")

        syn_lines = []
        for node, leafset, events in syn:
            if not events:
                continue
            taxa_str = ", ".join(sorted(leafset))
            syn_lines.append(f"Clade [{taxa_str}]:")
            for char, ps, cs in events:
                a, b = char
                syn_lines.append(
                    f"    {a}+{b} fusion: {STATE_NAME[ps]} -> {STATE_NAME[cs]}"
                )
        with open(os.path.join(args.outdir, f"synapomorphies_{fname_prefix}.txt"), "w") as f:
            f.write("\n".join(syn_lines) + "\n" if syn_lines else "(no state changes found)\n")

        if args.plot_trees:
            tree_fig_path = os.path.join(args.outdir, f"tree_{fname_prefix}.{ext}")
            ok = render_tree_figure(root, characters, matrix, mode, tree_fig_path,
                                     title=f"{name} (score={score})")
            if ok:
                print(f"[info] wrote {tree_fig_path}")

        report_lines.append("  synapomorphies (character-state changes mapped to clades):")
        if syn_lines:
            for line in syn_lines:
                report_lines.append("    " + line)
        else:
            report_lines.append("    (none)")

        # ---- Consistency Index / Retention Index ----
        per_char_ci_ri, ensemble = compute_ci_ri(root, characters, matrix, mode)
        report_lines.append(f"  ensemble CI: {ensemble['ci']:.3f}   ensemble RI: {ensemble['ri']:.3f}")
        with open(os.path.join(args.outdir, f"ci_ri_{fname_prefix}.tsv"), "w") as f:
            f.write("character\tm\ts\tg\tCI\tRI\n")
            for char in characters:
                a, b = char
                v = per_char_ci_ri[char]
                f.write(f"{a}+{b}\t{v['m']}\t{v['s']}\t{v['g']}\t{v['ci']:.4f}\t{v['ri']:.4f}\n")

        # ---- bootstrap support (already computed above; report it here) ----
        if support is not None:
            with open(os.path.join(args.outdir, f"bootstrap_{fname_prefix}.tsv"), "w") as f:
                f.write("clade\tsupport_fraction\n")
                for clade in target_clades:
                    f.write("{" + ", ".join(sorted(clade)) + "}\t" + f"{support[clade]:.3f}\n")
            report_lines.append(f"  bootstrap support ({args.bootstrap} replicates):")
            for clade in target_clades:
                taxa_str = ", ".join(sorted(clade))
                report_lines.append(f"    [{taxa_str}]: {support[clade] * 100:.0f}%")

        # ---- permutation / signal test ----
        if args.permutation_test > 0:
            print(f"[info] running {args.permutation_test} permutation-test replicates for '{name}' "
                  f"({'exhaustive' if exhaustive else 'heuristic'} search per replicate) ...")
            permuted_scores, p_value = permutation_test(
                species_list, characters, matrix, mode, outgroup, args.permutation_test, exhaustive,
                args.permutation_restarts, rng, constraint, score)
            mean_permuted = sum(permuted_scores) / len(permuted_scores)
            report_lines.append(f"  permutation test ({args.permutation_test} replicates): observed "
                                 f"score {score}, mean shuffled score {mean_permuted:.2f}, "
                                 f"p = {p_value:.4f} (fraction of shuffled replicates scoring <= observed)")
            with open(os.path.join(args.outdir, f"permutation_test_{fname_prefix}.tsv"), "w") as f:
                f.write("replicate\tpermuted_score\n")
                for i, s in enumerate(permuted_scores, start=1):
                    f.write(f"{i}\t{s}\n")

        report_lines.append("")

    emit_run("Reversible", "reversible", rev_score, rev_trees, "reversible")
    emit_run("Irreversible (fusion/mixing cannot revert)", "irreversible", irr_score, irr_trees,
              "irreversible")

    if args.plot_ribbon:
        if args.ribbon_order:
            ribbon_order = args.ribbon_order
        elif args.ribbon_order_from_tree:
            tree_root = roots_by_mode[args.ribbon_order_from_tree]
            ladderize(tree_root)
            ribbon_order = [leaf.name for leaf in collect_leaf_nodes(tree_root)]
            print(f"[info] taxon order for ribbon plot derived from the {args.ribbon_order_from_tree} "
                  f"tree (ladderized, rooted on {outgroup}): {', '.join(ribbon_order)}")
        else:
            ribbon_order = species_list
        out_path = os.path.join(args.outdir, f"ribbon.{ext}")
        if plot_ribbon(ribbon_order, species_data, alg_colors, out_path,
                        highlight_wgd=args.highlight_wgd):
            print(f"[info] wrote {out_path}")

        if args.highlight_wgd:
            wgd_taxa = [sp for sp in species_list
                        if "orthology_type" in species_data[sp].columns
                        and (species_data[sp]["orthology_type"] == "wgd_inferred").any()]
            if wgd_taxa:
                wgd_ribbon_dir = os.path.join(args.outdir, "ribbon_wgd")
                os.makedirs(wgd_ribbon_dir, exist_ok=True)
            for sp in wgd_taxa:
                # Same taxon set as the main ribbon plot - the reference genome's ALG
                # coordinate system underlies every track regardless of mode, and in
                # consensus mode the core/consensus species are literal taxa in
                # ribbon_order - plus anything else the user named there. `sp` itself
                # is guaranteed to be included even if --ribbon-order omitted it.
                sp_ribbon_order = ribbon_order if sp in ribbon_order else [sp] + list(ribbon_order)
                wgd_ribbon_path = os.path.join(wgd_ribbon_dir, f"{sp}.wgd.{ext}")
                if plot_ribbon(sp_ribbon_order, species_data, alg_colors, wgd_ribbon_path,
                                highlight_wgd=True,
                                title=f"{sp} vs. reference/consensus context "
                                      f"\u2014 WGD-rescued orthologs highlighted"):
                    print(f"[info] wrote {wgd_ribbon_path}")

    report_text = "\n".join(report_lines)
    with open(os.path.join(args.outdir, "report.txt"), "w") as f:
        f.write(report_text + "\n")

    print()
    print(report_text)
    print()

    # ---- optional TNT / PAUP* export and invocation ----
    if args.export_nexus or args.run_paup:
        nexus_path = os.path.join(args.outdir, "character_matrix.nex")
        write_nexus(nexus_path, species_list, characters, matrix, outgroup,
                    constraint_newick=constraint_newick_text)
        print(f"[info] wrote NEXUS file for PAUP*: {nexus_path}")
        if args.run_paup:
            try_run_external(args.paup_bin, os.path.basename(nexus_path), args.outdir,
                              "PAUP*", "paup_run.log")

    if args.export_tnt or args.run_tnt:
        tnt_path = os.path.join(args.outdir, "character_matrix.tnt")
        write_tnt(tnt_path, species_list, characters, matrix)
        print(f"[info] wrote TNT data/command file: {tnt_path} "
              f"(reversible search only - see file header)")
        if args.run_tnt:
            try_run_external(args.tnt_bin, os.path.basename(tnt_path), args.outdir,
                              "TNT", "tnt_run.log")

    print(f"[info] all output written to {args.outdir}/")


if __name__ == "__main__":
    main()
