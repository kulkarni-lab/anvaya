#!/usr/bin/env python3
"""
make_toy_dataset.py
====================

Generates a small, fully-synthetic set of genome assemblies (CDS + protein
+ GFF3 per taxon) with a KNOWN, designed macrosynteny history, for testing
genomes_to_rbh.py / alg_parsimony.py end to end - including against a real
`diamond` installation, not a stubbed one.

Unlike the various one-off test fixtures built during this project's own
development (which used exact-ID-matching tricks that only work against a
fake stand-in for diamond), every gene here gets a real, independently
mutated protein sequence descended from a shared "ancestral" sequence, so
that genuine reciprocal-best-hit search finds the intended orthology
relationships for the reasons a real ortholog search would - not because
the IDs happen to match. See "SEQUENCE MODEL" below for the substitution
rate this relies on, and --verify-only to re-run the sanity check that
justifies it without regenerating the dataset.

--------------------------------------------------------------------------
WHAT IT BUILDS
--------------------------------------------------------------------------
11 ancestral linkage groups (ALGs): A-H (the "core" set, 6 genes each),
I/J (dedicated to demonstrating the mixed-fusion state), and N (a single
gene dedicated to demonstrating --min-scaffold-length).

8 taxa, each a deliberate test case:

  T1            outgroup / natural --reference-species or core-set anchor;
                every ALG unfused, one gene of ALG N on its own tiny
                (~300 bp) scaffold
  T2            ALG G+H fused, clean            (independent of T5, below)
  T3            ALG A+B fused, clean; E+F fused, clean
  T4            same fusions as T3 (clade partner)
  T5            ALG A+B fused, clean; C+D fused, clean; G+H fused, clean
  T6            same A+B/C+D fusions as T5 (clade partner); also carries
                the ALG-N ortholog, on its own tiny scaffold, matched with
                T1's copy
  T7_wgd        every one of ALGs A-H present as TWO independently-mutated
                full copies on separate scaffolds (a clean, fusion-free
                whole-genome duplication) - for --detect-wgd
  T8_mixed      ALG I+J fused AND interleaved (mixed) - the one designed
                "mixed" character in the dataset, isolated from everything
                else so it doesn't disturb the T1-T6 signal

T1-T6 alone reproduce the true tree (T1,(T2,((T3,T4),(T5,T6)))) via three
clean, nested synapomorphies (A+B, E+F, C+D) plus one deliberately
homoplastic character (G+H, arising independently in T2 and T5) - this is
the same design used to validate CI/RI, bootstrap and the permutation test
during this project's own development, so the exact expected numbers are
already known (see MANIFEST.md, written alongside the data).

--------------------------------------------------------------------------
SEQUENCE MODEL
--------------------------------------------------------------------------
Each ancestral gene gets one randomly-generated 120-aa "master" sequence.
Every taxon's copy of that gene is the master sequence with each residue
independently mutated at rate --mutation-rate (default 0.06, i.e. ~94%
identity to the master per copy). Two taxa's independently-mutated copies
of the *same* ancestral gene are then ~88% identical to each other -
trivially detectable by a real ortholog search - while two genes from
*different* ancestral origins are ~5% identical (~random chance for a
20-letter alphabet), nowhere near enough to cross-hit. --verify-only
recomputes this separation without regenerating any files.

--------------------------------------------------------------------------
USAGE
--------------------------------------------------------------------------
    python3 make_toy_dataset.py --outdir toy_data

Then, e.g.:

    python3 genomes_to_rbh.py --input-dir toy_data --reference-species T1 \\
        --rbh-outdir toy_rbh --parsimony-outdir toy_parsimony \\
        --forward --min-genes 3 --alpha 0.05 --outgroup T1 \\
        --plot-dotplots --plot-trees

See MANIFEST.md (written into --outdir) for the full ground truth and more
example commands, including consensus mode, --detect-wgd and
--min-scaffold-length.
"""

import argparse
import random
import sys
from pathlib import Path

AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")
CODON_OF = {
    "A": "GCT", "R": "CGT", "N": "AAT", "D": "GAT", "C": "TGT", "Q": "CAA",
    "E": "GAA", "G": "GGT", "H": "CAT", "I": "ATT", "L": "CTT", "K": "AAA",
    "M": "ATG", "F": "TTT", "P": "CCT", "S": "TCT", "T": "ACT", "W": "TGG",
    "Y": "TAT", "V": "GTT",
}
STOP_CODON = "TAA"
PROTEIN_LEN = 120


def random_sequence(rng):
    return "".join(rng.choice(AMINO_ACIDS) for _ in range(PROTEIN_LEN))


def mutate(seq, rate, rng):
    out = []
    for aa in seq:
        if rng.random() < rate:
            choices = [a for a in AMINO_ACIDS if a != aa]
            out.append(rng.choice(choices))
        else:
            out.append(aa)
    return "".join(out)


def reverse_translate(aa_seq):
    return "".join(CODON_OF[aa] for aa in aa_seq) + STOP_CODON


def hamming_identity(a, b):
    matches = sum(1 for x, y in zip(a, b) if x == y)
    return matches / max(len(a), len(b))


# ==========================================================================
# ALG / taxon design
# ==========================================================================

CORE_ALGS = list("ABCDEFGH")
N_CORE_GENES = 6
N_MIXED_GENES = 6

# every ALG's total gene indices, keyed by ALG letter
ALG_GENE_COUNT = {alg: N_CORE_GENES for alg in CORE_ALGS}
ALG_GENE_COUNT["I"] = N_MIXED_GENES
ALG_GENE_COUNT["J"] = N_MIXED_GENES
ALG_GENE_COUNT["N"] = 1

SCAF_GENE_SPACING = 8000     # bp between consecutive genes on a normal scaffold
SCAF_START = 2000
TINY_SCAF_POS = 500          # the one gene on a "noise" scaffold


def unfused_layout(taxon, algs, start_scaf_suffix=""):
    """One scaffold per ALG, genes in index order - the 'ancestral', fully
    unfused state. Returns list of (alg, gene_idx, scaffold, pos)."""
    rows = []
    for alg in algs:
        scaf = f"{taxon}_{alg}{start_scaf_suffix}"
        pos = SCAF_START
        for g in range(ALG_GENE_COUNT[alg]):
            rows.append((alg, g, scaf, pos))
            pos += SCAF_GENE_SPACING
    return rows


def fused_layout(taxon, scaf_name, alg_blocks, mixed=False):
    """One scaffold carrying two or more ALGs. alg_blocks: list of ALG
    letters, each contributing all of its genes, in index order.
    mixed=True interleaves genes from the (exactly two) ALGs instead of
    placing them as contiguous blocks."""
    rows = []
    pos = SCAF_START
    if mixed:
        assert len(alg_blocks) == 2
        a, b = alg_blocks
        na, nb = ALG_GENE_COUNT[a], ALG_GENE_COUNT[b]
        assert na == nb
        for i in range(na):
            for alg in (a, b):
                rows.append((alg, i, scaf_name, pos))
                pos += SCAF_GENE_SPACING
    else:
        for alg in alg_blocks:
            for g in range(ALG_GENE_COUNT[alg]):
                rows.append((alg, g, scaf_name, pos))
                pos += SCAF_GENE_SPACING
    return rows


def build_taxon_layouts():
    """Returns dict taxon_name -> list of (alg, gene_idx, scaffold, pos)."""
    layouts = {}

    # T1: outgroup, everything unfused, + one tiny ALG-N scaffold
    rows = unfused_layout("T1", CORE_ALGS)
    rows.append(("N", 0, "T1_Nfrag", TINY_SCAF_POS))
    layouts["T1"] = rows

    # T2: G+H fused, clean
    rows = unfused_layout("T2", [a for a in CORE_ALGS if a not in ("G", "H")])
    rows += fused_layout("T2", "T2_GH", ["G", "H"])
    layouts["T2"] = rows

    # T3, T4: A+B fused clean, E+F fused clean (independent scaffolds per taxon)
    for taxon in ("T3", "T4"):
        rest = [a for a in CORE_ALGS if a not in ("A", "B", "E", "F")]
        rows = unfused_layout(taxon, rest)
        rows += fused_layout(taxon, f"{taxon}_AB", ["A", "B"])
        rows += fused_layout(taxon, f"{taxon}_EF", ["E", "F"])
        layouts[taxon] = rows

    # T5: A+B, C+D, G+H all fused clean
    rest = [a for a in CORE_ALGS if a not in ("A", "B", "C", "D", "G", "H")]
    rows = unfused_layout("T5", rest)
    rows += fused_layout("T5", "T5_AB", ["A", "B"])
    rows += fused_layout("T5", "T5_CD", ["C", "D"])
    rows += fused_layout("T5", "T5_GH", ["G", "H"])
    layouts["T5"] = rows

    # T6: A+B, C+D fused clean, + the ALG-N ortholog (tiny scaffold, matches T1)
    rest = [a for a in CORE_ALGS if a not in ("A", "B", "C", "D")]
    rows = unfused_layout("T6", rest)
    rows += fused_layout("T6", "T6_AB", ["A", "B"])
    rows += fused_layout("T6", "T6_CD", ["C", "D"])
    rows.append(("N", 0, "T6_Nfrag", TINY_SCAF_POS))
    layouts["T6"] = rows

    # T7_wgd: every core ALG present as two full, independent, unfused copies
    rows = []
    for alg in CORE_ALGS:
        for copy in ("cpA", "cpB"):
            scaf = f"T7_wgd_{alg}_{copy}"
            pos = SCAF_START
            for g in range(ALG_GENE_COUNT[alg]):
                rows.append((alg, g, scaf, pos))
                pos += SCAF_GENE_SPACING
    layouts["T7_wgd"] = rows

    # T8_mixed: core ALGs unfused (baseline), plus I+J fused AND mixed
    rows = unfused_layout("T8_mixed", CORE_ALGS, start_scaf_suffix="_b")
    rows += fused_layout("T8_mixed", "T8_mixed_IJ", ["I", "J"], mixed=True)
    layouts["T8_mixed"] = rows

    return layouts


# ==========================================================================
# Sequence generation + file writing
# ==========================================================================

def build_ancestral_sequences(rng):
    master = {}
    for alg, n in ALG_GENE_COUNT.items():
        for g in range(n):
            master[(alg, g)] = random_sequence(rng)
    return master


def write_taxon(taxon, layout, master_seqs, mutation_rate, rng, outdir):
    cds_path = outdir / f"{taxon}.cds"
    pep_path = outdir / f"{taxon}.pep"
    gff_path = outdir / f"{taxon}.gff3"
    with open(cds_path, "w") as fcds, open(pep_path, "w") as fpep, open(gff_path, "w") as fgff:
        fgff.write("##gff-version 3\n")
        for alg, gene_idx, scaf, pos in layout:
            gid = f"{taxon}_{alg}{gene_idx}"
            master = master_seqs[(alg, gene_idx)]
            protein = mutate(master, mutation_rate, rng)
            cds = reverse_translate(protein)
            fcds.write(f">{gid}\n{cds}\n")
            fpep.write(f">{gid}\n{protein}\n")
            end = pos + len(cds) - 1
            fgff.write(f"{scaf}\ttoy\tgene\t{pos}\t{end}\t.\t+\t.\tID={gid};Name={gid}\n")
    return cds_path, pep_path, gff_path


# ==========================================================================
# Diamond-compatibility self-check (no diamond required)
# ==========================================================================

def verify_sequence_separation(master_seqs, mutation_rate, rng, n_samples=200):
    """Confirms, without needing diamond installed, that the mutation rate
    gives real orthologs (independently-mutated copies of the same
    ancestral sequence) much higher identity than unrelated genes
    (independently-generated ancestral sequences) - i.e. that a real
    reciprocal-best-hit search run on this data will behave as intended."""
    keys = list(master_seqs.keys())

    ortholog_identities = []
    for _ in range(n_samples):
        k = rng.choice(keys)
        a = mutate(master_seqs[k], mutation_rate, rng)
        b = mutate(master_seqs[k], mutation_rate, rng)
        ortholog_identities.append(hamming_identity(a, b))

    unrelated_identities = []
    for _ in range(n_samples):
        k1, k2 = rng.sample(keys, 2)
        a = mutate(master_seqs[k1], mutation_rate, rng)
        b = mutate(master_seqs[k2], mutation_rate, rng)
        unrelated_identities.append(hamming_identity(a, b))

    def stats(xs):
        return min(xs), sum(xs) / len(xs), max(xs)

    o_min, o_mean, o_max = stats(ortholog_identities)
    u_min, u_mean, u_max = stats(unrelated_identities)

    print(f"[verify] true-ortholog pair identity:  min={o_min:.2f}  mean={o_mean:.2f}  max={o_max:.2f}")
    print(f"[verify] unrelated-gene pair identity: min={u_min:.2f}  mean={u_mean:.2f}  max={u_max:.2f}")

    ok = o_min > u_max and o_min > 0.5 and u_max < 0.5
    if ok:
        print(f"[verify] OK: every sampled ortholog pair (>= {o_min:.2f} identity) is clearly more "
              f"similar than every sampled unrelated pair (<= {u_max:.2f} identity) - real diamond "
              f"reciprocal-best-hit search should recover this dataset's designed orthology correctly.")
    else:
        print("[verify] WARNING: ortholog and unrelated-gene identity distributions overlap - "
              "lower --mutation-rate (less divergence) and/or increase PROTEIN_LEN before trusting "
              "this dataset's real-diamond behaviour.")
    return ok


# ==========================================================================
# Manifest
# ==========================================================================

MANIFEST_TEMPLATE = """\
# Toy dataset ground truth

Generated by make_toy_dataset.py (seed={seed}, mutation rate={mutation_rate}).

## Linkage groups

- **A-H** (the "core" set, {n_core} genes each): the fusion/mixing history below.
- **I, J** ({n_mixed} genes each): used only by T8_mixed, to demonstrate the
  "fused, mixed" (state 2) character call in isolation.
- **N** (1 gene): a single ortholog pair (T1 <-> T6) each sitting alone on
  its own ~300 bp scaffold, to demonstrate --min-scaffold-length.

## Taxa and designed states

| Taxon | A+B | C+D | E+F | G+H | I+J | Notes |
|---|---|---|---|---|---|---|
| T1 | 0 | 0 | 0 | 0 | - | outgroup; carries ALG N on a tiny scaffold |
| T2 | 0 | 0 | 0 | 1 (clean) | - | |
| T3 | 1 (clean) | 0 | 1 (clean) | 0 | - | clade partner of T4 |
| T4 | 1 (clean) | 0 | 1 (clean) | 0 | - | clade partner of T3 |
| T5 | 1 (clean) | 1 (clean) | 0 | 1 (clean) | - | clade partner of T6 |
| T6 | 1 (clean) | 1 (clean) | 0 | 0 | - | clade partner of T5; carries ALG N |
| T7_wgd | 0 | 0 | 0 | 0 | - | every core ALG duplicated (2 full copies each) |
| T8_mixed | 0 | 0 | 0 | 0 | 2 (mixed) | isolated mixing-state demonstration |

States: 0 = unfused, 1 = fused-clean, 2 = fused-mixed, "-" = ALG not present in this taxon.

## True tree for T1-T6

    (T1,(T2,((T3,T4),(T5,T6))))

Three clean, nested synapomorphies (A+B unites T3-T6; E+F unites T3+T4;
C+D unites T5+T6) plus one deliberately homoplastic character (G+H, arising
independently in T2 and T5) - the same design used to validate CI/RI,
bootstrap support and the permutation test during this project's own
development. Exact CI/RI values will depend on which other taxa (T7_wgd,
T8_mixed) are included in a given run, since they contribute additional
(uninformative, all-0) data points to the A-H characters' m/s/g counts;
run with just T1-T6 to reproduce the original validated numbers exactly
(ensemble CI = 0.800, RI = 0.750, per-character CI/RI in the manuscript).

## What each taxon is *for*

- **T1-T6**: core parsimony validation - character calling, tree search
  (reversible and irreversible), synapomorphy mapping, CI/RI, bootstrap,
  permutation test, backbone constraints.
- **T7_wgd**: --detect-wgd. Without it, plain reciprocal-best-hit should
  under-recover T7's orthologs for every core ALG (each duplicated ALG's
  hits split across its two scaffold copies); with --detect-wgd
  --wgd-species T7_wgd, both copies should be fully recovered and
  duplication_blocks.tsv should report 8 confirmed ohnologous scaffold
  pairs (one per core ALG).
- **T8_mixed**: mixing-state detection (state 2) and the ideogram/dot-plot
  visual signature of interleaving, in isolation from the rest of the
  dataset's signal.
- **T1 + T6's ALG-N genes**: --min-scaffold-length. Both sit on ~300 bp
  scaffolds; a threshold anywhere above that (e.g. 1000) should exclude
  them, logged as "excluded 1 gene(s) on 1 scaffold(s)" for each taxon.

## Example commands

```bash
# single-reference mode, T1 as reference, core validation set only
python3 genomes_to_rbh.py --input-dir {outdir} --reference-species T1 \\
    --species T1 T2 T3 T4 T5 T6 \\
    --rbh-outdir toy_rbh --parsimony-outdir toy_parsimony \\
    --forward --min-genes 3 --alpha 0.05 --outgroup T1 \\
    --plot-dotplots --plot-trees --bootstrap 100 --permutation-test 100

# consensus mode: T1 (outgroup) keeps the reconstruction honest even
# though T5 privately disagrees about C+D/G+H
python3 genomes_to_rbh.py --input-dir {outdir} \\
    --consensus-species T1 T3 T5 --species T1 T2 T3 T4 T5 T6 \\
    --rbh-outdir toy_rbh_consensus --parsimony-outdir toy_parsimony_consensus \\
    --forward --min-genes 3 --alpha 0.05 --outgroup T1

# whole-genome duplication rescue
python3 genomes_to_rbh.py --input-dir {outdir} --reference-species T1 \\
    --species T1 T7_wgd --detect-wgd --wgd-species T7_wgd \\
    --rbh-outdir toy_rbh_wgd --rbh-only
cat toy_rbh_wgd/duplication_blocks.tsv

# minimum scaffold length filter (excludes ALG N's tiny scaffolds)
python3 genomes_to_rbh.py --input-dir {outdir} --reference-species T1 \\
    --species T1 T6 --min-scaffold-length 1000 \\
    --rbh-outdir toy_rbh_minscaf --rbh-only

# mixed-fusion / ideogram demonstration
python3 genomes_to_rbh.py --input-dir {outdir} --reference-species T1 \\
    --species T1 T8_mixed \\
    --rbh-outdir toy_rbh_mixed --parsimony-outdir toy_parsimony_mixed \\
    --forward --min-genes 3 --alpha 0.05 --plot-dotplots --plot-ideograms
```
"""


def write_manifest(outdir, seed, mutation_rate):
    text = MANIFEST_TEMPLATE.format(
        seed=seed, mutation_rate=mutation_rate, n_core=N_CORE_GENES,
        n_mixed=N_MIXED_GENES, outdir=outdir,
    )
    path = outdir / "MANIFEST.md"
    with open(path, "w") as f:
        f.write(text)
    return path


# ==========================================================================
# Main
# ==========================================================================

def main():
    ap = argparse.ArgumentParser(
        description="Generate a small, fully-synthetic toy dataset (real, diamond-alignable "
                    "sequences) with known ground truth for testing genomes_to_rbh.py / "
                    "alg_parsimony.py.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument("--outdir", default="toy_data", help="output directory (default: toy_data)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--mutation-rate", type=float, default=0.06,
                     help="per-residue substitution probability applied independently to every "
                          "taxon's copy of every gene (default: 0.06)")
    ap.add_argument("--verify-only", action="store_true",
                     help="don't write any files - just run the sequence-separation sanity check "
                          "for the given --seed/--mutation-rate and report the result")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    master_seqs = build_ancestral_sequences(rng)

    ok = verify_sequence_separation(master_seqs, args.mutation_rate, rng)
    if args.verify_only:
        sys.exit(0 if ok else 1)
    if not ok:
        print("[warn] proceeding to write files despite the warning above - inspect the identity "
              "values before relying on a real diamond run against this data.")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    layouts = build_taxon_layouts()
    for taxon, layout in layouts.items():
        cds_path, pep_path, gff_path = write_taxon(taxon, layout, master_seqs, args.mutation_rate,
                                                     rng, outdir)
        print(f"[info] wrote {taxon}: {len(layout)} genes -> {cds_path.name}, {pep_path.name}, "
              f"{gff_path.name}")

    manifest_path = write_manifest(outdir, args.seed, args.mutation_rate)
    print(f"[info] wrote {manifest_path}")
    print(f"[info] {len(layouts)} taxa written to {outdir}/ - see MANIFEST.md for ground truth "
          f"and example commands")


if __name__ == "__main__":
    main()
