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
from collections import Counter

__version__ = "1.0.1b"
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

    # I and J are carried unfused by T1-T6 so that the mixed fusion in T8_mixed
    # can be called against T1 as reference and scored against other taxa.
    # (Before version 1.0 only T8_mixed carried them, so the mixed-fusion
    # demonstration had nothing to compare against.) Being unfused everywhere
    # in T1-T6, they add no character to the T1-T6 validation.
    IJ = ["I", "J"]

    # T1: outgroup, everything unfused, + one tiny ALG-N scaffold
    rows = unfused_layout("T1", CORE_ALGS + IJ)
    rows.append(("N", 0, "T1_Nfrag", TINY_SCAF_POS))
    layouts["T1"] = rows

    # T2: G+H fused, clean
    rows = unfused_layout("T2", [a for a in CORE_ALGS if a not in ("G", "H")] + IJ)
    rows += fused_layout("T2", "T2_GH", ["G", "H"])
    layouts["T2"] = rows

    # T3, T4: A+B fused clean, E+F fused clean (independent scaffolds per taxon)
    for taxon in ("T3", "T4"):
        rest = [a for a in CORE_ALGS if a not in ("A", "B", "E", "F")] + IJ
        rows = unfused_layout(taxon, rest)
        rows += fused_layout(taxon, f"{taxon}_AB", ["A", "B"])
        rows += fused_layout(taxon, f"{taxon}_EF", ["E", "F"])
        layouts[taxon] = rows

    # T5: A+B, C+D, G+H all fused clean
    rest = [a for a in CORE_ALGS if a not in ("A", "B", "C", "D", "G", "H")] + IJ
    rows = unfused_layout("T5", rest)
    rows += fused_layout("T5", "T5_AB", ["A", "B"])
    rows += fused_layout("T5", "T5_CD", ["C", "D"])
    rows += fused_layout("T5", "T5_GH", ["G", "H"])
    layouts["T5"] = rows

    # T6: A+B, C+D fused clean, + the ALG-N ortholog (tiny scaffold, matches T1)
    rest = [a for a in CORE_ALGS if a not in ("A", "B", "C", "D")] + IJ
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
        seen = Counter()
        for alg, gene_idx, scaf, pos in layout:
            # Every gene needs its own identifier. Before version 1.0 both
            # copies of a duplicated gene in T7_wgd received the same ID, so
            # the self-comparison could never pair them and --detect-wgd found
            # zero ohnologs instead of the eight designed.
            seen[(alg, gene_idx)] += 1
            copy = seen[(alg, gene_idx)]
            gid = f"{taxon}_{alg}{gene_idx}" + (f"_copy{copy}" if copy > 1 else "")
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

Generated by make_toy_dataset.py {version} (seed={seed}, mutation rate={mutation_rate}).
Every expectation below is checked automatically by `tests/validate_toy.sh`,
which runs the five validation commands and exits non-zero if any result
differs from what is stated here.

## Linkage groups

- **A-H** ({n_core} genes each): the fusion history below.
- **I, J** ({n_mixed} genes each): unfused in T1-T6; fused AND interleaved in
  T8_mixed, the one designed "fused-and-mixed" (state 2) character.
- **N** (1 gene): one ortholog pair (T1, T6), each copy alone on a ~300 bp
  scaffold, to demonstrate --min-scaffold-length.

## Taxa and designed states

| Taxon | A+B | C+D | E+F | G+H | I+J | Notes |
|---|---|---|---|---|---|---|
| T1 | 0 | 0 | 0 | 0 | 0 | outgroup / reference; carries ALG N on a tiny scaffold |
| T2 | 0 | 0 | 0 | 1 | 0 | |
| T3 | 1 | 0 | 1 | 0 | 0 | clade partner of T4 |
| T4 | 1 | 0 | 1 | 0 | 0 | clade partner of T3 |
| T5 | 1 | 1 | 0 | 1 | 0 | clade partner of T6 |
| T6 | 1 | 1 | 0 | 0 | 0 | clade partner of T5; carries ALG N |
| T7_wgd | 0 | 0 | 0 | 0 | - | every core ALG duplicated: two full copies on separate scaffolds |
| T8_mixed | 0 | 0 | 0 | 0 | 2 | I+J fused and interleaved |

States: 0 unfused, 1 fused-clean, 2 fused-and-mixed, "-" not present.

## True tree for T1-T6

    (T1,(T2,((T3,T4),(T5,T6))))

Three clean, nested synapomorphies (A+B unites T3-T6; E+F unites T3+T4; C+D
unites T5+T6) and one deliberately homoplastic character (G+H, gained
independently in T2 and T5).

## The five validation commands and their expected results

Run all five with `bash tests/validate_toy.sh` (needs DIAMOND on PATH).

**V1. Single-reference mode (T1 is the reference, so not a taxon).**
Taxa T2-T6; characters A+B, C+D, E+F, G+H; outgroup T2.
Expected: length 5 under both criteria. Irreversible: one most-parsimonious
tree, (T2,((T3,T4),(T5,T6))). Reversible: two, that tree and
(T2,(T5,(T6,(T3,T4)))), because with T2 as root the G+H fusion shared by T2
and T5 can be explained by one gain at the root and one reversal - exactly
what the irreversible criterion forbids.

**V2. Consensus mode (core genomes T1, T3, T5; T1 is a taxon).**
Taxa T1-T6; outgroup T1. Expected: one most-parsimonious tree under each
criterion, (T1,(T2,((T3,T4),(T5,T6)))); length 5 under both criteria;
ensemble CI 0.800 and RI 0.750; the G+H character alone has CI 0.5, RI 0.0.
T5 privately carries C+D and G+H fusions, but they are not shared by T1 and
T3, so the consensus keeps C, D, G and H as separate groups.

**V3. Whole-genome duplication rescue (T1 reference, T7_wgd compared).**
Expected: 48 intragenomic reciprocal-best-hit paralog pairs; 8 confirmed
ohnologous scaffold pairs (A_cpA-A_cpB ... H_cpA-H_cpB), each supported by 6
pairs, in duplication_blocks.tsv; 48 orthologs without --detect-wgd and 96
with it (48 flagged wgd_inferred).

**V4. Minimum scaffold length (T1 reference, T6 compared).**
Expected: with --min-scaffold-length 1000 the ALG-N gene is excluded in both
genomes, and the T6 table has one ortholog fewer than without the filter.

**V5. Mixed fusion (T1 reference; taxa T2-T6 and T8_mixed).**
Expected: the I+J character is coded 2 (fused-and-mixed) in T8_mixed and 0
in T2-T6; its mixing statistic (fraction of adjacent genes from different
linkage groups) is 1.0, against 1/11 for a clean two-block fusion.

## Commands

```bash
# V1 single-reference mode
python3 genomes_to_rbh.py --input-dir {outdir} --reference-species T1 \\
    --species T1 T2 T3 T4 T5 T6 \\
    --rbh-outdir v1_rbh --parsimony-outdir v1_out \\
    --forward --min-genes 3 --alpha 0.05 --outgroup T2 --seed 1

# V2 consensus mode
python3 genomes_to_rbh.py --input-dir {outdir} \\
    --consensus-species T1 T3 T5 --species T1 T2 T3 T4 T5 T6 \\
    --rbh-outdir v2_rbh --parsimony-outdir v2_out \\
    --forward --min-genes 3 --alpha 0.05 --outgroup T1 --seed 1 \\
    --bootstrap 100 --permutation-test 100

# V3 whole-genome duplication rescue (with and without)
python3 genomes_to_rbh.py --input-dir {outdir} --reference-species T1 \\
    --species T1 T7_wgd --rbh-outdir v3_rbh_plain --rbh-only
python3 genomes_to_rbh.py --input-dir {outdir} --reference-species T1 \\
    --species T1 T7_wgd --detect-wgd --wgd-species T7_wgd \\
    --rbh-outdir v3_rbh_wgd --rbh-only

# V4 minimum scaffold length (with and without)
python3 genomes_to_rbh.py --input-dir {outdir} --reference-species T1 \\
    --species T1 T6 --rbh-outdir v4_rbh_plain --rbh-only
python3 genomes_to_rbh.py --input-dir {outdir} --reference-species T1 \\
    --species T1 T6 --min-scaffold-length 1000 --rbh-outdir v4_rbh_minscaf --rbh-only

# V5 mixed fusion
python3 genomes_to_rbh.py --input-dir {outdir} --reference-species T1 \\
    --species T1 T2 T3 T4 T5 T6 T8_mixed \\
    --rbh-outdir v5_rbh --parsimony-outdir v5_out \\
    --forward --min-genes 3 --alpha 0.05 --outgroup T2 --seed 1 \\
    --plot-dotplots --plot-ideograms
```
"""


def write_manifest(outdir, seed, mutation_rate):
    text = MANIFEST_TEMPLATE.format(
        seed=seed, mutation_rate=mutation_rate, n_core=N_CORE_GENES,
        n_mixed=N_MIXED_GENES, outdir=outdir, version=__version__,
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
    ap.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
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
