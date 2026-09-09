# Anvaya

A toolkit for turning chromosome-scale synteny data into a parsimony-based
phylogenetic hypothesis, using ancestral-linkage-group (ALG) **fusion** and
post-fusion **mixing** as discrete characters.

| Script | Purpose |
|---|---|
| **`genomes_to_rbh.py`** | Raw genomes → `.rbh` ortholog tables. Builds reciprocal-best-hit orthology with DIAMOND, defines linkage groups from either one reference genome or a multi-genome consensus, optionally filters short scaffolds and rescues whole-genome-duplication ohnologs. |
| **`alg_parsimony.py`** | `.rbh` tables → character matrix, parsimony trees, synapomorphies, support statistics and publication figures. |
| **`make_toy_dataset.py`** | Generates a small synthetic dataset with known ground truth, for testing the pipeline end to end (including against a real DIAMOND install). |
| **`sweep_consensus_params.py`** | Sweeps the three consensus-mode thresholds and writes one TSV summarising the linkage-group counts, to guide parameter choice. |

**A full tutorial with eight worked examples is provided as
`Anvaya_tutorial.pdf`** — covering duplicated vs. unduplicated genomes, both
operating modes, threshold tuning, figure generation and troubleshooting, with
runnable scripts for cases across animals and plants.

Only `pandas` is a hard requirement. `matplotlib` (figures), `scipy`
(enrichment tests), `biopython` (CDS→protein translation) and `diamond`
(orthology search) are used when available and degrade gracefully — with a
warning, not a crash — when they are not.

---

## Before anything else: which of your genomes are duplicated?

Anvaya defines linkage groups from either one **reference** genome or a set of
**core** genomes. Those genomes must not carry a whole-genome duplication (WGD)
that the rest of your sample does not share.

Reciprocal best hit is a strict 1:1 mapping, so a duplicated genome breaks it:

- **As a reference** — each ancestral chromosome is present twice and gets
  registered as *two* linkage groups, making every non-duplicated taxon look
  like it is missing half of each.
- **As a core genome** — orthogroups containing two genes from the same core
  genome are discarded as ambiguous, so its ohnologous orthogroups vanish
  silently and gut the reconstruction.

> **Rule:** duplicated lineages belong on the *comparison* side, never in
> `--reference-species` or `--consensus-species`. Flag them with
> `--detect-wgd --wgd-species …` so their dropped ohnologs are recovered.

**"Duplicated" is relative to your sample.** All vertebrates share the ancient
2R rounds, so among vertebrates 2R is shared ancestry and causes no asymmetry.
Add an invertebrate outgroup and 2R becomes lineage-specific *relative to it*.
Always ask: duplicated relative to what?

| Lineage | Status |
|---|---|
| Nematodes, insects, ticks/mites, sea spiders | Unduplicated — good core/reference taxa |
| Amphioxus, spotted gar, *Xenopus tropicalis*, *Amborella* | Unduplicated relative to their usual comparators — ideal references |
| Most molluscs and annelids | Usually unduplicated — verify per species |
| Teleost fish (3R) | **Duplicated** — comparison side |
| Salmonids (3R + 4R) | **Duplicated**, multiple rounds — partial rescue only |
| *Xenopus laevis* (allotetraploid) | **Duplicated** — use *X. tropicalis* as reference |
| Spiders, scorpions, whip scorpions (Arachnopulmonata) | **Duplicated** — comparison side |
| Horseshoe crabs (Xiphosura) | **Duplicated**, multiple rounds — partial rescue only |
| Most crop plants (wheat, soy, *Brassica*) | **Duplicated** — verify per species |

`--detect-wgd` was designed and validated for a **single, clean, pairwise**
duplication. Lineages with multiple rounds, or heavily rearranged post-WGD
genomes, get partial rescue — a substantial improvement on the plain-RBH false
negative, but treat it as a lower bound and inspect `duplication_blocks.tsv`.

See §2 of the tutorial PDF for the full discussion.

## Which script do I start with?

```
Do you already have .rbh files built against a real ALG reference
(BCnSSimakov2022_*, odp_LGs_from_*, or similar)?
                    |
        ┌───────────┴───────────┐
       YES                      NO
        |                        |
alg_parsimony.py          genomes_to_rbh.py first
directly                  (builds .rbh from raw genomes,
                           then calls alg_parsimony.py)
```

---

## Flag reference

Every flag, one line each. Defaults in brackets.

### `genomes_to_rbh.py`

**Input / output**

| Flag | Meaning |
|---|---|
| `--input-dir` | Directory of genomes: `{name}.cds` + `{name}.gff3\|.gtf` (+ optional `{name}.pep`) per species. |
| `--rbh-outdir` | Where to write the `.rbh` tables (required). |
| `--species` | Restrict the run to these species, plus the reference/core set [all discovered]. |
| `--work-dir` | Scratch directory for DIAMOND databases and intermediates [inside `--rbh-outdir`]. |
| `--rbh-only` | Stop after writing `.rbh` files; don't run the parsimony step. |
| `--parsimony-script` | Path to `alg_parsimony.py` [auto-detected next to this script, in `$PWD`, or on `$PATH`]. |
| `--parsimony-outdir` | Where `alg_parsimony.py` should write its output. |
| `--forward` | Everything after this is passed straight through to `alg_parsimony.py`. |

**Linkage-group definition** (pick exactly one mode)

| Flag | Meaning |
|---|---|
| `--reference-species` | Single-reference mode: this genome's own chromosomes *are* the linkage groups. Fast; not itself a taxon in the tree. |
| `--consensus-species` | Consensus mode: reconstruct linkage groups de novo from these ≥2 "core" genomes. All genomes, core included, become taxa. |
| `--consensus-min-species` | An orthogroup must be found in at least this many core genomes to be used [all core species]. |
| `--consensus-min-frac` | Fraction of core genomes that must agree two orthogroups are co-located before merging them into one linkage group [1.0 = unanimous]. |
| `--consensus-min-alg-orthogroups` | Discard any reconstructed linkage group with fewer than this many member orthogroups [1 = keep everything]. |

**Filtering and duplication**

| Flag | Meaning |
|---|---|
| `--min-scaffold-length` | Exclude genes on scaffolds shorter than this many bp from the orthology search entirely [0 = off]. |
| `--detect-wgd` | Detect retained whole-genome/segmental duplication and rescue ohnologs that plain reciprocal-best-hit silently drops. |
| `--wgd-species` | Restrict `--detect-wgd` to these species [all non-reference/non-core species]. |
| `--wgd-min-paralog-pairs` | Minimum intragenomic paralog pairs to call two scaffolds ohnologous [5]. |
| `--wgd-alpha` | Significance cutoff for the ohnologous-scaffold enrichment test [0.05]. |
| `--wgd-output-dir` | Reuse already-translated `.faa` files from a prior `run_wgd_pipeline.py` run. |

**Search parameters**

| Flag | Meaning |
|---|---|
| `--evalue` | DIAMOND E-value cutoff [1e-3]. |
| `--max-target-seqs` | DIAMOND hits retained per query [5]. |
| `--threads` | DIAMOND threads [4]. |
| `--genetic-code` | NCBI translation table for CDS→protein [1]. |
| `--no-reuse-blast` | Recompute alignments even if a cached result with the same `--evalue`/`--max-target-seqs` exists in the work directory [off; cached results are reused]. |

### `alg_parsimony.py`

**Input / output**

| Flag | Meaning |
|---|---|
| `--indir` | Directory of `.rbh` files, one per taxon (required). |
| `--outdir` | Output directory [`alg_parsimony_out`]. |
| `--reference-prefix` | Override auto-detection of the shared ALG/reference column prefix. |

**Character calling**

| Flag | Meaning |
|---|---|
| `--alpha` | `whole_FET` significance cutoff for calling an ALG present on a scaffold [0.05]. |
| `--min-genes` | Minimum shared orthologs to call an ALG present [3]. |
| `--mixing-threshold` | Fraction of gene-order switches along a fused scaffold above which the fusion is "mixed" (state 2) rather than "clean" (state 1) [0.3]. |

**Tree search**

| Flag | Meaning |
|---|---|
| `--outgroup` | Taxon to root the irreversible analysis on [auto: fewest derived states — verify this choice]. |
| `--constraint` | Backbone constraint as a Newick string or file; may name any subset of taxa, polytomies allowed. |
| `--max-exhaustive` | Taxon-count ceiling for exact tree enumeration; above this a heuristic NNI search runs [8]. |
| `--nni-restarts` | Random restarts for the heuristic search [10]. |
| `--seed` | RNG seed, for reproducible heuristic searches and replicates [1]. |

**Support and diagnostics**

| Flag | Meaning |
|---|---|
| `--bootstrap N` | Nonparametric bootstrap replicates; support is written both into the `.nwk` file and a `.tsv` [0 = off]. |
| `--bootstrap-restarts` | NNI restarts per bootstrap replicate, heuristic regime only [3]. |
| `--permutation-test N` | Signal-test replicates: does the matrix beat randomly shuffled data at all? [0 = off]. |
| `--permutation-restarts` | NNI restarts per permutation replicate, heuristic regime only [3]. |

CI/RI homoplasy indices are always computed — no flag needed.

**Figures**

| Flag | Meaning |
|---|---|
| `--plot-dotplots` | Oxford-style synteny dot plot per taxon. |
| `--plot-ideograms` | Chromosome-painting bars per taxon — the collapsed version of the dot plot. |
| `--plot-ribbon` | Multi-track ribbon/riparian plot connecting orthologs between adjacent taxa. |
| `--ribbon-order` | Explicit top-to-bottom taxon order for the ribbon plot; any subset of ≥2 taxa. |
| `--ribbon-order-from-tree` | Order ribbon tracks by a ladderized traversal of the `reversible` or `irreversible` tree, so clades stay contiguous. |
| `--plot-matrix` | Taxa × characters heatmap of the character matrix. |
| `--plot-trees` | Rendered cladogram per criterion, with synapomorphies labelled on branches. |
| `--plot-format` | `pdf`, `png` or `svg` for every figure [pdf]. |
| `--highlight-wgd` | Extra dot plot + ribbon plot + summary table highlighting WGD-rescued orthologs. |

**Interoperability**

| Flag | Meaning |
|---|---|
| `--export-nexus` | Write a NEXUS matrix with a PAUP* command block. |
| `--export-tnt` | Write a Hennig86 matrix with a TNT script (reversible criterion only). |
| `--run-paup` / `--run-tnt` | Also execute PAUP*/TNT if present on the system. |
| `--paup-bin` / `--tnt-bin` | Executable name or path [`paup`, `tnt`]. |

---

## Choosing the consensus-mode thresholds

The three `--consensus-*` thresholds are the parameters most likely to
decide whether your reconstruction is usable, and none of them has a
universally correct value. They do genuinely different jobs, so tune them
in this order.

### The quick way: sweep all three at once

`sweep_consensus_params.py` runs the sweep for you and writes a single TSV:

```bash
python3 sweep_consensus_params.py \
    --input-dir ./Data \
    --consensus-species Ixod Varr Sarc Neoseiu \
    --outdir ./sweep \
    --min-species 3 4 \
    --min-frac 1.0 0.75 0.5 \
    --min-alg-orthogroups 1 5 10 15 20 30 \
    --threads 30
```

```
min_species  min_frac  min_alg_orthogroups  n_linkage_groups  n_orthogroups_kept  ...  largest_group  pct_orthogroups_retained
3            1.0       1                    9                 60                       14             100.0
3            1.0       15                   0                 60                       0              0.0
3            0.75      10                   4                 60                       15             100.0
```

Reading it:

- **`n_linkage_groups`** — the number you are trying to make biologically
  plausible; compare against known haploid chromosome counts for your clade.
- **`pct_orthogroups_retained`** — how much of your data ends up inside a kept
  group. A combination with a nice group count but tiny retention is a bad
  trade, and the row with `min_alg_orthogroups 15` above shows how a size
  filter set too high can silently discard everything.
- **`largest_group`** — if one group swallows most of the orthogroups,
  `--consensus-min-frac` is too permissive and unrelated groups are merging.

It is much faster than it looks. `--consensus-min-alg-orthogroups` is a pure
post-reconstruction size filter, so every value on that axis is derived
arithmetically from one reconstruction instead of rerunning it; and all
combinations share one work directory, so the DIAMOND alignments are computed
once for the whole sweep. Only `min-species × min-frac` actually triggers a
run.

Use `--extra` to pass anything else through, e.g.
`--extra --min-scaffold-length 100000`.

The rest of this section explains what each threshold does, for interpreting
that table.

### Step 0 — a diagnostic run with everything permissive

```bash
python3 genomes_to_rbh.py --input-dir ./Data \
    --consensus-species SpA SpB SpC SpD \
    --rbh-outdir ./rbh_diag --rbh-only --threads 30 \
    2>&1 | grep -E "orthogroups|reconstructed"
```

Two log lines matter:

```
built 7314 orthogroups from 10844 RBH-graph components (3530 dropped as ambiguous/multi-copy...)
2074/7314 orthogroups present in >= 6 core species        <- the presence filter
reconstructed 1896 consensus linkage groups (1896 before size filtering)
```

### `--consensus-min-species` — how many core genomes must contain a gene

Default is *all* core species. With a handful of closely related core
genomes that's fine. With many core genomes spanning deep divergences, it
gets strict fast: gene loss, elevated sequence divergence and annotation
gaps all compound, and a single miss in one lineage discards the orthogroup
entirely.

Watch the `N/M orthogroups present in >= K core species` line. If you are
keeping only a small fraction of orthogroups, relax it — allowing one or two
lineages to be missing (e.g. `4` with 6 core species) often recovers a large
amount of usable signal at little cost.

### `--consensus-min-frac` — how much disagreement to tolerate before merging

Default is `1.0`: every core genome that has both orthogroups must agree
they sit on the same chromosome. In clades with real rearrangement, one
lineage's private translocation is then enough to block a genuine ancestral
linkage group from being reconstructed at all — the reconstruction shatters
into many small groups.

Sweep it and watch how the *raw* group count (before size filtering) moves:

```bash
for frac in 1.0 0.8 0.7 0.6; do
  echo "--- min-frac $frac ---"
  python3 genomes_to_rbh.py --input-dir ./Data \
      --consensus-species SpA SpB SpC SpD \
      --consensus-min-species 3 --consensus-min-frac $frac \
      --rbh-outdir ./rbh_frac_$frac --rbh-only --threads 30 \
      2>&1 | grep reconstructed
done
```

You want a **smooth, monotonic decrease**. If the count collapses abruptly
toward 1, that means unrelated groups are being merged — back off to a
stricter value. (Merging uses agglomerative clustering with cluster-level
re-validation at every step, specifically so that weak pairwise links cannot
chain transitively into one giant group; a sudden collapse is nonetheless
worth treating as a warning that the threshold is too permissive for your
data.)

#### Watch out: only the fractions *k/n* are achievable

Agreement is measured as (genomes agreeing) / (genomes sharing both
orthogroups), so a threshold sitting just **above** one of those discrete
fractions silently behaves like the next one up. The classic trap:

> With 3 core genomes, `--consensus-min-frac 0.67` **rejects** 2-of-3
> agreement, because 2/3 = 0.6667 < 0.67 — making it indistinguishable from
> requiring unanimity. For "2 of 3 may disagree", use `0.66` or lower.

Every consensus run now prints what the threshold actually means, and warns
when the value sits just above a natural boundary:

```
--consensus-min-frac 0.67 means, for orthogroup pairs shared by 2/3/4 core
genomes respectively: 2/2, 3/3, 3/4 must agree
NOTE: --consensus-min-frac 0.67 is just above 2/3 (0.6667); if you intended
'2 of 3 genomes may agree', use 0.662 or lower - as given, 2-of-3 agreement
is REJECTED
```

If two different `--consensus-min-frac` values give byte-identical results,
this is almost certainly why.

#### Non-independent core taxa change what a threshold means

`--consensus-min-frac` counts genomes, not lineages. If two of your four
core genomes are sister taxa, then at `0.5` (2-of-4) that single lineage can
carry a merge over the objections of both other lineages. Set the threshold
high enough that at least one genome from outside any one lineage must
concur — with four core genomes, two of which are sisters, that means
`0.75`.

### `--consensus-min-alg-orthogroups` — discard small groups outright

The other two control *why* groups form. This one simply drops whatever
small groups remain, regardless of cause. Every run writes
`consensus_linkage_group_sizes.tsv`, ranking each candidate group by member
count with a `kept` column, so you can pick the cut-off from the data:

```bash
sort -k2,2 -rn rbh_diag/consensus_linkage_group_sizes.tsv | column -t | head -40
```

**If there is a clear cliff** — a cluster of large groups, then a sharp drop
to a tail of small ones — set the threshold just above the drop.

**If the distribution decays smoothly**, which is common for
phylogenetically broad samples, then no threshold is objectively correct.
Don't keep hunting for a knee that isn't there. Instead:

1. Pick a value giving a biologically plausible number of linkage groups
   (compare against known haploid chromosome counts for your clade).
2. State it in your write-up as a **chosen analytical cut-off**, not a
   data-derived one.
3. **Run the downstream analysis at several thresholds** and report whether
   your conclusions — topology, supported clades, inferred events — are
   stable across them.

If the conclusions are stable, the exact cut-off doesn't matter. If they
aren't, that instability *is* the finding, and belongs in the paper rather
than hidden behind a single flattering threshold. Because
`--consensus-min-frac` and `--consensus-min-alg-orthogroups` interact, vary
both in that sensitivity analysis — not just one.

### Choosing the core set itself

- **Never put a whole-genome-duplicated lineage in the core set.** Its
  ohnologs are dropped as ambiguous multi-copy components, silently gutting
  the reconstruction. Put duplicated genomes on the comparison side and use
  `--detect-wgd --wgd-species …` instead.
- **Prefer phylogenetically independent core lineages.** Several close
  relatives will agree with each other for reasons of recent shared
  ancestry, which makes the agreement check circular rather than
  informative.
- **Check assembly contiguity before committing.** A fragmented core genome
  inflates the linkage-group count no matter how the thresholds are set:

  ```bash
  for sp in SpA SpB SpC SpD; do
    echo "=== $sp ==="
    f=$(ls Data/${sp}.gff3 Data/${sp}.gtf 2>/dev/null | head -1)
    awk -F'\t' '!/^#/{if($5>max[$1]) max[$1]=$5} END{for (s in max) print max[s], s}' "$f" \
      | sort -rn | head -5
  done
  ```

  (Note the `.gff3`-or-`.gtf` handling — assuming one extension silently
  produces empty output for species using the other.)

---

## Detecting whole-genome duplication

### Which genomes go where

Duplication-carrying lineages must **never** be core/reference genomes. In
consensus mode their ohnologs are discarded as ambiguous multi-copy
components, silently gutting the reconstruction; in single-reference mode
each ancestral chromosome would define two "linkage groups", one per
retained copy. Put them on the comparison side:

```
--consensus-species <non-duplicated taxa> --detect-wgd --wgd-species <duplicated taxa>
```

### Calibrating `--wgd-min-paralog-pairs`

Same two-step approach as the consensus thresholds. Run with the count
filter effectively off so the report lists every candidate:

```bash
python3 genomes_to_rbh.py --input-dir ./Data \
    --consensus-species SpA SpB SpC \
    --detect-wgd --wgd-species DupTaxon \
    --wgd-min-paralog-pairs 1 \
    --rbh-outdir ./rbh_wgd_diag --rbh-only --threads 30

# NB: duplication_blocks.tsv starts with a taxon column, so the count is $4
sort -k4,4 -rn rbh_wgd_diag/duplication_blocks.tsv | head -40
awk 'NR>1{print $4}' rbh_wgd_diag/duplication_blocks.tsv | sort -rn | uniq -c
```

Look for the break between real ohnologous chromosome pairs (typically
dozens to hundreds of supporting paralog pairs, with very small p-values)
and incidental paralogy (a large mass of pairs at 1–5). Set the threshold
just above that noise floor.

Sanity checks:

- **Confirmed pairs vs. chromosome number.** One clean WGD should give
  roughly one ohnologous partner per ancestral chromosome. Many more than
  that means the threshold is still too permissive.
- **Scaffolds with several partners.** After an ancient WGD, post-duplication
  fusion and fission mean a modern chromosome often carries material from
  more than one ancestral ohnolog, so a *strong* partner plus one or two
  weaker links is expected and interpretable. A dense, evenly-connected
  network is not, and suggests either a too-low threshold or a duplication
  history beyond the pairwise model this implementation assumes.

Then verify visually rather than trusting counts alone:

```bash
--forward ... --plot-dotplots --plot-ribbon --highlight-wgd
```

`dotplots/<taxon>.wgd.pdf` should show rescued orthologs (drawn at double
diameter and ringed) forming paired diagonals across ohnologous chromosomes.

## Example scripts

### 1. Quick start from existing `.rbh` files

```bash
python3 alg_parsimony.py \
    --indir ./rbh_files \
    --outdir ./results \
    --min-genes 3 --alpha 0.05 \
    --outgroup MyOutgroup \
    --plot-dotplots --plot-trees --plot-matrix
```

### 2. Full run from raw genomes, single-reference mode

```bash
#!/usr/bin/env bash
set -euo pipefail

python3 genomes_to_rbh.py \
    --input-dir ./Data \
    --reference-species Ixod \
    --min-scaffold-length 100000 \
    --threads 30 \
    --rbh-outdir ./rbh_out \
    --parsimony-script "$(pwd)/alg_parsimony.py" \
    --parsimony-outdir ./parsimony_out \
    --forward \
      --min-genes 3 --alpha 0.05 \
      --outgroup Ixod \
      --plot-dotplots --plot-trees --plot-ideograms \
      --bootstrap 100 --permutation-test 100
```

### 3. Consensus mode with WGD lineages, and a tuning sweep

```bash
#!/usr/bin/env bash
# Core set: non-duplicated, well-assembled, phylogenetically spread taxa.
# WGD lineages go on the comparison side with --detect-wgd, never in the core.
set -euo pipefail

CORE="SeaSpider Tick Mite Harvestman"
WGD="Spider Scorpion HorseshoeCrab"

# --- Step 1: diagnostic, everything permissive -----------------------------
python3 genomes_to_rbh.py \
    --input-dir ./Data --consensus-species $CORE \
    --rbh-outdir ./rbh_diag --rbh-only --threads 30 \
    2>&1 | grep -E "orthogroups|reconstructed"

echo "--- linkage-group size distribution ---"
sort -k2,2 -rn ./rbh_diag/consensus_linkage_group_sizes.tsv | column -t | head -40

# --- Step 2: sweep the thresholds that interact ----------------------------
for frac in 1.0 0.8 0.7; do
  for minalg in 10 20; do
    echo "=== min-frac=$frac min-alg-orthogroups=$minalg ==="
    python3 genomes_to_rbh.py \
        --input-dir ./Data --consensus-species $CORE \
        --consensus-min-species 3 \
        --consensus-min-frac "$frac" \
        --consensus-min-alg-orthogroups "$minalg" \
        --rbh-outdir "./rbh_f${frac}_m${minalg}" --rbh-only --threads 30 \
        2>&1 | grep reconstructed
  done
done

# --- Step 3: full analysis at a chosen setting -----------------------------
python3 genomes_to_rbh.py \
    --input-dir ./Data \
    --consensus-species $CORE \
    --consensus-min-species 3 \
    --consensus-min-frac 0.7 \
    --consensus-min-alg-orthogroups 15 \
    --detect-wgd --wgd-species $WGD --wgd-min-paralog-pairs 30 \
    --min-scaffold-length 100000 \
    --threads 30 \
    --rbh-outdir ./rbh_out \
    --parsimony-script "$(pwd)/alg_parsimony.py" \
    --parsimony-outdir ./parsimony_out \
    --forward \
      --min-genes 3 --alpha 0.05 \
      --outgroup SeaSpider \
      --constraint constraint.txt \
      --plot-dotplots --plot-trees --plot-matrix --plot-ideograms \
      --plot-ribbon --ribbon-order-from-tree reversible \
      --highlight-wgd \
      --bootstrap 100 --permutation-test 100 \
      --export-nexus --export-tnt
```

### 4. Sensitivity analysis on the downstream conclusions

```bash
#!/usr/bin/env bash
# Are the biological conclusions stable across plausible thresholds?
set -euo pipefail
CORE="SeaSpider Tick Mite Harvestman"

for minalg in 10 15 20 30; do
  python3 genomes_to_rbh.py \
      --input-dir ./Data --consensus-species $CORE \
      --consensus-min-species 3 --consensus-min-frac 0.7 \
      --consensus-min-alg-orthogroups "$minalg" \
      --threads 30 \
      --rbh-outdir "./rbh_t${minalg}" \
      --parsimony-script "$(pwd)/alg_parsimony.py" \
      --parsimony-outdir "./parsimony_t${minalg}" \
      --forward --min-genes 3 --alpha 0.05 \
        --outgroup SeaSpider --bootstrap 100
done

echo "=== recovered topologies across thresholds ==="
for minalg in 10 15 20 30; do
  printf "%-6s " "t=$minalg"
  cat "./parsimony_t${minalg}/tree_reversible.nwk"
done
```

### 5. Test the pipeline end to end on synthetic data

```bash
python3 make_toy_dataset.py --outdir toy_data      # see toy_data/MANIFEST.md
python3 make_toy_dataset.py --verify-only          # sanity-check sequences only

# Note: in single-reference mode the reference genome (T1) is NOT a taxon in
# the tree, so it cannot also be the --outgroup. T2 is the correct outgroup
# for this dataset once T1 is acting as the reference.
python3 genomes_to_rbh.py \
    --input-dir toy_data --reference-species T1 \
    --species T1 T2 T3 T4 T5 T6 \
    --rbh-outdir toy_rbh \
    --parsimony-script "$(pwd)/alg_parsimony.py" \
    --parsimony-outdir toy_parsimony \
    --forward --min-genes 3 --alpha 0.05 --outgroup T2 \
      --bootstrap 100 --permutation-test 100 --plot-trees

# Or, in consensus mode, T1 IS a taxon and can serve as the outgroup:
python3 genomes_to_rbh.py \
    --input-dir toy_data --consensus-species T1 T3 T5 \
    --species T1 T2 T3 T4 T5 T6 \
    --rbh-outdir toy_rbh_consensus \
    --parsimony-script "$(pwd)/alg_parsimony.py" \
    --parsimony-outdir toy_parsimony_consensus \
    --forward --min-genes 3 --alpha 0.05 --outgroup T1
```

---

## Output files

Written to `alg_parsimony.py`'s `--outdir` (figures in `--plot-format`,
PDF by default):

| File | Contents |
|---|---|
| `character_matrix.tsv` | Taxa × characters, states 0/1/2/`?`. |
| `alg_scaffold_calls.tsv` | Every ALG-pair fusion/mixing call per taxon. |
| `tree_reversible.nwk` / `.txt` | Best tree under ordinary (Fitch/Sankoff) parsimony; bootstrap support embedded as internal-node labels when `--bootstrap` is used. |
| `tree_irreversible.nwk` / `.txt` | Same, forbidding fusion/mixing reversal (Camin–Sokal). |
| `synapomorphies_*.txt` | Which fusion/mixing events unite which clades. |
| `ci_ri_*.tsv` | Per-character consistency and retention indices (always written). |
| `bootstrap_*.tsv` | Per-clade bootstrap support. |
| `permutation_test_*.tsv` | Per-replicate scores from the signal test. |
| `report.txt` | Everything above, human-readable. |
| `dotplots/`, `ideograms/`, `ribbon.*`, `character_matrix_heatmap.*`, `tree_*.pdf` | Figures. |
| `ohnolog_summary.tsv`, `ribbon_wgd/` | With `--highlight-wgd`: rescued-vs-direct ortholog counts and duplication-highlighted figures. |

From `genomes_to_rbh.py`'s `--rbh-outdir`:

| File | Contents |
|---|---|
| `<taxon>.rbh` | Ortholog table per taxon, including an `orthology_type` column (`direct` / `wgd_inferred`). |
| `consensus_linkage_group_sizes.tsv` | Consensus mode: every candidate linkage group ranked by member-orthogroup count, with a `kept` column. |
| `duplication_blocks.tsv` | With `--detect-wgd`: confirmed ohnologous scaffold pairs and their support. |

---

## Character model

For every pair of linkage groups co-occurring on a scaffold in some taxon:

- **0 — unfused**: the two groups sit on different scaffolds.
- **1 — fused, clean**: co-located, still two contiguous blocks.
- **2 — fused, mixed**: co-located and interleaved.

Two searches run over the same characters:

- **Reversible** — ordinary Fitch/Sankoff parsimony; changes may go either
  direction.
- **Irreversible** — states may only increase (0→1→2), modelling the claim
  that fusion and especially post-fusion mixing do not spontaneously revert.
  Requires a rooted tree via `--outgroup`.

Ancestral states are reconstructed for both, and shared derived states are
reported as synapomorphies on the clades they support.

---

## Caveats

- Fusion calls use whole-scaffold association only; sub-scaffold breakpoint
  (`break_FET`) signal is not modelled.
- A scaffold carrying 3+ linkage groups is decomposed into all pairwise
  characters, following convention — but those characters are not strictly
  independent.
- The heuristic search (above `--max-exhaustive` taxa) is an NNI
  hill-climber, not a validated general-purpose search. Cross-check with
  `--export-nexus` / `--export-tnt` in PAUP*/TNT for anything that matters.
  The same caveat applies to bootstrap and permutation replicates run in
  that regime.
- Single-reference mode inherits the reference genome's own rearrangements
  and assembly errors. Consensus mode is more robust but only as good as its
  core set. Either way, describe results as fusions *relative to* that
  reference or consensus set.
- `--detect-wgd` rescues ohnologs only on the compared side, and assumes
  duplication resolves into clean paralog **pairs**. Lineages with multiple
  rounds of duplication (e.g. horseshoe crabs) will be partially, not fully,
  rescued — treat those results as a lower bound and inspect
  `duplication_blocks.tsv` directly.
- Scaffold lengths come from GFF3 `##sequence-region` pragmas where present,
  otherwise from the furthest annotated feature end — a lower bound. GTF
  never carries such pragmas.

---

## Requirements

- **Python 3** with `pandas` (both scripts)
- **`diamond`** on `PATH` — required by `genomes_to_rbh.py`
- **`matplotlib`** — optional, for all `--plot-*` output
- **`scipy`** — optional, for Fisher's exact tests (falls back to count-only
  thresholds)
- **`biopython`** — optional, only when a species has no `.pep` and needs
  CDS→protein translation
- **`paup`** / **`tnt`** on `PATH` — optional, only for `--run-paup` /
  `--run-tnt`
