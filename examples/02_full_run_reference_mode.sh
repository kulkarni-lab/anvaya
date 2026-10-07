#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 2. Full run from raw genomes, single-reference mode
# ---------------------------------------------------------------------------
# Linkage groups are read off ONE reference genome. Simple and fast, but the
# results inherit that genome's own rearrangements and assembly errors, so
# describe every fusion as relative to it.
#
# The reference genome must NOT carry a whole-genome duplication the rest of
# your sample lacks - see the README's first section.
#
# In single-reference mode the reference is not a taxon in the tree, so if you
# want it rooted on the reference you must also pass it on the comparison side.
# ---------------------------------------------------------------------------
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
