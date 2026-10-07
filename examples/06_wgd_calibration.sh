#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 6. Calibrating --wgd-min-paralog-pairs, per taxon
# ---------------------------------------------------------------------------
# Paralog-pair counts scale with genome size and gene content, so one threshold
# across taxa is dominated by whichever genome is largest. Calibrate each taxon
# separately.
#
# Note what calibrate_wgd.py does NOT do: it does not threshold on statistical
# significance. Fisher's exact test asks whether a chromosome pair is enriched
# for paralogy against the genome-wide background, which is not the same
# question as whether the two chromosomes are ohnologous - two small
# chromosomes sharing three paralogs can be significantly enriched. Significance
# is necessary and nowhere near sufficient.
#
# The usable constraint is karyotype expectation: after ONE duplication every
# chromosome has exactly one ohnologous partner, so roughly n/2 pairs are
# expected, rising toward n with post-duplication rearrangement or further
# rounds. A count several times larger than n is noise, not ohnology.
# ---------------------------------------------------------------------------
set -euo pipefail

WGD="Spider Scorpion HorseshoeCrab"

# --- Step 1: diagnostic run with the count filter switched OFF -------------
# --wgd-min-paralog-pairs 1 reports every candidate chromosome pair, which is
# what the calibration needs to see.
python3 genomes_to_rbh.py \
    --input-dir ./Data \
    --consensus-species SeaSpider Tick Mite Harvestman \
    --consensus-min-species 3 --consensus-min-frac 0.7 \
    --consensus-min-alg-orthogroups 15 \
    --detect-wgd --wgd-species $WGD --wgd-min-paralog-pairs 1 \
    --rbh-outdir ./rbh_wgd_diag --rbh-only --threads 30

# --- Step 2: calibrate, and look at the figure before deciding -------------
python3 calibrate_wgd.py ./rbh_wgd_diag/duplication_blocks.tsv \
    --plot --plot-out FigureS_wgd_thresholds.pdf

# Panel A is the rank-ordered support distribution per taxon, panel B support
# against significance, panel C what the threshold costs you. Pick from those,
# then re-run marking your choices so the figure shows them:
#
#   python3 calibrate_wgd.py ./rbh_wgd_diag/duplication_blocks.tsv \
#       --plot --plot-out FigureS_wgd_thresholds.pdf \
#       --chosen HorseshoeCrab=70 Spider=15 Scorpion=30

# --- Step 3: the real run, with per-species thresholds ---------------------
# SPECIES=N syntax; a bare number still sets the default for everything else.
#
#   --wgd-min-paralog-pairs HorseshoeCrab=70 Spider=15 Scorpion=30
#
# Horseshoe crabs carry multiple rounds of duplication, so they need a far
# higher threshold and are still only partially rescued. Treat their result as
# a lower bound and read duplication_blocks.tsv directly.

echo
echo "Chromosome counts with and without rescue are the check that matters:"
echo "    grep -c '^' rbh_out/dotplots_manifest.tsv   # if present"
echo "and compare dotplots/<taxon>.pdf against dotplots/<taxon>.wgd.pdf -"
echo "rescue should ADD chromosomes that fell below the threshold unrescued."
