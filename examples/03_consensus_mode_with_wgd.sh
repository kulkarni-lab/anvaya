#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 3. Consensus mode with WGD lineages, in three steps
# ---------------------------------------------------------------------------
# The recommended route for anything but a quick look. Linkage groups are
# reconstructed from agreement across several core genomes rather than taken
# from one, so no single assembly's rearrangements propagate into every
# character.
#
# Two rules decide the taxon sets:
#   * the CORE set must be non-duplicated, well assembled and phylogenetically
#     spread;
#   * duplicated lineages go on the COMPARISON side with --detect-wgd, never in
#     the core - as core genomes their ohnologous orthogroups are silently
#     discarded as ambiguous.
#
# Step 1 is a permissive diagnostic, step 2 sweeps the thresholds that
# interact, step 3 is the real run at the setting you chose. Do not skip
# step 2: no threshold is correct a priori.
# ---------------------------------------------------------------------------
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
# sweep_consensus_params.py does this far faster than a nested loop: it caches
# the diamond alignments across the whole sweep and derives the
# --consensus-min-alg-orthogroups axis arithmetically from one reconstruction.
python3 sweep_consensus_params.py \
    --input-dir ./Data \
    --consensus-species $CORE \
    --outdir ./sweep \
    --min-species 3 4 \
    --min-frac 1.0 0.8 0.7 0.5 \
    --min-alg-orthogroups 1 5 10 15 20 30 \
    --threads 30 \
    --plot --plot-out ./Figure1_thresholds.pdf \
    --chosen-min-species 3 --chosen-min-frac 0.7 --chosen-min-alg-orthogroups 15

column -t ./sweep/consensus_sweep_summary.tsv

# --- Step 3: full analysis at the chosen setting ---------------------------
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

echo
echo "Now cross-check the search in TNT or PAUP* using the exported matrix -"
echo "see examples/tnt/ and the Caveats section of the README."
