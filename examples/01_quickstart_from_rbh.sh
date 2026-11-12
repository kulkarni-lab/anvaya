#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 1. Quick start: you already have .rbh tables
# ---------------------------------------------------------------------------
# The fastest way to see what Anvaya does. Needs no diamond and no genomes -
# just a directory of .rbh files, one per compared taxon.
#
# Run from the repository root:
#     bash examples/01_quickstart_from_rbh.sh ./rbh_files MyOutgroup
# ---------------------------------------------------------------------------
set -euo pipefail

INDIR="${1:-./rbh_files}"
OUTGROUP="${2:?usage: $0 <rbh_dir> <outgroup_taxon>}"

python3 alg_parsimony.py \
    --indir "$INDIR" \
    --outdir ./results \
    --min-genes 3 --alpha 0.05 \
    --outgroup "$OUTGROUP" \
    --plot-dotplots --plot-trees --plot-matrix

echo
echo "Read ./results/report.txt first - it carries the scores, the"
echo "synapomorphies and every warning the run produced."
