#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 4. Sensitivity analysis: are the conclusions stable across thresholds?
# ---------------------------------------------------------------------------
# A linkage-group count that looks plausible is not evidence that the threshold
# was right. What matters is whether the TOPOLOGY and the supported clades hold
# across the plausible range. Run this, and report the range you tested
# alongside the result.
#
# This loop varies --consensus-min-alg-orthogroups. Repeat it for --min-genes
# and --alpha, which control the fusion CALLS rather than the linkage groups,
# by re-running alg_parsimony.py on one .rbh set at several values.
# ---------------------------------------------------------------------------
set -euo pipefail

CORE="SeaSpider Tick Mite Harvestman"
THRESHOLDS="10 15 20 30"

for minalg in $THRESHOLDS; do
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
for minalg in $THRESHOLDS; do
  printf "%-8s " "t=$minalg"
  cat "./parsimony_t${minalg}/tree_reversible.nwk"
done

echo
echo "=== character counts and ensemble fit across thresholds ==="
printf "%-8s %-10s %s\n" "thresh" "chars" "ensemble CI / RI"
for minalg in $THRESHOLDS; do
  d="./parsimony_t${minalg}"
  nchar=$(($(head -1 "$d/character_matrix.tsv" | tr '\t' '\n' | wc -l) - 1))
  fit=$(grep -E "ensemble CI" "$d/report.txt" | head -1)
  printf "%-8s %-10s %s\n" "t=$minalg" "$nchar" "$fit"
done

echo
echo "Columns 5 and 6 of ci_ri_*.tsv are CI and RI respectively - a frequent"
echo "off-by-one. Characters with CI < 0.5 are the homoplastic ones:"
echo "    sort -k5,5g parsimony_t15/ci_ri_reversible.tsv | head"
