#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 5. Test the installation end to end on synthetic data
# ---------------------------------------------------------------------------
# make_toy_dataset.py writes eight small genomes whose true tree, fusions and
# one whole-genome duplication are known and recorded in toy_data/MANIFEST.md.
# If this script reproduces the manifest's tree, your install - diamond
# included - is working.
#
#     bash examples/05_toy_dataset_test.sh
# ---------------------------------------------------------------------------
set -euo pipefail

python3 make_toy_dataset.py --outdir toy_data      # see toy_data/MANIFEST.md
python3 make_toy_dataset.py --verify-only          # sanity-check sequences only

# --- single-reference mode --------------------------------------------------
# In this mode the reference genome (T1) is NOT a taxon in the tree, so it
# cannot also be the --outgroup. T2 is the correct outgroup here.
python3 genomes_to_rbh.py \
    --input-dir toy_data --reference-species T1 \
    --species T1 T2 T3 T4 T5 T6 \
    --rbh-outdir toy_rbh \
    --parsimony-script "$(pwd)/alg_parsimony.py" \
    --parsimony-outdir toy_parsimony \
    --forward --min-genes 3 --alpha 0.05 --outgroup T2 \
      --bootstrap 100 --permutation-test 100 --plot-trees

# --- consensus mode --------------------------------------------------------
# Here T1 IS a taxon, so it can serve as the outgroup.
python3 genomes_to_rbh.py \
    --input-dir toy_data --consensus-species T1 T3 T5 \
    --species T1 T2 T3 T4 T5 T6 \
    --rbh-outdir toy_rbh_consensus \
    --parsimony-script "$(pwd)/alg_parsimony.py" \
    --parsimony-outdir toy_parsimony_consensus \
    --forward --min-genes 3 --alpha 0.05 --outgroup T1

echo
echo "=== recovered vs. designed ==="
echo "single-reference: $(cat toy_parsimony/tree_reversible.nwk)"
echo "consensus:        $(cat toy_parsimony_consensus/tree_reversible.nwk)"
echo "designed:         see the 'True tree' line in toy_data/MANIFEST.md"
