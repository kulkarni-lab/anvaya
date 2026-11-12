#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Same --seed  =>  identical results, whatever Python's string-hash seed.
# ---------------------------------------------------------------------------
# Runs the same analysis twice (and with 1 and 4 worker processes), in separate Python processes with different
# PYTHONHASHSEED values, once under exhaustive search and once under the
# heuristic, with bootstrap and permutation replicates, and requires every
# result file to be byte-identical. (run_info.json and the report header carry
# a timestamp and are compared after removing it.)
#
#     bash tests/test_determinism.sh
# ---------------------------------------------------------------------------
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
cd "$ROOT"
cp mock_figures/make_mock_rbh.py "$WORK/"
( cd "$WORK" && python3 make_mock_rbh.py >/dev/null )
RBH="$WORK/mockfigs/rbh_main"

for search in exhaustive tbr; do
  for hs in 1 4242; do
    PYTHONHASHSEED=$hs python3 alg_parsimony.py --indir "$RBH" --outdir "$WORK/$search.$hs" \
        --min-genes 3 --outgroup Sarc --search "$search" --replicates 5 --ratchet 3 \
        --bootstrap 20 --permutation-test 20 --seed 7 >/dev/null 2>&1
  done
  for f in tree_reversible.nwk tree_irreversible.nwk bootstrap_reversible.tsv \
           bootstrap_irreversible.tsv permutation_test_reversible.tsv \
           permutation_test_irreversible.tsv synapomorphies_reversible.txt \
           character_matrix.tsv; do
    cmp -s "$WORK/$search.1/$f" "$WORK/$search.4242/$f" \
      || { echo "FAIL: $search/$f differs between PYTHONHASHSEED values"; exit 1; }
  done
  grep -v -e '^Started:' -e '^Command:' -e '--outdir ' "$WORK/$search.1/report.txt" > "$WORK/a"
  grep -v -e '^Started:' -e '^Command:' -e '--outdir ' "$WORK/$search.4242/report.txt" > "$WORK/b"
  cmp -s "$WORK/a" "$WORK/b" || { echo "FAIL: $search report differs"; diff "$WORK/a" "$WORK/b" | head; exit 1; }
  echo "ok   $search search: identical outputs under PYTHONHASHSEED=1 and 4242"
done

# ---- same --seed, different --threads => identical results ----
for threads in 1 4; do
  python3 alg_parsimony.py --indir "$RBH" --outdir "$WORK/threads.$threads" \
      --min-genes 3 --outgroup Sarc --search tbr --replicates 5 --ratchet 3 \
      --bootstrap 20 --permutation-test 20 --seed 7 --threads $threads >/dev/null 2>&1
done
for f in tree_reversible.nwk tree_irreversible.nwk bootstrap_reversible.tsv \
         bootstrap_irreversible.tsv permutation_test_reversible.tsv \
         permutation_test_irreversible.tsv; do
  cmp -s "$WORK/threads.1/$f" "$WORK/threads.4/$f" \
    || { echo "FAIL: $f differs between --threads 1 and --threads 4"; exit 1; }
done
echo "ok   tbr search: identical outputs with --threads 1 and --threads 4"
