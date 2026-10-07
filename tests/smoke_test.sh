#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Smoke test: does this checkout work at all?
# ---------------------------------------------------------------------------
# Needs pandas and matplotlib. Does NOT need diamond - it builds synthetic
# .rbh tables directly with a designed, known structure, so the whole
# parsimony side is exercised without an orthology search.
#
#     bash tests/smoke_test.sh
#
# Exits non-zero on the first failure. This is what CI runs.
# ---------------------------------------------------------------------------
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
cd "$ROOT"

pass() { printf '  \033[32mok\033[0m   %s\n' "$1"; }
fail() { printf '  \033[31mFAIL\033[0m %s\n' "$1"; exit 1; }

echo "== 1. every script parses and prints usage =="
for f in alg_parsimony.py genomes_to_rbh.py sweep_consensus_params.py \
         calibrate_wgd.py make_toy_dataset.py; do
    python3 "$f" --help >/dev/null 2>&1 && pass "$f --help" || fail "$f --help"
done

echo "== 2. synthetic sequence data is internally consistent =="
python3 make_toy_dataset.py --verify-only >/dev/null 2>&1 \
    && pass "make_toy_dataset.py --verify-only" \
    || fail "make_toy_dataset.py --verify-only"

echo "== 3. end-to-end parsimony run on designed .rbh tables =="
cp mock_figures/make_mock_rbh.py "$WORK/"
( cd "$WORK" && python3 make_mock_rbh.py >/dev/null )
RBH="$WORK/mockfigs/rbh_main"
[ -s "$RBH/Ixod.rbh" ] || fail "mock .rbh tables were not written"
pass "mock .rbh tables written"

python3 alg_parsimony.py \
    --indir "$RBH" --outdir "$WORK/out" \
    --min-genes 3 --alpha 0.05 --outgroup Sarc \
    --plot-trees --plot-matrix --plot-ribbon \
    --bootstrap 50 --permutation-test 50 \
    --export-nexus --export-tnt >"$WORK/run.log" 2>&1 \
    || { cat "$WORK/run.log"; fail "alg_parsimony.py run"; }
pass "alg_parsimony.py completed"

echo "== 4. expected outputs exist =="
for f in report.txt character_matrix.tsv tree_reversible.nwk \
         tree_irreversible.nwk ci_ri_reversible.tsv \
         bootstrap_reversible.tsv character_matrix.nex character_matrix.tnt \
         ribbon.pdf tree_reversible.pdf character_matrix_heatmap.pdf; do
    [ -s "$WORK/out/$f" ] && pass "$f" || fail "$f missing or empty"
done

echo "== 5. figure text is editable, not outlined =="
# pdf.fonttype = 42 embeds TrueType subsets (/FontFile2). Type 3 would mean
# the text came out as unselectable outline shapes, which breaks editing in
# Illustrator and is a regression we have hit before.
python3 - "$WORK/out/ribbon.pdf" "$WORK/out/tree_reversible.pdf" <<'PY'
import sys
bad = []
for p in sys.argv[1:]:
    d = open(p, "rb").read()
    if b"/FontFile2" not in d or b"/Type3" in d:
        bad.append(p)
if bad:
    sys.exit("outlined or Type 3 text in: " + ", ".join(bad))
PY
pass "TrueType subsets embedded (pdf.fonttype = 42)"

echo "== 6. the designed phylogenetic structure is recovered =="
# make_mock_rbh.py designs cALG1+cALG2 as a synapomorphy of (Ixod,Varr,Neoseiu)
# and cALG3+cALG4 of (Varr,Neoseiu), with cALG7+cALG8 as deliberate homoplasy.
# The run must recover the first clade and report the homoplasy as homoplasy.
grep -q "Ixod, Neoseiu, Varr" "$WORK/out/report.txt" \
    && pass "clade (Ixod,Varr,Neoseiu) recovered" \
    || fail "designed clade (Ixod,Varr,Neoseiu) not recovered"
grep -q "cALG1+cALG2" "$WORK/out/synapomorphies_reversible.txt" \
    && pass "cALG1+cALG2 mapped as a synapomorphy" \
    || fail "cALG1+cALG2 not mapped as a synapomorphy"
awk -F'\t' 'NR>1 && $5+0 < 1 {n++} END {exit !(n>0)}' \
    "$WORK/out/ci_ri_reversible.tsv" \
    && pass "homoplasy detected (CI < 1 for at least one character)" \
    || fail "no homoplasy detected, but the dataset contains a designed one"

echo
printf '\033[32mall checks passed\033[0m\n'
