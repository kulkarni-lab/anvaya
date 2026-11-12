#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Run the five toy-dataset validations (V1-V5 in toy_data/MANIFEST.md) and
# check every stated expectation. Needs DIAMOND on PATH.
#
#     bash tests/validate_toy.sh [workdir]
#
# Exits non-zero at the first expectation that is not met.
# ---------------------------------------------------------------------------
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK="${1:-$(mktemp -d)}"
mkdir -p "$WORK"
cd "$WORK"
command -v diamond >/dev/null || { echo "diamond not found on PATH"; exit 2; }

G="python3 $ROOT/genomes_to_rbh.py"
P="--parsimony-script $ROOT/alg_parsimony.py"
python3 "$ROOT/make_toy_dataset.py" --outdir toy_data >/dev/null

echo "== V1 single-reference mode =="
$G --input-dir toy_data --reference-species T1 --species T1 T2 T3 T4 T5 T6 \
   --rbh-outdir v1_rbh --parsimony-outdir v1_out $P \
   --forward --min-genes 3 --alpha 0.05 --outgroup T2 --seed 1 > v1.log 2>&1
echo "== V2 consensus mode =="
$G --input-dir toy_data --consensus-species T1 T3 T5 --species T1 T2 T3 T4 T5 T6 \
   --rbh-outdir v2_rbh --parsimony-outdir v2_out $P \
   --forward --min-genes 3 --alpha 0.05 --outgroup T1 --seed 1 \
   --bootstrap 100 --permutation-test 100 > v2.log 2>&1
echo "== V3 whole-genome duplication rescue =="
$G --input-dir toy_data --reference-species T1 --species T1 T7_wgd \
   --rbh-outdir v3_rbh_plain --rbh-only > v3a.log 2>&1
$G --input-dir toy_data --reference-species T1 --species T1 T7_wgd \
   --detect-wgd --wgd-species T7_wgd --rbh-outdir v3_rbh_wgd --rbh-only > v3b.log 2>&1
echo "== V4 minimum scaffold length =="
$G --input-dir toy_data --reference-species T1 --species T1 T6 \
   --rbh-outdir v4_rbh_plain --rbh-only > v4a.log 2>&1
$G --input-dir toy_data --reference-species T1 --species T1 T6 \
   --min-scaffold-length 1000 --rbh-outdir v4_rbh_minscaf --rbh-only > v4b.log 2>&1
echo "== V5 mixed fusion =="
$G --input-dir toy_data --reference-species T1 --species T1 T2 T3 T4 T5 T6 T8_mixed \
   --rbh-outdir v5_rbh --parsimony-outdir v5_out $P \
   --forward --min-genes 3 --alpha 0.05 --outgroup T2 --seed 1 \
   --plot-dotplots --plot-ideograms > v5.log 2>&1

python3 - "$ROOT" <<'PY'
import csv, glob, os, re, sys
sys.path.insert(0, sys.argv[1])
import alg_parsimony as ap

def check(cond, msg):
    print(("  ok   " if cond else "  FAIL ") + msg)
    if not cond:
        sys.exit(1)

def has_clade(nwk_path, clade):
    g = ap.node_tree_to_graph(ap.parse_newick(open(nwk_path).read().strip()))
    return any(frozenset(clade) in bip for bip in ap.bipartitions(g))

def report(path):
    return open(path).read()

def matrix(path):
    rows = list(csv.reader(open(path), delimiter="\t"))
    head = rows[0][1:]
    return {r[0]: dict(zip(head, r[1:])) for r in rows[1:]}

def ci_ri(path):
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    return {r["character"]: (float(r["CI"]), float(r["RI"])) for r in rows}

print("== checking expectations ==")
# V1
r = report("v1_out/report.txt")
n_mpt = re.findall(r"number of equally parsimonious trees found: (\d+)", r)
check(n_mpt == ["2", "1"], f"V1 two MPTs (reversible), one (irreversible): {n_mpt}")
check(r.count("best score: 5") == 2, "V1 length 5 under both criteria")
g = ap.node_tree_to_graph(ap.parse_newick(open("v1_out/tree_irreversible.nwk").read().strip()))
true = ap.node_tree_to_graph(ap.parse_newick("(T2,((T3,T4),(T5,T6)));"))
check(ap.bipartitions(g) == ap.bipartitions(true), "V1 irreversible tree = (T2,((T3,T4),(T5,T6)))")
check(has_clade("v1_out/tree_reversible.nwk", {"T3", "T4"}), "V1 reversible tree has (T3,T4)")

# V2
r = report("v2_out/report.txt")
check(r.count("number of equally parsimonious trees found: 1") == 2, "V2 one MPT per criterion")
check(r.count("best score: 5") == 2, "V2 length 5 under both criteria")
check(r.count("ensemble CI: 0.800   ensemble RI: 0.750") == 2, "V2 ensemble CI 0.800, RI 0.750")
for crit in ("reversible", "irreversible"):
    nwk = open(f"v2_out/tree_{crit}.nwk").read()
    g = ap.node_tree_to_graph(ap.parse_newick(nwk.strip()))
    true = ap.node_tree_to_graph(ap.parse_newick("(T1,(T2,((T3,T4),(T5,T6))));"))
    check(ap.bipartitions(g) == ap.bipartitions(true), f"V2 {crit} tree = (T1,(T2,((T3,T4),(T5,T6))))")
    cr = ci_ri(f"v2_out/ci_ri_{crit}.tsv")
    homopl = [k for k, v in cr.items() if v[0] < 1]
    check(len(homopl) == 1 and cr[homopl[0]] == (0.5, 0.0),
          f"V2 {crit}: exactly one homoplastic character (CI 0.5, RI 0.0): {homopl}")

# V3
def nrows(p):
    return sum(1 for _ in open(p)) - 1
check(nrows("v3_rbh_plain/T7_wgd.rbh") == 48, "V3 48 orthologs without --detect-wgd")
check(nrows("v3_rbh_wgd/T7_wgd.rbh") == 96, "V3 96 orthologs with --detect-wgd")
blocks = list(csv.DictReader(open("v3_rbh_wgd/duplication_blocks.tsv"), delimiter="\t"))
check(len(blocks) == 8 and all(int(b["n_paralog_pairs"]) == 6 for b in blocks),
      "V3 8 ohnologous scaffold pairs, 6 paralog pairs each")
types = [row["orthology_type"] for row in csv.DictReader(open("v3_rbh_wgd/T7_wgd.rbh"), delimiter="\t")]
check(types.count("wgd_inferred") == 48, "V3 48 orthologs flagged wgd_inferred")

# V4
check(nrows("v4_rbh_minscaf/T6.rbh") == nrows("v4_rbh_plain/T6.rbh") - 1,
      "V4 --min-scaffold-length 1000 removes the ALG-N ortholog")

# V5
m = matrix("v5_out/character_matrix.tsv")
ij = [k for k in next(iter(m.values())) if k.split("+")[0].endswith("I") and k.split("+")[1].endswith("J")]
check(len(ij) == 1, f"V5 one I+J character ({ij})")
check(m["T8_mixed"][ij[0]] == "2", "V5 I+J coded 2 (fused-and-mixed) in T8_mixed")
check(all(m[t][ij[0]] == "0" for t in ("T2", "T3", "T4", "T5", "T6")), "V5 I+J coded 0 in T2-T6")
print("all five toy validations passed")
PY
