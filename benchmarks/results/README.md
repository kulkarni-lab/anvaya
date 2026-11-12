# Benchmark results reported in the paper (Anvaya v1.0.1b, Table 2)

| File | Command |
|---|---|
| `taxa12_exact.tsv` | `python3 benchmarks/search_benchmark.py --datasets 20 --taxa 12 --chars 20 --legacy-restarts 10 --threads 16` (seed 2026) |
| `taxa16_tnt-design.tsv` | `python3 benchmarks/search_benchmark.py --datasets 4 --taxa 16 --chars 20 --no-exact --long 100 --legacy-restarts 10 500 --criterion reversible --threads 16 --tnt-bin tnt` (seed 2026) |

Columns: `tbr_*` default search (20 replicates, 10 ratchet iterations); `exact_score` branch and bound;
`long_tbr_*` independent search with 100 replicates and 20 ratchet iterations; `nniN_*` the
pre-1.0 nearest-neighbour search with N starts; `tnt_score` TNT v1.6 implicit enumeration (Linux).
The matrices depend only on the benchmark seed, so the TNT lengths, obtained with the same command,
apply to these matrices unchanged.
