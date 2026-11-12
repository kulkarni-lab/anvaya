# Changelog

## 1.0.1b (November 2026)

- New `--threads N` (0 = all CPUs): independent TBR replicates, bootstrap pseudoreplicates and
  permutation replicates run in parallel worker processes. Each unit of work draws from its own
  seed, taken from `--seed` before any work starts, and results are merged in order, so the
  output is identical for any number of threads (`tests/test_determinism.sh`).
- Because of this per-replicate seeding, a given `--seed` gives different (equally valid)
  bootstrap and permutation draws than in 1.0.1. Optimal tree lengths found by exhaustive search
  or branch and bound are unaffected.

## 1.0.1 (November 2026)

First public release. The repository history begins here; the changes below are relative to
the unversioned code described in the preprint (commit 5ba7039, referred to as v0.9).

### Tree search
- New default heuristic: random-addition Wagner trees, TBR branch swapping and the
  parsimony ratchet (`--replicates`, `--ratchet`, `--hold`), on a bit-parallel scorer that
  is checked against the reference Sankoff implementation (`tests/test_engine.py`).
- New `--search exact`: branch and bound seeded with the heuristic score (practical to
  about 14 taxa in pure Python).
- The report now states how many independent replicates reached the best score and warns
  when only one did.
- The pre-1.0 NNI hill-climber is kept as `--search nni` for reproducing old results only.
  On random 16-taxon matrices it stopped several steps above the optimum
  (`benchmarks/search_benchmark.py`).
- Bootstrap and permutation replicates use the main search settings by default
  (`--bootstrap-restarts`, `--permutation-restarts` and the `*-ratchet` options lighten them).
- Backbone constraints: terminal branches now count when deciding whether a required clade
  is present in the restricted tree, so exhaustive and heuristic searches share the documented
  semantics.

### Reproducibility
- `--seed` now fully determines results. Tree adjacency was held in Python sets of taxon
  names, whose order depends on `PYTHONHASHSEED`; it is now insertion-ordered.
  `tests/test_determinism.sh` checks identical output under different hash seeds.
- Every run writes `run_info.json` (version, git commit, exact command, every parameter,
  seed, package versions, SHA-256 of every input) and the same header in `report.txt`.
- `--version` on every script.

### TNT / PAUP* export
- Tree lengths are read from TNT's `length` table and search summary line, tested against a
  real TNT 1.6 log (`tests/data/tnt_ienum_log.txt`). On the Acari dataset the exported scripts
  give 459 (reversible) and 468 (irreversible) steps by implicit enumeration, as Anvaya does.
- TNT: valid search commands (`ienum` up to 20 taxa, otherwise `xmult`), outgroup, backbone,
  `length`, tree files and log; a second script for the irreversible criterion (Sankoff step
  matrix with an all-zero ancestor). `--run-tnt` invokes `tnt "proc <file>;"` and reports
  failure unless trees and a length were actually written, then compares lengths in the report.
- PAUP*: the irreversible search uses a step matrix that matches Anvaya (0->2 costs one
  step); the built-in `irrev.up` charged two. Backbone constraints use PAUP*'s BACKBONE type.
  `--run-paup` runs `paup -n` and checks its output the same way.

### Toy dataset
- Duplicated genes in `T7_wgd` get distinct identifiers (the detector previously found zero
  ohnologs instead of eight).
- Linkage groups I and J are present, unfused, in T1-T6, so the mixed fusion in `T8_mixed`
  is scored against a reference; the mixed-fusion command now has six terminals.
- `MANIFEST.md` states five validations (V1-V5) with explicit expected results;
  `tests/validate_toy.sh` checks every one.

### Packaging and documentation
- `pyproject.toml` (pip-installable; console scripts), `CITATION.cff`, `.zenodo.json`,
  this changelog, and a GitHub Actions workflow that runs the tests on Python 3.9, 3.11 and 3.13.
- README: real clone URL, no reference to the removed tutorial PDF or to scripts that are not
  in the repository, current flag names.
