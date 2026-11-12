#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 5. Test the installation end to end on synthetic data (needs DIAMOND)
# ---------------------------------------------------------------------------
# make_toy_dataset.py writes eight small genomes whose tree, fusions, mixed
# fusion and whole-genome duplication are known; toy_data/MANIFEST.md lists
# five validations (V1-V5), the command for each and the result it must give.
# tests/validate_toy.sh runs all five and checks every expectation.
#
#     bash examples/05_toy_dataset_test.sh
# ---------------------------------------------------------------------------
set -euo pipefail
python3 make_toy_dataset.py --verify-only          # sequences behave as designed
bash tests/validate_toy.sh ./toy_validation       # V1-V5; non-zero exit on any failure
