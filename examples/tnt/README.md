# Cross-checking the search in TNT

From version 1.0, `alg_parsimony.py --export-tnt` writes two self-contained TNT
scripts and `--run-tnt` runs them and compares the lengths with Anvaya's:

```bash
python3 alg_parsimony.py --indir ./rbh_out --outdir ./parsimony_out \
    --outgroup MyOutgroup --constraint backbone.nwk --run-tnt
```

| File | What it does |
|---|---|
| `character_matrix.tnt` | Reversible criterion: all characters non-additive (`cc - .`), outgroup, backbone constraint (`force` + `constrain =`), exact search by implicit enumeration (`ienum`) up to 20 taxa, otherwise `xmult` + TBR. Reports `length`, saves trees in parenthetical form and writes a log. |
| `character_matrix_irreversible.tnt` | Irreversible criterion: Sankoff step matrix (any increase 1 step, any decrease 1000), an all-zero `Ancestor_allzero` taxon as outgroup and the ingroup constrained monophyletic. |

Run either by hand with

```bash
cd parsimony_out
tnt "proc character_matrix.tnt;"
```

`--run-tnt` does exactly this, and counts a run as successful only if the tree
file and a tree length were actually written; it then adds an "External
cross-check" section to `report.txt`.

## Comparing irreversible scores

TNT optimises the root state, so Camin-Sokal-type scoring is emulated with the
all-zero ancestor taxon. Anvaya additionally fixes the common ancestor of the
outgroup and the ingroup at state 0. For every character in which the outgroup
has the lowest ingroup state and that state is above 0, the TNT emulation is one
step shorter. This difference does not depend on the topology, so both programs
rank trees identically; the offset is printed in the script (`quote` line) and
added back by `--run-tnt` before comparing.

## The hand-written scripts in this directory

`script_anvaya.run` is the script used for the constrained bootstrap reported in
the paper. It reads the data block of an exported matrix:

```bash
sed -n '/^xread/,/^;$/p' character_matrix.tnt > examples/tnt/matrix_data.tnt
```

Edit its `force =` line to your own backbone constraint.

## Things that will bite you

- **Comment syntax depends on your TNT build.** Some versions take `/* … */`,
  others only `[ … ]`. If the parser dies on line 1, that is why.
- **`nstates` must precede the matrix.**
- **`tsave` writes binary.** `tsave *file.tre` writes the parenthetical form;
  without the asterisk you get an unreadable file.
- **`force` only declares the constraint; `constrain =;` activates it.** A
  script with `force` and no `constrain =;` runs unconstrained and returns a
  shorter tree.
- **`outgroup` must come before `force`.**
- **Compare constrained against constrained.** A constrained Anvaya score and
  an unconstrained TNT score differ by the cost of the constraint.
- **TNT exits with code 0 even when a script fails.** Check the log, or use
  `--run-tnt`, which does.
