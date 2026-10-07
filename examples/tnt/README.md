# Cross-checking the search in TNT

Anvaya's heuristic search above `--max-exhaustive` taxa is an NNI
hill-climber, not a validated general-purpose search. For anything that
matters, export the matrix and verify the score in TNT or PAUP\*:

```bash
python3 alg_parsimony.py --indir ./rbh_out --outdir ./parsimony_out \
    --outgroup MyOutgroup --export-tnt --export-nexus
```

That writes `character_matrix.tnt` and `character_matrix.nex`.

## Preparing the matrix files

`character_matrix.tnt` is self-contained: it carries its own search commands
after the matrix so it runs on its own. These scripts drive the search
themselves, so they need the **data block only** — everything from `xread`
down to the `;` that closes it, and nothing after.

```bash
# matrix_data.tnt  - for script_anvaya.run
sed '/^;$/q' character_matrix.tnt > examples/tnt/matrix_data.tnt
```

For the irreversible script, additionally append an all-zero `Ancestor` taxon
and bump the taxon count by one:

```bash
# matrix_irrev.tnt - for script_anvaya_irreversible.run
python3 - <<'PY'
lines = open("character_matrix.tnt").read().split("\n")
end = lines.index(";")                       # first bare ';' closes the matrix
nchar, ntax = lines[2].split()
lines[2] = f"{nchar} {int(ntax) + 1}"
lines.insert(end, "Ancestor " + "0" * int(nchar))   # the ';' shifts to end + 1
open("examples/tnt/matrix_irrev.tnt", "w").write("\n".join(lines[:end + 2]) + "\n")
PY
```

Then edit the `force =` line in each script to your own backbone constraint,
and the `0.324` character range in the irreversible script to `nchar - 1`.

| Script | What it does |
|---|---|
| `script_anvaya.run` | Reversible (Fitch) implied-enumeration search under a backbone constraint, then a constrained bootstrap of 1000 replicates. |
| `script_anvaya_irreversible.run` | The same under a Camin-Sokal (irreversible) model, which TNT has no native type for — see below. |

Run either with `tnt proc script_anvaya.run;` or by typing `proc
script_anvaya.run;` at the TNT prompt.

## Things that will bite you

These are all mistakes made while writing these scripts, in the order they
surfaced:

- **Comment syntax depends on your TNT build.** Some versions take
  `/* … */`, others only `[ … ]`. If the parser dies on line 1, that is why.
- **`nstates` must precede `proc`.** Set it after the matrix is read and it is
  ignored.
- **`tsave` writes binary.** `tsave *file.tre` writes the parenthetical form;
  without the asterisk you get an unreadable file.
- **`force` only declares the constraint; `constrain=;` activates it.** A
  script with `force` and no `constrain=;` runs completely unconstrained and
  silently returns a better score — which is exactly how a "the search is
  broken" conclusion gets manufactured.
- **`outgroup` must come before `force`.** TNT roots on its first taxon
  otherwise, and if that taxon sits inside a constrained group the constraint
  is unsatisfiable.
- **Compare constrained against constrained.** A constrained Anvaya score and
  an unconstrained TNT score are not comparable, and the difference between
  them is the cost of the constraint, not a deficiency in the search.

## The irreversible model

TNT has no native irreversible character type and optimises the root state, so
Camin-Sokal cannot be requested directly. Two things together emulate it:

1. a **step matrix** making every 1→0 and 2→1 reversal prohibitively
   expensive, and
2. a dummy **`Ancestor`** taxon scored all-zero, which pins the root to the
   unfused state.

Score the result against Anvaya's `tree_irreversible.*` output, not against
`tree_reversible.*`. The `Ancestor` taxon is an artefact of the encoding and is
not part of the result — prune it before reading the topology.

## Reading the result

TNT reports the score; Anvaya reports its own in `report.txt`. They should
match exactly when the constraint and the model match. If TNT's is lower, the
topology it found is the better one — check whether Anvaya emitted a
`[warn] the … search returned score …` line, which is its own cross-check
between the reversible and irreversible searches firing.
