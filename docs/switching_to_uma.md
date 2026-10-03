# Switching to UMA

`chemsmart guess` is experimental. It is not yet better than hand-built
guesses. On the xTB benchmark, the hand-built guesses beat xTB on all
3 cases. The chemistry comparison is deferred until this calculator
is available.

`chemsmart guess --calc uma` is the calculator this command is built for.
xTB (`--calc xtb`) is only there so the scan can be run before a Hugging
Face token exists. It is not a claim about speed or about agreement
with MN15.

UMA is the `facebook/UMA` checkpoint (default `uma-s-1p1`), loaded by
fairchem as an ASE `FAIRChemCalculator` with task `omol`. Charge comes
from the case YAML. Spin is the YAML multiplicity (1 for these
singlets), written to `atoms.info["spin"]`.

## Steps

1. Accept the licence at <https://huggingface.co/facebook/UMA>.
2. Create a token at <https://huggingface.co/settings/tokens>.
3. Export it in the shell you use for chemsmart. The process reads
   `HF_TOKEN` from the environment only. Do not put the token in a
   YAML file, a command line, or a commit. The commands never print it.

   ```bash
   export HF_TOKEN=hf_...
   ```

4. Install the optional extra. It is not part of the core dependencies.

   ```bash
   pip install 'chemsmart[mlip]'
   ```

5. Check the install and the gated checkpoint. This does not evaluate
   an energy and does not start a scan.

   ```bash
   chemsmart guess --calc uma --check-access
   ```

   A good check ends with `UMA access ok`, `task: omol`, and
   `HF_TOKEN: set`.

6. Run a guess from the reactant or precomplex and the case YAML.
   `scan.to` is the product-side length of the primary bond. Chi must
   set that number. There is no default. Templates leave it `unset`.

   ```bash
   chemsmart guess -f reactant.xyz --spec case.yaml --calc uma
   ```

   Add `--sella` only when you want a report-only Sella order-1 run.
   That writes `guess_sella.xyz` and does not replace `guess.xyz`.
   The default guess is the converged scan maximum itself.

If fairchem is missing, `HF_TOKEN` is unset, or Hugging Face rejects
the gated repo, the command stops with an install or access message
and does not run the scan.

xTB, when you still want the plumbing check:

```bash
chemsmart guess -f reactant.xyz --spec case.yaml --calc xtb
```

Charge and uhf come from the case YAML (`uhf = multiplicity - 1`, so
0/1 is charge 0, uhf 0) or from `-c` / `--uhf`.
