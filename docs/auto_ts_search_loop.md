# Loop.md output (folded into ts-scout pass 2)

The design is the critic-accepted **PASS 2 FINAL** in the ts-scout domain-bot
spec. This file records how that spec was implemented in the chemsmart fork.
It is not a parallel design.

## PROBLEM

Hands-on time per verified TS is too high. A verified TS is: nimag==1, hard
mode checks pass, both endpoints (IRC, or QRC as a fallback) optimized to
nimag=0 and meeting the case YAML endpoint bonds, and Chi gives the final OK.
A freq PASS means "candidate, queue IRC", never "verified".

## DESIGN (PASS 2 FINAL)

Two small pieces in the chemsmart fork. **No new job types.**

1. `chemsmart check <log> --spec case.yaml`
2. `chemsmart queue ts|irc|opt` writes a shell file and never submits.

Case YAML is written by Chi. QRC is a printed suggestion only.

## METRIC

Count Chi's manual actions per verified TS.

- Baseline (estimate): 12 actions without QRC, 15 with QRC.
- Measured replay (pass 3; those logs are not in this repository): A is 6 manual touchpoints against a
  baseline of 12, or 7–8 counting the case YAML. C is 7 against a baseline
  of 15, or 8 counting the YAML.
- Extra that the raw count leaves out: writing the case YAML is +1 per case;
  each command is still typed by Chi; nothing chains `check` into the next
  queue file. For A, a template product cutoff of < 1.5 Å false-flagged the
  TS1b IRC end (C8–N29 = 1.574 Å, which is INT2). That cutoff is removed.
  Until Chi sets thresholds, `check` prints `endpoint criteria unset` and
  does not add a connectivity-mismatch touchpoint.

## ACCEPTANCE CRITERIA (from PASS 2 FINAL / user)

1. `check` prints nimag, imag freq, charge/mult, top-8 pair projections with
   distances, spec-bond ranks, Zmax (info only). The table header notes
   `projections from 2-decimal displacements`.
2. Hard TS rejects that `--force` does not override: nimag=0, and nimag≥2
   when an extra imaginary mode has |freq| ≥ 15 cm⁻¹. Extra modes that are
   all under 15 cm⁻¹ print `small imag (<15): judge`, list each small mode,
   and stay `CANDIDATE: judge mode`. Primary spec bond not #1 within 3.2 Å
   is a reject `--force` overrides; rank printed on that path. If every
   spec-bond projection is about 0, the reject is `spec bonds not displaced
   in mode` and does not quote a rank. `--force` overrides that reject too.
3. Route and 0/1 charge/mult checked against the three expected MN15/def2svp
   routes. Neg1-style routes are flagged. Flags do not turn a CANDIDATE
   into REJECT.
4. Endpoint-opt: nimag≠0 prints `ENDPOINT NOT MIN` (not REJECT). Bond-criteria
   failure prints `connectivity mismatch: judge` (not REJECT) and one line
   per failing bond (`C8–N29: 1.574 Å, needs < 1.5 Å (product)`). While the
   YAML placeholder is unset, `check` prints `endpoint criteria unset` and
   skips the bond verdict. Case 1 thresholds stay unset.
5. table printed to stdout
6. No new job types. No MLIP code. No early-kill. No custom Gaussian
   templates. Do not edit ts-scout `calc/chemsmart/zn5.yaml`. `queue` still
   writes, never submits (see DESIGN). QRC is printed only.
7. Tests: synthetic Gaussian snippets, including
   `tests/data/tscheck/synthetic_irc_minimum.log` (invented points and a
   PES-minimum line), plus the public chemsmart fixture
   `pd_genecp_ts.log`. No research logs are committed.
8. Draft PR, author AsterGemini only.

## Passes 1–2 (ts-scout; accepted)

Evidence: pass2_evidence.md (A/B/C, Neg1, Neg2, IRC routes, early-kill false
flags, UMA gated, zn5.yaml stale vs live cluster yaml).

Deleted: early-kill, auto-reject on a fragile Zmax rule, custom templates,
custom QRC, MLIP on the critical path.

Kept manual: case YAML, mode judgment, every submit, final OK.

## Pass 3 (this builder implementation)

Folded the spec into `chemsmart check` / `chemsmart queue`. Did not re-open
wrap-vs-rebuild of autodE/pysisyphus/Sella/UMA.

### What was automated

| Automated | Why it survived |
| --- | --- |
| Parse nimag / imag / charge/mult / top-8 / spec ranks / Zmax | Mechanical; Zmax is printed only |
| Hard reject nimag=0, or nimag≥2 with an extra mode \|freq\| ≥ 15 cm⁻¹ | `--force` does not override |
| `small imag (<15): judge` when every extra mode is under 15 cm⁻¹ | Not a hard reject; verdict stays CANDIDATE |
| Reject primary bond not #1 within 3.2 Å (`--force` override) | Rank printed on this path only |
| Reject `spec bonds not displaced in mode` when every spec projection is ~0 (`--force` override) | No meaningful rank; override does not quote one |
| Route and 0/1 charge/mult flag | Catches Neg1's old route without auto-rejecting |
| `ENDPOINT NOT MIN` when an endpoint opt has nimag≠0 | Judge state, not a hard reject |
| `connectivity mismatch: judge` plus one line per failing bond | Judge state, not a hard reject |
| `endpoint criteria unset` while the YAML placeholder is unset | Skips a false bond verdict |
| Write (not submit) ts / irc / opt command files | File shuffling only |
| Print QRC suggestion when IRC points < 5 | Diagnosis only; not queued. Point count is not a minimum |
| `IRC did not reach a minimum` when ≥5 points and no PES-minimum string | Matches `PES minimum detected` and the older `minimum found` |

### What was deliberately left manual

| Manual | Reason |
| --- | --- |
| Case YAML (which bonds, endpoint criteria) | Human judgment of the reaction coordinate |
| Mode table judgment | Evidence: no robust auto-rule separates A from Neg1 |
| `--force` | Human override of the 3.2 Å rank reject |
| Every submit / running the queue file | Real compute, opt-in |
| Final OK | Definition of a verified TS |
| Cluster zn5.yaml | Flag for Chi; not edited here |
| QRC run | Printed suggestion only |
| UMA/MLIP | Gated HF repo; deferred |
| Early-kill | Would false-flag B and C |

## Benchmark in this environment

Pass 3 was measured on logs that are not in this repository. Functional
results:

- A TS1b −175.16, case1.yaml: CANDIDATE, C8–N29 #1, Zmax 0.1623 Zn21–N29.
- B TS2 −170.46, case1_ts2.yaml: CANDIDATE, O1–C8 #1, Zmax 0.3338 C8–Zn21.
- C case3 TS2 −84.61: CANDIDATE, C113–O114 #1, Zmax 0.2768 Zn5–C113.
- Neg1 −219.14: CANDIDATE with ROUTE FLAG (old `#p … MN15/Def2SVP`).
- Neg2 rotor −89.82: spec-bond projections are 0. The reject is
  `spec bonds not displaced in mode` with no rank (a distance tie-break
  previously printed rank 76).
- TS1b forward IRC (not committed): 23 points and `PES minimum
  detected`. With case1 endpoints unset, `check` prints `endpoint
  criteria unset`, no QRC, not REJECT. The old product cutoff < 1.5 Å
  was wrong: the IRC end is C8–N29 = 1.574 Å and that geometry is INT2.
  The file in this repo, `tests/data/tscheck/synthetic_irc_minimum.log`,
  is an invented fixture, not that log.
- Case 3 one-point IRC (also prints `PES minimum detected`): still
  `suggest qrc ±` because the point count is under 5.

Measured manual touchpoints: A 6 vs baseline 12 (7–8 counting the YAML);
C 7 vs baseline 15 (or 8).

Tests use synthetic Gaussian snippets for the small-imaginary count,
the stalled IRC (no minimum string), endpoint nimag≠0, and a short
invented IRC (`tests/data/tscheck/synthetic_irc_minimum.log`) that
reaches a PES minimum. The public chemsmart fixture `pd_genecp_ts.log`
is also used. No research logs are committed.

## Deliberate deletion

`ts_auto.csv` (and a matching INGEST note) is **not** part of this MVP.
That is a deliberate deletion, not an access limitation. The mode table is
printed to stdout (acceptance criterion 5). Nothing writes a results CSV.

## Disagreements with the spec

None on intended behavior. Implementation notes:

- Commands are top-level `chemsmart check` and `chemsmart queue`, not
  `chemsmart run …`, so `sub` cannot HPC-submit them.
- UMA is absent from the code, as required.
- `queue` still writes a shell file and never submits: TS with
  `--additional-opt-options maxstep=5`, one IRC command (forward and
  reverse together), endpoint opts batched. QRC is never queued.

## Related work (from pass 1–2; not wrapped)

chemsmart already covers ts, irc, qrc, opt. autodE, pysisyphus, and Sella
were reuse candidates for a later guess stage. That stage is the guess
section below. The packages are still not vendored.

## Guess generation (critic pass 1)

### PROBLEM

`check` and `queue` only look at a finished Gaussian TS. The missing
step is a structure to submit. The calculator that step is built for
is UMA (`facebook/UMA`, task `omol`). Hugging Face access is still
rejected, so the same command has to run today on GFN2-xTB.

### METRIC

Plumbing gate for xTB only, not a speedup. On at least 2 of 3 cases,
the non-driven distances (secondary spec bonds, plus Zn–O and Zn–N
pairs listed under `bonds` or `contacts`) satisfy `|Δd| ≤ 0.15 Å`
versus the DFT TS, and the guess's imaginary mode hits the primary
bond (rank 1 within 3.2 Å). Primary-bond `|Δd|` and the heavy-atom
Kabsch RMSD are printed as information. The baseline is Chi's
hand-built guess, read from paths supplied to `scripts/score_guess.py`
(`TS_GUESS_OUT` and `TS_GUESS_<CASE>_{REACTANT,TS,BASELINE,SPEC}`;
no path is hardcoded). This gate is not a claim that the pipeline is
faster. Pass 2 on xTB scored 1 of 3 and the command is not yet better
than hand-built guesses. The chemistry comparison is deferred to UMA.
Pass 3 checks plumbing only.

### ACCEPTANCE

1. `chemsmart guess -f reactant --spec case.yaml` runs a 1D
   constrained scan of the primary bond to YAML `scan.to`, 10–15
   points, everything else relaxed, and prints the energy profile.
2. A maximum on the first or last converged frame prints
   `FAIL (max at scan edge)` and writes no `guess.xyz`. The run still
   writes `scan.xyz` (every frame), `profile.txt`, and `stdout.txt`.
3. Otherwise the converged scan-maximum frame is `guess.xyz` and
   `guess.gjf`, with the tscheck mode table (primary rank, top pairs,
   Zmax as information) and a `queue` shell file. `--sella` is off
   unless asked. Gaussian is not submitted.
4. `--calc uma` is `FAIRChemCalculator`, task `omol`, charge and spin
   (multiplicity) from the YAML. `HF_TOKEN` comes from the environment
   and is never printed. `--check-access` does not evaluate an energy.
   `--calc xtb` is GFN2 via tblite or xtb-python, charge and uhf from
   the YAML or the CLI. Both live in `chemsmart[mlip]`.
5. The product is optional and is not a scan input. NEB is a hook
   (`neb_fallback`), not an implementation. Tests use synthetic
   molecules. Private geometries are not committed.

### QUESTION

- Two endpoints and a CI-NEB: the real case (Case 4) has no product.
  Product-required input is an assumption. Dropped.
- Ranking guesses by barrier: the scan already defines one frame, the
  energy maximum. A second ranking is not a separate decision. Dropped.
- Sella on every guess: extra Hessian iterations before the scan is
  known to be interior. Off unless `--sella`.
- Driving every spec bond: one coordinate matches the case YAML Chi
  already writes. The scan is the primary bond only.
- Atom remapping: same ordering is required. Fail, do not remap.
- Putting fairchem, Sella, and tblite in core dependencies: they are
  not needed for `check` or `queue`. Optional extra.

### DELETE

- Climbing-image NEB (hook left, code not built). Confirmed by the
  acceptance list: a 1D scan still proposes a guess.
- Product as a required input. Reporting against `--product` remains.
- Barrier ranking of several candidates.
- Sella as a default.
- Vendoring autodE, pysisyphus, or React-OT.
- Early-kill and custom route templates, unchanged from pass 2.

### SIMPLIFY

One command, one bond, one profile. Mode projection calls tscheck
(`rank_pairs`, `zmax_row`, the same table). The shell file calls the
existing queue writer. Structures go through `Molecule.from_filepath`.

### ACCELERATE

`--check-access` answers "can UMA run?" before a scan. The default
path does not call Sella. An analytical Hessian is used when the
calculator exposes `get_hessian`; otherwise a finite difference is
the calculator Hessian the mode test asked for.

### AUTOMATE

The scan, the edge check, the mode table, and writing xyz, gjf, and
the queue file. Still manual: `scan.to` and the bond list, judging
the mode, every submit, and the final OK.

## Guess generation (critic pass 3)

`chemsmart guess` is experimental. The pass-2 xTB plumbing gate was
1 of 3: only one case produced a guess, and the other two stopped on
`FAIL (max at scan edge)`. On the xTB benchmark, the hand-built
guesses beat xTB on all 3 cases. The chemistry gate waits for UMA.
This pass is plumbing.

- `scripts/score_guess.py` reads every path from `TS_GUESS_*` env vars.
- `--sella` stays opt-in and report-only. A converged Sella run writes
  `guess_sella.xyz` and does not replace `guess.xyz` or rerun the scan.
- The mode Hessian projects out translation and rotation, runs only on
  a converged frame, and reports how many spurious imaginary modes
  were removed.
- A frame that hits the relaxation step cap is flagged and left out of
  the maximum. If the highest frame did not converge, the table says
  `low confidence`.
- A failed run still writes `scan.xyz`, `profile.txt`, and `stdout.txt`.
  Calculator noise goes to `guess.log`. Stdout is the table.
- Templates set `scan.to: unset`. Chi must set it. There is no default
  and no second driven coordinate. Zn contacts are scoring-only.

### Related methods (not used)

- autodE (Duarte group) builds reaction profiles with constrained
  scans and NEB: <https://github.com/duartegroup/autodE>
- pysisyphus (Steinmetzer, Grimme, and co-workers) optimizes reaction
  paths, including scans and NEB: <https://github.com/eljost/pysisyphus>
- React-OT generates a TS from a reactant and a product by optimal
  transport (Duan, Liu, Du, et al., Nat. Mach. Intell. 7, 615–626,
  2025, <https://doi.org/10.1038/s42256-025-01010-0>). This command
  has no product, so that model is not applicable here.
- Sella is the optional saddle optimizer (Hermes, Sargsyan, Najm,
  Zádor, J. Chem. Theory Comput. 2019, 15, 6536). ASE's
  climbing-image NEB is the unused fallback hook.
