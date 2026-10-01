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
- Pipeline (predicted): 6–7 (build guess, write case YAML once per case, run
  TS queue, judge mode from the table, run IRC queue, run endpoint/QRC queue,
  final OK).

Pass 3 was to measure this on A (normal IRC) and C (QRC path) with real code.
Those ts-scout logs were **not readable** from this environment
(`AsterGemini/ts-scout` GitHub 404). The metric is not re-counted here.

## ACCEPTANCE CRITERIA (from PASS 2 FINAL / user)

1. `check` prints nimag, imag freq, charge/mult, top-8 pair projections with
   distances, spec-bond ranks, Zmax (info only).
2. nimag≠1 hard reject. Primary spec bond not #1 within 3.2 Å is a reject
   `--force` can override; rank always printed. Else CANDIDATE: judge mode.
3. Route and 0/1 charge/mult checked against the three expected MN15/def2svp
   routes. Neg1-style routes are flagged.
4. Endpoint-opt: nimag=0 and YAML endpoint bonds.
5. `queue` writes, never submits, ts → irc fwd+rev batched → endpoint opt
   batched. QRC is printed only: <5 IRC points → `suggest qrc ±` plus
   commands (amp 0.5, existing qrc job). Minima but failed bonds →
   `connectivity mismatch: judge`, no QRC.
6. No new job types. No MLIP code. No early-kill. No custom Gaussian
   templates. Do not edit ts-scout `calc/chemsmart/zn5.yaml`.
7. Tests: ts-scout logs if accessible, else chemsmart test data (and say so).
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
| Hard reject nimag≠1 | Robust on the benchmark set |
| Reject primary bond not #1 within 3.2 Å (`--force` override) | Catches the rotor (Neg2); rank always printed |
| Route and 0/1 charge/mult flag | Catches Neg1's old route without auto-rejecting |
| Endpoint nimag=0 + YAML bond criteria | Survived pass 1–2 |
| Write (not submit) ts / irc / opt command files | File shuffling only |
| Print QRC suggestion when IRC points < 5 | Diagnosis only; not queued |
| Print connectivity mismatch, no QRC | Critic fix 1 |

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

`AsterGemini/ts-scout` was not readable (GitHub 404). No A/B/C/Neg1/Neg2 or
Case 3 IRC logs on this box.

Tests use:

- Synthetic Gaussian logs that reproduce the spec outcomes: A/B/C-like
  CANDIDATE, Neg2-like REJECT (primary not #1), Neg1-like CANDIDATE with
  ROUTE FLAG, one-point IRC → `suggest qrc ±`, long IRC with failed endpoint
  bonds → `connectivity mismatch: judge` and no QRC.
- Public chemsmart fixture `pd_genecp_ts.log`: nimag=1 CANDIDATE with a route
  flag (not MN15/def2svp/maxstep=5); quiet-bond spec → REJECT.

Pass-3 metric replay on A and C is **blocked** without ts-scout logs.

## Disagreements with the spec

None on intended behavior. Implementation notes:

- Commands are top-level `chemsmart check` and `chemsmart queue`, not
  `chemsmart run …`, so `sub` cannot HPC-submit them.
- ts-scout CSV `calc/results/ts_auto.csv` / INGEST.md lives in ts-scout and
  was not edited (repo not readable; also out of this fork).
- UMA is absent from the code, as required.

## Related work (from pass 1–2; not wrapped)

chemsmart already covers ts, irc, qrc, opt. autodE, pysisyphus, and Sella
are reuse candidates for a later MLIP stage only after HF access and a 3-TS
mode-agreement test. They are not used in this MVP.
