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
are reuse candidates for a later MLIP stage only after HF access and a 3-TS
mode-agreement test. They are not used in this MVP.
