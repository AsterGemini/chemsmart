"""Score ``chemsmart guess`` outputs. Paths come from the environment.

Required:

- ``TS_GUESS_OUT`` — directory of ``<CASE>_nosella`` and ``<CASE>_sella``
- ``TS_GUESS_CASES`` — comma-separated case labels (for example ``A,B,C``)
- ``TS_GUESS_<CASE>_REACTANT``
- ``TS_GUESS_<CASE>_TS``
- ``TS_GUESS_<CASE>_BASELINE``
- ``TS_GUESS_<CASE>_SPEC``

No path is hardcoded. Writes ``$TS_GUESS_OUT/scores.json``.

``python scripts/score_guess.py baseline-mode-only CASE`` caches the
baseline mode Hessian and prints that record.
"""

import json
import os
import re
import sys
import time

import numpy as np

from chemsmart.analysis.calculators import resolve_charge_uhf, xtb_factory
from chemsmart.analysis.tscheck import load_case_spec
from chemsmart.analysis.tsguess import (
    IMAG_FLOOR_CM,
    assert_same_ordering,
    bond_distance,
    heavy_kabsch_rmsd,
    load_structure,
    lowest_mode,
    mode_summary,
    non_driven_pairs,
)


def require_out():
    out = os.environ.get("TS_GUESS_OUT")
    if not out:
        raise SystemExit(
            "TS_GUESS_OUT is unset. Point it at the directory of "
            "<CASE>_nosella and <CASE>_sella runs."
        )
    return out


def env(case, key):
    name = f"TS_GUESS_{case}_{key}"
    if name not in os.environ:
        raise SystemExit(f"{name} is unset.")
    return os.environ[name]


def geom_metrics(atoms, ts, spec):
    primary = spec["primary_bond"]
    pairs = non_driven_pairs(spec)
    deltas = {}
    for atom_i, atom_j in pairs:
        deltas[f"{atom_i}-{atom_j}"] = abs(
            bond_distance(atoms, (atom_i, atom_j))
            - bond_distance(ts, (atom_i, atom_j))
        )
    values = list(deltas.values())
    return {
        "primary_d": bond_distance(atoms, primary),
        "primary_dd": bond_distance(atoms, primary)
        - bond_distance(ts, primary),
        "nd_mean": float(np.mean(values)) if values else None,
        "nd_max": float(np.max(values)) if values else None,
        "nd": deltas,
        "rmsd_heavy": heavy_kabsch_rmsd(ts, atoms),
    }


def baseline_mode(case, out):
    path = os.path.join(out, f"{case}_baseline_mode.json")
    if os.path.isfile(path):
        with open(path, encoding="utf-8") as handle:
            return json.load(handle)
    spec = load_case_spec(env(case, "SPEC"))
    charge, uhf, _mult = resolve_charge_uhf(spec)
    atoms = load_structure(env(case, "BASELINE"))
    atoms.calc = xtb_factory(charge, uhf)()
    started = time.time()
    signed, mode = lowest_mode(atoms)[:2]
    lowest = float(np.min(signed))
    lines, info = mode_summary(
        list(atoms.get_chemical_symbols()),
        atoms.positions.copy(),
        mode,
        spec,
    )
    zmax = info["zmax"]
    result = {
        "lowest_cm": lowest,
        "nimag": int(np.sum(signed < -IMAG_FLOOR_CM)),
        "primary_rank_32": info["primary_rank_32"],
        "hit": bool(info["hit"] and lowest < -IMAG_FLOOR_CM),
        "zmax": (
            None
            if not zmax
            else {
                "label": zmax["label"],
                "abs_proj": float(zmax["abs_proj"]),
            }
        ),
        "lines": lines,
        "hessian_s": time.time() - started,
    }
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=1)
    return result


def _group(pattern, text):
    match = re.search(pattern, text)
    if match is None:
        return None
    return match.group(1)


def parse_stdout(directory):
    stdout_path = os.path.join(directory, "stdout.txt")
    stderr_path = os.path.join(directory, "stderr.txt")
    timing_path = os.path.join(directory, "timing.txt")
    text = ""
    err = ""
    timing = ""
    if os.path.isfile(stdout_path):
        with open(stdout_path, encoding="utf-8") as handle:
            text = handle.read()
    if os.path.isfile(stderr_path):
        with open(stderr_path, encoding="utf-8") as handle:
            err = handle.read()
    if os.path.isfile(timing_path):
        with open(timing_path, encoding="utf-8") as handle:
            timing = handle.read()
    tail = ""
    if err.strip():
        tail = err.strip().splitlines()[-1]
    return {
        "edge_fail": "FAIL (max at scan edge)" in text,
        "scan_max": _group(r"scan maximum: (frame \d+ at [\d.]+ Å)", text),
        "lowest_cm": _group(r"lowest frequency: (-?[\d.]+)", text),
        "nimag": _group(r"nimag \(freq < -1 cm\^-1\): (\d+)", text),
        "rank": _group(r"rank within 3\.2 Å: (\S+)", text),
        "hit": _group(r"mode hit: (\w+)", text) == "yes",
        "zmax": _group(r"Zmax \(Zn–X, r≤[\d.]+ Å\): ([^(]+)", text),
        "sella": _group(
            r"(Sella refined[^\n]*|Sella did not finish[^\n]*)", text
        ),
        "wall_s": _group(r"wall_s=(\d+)", timing),
        "rc": _group(r"rc=(\d+)", timing),
        "stderr_tail": tail,
    }


def main():
    out = require_out()
    os.makedirs(out, exist_ok=True)
    if len(sys.argv) == 3 and sys.argv[1] == "baseline-mode-only":
        print(json.dumps(baseline_mode(sys.argv[2], out), indent=1))
        return
    raw_cases = os.environ.get("TS_GUESS_CASES")
    if not raw_cases:
        raise SystemExit("TS_GUESS_CASES is unset.")
    results = {}
    for case in raw_cases.split(","):
        case = case.strip()
        if not case:
            continue
        spec = load_case_spec(env(case, "SPEC"))
        ts = load_structure(env(case, "TS"))
        rows = {}
        reactant = load_structure(env(case, "REACTANT"))
        assert_same_ordering(ts, reactant, "TS", "reactant")
        rows["reactant"] = geom_metrics(reactant, ts, spec)
        for variant in ("nosella", "sella"):
            directory = os.path.join(out, f"{case}_{variant}")
            info = parse_stdout(directory)
            guess_xyz = os.path.join(directory, "guess.xyz")
            if os.path.isfile(guess_xyz):
                guess = load_structure(guess_xyz)
                assert_same_ordering(ts, guess, "TS", "guess")
                info.update(geom_metrics(guess, ts, spec))
            rows[variant] = info
        baseline = load_structure(env(case, "BASELINE"))
        assert_same_ordering(ts, baseline, "TS", "baseline")
        mode = baseline_mode(case, out)
        rows["baseline"] = {
            **geom_metrics(baseline, ts, spec),
            **{key: value for key, value in mode.items() if key != "lines"},
        }
        results[case] = rows
    destination = os.path.join(out, "scores.json")
    with open(destination, "w", encoding="utf-8") as handle:
        json.dump(results, handle, indent=1, default=str)
    print(json.dumps(results, indent=1, default=str))


if __name__ == "__main__":
    main()
