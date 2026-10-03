"""1D constrained scan that proposes a TS guess.

Experimental. The pass-2 xTB plumbing gate did not beat hand-built
guesses. The scan drives the primary spec bond from the reactant
distance to ``scan.to`` (Chi sets that number; there is no default),
relaxes every other degree of freedom, and keeps the energy profile.
The maximum is taken only from frames that converged. An energy
maximum on the first or last of those frames is ``FAIL (max at scan
edge)``. A failed run still writes ``scan.xyz``, the profile, and the
result table, and does not write ``guess.xyz``. ``--sella`` is
opt-in and report-only: it writes ``guess_sella.xyz`` and does not
replace ``guess.xyz`` or rerun the scan. The mode table imports the
tscheck projection helpers and is run only on a converged frame, after
translation and rotation are projected out. Nothing here submits
Gaussian.

Climbing-image NEB is a later fallback (``neb_fallback``). It is not
built. autodE, pysisyphus, and React-OT are cited in the docs and are
not vendored.
"""

import hashlib
import json
import logging
import os
import sys

import numpy as np

from chemsmart.analysis.calculators import CalculatorError, scrub_secret
from chemsmart.analysis.tscheck import (
    CUTOFF_PRIMARY,
    CUTOFF_ZMAX,
    SPEC_PROJ_ZERO,
    _append_table,
    pair_label,
    rank_pairs,
    zmax_row,
)

# Same route `check` expects for a TS log. Casing is cosmetic; the
# canonical form is EXPECTED_TS_ROUTE.
GJF_TS_ROUTE = "# opt=(ts,calcfc,noeigentest,maxstep=5) freq MN15 def2svp"
SCAN_POINTS_MIN = 10
SCAN_POINTS_MAX = 15
HEAVY_SKIP = {"H", "D", "T"}
IMAG_FLOOR_CM = 1.0


def neb_fallback(*_args, **_kwargs):
    """Hook for a later climbing-image NEB. Not implemented.

    ASE's climbing-image NEB (and the NEB stages in autodE and
    pysisyphus) can be called from here if a 1D scan is the wrong
    coordinate. This MVP does not build that path.
    """
    raise NotImplementedError(
        "Climbing-image NEB is not in this MVP. "
        "chemsmart guess drives a 1D constrained scan of the primary bond."
    )


def non_driven_pairs(spec):
    """Secondary spec bonds plus YAML ``contacts``, excluding the driven bond.

    The driven bond is ``primary_bond``. Pairs are returned 1-indexed
    and de-duplicated. The benchmark scores |Δd| on this list.
    """
    driven = spec.get("primary_bond")
    driven_key = tuple(sorted(driven)) if driven else None
    pairs = []
    seen = set()
    for source in (spec.get("bonds") or [], spec.get("contacts") or []):
        for pair in source:
            key = tuple(sorted(pair))
            if key == driven_key or key in seen or key[0] == key[1]:
                continue
            seen.add(key)
            pairs.append(tuple(pair))
    return pairs


def load_structure(path):
    """Last geometry in an xyz, gjf/com, or Gaussian log (IRC endpoint)."""
    from chemsmart.io.molecules.structure import Molecule

    molecule = Molecule.from_filepath(path)
    if isinstance(molecule, list):
        molecule = molecule[-1] if molecule else None
    if molecule is None:
        raise ValueError(f"No structure in {path}")
    converted = molecule.to_ase()
    # Plain ASE Atoms. The project subclass does not round-trip through
    # Atoms.copy(), which the scan and the optimizers call.
    from ase import Atoms

    atoms = Atoms(
        symbols=list(converted.get_chemical_symbols()),
        positions=np.asarray(converted.positions, dtype=float),
    )
    return atoms


def assert_same_ordering(left, right, left_name, right_name):
    """Fail if element sequences differ. Guess does not remap atoms."""
    symbols_left = list(left.get_chemical_symbols())
    symbols_right = list(right.get_chemical_symbols())
    if len(symbols_left) != len(symbols_right):
        raise ValueError(
            f"Atom counts differ: {left_name} has {len(symbols_left)}, "
            f"{right_name} has {len(symbols_right)}. "
            "Same atom ordering is required; guess does not remap atoms."
        )
    mismatches = [
        (index + 1, symbols_left[index], symbols_right[index])
        for index in range(len(symbols_left))
        if symbols_left[index] != symbols_right[index]
    ]
    if mismatches:
        shown = ", ".join(
            f"#{index} {old}->{new}" for index, old, new in mismatches[:8]
        )
        raise ValueError(
            f"Atom ordering differs at {len(mismatches)} site(s) "
            f"({shown}). Same atom ordering is required; "
            "guess does not remap atoms."
        )


def bond_distance(atoms, pair):
    """Distance in Å for a 1-indexed pair."""
    atom_i, atom_j = pair
    return float(
        np.linalg.norm(
            atoms.positions[atom_j - 1] - atoms.positions[atom_i - 1]
        )
    )


def _set_bond_length(atoms, pair, target):
    atom_i, atom_j = pair[0] - 1, pair[1] - 1
    vector = atoms.positions[atom_j] - atoms.positions[atom_i]
    length = float(np.linalg.norm(vector))
    if length < 1e-8:
        raise ValueError(f"Scan bond atoms {pair[0]} and {pair[1]} coincide.")
    unit = vector / length
    delta = float(target) - length
    atoms.positions[atom_i] = atoms.positions[atom_i] - 0.5 * delta * unit
    atoms.positions[atom_j] = atoms.positions[atom_j] + 0.5 * delta * unit


def _scan_setup(spec, points_override):
    scan = spec.get("scan") or {}
    primary = spec.get("primary_bond")
    if primary is None:
        raise ValueError("case YAML needs primary_bond for the scan.")
    bond = scan.get("bond") or primary
    if tuple(sorted(bond)) != tuple(sorted(primary)):
        raise ValueError(
            "The scan is 1D on the primary spec bond only. "
            f"scan.bond {list(bond)} is not primary_bond {list(primary)}."
        )
    if scan.get("to") is None:
        raise ValueError(
            "Chi must set scan.to; there is no default. "
            "scan.to is the product-side distance in Å for this case."
        )
    if points_override is None:
        count = int(scan.get("points") or 12)
    else:
        count = int(points_override)
    if count < SCAN_POINTS_MIN or count > SCAN_POINTS_MAX:
        raise ValueError(
            f"The scan uses {SCAN_POINTS_MIN}–{SCAN_POINTS_MAX} points "
            f"(got {count})."
        )
    atom_i, atom_j = bond
    if atom_i < 1 or atom_j < 1 or atom_i == atom_j:
        raise ValueError(f"Scan bond is not a 1-indexed pair: {bond}")
    return bond, scan.get("from", "auto"), float(scan["to"]), count


def _start_distance(atoms, bond, start):
    measured = bond_distance(atoms, bond)
    if isinstance(start, str) and start.strip().lower() == "auto":
        distance = measured
    else:
        distance = float(start)
    return distance, measured


def _relax_fixed_bond(atoms, bond, target, factory, fmax, steps):
    from ase.constraints import FixBondLength
    from ase.optimize import BFGS

    work = atoms.copy()
    work.constraints = []
    work.calc = None
    _set_bond_length(work, bond, target)
    work.calc = factory()
    work.set_constraint(FixBondLength(bond[0] - 1, bond[1] - 1))
    # BFGS holds a fixed bond more steadily than FIRE on a stiff scan.
    opt = BFGS(work, logfile=None, maxstep=0.05)
    converged = opt.run(fmax=fmax, steps=steps)
    energy = float(work.get_potential_energy())
    if converged is None:
        forces = work.get_forces()
        converged = float(np.max(np.linalg.norm(forces, axis=1))) <= fmax
    clean = work.copy()
    clean.constraints = []
    clean.calc = None
    return clean, energy, bool(converged)


def constrained_scan(atoms, bond, distances, factory, fmax, steps):
    """Relax every distance with the primary bond fixed. Keep the profile."""
    frames = []
    current = atoms
    for index, distance in enumerate(distances, start=1):
        relaxed, energy, converged = _relax_fixed_bond(
            current, bond, distance, factory, fmax, steps
        )
        actual = bond_distance(relaxed, bond)
        frames.append(
            {
                "index": index,
                "target": float(distance),
                "distance": actual,
                "energy": energy,
                "converged": converged,
                "atoms": relaxed,
            }
        )
        current = relaxed
    return frames


def _as_square_hessian(raw, natoms):
    array = np.asarray(raw, dtype=float)
    dim = 3 * natoms
    if array.shape == (dim, dim):
        return array
    if array.shape == (natoms, 3, natoms, 3):
        return array.reshape(dim, dim)
    return None


def cartesian_hessian(atoms, step=0.002):
    """Cartesian Hessian in eV/Å².

    Uses ``get_hessian`` when the calculator provides it.
    """
    calc = atoms.calc
    if calc is None:
        raise CalculatorError(
            "No calculator is attached; cannot form a Hessian."
        )
    getter = getattr(calc, "get_hessian", None)
    if callable(getter):
        try:
            raw = getter(atoms)
        except TypeError:
            raw = getter()
        square = _as_square_hessian(raw, len(atoms))
        if square is not None:
            return 0.5 * (square + square.T)
    natoms = len(atoms)
    dim = 3 * natoms
    hessian = np.zeros((dim, dim))
    origin = atoms.get_positions().copy()
    for column in range(dim):
        disp = np.zeros(dim)
        disp[column] = step
        atoms.set_positions(origin + disp.reshape(natoms, 3))
        plus = atoms.get_forces().ravel()
        atoms.set_positions(origin - disp.reshape(natoms, 3))
        minus = atoms.get_forces().ravel()
        hessian[:, column] = -(plus - minus) / (2.0 * step)
    atoms.set_positions(origin)
    return 0.5 * (hessian + hessian.T)


def _signed_frequencies(freqs):
    signed = []
    for value in np.asarray(freqs):
        number = complex(value)
        if abs(number.imag) > 1e-6:
            signed.append(-abs(number))
        else:
            signed.append(float(number.real))
    return np.asarray(signed, dtype=float)


def _diagonalize_hessian(atoms, hessian):
    from ase.vibrations import VibrationsData

    data = VibrationsData.from_2d(atoms, hessian)
    signed = _signed_frequencies(data.get_frequencies())
    modes = np.real(np.asarray(data.get_modes(), dtype=complex))
    return signed, modes


def _tr_basis(atoms, tol=1e-8):
    """Mass-weighted translation and rotation vectors. Linear drops one."""
    masses = np.asarray(atoms.get_masses(), dtype=float)
    sqrt_m = np.sqrt(masses)
    natoms = len(atoms)
    dim = 3 * natoms
    com = np.asarray(atoms.get_center_of_mass(), dtype=float)
    pos = np.asarray(atoms.positions, dtype=float) - com
    raw = []
    for axis in range(3):
        vec = np.zeros(dim)
        vec[axis::3] = sqrt_m
        raw.append(vec)
    for axis in range(3):
        unit = np.zeros(3)
        unit[axis] = 1.0
        vec = np.zeros(dim)
        for index in range(natoms):
            vec[3 * index : 3 * index + 3] = sqrt_m[index] * np.cross(
                unit, pos[index]
            )
        raw.append(vec)
    basis = []
    for vec in raw:
        work = np.asarray(vec, dtype=float).copy()
        for kept in basis:
            work = work - np.dot(work, kept) * kept
        norm = float(np.linalg.norm(work))
        if norm > tol:
            basis.append(work / norm)
    return basis


def project_translation_rotation(hessian, atoms):
    """Return a Cartesian Hessian with translation and rotation removed.

    Projection is done in the mass-weighted Hessian. The result is
    transformed back so ``VibrationsData`` can mass-weight it again.
    """
    array = np.asarray(hessian, dtype=float)
    masses = np.asarray(atoms.get_masses(), dtype=float)
    weights = np.repeat(masses**-0.5, 3)
    weighted = array * weights[:, None] * weights[None, :]
    projector = np.eye(weighted.shape[0])
    for vec in _tr_basis(atoms):
        projector = projector - np.outer(vec, vec)
    projected = projector @ weighted @ projector
    inverse = 1.0 / weights
    cartesian = projected * inverse[:, None] * inverse[None, :]
    return 0.5 * (cartesian + cartesian.T)


def lowest_mode(atoms):
    """Lowest mode after translation and rotation are projected out.

    Returns ``(signed_cm, cartesian_mode, n_removed)``. ``n_removed``
    is how many spurious imaginary modes (frequency < -1 cm⁻¹)
    disappeared once translation and rotation were removed. The mode
    is a unit Cartesian displacement.
    """
    hessian = cartesian_hessian(atoms)
    raw_signed, _ = _diagonalize_hessian(atoms, hessian)
    projected = project_translation_rotation(hessian, atoms)
    signed, modes = _diagonalize_hessian(atoms, projected)
    n_raw = int(np.sum(raw_signed < -IMAG_FLOOR_CM))
    n_proj = int(np.sum(signed < -IMAG_FLOOR_CM))
    n_removed = max(0, n_raw - n_proj)
    index = int(np.argmin(signed))
    mode = np.asarray(modes[index], dtype=float)
    norm = float(np.linalg.norm(mode))
    if norm > 0.0:
        mode = mode / norm
    return signed, mode, n_removed


def _rank_of(pair, rows):
    for row in rows:
        atoms = row["atoms"]
        if atoms == tuple(pair) or atoms == (pair[1], pair[0]):
            return row
    return None


def mode_summary(symbols, positions, mode, spec):
    """tscheck projection of one Cartesian mode. Does not reject the guess.

    Returns ``(lines, info)``. ``info['hit']`` is true only when the
    lowest mode is imaginary, a spec bond moves, and the primary bond
    is rank 1 within 3.2 Å. Zmax is information. Judgment stays manual.
    """
    lines = []
    all_rows = rank_pairs(symbols, positions, mode, cutoff=None)
    rows_32 = rank_pairs(symbols, positions, mode, cutoff=CUTOFF_PRIMARY)
    _append_table(lines, all_rows, n=8)
    zmax = zmax_row(symbols, positions, mode)
    if zmax:
        lines.append(
            f"Zmax (Zn–X, r≤{CUTOFF_ZMAX} Å): {zmax['abs_proj']:.4f} "
            f"{zmax['label']} r={zmax['distance']:.3f} Å  "
            "(information only)"
        )
    else:
        lines.append(
            f"Zmax (Zn–X, r≤{CUTOFF_ZMAX} Å): none  (information only)"
        )
    primary = spec.get("primary_bond")
    spec_bonds = list(spec.get("bonds") or [])
    if primary and primary not in spec_bonds:
        spec_bonds = [primary] + spec_bonds
    hits = []
    for pair in spec_bonds:
        row = _rank_of(pair, all_rows)
        if row is None:
            continue
        hits.append((pair, row, _rank_of(pair, rows_32)))
    undisplaced = bool(hits) and all(
        abs(row["proj"]) < SPEC_PROJ_ZERO for _, row, _ in hits
    )
    primary_row = _rank_of(primary, rows_32) if primary else None
    rank = primary_row["rank"] if primary_row else None
    hit = False
    if undisplaced:
        lines.append("spec bonds not displaced in mode")
    elif primary:
        label = pair_label(symbols, primary[0], primary[1])
        lines.append("spec-bond ranks:")
        for pair, all_hit, cut_hit in hits:
            pair_name = pair_label(symbols, pair[0], pair[1])
            rank_32 = cut_hit["rank"] if cut_hit else None
            rank_txt = str(rank_32) if rank_32 is not None else "n/a"
            lines.append(
                f"  {pair_name}: rank_all={all_hit['rank']} "
                f"rank_3.2Å={rank_txt} proj={all_hit['proj']:.4f} "
                f"r={all_hit['distance']:.3f}"
            )
        lines.append(
            f"primary spec bond {label} rank within {CUTOFF_PRIMARY} Å: "
            f"{rank if rank is not None else 'n/a'}"
        )
        hit = rank == 1
        if rank != 1:
            lines.append(
                f"Primary spec bond not #1 within {CUTOFF_PRIMARY} Å "
                f"(rank {rank if rank is not None else 'n/a'})"
            )
    lines.append("Mode judgment stays manual. This is not a verified TS.")
    return lines, {
        "hit": hit,
        "undisplaced": undisplaced,
        "primary_rank_32": rank,
        "zmax": zmax,
    }


def heavy_kabsch_rmsd(reference, mobile):
    """Heavy-atom Kabsch RMSD in Å. Same atom ordering required."""
    assert_same_ordering(reference, mobile, "reference", "mobile")
    symbols = list(reference.get_chemical_symbols())
    indices = [
        index
        for index, symbol in enumerate(symbols)
        if symbol not in HEAVY_SKIP
    ]
    if len(indices) < 3:
        indices = list(range(len(symbols)))
    return _kabsch_rmsd(
        np.asarray(reference.positions, dtype=float)[indices],
        np.asarray(mobile.positions, dtype=float)[indices],
    )


def _kabsch_rmsd(reference, mobile):
    left = np.asarray(reference, dtype=float)
    right = np.asarray(mobile, dtype=float)
    left_c = left - left.mean(axis=0)
    right_c = right - right.mean(axis=0)
    covariance = right_c.T @ left_c
    rotation_v, _, rotation_w = np.linalg.svd(covariance)
    sign = np.sign(np.linalg.det(rotation_v @ rotation_w))
    if sign == 0.0:
        sign = 1.0
    rotation = rotation_v @ np.diag([1.0, 1.0, sign]) @ rotation_w
    aligned = right_c @ rotation
    delta = left_c - aligned
    return float(np.sqrt(np.mean(np.sum(delta * delta, axis=1))))


def _run_sella(atoms, factory, fmax, steps):
    """One Sella order-1 run. Returns ``(atoms, converged)``.

    ``dyn.run`` is honored. A step-cap stop is not convergence.
    Cartesian internals are tried only when the first run raises.
    """
    try:
        from sella import Sella
    except ImportError as exc:
        raise CalculatorError(
            "Sella is not installed. pip install 'chemsmart[mlip]'"
        ) from exc
    errors = []
    for internal in (True, False):
        trial = atoms.copy()
        trial.constraints = []
        trial.calc = factory()
        try:
            dyn = Sella(
                trial,
                order=1,
                internal=internal,
                logfile=None,
            )
            finished = dyn.run(fmax=fmax, steps=steps)
        except Exception as exc:  # noqa: BLE001 — try Cartesian next
            errors.append(f"internal={internal}: {type(exc).__name__}: {exc}")
            continue
        trial.constraints = []
        trial.calc = None
        return trial, bool(finished)
    raise CalculatorError("Sella failed. " + " | ".join(errors))


def write_gjf(path, atoms, charge, multiplicity, title):
    """Gaussian input with the TS route ``check`` already expects."""
    symbols = atoms.get_chemical_symbols()
    lines = [
        GJF_TS_ROUTE,
        "",
        title,
        "",
        f"{int(charge)} {int(multiplicity)}",
    ]
    for symbol, xyz in zip(symbols, atoms.positions):
        lines.append(
            f" {symbol:<2} {xyz[0]:14.8f} {xyz[1]:14.8f} {xyz[2]:14.8f}"
        )
    lines.append("")
    folder = os.path.dirname(os.path.abspath(path))
    if folder:
        os.makedirs(folder, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))


def write_xyz(path, atoms, comment):
    from ase.io import write

    folder = os.path.dirname(os.path.abspath(path))
    if folder:
        os.makedirs(folder, exist_ok=True)
    tagged = atoms.copy()
    tagged.info["comment"] = comment
    write(path, tagged, format="xyz", comment=comment)


def _profile_lines(frames):
    lines = [
        f"{'frame':>5}  {'r/A':>8}  {'E/eV':>14}  {'rel/eV':>10}  relaxed"
    ]
    origin = frames[0]["energy"]
    for frame in frames:
        flag = "yes" if frame["converged"] else "no (step cap)"
        lines.append(
            f"{frame['index']:5d}  {frame['distance']:8.4f}  "
            f"{frame['energy']:14.6f}  "
            f"{frame['energy'] - origin:10.4f}  {flag}"
        )
    return lines


def _clear_guess_files(outdir):
    """Remove a previous guess. Scan frames and the table stay."""
    for name in (
        "guess.xyz",
        "guess.gjf",
        "guess_sella.xyz",
        "guess_sella.gjf",
        "scan_max.xyz",
        "scan_max.gjf",
        "guess.sh",
    ):
        path = os.path.join(outdir, name)
        if os.path.isfile(path):
            os.remove(path)


def write_scan_xyz(path, frames):
    """Every scan frame, including step-cap frames."""
    lines = []
    for frame in frames:
        atoms = frame["atoms"]
        state = "converged" if frame["converged"] else "step-cap"
        lines.append(str(len(atoms)))
        lines.append(
            f"frame {frame['index']} r={frame['distance']:.4f} "
            f"E={frame['energy']:.6f} {state}"
        )
        for symbol, xyz in zip(atoms.get_chemical_symbols(), atoms.positions):
            lines.append(
                f"{symbol:<2} {xyz[0]:14.8f} {xyz[1]:14.8f} {xyz[2]:14.8f}"
            )
    folder = os.path.dirname(os.path.abspath(path))
    if folder:
        os.makedirs(folder, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")


def _reactant_key(atoms):
    symbols = ",".join(atoms.get_chemical_symbols()).encode()
    coords = np.round(np.asarray(atoms.positions, dtype=float), 8).tobytes()
    return hashlib.sha256(symbols + b"\n" + coords).hexdigest()


def _cache_path(outdir):
    return os.path.join(outdir, "scan_cache.json")


def _cache_matches(payload, bond, distances, fmax, steps, reactant):
    if payload.get("bond") != [int(bond[0]), int(bond[1])]:
        return False
    cached = payload.get("distances") or []
    if len(cached) != len(distances):
        return False
    if not np.allclose(cached, np.asarray(distances, dtype=float), atol=1e-8):
        return False
    if abs(float(payload.get("fmax", -1)) - float(fmax)) > 1e-12:
        return False
    if int(payload.get("relax_steps", -1)) != int(steps):
        return False
    return payload.get("reactant") == _reactant_key(reactant)


def _frames_from_cache(payload):
    from ase import Atoms

    frames = []
    for row in payload["frames"]:
        frames.append(
            {
                "index": int(row["index"]),
                "target": float(row["target"]),
                "distance": float(row["distance"]),
                "energy": float(row["energy"]),
                "converged": bool(row["converged"]),
                "atoms": Atoms(
                    symbols=list(row["symbols"]),
                    positions=np.asarray(row["positions"], dtype=float),
                ),
            }
        )
    return frames


def _save_scan_cache(outdir, bond, distances, fmax, steps, reactant, frames):
    payload = {
        "bond": [int(bond[0]), int(bond[1])],
        "distances": [float(value) for value in distances],
        "fmax": float(fmax),
        "relax_steps": int(steps),
        "reactant": _reactant_key(reactant),
        "frames": [
            {
                "index": int(frame["index"]),
                "target": float(frame["target"]),
                "distance": float(frame["distance"]),
                "energy": float(frame["energy"]),
                "converged": bool(frame["converged"]),
                "symbols": list(frame["atoms"].get_chemical_symbols()),
                "positions": np.asarray(
                    frame["atoms"].positions, dtype=float
                ).tolist(),
            }
            for frame in frames
        ],
    }
    path = _cache_path(outdir)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle)


def _load_scan_cache(outdir, bond, distances, fmax, steps, reactant):
    path = _cache_path(outdir)
    if not os.path.isfile(path):
        return None
    with open(path, encoding="utf-8") as handle:
        payload = json.load(handle)
    if not _cache_matches(payload, bond, distances, fmax, steps, reactant):
        return None
    return _frames_from_cache(payload)


def _choose_maximum(frames):
    """Maximum among converged frames.

    Returns ``(chosen, low_confidence, edge)``. ``chosen`` is None when
    every frame hit the step cap. ``edge`` is true when the frame used
    for the decision is the first or last scan frame. ``low_confidence``
    is true when that decision is not the highest frame in the profile,
    or when the highest frame did not converge.
    """
    energies = [frame["energy"] for frame in frames]
    global_index = int(np.argmax(energies))
    global_max = frames[global_index]
    converged = [frame for frame in frames if frame["converged"]]
    if not converged:
        edge = global_index in (0, len(frames) - 1)
        return None, True, edge, global_max
    conv_energy = [
        frame["energy"] if frame["converged"] else -np.inf for frame in frames
    ]
    chosen_index = int(np.argmax(conv_energy))
    chosen = frames[chosen_index]
    edge = chosen_index in (0, len(frames) - 1)
    low_confidence = (not global_max["converged"]) or (
        chosen_index != global_index
    )
    return chosen, low_confidence, edge, global_max


class _SecretFilter(logging.Filter):
    """Drop HF_TOKEN from any log record before it is written."""

    def filter(self, record):
        token = os.environ.get("HF_TOKEN")
        if not token:
            return True
        message = record.getMessage()
        if token not in message and quote_token(token) not in message:
            return True
        record.msg = scrub_secret(message, token)
        record.args = ()
        return True


def quote_token(token):
    from urllib.parse import quote

    return quote(token)


def _route_guess_logs(outdir):
    """Point root-logger stdout noise at ``guess.log`` for this call."""
    root = logging.getLogger()
    removed = []
    for handler in list(root.handlers):
        stream = getattr(handler, "stream", None)
        is_file = isinstance(handler, logging.FileHandler)
        if isinstance(handler, logging.StreamHandler) and not is_file:
            if stream is sys.stdout:
                root.removeHandler(handler)
                removed.append(handler)
    path = os.path.join(outdir, "guess.log")
    handler = logging.FileHandler(path, mode="w", encoding="utf-8")
    handler.setFormatter(
        logging.Formatter(
            "{asctime} - {levelname:6s} - [{name}] {message}",
            style="{",
        )
    )
    handler.addFilter(_SecretFilter())
    root.addHandler(handler)
    return handler, removed


def _restore_guess_logs(state):
    handler, removed = state
    root = logging.getLogger()
    root.removeHandler(handler)
    handler.close()
    for item in removed:
        root.addHandler(item)


def _scrub_log_file(path):
    token = os.environ.get("HF_TOKEN")
    if not token or not os.path.isfile(path):
        return
    with open(path, encoding="utf-8") as handle:
        text = handle.read()
    cleaned = scrub_secret(text, token)
    if cleaned != text:
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(cleaned)


def _write_text(path, text):
    folder = os.path.dirname(os.path.abspath(path))
    if folder:
        os.makedirs(folder, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(text)


def _maximum_line(frame):
    return f"scan maximum: frame {frame['index']} at {frame['distance']:.4f} Å"


def _append_mode(lines, atoms, factory, spec):
    """Mode table for one converged frame. Returns the hit flag."""
    mode_atoms = atoms.copy()
    mode_atoms.calc = factory()
    signed, mode, n_removed = lowest_mode(mode_atoms)
    lowest = float(signed[int(np.argmin(signed))])
    n_imag = int(np.sum(signed < -IMAG_FLOOR_CM))
    lines.append(f"lowest frequency: {lowest:.2f} cm^-1")
    lines.append(f"nimag (freq < -{IMAG_FLOOR_CM:.0f} cm^-1): {n_imag}")
    lines.append(f"spurious imaginary modes removed: {n_removed}")
    summary, info = mode_summary(
        list(mode_atoms.get_chemical_symbols()),
        np.asarray(mode_atoms.positions, dtype=float),
        mode,
        spec,
    )
    lines.extend(summary)
    hit = bool(info["hit"] and lowest < -IMAG_FLOOR_CM)
    lines.append(f"mode hit: {'yes' if hit else 'no'}")
    return hit


def _finish(outdir, lines, ok):
    text = "\n".join(lines) + "\n"
    _write_text(os.path.join(outdir, "stdout.txt"), text)
    return text, ok


def run_guess(
    reactant_path,
    spec,
    factory,
    *,
    charge,
    uhf,
    multiplicity,
    outdir,
    calc_name,
    project=None,
    label=None,
    sella=False,
    sella_steps=40,
    fmax=0.05,
    relax_steps=80,
    points=None,
    product_path=None,
):
    """Scan, optionally report Sella, write the table.

    Returns ``(text, ok)``. ``ok`` is false for ``FAIL (max at scan
    edge)`` and when no converged frame can be a guess. Those runs
    still write ``scan.xyz``, ``profile.txt``, and ``stdout.txt``.
    ``guess.xyz`` is the converged scan maximum. ``--sella`` writes
    ``guess_sella.xyz`` and does not replace it. Gaussian is not
    submitted.
    """
    bond, start, target, count = _scan_setup(spec, points)
    os.makedirs(outdir, exist_ok=True)
    log_state = _route_guess_logs(outdir)
    try:
        text, ok = _run_guess_logged(
            reactant_path,
            spec,
            factory,
            bond=bond,
            start=start,
            target=target,
            count=count,
            charge=charge,
            uhf=uhf,
            multiplicity=multiplicity,
            outdir=outdir,
            calc_name=calc_name,
            project=project,
            label=label,
            sella=sella,
            sella_steps=sella_steps,
            fmax=fmax,
            relax_steps=relax_steps,
            product_path=product_path,
        )
    finally:
        _restore_guess_logs(log_state)
        _scrub_log_file(os.path.join(outdir, "guess.log"))
    return text, ok


def _run_guess_logged(
    reactant_path,
    spec,
    factory,
    *,
    bond,
    start,
    target,
    count,
    charge,
    uhf,
    multiplicity,
    outdir,
    calc_name,
    project,
    label,
    sella,
    sella_steps,
    fmax,
    relax_steps,
    product_path,
):
    reactant = load_structure(reactant_path)
    natoms = len(reactant)
    atom_i, atom_j = bond
    if atom_i > natoms or atom_j > natoms:
        raise ValueError(
            f"Scan bond {list(bond)} is outside the {natoms}-atom structure."
        )
    d_from, measured = _start_distance(reactant, bond, start)
    if abs(d_from - target) < 1e-6:
        raise ValueError(
            f"scan from ({d_from:.4f} Å) and scan.to ({target:.4f} Å) match, "
            "so there is nothing to drive."
        )
    product = None
    if product_path:
        product = load_structure(product_path)
        assert_same_ordering(reactant, product, "reactant", "product")

    symbols = list(reactant.get_chemical_symbols())
    bond_name = pair_label(symbols, bond[0], bond[1])
    distances = np.linspace(d_from, target, count)
    lines = [
        f"calculator: {calc_name} charge={int(charge)} uhf={int(uhf)} "
        f"(UMA spin={int(multiplicity)} when --calc uma)",
        f"scan: primary {bond_name} from {d_from:.4f} Å to {target:.4f} Å, "
        f"{count} points"
        + (
            f" (reactant r={measured:.4f} Å)"
            if abs(measured - d_from) > 1e-4
            else ""
        ),
        "Sella: on (order=1, report only)" if sella else "Sella: off",
        "energy profile (constrained primary bond, other atoms relaxed):",
    ]
    cached = _load_scan_cache(
        outdir, bond, distances, fmax, relax_steps, reactant
    )
    if cached is None:
        frames = constrained_scan(
            reactant, bond, distances, factory, fmax, relax_steps
        )
        _save_scan_cache(
            outdir, bond, distances, fmax, relax_steps, reactant, frames
        )
    else:
        frames = cached
        lines.append("reused saved scan (not rerun)")
    profile = _profile_lines(frames)
    lines.extend(profile)
    _write_text(os.path.join(outdir, "profile.txt"), "\n".join(profile) + "\n")
    scan_path = os.path.abspath(os.path.join(outdir, "scan.xyz"))
    write_scan_xyz(scan_path, frames)
    lines.append(f"wrote {scan_path}")

    chosen, low_confidence, edge, global_max = _choose_maximum(frames)
    if chosen is None:
        lines.append(_maximum_line(global_max))
        lines.append(
            "low confidence: the energy maximum frame did not converge"
        )
        if edge:
            lines.append("FAIL (max at scan edge)")
        else:
            lines.append("no converged scan frame; no guess written")
        _clear_guess_files(outdir)
        return _finish(outdir, lines, False)

    lines.append(_maximum_line(chosen))
    if low_confidence:
        lines.append(
            "low confidence: the energy maximum frame did not converge "
            f"(highest frame {global_max['index']})"
        )
    if edge:
        lines.append("FAIL (max at scan edge)")
        _clear_guess_files(outdir)
        return _finish(outdir, lines, False)

    guess_atoms = chosen["atoms"]
    from chemsmart.cli.queue import _ts_cmd, _write

    project = project or spec.get("project") or "zn5"
    label = label or spec.get("name") or "guess"
    title = f"chemsmart guess {label}"
    guess_xyz = os.path.abspath(os.path.join(outdir, "guess.xyz"))
    guess_gjf = os.path.abspath(os.path.join(outdir, "guess.gjf"))
    write_xyz(guess_xyz, guess_atoms, title)
    write_gjf(guess_gjf, guess_atoms, charge, multiplicity, title)
    for name in ("scan_max.xyz", "scan_max.gjf", "guess_sella.gjf"):
        stale = os.path.join(outdir, name)
        if os.path.isfile(stale):
            os.remove(stale)
    written = [guess_xyz, guess_gjf]
    if sella:
        sella_xyz = os.path.abspath(os.path.join(outdir, "guess_sella.xyz"))
        try:
            sella_atoms, converged = _run_sella(
                chosen["atoms"], factory, fmax, sella_steps
            )
        except CalculatorError as exc:
            lines.append(f"Sella did not finish: {exc}")
            lines.append("guess.xyz stays the scan-maximum frame.")
            if os.path.isfile(sella_xyz):
                os.remove(sella_xyz)
        else:
            write_xyz(sella_xyz, sella_atoms, f"{title} sella report")
            written.append(sella_xyz)
            if converged:
                lines.append(
                    "Sella refined the scan maximum "
                    f"(step cap {int(sella_steps)}). "
                    "Report only; guess.xyz is the scan maximum."
                )
            else:
                lines.append(
                    "Sella did not finish "
                    f"(step cap {int(sella_steps)}). "
                    "guess.xyz stays the scan-maximum frame."
                )
    else:
        stale = os.path.join(outdir, "guess_sella.xyz")
        if os.path.isfile(stale):
            os.remove(stale)

    if chosen["converged"]:
        _append_mode(lines, guess_atoms, factory, spec)
    else:
        lines.append("mode check skipped: frame did not converge")
        lines.append("mode hit: no")
    if product is not None:
        lines.extend(_product_report(guess_atoms, product, spec, bond))
    commands = [
        _ts_cmd(project, guess_xyz, charge, multiplicity, f"{label}_guess")
    ]
    shell = os.path.abspath(os.path.join(outdir, "guess.sh"))
    _write(shell, commands)
    written.append(shell)
    for path in written:
        lines.append(f"wrote {path}")
    lines.append("Wrote guess.sh (not submitted). Run it yourself.")
    return _finish(outdir, lines, True)


def _product_report(guess, product, spec, driven):
    """Optional product is reporting only. It does not drive the scan."""
    lines = ["product (reporting only):"]
    try:
        rmsd = heavy_kabsch_rmsd(product, guess)
    except ValueError as exc:
        lines.append(f"  {exc}")
        return lines
    lines.append(f"  heavy-atom Kabsch RMSD vs product: {rmsd:.4f} Å")
    lines.append(
        f"  primary |Δd| vs product: {_delta(guess, product, driven):.4f} Å"
    )
    for pair in non_driven_pairs(spec):
        lines.append(
            f"  non-driven |Δd| {pair[0]}-{pair[1]}: "
            f"{_delta(guess, product, pair):.4f} Å"
        )
    return lines


def _delta(left, right, pair):
    return abs(bond_distance(left, pair) - bond_distance(right, pair))
