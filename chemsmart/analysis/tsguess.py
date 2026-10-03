"""1D constrained scan that proposes a TS guess.

The scan drives the primary spec bond from the reactant distance to
``scan.to`` in the case YAML, relaxes every other degree of freedom,
and keeps the energy profile. An energy maximum on the first or last
frame is ``FAIL (max at scan edge)`` and writes no guess. Otherwise
the maximum frame is the guess. ``--sella`` optionally refines that
frame (order 1, step cap). The mode table imports the tscheck
projection helpers. Nothing here submits Gaussian.

Climbing-image NEB is a later fallback (``neb_fallback``). It is not
built. autodE, pysisyphus, and React-OT are cited in the docs and are
not vendored.
"""

import os

import numpy as np

from chemsmart.analysis.calculators import CalculatorError
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
            "case YAML needs scan.to, the product-side distance in Å "
            "(for example forming ~1.5 or breaking ~2.4)."
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


def lowest_mode(atoms):
    """Lowest Hessian mode as a unit Cartesian displacement, plus cm⁻¹."""
    from ase.vibrations import VibrationsData

    hessian = cartesian_hessian(atoms)
    data = VibrationsData.from_2d(atoms, hessian)
    signed = _signed_frequencies(data.get_frequencies())
    modes = np.real(np.asarray(data.get_modes(), dtype=complex))
    index = int(np.argmin(signed))
    mode = np.asarray(modes[index], dtype=float)
    norm = float(np.linalg.norm(mode))
    if norm > 0.0:
        mode = mode / norm
    return signed, mode


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
            dyn.run(fmax=fmax, steps=steps)
        except Exception as exc:  # noqa: BLE001 — try Cartesian next
            errors.append(f"internal={internal}: {type(exc).__name__}: {exc}")
            continue
        trial.constraints = []
        trial.calc = None
        return trial
    raise CalculatorError("Sella failed. " + " | ".join(errors))


def _geometry_rmsd(left, right):
    delta = np.asarray(left.positions) - np.asarray(right.positions)
    return float(np.sqrt(np.mean(np.sum(delta * delta, axis=1))))


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
    lines = [f"{'frame':>5}  {'r/A':>8}  {'E/eV':>14}  {'rel/eV':>10}  ok"]
    origin = frames[0]["energy"]
    for frame in frames:
        flag = "yes" if frame["converged"] else "no"
        lines.append(
            f"{frame['index']:5d}  {frame['distance']:8.4f}  "
            f"{frame['energy']:14.6f}  "
            f"{frame['energy'] - origin:10.4f}  {flag}"
        )
    return lines


def _clear_guess_files(outdir):
    for name in (
        "guess.xyz",
        "guess.gjf",
        "scan_max.xyz",
        "scan_max.gjf",
        "guess.sh",
    ):
        path = os.path.join(outdir, name)
        if os.path.isfile(path):
            os.remove(path)


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
    """Scan, optionally refine, write guess files and a queue script.

    Returns ``(text, ok)``. ``ok`` is false for ``FAIL (max at scan edge)``
    and no guess files are left behind. Gaussian is not submitted.
    """
    bond, start, target, count = _scan_setup(spec, points)
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
        "Sella: on (order=1)" if sella else "Sella: off",
        "energy profile (constrained primary bond, other atoms relaxed):",
    ]
    frames = constrained_scan(
        reactant, bond, distances, factory, fmax, relax_steps
    )
    lines.extend(_profile_lines(frames))
    energies = [frame["energy"] for frame in frames]
    imax = int(np.argmax(energies))
    chosen = frames[imax]
    lines.append(
        f"scan maximum: frame {chosen['index']} at "
        f"{chosen['distance']:.4f} Å"
    )
    os.makedirs(outdir, exist_ok=True)
    if imax == 0 or imax == len(frames) - 1:
        lines.append("FAIL (max at scan edge)")
        _clear_guess_files(outdir)
        return "\n".join(lines) + "\n", False

    guess_atoms = chosen["atoms"]
    sella_note = None
    if sella:
        try:
            guess_atoms = _run_sella(
                chosen["atoms"], factory, fmax, sella_steps
            )
            lines.append(
                "Sella refined the scan maximum "
                f"(step cap {int(sella_steps)})."
            )
        except CalculatorError as exc:
            sella_note = str(exc)
            lines.append(f"Sella did not finish: {sella_note}")
            lines.append("Guess stays the scan-maximum frame.")
            guess_atoms = chosen["atoms"]

    differs = _geometry_rmsd(guess_atoms, chosen["atoms"]) > 1e-3
    from chemsmart.cli.queue import _ts_cmd, _write

    project = project or spec.get("project") or "zn5"
    label = label or spec.get("name") or "guess"
    title = f"chemsmart guess {label}"
    guess_xyz = os.path.abspath(os.path.join(outdir, "guess.xyz"))
    guess_gjf = os.path.abspath(os.path.join(outdir, "guess.gjf"))
    write_xyz(guess_xyz, guess_atoms, title)
    write_gjf(guess_gjf, guess_atoms, charge, multiplicity, title)
    commands = [
        _ts_cmd(project, guess_xyz, charge, multiplicity, f"{label}_guess")
    ]
    written = [guess_xyz, guess_gjf]
    if differs:
        scan_xyz = os.path.abspath(os.path.join(outdir, "scan_max.xyz"))
        scan_gjf = os.path.abspath(os.path.join(outdir, "scan_max.gjf"))
        write_xyz(scan_xyz, chosen["atoms"], f"{title} scan maximum")
        write_gjf(
            scan_gjf,
            chosen["atoms"],
            charge,
            multiplicity,
            f"{title} scan maximum",
        )
        commands.append(
            _ts_cmd(
                project,
                scan_xyz,
                charge,
                multiplicity,
                f"{label}_scan_max",
            )
        )
        written.extend([scan_xyz, scan_gjf])
        lines.append(
            "scan-max frame differs from the Sella guess; both written."
        )
    else:
        for name in ("scan_max.xyz", "scan_max.gjf"):
            stale = os.path.join(outdir, name)
            if os.path.isfile(stale):
                os.remove(stale)

    shell = os.path.abspath(os.path.join(outdir, "guess.sh"))
    _write(shell, commands)
    written.append(shell)

    mode_atoms = guess_atoms.copy()
    mode_atoms.calc = factory()
    signed, mode = lowest_mode(mode_atoms)
    lowest = float(signed[int(np.argmin(signed))])
    n_imag = int(np.sum(signed < -IMAG_FLOOR_CM))
    lines.append(f"lowest frequency: {lowest:.2f} cm^-1")
    lines.append(f"nimag (freq < -{IMAG_FLOOR_CM:.0f} cm^-1): {n_imag}")
    summary, info = mode_summary(
        list(mode_atoms.get_chemical_symbols()),
        np.asarray(mode_atoms.positions, dtype=float),
        mode,
        spec,
    )
    lines.extend(summary)
    hit = bool(info["hit"] and lowest < -IMAG_FLOOR_CM)
    lines.append(f"mode hit: {'yes' if hit else 'no'}")
    if product is not None:
        lines.extend(_product_report(guess_atoms, product, spec, bond))
    for path in written:
        lines.append(f"wrote {path}")
    lines.append("Wrote guess.sh (not submitted). Run it yourself.")
    return "\n".join(lines) + "\n", True


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
