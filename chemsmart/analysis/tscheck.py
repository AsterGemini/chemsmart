"""TS candidate check: mode table, hard rejects, route flags.

This is analysis only. It does not submit Gaussian or HPC jobs and does
not add a new job type. The case YAML (reaction-coordinate bonds and
endpoint criteria) is written by the user.

Projections use Gaussian Cartesian displacements
``proj = (d_a - d_b) · u_ab`` (not mass-weighted). Rankings among
heavy-atom pairs follow the ts-scout pass-2 evidence method.
"""

import re
from dataclasses import dataclass, field

import numpy as np
import yaml

from chemsmart.io.gaussian.output import Gaussian16Output

HEAVY_SKIP = {"H", "D", "T"}
CUTOFF_PRIMARY = 3.2
CUTOFF_ZMAX = 3.6
IRC_QRC_MIN_POINTS = 5
QRC_AMP = 0.5

EXPECTED_CHARGE = 0
EXPECTED_MULTIPLICITY = 1

EXPECTED_TS_ROUTE = "# opt=(ts,calcfc,noeigentest,maxstep=5) freq mn15 def2svp"
EXPECTED_OPT_ROUTE = "# opt freq mn15 def2svp"
EXPECTED_IRC_RE = re.compile(
    r"^# mn15 def2svp irc\(calcfc,recalc=6,"
    r"(forward|reverse),maxpoints=512,maxcycle=128\)$"
)


def canonicalize_route(route):
    if not route:
        return ""
    text = str(route).strip().lower()
    text = text.replace("\n", " ")
    return re.sub(r"\s+", " ", text)


def parse_atom_pair(value):
    """Return a 1-indexed atom pair from [i, j], 'i,j', or 'C8-N29'."""
    if value is None:
        return None
    if isinstance(value, (list, tuple)) and len(value) == 2:
        return int(value[0]), int(value[1])
    text = str(value).strip()
    nums = [int(part) for part in re.findall(r"\d+", text)]
    if len(nums) != 2:
        raise ValueError(f"Bond must be two atom indices: {value!r}")
    if nums[0] == nums[1] or nums[0] < 1 or nums[1] < 1:
        raise ValueError(f"Bond atoms must be distinct and 1-indexed: {value}")
    return nums[0], nums[1]


def pair_label(symbols, atom_i, atom_j):
    return f"{symbols[atom_i - 1]}{atom_i}-{symbols[atom_j - 1]}{atom_j}"


def is_heavy(symbol):
    return str(symbol) not in HEAVY_SKIP


def pair_projection(positions, mode, atom_i, atom_j):
    """(d_a - d_b) · u_ab for 1-indexed atoms a, b."""
    a = atom_i - 1
    b = atom_j - 1
    bond = np.asarray(positions[b], dtype=float) - np.asarray(
        positions[a], dtype=float
    )
    norm = np.linalg.norm(bond)
    if norm < 1e-8:
        return 0.0, 0.0
    u_ab = bond / norm
    delta = np.asarray(mode[a], dtype=float) - np.asarray(mode[b], dtype=float)
    return float(np.dot(delta, u_ab)), float(norm)


def heavy_pairs(symbols):
    idxs = [i + 1 for i, sym in enumerate(symbols) if is_heavy(sym)]
    pairs = []
    for ii, atom_i in enumerate(idxs):
        for atom_j in idxs[ii + 1 :]:
            pairs.append((atom_i, atom_j))
    return pairs


def rank_pairs(symbols, positions, mode, cutoff=None):
    """Pairs sorted by |proj| descending. Optional distance cutoff in Å."""
    rows = []
    for atom_i, atom_j in heavy_pairs(symbols):
        proj, dist = pair_projection(positions, mode, atom_i, atom_j)
        if cutoff is not None and dist > cutoff:
            continue
        rows.append(
            {
                "atoms": (atom_i, atom_j),
                "label": pair_label(symbols, atom_i, atom_j),
                "proj": proj,
                "abs_proj": abs(proj),
                "distance": dist,
            }
        )
    rows.sort(key=lambda row: (-row["abs_proj"], row["distance"]))
    for rank, row in enumerate(rows, start=1):
        row["rank"] = rank
    return rows


def zmax_row(symbols, positions, mode, cutoff=CUTOFF_ZMAX):
    best = None
    for atom_i, atom_j in heavy_pairs(symbols):
        si = symbols[atom_i - 1]
        sj = symbols[atom_j - 1]
        zn_x = (si == "Zn" and sj != "Zn") or (sj == "Zn" and si != "Zn")
        if not zn_x:
            continue
        proj, dist = pair_projection(positions, mode, atom_i, atom_j)
        if dist > cutoff:
            continue
        row = {
            "atoms": (atom_i, atom_j),
            "label": pair_label(symbols, atom_i, atom_j),
            "proj": proj,
            "abs_proj": abs(proj),
            "distance": dist,
        }
        if best is None or row["abs_proj"] > best["abs_proj"]:
            best = row
    return best


def load_case_spec(path):
    with open(path, encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    bonds = []
    for item in data.get("bonds") or []:
        bonds.append(parse_atom_pair(item))
    primary = data.get("primary_bond")
    if primary is None and bonds:
        primary = bonds[0]
    elif primary is not None:
        primary = parse_atom_pair(primary)
    endpoints = data.get("endpoints") or {}
    return {
        "name": data.get("name"),
        "charge": data.get("charge", EXPECTED_CHARGE),
        "multiplicity": data.get("multiplicity", EXPECTED_MULTIPLICITY),
        "primary_bond": primary,
        "bonds": bonds,
        "endpoints": endpoints,
        "project": data.get("project", "zn5"),
        "raw": data,
    }


def _eval_bond_criterion(positions, atoms, op, value):
    atom_i, atom_j = parse_atom_pair(atoms)
    dist = float(
        np.linalg.norm(
            np.asarray(positions[atom_j - 1])
            - np.asarray(positions[atom_i - 1])
        )
    )
    ops = {
        "<": dist < value,
        "<=": dist <= value,
        ">": dist > value,
        ">=": dist >= value,
    }
    if op not in ops:
        raise ValueError(f"Unsupported endpoint op {op!r}")
    return ops[op], dist


def evaluate_endpoints(spec, positions, symbols):
    reports = []
    matched = []
    for name, body in (spec.get("endpoints") or {}).items():
        criteria = []
        ok = True
        for rule in body.get("bonds") or []:
            passed, dist = _eval_bond_criterion(
                positions,
                rule.get("atoms"),
                rule.get("op"),
                float(rule.get("value")),
            )
            atom_i, atom_j = parse_atom_pair(rule.get("atoms"))
            criteria.append(
                {
                    "label": pair_label(symbols, atom_i, atom_j),
                    "op": rule.get("op"),
                    "value": float(rule.get("value")),
                    "distance": dist,
                    "passed": passed,
                }
            )
            ok = ok and passed
        reports.append({"name": name, "passed": ok, "criteria": criteria})
        if ok:
            matched.append(name)
    return reports, matched


def last_freq_block(output):
    """Last harmonic-frequency list and matching modes."""
    freqs = list(output.vibrational_frequencies or [])
    modes = list(output.vibrational_modes or [])
    natoms = output.num_atoms
    usable = []
    for mode in modes:
        arr = np.asarray(mode)
        if arr.ndim == 2 and arr.shape[0] == natoms and arr.shape[1] == 3:
            usable.append(arr)
    n = min(len(freqs), len(usable))
    return freqs[:n], usable[:n]


def count_imaginary(freqs, ignore_threshold=-15.0):
    imag = [freq for freq in freqs if freq < ignore_threshold]
    return len(imag), imag


def first_imaginary_mode(freqs, modes, ignore_threshold=-15.0):
    for freq, mode in zip(freqs, modes):
        if freq < ignore_threshold:
            return freq, mode
    return None, None


def detect_kind(route, freqs=None):
    text = canonicalize_route(route)
    if "irc(" in text or " irc" in f" {text}":
        return "irc"
    if "opt=(ts" in text or "opt=( ts" in text:
        return "ts"
    if re.search(r"\bopt\b", text) and "freq" in text:
        return "opt"
    if re.search(r"\bopt\b", text):
        return "opt"
    nimag, _ = count_imaginary(freqs or [])
    if nimag:
        return "ts"
    return "unknown"


def check_route(kind, route):
    canon = canonicalize_route(route)
    if kind == "ts":
        match = canon == EXPECTED_TS_ROUTE
        expected = EXPECTED_TS_ROUTE
    elif kind == "opt":
        match = canon == EXPECTED_OPT_ROUTE
        expected = EXPECTED_OPT_ROUTE
    elif kind == "irc":
        match = bool(EXPECTED_IRC_RE.match(canon))
        expected = (
            "# mn15 def2svp irc(calcfc,recalc=6,"
            "forward|reverse,maxpoints=512,maxcycle=128)"
        )
    else:
        match = False
        expected = "ts, irc, or opt route"
    return {
        "route": route,
        "canonical": canon,
        "expected": expected,
        "match": match,
    }


def irc_point_count(filename, output=None):
    """Max Gaussian IRC Point Number, else structure count."""
    max_point = 0
    found = False
    with open(filename, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if "Point Number" in line:
                nums = re.findall(r"\d+", line)
                if nums:
                    found = True
                    max_point = max(max_point, int(nums[0]))
            if "of points along the reaction path" in line.lower():
                nums = re.findall(r"\d+", line)
                if nums:
                    found = True
                    max_point = max(max_point, int(nums[0]))
    if found:
        return max_point
    if output is not None and output.all_structures:
        return len(output.all_structures)
    return 0


def irc_reached_minimum(filename):
    with open(filename, encoding="utf-8", errors="replace") as handle:
        text = handle.read()
    return "minimum found" in text.lower()


def qrc_commands(project, ts_file, charge, multiplicity, label=None):
    bits = [
        "chemsmart sub gaussian",
        f"-p {project}",
        f"-f {ts_file}",
        f"-c {charge}",
        f"-m {multiplicity}",
    ]
    if label:
        bits.append(f"-l {label}")
    bits.extend(["qrc", "-a", str(QRC_AMP)])
    return " ".join(bits)


@dataclass
class CheckReport:
    kind: str
    verdict: str
    reasons: list = field(default_factory=list)
    nimag: int = None
    imag_freq: float = None
    charge: int = None
    multiplicity: int = None
    route_flag: bool = False
    force_used: bool = False
    top8: list = field(default_factory=list)
    spec_bonds: list = field(default_factory=list)
    primary_rank_32: int = None
    zmax: dict = None
    endpoints: list = field(default_factory=list)
    irc_points: int = None
    qrc_suggestion: str = None
    qrc_commands: str = None
    lines: list = field(default_factory=list)

    def text(self):
        return "\n".join(self.lines) + ("\n" if self.lines else "")


def _append_table(lines, rows, n=8):
    lines.append("top-8 heavy-atom pair projections (no cutoff):")
    lines.append(f"{'#':>3}  {'pair':<16} {'r/A':>8} {'proj':>10}")
    for row in rows[:n]:
        lines.append(
            f"{row['rank']:3d}  {row['label']:<16} "
            f"{row['distance']:8.3f} {row['proj']:10.4f}"
        )
    if not rows:
        lines.append("  (no heavy-atom pairs)")


def check_ts_log(
    filename,
    spec,
    force=False,
    project=None,
    ts_file=None,
    kind=None,
):
    output = Gaussian16Output(filename=filename)
    freqs, modes = last_freq_block(output)
    route = output.route_string
    kind = kind or detect_kind(route, freqs)
    structure = output.last_structure
    symbols = list(structure.chemical_symbols)
    positions = np.asarray(structure.positions, dtype=float)
    charge = output.charge
    multiplicity = output.multiplicity
    nimag, imag = count_imaginary(freqs)
    imag_freq = imag[0] if imag else None
    freq, mode = first_imaginary_mode(freqs, modes)

    report = CheckReport(
        kind=kind,
        verdict="CANDIDATE: judge mode",
        nimag=nimag,
        imag_freq=imag_freq,
        charge=charge,
        multiplicity=multiplicity,
    )
    lines = report.lines
    lines.append(f"log: {filename}")
    lines.append(f"kind: {kind}")
    lines.append(f"nimag: {nimag}")
    if imag_freq is not None:
        lines.append(f"imaginary frequency: {imag_freq:.4f} cm^-1")
    lines.append(f"charge/mult: {charge} {multiplicity}")

    route_info = check_route(kind, route)
    report.route_flag = not route_info["match"]
    lines.append(f"route: {route_info['canonical']}")
    exp_c = spec.get("charge", EXPECTED_CHARGE)
    exp_m = spec.get("multiplicity", EXPECTED_MULTIPLICITY)
    cm_ok = charge == exp_c and multiplicity == exp_m
    if report.route_flag:
        lines.append(f"ROUTE FLAG: expected {route_info['expected']}")
    if not cm_ok:
        report.route_flag = True
        lines.append(
            f"CHARGE/MULT FLAG: expected {exp_c}/{exp_m}, "
            f"found {charge}/{multiplicity}"
        )

    if kind == "irc":
        return _fill_irc(
            report,
            filename,
            output,
            spec,
            structure,
            symbols,
            positions,
            project=project or spec.get("project", "zn5"),
            ts_file=ts_file,
            charge=charge if charge is not None else exp_c,
            multiplicity=(multiplicity if multiplicity is not None else exp_m),
        )

    if kind == "opt":
        return _fill_opt(report, spec, symbols, positions, nimag)

    if nimag != 1:
        report.verdict = "REJECT"
        report.reasons.append(f"nimag≠1 is a hard reject (found {nimag})")
        lines.append(f"verdict: {report.verdict}")
        lines.append(report.reasons[-1])
        return report

    if mode is None:
        report.verdict = "REJECT"
        report.reasons.append("imaginary mode displacements were not parsed")
        lines.append(f"verdict: {report.verdict}")
        lines.append(report.reasons[-1])
        return report

    all_rows = rank_pairs(symbols, positions, mode, cutoff=None)
    rows_32 = rank_pairs(symbols, positions, mode, cutoff=CUTOFF_PRIMARY)
    report.top8 = all_rows[:8]
    _append_table(lines, all_rows, n=8)

    zmax = zmax_row(symbols, positions, mode)
    report.zmax = zmax
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
        spec_bonds.insert(0, primary)

    def _rank_of(pair, rows):
        for row in rows:
            atoms = row["atoms"]
            if atoms == pair or atoms == (pair[1], pair[0]):
                return row
        return None

    lines.append("spec-bond ranks:")
    for pair in spec_bonds:
        all_hit = _rank_of(pair, all_rows)
        cut_hit = _rank_of(pair, rows_32)
        label = pair_label(symbols, pair[0], pair[1])
        if all_hit is None:
            lines.append(f"  {label}: not a heavy-atom pair")
            continue
        rank_32 = cut_hit["rank"] if cut_hit else None
        report.spec_bonds.append(
            {
                "label": label,
                "atoms": pair,
                "rank_all": all_hit["rank"],
                "rank_32": rank_32,
                "proj": all_hit["proj"],
                "distance": all_hit["distance"],
            }
        )
        rank_32_txt = str(rank_32) if rank_32 is not None else "n/a"
        lines.append(
            f"  {label}: rank_all={all_hit['rank']} "
            f"rank_3.2Å={rank_32_txt} proj={all_hit['proj']:.4f} "
            f"r={all_hit['distance']:.3f}"
        )

    if primary:
        cut_hit = _rank_of(primary, rows_32)
        report.primary_rank_32 = cut_hit["rank"] if cut_hit else None
        label = pair_label(symbols, primary[0], primary[1])
        rank = report.primary_rank_32
        lines.append(
            f"primary spec bond {label} rank within {CUTOFF_PRIMARY} Å: "
            f"{rank if rank is not None else 'n/a'}"
        )
        if rank != 1:
            msg = (
                f"Primary spec bond not #1 within {CUTOFF_PRIMARY} Å "
                f"(rank {rank if rank is not None else 'n/a'})"
            )
            if force:
                report.force_used = True
                report.reasons.append(msg + " [--force override]")
                lines.append(msg + " [--force override]")
            else:
                report.verdict = "REJECT"
                report.reasons.append(msg)
                lines.append(f"verdict: {report.verdict}")
                lines.append(msg)
                lines.append("use --force to keep as CANDIDATE: judge mode")
                return report

    if report.route_flag:
        report.reasons.append("route or charge/mult mismatch")
    report.verdict = "CANDIDATE: judge mode"
    lines.append(f"verdict: {report.verdict}")
    lines.append("Mode judgment stays manual. This is not a verified TS.")
    return report


def _fill_opt(report, spec, symbols, positions, nimag):
    lines = report.lines
    if nimag != 0:
        report.verdict = "REJECT"
        report.reasons.append(f"endpoint opt requires nimag=0 (found {nimag})")
        lines.append(f"verdict: {report.verdict}")
        lines.append(report.reasons[-1])
        return report
    endpoints, matched = evaluate_endpoints(spec, positions, symbols)
    report.endpoints = endpoints
    for item in endpoints:
        status = "PASS" if item["passed"] else "FAIL"
        lines.append(f"endpoint {item['name']}: {status}")
        for rule in item["criteria"]:
            lines.append(
                f"  {rule['label']} {rule['op']} {rule['value']} "
                f"(r={rule['distance']:.3f} Å) "
                f"{'ok' if rule['passed'] else 'no'}"
            )
    if not endpoints:
        lines.append("no endpoint bond criteria in spec YAML")
        report.verdict = "CANDIDATE: judge mode"
    elif matched:
        report.verdict = "CANDIDATE: judge mode"
        lines.append(f"matched endpoints: {', '.join(matched)}")
    else:
        report.verdict = "REJECT"
        report.reasons.append("YAML endpoint bond criteria not met")
        lines.append(f"verdict: {report.verdict}")
        lines.append(report.reasons[-1])
        return report
    lines.append(f"verdict: {report.verdict}")
    return report


def _fill_irc(
    report,
    filename,
    output,
    spec,
    structure,
    symbols,
    positions,
    project,
    ts_file,
    charge,
    multiplicity,
):
    lines = report.lines
    npoints = irc_point_count(filename, output=output)
    report.irc_points = npoints
    lines.append(f"IRC points: {npoints}")
    endpoints, matched = evaluate_endpoints(spec, positions, symbols)
    report.endpoints = endpoints
    for item in endpoints:
        status = "PASS" if item["passed"] else "FAIL"
        lines.append(f"IRC end vs {item['name']}: {status}")
        for rule in item["criteria"]:
            lines.append(
                f"  {rule['label']} {rule['op']} {rule['value']} "
                f"(r={rule['distance']:.3f} Å)"
            )

    reached = npoints >= IRC_QRC_MIN_POINTS or irc_reached_minimum(filename)
    ts_for_qrc = ts_file or filename
    if npoints < IRC_QRC_MIN_POINTS:
        report.qrc_suggestion = "suggest qrc ±"
        report.qrc_commands = qrc_commands(
            project, ts_for_qrc, charge, multiplicity
        )
        report.verdict = "CANDIDATE: judge mode"
        lines.append("suggest qrc ±")
        lines.append(report.qrc_commands)
        lines.append(
            "QRC is a printed suggestion only; not written to a queue "
            "file and not submitted."
        )
    elif reached and endpoints and not matched:
        report.qrc_suggestion = None
        report.verdict = "CANDIDATE: judge mode"
        report.reasons.append("connectivity mismatch: judge")
        lines.append("connectivity mismatch: judge")
        lines.append("suggest no QRC")
    else:
        report.verdict = "CANDIDATE: judge mode"
        lines.append("no QRC suggestion")
    lines.append(f"verdict: {report.verdict}")
    return report
