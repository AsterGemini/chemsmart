"""Tests for chemsmart check / queue (no new job types, no MLIP)."""

import numpy as np
import pytest
import yaml
from click.testing import CliRunner

from chemsmart.analysis.tscheck import (
    check_ts_log,
    load_case_spec,
    pair_projection,
    rank_pairs,
)
from chemsmart.cli.main import entry_point
from chemsmart.io.gaussian.output import Gaussian16Output


def _write_spec(path, primary, endpoints=None, bonds=None):
    data = {
        "name": "test",
        "charge": 0,
        "multiplicity": 1,
        "primary_bond": list(primary),
        "bonds": [list(p) for p in (bonds or [primary])],
        "endpoints": endpoints or {},
        "project": "zn5",
    }
    path.write_text(yaml.safe_dump(data), encoding="utf-8")
    return str(path)


def _mode_block(freqs, modes, atomic_numbers):
    """One Gaussian 3-column frequency block."""
    natoms = len(atomic_numbers)
    padded = list(freqs) + [20.0, 30.0]
    f1, f2, f3 = padded[:3]
    cols = []
    for k in range(3):
        if k < len(modes):
            cols.append(np.asarray(modes[k], dtype=float))
        else:
            cols.append(np.zeros((natoms, 3)))
    lines = [
        " Harmonic frequencies (cm**-1), IR intensities (KM/Mole), "
        "Raman scattering",
        " activities (A**4/AMU), depolarization ratios for plane and "
        "unpolarized",
        " incident light, reduced masses (AMU), force constants " "(mDyne/A),",
        " and normal coordinates:",
        "                      1                      2"
        "                      3",
        "                      A                      A"
        "                      A",
        f" Frequencies --  {f1:10.4f}            {f2:10.4f}"
        f"            {f3:10.4f}",
        " Red. masses --      1.0000                 1.0000"
        "                 1.0000",
        " Frc consts  --      1.0000                 1.0000"
        "                 1.0000",
        " IR Inten    --      1.0000                 1.0000"
        "                 1.0000",
        "  Atom  AN      X      Y      Z        X      Y      Z"
        "        X      Y      Z",
    ]
    for i, z in enumerate(atomic_numbers):
        a, b, c = cols[0][i], cols[1][i], cols[2][i]
        lines.append(
            f"    {i + 1:2d}  {z:2d}   {a[0]:6.2f} {a[1]:6.2f} "
            f"{a[2]:6.2f}   {b[0]:6.2f} {b[1]:6.2f} {b[2]:6.2f}   "
            f"{c[0]:6.2f} {c[1]:6.2f} {c[2]:6.2f}"
        )
    return "\n".join(lines)


def write_log(
    path,
    *,
    route,
    symbols,
    numbers,
    positions,
    freqs,
    modes,
    charge=0,
    multiplicity=1,
    extra="",
):
    natoms = len(symbols)
    orient = [
        "                         Standard orientation:                         ",
        " ---------------------------------------------------------------------",
        " Center     Atomic      Atomic             Coordinates (Angstroms)",
        " Number     Number       Type             X           Y           Z",
        " ---------------------------------------------------------------------",
    ]
    for i, (z, xyz) in enumerate(zip(numbers, positions), start=1):
        orient.append(
            f"      {i:2d}         {z:3d}           0"
            f"       {xyz[0]:10.6f}   {xyz[1]:10.6f}   {xyz[2]:10.6f}"
        )
    orient.append(
        " ---------------------------------------------------------------------"
    )
    zmat = [
        "Symbolic Z-matrix:",
        f"Charge =  {charge} Multiplicity = {multiplicity}",
    ]
    for sym, xyz in zip(symbols, positions):
        zmat.append(f"{sym}    {xyz[0]:10.6f}  {xyz[1]:10.6f}  {xyz[2]:10.6f}")
    zmat.append("")
    chunks = [
        f" {route}",
        " ---------------------------------------------------------------------",
        f" Charge =  {charge} Multiplicity = {multiplicity}",
        f" NAtoms=   {natoms}",
        *zmat,
        *orient,
    ]
    if extra:
        chunks.append(extra.rstrip())
    if freqs is not None:
        chunks.append(_mode_block(freqs, modes, numbers))
    chunks.append(" Normal termination of Gaussian 16 at test.")
    path.write_text("\n".join(chunks) + "\n", encoding="utf-8")
    return str(path)


SYMBOLS = ["C", "O", "Zn", "H"]
NUMBERS = [6, 8, 30, 1]
POS = np.array(
    [
        [0.00, 0.00, 0.00],
        [1.40, 0.00, 0.00],
        [0.00, 2.00, 0.00],
        [0.00, 0.00, 1.00],
    ]
)
MODE_TS = np.array(
    [
        [0.50, 0.00, 0.00],
        [-0.50, 0.00, 0.00],
        [0.12, 0.00, 0.00],
        [0.00, 0.00, 0.00],
    ]
)
MODE_ROTOR = np.array(
    [
        [0.00, 0.00, 0.00],
        [0.00, 0.00, 0.00],
        [0.00, 0.55, 0.00],
        [0.80, 0.20, 0.00],
    ]
)


@pytest.fixture()
def ts_log(tmp_path):
    return write_log(
        tmp_path / "ts.log",
        route="# opt=(ts,calcfc,noeigentest,maxstep=5) freq MN15 def2svp",
        symbols=SYMBOLS,
        numbers=NUMBERS,
        positions=POS,
        freqs=[-175.16, 20.0, 30.0],
        modes=[MODE_TS],
    )


@pytest.fixture()
def neg1_route_log(tmp_path):
    return write_log(
        tmp_path / "neg1.log",
        route=(
            "#p opt=(ts,calcfc,noeigentest) freq MN15/Def2SVP "
            "int=ultrafine Temperature=383.15"
        ),
        symbols=SYMBOLS,
        numbers=NUMBERS,
        positions=POS,
        freqs=[-219.14, 20.0, 30.0],
        modes=[MODE_TS],
    )


@pytest.fixture()
def rotor_log(tmp_path):
    return write_log(
        tmp_path / "rotor.log",
        route="# opt=(ts,calcfc,noeigentest,maxstep=5) freq MN15 def2svp",
        symbols=SYMBOLS,
        numbers=NUMBERS,
        positions=POS,
        freqs=[-89.82, 20.0, 30.0],
        modes=[MODE_ROTOR],
    )


@pytest.fixture()
def opt_log(tmp_path):
    return write_log(
        tmp_path / "opt.log",
        route="# opt freq MN15 def2svp",
        symbols=SYMBOLS,
        numbers=NUMBERS,
        positions=np.array(
            [
                [0.00, 0.00, 0.00],
                [1.20, 0.00, 0.00],
                [0.00, 2.00, 0.00],
                [0.00, 0.00, 1.00],
            ]
        ),
        freqs=[20.0, 30.0, 40.0],
        modes=[np.zeros((4, 3))],
    )


@pytest.fixture()
def irc_one_point(tmp_path):
    return write_log(
        tmp_path / "irc1.log",
        route=(
            "# MN15 def2svp irc(calcfc,recalc=6,forward,"
            "maxpoints=512,maxcycle=128)"
        ),
        symbols=SYMBOLS,
        numbers=NUMBERS,
        positions=POS,
        freqs=None,
        modes=None,
        extra="Point Number 1 in 1st direction.\n",
    )


@pytest.fixture()
def irc_long(tmp_path):
    extra = "\n".join(
        f"Point Number {i} in 1st direction." for i in range(1, 9)
    )
    extra += "\n Minimum found on this side of the path.\n"
    mismatch_pos = np.array(
        [
            [0.00, 0.00, 0.00],
            [1.90, 0.00, 0.00],
            [0.00, 2.00, 0.00],
            [0.00, 0.00, 1.00],
        ]
    )
    return write_log(
        tmp_path / "irc8.log",
        route=(
            "# MN15 def2svp irc(calcfc,recalc=6,reverse,"
            "maxpoints=512,maxcycle=128)"
        ),
        symbols=SYMBOLS,
        numbers=NUMBERS,
        positions=mismatch_pos,
        freqs=None,
        modes=None,
        extra=extra + "\n",
    )


@pytest.fixture()
def case_yaml(tmp_path):
    return _write_spec(
        tmp_path / "case.yaml",
        (1, 2),
        endpoints={
            "INT2": {"bonds": [{"atoms": [1, 2], "op": "<", "value": 1.5}]},
            "product": {"bonds": [{"atoms": [1, 2], "op": ">", "value": 2.3}]},
        },
    )


class TestPairProjection:
    def test_pure_stretch_sign(self):
        proj, dist = pair_projection(POS, MODE_TS, 1, 2)
        assert dist == pytest.approx(1.4)
        assert proj == pytest.approx(1.0)


class TestCheckTS:
    def test_candidate_abc_like(self, ts_log, case_yaml):
        spec = load_case_spec(case_yaml)
        report = check_ts_log(ts_log, spec)
        assert report.nimag == 1
        assert report.verdict == "CANDIDATE: judge mode"
        assert report.primary_rank_32 == 1
        assert report.route_flag is False
        assert report.zmax is not None
        text = report.text()
        assert "CANDIDATE: judge mode" in text
        assert "Zmax" in text
        assert "information only" in text

    def test_neg1_route_flag_still_candidate(self, neg1_route_log, case_yaml):
        spec = load_case_spec(case_yaml)
        report = check_ts_log(neg1_route_log, spec)
        assert report.verdict == "CANDIDATE: judge mode"
        assert report.route_flag is True
        assert "ROUTE FLAG" in report.text()

    def test_neg2_primary_not_first_rejects(self, rotor_log, case_yaml):
        spec = load_case_spec(case_yaml)
        report = check_ts_log(rotor_log, spec)
        assert report.verdict == "REJECT"
        assert report.primary_rank_32 != 1
        assert "Primary spec bond not #1" in report.text()

    def test_force_overrides_primary_rank(self, rotor_log, case_yaml):
        spec = load_case_spec(case_yaml)
        report = check_ts_log(rotor_log, spec, force=True)
        assert report.verdict == "CANDIDATE: judge mode"
        assert report.force_used is True

    def test_nimag_zero_hard_reject(self, opt_log, case_yaml):
        spec = load_case_spec(case_yaml)
        report = check_ts_log(opt_log, spec, kind="ts")
        assert report.nimag == 0
        assert report.verdict == "REJECT"
        assert "nimag" in report.text()


class TestCheckOptAndIRC:
    def test_opt_endpoint_criteria(self, opt_log, case_yaml):
        spec = load_case_spec(case_yaml)
        report = check_ts_log(opt_log, spec, kind="opt")
        assert report.nimag == 0
        assert report.verdict == "CANDIDATE: judge mode"
        assert any(
            e["name"] == "INT2" and e["passed"] for e in report.endpoints
        )

    def test_irc_one_point_suggests_qrc(self, irc_one_point, case_yaml):
        spec = load_case_spec(case_yaml)
        report = check_ts_log(
            irc_one_point, spec, ts_file="ts.log", project="zn5"
        )
        assert report.irc_points == 1
        assert report.qrc_suggestion == "suggest qrc ±"
        text = report.text()
        assert "suggest qrc ±" in text
        assert "qrc -a 0.5" in text

    def test_irc_minima_connectivity_mismatch_no_qrc(
        self, irc_long, case_yaml
    ):
        spec = load_case_spec(case_yaml)
        report = check_ts_log(irc_long, spec)
        assert report.irc_points >= 5
        assert report.qrc_suggestion is None
        assert "connectivity mismatch: judge" in report.text()
        assert "suggest no QRC" in report.text()


class TestCLI:
    def test_entry_point_lists_check_and_queue(self):
        result = CliRunner().invoke(entry_point, ["--help"])
        assert result.exit_code == 0, result.output
        assert "check" in result.output
        assert "queue" in result.output

    def test_check_cli_candidate(self, ts_log, case_yaml):
        result = CliRunner().invoke(
            entry_point,
            ["check", ts_log, "--spec", case_yaml],
        )
        assert result.exit_code == 0, result.output
        assert "CANDIDATE: judge mode" in result.output

    def test_check_cli_reject_exit(self, rotor_log, case_yaml):
        result = CliRunner().invoke(
            entry_point,
            ["check", rotor_log, "--spec", case_yaml],
        )
        assert result.exit_code == 1
        assert "REJECT" in result.output

    def test_queue_ts_writes_not_submit(self, tmp_path, ts_log):
        out = tmp_path / "ts.sh"
        result = CliRunner().invoke(
            entry_point,
            [
                "queue",
                "-p",
                "zn5",
                "-c",
                "0",
                "-m",
                "1",
                "-o",
                str(out),
                "ts",
                "-f",
                ts_log,
            ],
        )
        assert result.exit_code == 0, result.output
        assert "not submitted" in result.output.lower()
        text = out.read_text()
        assert "chemsmart sub gaussian" in text
        assert "--additional-opt-options maxstep=5 ts" in text
        assert "Not submitted" in text

    def test_queue_irc_and_opt_batched(self, tmp_path, ts_log, irc_one_point):
        irc_out = tmp_path / "irc.sh"
        result = CliRunner().invoke(
            entry_point,
            ["queue", "-p", "zn5", "-o", str(irc_out), "irc", "-f", ts_log],
        )
        assert result.exit_code == 0, result.output
        assert irc_out.read_text().count("chemsmart sub gaussian") == 1
        assert irc_out.read_text().rstrip().endswith("irc")

        opt_out = tmp_path / "opt.sh"
        result = CliRunner().invoke(
            entry_point,
            [
                "queue",
                "-p",
                "zn5",
                "-o",
                str(opt_out),
                "opt",
                "--ircf",
                irc_one_point,
                "--ircr",
                irc_one_point,
            ],
        )
        assert result.exit_code == 0, result.output
        text = opt_out.read_text()
        assert text.count("chemsmart sub gaussian") == 2
        assert text.count(" opt") == 2


class TestRealChemsmartLog:
    """Public chemsmart fixture. ts-scout Case 1/3 logs were not readable."""

    def test_pd_ts_is_candidate_with_route_flag(
        self, gaussian_ts_genecp_outfile, tmp_path
    ):
        output = Gaussian16Output(filename=gaussian_ts_genecp_outfile)
        mode = output.vibrational_modes[0]
        symbols = list(output.last_structure.chemical_symbols)
        pos = np.asarray(output.last_structure.positions)
        rows = rank_pairs(symbols, pos, mode, cutoff=3.2)
        assert rows
        primary = rows[0]["atoms"]
        spec_path = _write_spec(tmp_path / "pd.yaml", primary)
        spec = load_case_spec(spec_path)
        report = check_ts_log(gaussian_ts_genecp_outfile, spec)
        assert report.nimag == 1
        assert report.verdict == "CANDIDATE: judge mode"
        assert report.route_flag is True
        assert report.primary_rank_32 == 1

    def test_pd_ts_quiet_bond_rejects(
        self, gaussian_ts_genecp_outfile, tmp_path
    ):
        output = Gaussian16Output(filename=gaussian_ts_genecp_outfile)
        mode = output.vibrational_modes[0]
        symbols = list(output.last_structure.chemical_symbols)
        pos = np.asarray(output.last_structure.positions)
        rows = rank_pairs(symbols, pos, mode, cutoff=3.2)
        quiet = rows[-1]["atoms"]
        spec_path = _write_spec(tmp_path / "quiet.yaml", quiet)
        spec = load_case_spec(spec_path)
        report = check_ts_log(gaussian_ts_genecp_outfile, spec)
        assert report.verdict == "REJECT"
        assert report.primary_rank_32 != 1
