"""TS guess scan. Synthetic molecules only. No research geometries."""

import sys
import types

import numpy as np
import pytest
import yaml
from ase import Atoms
from ase.calculators.calculator import Calculator
from click.testing import CliRunner

from chemsmart.analysis.calculators import (
    check_uma_access,
    resolve_charge_uhf,
    uma_factory,
    xtb_factory,
)
from chemsmart.analysis.tscheck import EXPECTED_TS_ROUTE, canonicalize_route
from chemsmart.analysis.tsguess import (
    GJF_TS_ROUTE,
    heavy_kabsch_rmsd,
    neb_fallback,
    run_guess,
)
from chemsmart.cli.main import entry_point


def _spec(path, target=1.2, points=10, bonds=None):
    data = {
        "name": "synth",
        "project": "zn5",
        "charge": 0,
        "multiplicity": 1,
        "primary_bond": [1, 2],
        "bonds": bonds or [[1, 2], [2, 3]],
        "endpoints": "unset",
        "scan": {
            "bond": [1, 2],
            "from": "auto",
            "to": target,
            "points": points,
        },
    }
    path.write_text(yaml.safe_dump(data), encoding="utf-8")
    return str(path)


def _write_xyz(path, symbols, positions):
    lines = [str(len(symbols)), "synthetic"]
    for symbol, xyz in zip(symbols, positions):
        lines.append(f"{symbol} {xyz[0]:.6f} {xyz[1]:.6f} {xyz[2]:.6f}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return str(path)


def _reactant(path):
    # Public-style triatomic. C–N is long; the scan drives it shorter.
    symbols = ["C", "N", "O"]
    positions = [[0.0, 0.0, 0.0], [2.2, 0.0, 0.0], [2.2, 1.3, 0.0]]
    return _write_xyz(path, symbols, positions), symbols, positions


class _ProfileCalc(Calculator):
    """Bond maximum at 1.6 Å plus springs on the other coordinates."""

    implemented_properties = ["energy", "forces"]

    def __init__(self, monotone=False, **kwargs):
        super().__init__(**kwargs)
        self.monotone = monotone

    def _energy(self, positions):
        pos = np.asarray(positions, dtype=float)
        bond = pos[1] - pos[0]
        length = np.linalg.norm(bond)
        if self.monotone:
            energy = 5.0 * length
        else:
            delta = length - 1.6
            energy = -40.0 * delta**2 + 8.0 * delta**4
        no = pos[2] - pos[1]
        energy += 80.0 * (np.linalg.norm(no) - 1.3) ** 2
        energy += 30.0 * (pos[0, 1] ** 2 + pos[0, 2] ** 2)
        energy += 30.0 * (pos[1, 1] ** 2 + pos[1, 2] ** 2)
        energy += 30.0 * pos[2, 0] ** 2
        return float(energy)

    def calculate(self, atoms, properties, system_changes):
        super().calculate(atoms, properties, system_changes)
        positions = np.array(atoms.positions, dtype=float, copy=True)
        energy = self._energy(positions)
        forces = np.zeros_like(positions)
        step = 1e-5
        for index in range(len(positions)):
            for coord in range(3):
                positions[index, coord] += step
                plus = self._energy(positions)
                positions[index, coord] -= 2 * step
                minus = self._energy(positions)
                positions[index, coord] += step
                forces[index, coord] = -(plus - minus) / (2 * step)
        self.results = {"energy": energy, "forces": forces}

    def get_hessian(self, atoms):
        """Negative curvature along C–N, positive everywhere else."""
        dim = 3 * len(atoms)
        hessian = np.eye(dim) * 5.0
        vector = atoms.positions[1] - atoms.positions[0]
        vector = vector / np.linalg.norm(vector)
        direction = np.zeros(dim)
        direction[0:3] = -vector
        direction[3:6] = vector
        direction /= np.linalg.norm(direction)
        hessian -= 7.0 * np.outer(direction, direction)
        return hessian


def _factory(monotone=False):
    def factory():
        return _ProfileCalc(monotone=monotone)

    return factory


def _case(path):
    from chemsmart.analysis.tscheck import load_case_spec

    return load_case_spec(path)


class TestChargeAndOrdering:
    def test_uhf_from_multiplicity(self):
        charge, uhf, mult = resolve_charge_uhf(
            {"charge": 0, "multiplicity": 1}
        )
        assert (charge, uhf, mult) == (0, 0, 1)

    def test_cli_uhf_overrides(self):
        charge, uhf, mult = resolve_charge_uhf(
            {"charge": 0, "multiplicity": 1}, uhf=2
        )
        assert (charge, uhf, mult) == (0, 2, 3)

    def test_yaml_uhf(self):
        charge, uhf, mult = resolve_charge_uhf(
            {"charge": -1, "multiplicity": 1, "uhf": 0}
        )
        assert (charge, uhf, mult) == (-1, 0, 1)

    def test_product_ordering_must_match(self, tmp_path):
        reactant, _, _ = _reactant(tmp_path / "r.xyz")
        product = _write_xyz(
            tmp_path / "p.xyz",
            ["N", "C", "O"],
            [[0.0, 0.0, 0.0], [1.2, 0.0, 0.0], [1.2, 1.3, 0.0]],
        )
        spec = _case(_spec(tmp_path / "case.yaml"))
        with pytest.raises(ValueError, match="does not remap atoms"):
            run_guess(
                reactant,
                spec,
                _factory(),
                charge=0,
                uhf=0,
                multiplicity=1,
                outdir=str(tmp_path / "out"),
                calc_name="fake",
                product_path=product,
                relax_steps=5,
            )

    def test_points_outside_10_15(self, tmp_path):
        reactant, _, _ = _reactant(tmp_path / "r.xyz")
        spec = _case(_spec(tmp_path / "case.yaml", points=9))
        with pytest.raises(ValueError, match="10"):
            run_guess(
                reactant,
                spec,
                _factory(),
                charge=0,
                uhf=0,
                multiplicity=1,
                outdir=str(tmp_path / "out"),
                calc_name="fake",
            )


class TestScan:
    def test_edge_maximum_writes_nothing(self, tmp_path):
        reactant, _, _ = _reactant(tmp_path / "r.xyz")
        spec = _case(_spec(tmp_path / "case.yaml", target=1.2))
        text, ok = run_guess(
            reactant,
            spec,
            _factory(monotone=True),
            charge=0,
            uhf=0,
            multiplicity=1,
            outdir=str(tmp_path / "out"),
            calc_name="fake",
            relax_steps=15,
            fmax=0.2,
        )
        assert ok is False
        assert "FAIL (max at scan edge)" in text
        assert not (tmp_path / "out" / "guess.xyz").exists()
        assert not (tmp_path / "out" / "guess.sh").exists()

    def test_interior_maximum_writes_guess_and_queue(self, tmp_path):
        reactant, _, _ = _reactant(tmp_path / "r.xyz")
        spec = _case(_spec(tmp_path / "case.yaml", target=1.2))
        out = tmp_path / "out"
        text, ok = run_guess(
            reactant,
            spec,
            _factory(),
            charge=0,
            uhf=0,
            multiplicity=1,
            outdir=str(out),
            calc_name="fake",
            relax_steps=40,
            fmax=0.1,
        )
        assert ok, text
        assert "Sella: off" in text
        assert "mode hit: yes" in text
        assert "primary spec bond C1-N2 rank within 3.2" in text
        assert "Zmax" in text
        assert "information only" in text
        assert "not submitted" in text
        assert (out / "guess.xyz").is_file()
        gjf = (out / "guess.gjf").read_text(encoding="utf-8")
        assert GJF_TS_ROUTE in gjf
        assert canonicalize_route(GJF_TS_ROUTE) == EXPECTED_TS_ROUTE
        shell = (out / "guess.sh").read_text(encoding="utf-8")
        assert "chemsmart sub gaussian" in shell
        assert "--additional-opt-options maxstep=5" in shell
        assert shell.strip().endswith("ts")
        assert "Not submitted" in shell
        assert not (out / "scan_max.xyz").exists()

    def test_sella_optional_still_writes_guess(self, tmp_path):
        pytest.importorskip("sella")
        reactant, _, _ = _reactant(tmp_path / "r.xyz")
        spec = _case(_spec(tmp_path / "case.yaml", target=1.2))
        out = tmp_path / "out"
        text, ok = run_guess(
            reactant,
            spec,
            _factory(),
            charge=0,
            uhf=0,
            multiplicity=1,
            outdir=str(out),
            calc_name="fake",
            relax_steps=40,
            fmax=0.1,
            sella=True,
            sella_steps=3,
        )
        assert ok, text
        assert "Sella:" in text
        assert (out / "guess.xyz").is_file()
        assert (out / "guess.sh").is_file()

    def test_neb_is_only_a_hook(self):
        with pytest.raises(NotImplementedError, match="not in this MVP"):
            neb_fallback()

    def test_kabsch_ignores_a_rotation(self):
        symbols = ["C", "N", "O"]
        positions = np.array(
            [[0.0, 0.0, 0.0], [1.4, 0.0, 0.0], [1.4, 1.2, 0.1]]
        )
        angle = 0.4
        rotation = np.array(
            [
                [np.cos(angle), -np.sin(angle), 0.0],
                [np.sin(angle), np.cos(angle), 0.0],
                [0.0, 0.0, 1.0],
            ]
        )
        moved = positions @ rotation.T + np.array([1.0, -2.0, 0.5])
        left = Atoms(symbols=symbols, positions=positions)
        right = Atoms(symbols=symbols, positions=moved)
        assert heavy_kabsch_rmsd(left, right) == pytest.approx(0.0, abs=1e-8)


class TestCLI:
    def test_help_lists_guess(self):
        result = CliRunner().invoke(entry_point, ["--help"])
        assert result.exit_code == 0, result.output
        assert "guess" in result.output

    def test_check_access_uma_missing_fairchem(self, monkeypatch):
        monkeypatch.setenv("HF_TOKEN", "hf_SUPER_SECRET_TOKEN")
        monkeypatch.setitem(sys.modules, "fairchem", None)
        monkeypatch.setitem(sys.modules, "fairchem.core", None)
        result = CliRunner().invoke(
            entry_point,
            ["guess", "--calc", "uma", "--check-access"],
        )
        assert result.exit_code == 1
        combined = result.output
        assert "chemsmart[mlip]" in combined
        assert "hf_SUPER_SECRET_TOKEN" not in combined

    def test_check_access_requires_token(self, monkeypatch):
        monkeypatch.delenv("HF_TOKEN", raising=False)
        _install_fake_fairchem(monkeypatch)
        with pytest.raises(Exception) as caught:
            check_uma_access()
        assert "HF_TOKEN is unset" in str(caught.value)
        assert "hf_" not in str(caught.value)

    def test_token_is_not_in_the_access_error(self, monkeypatch):
        token = "hf_SUPER_SECRET_TOKEN"
        monkeypatch.setenv("HF_TOKEN", token)
        _install_fake_fairchem(monkeypatch)
        _install_fake_hub(monkeypatch, fail=token)
        with pytest.raises(Exception) as caught:
            check_uma_access()
        message = str(caught.value)
        assert token not in message
        assert "facebook/UMA" in message

    def test_check_access_ok_does_not_print_token(self, monkeypatch):
        token = "hf_SUPER_SECRET_TOKEN"
        monkeypatch.setenv("HF_TOKEN", token)
        _install_fake_fairchem(monkeypatch)
        _install_fake_hub(monkeypatch, fail=None)
        lines = check_uma_access()
        text = "\n".join(lines)
        assert "UMA access ok" in text
        assert "omol" in text
        assert token not in text
        assert "HF_TOKEN: set" in text

    def test_uma_sets_charge_and_spin(self, monkeypatch):
        monkeypatch.setenv("HF_TOKEN", "hf_SUPER_SECRET_TOKEN")
        seen = {}

        class _Inner(Calculator):
            implemented_properties = ["energy", "forces"]

            def calculate(self, atoms, properties, system_changes):
                seen["charge"] = atoms.info.get("charge")
                seen["spin"] = atoms.info.get("spin")
                self.results = {
                    "energy": -1.0,
                    "forces": np.zeros((len(atoms), 3)),
                }

        def _predict(model, device="cpu"):
            seen["model"] = model
            seen["device"] = device
            return object()

        core = _install_fake_fairchem(monkeypatch)
        core.FAIRChemCalculator = lambda predictor, task_name: _Inner()
        core.pretrained_mlip.get_predict_unit = _predict
        factory = uma_factory(-1, 3, model="uma-s-1p1", device="cpu")
        atoms = Atoms("H", positions=[[0.0, 0.0, 0.0]])
        atoms.calc = factory()
        assert atoms.get_potential_energy() == -1.0
        assert seen["charge"] == -1
        assert seen["spin"] == 3
        assert seen["model"] == "uma-s-1p1"
        assert "hf_SUPER_SECRET_TOKEN" not in repr(seen)

    def test_cli_edge_does_not_submit(self, tmp_path, monkeypatch):
        def _refuse(*_args, **_kwargs):
            raise AssertionError("a job was submitted")

        monkeypatch.setattr(
            "chemsmart.cli.guess.make_factory",
            lambda *args, **kwargs: _factory(monotone=True),
        )
        monkeypatch.setattr("subprocess.run", _refuse)
        monkeypatch.setattr("subprocess.Popen", _refuse)
        reactant, _, _ = _reactant(tmp_path / "r.xyz")
        spec = _spec(tmp_path / "case.yaml", target=1.2)
        result = CliRunner().invoke(
            entry_point,
            [
                "guess",
                "-f",
                reactant,
                "--spec",
                spec,
                "-o",
                str(tmp_path / "out"),
                "--relax-steps",
                "15",
                "--fmax",
                "0.2",
            ],
        )
        assert result.exit_code == 1, result.output
        assert "FAIL (max at scan edge)" in result.output
        assert not (tmp_path / "out" / "guess.xyz").exists()
        assert "hf_" not in result.output


def _install_fake_fairchem(monkeypatch):
    core = types.ModuleType("fairchem.core")

    class FAIRChemCalculator:
        pass

    core.FAIRChemCalculator = FAIRChemCalculator
    core.pretrained_mlip = types.SimpleNamespace(
        get_predict_unit=lambda *args, **kwargs: object()
    )
    fake = types.ModuleType("fairchem")
    fake.core = core
    monkeypatch.setitem(sys.modules, "fairchem", fake)
    monkeypatch.setitem(sys.modules, "fairchem.core", core)
    return core


def _install_fake_hub(monkeypatch, fail):
    hub = types.ModuleType("huggingface_hub")
    utils = types.ModuleType("huggingface_hub.utils")

    def hf_hub_url(repo_id, filename):
        return f"https://example.invalid/{repo_id}/{filename}"

    def get_hf_file_metadata(url, token=None):
        if fail:
            raise RuntimeError(f"401 {token}")
        return types.SimpleNamespace(size=1)

    hub.hf_hub_url = hf_hub_url
    utils.get_hf_file_metadata = get_hf_file_metadata
    monkeypatch.setitem(sys.modules, "huggingface_hub", hub)
    monkeypatch.setitem(sys.modules, "huggingface_hub.utils", utils)


def test_gfn2_scan_on_water(tmp_path):
    """Plumbing: GFN2 can drive one O–H of water. Public geometry."""
    pytest.importorskip("tblite")
    from chemsmart.analysis.calculators import xtb_factory as build_xtb

    xyz = _write_xyz(
        tmp_path / "water.xyz",
        ["O", "H", "H"],
        [[0.000, 0.000, 0.000], [0.757, 0.586, 0.000], [-0.757, 0.586, 0.000]],
    )
    spec = _case(_spec(tmp_path / "case.yaml", target=1.6, points=10))
    text, ok = run_guess(
        xyz,
        spec,
        build_xtb(0, 0),
        charge=0,
        uhf=0,
        multiplicity=1,
        outdir=str(tmp_path / "out"),
        calc_name="xtb",
        relax_steps=25,
        fmax=0.1,
        sella=False,
    )
    assert "energy profile" in text
    assert "charge=0 uhf=0" in text
    if ok:
        assert "not submitted" in text
        assert (tmp_path / "out" / "guess.xyz").is_file()
    else:
        assert "FAIL (max at scan edge)" in text
        assert not (tmp_path / "out" / "guess.xyz").exists()


def test_xtb_factory_passes_charge_and_multiplicity(monkeypatch):
    tblite = pytest.importorskip("tblite")
    captured = {}

    class _Fake:
        def __init__(self, **kwargs):
            captured.update(kwargs)

    monkeypatch.setattr(tblite.ase, "TBLite", _Fake)
    factory = xtb_factory(0, 0)
    factory()
    assert captured["method"] == "GFN2-xTB"
    assert captured["charge"] == 0
    assert captured["multiplicity"] == 1
