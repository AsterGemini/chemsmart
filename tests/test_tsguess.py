"""TS guess scan. Synthetic molecules only. No research geometries."""

import importlib.util
import json
import os
import subprocess
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
from chemsmart.analysis.tscheck import (
    EXPECTED_TS_ROUTE,
    canonicalize_route,
    load_case_spec,
)
from chemsmart.analysis.tsguess import (
    GJF_TS_ROUTE,
    heavy_kabsch_rmsd,
    lowest_mode,
    neb_fallback,
    non_driven_pairs,
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


class _TranslatingHessian(_ProfileCalc):
    """Bond curvature plus a spurious imaginary translation."""

    def get_hessian(self, atoms):
        hessian = super().get_hessian(atoms)
        trans = np.zeros(hessian.shape[0])
        trans[0::3] = 1.0
        trans /= np.linalg.norm(trans)
        return hessian - 9.0 * np.outer(trans, trans)


class TestScan:
    def test_edge_maximum_still_writes_frames_and_table(self, tmp_path):
        reactant, _, _ = _reactant(tmp_path / "r.xyz")
        spec = _case(_spec(tmp_path / "case.yaml", target=1.2))
        out = tmp_path / "out"
        text, ok = run_guess(
            reactant,
            spec,
            _factory(monotone=True),
            charge=0,
            uhf=0,
            multiplicity=1,
            outdir=str(out),
            calc_name="fake",
            relax_steps=80,
            fmax=0.05,
        )
        assert ok is False
        assert "FAIL (max at scan edge)" in text
        assert not (out / "guess.xyz").exists()
        assert not (out / "guess.sh").exists()
        assert (out / "scan.xyz").is_file()
        assert (out / "profile.txt").is_file()
        table = (out / "stdout.txt").read_text(encoding="utf-8")
        assert table == text
        assert "FAIL (max at scan edge)" in table
        assert "frame" in (out / "profile.txt").read_text(encoding="utf-8")
        from ase.io import read

        frames = read(out / "scan.xyz", index=":")
        assert len(frames) == 10

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
        assert "spurious imaginary modes removed:" in text
        assert (out / "scan.xyz").is_file()
        assert (out / "stdout.txt").is_file()
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
            sella_steps=1,
        )
        assert ok, text
        assert "Sella did not finish" in text or "Sella refined" in text
        assert (out / "guess.xyz").is_file()
        assert (out / "guess_sella.xyz").is_file()
        assert (out / "guess.sh").is_file()
        assert "guess.xyz" in (out / "guess.sh").read_text(encoding="utf-8")
        assert "guess_sella" not in (out / "guess.sh").read_text(
            encoding="utf-8"
        )

    def test_sella_convergence_keeps_scan_maximum(self, tmp_path, monkeypatch):
        from chemsmart.analysis import tsguess as tsguess_mod

        calls = {"n": 0}
        real_scan = tsguess_mod.constrained_scan

        def _counted(*args, **kwargs):
            calls["n"] += 1
            return real_scan(*args, **kwargs)

        monkeypatch.setattr(tsguess_mod, "constrained_scan", _counted)

        class _Done:
            def __init__(self, atoms, order=1, internal=True, logfile=None):
                self.atoms = atoms
                assert order == 1

            def run(self, fmax, steps):
                self.atoms.positions[0, 0] += 0.35
                return True

        fake = types.ModuleType("sella")
        fake.Sella = _Done
        monkeypatch.setitem(sys.modules, "sella", fake)

        reactant, _, _ = _reactant(tmp_path / "r.xyz")
        spec = _case(_spec(tmp_path / "case.yaml", target=1.2))
        out = tmp_path / "out"
        kwargs = dict(
            charge=0,
            uhf=0,
            multiplicity=1,
            outdir=str(out),
            calc_name="fake",
            relax_steps=40,
            fmax=0.1,
        )
        first, ok = run_guess(
            reactant, spec, _factory(), sella=False, **kwargs
        )
        assert ok, first
        assert calls["n"] == 1
        guess_before = (out / "guess.xyz").read_text(encoding="utf-8")
        second, ok = run_guess(
            reactant, spec, _factory(), sella=True, sella_steps=5, **kwargs
        )
        assert ok, second
        assert calls["n"] == 1
        assert "reused saved scan" in second
        assert "Sella refined" in second
        assert "Sella did not finish" not in second
        assert (out / "guess.xyz").read_text(encoding="utf-8") == guess_before
        sella_text = (out / "guess_sella.xyz").read_text(encoding="utf-8")
        assert sella_text != guess_before

    def test_unconverged_sella_is_not_called_refined(
        self, tmp_path, monkeypatch
    ):
        class _Stopped:
            def __init__(self, atoms, order=1, internal=True, logfile=None):
                self.atoms = atoms

            def run(self, fmax, steps):
                self.atoms.positions[1, 1] += 0.2
                return False

        fake = types.ModuleType("sella")
        fake.Sella = _Stopped
        monkeypatch.setitem(sys.modules, "sella", fake)
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
            sella_steps=2,
        )
        assert ok, text
        assert "Sella did not finish" in text
        assert "Sella refined" not in text
        assert (out / "guess_sella.xyz").is_file()

    def test_step_cap_maximum_is_low_confidence(self, tmp_path, monkeypatch):
        from chemsmart.analysis import tsguess as tsguess_mod

        real_relax = tsguess_mod._relax_fixed_bond

        def _flag_long(atoms, bond, target, factory, fmax, steps):
            clean, energy, converged = real_relax(
                atoms, bond, target, factory, fmax, steps
            )
            if float(target) > 2.05:
                return clean, energy + 50.0, False
            return clean, energy, converged

        monkeypatch.setattr(tsguess_mod, "_relax_fixed_bond", _flag_long)
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
        assert "low confidence" in text
        assert "step cap" in text
        assert "FAIL (max at scan edge)" not in text
        assert (out / "guess.xyz").is_file()
        assert "mode hit:" in text

    def test_step_capped_last_frame_is_still_an_edge(
        self, tmp_path, monkeypatch
    ):
        """Energies 0,1,2,3 with only the last frame step-capped."""
        from chemsmart.analysis import tsguess as tsguess_mod

        def _four(_atoms, _bond, _distances, _factory, _fmax, _steps):
            frames = []
            for index, energy in enumerate((0.0, 1.0, 2.0, 3.0)):
                frames.append(
                    {
                        "index": index + 1,
                        "target": 2.0 - 0.2 * index,
                        "distance": 2.0 - 0.2 * index,
                        "energy": energy,
                        "converged": index < 3,
                        "atoms": _atoms.copy(),
                    }
                )
            return frames

        monkeypatch.setattr(tsguess_mod, "constrained_scan", _four)
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
            relax_steps=5,
            fmax=0.2,
        )
        assert ok is False
        assert "FAIL (max at scan edge)" in text
        assert not (out / "guess.xyz").exists()
        assert (out / "scan.xyz").is_file()

    def test_scan_cache_rejects_a_different_charge_or_calc(
        self, tmp_path, monkeypatch
    ):
        from chemsmart.analysis import tsguess as tsguess_mod

        calls = {"n": 0}

        def _rising(atoms, _bond, distances, _factory, _fmax, _steps):
            calls["n"] += 1
            frames = []
            for index, distance in enumerate(distances):
                frames.append(
                    {
                        "index": index + 1,
                        "target": float(distance),
                        "distance": float(distance),
                        "energy": float(index),
                        "converged": True,
                        "atoms": atoms.copy(),
                    }
                )
            return frames

        monkeypatch.setattr(tsguess_mod, "constrained_scan", _rising)
        reactant, _, _ = _reactant(tmp_path / "r.xyz")
        spec = _case(_spec(tmp_path / "case.yaml", target=1.2))
        out = str(tmp_path / "out")
        common = dict(
            uhf=0,
            multiplicity=1,
            outdir=out,
            relax_steps=5,
            fmax=0.2,
        )
        run_guess(
            reactant, spec, _factory(), charge=0, calc_name="fake", **common
        )
        assert calls["n"] == 1
        text = run_guess(
            reactant, spec, _factory(), charge=2, calc_name="fake", **common
        )[0]
        assert calls["n"] == 2
        assert "reused saved scan" not in text
        text = run_guess(
            reactant, spec, _factory(), charge=2, calc_name="xtb", **common
        )[0]
        assert calls["n"] == 3
        assert "reused saved scan" not in text
        text = run_guess(
            reactant, spec, _factory(), charge=2, calc_name="xtb", **common
        )[0]
        assert calls["n"] == 3
        assert "reused saved scan" in text

    def test_translation_is_projected_out_of_the_mode(self):
        atoms = Atoms(
            symbols=["C", "N", "O"],
            positions=[[0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [1.5, 1.2, 0.1]],
        )
        atoms.calc = _TranslatingHessian()
        _signed, mode, n_removed = lowest_mode(atoms)
        assert n_removed >= 1
        stretch = float(mode[1, 0] - mode[0, 0])
        assert abs(stretch) > 0.2

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

    def test_check_access_uses_hub_toplevel_metadata(self, monkeypatch):
        """hub 2.x exports get_hf_file_metadata from huggingface_hub."""
        token = "hf_SUPER_SECRET_TOKEN"
        monkeypatch.setenv("HF_TOKEN", token)
        _install_fake_fairchem(monkeypatch)
        hub = types.ModuleType("huggingface_hub")
        utils = types.ModuleType("huggingface_hub.utils")

        def hf_hub_url(repo_id, filename):
            return f"https://example.invalid/{repo_id}/{filename}"

        def get_hf_file_metadata(url, token=None):
            return types.SimpleNamespace(size=1)

        hub.hf_hub_url = hf_hub_url
        hub.get_hf_file_metadata = get_hf_file_metadata
        monkeypatch.setitem(sys.modules, "huggingface_hub", hub)
        monkeypatch.setitem(sys.modules, "huggingface_hub.utils", utils)
        lines = check_uma_access()
        text = "\n".join(lines)
        assert "UMA access ok" in text
        assert "huggingface_hub is not installed" not in text
        assert token not in text
        assert not hasattr(utils, "get_hf_file_metadata")

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

        monkeypatch.setenv("HF_TOKEN", "hf_SUPER_SECRET_TOKEN")
        monkeypatch.setattr(
            "chemsmart.cli.guess.make_factory",
            lambda *args, **kwargs: _factory(monotone=True),
        )
        monkeypatch.setattr("subprocess.run", _refuse)
        monkeypatch.setattr("subprocess.Popen", _refuse)
        reactant, _, _ = _reactant(tmp_path / "r.xyz")
        spec = _spec(tmp_path / "case.yaml", target=1.2)
        out = tmp_path / "out"
        result = CliRunner().invoke(
            entry_point,
            [
                "guess",
                "-f",
                reactant,
                "--spec",
                spec,
                "-o",
                str(out),
                "--relax-steps",
                "80",
                "--fmax",
                "0.05",
            ],
        )
        assert result.exit_code == 1, result.output
        assert "FAIL (max at scan edge)" in result.output
        assert "____ _   _ _____" not in result.output
        assert not (out / "guess.xyz").exists()
        assert (out / "scan.xyz").is_file()
        assert (out / "stdout.txt").is_file()
        assert "FAIL (max at scan edge)" in (out / "stdout.txt").read_text(
            encoding="utf-8"
        )
        log = (out / "guess.log").read_text(encoding="utf-8")
        assert "hf_SUPER_SECRET_TOKEN" not in log
        assert "hf_SUPER_SECRET_TOKEN" not in result.output
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
    hub.get_hf_file_metadata = get_hf_file_metadata
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
        assert (tmp_path / "out" / "scan.xyz").is_file()
        assert (tmp_path / "out" / "stdout.txt").is_file()


def test_templates_leave_scan_to_for_chi():
    root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    expected = {
        "case1.yaml": (19, 27),
        "case1_ts2.yaml": (19, 1),
        "case3.yaml": (5, 114),
    }
    for name, contact in expected.items():
        path = os.path.join(
            root, "chemsmart", "settings", "templates", "tscheck", name
        )
        spec = load_case_spec(path)
        text = open(path, encoding="utf-8").read()
        assert spec["scan"]["to"] is None
        assert "Chi must set" in text
        assert "to: unset" in text
        assert "to: 1.5" not in text
        assert "to: 2.4" not in text
        pairs = {tuple(sorted(pair)) for pair in non_driven_pairs(spec)}
        assert tuple(sorted(contact)) in pairs


def test_unset_scan_to_is_rejected(tmp_path):
    reactant, _, _ = _reactant(tmp_path / "r.xyz")
    spec = _case(_spec(tmp_path / "case.yaml"))
    raw = yaml.safe_load(open(tmp_path / "case.yaml", encoding="utf-8"))
    raw["scan"]["to"] = "unset"
    (tmp_path / "case.yaml").write_text(yaml.safe_dump(raw), encoding="utf-8")
    spec = _case(str(tmp_path / "case.yaml"))
    assert spec["scan"]["to"] is None
    with pytest.raises(ValueError, match="Chi must set scan.to"):
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


def _load_scorer():
    path = os.path.join(
        os.path.abspath(os.path.join(os.path.dirname(__file__), "..")),
        "scripts",
        "score_guess.py",
    )
    spec = importlib.util.spec_from_file_location("score_guess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_scorer_reads_paths_from_the_environment(tmp_path):
    source_path = os.path.join(
        os.path.abspath(os.path.join(os.path.dirname(__file__), "..")),
        "scripts",
        "score_guess.py",
    )
    source = open(source_path, encoding="utf-8").read()
    assert "/workspace/tsauto" not in source
    assert "TS_GUESS_OUT" in source
    env = os.environ.copy()
    env.pop("TS_GUESS_OUT", None)
    result = subprocess.run(
        [sys.executable, source_path],
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode != 0
    assert "TS_GUESS_OUT" in result.stderr
    assert "hf_" not in result.stderr
    folder = tmp_path / "A_nosella"
    folder.mkdir()
    (folder / "stdout.txt").write_text(
        "scan maximum: frame 12 at 2.2000 Å\nFAIL (max at scan edge)\n",
        encoding="utf-8",
    )
    parsed = _load_scorer().parse_stdout(str(folder))
    assert parsed["edge_fail"] is True
    assert parsed["scan_max"] == "frame 12 at 2.2000 Å"


def test_sella_row_scores_guess_sella_and_stale_mode_cache_is_dropped(
    tmp_path,
):
    module = _load_scorer()
    folder = tmp_path / "A_sella"
    folder.mkdir()
    _write_xyz(
        folder / "guess.xyz",
        ["C", "N", "O"],
        [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [2.0, 1.3, 0.0]],
    )
    _write_xyz(
        folder / "guess_sella.xyz",
        ["C", "N", "O"],
        [[0.0, 0.0, 0.0], [1.2, 0.0, 0.0], [1.2, 1.3, 0.0]],
    )
    sella_path = module.variant_guess_path(str(folder), "sella")
    plain_path = module.variant_guess_path(str(folder), "nosella")
    assert sella_path.endswith("guess_sella.xyz")
    assert plain_path.endswith("guess.xyz")
    from chemsmart.analysis.tsguess import bond_distance, load_structure

    assert bond_distance(load_structure(sella_path), (1, 2)) == pytest.approx(
        1.2
    )
    assert bond_distance(load_structure(plain_path), (1, 2)) == pytest.approx(
        2.0
    )
    stale = tmp_path / "A_baseline_mode.json"
    stale.write_text(json.dumps({"nimag": 99, "hit": False}), encoding="utf-8")
    assert module.load_cached_baseline_mode(str(stale)) is None
    stale.write_text(
        json.dumps(
            {
                "nimag": 1,
                "hit": True,
                "version": module.baseline_code_version(),
            }
        ),
        encoding="utf-8",
    )
    loaded = module.load_cached_baseline_mode(str(stale))
    assert loaded["nimag"] == 1
    assert loaded["version"] == module.baseline_code_version()


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
