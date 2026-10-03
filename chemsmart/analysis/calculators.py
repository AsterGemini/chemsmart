"""ASE calculators for TS guess generation.

UMA (fairchem ``FAIRChemCalculator``, task ``omol``) is the calculator
this command is built for. xTB (GFN2 via tblite, else xtb-python) is
plumbing so the scan can be exercised before a Hugging Face token
exists. Neither stack is imported unless that calculator is selected,
and neither is a core dependency. Install with
``pip install 'chemsmart[mlip]'``.
"""

import logging
import os
from urllib.parse import quote

UMA_REPO = "facebook/UMA"
DEFAULT_UMA_MODEL = "uma-s-1p1"
UMA_TASK = "omol"


class CalculatorError(RuntimeError):
    """Missing optional calculator, token, or gated-model access."""


def resolve_charge_uhf(spec, charge=None, uhf=None, multiplicity=None):
    """Charge and unpaired electrons from the case YAML or CLI overrides.

    CLI ``--uhf`` wins, then CLI ``--multiplicity``, then YAML ``uhf``,
    then YAML ``multiplicity - 1``. UMA spin is this multiplicity
    (``uhf + 1``). These cases are charge 0, uhf 0.
    """
    if charge is None:
        charge = int(spec.get("charge", 0))
    else:
        charge = int(charge)
    if uhf is not None:
        uhf = int(uhf)
        mult = uhf + 1
    elif multiplicity is not None:
        mult = int(multiplicity)
        uhf = mult - 1
    elif spec.get("uhf") is not None:
        uhf = int(spec["uhf"])
        mult = uhf + 1
    else:
        mult = int(spec.get("multiplicity", 1))
        uhf = mult - 1
    if uhf < 0:
        raise CalculatorError(
            "uhf (number of unpaired electrons) cannot be negative."
        )
    return charge, uhf, mult


def scrub_secret(text, secret):
    """Remove a token from text. The token is never written back."""
    if not text:
        return ""
    if not secret:
        return text
    cleaned = text.replace(secret, "<redacted>")
    encoded = quote(secret)
    if encoded and encoded != secret:
        cleaned = cleaned.replace(encoded, "<redacted>")
    if secret in cleaned:
        return "Hugging Face access check failed."
    return cleaned


def _silence_hub_logs():
    logging.getLogger("huggingface_hub").setLevel(logging.ERROR)
    logging.getLogger("filelock").setLevel(logging.ERROR)


def _require_token():
    token = os.environ.get("HF_TOKEN")
    if token is None or not str(token).strip():
        raise CalculatorError(
            "HF_TOKEN is unset. Accept the facebook/UMA licence, create "
            "a token, and export HF_TOKEN. The token is read from the "
            "environment only and is not printed."
        )
    return str(token).strip()


def _require_fairchem():
    try:
        from fairchem.core import FAIRChemCalculator, pretrained_mlip
    except ImportError as exc:
        raise CalculatorError(
            "fairchem is not installed, so --calc uma cannot run. "
            "Install the optional extra: pip install 'chemsmart[mlip]'"
        ) from exc
    return FAIRChemCalculator, pretrained_mlip


def _checkpoint_filenames(model):
    return (f"checkpoints/{model}.pt", f"{model}.pt")


def _access_failure(exc, token):
    text = scrub_secret(str(exc), token).lower()
    name = type(exc).__name__.lower()
    gated = (
        "gated" in text
        or "gated" in name
        or "401" in text
        or "403" in text
        or "forbidden" in text
        or "unauthorized" in text
    )
    if gated:
        return (
            "Hugging Face rejected access to facebook/UMA. Accept the "
            "licence at https://huggingface.co/facebook/UMA and export "
            "HF_TOKEN. No energy was evaluated."
        )
    return (
        "Could not confirm access to facebook/UMA "
        f"({type(exc).__name__}). No energy was evaluated."
    )


def check_uma_access(model=DEFAULT_UMA_MODEL):
    """Confirm fairchem and gated facebook/UMA access. No energy call.

    The checkpoint is checked with a metadata request, not a download
    and not an ASE calculation. ``HF_TOKEN`` is read from the
    environment and is not included in the returned lines.
    """
    _require_fairchem()
    token = _require_token()
    _silence_hub_logs()
    try:
        # hub 2.x exports this from the top-level package. It is not
        # in huggingface_hub.utils, so that import falsely looks missing.
        from huggingface_hub import get_hf_file_metadata, hf_hub_url
    except ImportError as exc:
        raise CalculatorError(
            "huggingface_hub is not installed. It is pulled in by "
            "pip install 'chemsmart[mlip]'."
        ) from exc

    last_error = None
    seen = None
    for filename in _checkpoint_filenames(model):
        url = hf_hub_url(repo_id=UMA_REPO, filename=filename)
        try:
            seen = get_hf_file_metadata(url, token=token)
            break
        except Exception as exc:  # noqa: BLE001 — hub errors vary by version
            last_error = exc
            status = getattr(
                getattr(exc, "response", None), "status_code", None
            )
            if status == 404 or "404" in scrub_secret(str(exc), token):
                continue
            raise CalculatorError(_access_failure(exc, token)) from exc
    if seen is None:
        detail = _access_failure(last_error, token) if last_error else ""
        raise CalculatorError(
            f"Checkpoint for model {model} was not readable on {UMA_REPO}. "
            + detail
        )
    return [
        "UMA access ok",
        f"repo: {UMA_REPO}",
        f"model: {model}",
        f"task: {UMA_TASK}",
        "HF_TOKEN: set",
        "fairchem: installed",
        "gated checkpoint: readable",
        "no energy was evaluated",
    ]


def _uma_calculator_class(predictor, charge, spin):
    from ase.calculators.calculator import Calculator
    from fairchem.core import FAIRChemCalculator

    class _UMA(Calculator):
        """Set omol charge and spin, then call FAIRChemCalculator."""

        implemented_properties = ["energy", "forces"]

        def __init__(self):
            super().__init__()
            self._inner = FAIRChemCalculator(predictor, task_name=UMA_TASK)
            self._charge = int(charge)
            self._spin = int(spin)

        def calculate(self, atoms, properties, system_changes):
            atoms.info["charge"] = self._charge
            atoms.info["spin"] = self._spin
            self._inner.calculate(atoms, properties, system_changes)
            self.results = dict(self._inner.results)

    return _UMA


def uma_factory(charge, multiplicity, model=DEFAULT_UMA_MODEL, device="cpu"):
    """Factory of omol FAIRChemCalculator wrappers sharing one predict unit.

    ``multiplicity`` is written to ``atoms.info['spin']`` (UMA's spin
    is the multiplicity). The predict unit is loaded once.
    """
    _, pretrained_mlip = _require_fairchem()
    token = _require_token()
    _silence_hub_logs()
    try:
        predictor = pretrained_mlip.get_predict_unit(model, device=device)
    except Exception as exc:  # noqa: BLE001 — gated, 401, missing weights
        raise CalculatorError(_access_failure(exc, token)) from exc
    cls = _uma_calculator_class(predictor, charge, multiplicity)
    return cls


def _xtb_backend():
    try:
        from tblite.ase import TBLite  # noqa: F401

        return "tblite"
    except ImportError:
        pass
    try:
        from xtb.ase.calculator import XTB  # noqa: F401

        return "xtb-python"
    except ImportError:
        return None


def xtb_factory(charge, uhf):
    """GFN2-xTB via tblite (charge, multiplicity) or xtb-python (uhf)."""
    backend = _xtb_backend()
    if backend is None:
        raise CalculatorError(
            "xTB calculator needs tblite or xtb-python. "
            "Install the optional extra: pip install 'chemsmart[mlip]'"
        )
    charge = int(charge)
    uhf = int(uhf)
    multiplicity = uhf + 1

    def factory():
        if backend == "tblite":
            from tblite.ase import TBLite

            return TBLite(
                method="GFN2-xTB",
                charge=charge,
                multiplicity=multiplicity,
                verbosity=0,
            )
        from xtb.ase.calculator import XTB

        return XTB(method="GFN2-xTB", charge=charge, uhf=uhf)

    return factory


def check_xtb_import():
    """Report whether the xTB plumbing imports. No energy call."""
    backend = _xtb_backend()
    if backend is None:
        raise CalculatorError(
            "xTB calculator needs tblite or xtb-python. "
            "Install the optional extra: pip install 'chemsmart[mlip]'"
        )
    return [
        f"xTB import ok ({backend}, GFN2-xTB)",
        "xTB does not use Hugging Face.",
        "no energy was evaluated",
    ]


def make_factory(kind, charge, uhf, multiplicity, model, device):
    """Return a zero-arg ASE calculator factory for ``xtb`` or ``uma``."""
    if kind == "xtb":
        return xtb_factory(charge, uhf)
    if kind == "uma":
        return uma_factory(charge, multiplicity, model=model, device=device)
    raise CalculatorError(f"Unknown calculator {kind!r}. Use xtb or uma.")
