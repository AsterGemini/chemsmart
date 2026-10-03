"""CLI: propose a TS guess from a reactant and a case YAML.

Writes xyz, gjf, and a ``chemsmart queue`` shell file. Does not submit
Gaussian or any other job.
"""

import sys

import click

from chemsmart.analysis.calculators import (
    CalculatorError,
    check_uma_access,
    check_xtb_import,
    make_factory,
    resolve_charge_uhf,
)
from chemsmart.analysis.tscheck import load_case_spec
from chemsmart.analysis.tsguess import run_guess


@click.command("guess")
@click.option(
    "-f",
    "--filename",
    type=click.Path(exists=True, dir_okay=False),
    default=None,
    help="Reactant or precomplex (xyz, gjf/com, or Gaussian log).",
)
@click.option(
    "--spec",
    type=click.Path(exists=True, dir_okay=False),
    default=None,
    help="Case YAML: primary bond and scan.to (product-side distance).",
)
@click.option(
    "--product",
    type=click.Path(exists=True, dir_okay=False),
    default=None,
    help="Optional. Reporting only. Same atom ordering required.",
)
@click.option(
    "--calc",
    type=click.Choice(["uma", "xtb"], case_sensitive=False),
    default="xtb",
    show_default=True,
    help="uma is the intended calculator. xtb is plumbing.",
)
@click.option(
    "--check-access",
    is_flag=True,
    default=False,
    help="With --calc uma, check fairchem and gated-model access. No job.",
)
@click.option(
    "--model",
    default="uma-s-1p1",
    show_default=True,
    help="UMA checkpoint name on facebook/UMA.",
)
@click.option(
    "--device",
    default="cpu",
    show_default=True,
    help="Device passed to fairchem get_predict_unit.",
)
@click.option(
    "--sella",
    is_flag=True,
    default=False,
    help=(
        "Report-only Sella order=1 run. Writes guess_sella.xyz and "
        "does not replace guess.xyz. Off by default."
    ),
)
@click.option(
    "--sella-steps",
    type=int,
    default=40,
    show_default=True,
    help="Step cap for --sella.",
)
@click.option(
    "--points",
    type=int,
    default=None,
    help="Scan points (10–15). Default: YAML scan.points or 12.",
)
@click.option(
    "--fmax",
    type=float,
    default=0.05,
    show_default=True,
    help="Relaxation force threshold in eV/Å.",
)
@click.option(
    "--relax-steps",
    type=int,
    default=80,
    show_default=True,
    help="Step cap for each constrained relaxation.",
)
@click.option(
    "-c",
    "--charge",
    type=int,
    default=None,
    help="Override the case YAML charge.",
)
@click.option(
    "--uhf",
    type=int,
    default=None,
    help="Unpaired electrons. Overrides multiplicity. 0 for these singlets.",
)
@click.option(
    "-m",
    "--multiplicity",
    type=int,
    default=None,
    help="Override the case YAML multiplicity. uhf = multiplicity - 1.",
)
@click.option(
    "-p",
    "--project",
    type=str,
    default=None,
    help="Project name for the queue file. Default: YAML project.",
)
@click.option(
    "-l",
    "--label",
    type=str,
    default=None,
    help="Queue label. Default: YAML name.",
)
@click.option(
    "-o",
    "--outdir",
    type=click.Path(file_okay=False),
    default="guess_out",
    show_default=True,
    help="Directory for guess.xyz, guess.gjf, and guess.sh.",
)
def guess(
    filename,
    spec,
    product,
    calc,
    check_access,
    model,
    device,
    sella,
    sella_steps,
    points,
    fmax,
    relax_steps,
    charge,
    uhf,
    multiplicity,
    project,
    label,
    outdir,
):
    """Drive the primary spec bond and write a TS guess.

    Experimental: not yet better than hand-built guesses (pass-2 xTB
    result). Reads the reactant with -f and the case YAML scan block.
    The scan is 1D: the primary bond is constrained, everything else is
    relaxed, 10–15 points out to scan.to. Chi sets scan.to; there is
    no default. The maximum is taken only from converged frames. If
    that frame is the first or last point, prints ``FAIL (max at scan
    edge)``. A failed run still writes scan.xyz, the profile, and the
    result table, and does not write guess.xyz.

    The default guess is that converged maximum frame. Its Hessian mode
    is scored only when the frame converged, after translation and
    rotation are projected out, with the same pair projection as
    ``chemsmart check`` (primary rank within 3.2 Å, top pairs, Zmax as
    information). ``--sella`` is report-only: it writes guess_sella.xyz
    and does not replace guess.xyz or rerun the scan. Also writes a
    shell file of ``chemsmart queue`` TS lines. Does not submit.
    Stdout is the table; calculator noise goes to guess.log.

    ``--calc uma`` reads HF_TOKEN from the environment and does not
    print it. ``--check-access`` stops after the access check.

    Examples:

      chemsmart guess --calc uma --check-access

      chemsmart guess -f reactant.xyz --spec case.yaml --calc xtb
    """
    kind = calc.lower()
    try:
        if check_access:
            if kind == "uma":
                lines = check_uma_access(model=model)
            else:
                lines = check_xtb_import()
            click.echo("\n".join(lines))
            return
        if not filename or not spec:
            raise click.UsageError(
                "guess needs -f and --spec (or --check-access)."
            )
        case = load_case_spec(spec)
        charge_v, uhf_v, mult_v = resolve_charge_uhf(
            case,
            charge=charge,
            uhf=uhf,
            multiplicity=multiplicity,
        )
        factory = make_factory(
            kind,
            charge_v,
            uhf_v,
            mult_v,
            model=model,
            device=device,
        )
        text, ok = run_guess(
            filename,
            case,
            factory,
            charge=charge_v,
            uhf=uhf_v,
            multiplicity=mult_v,
            outdir=outdir,
            calc_name=kind,
            project=project,
            label=label,
            sella=sella,
            sella_steps=sella_steps,
            fmax=fmax,
            relax_steps=relax_steps,
            points=points,
            product_path=product,
        )
    except (CalculatorError, ValueError, NotImplementedError) as exc:
        click.echo(str(exc), err=True)
        sys.exit(1)
    click.echo(text, nl=False)
    if not ok:
        sys.exit(1)
