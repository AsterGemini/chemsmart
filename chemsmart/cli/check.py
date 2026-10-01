"""CLI: check a Gaussian log against a case spec YAML."""

import sys

import click

from chemsmart.analysis.tscheck import check_ts_log, load_case_spec


@click.command("check")
@click.argument("log", type=click.Path(exists=True, dir_okay=False))
@click.option(
    "--spec",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Case YAML with reaction-coordinate bonds and endpoint criteria.",
)
@click.option(
    "--force",
    is_flag=True,
    default=False,
    help=(
        "Override a displaced-bond or zero-projection reject "
        "to CANDIDATE. Does not override nimag≠1."
    ),
)
@click.option(
    "--kind",
    type=click.Choice(["ts", "irc", "opt"], case_sensitive=False),
    default=None,
    help="Log kind. Default: detect from the route.",
)
@click.option(
    "-p",
    "--project",
    type=str,
    default=None,
    help="Project name used in printed QRC suggestion commands.",
)
@click.option(
    "--ts-file",
    type=click.Path(exists=True, dir_okay=False),
    default=None,
    help="TS log for printed QRC commands (IRC checks).",
)
def check(log, spec, force, kind, project, ts_file):
    """Check a Gaussian TS, IRC, or endpoint-opt log.

    Prints nimag, imaginary frequency, charge/mult, the top-8 pair
    projections, spec-bond ranks, and Zmax (information only).

    Hard TS reject that --force does not override: nimag≠1. Rejects
    that --force overrides: primary spec bond not #1 within 3.2 Å, and
    spec bonds with ~0 projection ('spec bonds not displaced in mode',
    no rank). A small extra imaginary mode still counts in nimag and
    prints 'small imag present (<15 cm⁻¹); counted in nimag'.
    Endpoint nimag≠0 is 'ENDPOINT NOT MIN'. Failed endpoint bonds are
    'connectivity mismatch: judge'. Unset YAML criteria print
    'endpoint criteria unset'. Does not submit jobs.

    Examples:

      chemsmart check ts.log --spec case1.yaml

      chemsmart check ircf.log --spec case3.yaml --ts-file ts.log
    """
    case = load_case_spec(spec)
    report = check_ts_log(
        log,
        case,
        force=force,
        project=project,
        ts_file=ts_file,
        kind=kind.lower() if kind else None,
    )
    click.echo(report.text(), nl=False)
    if report.verdict == "REJECT" and not report.force_used:
        sys.exit(1)
