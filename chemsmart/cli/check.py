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
    help="Override 'primary spec bond not #1 within 3.2 Å' to CANDIDATE.",
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

    Hard reject: nimag≠1. Reject (overridable with --force): primary
    spec bond not #1 within 3.2 Å. Everything else is
    'CANDIDATE: judge mode'. Does not submit jobs.

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
