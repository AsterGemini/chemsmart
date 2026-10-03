"""CLI: write, never submit, the next-stage chemsmart command file."""

import os

import click

from chemsmart.utils.cli import MyCommand, MyGroup

HEADER = """#!/bin/sh
# Written by `chemsmart queue`. Not submitted.
# Review, then run this file yourself (every submit stays manual).
set -e
"""


def _write(path, commands):
    text = HEADER + "\n".join(commands) + "\n"
    folder = os.path.dirname(os.path.abspath(path))
    if folder:
        os.makedirs(folder, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(text)
    try:
        os.chmod(path, 0o755)
    except OSError:
        pass
    return path


def _ts_cmd(project, filename, charge, multiplicity, label):
    bits = [
        "chemsmart sub gaussian",
        f"-p {project}",
        f"-f {filename}",
        f"-c {charge}",
        f"-m {multiplicity}",
    ]
    if label:
        bits.extend(["-l", label])
    bits.extend(["--additional-opt-options", "maxstep=5", "ts"])
    return " ".join(bits)


def _irc_cmd(project, filename, charge, multiplicity, label):
    bits = [
        "chemsmart sub gaussian",
        f"-p {project}",
        f"-f {filename}",
        f"-c {charge}",
        f"-m {multiplicity}",
    ]
    if label:
        bits.extend(["-l", label])
    bits.append("irc")
    return " ".join(bits)


def _opt_cmd(project, filename, charge, multiplicity, label):
    bits = [
        "chemsmart sub gaussian",
        f"-p {project}",
        f"-f {filename}",
        f"-c {charge}",
        f"-m {multiplicity}",
    ]
    if label:
        bits.extend(["-l", label])
    bits.append("opt")
    return " ".join(bits)


@click.group("queue", cls=MyGroup)
@click.option(
    "-p",
    "--project",
    type=str,
    default="zn5",
    show_default=True,
    help="Gaussian project name (cluster yaml, not edited here).",
)
@click.option(
    "-c",
    "--charge",
    type=int,
    default=0,
    show_default=True,
)
@click.option(
    "-m",
    "--multiplicity",
    type=int,
    default=1,
    show_default=True,
)
@click.option(
    "-l",
    "--label",
    type=str,
    default=None,
    help="Job label forwarded to chemsmart sub.",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(dir_okay=False),
    default="queue.sh",
    show_default=True,
    help="Shell file to write. Never submitted by this command.",
)
@click.pass_context
def queue(ctx, project, charge, multiplicity, label, output):
    """Write the next-stage chemsmart commands to a shell file.

    Stages: ts -> irc (fwd+rev in one chemsmart irc job) -> opt of
    the IRC endpoints (two commands in one file). QRC is not queued;
    `chemsmart check` may print a qrc ± suggestion.

    Does not submit Gaussian or HPC jobs.
    """
    ctx.ensure_object(dict)
    ctx.obj["project"] = project
    ctx.obj["charge"] = charge
    ctx.obj["multiplicity"] = multiplicity
    ctx.obj["label"] = label
    ctx.obj["output"] = output


@queue.command("ts", cls=MyCommand)
@click.option(
    "-f",
    "--filename",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="TS guess geometry (xyz/log/com).",
)
@click.pass_context
def queue_ts(ctx, filename):
    """Write a Gaussian TS/freq submit command (maxstep=5)."""
    path = _write(
        ctx.obj["output"],
        [
            _ts_cmd(
                ctx.obj["project"],
                filename,
                ctx.obj["charge"],
                ctx.obj["multiplicity"],
                ctx.obj["label"],
            )
        ],
    )
    click.echo(f"Wrote {path} (not submitted)")


@queue.command("irc", cls=MyCommand)
@click.option(
    "-f",
    "--filename",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="TS/freq log used as the IRC start.",
)
@click.pass_context
def queue_irc(ctx, filename):
    """Write one IRC command; chemsmart runs forward and reverse."""
    path = _write(
        ctx.obj["output"],
        [
            _irc_cmd(
                ctx.obj["project"],
                filename,
                ctx.obj["charge"],
                ctx.obj["multiplicity"],
                ctx.obj["label"],
            )
        ],
    )
    click.echo(f"Wrote {path} (not submitted)")
    click.echo(
        "IRC fwd+rev are one chemsmart irc job (batched). Run the file."
    )


@queue.command("opt", cls=MyCommand)
@click.option(
    "--ircf",
    type=click.Path(exists=True, dir_okay=False),
    default=None,
    help="Forward IRC log.",
)
@click.option(
    "--ircr",
    type=click.Path(exists=True, dir_okay=False),
    default=None,
    help="Reverse IRC log.",
)
@click.pass_context
def queue_opt(ctx, ircf, ircr):
    """Write batched endpoint-opt commands for IRC ends."""
    if not ircf and not ircr:
        raise click.UsageError("Provide --ircf and/or --ircr.")
    cmds = []
    label = ctx.obj["label"] or "irc"
    if ircf:
        cmds.append(
            _opt_cmd(
                ctx.obj["project"],
                ircf,
                ctx.obj["charge"],
                ctx.obj["multiplicity"],
                f"{label}_ircf_opt",
            )
        )
    if ircr:
        cmds.append(
            _opt_cmd(
                ctx.obj["project"],
                ircr,
                ctx.obj["charge"],
                ctx.obj["multiplicity"],
                f"{label}_ircr_opt",
            )
        )
    path = _write(ctx.obj["output"], cmds)
    click.echo(f"Wrote {path} (not submitted)")
    click.echo("Endpoint opts are batched in this file. Run it yourself.")
