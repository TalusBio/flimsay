import rich_click as click

from .blib import fill_blib


@click.group()
def cli():
    pass


@cli.command(name="fill_blib")
@click.argument("blib")
@click.argument("out_blib")
@click.option(
    "--overwrite",
    is_flag=True,
    help="Whether to overwrite output file, if it exists",
)
def fill_blib_cli(blib, out_blib, overwrite):
    """Add ion mobility prediction to a .blib file."""
    fill_blib(blib, out_blib, overwrite)


if __name__ == "__main__":
    cli()
