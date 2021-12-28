import click
import pandas as pd
from tabulate import tabulate
from rna_lib_design import logger, settings

log = logger.setup_applevel_logger()

@click.group()
def cli():
    pass


@cli.command()
@click.argument('rtype')
@click.option('-lmin', type=int)
@click.option('-lmax', type=int)
@click.option('-smin', type=int)
@click.option('-smax', type=int)
def maxdiff(rtype, lmin, lmax, smin, smax):
    rtype = rtype.lower()
    if rtype == 'helix':
        path = settings.RESOURCES_PATH + "barcodes/helices.csv"
    else:
        path = ""
        log.error(f"{rtype} type is not supported")
        exit()
    df = pd.read_csv(path)
    if lmin is None:
        lmin = df['length'].min()
    if lmax is None:
        lmax = df['length'].max()
    if smin is None:
        smin = df['size'].min()
    if smax is None:
        smax = df['size'].max()
    df = df[(df['length'] <= lmax) & (df['length'] >= lmin)]
    df = df[(df['size'] <= smax) & (df['size'] >= smin)]
    df = df.sort_values(by='diff', ascending=False)
    df_sub = df.iloc[0:5]
    print(tabulate(df_sub, headers='keys', tablefmt='psql'))


if __name__ == "__main__":
    cli()
