import os
import cloup
import numpy as np
from cloup import option_group, option
import pandas as pd
from pathlib import Path
from typing import Dict, Tuple, Optional
from dataclasses import dataclass

import vienna
from seq_tools import fold, trim

from rna_lib_design import logger, settings, util
from rna_lib_design.design import (
    Designer,
    DesignOpts,
)

# TODO add cli tool "list" to list all the available p5, p3, and common sequences
# TODO add cli tool to list which barcode sets are available

def validate_barcode_type(btype : str):
    pass


@cloup.group()
def cli():
    pass

@cloup.argument("csv", type=cloup.Path(exists=True))
@cloup.option("-t", "--btype", type=str)
@cloup.option("--param-file", type=cloup.Path(exists=True))
def barcode(csv, btype, param_file):
    pass

@cli.command()
@cloup.argument("csv", type=cloup.Path(exists=True))
def editdistance(csv):
    df = pd.read_csv(csv)
    print(util.compute_edit_distance(df))


if __name__ == "__main__":
    cli()
