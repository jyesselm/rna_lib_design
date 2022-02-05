import os
import shutil
from pathlib import Path
from click.testing import CliRunner

from rna_lib_design.barcode import *
from rna_lib_design import testing
from rna_lib_design.cli import CLIParser

TEST_RESOURCES = testing.TEST_RESOURCES


def test_single_barcode():
    path = TEST_RESOURCES + "opool_final_2.csv"
    df = pd.read_csv(path)
    bcoder = SingleBarcoder("helix", 6)
    df = bcoder.barcode(df, DesignOptions())
    assert "org_sequence" in df
    assert len(df["sequence"][0]) == 87
    # now add p5 and p3
    args = {
        "no_p5": False,
        "p5_sequence": None,
        "p5_structure": None,
        "p5_name": None,
        "no_p3": False,
        "p3_sequence": None,
        "p3_structure": None,
        "p3_name": None,
    }
    p5 = CLIParser.get_p5(args)
    p3 = CLIParser.get_p3(args)
    df = pd.read_csv(path)
    bcoder = SingleBarcoder("helix", 6)
    bcoder.set_p5_and_p3(p5, p3)
    df = bcoder.barcode(df, DesignOptions())
    assert len(df["sequence"][0]) == 127


def test_double_barcode():
    path = TEST_RESOURCES + "opool_final_2.csv"
    df = pd.read_csv(path)
    args = {"loop_sequence": None, "loop_structure": None, "loop_name": None}
    loop = CLIParser.get_loop(args)
    bcoder = DoubleBarcode("hairpin_helix", (6, 6), loop)
    df = bcoder.barcode(df, DesignOptions())
    assert "org_sequence" in df
