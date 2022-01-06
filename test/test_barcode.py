import os
from pathlib import Path
from click.testing import CliRunner

from rna_lib_design.barcode import *
from rna_lib_design import testing

TEST_RESOURCES = testing.TEST_RESOURCES


def remove_output(name):
    os.remove(f"{name}-all.csv")
    os.remove(f"{name}-rna.csv")
    os.remove(f"{name}-dna.csv")
    os.remove(f"{name}-opool.xlsx")


def test_setup_dataframe():
    path = TEST_RESOURCES + "opool_final_2.csv"
    args = {"csv": path, "trim_5p": 0, "trim_3p": 0}
    df = setup_dataframe_from_cli(args)
    assert len(df) == 2
    assert (
        df.loc[0]["sequence"]
        == "GAUAUGGAUAGAGUAAGAGAGAUGGAAGUCUCAGGGGAAACUUUGAGAUGGACGGUUUACAAGUUGUCCUAAGUC"
    )


def test_setup_dataframe_2():
    path = TEST_RESOURCES + "opool_final.csv"
    args = {"csv": path, "trim_5p": 0, "trim_3p": 0}
    df = setup_dataframe_from_cli(args)
    assert "structure" in df
    assert "mfe" in df


def test_setup_dataframe_3():
    path = TEST_RESOURCES + "opool_final_2.csv"
    args = {"csv": path, "trim_5p": 2, "trim_3p": 2}
    df = setup_dataframe_from_cli(args)
    assert (
        df.loc[0]["sequence"]
        == "UAUGGAUAGAGUAAGAGAGAUGGAAGUCUCAGGGGAAACUUUGAGAUGGACGGUUUACAAGUUGUCCUAAG"
    )


def test_single_barcode():
    path = TEST_RESOURCES + "opool_final_2.csv"
    df = pd.read_csv(path)
    p5 = defaults.get_p5_from_str(None)
    p3 = defaults.get_p3_from_str(None)
    add_type = structure_set.AddType.RIGHT
    df_results = single_barcode(df, "helix", 6, p5, p3, DesignOptions())
    assert len(df) == 2
    assert "org_sequence" in df


def test_helix_barcode_cli():
    runner = CliRunner()
    args = [
        "barcode",
        TEST_RESOURCES + "opool_final_2.csv",
    ]
    result = runner.invoke(cli, args, prog_name="helix_barcode")
    assert result.exit_code == 0
    fname = Path("opool_final_2.csv").stem
    assert os.path.exists(f"{fname}-all.csv")
    assert os.path.exists(f"{fname}-opool.xlsx")
    df = pd.read_csv(f"{fname}-all.csv")
    assert "org_sequence" in df
    remove_output(fname)


def test_hairpin_barcode_cli():
    runner = CliRunner()
    args = [
        "barcode",
        "--btype",
        "hairpin",
        TEST_RESOURCES + "opool_final_2.csv",
    ]
    result = runner.invoke(cli, args, prog_name="helix_barcode")
    assert result.exit_code == 0
    fname = Path("opool_final_2.csv").stem
    assert os.path.exists(f"{fname}-all.csv")
    assert os.path.exists(f"{fname}-opool.xlsx")
    df = pd.read_csv(f"{fname}-all.csv")
    ss = df.loc[0]["structure"]
    assert ss.find("(((((((....)))))))") != -1
    remove_output(fname)
