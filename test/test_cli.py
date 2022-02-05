import pytest
import shutil
import pandas as pd
from click.testing import CliRunner
from rna_lib_design.cli import CLIParser
from rna_lib_design import cli
from rna_lib_design import testing

TEST_RESOURCES = testing.TEST_RESOURCES


def test_setup_dataframe():
    path = TEST_RESOURCES + "opool_final_2.csv"
    args = {"csv": path, "trim_5p": 0, "trim_3p": 0}
    df = CLIParser.get_dataframe(args)
    assert len(df) == 2
    assert (
        df.loc[0]["sequence"]
        == "GAUAUGGAUAGAGUAAGAGAGAUGGAAGUCUCAGGGGAAACUUUGAGAUGGACGGUUUACAAGUUGUCCUAAGUC"
    )


def test_setup_dataframe_2():
    path = TEST_RESOURCES + "opool_final.csv"
    args = {"csv": path, "trim_5p": 0, "trim_3p": 0}
    df = CLIParser.get_dataframe(args)
    assert "structure" in df
    assert "mfe" in df


def test_setup_dataframe_3():
    path = TEST_RESOURCES + "opool_final_2.csv"
    args = {"csv": path, "trim_5p": 2, "trim_3p": 2}
    df = CLIParser.get_dataframe(args)
    assert (
        df.loc[0]["sequence"]
        == "UAUGGAUAGAGUAAGAGAGAUGGAAGUCUCAGGGGAAACUUUGAGAUGGACGGUUUACAAGUUGUCCUAAG"
    )


def test_get_p5_default():
    args = {
        "no_p5": False,
        "p5_sequence": None,
        "p5_structure": None,
        "p5_name": None,
    }
    p5 = CLIParser.get_p5(args)
    struct = p5.get_random()[0]
    assert struct.sequence == "GGAAGAUCGAGUAGAUCAAA"
    assert struct.dot_bracket == "....((((.....))))..."


def test_get_p5_by_name():
    args = {
        "no_p5": False,
        "p5_sequence": None,
        "p5_structure": None,
        "p5_name": "uucg_p5_rev_primer",
    }
    p5 = CLIParser.get_p5(args)
    struct = p5.get_random()[0]
    assert struct.sequence == "GGAACAGCACUUCGGUGCAAA"
    assert struct.dot_bracket == "......((((....))))..."


def test_get_p5_by_seq():
    args = {
        "no_p5": False,
        "p5_sequence": "GGAACAGCACUUCGGUGCAAA",
        "p5_structure": None,
        "p5_name": None,
    }
    p5 = CLIParser.get_p5(args)
    struct = p5.get_random()[0]
    assert struct.sequence == "GGAACAGCACUUCGGUGCAAA"
    assert struct.dot_bracket == "......((((....))))..."


def test_get_p5_raises_supplied_name_and_seq():
    args = {
        "no_p5": False,
        "p5_sequence": "GGAACAGCACUUCGGUGCAAA",
        "p5_structure": None,
        "p5_name": "uucg_p5_rev_primer",
    }
    with pytest.raises(ValueError):
        CLIParser.get_p5(args)


def test_get_p3_default():
    args = {
        "no_p3": False,
        "p3_sequence": None,
        "p3_structure": None,
        "p3_name": None,
    }
    p3 = CLIParser.get_p3(args)
    struct = p3.get_random()[0]
    assert struct.sequence == "AAAGAAACAACAACAACAAC"


def test_get_p5_buffer_default():
    args = {"p5b_sequence": None, "p5b_structure": None}
    p5b = CLIParser.get_p5_buffer(args)
    assert p5b is None


def test_get_p5_buffer_sequence():
    args = {"p5b_sequence": "AC", "p5b_structure": None}
    p5b = CLIParser.get_p5_buffer(args)
    struct = p5b.get_random()[0]
    assert struct.sequence == "AC"
    assert struct.dot_bracket == ".."


def test_get_p5_buffer_structure():
    args = {"p5b_sequence": "AC", "p5b_structure": "()"}
    p5b = CLIParser.get_p5_buffer(args)
    struct = p5b.get_random()[0]
    assert struct.sequence == "AC"
    assert struct.dot_bracket == "()"


def test_get_loop_default():
    args = {"loop_sequence": None, "loop_structure": None, "loop_name": None}
    loop = CLIParser.get_loop(args)
    assert loop.sequence == "CUUCGG"
    assert loop.dot_bracket == "(....)"


def test_get_loop_sequence():
    args = {
        "loop_sequence": "CGCGAGUAGCG",
        "loop_structure": None,
        "loop_name": None,
    }
    loop = CLIParser.get_loop(args)
    assert loop.sequence == "CGCGAGUAGCG"
    assert loop.dot_bracket == "(((.....)))"


def test_get_loop_structure():
    args = {
        "loop_sequence": "CGAGUAG",
        "loop_structure": "(.....)",
        "loop_name": None,
    }
    loop = CLIParser.get_loop(args)
    assert loop.sequence == "CGAGUAG"
    assert loop.dot_bracket == "(.....)"


def test_barcode_default():
    runner = CliRunner()
    args = [
        "barcode",
        TEST_RESOURCES + "opool_final_2.csv",
    ]
    result = runner.invoke(cli.cli, args, prog_name="helix_barcode")
    assert result.exit_code == 0
    df = pd.read_csv("results/results-opool.csv")
    assert len(df) == 2
    shutil.rmtree("results")


def test_barcode_hairpin():
    runner = CliRunner()
    args = [
        "barcode",
        "--btype",
        "hairpin",
        TEST_RESOURCES + "opool_final_2.csv",
    ]
    result = runner.invoke(cli.cli, args, prog_name="hairpin_barcode")
    assert result.exit_code == 0
    df = pd.read_csv("results/results-rna.csv")
    assert df.loc[0]["sequence"].find('CUUCGG') != -1
    shutil.rmtree("results")

