import pytest
import shutil
from pathlib import Path
import pandas as pd

from click.testing import CliRunner
from seq_tools import has_5p_sequence, has_3p_sequence
from rna_lib_design import cli
from rna_lib_design.settings import get_resources_path, get_test_path


TEST_RESOURCES = get_test_path() / "resources"


def test_cli():
    runner = CliRunner()
    result = runner.invoke(cli.cli, ["--help"])
    assert result.exit_code == 0
    assert "Show this message and exit." in result.output


def test_edit_distance():
    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "edit-distance",
            str(TEST_RESOURCES / "libs/minittr2.csv"),
        ],
    )
    assert result.exit_code == 0


def test_add_common():
    runner = CliRunner()
    result = runner.invoke(
        cli.cli,
        [
            "add-common",
            str(TEST_RESOURCES / "libs/minittr2.csv"),
        ],
    )
    assert result.exit_code == 0
    assert Path("results").is_dir()
    assert Path("results/results-rna.csv").is_file()
    assert Path("results/results-dna.csv").is_file()
    df_rna = pd.read_csv("results/results-rna.csv")
    assert has_5p_sequence(df_rna, "GGAACAGCACUUCGGUGCAAA")
    assert has_3p_sequence(df_rna, "AAAGAAACAACAACAACAAC")
    df_dna = pd.read_csv("results/results-dna.csv")
    assert has_5p_sequence(df_dna, "TTCTAATACGACTCACTATA")
    shutil.rmtree("results")


class TestBarcode:
    def test_standard(self):
        runner = CliRunner()
        result = runner.invoke(
            cli.cli,
            [
                "barcode",
                str(TEST_RESOURCES / "libs/minittr2.csv"),
            ],
        )
        assert result.exit_code == 0
        assert Path("results").is_dir()
        shutil.rmtree("results")
    
    def test_helix(self):
        runner = CliRunner()
        result = runner.invoke(
            cli.cli,
            [
                "barcode",
                "--btype", "helix",
                str(TEST_RESOURCES / "libs/minittr2.csv"),
            ],
        )
        assert result.exit_code == 0
        assert Path("results").is_dir()
        shutil.rmtree("results")
    
    def test_5p_hairpin(self):
        runner = CliRunner()
        result = runner.invoke(
            cli.cli,
            [
                "barcode",
                "--btype", "5p_hairpin",
                str(TEST_RESOURCES / "libs/minittr2.csv"),
            ],
        )
        assert result.exit_code == 0
        assert Path("results").is_dir()
        shutil.rmtree("results")
 