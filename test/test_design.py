import pandas as pd
import os
import pytest
from rna_lib_design import design, structure, structure_set, settings
from click.testing import CliRunner


def get_test_helices():
    df_path = settings.TEST_PATH + "/resources/helix_barcodes.csv"
    df = pd.read_csv(df_path)
    struct_set = structure_set.StructureSet(df, structure_set.AddType.HELIX)
    return struct_set


def get_test_sstrand_right():
    df_path = settings.TEST_PATH + "/resources/sstrand_barcodes.csv"
    df = pd.read_csv(df_path)
    struct_set = structure_set.StructureSet(df, structure_set.AddType.RIGHT)
    return struct_set


def get_test_sstrand_left():
    df_path = settings.TEST_PATH + "/resources/sstrand_barcodes.csv"
    df = pd.read_csv(df_path)
    struct_set = structure_set.StructureSet(df, structure_set.AddType.LEFT)
    return struct_set


def get_test_hairpin_right():
    df_path = settings.TEST_PATH + "/resources/helix_barcodes.csv"
    df = pd.read_csv(df_path)
    loop = structure.rna_structure("GGAAAC", "(....)")
    struct_set = structure_set.HairpinStructureSet(
        loop, df, structure_set.AddType.RIGHT
    )
    return struct_set


def get_test_hairpin_left():
    df_path = settings.TEST_PATH + "/resources/helix_barcodes.csv"
    df = pd.read_csv(df_path)
    loop = structure.rna_structure("GGAAAC", "(....)")
    struct_set = structure_set.HairpinStructureSet(loop, df, structure_set.AddType.LEFT)
    return struct_set


def test_design():
    rna = structure.rna_structure("GGGGAAAACCCC", "((((....))))")
    set_1 = get_test_helices()
    sol = design.get_best_design([set_1], rna, design.DesignOptions())
    assert sol.design.dot_bracket == "((((((((((....))))))))))"
    assert sol.ens_defect < 1.0


def test_design_2():
    rna = structure.rna_structure("GGGGAAAACCCC", "((((....))))")
    set_1 = get_test_helices()
    set_2 = structure_set.get_common_seq_structure_set("ref_hairpin_5prime")
    set_3 = structure_set.get_tail_structure_set()
    sol = design.get_best_design([set_1, set_2, set_3], rna, design.DesignOptions())
    assert (
        sol.design.dot_bracket
        == "....((((.....))))...((((((((((....))))))))))...................."
    )
    assert sol.ens_defect < 1.0


def test_helix_barcode():
    path = settings.TEST_PATH + "/resources/opool_final_2.csv"
    df = pd.read_csv(path)
    p5 = structure.get_common_struct("uucg_5prime")
    p3 = structure.get_common_struct("rt_tail")
    # df_new = design.helix_barcode(df, 6, p5, p3, design.DesignOptions())


def test_haiohttpelix_barcode_cli():
    runner = CliRunner()
    args = [
        "barcode",
        "--type",
        "helix",
        "--length",
        6,
        settings.TEST_PATH + "/resources/opool_final_2.csv",
    ]
    result = runner.invoke(design.cli, args, prog_name="helix_barcode")
    assert result.exit_code == 0
    assert os.path.exists("out.csv")
