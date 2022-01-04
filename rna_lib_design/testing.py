import pandas as pd
from rna_lib_design import settings, structure
from rna_lib_design.structure_set import (
    AddType,
    StructureSet,
    HairpinStructureSet,
)

TEST_RESOURCES = settings.RESOURCES_PATH + "/testing/"


def get_test_rna_structure():
    return structure.rna_structure("GGGGAAAACCCC", "((((....))))")


def get_test_helices():
    df_path = TEST_RESOURCES + "helix_barcodes.csv"
    df = pd.read_csv(df_path)
    struct_set = StructureSet(df, AddType.HELIX)
    return struct_set


def get_test_sstrand(add_type=AddType.RIGHT):
    df_path = TEST_RESOURCES + "sstrand_barcodes.csv"
    df = pd.read_csv(df_path)
    struct_set = StructureSet(df, add_type)
    return struct_set


def get_test_hairpin(add_type=AddType.RIGHT):
    df_path = TEST_RESOURCES + "helix_barcodes.csv"
    df = pd.read_csv(df_path)
    loop = structure.rna_structure("GGAAAC", "(....)")
    struct_set = HairpinStructureSet(loop, df, add_type)
    return struct_set
