import os
import numpy as np
from seq_tools import clt


def get_base_dir():
    file_path = os.path.realpath(__file__)
    spl = file_path.split("/")
    base_dir = "/".join(spl[:-1])
    return base_dir


def get_default_args():
    p = {
        "ss"       : None,
        "trim5"    : None,
        "trim3"    : None,
        "add5"     : None,
        "add3"     : None,
        "output"   : "output.csv",
        "to_rna"   : False,
        "to_dna"   : False,
        "remove_t7": False,
        "add_t7"   : False,
        "ds"       : False,
        "fold"     : False,
        "calc"     : None
    }
    return p


# tests ##############################################################################

def test_seq_len():
    p = get_default_args()
    p["input"] = "GGGAAACCC"
    p["calc"] = "len"
    df = clt.run_seq_tools(p)
    row = df.loc[1]
    assert "len" in df.columns
    assert row["type"] == "DNA"
    assert row["len"] == 9


def test_csv_mw():
    p = get_default_args()
    base_dir = get_base_dir()
    csv_path = base_dir + "/resources/opool_final.csv"
    p["input"] = csv_path
    p["calc"] = "mw"
    df = clt.run_seq_tools(p)
    assert "molecular weight" in df.columns
    row = df.loc[1]
    assert row["ds"] == False
    assert row["molecular weight"] == 24537


def test_mult_csv():
    p = get_default_args()
    base_dir = get_base_dir()
    csv_path = base_dir + "/resources/opool_final.csv"
    p["input"] = csv_path
    p["calc"] = "mw,len"
    df = clt.run_seq_tools(p)
    assert "molecular weight" in df.columns
    assert "len" in df.columns


def test_csv_convert_to_dna():
    p = get_default_args()
    base_dir = get_base_dir()
    csv_path = base_dir + "/resources/opool_final.csv"
    p["input"] = csv_path
    p["add_t7"] = True
    p["to_dna"] = True
    df = clt.run_seq_tools(p)
    row = df.loc[1]
    assert (
            row["sequence"]
            == "TTCTAATACGACTCACTATAGATATGGATGATTAGGACATGCATTGCTGAGGGGAAACTTTTTGCAATGCAACAG"
               "CCAAATCGTCCTAAGTC"
    )


def test_trim():
    p = get_default_args()
    base_dir = get_base_dir()
    csv_path = base_dir + "/resources/opool_final.csv"
    p["input"] = csv_path
    p["calc"] = "l"
    df = clt.run_seq_tools(p)
    total = np.sum(df["len"])
    assert total == 430
    p["trim5"] = 5
    df = clt.run_seq_tools(p)
    total = np.sum(df["len"])
    assert total == 400


def test_ec():
    p = get_default_args()
    base_dir = get_base_dir()
    p["input"] = (
        "GGAAGAUCGAGUAGAUCAAAGAGCCUAUGGCUGCCACCCGAGCCCUUGAACUACAGGGAACACUGGAAACAGUACCCCCU"
        "GCAAGGGCGUUUGACGGUGGCAGCCUAAGGGCUCAAAGAAACAACAACAACAAC"
    )
    p["calc"] = "ec"
    df = clt.run_seq_tools(p)
    assert "extinction coeff" in df
    row = df.loc[1]
    assert row["extinction coeff"] == 1340200

    p["ss"] = (
        "....((((.....))))...((((((..((((((((((((((((((((.....(((((...((((....))))...))"
        "))))))))))..)))..))))))))))...))))))...................."
    )

    df = clt.run_seq_tools(p)
    row = df.loc[1]
    assert row["extinction coeff"] == 1252307
