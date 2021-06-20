import click
import sys
import pandas as pd
import numpy as np
import xlsxwriter

from seq_tools import structure, vienna, extinction_coeff

mw_rna = {"A": 347.2, "C": 323.2, "G": 363.2, "U": 324.2}

mw_dna = {"A": 331.2, "C": 307.2, "G": 347.2, "T": 322.2}


def get_type(seq):
    t = "DNA"
    for e in seq:
        if e == "U":
            t = "RNA"
    return t


def get_mw(seq, t, ds):
    total = 0
    for e in seq:
        if t == "RNA":
            total += mw_rna[e]
        else:
            total += mw_dna[e]
    if ds:
        rc = get_reverse_complement(seq, t)
        for e in rc:
            if t == "RNA":
                total += mw_rna[e]
            else:
                total += mw_dna[e]
    return total



def convert_to_rna(seq):
    new_seq = ""
    for e in seq:
        if e == "T":
            new_seq += "U"
        else:
            new_seq += e
    return new_seq


def get_reverse_complement(seq, type):
    complement = ""
    for e in seq:
        if e == "A" and type == "RNA":
            complement += "U"
        elif e == "C":
            complement += "G"
        elif e == "G":
            complement += "C"
        elif e == "U":
            complement += "A"
        elif e == "T":
            complement += "A"
    return complement[::-1]


def get_df_from_seq_and_ss(seq, ss):
    df = pd.DataFrame(columns="name,sequence,structure".split(","))
    df.loc[1] = ["seq", seq, ss]
    return df


# apply to dataframe
def update_sequence_to_dna(df):
    seqs = []
    for i, row in df.iterrows():
        seqs.append(convert_to_dna(row["sequence"]))
    df["sequence"] = seqs
    df["type"] = "DNA"


def update_sequence_to_rna(df):
    seqs = []
    for i, row in df.iterrows():
        seqs.append(convert_to_rna(row["sequence"]))
    df["sequence"] = seqs
    df["type"] = "RNA"


def fold_rna(df):
    ss = []
    for i, row in df.iterrows():
        if row["type"] != "RNA":
            print("only RNA can be folded, sequence is DNA!")
            exit()
        ss.append(str(vienna.fold(row["sequence"]).dot_bracket))
    df["structure"] = ss


def add_mw_to_df(df):
    mws = []
    for i, row in df.iterrows():
        mws.append(round(get_mw(row["sequence"], row["type"], row["ds"])))
    df["molecular weight"] = mws


def add_ec_to_df(df):
    ecs = []
    for i, row in df.iterrows():
        if row["type"] == "RNA":
            ec = extinction_coeff.get_coefficient_rna(row["sequence"], row["structure"])
        else:
            ec = extinction_coeff.get_coefficient_dna(row["sequence"], row["ds"])
        ecs.append(ec)
    df["extinction coeff"] = ecs


def add_rc_to_df(df):
    rcs = []
    for i, row in df.iterrows():
        rcs.append(get_reverse_complement(row["sequence"], row["type"]))
    df["reverse complement"] = rcs


def trim_5p(df, length):
    seqs = []
    ss = []
    include_ss = False
    if df.loc[1]["structure"] is not None:
        include_ss = True
    for i, row in df.iterrows():
        seqs.append(row["sequence"][length:])
        if include_ss:
            ss.append(row["structure"][length:])
    df["sequence"] = seqs
    if include_ss:
        df["structure"] = ss


def trim_3p(df, length):
    seqs = []
    ss = []
    include_ss = False
    if df.loc[1]["structure"] is not None:
        include_ss = True
    for i, row in df.iterrows():
        seqs.append(row["sequence"][:-length])
        if include_ss:
            ss.append(row["structure"][:-length])
    df["sequence"] = seqs
    if include_ss:
        df["structure"] = ss


def add_5p(df, p5):
    spl = p5.split(",")
    seq_p5 = p5
    ss_p5 = ""
    include_ss = False
    if len(spl) == 2:
        include_ss = True
        seq_p5 = p5[0]
        ss_p5 = p5[1]
    seqs = []
    ss = []
    if include_ss:
        if "structure" not in df.columns:
            print("cannot add to sequence and structure, structure is not defined!")
            exit()
        if df.loc[1]["structure"] is None:
            print("cannot add to sequence and structure, structure is not defined!")
            exit()
    for i, row in df.iterrows():
        seqs.append(seq_p5 + row["sequence"])
        if include_ss:
            ss.append(ss_p5 + row["structure"])
    df["sequence"] = seqs
    if include_ss:
        df["structure"] = ss


def add_3p(df, p3):
    spl = p3.split(",")
    seq_p3 = p3
    ss_p3 = ""
    include_ss = False
    if len(spl) == 2:
        include_ss = True
        seq_p3 = p3[0]
        ss_p3 = p3[1]
    seqs = []
    ss = []
    if include_ss:
        if "structure" not in df.columns:
            print("cannot add to sequence and structure, structure is not defined!")
            exit()
        if df.loc[1]["structure"] is None:
            print("cannot add to sequence and structure, structure is not defined!")
            exit()
    for i, row in df.iterrows():
        seqs.append(row["sequence"] + seq_p3)
        if include_ss:
            ss.append(row["structure"] + ss_p3)
    df["sequence"] = seqs
    if include_ss:
        df["structure"] = ss


def remove_t7_from_df(df):
    t7 = "TTCTAATACGACTCACTATA"
    type = df.loc[1]["type"]
    if type != "DNA":
        print("cannot remove t7 promoter to RNA!")
        exit()
    seqs = []
    for i, row in df.iterrows():
        if row["sequence"].find(t7) == -1:
            print("sequence does not contain the t7 promoter cannot remove it!")
            seqs.append(row["sequence"])
            continue
        seqs.append(row["sequence"][20:])
    df["sequence"] = seqs


def is_seq_or_csv(input):
    input = input.upper()
    for e in input:
        if e == 'A' or e == 'C' or e == 'G' or e == 'U' or e =='T' or e == 'N':
            continue
        else:
            return "CSV"
    return "SEQ"

@click.command()
@click.argument("input")
@click.option("-ss", required=False, default=None)
@click.option("-c", "--calc", required=False, default=None)
@click.option("-trim5", required=False, default=None)
@click.option("-trim3", required=False, default=None)
@click.option("-add5", required=False, default=None)
@click.option("-add3", required=False, default=None)
@click.option("-o", "--output", required=False, default="output.csv")
@click.option("-to_rna", required=False, is_flag=True, default=False)
@click.option("-to_dna", required=False, is_flag=True, default=False)
@click.option("-to_fasta", required=False, is_flag=True, default=False)
@click.option("-remove_t7", required=False, is_flag=True, default=False)
@click.option("-add_t7", required=False, is_flag=True, default=False)
@click.option("-ds", required=False, is_flag=True, default=False)
@click.option("-fold", required=False, is_flag=True, default=False)
@click.option("-avg", required=False, is_flag=True, default=False)
def main(**p):
    if p["to_rna"] and p["to_dna"]:
        print("cannot use flag -to_rna and -to_dna")
        exit()
    if p["remove_t7"] and p["add_t7"]:
        print("cannot use flag -remove_t7 and -add_t7")
        exit()
    if p["add5"] and p["trim5"]:
        print("cannot use flag -trim5 and -add5")
        exit()
    if p["add3"] and p["trim3"]:
        print("cannot use flag -trim3 and -add3")
        exit()

    df = run_seq_tools(p)
    df = df.drop(["type", "ds"], axis=1)
    if p["to_fasta"]:
        f = open("test.fasta", "w")
        for i, row in df.iterrows():
            f.write(f">{row['name']}\n")
            f.write(row["sequence"]+"\n")
        f.close()
    if len(df) == 1:
        df = df.drop(["name"], axis=1)
        d = df.loc[1].to_dict()
        for k, v in d.items():
            if v == False or v == None:
                continue
            print(f"{k} : {v}")
    else:
        if "structure" is df.columns:
            if df.loc[1]["structure"] is None:
                df = df.drop("structure", axis=1)
        print("outputing results to csv: output.csv")
        df.to_csv("output.csv", index=False)
        if p["avg"]:
            excess = "name,sequence,structure,rc".split(",")
            for c in excess:
                if c in df.columns:
                    df = df.drop(c, axis=1)
            cols = df.columns
            print("averages: ")
            for c in cols:
                try:
                    avg = np.mean(df[c])
                    print(c + " : " + str(round(avg, 2)))
                except:
                    pass


def run_seq_tools(p):
    input_type = is_seq_or_csv(p["input"])
    if input_type == "CSV":
        df = pd.read_csv(p["input"])
    else:
        df = get_df_from_seq_and_ss(p["input"], p["ss"])
    calcs = {}
    if p["calc"] is not None:
        spl = p["calc"].split(",")
        for e in spl:
            calcs[e] = True;
    if "name" not in df:
        names = []
        for i, row in df.iterrows():
            names.append(f"seq_{i}")
        df["name"] = names
    type = get_type(df.loc[1]["sequence"])
    df["type"] = type
    if "ds" not in df.columns:
        df["ds"] = p["ds"]
    if p["remove_t7"]:
        remove_t7_from_df(df)
    if p["to_rna"]:
        update_sequence_to_rna(df)
        type = "RNA"
    if p["to_dna"]:
        update_sequence_to_dna(df)
        type = "DNA"
    if p["add_t7"]:
        add_5p(df, "TTCTAATACGACTCACTATA")
    if "structure" not in df.columns and type == "RNA":
        df["structure"] = None
    if p["trim5"] is not None:
        trim_5p(df, int(p["trim5"]))
    if p["trim3"] is not None:
        trim_3p(df, int(p["trim3"]))
    if p["fold"]:
        fold_rna(df)

    # things to add to df
    if "rc" in calcs:
        add_rc_to_df(df)
    if "len" in calcs or "l" in calcs:
        lens = []
        for i, row in df.iterrows():
            lens.append(len(row["sequence"]))
        df["len"] = lens
    if "ec" in calcs:
        add_ec_to_df(df)
    if "mw" in calcs:
        add_mw_to_df(df)
    return df


@click.command()
@click.argument("csv")
@click.option("-n", "--name", required=False, default=None)
def gen_opool(csv, name):
    df = pd.read_csv(csv)
    if name is None:
        if "pool_name" not in df:
            print("must supply pool name wiht -n/--name or have pool_name in csv")
            exit()
    else:
        df["pool_name"] = name
    workbook = xlsxwriter.Workbook("test.xlsx")
    worksheet = workbook.add_worksheet()
    worksheet.write("A1", "Pool name")
    worksheet.write("B1", "Sequence")
    pos = 2
    for i, row in df.iterrows():
        worksheet.write(f"A{pos}", row["pool_name"])
        worksheet.write(f"B{pos}", row["sequence"])
        pos += 1
    workbook.close()



if __name__ == "__main__":
    main()
