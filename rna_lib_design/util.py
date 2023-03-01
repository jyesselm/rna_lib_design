import pandas as pd
import numpy as np
import random
import editdistance
from dataclasses import dataclass
from typing import List

from seq_tools.sequence import to_dna, get_reverse_complement
from seq_tools.dataframe import has_5p_sequence, has_3p_sequence
from rna_lib_design.logger import get_logger

BASEPAIRS = ["AU", "UA", "GC", "CG", "GU", "UG"]
BASEPAIRS_WC = ["AU", "UA", "GC", "CG"]
BASEPAIRS_GU = ["GU", "UG"]

log = get_logger("UTIL")


def get_seq_fwd_primer_code(df: pd.DataFrame) -> str:
    """
    gets the sequence forward primer code
    :param df: the dataframe with sequences
    """
    df = df.copy()
    df = to_dna(df)
    path = os.path.join(LIB_PATH, "resources", "p5_sequences.csv")
    df_p5 = pd.read_csv(path)
    for _, row in df_p5.iterrows():
        # if all sequences in df start with the p5 sequence then return the p5 code
        if all(df["sequence"].str.startswith(row["sequence"])):  # type: ignore
            return row["code"]
    return ""


def get_primer_dataframe(file_path: str) -> pd.DataFrame:
    df = pd.read_csv(file_path)
    df["len"] = [len(x) for x in df["sequence"]]
    df.sort_values(["len"], ascending=False, inplace=True)
    df.reset_index(inplace=True)
    return df


@dataclass(frozen=True, order=True)
class SequenceInfo:
    name: str
    sequence: str
    code: str


def find_common_subsequence(common_seqs: pd.DataFrame, seqs: List[str]) -> SequenceInfo:
    seqs = [to_dna(seq) for seq in seqs]
    saved_row = None
    for i, row in common_seqs.iterrows():
        fail = False
        for seq in seqs:
            if seq.find(row["sequence"]) == -1:
                fail = True
                break
        if fail:
            continue
        saved_row = row
        break
    if saved_row is None:
        return SequenceInfo("", "", "")
    else:
        return SequenceInfo(saved_row["name"], saved_row["sequence"], saved_row["code"])


def find_valid_subsequences(common_seqs: pd.DataFrame, seqs: List[str]) -> pd.DataFrame:
    seqs = [to_dna(seq) for seq in seqs]
    mask = [False for _ in range(len(common_seqs))]
    for i, row in common_seqs.iterrows():
        fail = False
        for seq in seqs:
            if seq.find(row["sequence"]) == -1:
                fail = True
                break
        if fail:
            continue
        mask[i] = True
    return common_seqs[mask]


def indentify_p5_sequence(seqs: List[str]) -> SequenceInfo:
    df = get_primer_dataframe(settings.RESOURCES_PATH + "p5_sequences.csv")
    return find_common_subsequence(df, seqs)


def indentify_fwd_primer(seqs: List[str]) -> SequenceInfo:
    df = get_primer_dataframe(settings.RESOURCES_PATH + "fwd_primers.csv")
    return find_common_subsequence(df, seqs)


def indentify_rev_primer(seqs: List[str]) -> SequenceInfo:
    df = get_primer_dataframe(settings.RESOURCES_PATH + "rev_primers.csv")
    seqs = [get_reverse_complement(seq, "DNA") for seq in seqs]
    return find_common_subsequence(df, seqs)


def find_valid_fwd_primers(seqs: List[str]) -> pd.DataFrame:
    df = get_primer_dataframe(settings.RESOURCES_PATH + "fwd_primers.csv")
    return find_valid_subsequences(df, seqs)


def find_valid_rev_primers(seqs: List[str]) -> pd.DataFrame:
    df = get_primer_dataframe(settings.RESOURCES_PATH + "rev_primers.csv")
    seqs = [get_reverse_complement(seq, "DNA") for seq in seqs]
    return find_valid_subsequences(df, seqs)


def get_p5_by_name(name: str):
    p5_sequences = pd.read_csv(settings.RESOURCES_PATH + "/p5_sequences.csv")
    row = p5_sequences[p5_sequences["name"] == name].iloc[0]
    p5 = structure.rna_structure(row["sequence"], row["structure"])
    return structure_set.get_single_struct_set(p5, structure_set.AddType.LEFT)


def get_p3_by_name(name: str):
    p3_sequences = pd.read_csv(settings.RESOURCES_PATH + "/p3_sequences.csv")
    row = p3_sequences[p3_sequences["name"] == name].iloc[0]
    p5 = structure.rna_structure(row["sequence"], row["structure"])
    return structure_set.get_single_struct_set(p5, structure_set.AddType.RIGHT)


def compute_edit_distance(df_result):
    scores = [100 for _ in range(len(df_result))]
    sequences = list(df_result["sequence"])
    for i, seq1 in enumerate(sequences):
        # if i % 10 == 0:
        #    print(i)
        for j, seq2 in enumerate(sequences):
            if i >= j:
                continue
            diff = editdistance.eval(seq1, seq2)
            if scores[i] > diff:
                scores[i] = diff
            if scores[j] > diff:
                scores[j] = diff
    avg = np.mean(scores)
    return avg


def random_wc_basepair():
    return random.choice(params.basepairs_wc)


def random_basepair():
    return random.choice(params.basepairs)


def random_gu_basepair():
    return random.choice(params.basepairs_gu)


def random_weighted_basepair():
    if random.randint(0, 1000) > 300:
        return random_wc_basepair()
    else:
        return random_gu_basepair()


def hamming(a, b):
    """hamming distance between two strings"""
    dist = 0
    for i, j in zip(a, b):
        if i != j:
            dist += 1
    return dist
