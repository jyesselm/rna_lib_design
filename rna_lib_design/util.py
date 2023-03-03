import random
import pandas as pd
import numpy as np
import editdistance
from dataclasses import dataclass
from typing import List, Optional
from pathlib import Path

from seq_tools.sequence import to_dna, get_reverse_complement
from seq_tools.dataframe import has_5p_sequence, has_3p_sequence

from rna_lib_design.logger import get_logger
from rna_lib_design.settings import get_resources_path

log = get_logger("UTIL")


BASEPAIRS = ["AU", "UA", "GC", "CG", "GU", "UG"]
BASEPAIRS_WC = ["AU", "UA", "GC", "CG"]
BASEPAIRS_GU = ["GU", "UG"]


@dataclass(frozen=True, order=True)
class SequenceInfo:
    name: str
    sequence: str
    code: str


def get_seq_fwd_primer(df: pd.DataFrame) -> Optional[SequenceInfo]:
    """
    gets the sequence forward primer information
    :param df: the dataframe with sequences
    :return: the sequence forward primer information
    """
    df = df.copy()
    df = to_dna(df)
    path = get_resources_path() / "named_seqs/p5_sequences.csv"
    df_p5 = pd.read_csv(path)
    for _, row in df_p5.iterrows():
        if has_5p_sequence(df, row["sequence"]):
            return SequenceInfo(row["name"], row["sequence"], row["code"])
    return None


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
