import pandas as pd
import numpy as np
import random
import editdistance
from dataclasses import dataclass
from typing import List

from seq_tools.sequence import convert_to_dna, get_reverse_complement
from rna_lib_design import params, logger, settings, structure, structure_set

log = logger.setup_applevel_logger()


@dataclass(frozen=True, order=True)
class SequenceInfo(object):
    name: str
    sequence: str
    code: str


def get_primer_dataframe(file_path: str) -> pd.DataFrame:
    df = pd.read_csv(file_path)
    df["len"] = [len(x) for x in df["sequence"]]
    df.sort_values(["len"], ascending=False, inplace=True)
    df.reset_index(inplace=True)
    return df


def find_common_subsequence(
    common_seqs: pd.DataFrame, seqs: List[str]
) -> SequenceInfo:
    seqs = [convert_to_dna(seq) for seq in seqs]
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
        return SequenceInfo(
            saved_row["name"], saved_row["sequence"], saved_row["code"]
        )


def find_valid_subsequences(
    common_seqs: pd.DataFrame, seqs: List[str]
) -> pd.DataFrame:
    seqs = [convert_to_dna(seq) for seq in seqs]
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


def max_stretch(s):
    """returns max stretch of the same letter in string"""
    max_stretch = 0
    last = None
    stretch = 1
    for n in s:
        if last == None:
            last = n
            stretch = 1
            continue
        if n == last:
            stretch += 1
            if stretch > max_stretch:
                max_stretch = stretch
        else:
            stretch = 1
        last = n
    if stretch > max_stretch:
        max_stretch = stretch
    return max_stretch


def max_gc_stretch(s1, s2):
    i = -1
    j = len(s2)
    max_count = 0
    count = 0
    while i < len(s1) - 1:
        i += 1
        j -= 1
        flag = 0
        if s1[i] == "G" and s2[j] == "C":
            flag = 1
        elif s1[i] == "C" and s2[j] == "G":
            flag = 1
        if flag:
            count += 1
            continue
        else:
            if count > max_count:
                max_count = count
            count = 0
    if count > max_count:
        max_count = count
    return max_count


def random_helix(length, gu=0):
    seq_1 = ""
    seq_2 = ""
    bps = []
    for i in range(0, gu):
        bps.append(random_gu_basepair())
    for i in range(0, length - gu):
        bps.append(random_wc_basepair())
    random.shuffle(bps)
    for bp in bps:
        seq_1 += bp[0]
        seq_2 = bp[1] + seq_2
    return [seq_1, seq_2]


def num_of_basepairs(self, ss):
    pass


@dataclass(frozen=True, order=True)
class Stretches:
    max_stretch_1: int
    max_stretch_2: int
    max_gc_stretch: int


def compute_stretches(seq1, seq2):
    return Stretches(
        max_stretch(seq1), max_stretch(seq2), max_gc_stretch(seq1, seq2)
    )
