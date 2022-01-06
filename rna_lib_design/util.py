import pandas as pd
import random

from rna_lib_design import params


def check_cols_in_dataframe(df, cols):
    for c in cols:
        if c not in df:
            raise ValueError("column {} must be present in dataframe".format(c))


def seq_to_dna(seq):
    new_seq = ""
    for e in seq:
        if e == "U":
            new_seq += "T"
        else:
            new_seq += e
    return new_seq


def seq_to_rna(seq):
    new_seq = ""
    for e in seq:
        if e == "T":
            new_seq += "U"
        else:
            new_seq += e
    return new_seq


def random_wc_basepair():
    return random.choice(params.basepairs_wc)


def random_basepair():
    return random.choice(params.basepairs)


def random_gu_basepair():
    return random.choice(params.basepairs_gu)


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
        bps.append(random_basepair())
    random.shuffle(bps)
    for bp in bps:
        seq_1 += bp[0]
        seq_2 = bp[1] + seq_2
    return [seq_1, seq_2]


def num_of_basepairs(self, ss):
    pass
