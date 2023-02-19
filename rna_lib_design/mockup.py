from typing import List
import random
import re

from seq_tools import SequenceStructure
import pandas as pd

from rna_lib_design.settings import get_resources_path
from rna_lib_design.logger import get_logger

log = get_logger(__file__)


def str_to_range(x):
    """
    Convert a string representation of a range of numbers to a list of integers.

    Given a string representation of a range of numbers, this function returns a
    list of integers corresponding to the numbers in the range. The string can
    contain single numbers separated by commas, and ranges of numbers separated by
    a hyphen.

    :param x: A string representation of a range of numbers.
    :return: A list of integers corresponding to the numbers in the range.
    """
    return sum(
        (
            i if len(i) == 1 else list(range(i[0], i[1] + 1))
            for i in (
                [int(j) for j in i if j]
                for i in re.findall(r"(\d+),?(?:-(\d+))?", x)
            )
        ),
        [],
    )


def generate_helix(length, symbol):
    seq = symbol * length
    return SequenceStructure(seq + "&" + seq, "1" * length + "&" + "1" * length)


def parse_sequence_string(seq_str):
    """
    Parses a string of the form
    "P5-HPBARCODE-HBARCODE6A-SOI-HBARCODE6B-AC-P3" to return a dictionary
    containing the order and direction in which the other segments appear
    in reference to SOI.

    :param seq_str: The sequence string to parse.
    :type seq_str: str
    :return: A dictionary containing the order and direction in
    which the other segments appear in reference to SOI.
    :rtype: dict[str, int]
    """
    segments = seq_str.split("-")
    soi_index = segments.index("SOI")
    order = {}
    for i, segment in enumerate(segments):
        if i == soi_index:
            continue
        distance = i - soi_index
        key = segment.rstrip("AB")
        if segment.endswith("A"):
            if key in order:
                order[key] = min(distance, order[key])
            else:
                order[key] = distance
        elif segment.endswith("B"):
            if key in order:
                order[key] = max(distance, order[key])
            else:
                order[key] = distance
        else:
            if key in order:
                order[key] = min(distance, order[key])
            else:
                order[key] = distance
    return order


class SequenceStructureSet:
    """
    A set of SequenceStructures that can be used to build up
    new sequences.
    """

    def __init__(self, seqstructs: List[SequenceStructure]):
        self.seqstructs = seqstructs
        self.used = [False] * len(seqstructs)
        self.allow_duplicates = False
        self.last = None

    @classmethod
    def from_csv(cls, csv_path: str):
        """
        Creates a SequenceStructureSet from a csv file.
        """
        df = pd.read_csv(csv_path)
        seqstructs = []
        for index, row in df.iterrows():
            seqstructs.append(
                SequenceStructure(row["sequence"], row["structure"])
            )
        return cls(seqstructs)

    @classmethod
    def from_single(cls, seqstruct: SequenceStructure):
        """
        Creates a SequenceStructureSet from a single SequenceStructure.
        """
        sss = cls([seqstruct])
        sss.allow_duplicates = True
        return sss

    def __len__(self):
        return len(self.seqstructs)

    def __add__(self, other):
        seq_struct_set = SequenceStructureSet(
            self.seqstructs + other.seqstructs
        )
        seq_struct_set.used = self.used + other.used
        return seq_struct_set

    def get_random(self) -> SequenceStructure:
        if all(self.used):
            raise Exception("All SequenceStructures have been used.")
        while True:
            index = random.randint(0, len(self.seqstructs) - 1)
            if not self.used[index]:
                self.last = index
                return self.seqstructs[index]

    def set_used(self, sec_struct) -> None:
        index = self.seqstructs.index(sec_struct)
        if not self.allow_duplicates:
            self.used[index] = True

    def set_last_used(self) -> None:
        if self.last is not None:
            self.used[self.last] = True


def sequence_structure_set_from_params(num_seqs, params: List[str]):
    """
    Creates a SequenceStructureSet from a list of parameters.
    """
    if "m_type" in params:
        log.info(f"structure type is helix")
        m_type = params["m_type"].upper()
        if "length" not in params:
            raise ValueError("length must be specifiied with m_type")
        lengths = str_to_range(str(params["length"]))
        if m_type == "HELIX":
            gu = True
            if "gu" in params:
                gu = params["gu"]
            sets = None
            for length in lengths:
                if sets is None:
                    sets = get_optimal_helix_set(length, num_seqs, gu=gu)
                else:
                    sets = sets + get_optimal_helix_set(length, num_seqs, gu=gu)
            return sets
        elif m_type == "SSTRAND":
            sets = None
            for length in lengths:
                if sets is None:
                    sets = get_optimal_sstrand_set(length, num_seqs)
                else:
                    sets = sets + get_optimal_sstrand_set(length, num_seqs)
        elif m_type == "HAIRPIN":
            pass
        else:
            raise ValueError(f"unknown structure type {m_type}")


# get sets from csv files #############################################################


def get_optimal_set(path, length, min_count, **kwargs) -> str:
    df = pd.read_csv(path)
    df = df[df["length"] == length]
    if len(df) == 0:
        raise ValueError(f"no available with length {length} in {path}")
    df = df[df["size"] > min_count]
    df = df.sort_values(["diff"], ascending=False)
    if "gu" in kwargs and kwargs["gu"] == False:
        df = df[df["gu"] == 0]
    if len(df) == 0:
        raise ValueError(
            f"no set available with length {length} with max_count {min_count}"
        )
    return df.iloc[0]["path"]


def get_optimal_helix_set(length, min_count, gu=True):
    fname = get_resources_path() / "barcodes/helices.csv"
    csv_path = get_optimal_set(fname, length, min_count, gu=gu)
    return SequenceStructureSet.from_csv(
        get_resources_path() / "barcodes" / csv_path
    )


def get_optimal_sstrand_set(length, min_count):
    fname = get_resources_path() / "barcodes/sstrands.csv"
    csv_path = get_optimal_set(fname, length, min_count)
    return SequenceStructureSet.from_csv(
        get_resources_path() / "barcodes" / csv_path
    )


""" 
def unsued():
    symbols = [
        "1",
        "2",
        "3",
        "4",
        "5",
        "6",
        "7",
        "8",
        "9",
        "0",
        "*",
        "#",
        "$",
        "%",
        "^",
        "*",
        "B",
        "D",
        "F",
    ]
    replaces = []
    pos = 0
    return seqstruct
    for info in params:
        if info.startswith("HELIX-"):
            helix_len = int(info.split("-")[1])
            helix_strands = generate_helix(
                helix_len, symbols[pos]
            ).split_strands()
            seqstruct = helix_strands[0] + seqstruct + helix_strands[1]
    print(seqstruct)

    return replaces, seqstruct
"""
