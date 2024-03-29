from typing import List, Dict
import re
import pandas as pd
import numpy as np
from numpy import random
from dataclasses import dataclass

from seq_tools import SequenceStructure

from rna_lib_design.logger import get_logger
from rna_lib_design.settings import get_resources_path

log = get_logger("SSET")


def split_into_n(lst, n):
    """
    Split a list into n sublists with sizes as close to equal as possible.

    Args:
    - lst (list): The list to be split.
    - n (int): The number of sublists.

    Returns:
    - list of lists: The split sublists.
    """
    avg = len(lst) // n
    remainder = len(lst) % n
    result = []
    idx = 0

    for i in range(n):
        size = avg + (1 if i < remainder else 0)
        result.append(lst[idx : idx + size])
        idx += size

    return result


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
                [int(j) for j in i if j] for i in re.findall(r"(\d+),?(?:-(\d+))?", x)
            )
        ),
        [],
    )


class SequenceStructureSet:
    """
    A set of SequenceStructures that can be used to build up
    new sequences.
    """

    def __init__(self, seqstructs: List[SequenceStructure]):
        if len(seqstructs) == 0:
            raise ValueError("seqstructs must have at least one SequenceStructure")
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
        for _, row in df.iterrows():
            seqstructs.append(SequenceStructure(row["sequence"], row["structure"]))
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
        seq_struct_set = SequenceStructureSet(self.seqstructs + other.seqstructs)
        seq_struct_set.used = self.used + other.used
        return seq_struct_set

    def get_random(self) -> SequenceStructure:
        if all(self.used):
            raise Exception("All SequenceStructures have been used.")
        if len(self.seqstructs) == 1:
            return self.seqstructs[0]
        count = 0
        while True:
            index = random.randint(0, len(self.seqstructs) - 1)
            if not self.used[index]:
                self.last = index
                return self.seqstructs[index]
            count += 1
            if count > 10000:
                raise ValueError("cannot find a random sequence structure")

    def set_used(self, sec_struct) -> None:
        index = self.seqstructs.index(sec_struct)
        if not self.allow_duplicates:
            self.used[index] = True

    def set_last_used(self) -> None:
        if self.last is not None:
            self.used[self.last] = True

    def split(self, num_sets: int):
        """
        Splits the SequenceStructureSet into a list of SequenceStructureSets.
        """

        # if just 1 then we need to just distribute to each split set
        if len(self.seqstructs) == 1:
            split_sets = []
            for i in range(num_sets):
                new_set = SequenceStructureSet([self.seqstructs[0]])
                new_set.allow_duplicates = True
                split_sets.append(new_set)
            return split_sets

        if num_sets > len(self.seqstructs):
            raise ValueError(
                "num_sets must be less than or equal to the number of sets"
            )
        random.shuffle(self.seqstructs)
        seq_struct_splits = split_into_n(self.seqstructs, num_sets)
        return [SequenceStructureSet(s) for s in seq_struct_splits]

    def num_used(self):
        return sum(self.used)

    def num_available(self):
        return len(self) - self.num_used()


class SequenceStructureSetParser:
    def __init__(self):
        pass

    def parse(self, num_seqs: int, params: Dict) -> Dict[str, SequenceStructureSet]:
        """
        Parses a dictionary of parameters into a dictionary of SequenceStructureSets.
        The name of the each dictionary key is the name of the SequenceStructureSet.
        :param num_seqs: The number of sequences of interest
        :param params: A dictionary of parameters
        :return: A dictionary of SequenceStructureSets
        """
        set_dict = {}
        for name, params in params.items():
            set_dict[name] = self.__parse_entry(num_seqs, name, params)
        return set_dict

    # TODO add other types of entries. Can do from a csv? Or list of sequences and
    # TODO and structures?
    def __parse_entry(
        self, num_seqs: int, name: str, params: Dict
    ) -> SequenceStructureSet:
        """
        Parses a single entry in the params dictionary.
        """
        if "m_type" in params:
            return self.__parse_by_type(num_seqs, name, params)
        elif "sequence" in params and "structure" in params:
            seq_struct = SequenceStructure(params["sequence"], params["structure"])
            return SequenceStructureSet.from_single(seq_struct)
        elif "name" in params:
            seq_struct = get_named_seq_struct(params["name"])
            log.info(f"{name} is using a named sequence/structure: {params['name']}")
            log.info(f"{name} -> {seq_struct}")
            return SequenceStructureSet.from_single(seq_struct)
        else:
            raise ValueError(
                "sequence and structure or name must be specified or m_type"
            )

    def __parse_by_type(self, num_seqs, name, params: Dict) -> SequenceStructureSet:
        m_type = params["m_type"].upper()
        if "length" not in params:
            raise ValueError("length must be specified with m_type")
        lengths = str_to_range(str(params["length"]))
        if m_type == "HELIX":
            return self.__parse_helix_type(name, num_seqs, lengths, params)
        elif m_type == "SSTRAND":
            return self.__parse_sstrand_type(name, num_seqs, lengths, params)
        elif m_type == "HAIRPIN":
            return self.__parse_hairpin_type(name, num_seqs, lengths, params)

    def __parse_helix_type(
        self, name, num_seqs, lengths, params: Dict
    ) -> SequenceStructureSet:
        log.info(f"{name} structure type is helix")
        gu = True
        if "gu" in params:
            gu = params["gu"]
        sets = None
        for length in lengths:
            if sets is None:
                sets = get_optimal_helix_set(length, num_seqs, gu=gu)
            else:
                sets = sets + get_optimal_helix_set(length, num_seqs, gu=gu)
        log.info(f"{name} has {len(sets)} helix structures")
        return sets

    def __parse_sstrand_type(
        self, name, num_seqs, lengths, params: Dict
    ) -> SequenceStructureSet:
        log.info(f"{name} structure type is sstrand")
        sets = None
        for length in lengths:
            if sets is None:
                sets = get_optimal_sstrand_set(length, num_seqs)
            else:
                sets += sets + get_optimal_sstrand_set(length, num_seqs)
        log.info(f"{name} has {len(sets)} sstrand structures")
        return sets

    def __parse_hairpin_type(
        self, name, num_seqs, lengths, params: Dict
    ) -> SequenceStructureSet:
        log.info(f"{name} structure type is hairpin")
        if "loop_seq" in params and "loop_ss" in params:
            seq_struct = SequenceStructure(params["loop_seq"], params["loop_ss"])
        elif "name" in params:
            seq_struct = get_named_seq_struct(params["name"])
        else:
            raise ValueError("loop_seq and loop_ss or name must be specified")
        buffer_5p = SequenceStructure("", "")
        buffer_3p = SequenceStructure("AAA", "...")
        if "buffer_5p_seq" in params:
            seq = params["buffer_5p_seq"]
            if "buffer_5p_ss" in params:
                struct = params["buffer_5p_struct"]
            else:
                struct = "." * len(seq)
            buffer_5p = SequenceStructure(seq, struct)
        if "buffer_3p_seq" in params:
            seq = params["buffer_3p_seq"]
            if "buffer_3p_ss" in params:
                struct = params["buffer_3p_ss"]
            else:
                struct = "." * len(seq)
            buffer_3p = SequenceStructure(seq, struct)
        gu = True
        if "gu" in params:
            gu = params["gu"]
        sets = None
        for length in lengths:
            if sets is None:
                sets = get_optimal_hairpin_set(
                    seq_struct,
                    length,
                    num_seqs,
                    gu=gu,
                    buffer_5p=buffer_5p,
                    buffer_3p=buffer_3p,
                )
            else:
                sets = sets + get_optimal_hairpin_set(
                    seq_struct,
                    length,
                    num_seqs,
                    gu=gu,
                    buffer_5p=buffer_5p,
                    buffer_3p=buffer_3p,
                )
        log.info(f"{name} has {len(sets)} hairpin structures")
        return sets


# get sets from csv files #############################################################


def get_optimal_set(path, length, min_count, **kwargs) -> str:
    df = pd.read_csv(path)
    df = df[df["length"] == length]
    if len(df) == 0:
        raise ValueError(f"no available with length {length} in {path}")
    df = df[df["size"] > min_count]
    df = df.sort_values(["diff"], ascending=False)
    if "gu" in kwargs and not kwargs["gu"]:
        df = df[df["gu"] == 0]
    if len(df) == 0:
        raise ValueError(
            f"no set available with length {length} with max_count {min_count}"
        )
    return df.iloc[0]["path"]


def get_optimal_helix_set(length, min_count, gu=True):
    fname = get_resources_path() / "barcodes/helices.csv"
    csv_path = get_optimal_set(fname, length, min_count, gu=gu)
    return SequenceStructureSet.from_csv(get_resources_path() / "barcodes" / csv_path)


def get_optimal_sstrand_set(length, min_count):
    fname = get_resources_path() / "barcodes/sstrand.csv"
    csv_path = get_optimal_set(fname, length, min_count)
    return SequenceStructureSet.from_csv(get_resources_path() / "barcodes" / csv_path)


def get_optimal_hairpin_set(
    seq_struct, length, min_count, gu=True, buffer_5p=None, buffer_3p=None
):
    if buffer_5p is None:
        buffer_5p = SequenceStructure("", "")
    if buffer_3p is None:
        buffer_3p = SequenceStructure("", "")
    fname = get_resources_path() / "barcodes/helices.csv"
    csv_path = get_optimal_set(fname, length, min_count, gu=gu)
    df = pd.read_csv(get_resources_path() / "barcodes" / csv_path)
    seqstructs = []
    for index, row in df.iterrows():
        h_seq_struct = SequenceStructure(row["sequence"], row["structure"])
        h_seq_structs = h_seq_struct.split_strands()
        seqstructs.append(
            buffer_5p + h_seq_structs[0] + seq_struct + h_seq_structs[1] + buffer_3p
        )
    return SequenceStructureSet(seqstructs)


# get seq_structs from dataframes #####################################################


def get_named_seq_structs():
    dir_path = get_resources_path() / "named_seqs/rna"
    csv_files = list(dir_path.glob("*.csv"))
    dfs = []
    for fname in csv_files:
        dfs.append(pd.read_csv(fname)[["name", "sequence", "structure"]])
    return pd.concat(dfs)


def get_named_seq_struct(name):
    df = get_named_seq_structs()
    df = df[df["name"] == name]
    if len(df) == 0:
        raise ValueError(f"no sequence structure with name {name}")
    return SequenceStructure(df.iloc[0]["sequence"], df.iloc[0]["structure"])
