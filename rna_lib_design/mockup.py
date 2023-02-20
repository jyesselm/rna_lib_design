from typing import List
import random
import re
import pandas as pd

from seq_tools import SequenceStructure
from vienna import fold

from rna_lib_design.structure_set import SequenceStructureSet
from rna_lib_design.settings import get_resources_path
from rna_lib_design.logger import get_logger

log = get_logger(__file__)

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


def parse_sequence_structure_sets(num_seqs, params):
    sets = {}
    for key, value in params.items():
        sets[key] = sequence_structure_set_from_params(num_seqs, value)
    return sets


def get_design_setup(num_seqs, build_str, params):
    sets = parse_sequence_structure_sets(num_seqs, params)
    segments = parse_sequence_string(build_str)
    pos = 1
    seen = []
    build_up = []
    while True:
        found = False
        for seg_name, seg_pos in segments.items():
            if pos != abs(seg_pos):
                continue
            if seg_name in seen:
                continue
            if seg_name not in sets:
                raise ValueError(f"no set for {seg_name}")
            found = True
            seen.append(seg_name)
            direction = "3PRIME"
            if seg_pos < 0:
                direction = "5PRIME"
            cur_set = sets[seg_name]
            # hacky way to check if helix
            if cur_set.get_random().sequence.count("&") > 0:
                direction = "HELIX"
            build_up.append((direction, seg_name, cur_set))
        if not found:
            pos += 1
        if len(seen) == len(segments):
            break

    return build_up


def design(soi_seq_struct, build_up):
    # build up template seq_struct
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
    seq_struct = soi_seq_struct
    iterating_sets = []
    for direction, seg_name, cur_set in build_up:
        symbol = symbols.pop(0)
        if direction == "HELIX":
            if len(cur_set) > 1:
                length = len(cur_set.get_random().split_strands()[0])
                seq = symbol * length + "&" + symbol * length
                h_seq_struct = SequenceStructure(seq, seq).split_strands()
                seq_struct = h_seq_struct[0] + seq_struct + h_seq_struct[1]
                iterating_sets.append(
                    [[symbol * length, symbol * length], cur_set]
                )
            else:
                h_seq_struct = cur_set.get_random().split_strands()
                seq_struct = h_seq_struct[0] + seq_struct + h_seq_struct[1]
        elif direction == "5PRIME":
            if len(cur_set) > 1:
                seq = symbol * len(cur_set.get_random())
                seq_struct = SequenceStructure(seq, seq) + seq_struct
                iterating_sets.append([[seq], cur_set])
            else:
                seq_struct = cur_set.get_random() + seq_struct
        elif direction == "3PRIME":
            if len(cur_set) > 1:
                seq = symbol * len(cur_set.get_random())
                seq_struct = seq_struct + SequenceStructure(seq, seq)
                iterating_sets.append([[seq], cur_set])
            else:
                seq_struct = seq_struct + cur_set.get_random()
        else:
            raise ValueError(f"unknown direction {direction}")

    attempts = 10
    best = []
    best_seq = ""
    best_ens_defect = 9999
    for i in range(0, attempts):
        sequence = seq_struct.sequence
        structure = seq_struct.structure
        used = []
        for replace_str, cur_set in iterating_sets:
            sol = cur_set.get_random()
            used.append(sol)
            strands = sol.split_strands()
            for rs, strand in zip(replace_str, strands):
                sequence = sequence.replace(rs, strand.sequence, 1)
                structure = structure.replace(rs, strand.structure, 1)
        r = fold(sequence)
        success = False
        if r.dot_bracket == structure:
            success = True
        if not success:
            # TODO should put in some logic here to not be too strict about barcode structure
            continue
        if r.ens_defect < best_ens_defect:
            best_ens_defect = r.ens_defect
            best_seq = sequence
            best = used
    if len(best) != 0:
        for sol, sss in zip(best, iterating_sets):
            sss[1].set_used(sol)
    print(best_seq)
    print(best_ens_defect)



# get sets from csv files #############################################################


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
