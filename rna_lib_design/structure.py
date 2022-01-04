from __future__ import annotations

import pandas as pd
import collections as col

from rna_lib_design import settings


class Sequence(object):
    def __init__(self, sequence, seq_type="RNA"):
        self.__seq_type = seq_type
        self.__seq = sequence
        if self.__seq_type == "RNA":
            self.__seq = self.__convert_to_rna(self.__seq)
        else:
            self.__seq = self.__convert_to_dna(self.__seq)

    def __str__(self):
        return "".join(self.__seq)

    def __len__(self):
        return len(self.__seq)

    def __eq__(self, other):
        if type(other) == Sequence:
            if self.__seq_type != other.__seq_type:
                return False
            return self.seq == other.seq
        else:
            return str(self) == other

    def __add__(self, other):
        if type(other) == str:
            return str(self) + other

        if self.__seq_type != other.__seq_type:
            raise ValueError(
                "cannot add {} to {} they are not of the same sequence type".format(
                    self, other
                )
            )
        return Sequence(self.__seq + other.__seq, self.__seq_type)

    def __radd__(self, other):
        if type(other) == str:
            return other + str(self)
        else:
            raise TypeError(f"cannot add {type(other)} to Sequence")

    def __getitem__(self, item):
        return Sequence(self.__seq[item], self.__seq_type)

    def remove_segment(self, pos1, pos2):
        new_seq = list(self.__seq)
        del new_seq[pos1:pos2]
        return Sequence(new_seq, self.__seq_type)

    def remove_segments(self, bounds1, bounds2):
        return self.remove_segment(bounds2[0], bounds2[1]).remove_segment(
            bounds1[0], bounds1[1]
        )

    def insert(self, seq, pos) -> Sequence:
        new_seq = list(seq)
        if pos == 0:
            return Sequence(new_seq + self.__seq, self.__seq_type)
        elif pos == len(self) - 1:
            return Sequence(self.__seq + new_seq, self.__seq_type)

        seqs = list(self.__seq[0:pos])
        seqs.extend(new_seq)
        seqs.extend(list(self.__seq[pos:]))

        return Sequence(seqs, self.__seq_type)

    def insert_segments(self, seq1, seq2, pos1, pos2) -> Sequence:
        if pos1 >= pos2:
            raise ValueError("pos2 must be larger than pos1")
        return self.insert(seq2, pos2).insert(seq1, pos1)

    def is_rna(self) -> bool:
        return self.__seq_type == "RNA"

    def is_dna(self) -> bool:
        return self.__seq_type == "DNA"

    def to_rna(self) -> Sequence:
        return Sequence(list(self.__seq), "RNA")

    def to_dna(self) -> Sequence:
        return Sequence(list(self.__seq), "DNA")

    def get_complement(self) -> Sequence:
        return Sequence(self.__get_complement(self.__seq), self.__seq_type)

    def get_dna_complement(self) -> Sequence:
        return Sequence(self.__get_complement(self.__seq), "DNA")

    def get_rna_complement(self) -> Sequence:
        return Sequence(self.__get_complement(self.__seq), "RNA")

    def get_reverse_complement(self) -> Sequence:
        return Sequence(self.__get_complement(self.__seq)[::-1], self.__seq_type)

    def get_dna_reverse_complement(self) -> Sequence:
        return Sequence(self.__get_complement(self.__seq)[::-1], "DNA")

    def get_rna_reverse_complement(self) -> Sequence:
        return Sequence(self.__get_complement(self.__seq)[::-1], "RNA")

    def str(self) -> str:
        return str(self)

    # private
    def __convert_to_dna(self, seq):
        new_seq = []
        for s in seq:
            if s == "U":
                new_seq.append("T")
            else:
                new_seq.append(s)
        return new_seq

    def __convert_to_rna(self, seq):
        new_seq = []
        for s in seq:
            if s == "T":
                new_seq.append("U")
            else:
                new_seq.append(s)
        return new_seq

    def __get_complement(self, seq):
        comp_list = ["A"] * len(seq)
        for i, e in enumerate(seq):
            if e == "C":
                comp_list[i] = "G"
            elif e == "G":
                comp_list[i] = "C"
            elif e == "A" and self.is_rna():
                comp_list[i] = "U"
            elif e == "A" and self.is_dna():
                comp_list[i] = "T"
            elif e == "U" and self.is_rna():
                comp_list[i] = "A"
            elif e == "T" and self.is_dna():
                comp_list[i] = "A"
            else:
                raise ValueError(
                    "invalid sequence character: {} in {}".format(e, self._seq_type)
                )
        return comp_list

    def __are_sequences_equal(self, seq1, seq2):
        if len(seq1) != len(seq2):
            return False
        for e1, e2 in zip(seq1, seq1):
            if e1 != e2:
                return False
        return True


def rna_sequence(sequence_str) -> Sequence:
    return Sequence(list(sequence_str), seq_type="RNA")


def dna_sequence(sequence_str) -> Sequence:
    return Sequence(list(sequence_str), seq_type="DNA")


bracket_left = {x: i for i, x in enumerate("([{<ABCDEFGHIJKLMNOPQRSTUVWXYZ")}
bracket_right = {x: i for i, x in enumerate(")]}>abcdefghijklmnopqrstuvwxyz")}


class DotBracket(object):
    def __init__(self, dot_bracket_str):
        self.__is_valid = 1
        self.__dot_bracket = list(dot_bracket_str)
        self.__pairs = self.__assign_pairs()

    def __str__(self):
        return "".join(self.__dot_bracket)

    def __len__(self):
        return len(self.__dot_bracket)

    def __add__(self, other):
        if type(other) == str:
            return str(self) + other
        return DotBracket(self.__dot_bracket + other.__dot_bracket)

    def __radd__(self, other):
        if type(other) == str:
            return other + str(self)
        else:
            raise TypeError(f"cannot add {type(other)} to DotBracket")

    def __eq__(self, other):
        if type(other) == DotBracket:
            return self.difference(other) == 0
        else:
            return str(self) == other

    def insert(self, db, pos) -> DotBracket:
        new_db = list(db)
        if pos == 0:
            return DotBracket(new_db + self.__dot_bracket)
        elif pos == len(self) - 1:
            return DotBracket(self.__dot_bracket + new_db)

        dbs = list(self.__dot_bracket[0:pos])
        dbs.extend(list(new_db))
        dbs.extend(list(self.__dot_bracket[pos:]))

        return DotBracket(dbs)

    def insert_segments(self, seq1, seq2, pos1, pos2) -> DotBracket:
        if pos1 >= pos2:
            raise ValueError("pos2 must be larger than pos1")
        return self.insert(seq2, pos2).insert(seq1, pos1)

    def __getitem__(self, item):
        return DotBracket(self.__dot_bracket[item])

    def remove_segment(self, pos1, pos2):
        new_db = list(self.__dot_bracket)
        del new_db[pos1:pos2]
        return DotBracket(new_db)

    def remove_segments(self, bounds1, bounds2):
        return self.remove_segment(bounds2[0], bounds2[1]).remove_segment(
            bounds1[0], bounds1[1]
        )

    def __assign_pairs(self):
        i = 0
        stack = col.defaultdict(list)
        pairs = {}
        a = ""
        for a in self.__dot_bracket:
            if a == "&":
                continue
            i += 1
            if a == ".":
                continue
            if a in bracket_left:
                stack[bracket_left[a]].append(i)
            else:
                if len(stack[bracket_right[a]]) == 0:
                    self.__is_valid = 0
                    return pairs
                j = stack[bracket_right[a]].pop()
                pairs[i - 1] = j - 1
                pairs[j - 1] = i - 1

        if a in bracket_left:
            if len(stack[bracket_left[a]]) != 0:
                self.__is_valid = 0
        return pairs

    def is_valid(self):
        return self.__is_valid

    def is_paired(self, pos) -> bool:
        if pos in self.__pairs:
            return True
        else:
            return False

    def get_pair_partner(self, pos) -> int:
        if pos not in self.__pairs:
            raise ValueError("cannot get pair partner, position is not paired")
        return self.__pairs[pos]

    def difference(self, other: DotBracket) -> int:
        diff = 0
        for e1, e2 in zip(self.__dot_bracket, other.__dot_bracket):
            if e1 != e2:
                diff += 1
        return diff

    def __search(self, pos, seen):
        pass

    def get_helical_segments(self):
        if not self.__is_valid:
            raise ValueError("cannot assign motifs with in a invalid dotbracket")

        current = [[], []]
        for i in range(0, len(self)):
            if self.__dot_bracket[i] == "(" and len(current[0]) == 0:
                if len(current[0]) == 0:
                    current[0].append(i)
                    current[1].append(self.get_pair_partner(i))
                else:
                    diff1 = current[0][-1] - i
                    diff2 = current[1][-1] - self.get_pair_partner(i)
                    if diff1 == -1 and diff2 == 1:
                        current[0].append(i)
                        current[1].append(self.get_pair_partner(i))
                        continue

        exit()
        all_helices = []
        for i in range(0, len(self)):
            if self.__dot_bracket[i] == "(":
                current[0].append(i)
                current[1].append(self.get_pair_partner(i))
            elif len(current[0]) > 0:
                all_helices.append(current)
                current = [[], []]

        if len(current[0]) > 0:
            all_helices.append(current)
        return all_helices


class Structure(object):
    def __init__(self, sequence: Sequence, dot_bracket: DotBracket):
        if len(sequence) != len(dot_bracket):
            raise ValueError(
                "sequence and dot bracket must be the same length: {} {}".format(
                    sequence, dot_bracket
                )
            )
        self.__sequence = sequence
        self.__dot_bracket = dot_bracket

    def __len__(self):
        return len(self.__sequence)

    def __str__(self):
        return str(self.__sequence) + " " + str(self.__dot_bracket)

    def __add__(self, other):
        return Structure(
            self.__sequence + other.__sequence,
            self.__dot_bracket + other.__dot_bracket,
        )

    def __eq__(self, other: Structure):
        return self.sequence == other.sequence and self.dot_bracket == other.dot_bracket

    def __ne__(self, other: Structure):
        return not self == other

    def insert(self, seq, ss, pos):
        return Structure(
            self.__sequence.insert(seq, pos), self.__dot_bracket.insert(ss, pos)
        )

    def insert_segments(self, seq1, ss1, seq2, ss2, pos1, pos2):
        return Structure(
            self.__sequence.insert_segments(seq1, seq2, pos1, pos2),
            self.__dot_bracket.insert_segments(ss1, ss2, pos1, pos2),
        )

    def insert_bps(self, bp_str, pos):
        if not self.is_paired(pos):
            raise ValueError("must insert basepairs at an existing helix")
        spl = bp_str.split("&")
        seg1, seg2 = spl[0], spl[1][::-1]
        pair_pos = self.get_pair_partner(pos)
        return self.insert_segments(
            seg1, "(" * len(seg1), seg2, ")" * len(seg2), pos, pair_pos
        )

    def __getitem__(self, item):
        return Structure(self.__sequence[item], self.__dot_bracket[item])

    def remove_segment(self, pos1, pos2):
        return Structure(
            self.__sequence.remove_segment(pos1, pos2),
            self.__dot_bracket.remove_segment(pos1, pos2),
        )

    def remove_segments(self, bounds1, bounds2):
        return self.remove_segment(bounds2[0], bounds2[1]).remove_segment(
            bounds1[0], bounds1[1]
        )

    def remove_bp(self, pos):
        if not self.is_paired(pos):
            raise ValueError("cannot delete a bp at an existing bp")
        pair_pos = self.get_pair_partner(pos)
        if pos < pair_pos:
            return self.remove_segments([pos, pos + 1], [pair_pos, pair_pos + 1])
        else:
            return self.remove_segments([pair_pos, pair_pos + 1], [pos, pos + 1])

    def is_paired(self, pos):
        return self.__dot_bracket.is_paired(pos)

    def get_pair_partner(self, pos):
        return self.__dot_bracket.get_pair_partner(pos)

    def dot_bracket_difference(self, db):
        return self.__dot_bracket.difference(db)

    def is_dot_bracket_different(self, db):
        return self.__dot_bracket == db

    @property
    def sequence(self):
        return self.__sequence

    @property
    def dot_bracket(self):
        return self.__dot_bracket


def rna_structure(sequence_str, dot_bracket_str) -> Structure:
    return Structure(rna_sequence(sequence_str), DotBracket(dot_bracket_str))


def rna_structure_unpaired(sequence_str) -> Structure:
    return rna_structure(sequence_str, len(sequence_str) * ".")


def rna_structure_break() -> Structure:
    return rna_structure("&", "&")


def dna_structure(sequence_str, dot_bracket_str) -> Structure:
    return Structure(dna_sequence(sequence_str), DotBracket(dot_bracket_str))


def dna_structure_left_helix(sequence_str) -> Structure:
    return dna_structure(sequence_str, len(sequence_str) * "(")


def dna_structure_right_helix(sequence_str) -> Structure:
    return dna_structure(sequence_str, len(sequence_str) * ")")


def rna_structure_from_row(row) -> Structure:
    return rna_structure(row["sequence"], row["structure"])


def common_structure_dataframe():
    path = settings.RESOURCES_PATH + "/common_seqs.csv"
    df = pd.read_csv(path)
    return df


def common_structures():
    structures = {}
    df = common_structure_dataframe()
    for i, row in df.iterrows():
        structures[row["name"]] = rna_structure(row["sequence"], row["structure"])
    return structures


def get_common_struct(name):
    common_structs = common_structures()
    return common_structs[name]
