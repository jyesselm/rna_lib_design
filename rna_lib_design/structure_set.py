import pandas as pd
import random
from enum import IntEnum

from rna_lib_design import structure, settings


class AddType(IntEnum):
    HELIX = 0
    LEFT = 1
    RIGHT = 2

    def to_str(self):
        names = ("HELIX", "LEFT", "RIGHT")
        return names[int(self)]


class StructureSet(object):
    def __init__(self, df, add_type):
        self.df = df.sample(frac=1).reset_index(drop=True)
        self.df['used'] = 0
        self.add_type = add_type
        self.current = -1

    def __len__(self):
        return len(self.df)

    def __get_return(self, pos):
        row = self.df.loc[pos]
        if self.add_type != AddType.HELIX:
            return [structure.Structure(row['seq'], row['ss'])]
        else:
            return [
                structure.Structure(row['seq_1'], row['ss_1']),
                structure.Structure(row['seq_2'], row['ss_2'])
            ]

    def get_random(self):
        while 1:
            df = self.df.sample()
            row = df.iloc[0]
            i = df.index[0]
            if row['used'] != 0:
                continue
            self.current = i
            return self.__get_return(i)

    def get(self, pos):
        return self.__get_return(pos)

    def get_current_pos(self):
        return self.current

    def set_used(self, pos=None):
        if pos is None:
            self.df.at[self.current, 'used'] = 1
        else:
            self.df.at[pos, 'used'] = 1

    def apply_random(self, struct) -> structure.Structure:
        if self.add_type == AddType.LEFT:
            return self.get_random()[0] + struct
        elif self.add_type == AddType.RIGHT:
            return struct + self.get_random()[0]
        elif self.add_type == AddType.HELIX:
            s1, s2 = self.get_random()
            return s1 + struct + s2

    def apply(self, struct, pos) -> structure.Structure:
        if self.add_type == AddType.LEFT:
            return self.get(pos)[0] + struct
        elif self.add_type == AddType.RIGHT:
            return struct + self.get(pos)[0]
        elif self.add_type == AddType.HELIX:
            s1, s2 = self.get(pos)
            return s1 + struct + s2


class HairpinStructureSet(StructureSet):
    def __init__(self, loop, df, add_type):
        super().__init__(df, add_type)
        self.loop = loop
        self.buffer = structure.rna_structure_unpaired('AAA')

    def set_buffer(self, struct):
        self.buffer = struct

    def get_random(self):
        s1, s2 = super().get_random()
        return s1 + self.loop + s2 + self.buffer


class SingleStructureSet(StructureSet):
    def __init__(self, df, add_type):
        super().__init__(df, add_type)

    # redefine set_used() so it never uses up the one entry
    def set_used(self, pos=None):
        pass


"""

def apply(struct_dicts, struct):
    struct_dicts = struct_dicts[::-1]
    temp_struct = struct
    for sd in struct_dicts:
        temp_struct = sd.apply_next(temp_struct)
    return temp_struct


def get_helices(length, max_count=200):
    fname = settings.RESOURCES_PATH + "/barcodes/helices.csv"
    df = pd.read_csv(fname)
    df = df[df["length"] == length]
    if len(df) == 0:
        raise ValueError(f"no helices available with length {length}")
    df = df.sort_values(["count"])
    data_fname = None
    for i, row in df.iterrows():
        if max_count < row["count"]:
            data_fname = row["path"]
            break
    if data_fname is None:
        raise ValueError(f"no helices available with length {length} with max_count {max_count}")
    return HelixDict(data_fname)


def get_hairpins(h_length, loop, max_count=200, type="LEFT"):
    helix_dict = get_helices(h_length, max_count)
    return HairpinDict(helix_dict, loop, type=type)


def get_sstrands(length, max_count=200, type="LEFT"):
    r_min, r_max = 5, 29
    fname = settings.RESOURCES_PATH + "/barcodes/sstrand.csv"
    df = pd.read_csv(fname)
    df = df[df["length"] == length]
    if len(df) == 0:
        raise ValueError(f"no helices available with length {length}")
    df = df.sort_values(["count"])
    data_fname = None
    for i, row in df.iterrows():
        if max_count < row["count"]:
            data_fname = row["path"]
            break
    if data_fname is None:
        raise ValueError(f"no helices available with length {length} with max_count {max_count}")
    return SStrandDict(data_fname, type=type)


def get_common_seq(name, direction="LEFT"):
    common_structs = structure.common_structures()
    return SingleDict(common_structs[name], type=direction)


def get_tail():
    return get_common_seq('rt_tail', direction='RIGHT')


def get_p5(name):
    return get_common_seq('5PRIME', name)
"""
