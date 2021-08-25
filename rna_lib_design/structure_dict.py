import pandas as pd
import random

from rna_lib_design import structure, settings


class StructureDict(object):
    def __init__(self):
        pass

    def get_next(self):
        pass

    def last_was_used(self):
        pass

    def set_used(self, i):
        pass

    def is_done(self):
        pass

    def list_used(self):
        pass

    def apply_next(self, struct) -> structure.Structure:
        pass


class SStrandDict(StructureDict):
    def __init__(self, fname, type="LEFT"):
        self.type = type
        self.df = pd.read_csv(fname)
        self.df = self.df.sample(frac=1).reset_index(drop=True)
        self.last = -1
        self.used = {}
        self.trim = 0

    def __len__(self):
        return len(self.df)

    def set_trim(self, trim):
        self.trim = trim

    def get_next(self) -> structure.Structure:
        count = 0
        while 1:
            count += 1
            if count > len(self.df) * 2:
                print("cannot get_next! exiting")
                exit()
            i = random.randrange(1, len(self.df))
            row = self.df.loc[i]
            if row["seq"] in self.used:
                continue
            seq, ss = row["seq"], row["ss"]
            if self.trim != 0:
                seq = seq[0:-self.trim]
                ss = ss[0:-self.trim]
            s = structure.rna_structure(seq, ss)
            self.last = row["seq"]
            break
        return s

    def last_was_used(self):
        self.used[self.last] = 1

    def set_used(self, i):
        self.used[i] = 1

    def list_used(self):
        s = ""
        for i in self.used:
            s += self.df.loc[i]["seq"] + ","
        return s

    def is_done(self):
        if len(self.used) == len(self.df):
            return 1

    def apply_next(self, struct) -> structure.Structure:
        if self.type == "LEFT":
            return self.get_next() + struct
        else:
            return struct + self.get_next()


class HelixDict(StructureDict):
    def __init__(self, fname):
        self.df = pd.read_csv(fname)
        self.df = self.df.sample(frac=1).reset_index(drop=True)
        self.last = -1
        self.used = {}
        self.trim_structs = False
        self.trim = None

    def __len__(self):
        return len(self.df)

    def set_trim(self, trim):
        if trim > 0:
            self.trim_structs = True
            self.trim = trim
        else:
            self.trim_structs = False

    def get_next(self):
        count = 0
        while 1:
            count += 1
            if count > len(self.df) * 2:
                print("cannot get_next! exiting")
                exit()
            i = random.randrange(1, len(self.df))
            if i in self.used:
                continue
            row = self.df.loc[i]
            if not self.trim_structs:
                s1 = structure.rna_structure(row["seq_1"], row["ss_1"])
                s2 = structure.rna_structure(row["seq_2"], row["ss_2"])
            else:
                new_len = len(row["seq_1"]) - self.trim
                s1 = structure.rna_structure(row["seq_1"][:new_len], row["ss_1"][:new_len])
                s2 = structure.rna_structure(row["seq_2"][self.trim:], row["ss_2"][self.trim:])
            self.last = i
            break
        return [s1, s2]

    def last_was_used(self):
        self.used[self.last] = 1

    def set_used(self, i):
        self.used[i] = 1

    def list_used(self):
        s = ""
        for i in self.used:
            row = self.df.loc[i]
            s += row["seq_1"] + "&" + row["seq_2"] + ","
        return s

    def is_done(self):
        if len(self.used) == len(self.df):
            return 1

    def apply_next(self, struct) -> structure.Structure:
        s1, s2 = self.get_next()
        return s1 + struct + s2


class HairpinDict(StructureDict):
    def __init__(
            self, helix_dict, loop, type="LEFT", buffer=structure.rna_structure_unpaired('AAA')):
        self.helix_dict = helix_dict
        self.buffer = buffer
        self.type = type
        self.loop = loop

    def set_trim(self, trim):
        self.helix_dict.set_trim(trim)

    def get_next(self):
        while 1:
            h1, h2 = self.helix_dict.get_next()
            hp = h1 + self.loop + h2
            if self.buffer is not None:
                if self.type == 'LEFT':
                    hp += self.buffer
                else:
                    hp = self.buffer + hp
            return hp

    def last_was_used(self):
        self.helix_dict.last_was_used()

    def set_used(self, i):
        self.helix_dict.set_used(i)

    def list_used(self):
        return self.helix_dict.list_used()

    def is_done(self):
        return self.helix_dict.is_done()

    def apply_next(self, struct) -> structure.Structure:
        if self.type == "LEFT":
            return self.get_next() + struct
        else:
            return struct + self.get_next()

    @property
    def last(self):
        return self.helix_dict.last


class SingleDict(StructureDict):
    def __init__(self, struct, type="LEFT"):
        self.struct = struct
        self.type = type
        self.last = 1

    def get_next(self) -> structure.Structure:
        return self.struct

    def apply_next(self, struct) -> structure.Structure:
        if self.type == "LEFT":
            return self.get_next() + struct
        else:
            return struct + self.get_next()


class SingleHelixDict(StructureDict):
    def __init__(self, struct1, struct2):
        self.struct1 = struct1
        self.struct2 = struct2
        self.trim_structs = False
        self.trim = None

    def set_trim(self, trim):
        self.trim_structs = True
        self.trim = trim

    def get_next(self):
        if not self.trim_structs:
            return self.struct1, self.struct2
        else:
            new_len = len(self.struct1) - self.trim
            s1 = self.struct1[:new_len]
            s2 = self.struct2[self.trim:]
            return s1, s2

    def apply_next(self, struct) -> structure.Structure:
        s1, s2 = self.get_next()
        return s1 + struct + s2


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
