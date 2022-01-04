import pandas as pd
from enum import IntEnum

from rna_lib_design import structure, settings, logger

log = logger.setup_applevel_logger()


class AddType(IntEnum):
    HELIX = 0
    LEFT = 1
    RIGHT = 2

    def to_str(self):
        names = ("HELIX", "LEFT", "RIGHT")
        return names[int(self)]


def str_to_add_type(s: str) -> AddType:
    s = s.upper()
    if s == "HELIX":
        return AddType.HELIX
    elif s == "LEFT":
        return AddType.LEFT
    elif s == "RIGHT":
        return AddType.RIGHT
    else:
        log.error(f"invalid AddType: {s}")
        exit()


class StructureSet(object):
    def __init__(self, df, add_type):
        self.df = df.sample(frac=1).reset_index(drop=True)
        self.df["used"] = 0
        self.add_type = add_type
        self.current = -1

    def __len__(self):
        return len(self.df)

    def __get_return(self, pos):
        row = self.df.loc[pos]
        if self.add_type != AddType.HELIX:
            return [structure.Structure(row["seq"], row["ss"])]
        else:
            return [
                structure.Structure(row["seq_1"], row["ss_1"]),
                structure.Structure(row["seq_2"], row["ss_2"]),
            ]

    def remove_dots(self):
        self.df = self.df[~self.df["ss_1"].str.contains("\.")]
        self.df = self.df.reset_index(drop=True)

    def get_random(self):
        while 1:
            df = self.df.sample()
            row = df.iloc[0]
            i = df.index[0]
            if row["used"] != 0:
                continue
            self.current = i
            return self.__get_return(i)

    def get(self, pos):
        return self.__get_return(pos)

    def get_current_pos(self):
        return self.current

    def set_used(self, pos=None):
        if pos is None:
            self.df.at[self.current, "used"] = 1
        else:
            self.df.at[pos, "used"] = 1

    def apply_random(self, struct) -> structure.Structure:
        if self.add_type == AddType.LEFT:
            return self.get_random()[0] + struct
        elif self.add_type == AddType.RIGHT:
            return struct + self.get_random()[0]
        elif self.add_type == AddType.HELIX:
            s1, s2 = self.get_random()
            if len(struct) > 0:
                return s1 + struct + s2
            else:
                return s1 + structure.rna_structure_break() + s2

    def apply(self, struct, pos) -> structure.Structure:
        if self.add_type == AddType.LEFT:
            return self.get(pos)[0] + struct
        elif self.add_type == AddType.RIGHT:
            return struct + self.get(pos)[0]
        elif self.add_type == AddType.HELIX:
            s1, s2 = self.get(pos)
            if len(struct) > 0:
                return s1 + struct + s2
            else:
                return s1 + structure.rna_structure_break() + s2


class HairpinStructureSet(object):
    def __init__(self, loop, df, add_type):
        if not (add_type == AddType.LEFT or add_type == AddType.RIGHT):
            raise ValueError(f"incorrect type {add_type}")
        self.add_type = add_type
        self.loop = loop
        self.helices = StructureSet(df, AddType.HELIX)
        self.buffer = structure.rna_structure_unpaired("AAA")

    def __len__(self):
        return len(self.helices)

    def set_buffer(self, struct):
        self.buffer = struct

    def remove_dots(self):
        self.helices.remove_dots()

    def get(self, pos):
        s1, s2 = self.helices.get(pos)
        if self.buffer is not None:
            if self.add_type == AddType.LEFT:
                return [s1 + self.loop + s2 + self.buffer]
            else:
                return [self.buffer + s1 + self.loop + s2]
        else:
            return [s1 + self.loop + s2]

    def get_current_pos(self):
        return self.helices.get_current_pos()

    def get_random(self):
        s1, s2 = self.helices.get_random()
        if self.buffer is not None:
            if self.add_type == AddType.LEFT:
                return [s1 + self.loop + s2 + self.buffer]
            else:
                return [self.buffer + s1 + self.loop + s2]
        else:
            return [s1 + self.loop + s2]

    def set_used(self, pos=None):
        self.helices.set_used(pos)

    def apply_random(self, struct) -> structure.Structure:
        if self.add_type == AddType.LEFT:
            return self.get_random()[0] + struct
        elif self.add_type == AddType.RIGHT:
            return struct + self.get_random()[0]

    def apply(self, struct, pos) -> structure.Structure:
        if self.add_type == AddType.LEFT:
            return self.get(pos)[0] + struct
        elif self.add_type == AddType.RIGHT:
            return struct + self.get(pos)[0]


class SingleStructureSet(StructureSet):
    def __init__(self, df, add_type):
        super().__init__(df, add_type)

    # redefine set_used() so it never uses up the one entry
    def set_used(self, pos=None):
        pass


def get_single_struct_set(struct, add_type):
    df = pd.DataFrame(columns="seq ss".split())
    df.loc[0] = [str(struct.sequence), str(struct.dot_bracket)]
    return SingleStructureSet(df, add_type)


def apply(struct_sets, struct):
    temp_struct = struct
    for sd in list(struct_sets):
        temp_struct = sd.apply_random(temp_struct)
    return temp_struct


def get_optimal_helix_set(length, min_count=10):
    fname = settings.RESOURCES_PATH + "/barcodes/helices.csv"
    df = pd.read_csv(fname)
    df = df[df["length"] == length]
    if len(df) == 0:
        raise ValueError(f"no helices available with length {length}")
    df = df[df["size"] > min_count]
    df = df.sort_values(["diff"], ascending=False)
    if len(df) == 0:
        raise ValueError(
            f"no helices available with length {length} with max_count {min_count}"
        )
    df_helix = pd.read_csv(
        settings.RESOURCES_PATH + "barcodes/" + df.iloc[0]["path"]
    )
    return StructureSet(df_helix, AddType.HELIX)


def get_optimal_hairpin_set(
    length, loop_struct, min_count=10, type=AddType.LEFT
):
    h_df = get_optimal_helix_set(length, min_count).df
    return HairpinStructureSet(loop_struct, h_df, type)


def get_optimal_sstrand_set(length, min_count=10, type=AddType.LEFT):
    fname = settings.RESOURCES_PATH + "/barcodes/sstrand.csv"
    df = pd.read_csv(fname)
    df = df[df["length"] == length]
    if len(df) == 0:
        raise ValueError(f"no sstrand available with length {length}")
    df = df[df["size"] > min_count]
    df = df.sort_values(["diff"], ascending=False)
    if len(df) == 0:
        raise ValueError(
            f"no sstrand available with length {length} with max_count {min_count}"
        )
    df_ss = pd.read_csv(
        settings.RESOURCES_PATH + "barcodes/" + df.iloc[0]["path"]
    )
    return StructureSet(df_ss, type)


def get_common_seq_structure_set(name, type=AddType.LEFT):
    struct = structure.get_common_struct(name)
    return get_single_struct_set(struct, type)


def get_tail_structure_set():
    return get_common_seq_structure_set("rt_tail", AddType.RIGHT)
