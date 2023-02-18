import pandas as pd
import numpy as np
from typing import Dict, Tuple, Optional
from pathlib import Path
from multiprocessing import Pool

from rna_lib_design.design import (
    get_best_designs_in_dataframe,
    write_results_to_file,
    DesignOptions,
)
from rna_lib_design import logger, structure_set, util, design
from rna_lib_design.structure_set import StructureSet, AddType

log = logger.setup_applevel_logger()


class Barcoder(object):
    def __init__(self, btype):
        self.btype: str = btype
        self.p5: Optional[StructureSet] = None
        self.p3: Optional[StructureSet] = None
        self.p5_buffer: Optional[StructureSet] = None
        self.p3_buffer: Optional[StructureSet] = None
        self.max: int = 1000000
        self.multiprocess = -1

    def set_p5_and_p3(
        self, p5: Optional[StructureSet], p3: Optional[StructureSet]
    ):
        self.p5, self.p3 = p5, p3

    def set_p5_buffer(self, p5_buffer: Optional[StructureSet]):
        self.p5_buffer = p5_buffer

    def set_p3_buffer(self, p3_buffer: Optional[StructureSet]):
        self.p3_buffer = p3_buffer

    def _barcode_w_sets(self, df, sets, design_opts) -> pd.DataFrame:
        add_ons = [self.p5_buffer, self.p3_buffer, self.p5, self.p3]
        for ao in add_ons:
            if ao is not None:
                sets.append(ao)
        if self.multiprocess == -1:
            return get_best_designs_in_dataframe(df, sets, design_opts)
        else:
            p = Pool(processes=self.multiprocess)
            spls = [e.split(self.multiprocess) for e in sets]
            set_spls = []
            for i in range(self.multiprocess):
                split_set = []
                for j in range(len(sets)):
                    split_set.append(spls[j][i])
                set_spls.append(split_set)
            dfs = np.array_split(df, self.multiprocess)
            all_design_opts = [design_opts for x in range(self.multiprocess)]
            output_dfs = p.starmap(
                get_best_designs_in_dataframe, zip(dfs, set_spls, all_design_opts)
            )
            return pd.concat(output_dfs)


class SingleBarcoder(Barcoder):
    def __init__(self, btype, length, add_type=AddType.RIGHT, loop=None):
        super().__init__(btype)
        self.length = length
        self.add_type = add_type
        self.loop: Optional[Structure] = loop

    def set_loop(self, loop: Structure):
        self.loop = loop

    def barcode(
        self, df: pd.DataFrame, design_opts: DesignOptions
    ) -> pd.DataFrame:
        if self.btype == "helix":
            v_set = structure_set.get_optimal_helix_set(
                self.length, len(df) * 1.2
            )
        elif self.btype == "sstrand":
            v_set = structure_set.get_optimal_sstrand_set(
                self.length, len(df) * 1.2, self.add_type
            )
        elif self.btype == "hairpin":
            if self.loop is None:
                log.error(
                    "cannot generate a hairpin barcode without supplying loop"
                    " sequence"
                )
                exit()
            # loop_struct = defaults.get_loop_from_str(self.loop)
            v_set = structure_set.get_optimal_hairpin_set(
                self.length, self.loop, len(df) * 1.2, self.add_type
            )
        else:
            log.error(f"{self.btype} is not a valid type")
            exit()
        return self._barcode_w_sets(df, [v_set], design_opts)


class CustomSingleBarcoder(Barcoder):
    def __init__(self, btype, bfile, add_type=AddType.RIGHT, loop=None):
        super().__init__(btype)
        self.barcode_file = bfile
        self.add_type = add_type
        self.loop: Optional[Structure] = loop

    def set_loop(self, loop: Structure):
        self.loop = loop

    def barcode(
        self, df: pd.DataFrame, design_opts: DesignOptions
    ) -> pd.DataFrame:
        df_set = pd.read_csv(self.barcode_file)
        if self.btype == "helix":
            v_set = structure_set.StructureSet(df_set, AddType.HELIX)
        elif self.btype == "sstrand":
            v_set = structure_set.StructureSet(df_set, self.add_type)
        elif self.btype == "hairpin":
            if self.loop is None:
                log.error(
                    "cannot generate a hairpin barcode without supplying loop"
                    " sequence"
                )
                exit()
            v_set = structure_set.HairpinStructureSet(
                self.loop, df_set, self.add_type
            )
        else:
            log.error(f"{self.btype} is not a valid type")
            exit()
        return self._barcode_w_sets(df, [v_set], design_opts)


class DoubleBarcode(Barcoder):
    def __init__(self, btype, lengths, loop):
        super().__init__(btype)
        self.lengths = lengths
        self.loop = loop

    def set_loop(self, loop: Structure):
        self.loop = loop

    def barcode(
        self, df: pd.DataFrame, design_opts: DesignOptions
    ) -> pd.DataFrame:
        sets = []
        if self.btype == "helix_hairpin":
            h_set = structure_set.get_optimal_helix_set(
                self.lengths[0], len(df) * 1.2
            )
            hp_set = structure_set.get_optimal_hairpin_set(
                self.lengths[1], self.loop, len(df) * 1.2, AddType.RIGHT
            )
            sets.extend([h_set, hp_set])
        elif self.btype == "hairpin_helix":
            h_set = structure_set.get_optimal_helix_set(
                self.lengths[1], len(df) * 1.2
            )
            hp_set = structure_set.get_optimal_hairpin_set(
                self.lengths[0],
                self.loop,
                len(df) * 1.2,
                AddType.LEFT,
            )
            sets.extend([h_set, hp_set])
        elif self.btype == "hairpin_hairpin":
            hp_set_1 = structure_set.get_optimal_hairpin_set(
                self.lengths[0],
                self.loop,
                len(df) * 1.2,
                AddType.LEFT,
            )
            hp_set_2 = structure_set.get_optimal_hairpin_set(
                self.lengths[1],
                self.loop,
                len(df) * 1.2,
                AddType.RIGHT,
            )
            sets.extend([hp_set_1, hp_set_2])
        else:
            log.error(f"unknown type {self.btype}")
            exit()
        return self._barcode_w_sets(df, sets, design_opts)


class CustomDoubleBarcode(Barcoder):
    def __init__(self, btype, bfiles, loop):
        super().__init__(btype)
        self.bfiles = bfiles
        self.loop = loop

    def set_loop(self, loop: Structure):
        self.loop = loop

    def barcode(
        self, df: pd.DataFrame, design_opts: DesignOptions
    ) -> pd.DataFrame:
        barcode_files = self.bfiles.split(",")
        if len(barcode_files) > 2:
            raise ValueError(
                f"too many barcode files were supplied {barcode_files}"
            )
        barcode_file_1, barcode_file_2 = "", ""
        if len(barcode_files) == 1:
            barcode_file_1 = barcode_files[0]
            barcode_file_2 = barcode_files[0]
        else:
            barcode_file_1 = barcode_files[0]
            barcode_file_2 = barcode_files[1]
        df_1 = pd.read_csv(barcode_file_1)
        df_2 = pd.read_csv(barcode_file_2)
        sets = []
        if self.btype == "helix_hairpin":
            h_set = structure_set.StructureSet(df_1, AddType.HELIX)
            hp_set = structure_set.HairpinStructureSet(
                self.loop, df_2, AddType.RIGHT
            )
            sets.extend([h_set, hp_set])
        elif self.btype == "hairpin_helix":
            h_set = structure_set.StructureSet(df_2, AddType.HELIX)
            hp_set = structure_set.HairpinStructureSet(
                self.loop, df_1, AddType.LEFT
            )
            sets.extend([h_set, hp_set])
        elif self.btype == "hairpin_hairpin":
            hp_set_1 = structure_set.HairpinStructureSet(
                self.loop, df_1, AddType.LEFT
            )
            hp_set_2 = structure_set.HairpinStructureSet(
                self.loop, df_1, AddType.RIGHT
            )
            sets.extend([hp_set_1, hp_set_2])
        else:
            log.error(f"unknown type {self.btype}")
            exit()
        return self._barcode_w_sets(df, sets, design_opts)


class TripleBarcode(Barcoder):
    def __init__(self, btype, lengths, loop):
        super().__init__(btype)
        self.lengths = lengths
        self.loop = loop

    def set_loop(self, loop: Structure):
        self.loop = loop

    def barcode(
        self, df: pd.DataFrame, design_opts: DesignOptions
    ) -> pd.DataFrame:
        sets = []
        if self.btype == "hairpin_helix_hairpin":
            hp_set_1 = structure_set.get_optimal_hairpin_set(
                self.lengths[0],
                self.loop,
                len(df) * 1.2,
                structure_set.AddType.LEFT,
            )
            h_set = structure_set.get_optimal_helix_set(
                self.lengths[1], len(df) * 1.2
            )
            hp_set_2 = structure_set.get_optimal_hairpin_set(
                self.lengths[2],
                self.loop,
                len(df) * 1.2,
                structure_set.AddType.RIGHT,
            )
            sets.extend([h_set, hp_set_1, hp_set_2])
        else:
            log.error(f"unknown type {self.btype}")
            exit()

        return self._barcode_w_sets(df, sets, design_opts)


def output_results(df_result, kwargs):
    org_total = len(df_result)
    df_result = df_result.dropna()
    diff = org_total - len(df_result)
    log.info(f"{diff} constructs were discarded as they were below cutoffs")
    util.compute_edit_distance(df_result)
    design.write_results_to_file(
        df_result, f"{kwargs['name']}/results", kwargs["name"]
    )
