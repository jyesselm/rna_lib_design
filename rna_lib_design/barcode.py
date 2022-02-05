import click
import cloup
import pandas as pd
import shutil
import os
from typing import Dict, Tuple, Optional
from pathlib import Path

from rna_lib_design.design import (
    get_best_designs_in_dataframe,
    write_results_to_file,
    DesignOptions,
)
from rna_lib_design import logger, structure_set, util, design
from rna_lib_design.structure_set import StructureSet, AddType
from rna_lib_design.structure import rna_structure, Structure

log = logger.setup_applevel_logger()


class Barcoder(object):
    def __init__(self, btype):
        self.btype: str = btype
        self.p5: Optional[StructureSet] = None
        self.p3: Optional[StructureSet] = None
        self.p5_buffer: Optional[StructureSet] = None
        self.p3_buffer: Optional[StructureSet] = None

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
        return get_best_designs_in_dataframe(df, sets, design_opts)


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
        df_result, f"{kwargs['name']}/results", kwargs["opool_name"]
    )
