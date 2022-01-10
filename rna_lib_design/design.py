import logging

import click
import pandas as pd
import numpy as np
from dataclasses import dataclass
from typing import List, Dict
import yaml

import vienna
import seq_tools.data_frame as stdf
from rna_lib_design import structure_set, logger, structure

log = logger.setup_applevel_logger()


################################################################################
# classes and functions                                                        #
################################################################################


@dataclass(frozen=True, order=True)
class DesignResults:
    ens_defect: float
    design: structure.Structure


@dataclass(frozen=True, order=True)
class DesignOptions:
    max_ens_defect: float = 2.0
    max_design_attempts: int = 100
    max_design_solutions: int = 10


def get_best_design(
    sd_dicts: List[structure_set.StructureSet],
    struct: structure.Structure,
    opts: DesignOptions,
    org_ens_defect: float = 999.0,
) -> DesignResults:
    count = 0
    best_score = 1000
    best = None
    best_uses = []
    if org_ens_defect > 100:
        if len(struct) > 5:
            org_ens_defect = vienna.fold(struct).ensemble_diversity
        else:
            org_ens_defect = 1
    while 1:
        count += 1
        final_struct = structure_set.apply(sd_dicts, struct)
        vr = vienna.fold(final_struct.sequence)
        if count > opts.max_design_attempts:
            return DesignResults(999, structure.rna_structure("A", "."))
        # score = final_struct.dot_bracket_difference(structure.DotBracket(vr.dot_bracket))
        # print(final_struct.dot_bracket, vr.ensemble_diversity, score)
        if final_struct.dot_bracket != vr.dot_bracket:
            continue
        if opts.max_ens_defect < vr.ensemble_diversity - org_ens_defect:
            continue
        if vr.ensemble_diversity < best_score:
            best_score = vr.ensemble_diversity
            best = final_struct
            best_uses = [x.get_current_pos() for x in sd_dicts]
        if count > opts.max_design_solutions:
            break
    for i, bu in enumerate(best_uses):
        sd_dicts[i].set_used(bu)
    return DesignResults(best_score, best)


def get_best_designs_in_dataframe(
    df: pd.DataFrame,
    sets: List[structure_set.StructureSet],
    opts: DesignOptions,
) -> pd.DataFrame:
    df["org_sequence"] = df["sequence"]
    df["org_structure"] = df["structure"]
    if "ens_defect" in df:
        df["org_ens_defect"] = df["ens_defect"]
    else:
        df["ens_defect"] = np.nan
    for i, row in df.iterrows():
        struct = structure.rna_structure(row["sequence"], row["structure"])
        if "org_ens_defect" in df:
            sol = get_best_design(sets, struct, opts, row["org_ens_defect"])
        else:
            sol = get_best_design(sets, struct, opts)
        if sol.ens_defect > 100:
            log.error(
                f"unable to find a design for {row['name']} below the cutoff "
                f"{opts.max_ens_defect}"
            )
            df.at[i, ["sequence", "structure", "ens_defect"]] = [
                np.nan,
                np.nan,
                np.nan,
            ]
        else:
            log.debug(f"solution found: name={row['name']}, sol={sol},")
            df.at[i, ["sequence", "structure", "ens_defect"]] = [
                str(sol.design.sequence),
                str(sol.design.dot_bracket),
                sol.ens_defect,
            ]
    return df


def write_results_to_file(
    df: pd.DataFrame, fname="results", opool_name="opool"
) -> None:
    df.to_csv(f"{fname}-all.csv", index=False)
    df_sub = df[["name", "sequence", "structure", "ens_defect"]]
    df_sub.to_csv(f"{fname}-rna.csv", index=False)

    df_sub = df[["name", "sequence"]].copy()
    df_sub["sequence"] = [
        "TTCTAATACGACTCACTATA" + seq for seq in df_sub["sequence"]
    ]
    stdf.convert_to_dna(df_sub)
    df_sub.to_csv(f"{fname}-dna.csv", index=False)
    df_sub = df_sub.rename(
        columns={"name": "Pool name", "sequence": "Sequence"}
    )
    df_sub["Pool name"] = opool_name
    df_sub.to_excel(f"{fname}-opool.xlsx", index=False)
