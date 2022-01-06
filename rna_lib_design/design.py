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


################################################################################
# cli functions                                                                #
################################################################################


@click.group()
def cli():
    pass


@cli.command()
@click.option(
    "-l", "--lengths", default=(6, 6), type=(int, int), help="lengths"
)
def barcode2(type, lengths, csv, p5_common, p3_common, output, **kwargs):
    __setup_logging(type, kwargs["log_level"].lower())
    df = __df_setup(csv, kwargs)
    opts = DesignOptions(
        kwargs["max_ens_defect"],
        kwargs["max_design_attempts"],
        kwargs["max_design_solutions"],
    )
    loop_struct = __parse_cli_loop(kwargs["loop"])
    sets = []
    type = type.lower()
    if type == "helix_hairpin":
        h_set = structure_set.get_optimal_helix_set(lengths[0], len(df) * 1.1)
        hp_set = structure_set.get_optimal_hairpin_set(
            lengths[1], loop_struct, len(df) * 1.1, structure_set.AddType.RIGHT
        )
        sets.extend([h_set, hp_set])
    elif type == "hairpin_helix":
        h_set = structure_set.get_optimal_helix_set(lengths[1], len(df) * 1.1)
        hp_set = structure_set.get_optimal_hairpin_set(
            lengths[0], loop_struct, len(df) * 1.1, structure_set.AddType.LEFT
        )
        sets.extend([h_set, hp_set])
    elif type == "hairpin_hairpin":
        hp_set_1 = structure_set.get_optimal_hairpin_set(
            lengths[0], loop_struct, len(df) * 1.1, structure_set.AddType.LEFT
        )
        hp_set_2 = structure_set.get_optimal_hairpin_set(
            lengths[1], loop_struct, len(df) * 1.1, structure_set.AddType.RIGHT
        )
        sets.extend([hp_set_1, hp_set_2])
    else:
        log.error(f"unknown type {type}")
        exit()

    sets.extend(__parse_cli_p5_and_p3(p5_common, p3_common))
    df_result = __generate_designs_for_dataframe(df, sets, opts)
    org_total = len(df_result)
    df_result = df_result.dropna()
    diff = org_total - len(df_result)
    log.info(f"{diff} constructs were discarded as they were below cutoffs")
    df_result.to_csv(output, index=False)


@cli.command()
@click.option(
    "-l", "--lengths", default=(6, 6, 6), type=(int, int, int), help="lengths"
)
def barcode3(type, lengths, csv, p5_common, p3_common, output, **kwargs):
    __setup_logging(type, kwargs["log_level"].lower())
    df = __df_setup(csv, kwargs)
    opts = DesignOptions(
        kwargs["max_ens_defect"],
        kwargs["max_design_attempts"],
        kwargs["max_design_solutions"],
    )
    loop_struct = __parse_cli_loop(kwargs["loop"])
    sets = []
    type = type.lower()
    if type == "hairpin_helix_hairpin":
        h_set = structure_set.get_optimal_helix_set(lengths[0], len(df) * 1.1)
        hp_set_1 = structure_set.get_optimal_hairpin_set(
            lengths[1], loop_struct, len(df) * 1.1, structure_set.AddType.RIGHT
        )
        hp_set_2 = structure_set.get_optimal_hairpin_set(
            lengths[1], loop_struct, len(df) * 1.1, structure_set.AddType.LEFT
        )
        sets.extend([h_set, hp_set_1, hp_set_2])
    sets.extend(__parse_cli_p5_and_p3(p5_common, p3_common))
    df_result = __generate_designs_for_dataframe(df, sets, opts)
    org_total = len(df_result)
    df_result = df_result.dropna()
    diff = org_total - len(df_result)
    log.info(f"{diff} constructs were discarded as they were below cutoffs")
    df_result.to_csv(output, index=False)


################################################################################
# assemble functions                                                           #
################################################################################


def __get_set_from_sequence(name: str, n: Dict):
    if "add" in n:
        add_type = structure_set.str_to_add_type(n["add"])
    else:
        add_type = None
    sequence = n["sequence"]
    if "structure" in n:
        struct = n["structure"]
    else:
        struct = vienna.fold(sequence).dot_bracket
    log.info(
        f"node: {name} is intepreted as sequence node with seq: "
        f"{sequence} and ss: {struct}"
    )
    struct = structure.rna_structure(sequence, struct)
    if add_type is None:
        add_type = structure_set.AddType.LEFT
        if sequence.find("&") > -1:
            add_type = structure_set.AddType.HELIX
    return structure_set.get_single_struct_set(struct, add_type)


def __get_set_from_library(name: str, n: Dict, num: int):
    log.info(
        f"node: {name} is intepreted as library node of type: "
        f"{n['library']}"
    )
    if "add" in n:
        add_type = structure_set.str_to_add_type(n["add"])
    else:
        add_type = None
    lib_name = n["library"].lower()
    length = 6
    size = 10
    if num is not None:
        size = num
    if "length" in n:
        length = n["length"]
    if "size" in n:
        size = n["size"]
    if lib_name == "helix":
        h_set = structure_set.get_optimal_helix_set(length, size)
        return h_set
    elif lib_name == "hairpin":
        loop_struct = None
        if "loop_name" in n:
            loop_struct = structure.get_common_struct(n["loop_name"])
        elif "loop_sequence" in n and "loop_structure" in n:
            loop_struct = structure.rna_structure(
                n["loop_sequence"], n["loop_structure"]
            )
        else:
            loop_struct = structure.get_common_struct("uucg_loop")
        hp_set = structure_set.get_optimal_hairpin_set(
            length, loop_struct, size, add_type
        )
        return hp_set
    else:
        log.error(f"{lib_name} not implemented yet!")
        exit()


@cli.command()
@click.argument("yml")
@click.option("-n", "--num")
def assemble(yml, num):
    f = open(yml)
    nodes = yaml.load(f, Loader=yaml.FullLoader)
    sets = []
    start_struct = structure.rna_structure("", "")
    mode = "num"
    i = 0
    for name, n in nodes.items():
        add_type = structure_set.AddType.LEFT
        if "add" in n:
            add_type = structure_set.str_to_add_type(n["add"])
        if "sequence" in n:
            sets.append(__get_set_from_sequence(name, n))
        elif "library" in n:
            sets.append(__get_set_from_library(name, n, num))
        elif "named_struct" in n:
            sets.append(
                structure_set.get_common_seq_structure_set(
                    n["named_struct"], add_type
                )
            )
        elif "5_and_3p" in n:
            p5_name, p3_name = n["5_and_3p"]
            sets.extend(
                [
                    structure_set.get_common_seq_structure_set(
                        p5_name, structure_set.AddType.LEFT
                    ),
                    structure_set.get_common_seq_structure_set(
                        p3_name, structure_set.AddType.RIGHT
                    ),
                ]
            )

            pass
        elif "file" in n:
            pass

    opts = DesignOptions()
    sol = get_best_design(sets, start_struct, opts)
    print(sol.design)
