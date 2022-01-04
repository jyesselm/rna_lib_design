import logging

import click
import pandas as pd
import numpy as np
from dataclasses import dataclass
from typing import List, Dict
from tabulate import tabulate
import textwrap
import yaml

import vienna
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


def write_results_to_csv(constructs, fname="final"):
    f_rna = open(fname + "_rna.csv", "w")
    f_rna.write("name,sequence,structure,ens_defect\n")
    f_dna = open(fname + "_dna.csv", "w")
    f_dna.write("name,sequence\n")
    for c in constructs:
        f_rna.write(f"{c[0]},{c[1].sequence},{c[1].dot_bracket},{c[2]}\n")
        f_dna.write(f"{c[0]},TTCTAATACGACTCACTATA{c[1].sequence.to_dna()}\n")
    f_rna.close()
    f_dna.close()


################################################################################
# helper functions                                                             #
################################################################################


def __print_design_solution(name: str, sol: DesignResults) -> None:
    data = {
        "construct": name,
        "sequence": textwrap.fill(str(sol.design.sequence), 80),
        "structure": textwrap.fill(str(sol.design.dot_bracket), 80),
        "ens_defect": sol.ens_defect,
    }
    df = pd.DataFrame(list(data.items()), columns="name value".split())
    log.debug("\n" + tabulate(df, headers="keys", tablefmt="psql", showindex=False))


def __generate_designs_for_dataframe(
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
            __print_design_solution(row["name"], sol)
            df.at[i, ["sequence", "structure", "ens_defect"]] = [
                str(sol.design.sequence),
                str(sol.design.dot_bracket),
                sol.ens_defect,
            ]
    return df


def __parse_cli_p5_and_p3(
    p5_common: str, p3_common: str
) -> List[structure_set.StructureSet]:
    common_structs = structure.common_structures()
    if p5_common is None:
        p5 = common_structs["ref_hairpin_5prime"]
        log.info(f"no p5 sequence supplied using: {p5.sequence}")
    else:
        r = vienna.fold(p5_common)
        log.info(f"p5 sequence supplied: {p5_common}")
        log.info(
            f"p5 sequence has a folded structure of {r.dot_bracket} with "
            f"ens_defect of {r.ensemble_diversity}"
        )
        p5 = structure.rna_structure(p5_common, r.dot_bracket)

    if p3_common is None:
        p3 = common_structs["rt_tail"]
        log.info(f"no p3 sequence supplied using: {p3.sequence}")
    else:
        r = vienna.fold(p3_common)
        log.info(f"p5 sequence supplied: {p3_common}")
        log.info(
            f"p5 sequence has a folded structure of {r.dot_bracket} with "
            f"ens_defect of {r.ensemble_diversity}"
        )
        p3 = structure.rna_structure(p3_common, r.dot_bracket)
    p5_set = structure_set.get_single_struct_set(p5, structure_set.AddType.LEFT)
    p3_set = structure_set.get_single_struct_set(p3, structure_set.AddType.RIGHT)
    return [p5_set, p3_set]


def __parse_cli_loop(loop: str) -> structure.Structure:
    if loop is None:
        loop_struct = structure.get_common_struct("uucg_loop")
    else:
        r = vienna.fold(loop)
        log.info(f"loop sequence supplied: {loop}")
        log.info(
            f"loop sequence has a folded structure of {r.dot_bracket} with "
            f"ens_defect of {r.ensemble_diversity}"
        )
        loop_struct = structure.rna_structure(loop, r.dot_bracket)
    return loop_struct


def __add_folded_structure(df):
    structures = []
    ens_defects = []
    for i, row in df.iterrows():
        vr = vienna.fold(row["sequence"])
        structures.append(vr.dot_bracket)
        ens_defects.append(vr.ensemble_diversity)
    df["structure"] = structures
    df["ens_defect"] = ens_defects


def __trim_sequence(df, trim_5p, trim_3p):
    if trim_5p == 0 and trim_3p == 0:
        return
    for i, row in df.iterrows():
        sequence = row["sequence"]
        structure = row["structure"]
        if trim_5p > 0:
            sequence = sequence[trim_5p:]
            structure = structure[trim_5p:]
        if trim_3p > 0:
            sequence = sequence[:-trim_3p]
            structure = structure[:-trim_3p]
        df.at[i, ["sequence", "structure"]] = [sequence, structure]


def __df_setup(csv, opts) -> pd.DataFrame:
    df_input = pd.read_csv(csv)
    if "sequence" not in df_input:
        log.error(f"csv must contain a `sequence` column")
        exit()
    if "structure" not in df_input:
        __add_folded_structure(df_input)
    __trim_sequence(df_input, opts["trim_5p"], opts["trim_3p"])
    return df_input


def __setup_logging(type, log_level):
    log.info(f"barcode type is {type}")
    if log_level == "debug":
        log.setLevel(logging.DEBUG)


################################################################################
# cli functions                                                                #
################################################################################


def common_options(function):
    function = click.argument("csv", type=click.Path(exists=True))(function)
    function = click.option("-t", "--type", default="helix", help="type")(function)
    function = click.option("-p5", "--p5-common", help="p5 sequence")(function)
    function = click.option("-p3", "--p3-common", help="p3 sequence")(function)
    function = click.option("-o", "--output", default="out.csv")(function)
    function = click.option("-ll", "--log-level", default="info")(function)
    function = click.option("--trim_5p", default=0, type=int)(function)
    function = click.option("--trim_3p", default=0, type=int)(function)
    function = click.option("-med", "--max_ens_defect", default=2, type=int)(function)
    function = click.option("-mda", "--max_design_attempts", default=100, type=int)(
        function
    )
    function = click.option("-mds", "--max_design_solutions", default=10, type=int)(
        function
    )
    function = click.option("--loop")(function)
    return function


@click.group()
def cli():
    pass


@cli.command()
@click.option("-l", "--length", default=6, type=int, help="length")
@click.option(
    "--add_3p", "add_type", flag_value=structure_set.AddType.RIGHT, default=True
)
@click.option("--add_5p", "add_type", flag_value=structure_set.AddType.LEFT)
@common_options
def barcode(type, length, csv, p5_common, p3_common, output, **kwargs):
    __setup_logging(type, kwargs["log_level"].lower())
    df = __df_setup(csv, kwargs)
    opts = DesignOptions(
        kwargs["max_ens_defect"],
        kwargs["max_design_attempts"],
        kwargs["max_design_solutions"],
    )
    if type.lower() == "helix":
        v_set = structure_set.get_optimal_helix_set(length, len(df) * 1.1)
    elif type.lower() == "sstrand":
        v_set = structure_set.get_optimal_sstrand_set(
            length, len(df) * 1.1, kwargs["add_type"]
        )
    elif type.lower() == "hairpin":
        loop_struct = __parse_cli_loop(kwargs["loop"])
        v_set = structure_set.get_optimal_hairpin_set(
            length, loop_struct, len(df) * 1.1, kwargs["add_type"]
        )
    else:
        log.error(f"{type} is not a valid type")
        exit()
    p5_set, p3_set = __parse_cli_p5_and_p3(p5_common, p3_common)
    sets = [v_set, p5_set, p3_set]
    df_result = __generate_designs_for_dataframe(df, sets, opts)
    org_total = len(df_result)
    df_result = df_result.dropna()
    diff = org_total - len(df_result)
    log.info(f"{diff} constructs were discarded as they were below cutoffs")
    df_result.to_csv(output, index=False)


@cli.command()
@click.option("-l", "--lengths", default=(6, 6), type=(int, int), help="lengths")
@common_options
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
@common_options
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
    log.info(f"node: {name} is intepreted as library node of type: " f"{n['library']}")
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
                structure_set.get_common_seq_structure_set(n["named_struct"], add_type)
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


if __name__ == "__main__":
    cli()
