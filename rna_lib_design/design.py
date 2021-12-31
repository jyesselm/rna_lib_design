import logging

import click
import pandas as pd
import numpy as np
from dataclasses import dataclass
from typing import List, Tuple
from tabulate import tabulate
import textwrap

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
    max_ens_defect: float = 10.0
    max_design_attempts: int = 100
    max_design_solutions: int = 10


def get_best_design(
    sd_dicts: List[structure_set.StructureSet],
    struct: structure.Structure,
    opts: DesignOptions,
) -> DesignResults:
    count = 0
    best_score = 1000
    best = None
    best_uses = []
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
        if opts.max_ens_defect < vr.ensemble_diversity:
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
    log.debug(
        "\n" + tabulate(df, headers="keys", tablefmt="psql", showindex=False)
    )


def __update_dataframe_with_designs(
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
    p3_set = structure_set.get_single_struct_set(
        p3, structure_set.AddType.RIGHT
    )
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
    for i, row in df.iterrows():
        vr = vienna.fold(row["sequence"])
        structures.append(vr.dot_bracket)
    df["structure"] = structures


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


################################################################################
# cli functions                                                                #
################################################################################


def helix_hairpin_barcodes(
    df: pd.DataFrame,
    lengths: Tuple[int, int],
    p5: structure.Structure,
    p3: structure.Structure,
    loop: str,
    opts: DesignOptions,
) -> pd.DataFrame:
    loop_struct = __parse_cli_loop(loop)
    hp_set = structure_set.get_optimal_hairpin_set(
        lengths[1], loop_struct, len(df) * 1.1, structure_set.AddType.RIGHT
    )
    h_set = structure_set.get_optimal_helix_set(lengths[0], len(df) * 1.1)
    p5_set = structure_set.get_single_struct_set(p5, structure_set.AddType.LEFT)
    p3_set = structure_set.get_single_struct_set(
        p3, structure_set.AddType.RIGHT
    )
    sets = [h_set, hp_set, p5_set, p3_set]
    return __update_dataframe_with_designs(df, sets, opts)


def common_options(function):
    function = click.argument("csv", type=click.Path(exists=True))(function)
    function = click.option("-t", "--type", default="helix", help="type")(
        function
    )
    function = click.option("-p5", "--p5-common", help="p5 sequence")(function)
    function = click.option("-p3", "--p3-common", help="p3 sequence")(function)
    function = click.option("-o", "--output", default="out.csv")(function)
    function = click.option("-ll", "--log-level", default="info")(function)
    function = click.option("--trim_5p", default=0, type=int)(function)
    function = click.option("--trim_3p", default=0, type=int)(function)
    function = click.option("-med", "--max_ens_defect", default=10, type=int)(
        function
    )
    function = click.option(
        "-mda", "--max_design_attempts", default=100, type=int
    )(function)
    function = click.option(
        "-mds", "--max_design_solutions", default=10, type=int
    )(function)
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
@click.option("--loop")
@common_options
def barcode(type, length, csv, p5_common, p3_common, output, **kwargs):
    log.info(f"barcode type is {type}")
    log_level = kwargs["log_level"].lower()
    if log_level == "debug":
        log.setLevel(logging.DEBUG)
    df = __df_setup(csv, kwargs)
    opts = DesignOptions(
        kwargs["max_ens_defect"],
        kwargs["max_design_attempts"],
        kwargs["max_design_solutions"],
    )
    df_result = pd.DataFrame()
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
    df_result = __update_dataframe_with_designs(df, sets, opts)
    org_total = len(df_result)
    df_result = df_result.dropna()
    diff = org_total - len(df_result)
    log.info(f"{diff} constructs were discarded as they were below cutoffs")
    df_result.to_csv(output, index=False)


@cli.command()
@click.option(
    "-l", "--lengths", default=(6, 6), type=(int, int), help="lengths"
)
@click.option("--loop")
@common_options
def barcode2(type, lengths, csv, p5_common, p3_common, output, **kwargs):
    log.info(f"barcode type is {type}")
    log_level = kwargs["log_level"].lower()
    if log_level == "debug":
        log.setLevel(logging.DEBUG)
    df_input = __df_setup(csv, kwargs)
    p5, p3 = __parse_cli_p5_and_p3(p5_common, p3_common)
    opts = DesignOptions(
        kwargs["max_ens_defect"],
        kwargs["max_design_attempts"],
        kwargs["max_design_solutions"],
    )
    df_result = pd.DataFrame()
    loop_struct = __parse_cli_loop(kwargs["loop"])
    if type.lower() == "helix_hairpin":
        df_result = helix_hairpin_barcodes(
            df_input, lengths, p5, p3, kwargs["loop"], opts
        )
    org_total = len(df_result)
    df_result = df_result.dropna()
    diff = org_total - len(df_result)
    log.info(f"{diff} constructs were discarded as they were below cutoffs")
    df_result.to_csv(output, index=False)


if __name__ == "__main__":
    cli()
