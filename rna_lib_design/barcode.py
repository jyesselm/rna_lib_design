import click
import pandas as pd
import numpy as np
import shutil
import os
import logging
from typing import Dict, Tuple
from pathlib import Path

import vienna
import seq_tools.data_frame as stdf
from rna_lib_design.design import (
    get_best_designs_in_dataframe,
    write_results_to_file,
    DesignOptions,
)
from rna_lib_design import (
    logger,
    structure_set,
    structure,
    defaults,
    util,
    design,
)

log = logger.setup_applevel_logger(file_name="barcode.log")


def single_barcode(
    df: pd.DataFrame,
    btype: str,
    length: int,
    p5: structure_set.StructureSet,
    p3: structure_set.StructureSet,
    design_opts: DesignOptions,
    add_type: structure_set.AddType = structure_set.AddType.RIGHT,
    loop: str = None,
    p5_buffer = None,
    p3_buffer = None
) -> pd.DataFrame:
    if btype == "helix":
        v_set = structure_set.get_optimal_helix_set(length, len(df) * 1.2)
    elif btype == "sstrand":
        v_set = structure_set.get_optimal_sstrand_set(
            length, len(df) * 1.2, add_type
        )
    elif btype == "hairpin":
        loop_struct = defaults.get_loop_from_str(loop)
        v_set = structure_set.get_optimal_hairpin_set(
            length, loop_struct, len(df) * 1.2, add_type
        )
    else:
        log.error(f"{btype} is not a valid type")
        exit()
    sets = [v_set]
    setup_buffer_sequences(p5_buffer, p3_buffer, sets)
    sets.extend([p5, p3])
    df_result = get_best_designs_in_dataframe(df, sets, design_opts)
    return df_result


def double_barcode(
    df: pd.DataFrame,
    btype: str,
    lengths: Tuple[int, int],
    p5: structure_set.StructureSet,
    p3: structure_set.StructureSet,
    design_opts: DesignOptions,
    loop: str = None,
    p5_buffer: str = None,
    p3_buffer: str = None,
):
    sets = []
    loop_struct = defaults.get_loop_from_str(loop)
    if btype == "helix_hairpin":
        h_set = structure_set.get_optimal_helix_set(lengths[0], len(df) * 1.2)
        hp_set = structure_set.get_optimal_hairpin_set(
            lengths[1], loop_struct, len(df) * 1.2, structure_set.AddType.RIGHT
        )
        sets.extend([h_set, hp_set])
    elif btype == "hairpin_helix":
        h_set = structure_set.get_optimal_helix_set(lengths[1], len(df) * 1.2)
        hp_set = structure_set.get_optimal_hairpin_set(
            lengths[0], loop_struct, len(df) * 1.2, structure_set.AddType.LEFT
        )
        sets.extend([h_set, hp_set])
    elif btype == "hairpin_hairpin":
        hp_set_1 = structure_set.get_optimal_hairpin_set(
            lengths[0], loop_struct, len(df) * 1.2, structure_set.AddType.LEFT
        )
        hp_set_2 = structure_set.get_optimal_hairpin_set(
            lengths[1], loop_struct, len(df) * 1.2, structure_set.AddType.RIGHT
        )
        sets.extend([hp_set_1, hp_set_2])
    else:
        log.error(f"unknown type {btype}")
        exit()
    setup_buffer_sequences(p5_buffer, p3_buffer, sets)
    sets.extend([p5, p3])
    df_result = get_best_designs_in_dataframe(df, sets, design_opts)
    return df_result


def triple_barcode(
    df: pd.DataFrame,
    btype: str,
    lengths: Tuple[int, int, int],
    p5: structure_set.StructureSet,
    p3: structure_set.StructureSet,
    design_opts: DesignOptions,
    loop: str = None,
):
    sets = []
    loop_struct = defaults.get_loop_from_str(loop)
    if btype == "hairpin_helix_hairpin":
        hp_set_1 = structure_set.get_optimal_hairpin_set(
            lengths[0], loop_struct, len(df) * 1.2, structure_set.AddType.LEFT
        )
        h_set = structure_set.get_optimal_helix_set(lengths[1], len(df) * 1.2)
        hp_set_2 = structure_set.get_optimal_hairpin_set(
            lengths[2], loop_struct, len(df) * 1.2, structure_set.AddType.RIGHT
        )
        sets.extend([h_set, hp_set_1, hp_set_2])
    sets.extend([p5, p3])
    df_result = get_best_designs_in_dataframe(df, sets, design_opts)
    return df_result


################################################################################
# cli functions                                                                #
################################################################################


def setup_buffer_sequences(p5_buffer, p3_buffer, sets):
    if p5_buffer is not None:
        log.info(f"p5 buffer sequence supplied {p5_buffer}")
        p5_buffer_struct = structure.rna_structure(
            p5_buffer, vienna.fold(p5_buffer).dot_bracket
        )
        sets.append(
            structure_set.get_single_struct_set(
                p5_buffer_struct, structure_set.AddType.LEFT
            )
        )
    if p3_buffer is not None:
        log.info(f"p3 buffer sequence supplied {p3_buffer}")
        p3_buffer_struct = structure.rna_structure(
            p3_buffer, vienna.fold(p3_buffer).dot_bracket
        )
        sets.append(
            structure_set.get_single_struct_set(
                p3_buffer_struct, structure_set.AddType.RIGHT
            )
        )


def setup_logging(type, log_level):
    log.info(f"barcode type is {type}")
    if log_level == "debug":
        log.setLevel(logging.DEBUG)


def setup_dataframe_from_cli(args: Dict) -> pd.DataFrame:
    df = pd.read_csv(args["csv"])
    if "name" not in df:
        df["name"] = [f"seq_{x}" for x in range(len(df))]
    if "sequence" not in df:
        log.error(f"csv must contain a `sequence` column")
        exit()
    if "structure" not in df:
        stdf.get_folded_structure(df)
    if args["trim_5p"] > 0:
        stdf.trim_5p(df, args["trim_5p"])
    if args["trim_3p"] > 0:
        stdf.trim_3p(df, args["trim_3p"])
    return df


def setup_design_opts_from_cli(args: Dict) -> DesignOptions:
    opts = DesignOptions(
        args["max_ens_defect"],
        args["max_design_attempts"],
        args["max_design_solutions"],
    )
    log.debug(f"design options are {opts}")
    return opts


def setup_from_cli(kwargs: Dict):
    setup_logging(kwargs["btype"], kwargs["log_level"])
    opts = setup_design_opts_from_cli(kwargs)
    df = setup_dataframe_from_cli(kwargs)
    p5 = defaults.get_p5_from_str(kwargs["p5_common"])
    p3 = defaults.get_p3_from_str(kwargs["p3_common"])
    return [df, p5, p3, opts]


def final_results(df_result, kwargs):
    org_total = len(df_result)
    df_result = df_result.dropna()
    diff = org_total - len(df_result)
    log.info(f"{diff} constructs were discarded as they were below cutoffs")
    return df_result


def common_options(function):
    function = click.argument("csv", type=click.Path(exists=True))(function)
    function = click.option("-p5", "--p5-common", help="p5 sequence")(function)
    function = click.option("-p3", "--p3-common", help="p3 sequence")(function)
    function = click.option("-o", "--output", default="out.csv")(function)
    function = click.option("-ll", "--log-level", default="info")(function)
    function = click.option("--trim_5p", default=0, type=int)(function)
    function = click.option("--trim_3p", default=0, type=int)(function)
    function = click.option("-med", "--max_ens_defect", default=2, type=int)(
        function
    )
    function = click.option(
        "-mda", "--max_design_attempts", default=100, type=int
    )(function)
    function = click.option(
        "-mds", "--max_design_solutions", default=10, type=int
    )(function)
    function = click.option("--loop")(function)
    function = click.option("-n", "--name", help="what to call the results")(
        function
    )
    function = click.option("-on", "--opool-name", default="opool")(function)
    function = click.option(
        "-n",
        "--name",
        default="barcode_results",
        help="what to call the results",
    )(function)
    function = click.option("-p3b", "--p3-buffer")(function)
    function = click.option("-p5b", "--p5-buffer")(function)
    return function


@click.group()
def cli():
    pass


@cli.command()
@click.option("-t", "--btype", default="helix", help="type")
@click.option("-l", "--length", default=6, type=int, help="length")
@click.option("--add_5p", is_flag=True)
@common_options
def barcode(**kwargs):
    df, p5, p3, opts = setup_from_cli(kwargs)
    add_type = structure_set.AddType.RIGHT
    if kwargs["add_5p"]:
        add_type = structure_set.AddType.LEFT
    btype = kwargs["btype"].lower()
    df_result = single_barcode(
        df,
        btype,
        kwargs["length"],
        p5,
        p3,
        opts,
        loop=kwargs["loop"],
        add_type=add_type,
        p5_buffer=kwargs["p5_buffer"],
        p3_buffer=kwargs["p3_buffer"]
    )
    final_results(df_result, kwargs)
    util.compute_edit_distance(df_result)


@cli.command()
@click.option("-t", "--btype", default="helix_hairpin", help="type")
@click.option(
    "-l", "--lengths", default=(6, 6), type=(int, int), help="lengths"
)
@common_options
def barcode2(**kwargs):
    df, p5, p3, opts = setup_from_cli(kwargs)
    log.info(f"lengths: {kwargs['lengths']}")
    btype = kwargs["btype"].lower()
    df_result = double_barcode(
        df,
        btype,
        kwargs["lengths"],
        p5,
        p3,
        opts,
        loop=kwargs["loop"],
        p5_buffer=kwargs["p5_buffer"],
        p3_buffer=kwargs["p3_buffer"],
    )
    os.makedirs(kwargs["name"], exist_ok=True)
    design.write_results_to_file(
        df, f"{kwargs['name']}/results", kwargs["name"]
    )
    if os.path.isfile(f"{kwargs['name']}/barcode.log"):
        os.remove(f"{kwargs['name']}/barcode.log")
    util.compute_edit_distance(df_result)
    shutil.move("barcode.log", kwargs["name"])


@cli.command()
@click.option("-t", "--btype", default="hairpin_helix_hairpin", help="type")
@click.option(
    "-l", "--lengths", default=(6, 6, 6), type=(int, int, int), help="lengths"
)
@common_options
def barcode3(**kwargs):
    df, p5, p3, opts = setup_from_cli(kwargs)
    btype = kwargs["btype"].lower()
    df_result = triple_barcode(
        df, btype, kwargs["lengths"], p5, p3, opts, loop=kwargs["loop"]
    )
    df_result = final_results(df_result, kwargs)
    util.compute_edit_distance(df_result)


if __name__ == "__main__":
    cli()
