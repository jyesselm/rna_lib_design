import click
import pandas as pd
import logging
from typing import Dict
from pathlib import Path

import seq_tools.data_frame as stdf
from rna_lib_design.design import (
    get_best_designs_in_dataframe,
    write_results_to_file,
    DesignOptions,
)
from rna_lib_design import logger, structure_set, structure, defaults

log = logger.setup_applevel_logger()


def single_barcode(
    df: pd.DataFrame,
    btype: str,
    length: int,
    p5: structure_set.StructureSet,
    p3: structure_set.StructureSet,
    design_opts: DesignOptions,
    add_type: structure_set.AddType = structure_set.AddType.RIGHT,
    loop: str = None,
) -> pd.DataFrame:
    if btype == "helix":
        v_set = structure_set.get_optimal_helix_set(length, len(df) * 1.1)
    elif btype == "sstrand":
        v_set = structure_set.get_optimal_sstrand_set(
            length, len(df) * 1.1, add_type
        )
    elif btype == "hairpin":
        loop_struct = defaults.get_loop_from_str(loop)
        v_set = structure_set.get_optimal_hairpin_set(
            length, loop_struct, len(df) * 1.1, add_type
        )
    else:
        log.error(f"{type} is not a valid type")
        exit()
    sets = [v_set, p5, p3]
    df_result = get_best_designs_in_dataframe(df, sets, design_opts)
    return df_result


################################################################################
# cli functions                                                                #
################################################################################


def setup_logging(type, log_level):
    log.info(f"barcode type is {type}")
    if log_level == "debug":
        log.setLevel(logging.DEBUG)


def setup_dataframe_from_cli(args: Dict) -> pd.DataFrame:
    df = pd.read_csv(args["csv"])
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
    fname = Path(kwargs["csv"]).stem
    write_results_to_file(df_result, fname, kwargs["opool_name"])


def common_options(function):
    function = click.argument("csv", type=click.Path(exists=True))(function)
    function = click.option("-t", "--btype", default="helix", help="type")(
        function
    )
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
    return function


@click.group()
def cli():
    pass


@cli.command()
@click.option("-l", "--length", default=6, type=int, help="length")
@click.option("--add_5p", is_flag=True)
@common_options
def barcode(**kwargs):
    df, p5, p3, opts = setup_from_cli(kwargs)
    add_type = structure_set.AddType.RIGHT
    if kwargs["add_5p"]:
        add_type = structure_set.AddType.LEFT
    btype = kwargs["btype"].lower()
    if btype != "hairpin":
        df_result = single_barcode(
            df, btype, kwargs["length"], p5, p3, opts, add_type=add_type
        )
    else:
        df_result = single_barcode(
            df,
            btype,
            kwargs["length"],
            p5,
            p3,
            opts,
            loop=kwargs["loop"],
            add_type=add_type,
        )
    final_results(df_result, kwargs)


if __name__ == "__main__":
    cli()
