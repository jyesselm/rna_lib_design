import json
import os
import shutil
import pandas as pd
import cloup
import yaml
from tabulate import tabulate
from cloup import option_group, option
from pathlib import Path

from rna_lib_design.design import (
    DesignOpts,
    design_and_save_output,
    write_output_dir,
    log_failed_design_sequences,
)
from seq_tools import calc_edit_distance, trim
from rna_lib_design.logger import get_logger, setup_applevel_logger
from rna_lib_design.parameters import (
    parse_parameters_from_file,
    combine_params,
    get_preset_parameters,
)
from rna_lib_design.settings import get_resources_path

log = get_logger("CLI")


# TODO add triming of sequences p5 and p3
# TODO validate build str does it include everythingp ?


def setup_log_and_log_inputs(csv, preset, param_file, output, debug):
    setup_applevel_logger(is_debug=debug, file_name=f"{output}/log.txt")
    log.info(f"Using csv: {csv}")
    log.info(f"Using output dir: {output}")
    log.info(f"Copying {csv} to {output}/input.csv")
    try:
        shutil.copy(csv, f"{output}/input.csv")
    except:
        pass
    if debug:
        log.info("Debug mode enabled")
    if preset is not None:
        log.info(f"Using preset: {preset}")
    if param_file is not None:
        log.info(f"Using parameter file: {param_file}")


def is_valid_method(method_name):
    methods = ["add_common", "single_barcode", "double_barcode", "triple_barcode"]
    if method_name not in methods:
        raise ValueError(f"Invalid method: {method_name}")


# TODO check for edit distance of library?
def validate_initial_library(csv):
    df = pd.read_csv(csv)
    log.info(f"csv has {len(df)} sequences")
    min_len = df["sequence"].str.len().min()
    max_len = df["sequence"].str.len().max()
    avg_len = df["sequence"].str.len().mean()
    # max_len - min_len is greater than 10% of avg length error
    if max_len - min_len > avg_len * 0.1:
        raise ValueError(
            "The library size difference is too large must be under 10% of the"
            "average length, can turn this off with --skip-length-check"
        )


def get_method_params(method_name, btype, param_file, args):
    is_valid_method(method_name)
    schema_file = get_resources_path() / "schemas" / f"{method_name}.json"
    log.info(f"Using schema file: {schema_file}")
    params = {}
    if btype is not None:
        params = get_preset_parameters(btype.lower(), method_name)
    if param_file is not None:
        if len(params) == 0:
            params = parse_parameters_from_file(param_file, schema_file)
        else:
            user_params = parse_parameters_from_file(param_file, schema_file)
            combine_params(params, user_params)
    else:
        log.info("No preset or param file supplied, using standard preset")
        preset_file = get_resources_path() / "presets" / f"{method_name}_standard.yml"
        log.info(f"using file {preset_file}")
        params = parse_parameters_from_file(preset_file, schema_file)
    if len(params) == 0:
        raise ValueError("No parameters supplied")
    # hacky bring back the script that updates this automatically
    params["debug"] = args["debug"]
    params["num_of_processes"] = args["num_processes"]
    params["preprocess"]["trim_p5"] = args["trim_p5"]
    params["preprocess"]["trim_p3"] = args["trim_p3"]
    params["preprocess"]["skip_length_check"] = args["skip_length_check"]
    # params["preprocess"]["skip_edit_distance_check"] = args["skip_edit_distance_check"]
    params["postprocess"]["skip_edit_distance"] = args["skip_edit_dist"]
    return params


def setup_method(method_name, csv, btype, param_file, output, args):
    is_valid_method(method_name)
    os.makedirs(output, exist_ok=True)
    setup_log_and_log_inputs(csv, btype, param_file, output, args["debug"])
    params = get_method_params(method_name, btype, param_file, args)
    df_seqs = pd.read_csv(csv)
    if params["preprocess"]["trim_p5"] != 0 or params["preprocess"]["trim_p3"] != 0:
        log.info(
            f"trimming sequences by {params['preprocess']['trim_p5']} at 5' and"
            f"{params['preprocess']['trim_p3']} at 3'"
        )
        df_seqs = trim(
            df_seqs, params["preprocess"]["trim_p5"], params["preprocess"]["trim_p3"]
        )
    return params, df_seqs


# cli commands ########################################################################


def main_options():
    return option_group(
        "Main options",
        "These are the main options for generating a library",
        option(
            "-t",
            "--btype",
            type=str,
            default=None,
            help="what type of barcode to use see full list in resources/presets",
        ),
        option(
            "--param-file",
            type=cloup.Path(exists=True),
            default=None,
            help=(
                "supply a new param file to override specific present or to manually"
                "determine each option"
            ),
        ),
        option("-o", "--output", default="results", help="the path to save results to"),
        option(
            "-p",
            "--num-processes",
            type=int,
            default=1,
            help="number of processes to run simultaneously",
        ),
        option(
            "--debug", is_flag=True, help="turn on debug logging for the application"
        ),
        option(
            "--skip-edit-dist", is_flag=True, help="skip the edit distance calculation"
        ),
        option(
            "--skip-length-check",
            is_flag=True,
            help="skip the check for the length variation of the library",
        ),
        option(
            "--trim-p5",
            type=int,
            default=0,
            help="trim sequence at 5' end by this length",
        ),
        option(
            "--trim-p3",
            type=int,
            default=0,
            help="trim sequence at 3' end by this length",
        ),
    )


@cloup.group()
def cli():
    """
    This is a simple tool to generate barcodes for a library of RNAs to be ordered
    as an DNA oligo library.
    """
    pass


@cli.command()
@cloup.argument("csv", type=cloup.Path(exists=True))
@main_options()
def add_common(csv, btype, param_file, output, **args):
    """
    add common p5/p3 sequences
    """
    params, df_seqs = setup_method("add_common", csv, btype, param_file, output, args)
    design_and_save_output(df_seqs, output, params)


@cli.command()
@cloup.argument("csv", type=cloup.Path(exists=True))
@main_options()
def barcode(csv, btype, param_file, output, **args):
    """
    adds a single barcode
    """
    params, df_seqs = setup_method(
        "single_barcode", csv, btype, param_file, output, args
    )
    design_and_save_output(df_seqs, output, params)


@cli.command()
@cloup.argument("csv", type=cloup.Path(exists=True))
@main_options()
def barcode2(csv, btype, param_file, output, **args):
    """
    adds two barcodes
    """
    params, df_seqs = setup_method(
        "double_barcode", csv, btype, param_file, output, args
    )
    design_and_save_output(df_seqs, output, params)


@cli.command()
@cloup.argument("csv", type=cloup.Path(exists=True))
def edit_distance(csv):
    """
    compute edit distance of library
    """
    setup_applevel_logger()
    log.info(f"Using csv: {csv}")
    df = pd.read_csv(csv)
    log.info("edit distance:" + str(calc_edit_distance(df)))


@cli.command()
@cloup.argument("name", type=str)
def list(name):
    """
    lists resources available
    """
    setup_applevel_logger()
    seqs = {"p5": "p5_sequences", "p3": "p3_sequences", "loops": "loops"}
    barcodes = {"helix_barcodes": "helices", "sstrand_barcodes": "sstrands"}
    if name in seqs:
        dir_path = get_resources_path() / f"named_seqs/rna/{seqs[name]}.csv"
        df = pd.read_csv(dir_path)
        log.info(
            f"{name} sequences available\n"
            + tabulate(df, headers="keys", tablefmt="psql", showindex=False)
        )
    elif name in barcodes:
        dir_path = get_resources_path() / f"barcodes/{barcodes[name]}.csv"
        df = pd.read_csv(dir_path)
        log.info(
            f"{name} sets available\n"
            + tabulate(df, headers="keys", tablefmt="psql", showindex=False)
        )
    else:
        raise ValueError(f"unknown name: {name}")


if __name__ == "__main__":
    cli()
