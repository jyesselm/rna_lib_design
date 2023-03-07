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
    design,
    write_output_dir,
)
from seq_tools import calc_edit_distance
from rna_lib_design.logger import get_logger, setup_applevel_logger
from rna_lib_design.parameters import parse_parameters_from_file
from rna_lib_design.settings import get_resources_path

log = get_logger("CLI")


# TODO validate build str does it include everythingp ?


def get_preset_parameters(btype: str, barcode_name: str):
    full_name = f"{barcode_name}_{btype}.yml"
    if not (get_resources_path() / "presets" / full_name).exists():
        presets = get_resources_path() / "presets").glob(f"{barcode_name}*.yml")
        log.error(presets)
        raise ValueError(f"Invalid barcode type: {btype}")
    return get_resources_path() / "presets" / full_name


def merge_dict(original_dict, new_dict):
    for key, value in new_dict.items():
        if isinstance(value, dict):
            original_dict[key] = merge_dict(original_dict.get(key, {}), value)
        else:
            # safe guard since this might really miss stuff up
            if key == "m_type":
                if original_dict[key] != value:
                    raise ValueError(
                        f"Cannot change m_type from {original_dict[key]} to {value}"
                    )
            original_dict[key] = value
    return original_dict


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
            "average length"
        )



def log_failed_sequences(results):
    df_results = results.df_results
    failures = results.failures
    table = []
    total = 0
    for key, value in failures.items():
        if value > 0:
            table.append([key, value])
        total += value
    if len(table) > 0:
        log.info(
            "summary of sequences discarded\n"
            + tabulate(table, headers=["key", "value"], tablefmt="psql")
        )
        log.info(f"total sequences discarded: {total}")
        log.info(f"total remaining sequences: {len(df_results)}")
    else:
        log.info("no sequences discarded")


def run_method(method_name, csv, btype, param_file, output, args):
    is_valid_method(method_name)
    os.makedirs(output, exist_ok=True)
    setup_log_and_log_inputs(csv, btype, param_file, output, args["debug"])
    validate_initial_library(csv)
    schema_file = get_resources_path() / "schemas" / f"{method_name}.json"
    # setup parameters
    params = {}
    if btype is None:
        preset_file = None
    else:
        preset_file = get_preset_parameters(btype.lower(), method_name)
        params = parse_parameters_from_file(preset_file, schema_file)
    if param_file is not None:
        if len(params) == 0:
            params = parse_parameters_from_file(param_file, schema_file)
        else:
            user_params = parse_parameters_from_file(param_file, schema_file)
            merge_dict(params, user_params)
    log.info(f"Writing parameters to {output}/params.yml")
    yaml.dump(params, open(f"{output}/params.yml", "w"))
    # default to standard preset if no preset or param file is supplied
    if preset_file is None and param_file is None:
        log.info("No preset or param file supplied, using standard preset")
        preset_file = get_resources_path() / "presets" / f"{method_name}_standard.yml"
        params = parse_parameters_from_file(preset_file, schema_file)
    if len(params) == 0:
        raise ValueError("No parameters supplied")
    log.info(f"Using parameters:\n{json.dumps(params, indent=4)}")
    design_opts = DesignOpts(**params["design_opts"])
    df_seqs = pd.read_csv(csv)
    results = design(
        args["num_processes"],
        df_seqs,
        params["build_str"],
        params["segments"],
        design_opts,
    )
    log_failed_sequences(results)
    df_results = results.df_results
    write_output_dir(df_results, Path(output))
    if not args["skip_edit_dist"]:
        edit_dist = calc_edit_distance(df_results)
        log.info(f"the edit distance of lib is: {edit_dist}")
    else:
        log.info("skipping edit distance calculation")


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
    run_method("add_common", csv, btype, param_file, output, args)


@cli.command()
@cloup.argument("csv", type=cloup.Path(exists=True))
@main_options()
def barcode(csv, btype, param_file, output, **args):
    """
    adds a single barcode
    """
    run_method("single_barcode", csv, btype, param_file, output, args)


@cli.command()
@cloup.argument("csv", type=cloup.Path(exists=True))
@main_options()
def barcode2(csv, btype, param_file, output, **args):
    """
    adds two barcodes
    """
    run_method("double_barcode", csv, btype, param_file, output, args)


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
