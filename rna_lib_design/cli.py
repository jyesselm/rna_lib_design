import json
import os
import pandas as pd
import cloup
from cloup import option_group, option
from pathlib import Path

from rna_lib_design.design import (
    DesignOpts,
    design,
    write_results_to_file,
)
from seq_tools import calc_edit_distance
from rna_lib_design.logger import get_logger, setup_applevel_logger
from rna_lib_design.parameters import parse_parameters_from_file
from rna_lib_design.settings import get_resources_path

log = get_logger("CLI")


# TODO add cli tool "list" to list all the available p5, p3, and common sequences
# TODO add cli tool to list which barcode sets are available
# TODO have checks for the library like diversity and size difference between sequences
# TODO validate build str does it include everythingp ?
# TODO list and determine the P5 and P3 code if used
# TODO need to update primers
# TODO add standard option if nothing is supplied for btype


def get_preset_parameters(btype: str, barcode_name: str):
    full_name = f"{barcode_name}_{btype}.yml"
    if not (get_resources_path() / "presets" / full_name).exists():
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


def setup_log_and_log_inputs(csv, preset, param_file):
    setup_applevel_logger()
    log.info(f"Using csv: {csv}")
    if preset is not None:
        log.info(f"Using preset: {preset}")
    if param_file is not None:
        log.info(f"Using parameter file: {param_file}")


def run_design(csv, schema_file, preset_file, param_file, num_processes):
    params = {}
    if preset_file is not None:
        params = parse_parameters_from_file(preset_file, schema_file)
    if param_file is not None:
        if len(params) == 0:
            params = parse_parameters_from_file(param_file, schema_file)
        else:
            user_params = parse_parameters_from_file(param_file, schema_file)
            merge_dict(params, user_params)
    if len(params) == 0:
        raise ValueError("No parameters supplied")
    log.info(f"Using parameters:\n{json.dumps(params, indent=4)}")
    design_opts = DesignOpts(**params["design_opts"])
    df_seqs = pd.read_csv(csv)
    df_results = design(
        num_processes,
        df_seqs,
        params["build_str"],
        params["segments"],
        design_opts,
    )
    return df_results


def is_valid_method(method_name):
    methods = ["single_barcode", "double_barcode", "triple_barcode"]
    if method_name not in methods:
        raise ValueError(f"Invalid method: {method_name}")


def run_method(method_name, csv, btype, param_file, output, args):
    setup_log_and_log_inputs(csv, btype, param_file)
    schema_file = get_resources_path() / "schemas" / f"{method_name}.json"
    if btype is None:
        preset_file = None
    else:
        preset_file = get_preset_parameters(btype.lower(), method_name)
    if preset_file is None and param_file is None:
        preset_file = get_resources_path() / "presets" / f"{method_name}_standard.yml"
    results = run_design(
        csv, schema_file, preset_file, param_file, args["num_processes"]
    )
    df_results = results.df_results
    write_results_to_file(df_results, Path("results"))
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
            "--skip-edit-dist", is_flag=True, help="skip the edit distance calculation"
        ),
    )


@cloup.group()
def cli():
    pass


@cli.command()
@cloup.argument("csv", type=cloup.Path(exists=True))
@main_options()
def barcode(csv, btype, param_file, output, **args):
    run_method("single_barcode", csv, btype, param_file, output, args)


@cli.command()
@cloup.argument("csv", type=cloup.Path(exists=True))
@main_options()
def add_common(csv, btype, param_file, output, **args):
    run_method("add_common", csv, btype, param_file, output, args)


@cli.command()
@cloup.argument("csv", type=cloup.Path(exists=True))
def edit_distance(csv):
    setup_log_and_log_inputs(csv, None, None)
    df = pd.read_csv(csv)
    log.info("edit distance:" + str(calc_edit_distance(df)))


if __name__ == "__main__":
    cli()
