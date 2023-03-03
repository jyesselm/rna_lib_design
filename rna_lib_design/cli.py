import json
import os
import pandas as pd
import cloup

from rna_lib_design.design import (
    Designer,
    DesignOpts,
    design_multiprocess,
    write_results_to_file,
)
from rna_lib_design.logger import get_logger, setup_applevel_logger
from rna_lib_design.parameters import parse_parameters_from_file
from rna_lib_design.settings import get_resources_path
from rna_lib_design.util import compute_edit_distance

log = get_logger("CLI")


# TODO add cli tool "list" to list all the available p5, p3, and common sequences
# TODO add cli tool to list which barcode sets are available
# TODO have checks for the library like diversity and size difference between sequences
# TODO validate build str does it include everything?
# TODO list and determine the P5 and P3 code if used
# TODO need to update primers


def get_barcode_preset_parameters(btype: str, barcode_name: str):
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
    log.info(f"Using parameters:\n{json.dumps(params, indent=4)}")
    design_opts = DesignOpts(**params["design_opts"])
    df_seqs = pd.read_csv(csv)
    if num_processes == 1:
        designer = Designer()
        designer.setup(design_opts)
        df_results = designer.design(df_seqs, params["build_str"], params["segments"])
    else:
        df_results = design_multiprocess(
            num_processes,
            df_seqs,
            params["build_str"],
            params["segments"],
            design_opts,
        )
    return df_results


@cloup.group()
def cli():
    pass


@cli.command()
@cloup.argument("csv", type=cloup.Path(exists=True))
@cloup.option("-t", "--btype", type=str, default=None)
@cloup.option("--param-file", type=cloup.Path(exists=True), default=None)
@cloup.option("-o", "--output", default="results")
def barcode(csv, btype, param_file, output):
    setup_log_and_log_inputs(csv, btype, param_file)
    schema_file = get_resources_path() / "schemas" / "single_barcode.json"
    if btype is None:
        preset_file = None
    else:
        preset_file = get_barcode_preset_parameters(btype.lower(), "single_barcode")
    df_results = run_design(csv, schema_file, preset_file, param_file, 8)
    os.makedirs(output, exist_ok=True)
    write_results_to_file(df_results)
    edit_dist = compute_edit_distance(df_results)
    log.info(f"the edit distance of lib is: {edit_dist}")


@cli.command()
@cloup.argument("csv", type=cloup.Path(exists=True))
def edit_distance(csv):
    setup_log_and_log_inputs(csv, None, None)
    df = pd.read_csv(csv)
    log.info("edit distance:" + str(compute_edit_distance(df)))


if __name__ == "__main__":
    cli()
