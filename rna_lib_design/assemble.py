import click
import sys
import os
import yaml
import pandas as pd
from typing import Dict
import logging
import shutil

import vienna
from rna_lib_design import structure, structure_set, logger, design, defaults
from rna_lib_design.design import DesignOptions

log = logger.setup_applevel_logger(file_name="assemble.log")


def setup_design_opts_from_cli(args: Dict) -> DesignOptions:
    opts = DesignOptions(
        args["max_ens_defect"],
        args["max_design_attempts"],
        args["max_design_solutions"],
    )
    log.debug(f"design options are {opts}")
    return opts


def setup_from_cli(kwargs: Dict):
    if kwargs["log_level"] == "debug":
        log.setLevel(logging.DEBUG)
    opts = setup_design_opts_from_cli(kwargs)
    df = pd.read_csv(kwargs["csv"])
    p5 = defaults.get_p5_from_str(kwargs["p5_common"])
    p3 = defaults.get_p3_from_str(kwargs["p3_common"])
    return [df, p5, p3, opts]


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


@click.group()
def cli():
    pass


@cli.command()
@click.argument("yml")
@click.option("-n", "--num", default=1, type=int)
def assemble(yml, num):
    f = open(yml)
    nodes = yaml.load(f, Loader=yaml.FullLoader)
    sets = []
    start_struct = structure.rna_structure("", "")
    mode = "num"
    i = 0
    for name, n in nodes.items():
        if "file" in n:
            mode = "file"
            print(n)
    exit()

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

    # opts = DesignOptions()
    # sol = get_best_design(sets, start_struct, opts)
    # print(sol.design)


def common_options(function):
    function = click.argument("csv", type=click.Path(exists=True))(function)
    function = click.option("-p5", "--p5-common", help="p5 sequence")(function)
    function = click.option("-p3", "--p3-common", help="p3 sequence")(function)
    function = click.option("-ll", "--log-level", default="info")(function)
    function = click.option("-med", "--max_ens_defect", default=2, type=int)(
        function
    )
    function = click.option(
        "-mda", "--max_design_attempts", default=100, type=int
    )(function)
    function = click.option(
        "-mds", "--max_design_solutions", default=10, type=int
    )(function)
    function = click.option(
        "-n",
        "--name",
        default="assemble_results",
        help="what to call the results",
    )(function)
    return function


def hairpin_options(function):
    function = click.option("-urh", "--use-ref-hairpin", is_flag=True)(function)
    function = click.option("--ref-hairpin-length", default=5, type=int)(
        function
    )
    function = click.option(
        "--ref-hairpin-orientation", default="LEFT", type=str
    )(function)
    function = click.option("--ref-hairpin-loop")(function)
    return function


@cli.command()
@click.option("-hpl", "--hairpin-length", default=6, type=int)
@click.option("-hl", "--helix-length", default=6, type=int)
@click.option("-p3b", "--p3-buffer")
@click.option("--add_5p")
@common_options
@hairpin_options
def twoways(**kwargs):
    log.info(" ".join(sys.argv))
    start_struct = structure.rna_structure("", "")
    df, p5, p3, opts = setup_from_cli(kwargs)
    if "name" not in df:
        df["name"] = [f"seq_{x}" for x in range(len(df))]
    buffer_set = None
    os.makedirs(kwargs["name"], exist_ok=True)
    shutil.copy(kwargs["csv"], kwargs["name"] + "/input.csv")
    if kwargs["p3_buffer"] is not None:
        buff_seq = kwargs["p3_buffer"]
        buff_ss = vienna.folded_structure(buff_seq)
        buffer_set = structure_set.get_single_struct_set(
            structure.rna_structure(buff_seq, buff_ss),
            structure_set.AddType.RIGHT,
        )
    loop_struct = structure.get_common_struct("uucg_loop")
    hp_set = structure_set.get_optimal_hairpin_set(
        kwargs["hairpin_length"],
        loop_struct,
        len(df) * 2,
        structure_set.AddType.RIGHT,
    )
    hp_set.remove_dots()
    hp_set.set_buffer(None)
    h_set = structure_set.get_optimal_helix_set(
        kwargs["helix_length"], len(df) * 2
    )
    extras = [p5]
    if buffer_set is not None:
        extras.append(buffer_set)
    if kwargs["use_ref_hairpin"]:
        if kwargs["ref_hairpin_loop"]:
            hp_loop = structure.get_common_struct(kwargs["ref_hairpin_loop"])
        else:
            hp_loop = defaults.get_loop_from_str(None)
        add_type = structure_set.AddType.RIGHT
        if kwargs["ref_hairpin_orientation"] == "LEFT":
            add_type = structure_set.AddType.LEFT
        ref_hairpin_set = structure_set.get_optimal_hairpin_set(
            kwargs["ref_hairpin_length"], hp_loop, len(df) * 2, add_type
        )
        if kwargs["ref_hairpin_orientation"] == "LEFT":
            extras.insert(0, ref_hairpin_set)
        else:
            extras.append(ref_hairpin_set)
    extras.append(p3)
    df.insert(0, "ens_defect", "")
    df.insert(0, "structure", "")
    df.insert(0, "sequence", "")
    for i, row in df.iterrows():
        df_sub = df[df.index == i]
        df_set = structure_set.SingleStructureSet(
            df_sub, structure_set.AddType.HELIX
        )
        sets = [hp_set, df_set, h_set]
        sets.extend(extras)
        r = design.get_best_design(sets, start_struct, opts)
        if len(r.design) == 1:
            log.info(f"could not find solution to {row['name']}")
        else:
            log.info(r)
        df.at[i, ["sequence", "structure", "ens_defect"]] = [
            str(r.design.sequence),
            str(r.design.dot_bracket),
            r.ens_defect,
        ]
    design.write_results_to_file(
        df, f"{kwargs['name']}/results", kwargs["name"]
    )
    if os.path.isfile(f"{kwargs['name']}/assemble.log"):
        os.remove(f"{kwargs['name']}/assemble.log")
    shutil.move("assemble.log", kwargs["name"])
    # print(r)
    exit()


if __name__ == "__main__":
    cli()
