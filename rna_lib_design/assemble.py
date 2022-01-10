import click
import yaml
from typing import Dict

import vienna
from rna_lib_design import structure, structure_set, logger

log = logger.setup_applevel_logger()


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


@click.command()
@click.argument("yml")
@click.option("-n", "--num")
def cli(yml, num):
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


if __name__ == "__main__":
    cli()
