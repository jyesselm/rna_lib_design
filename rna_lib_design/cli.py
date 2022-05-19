import os
import cloup
import numpy as np
from cloup import option_group, option
import pandas as pd
from typing import Dict, Tuple, Optional
from dataclasses import dataclass

import vienna
from seq_tools.data_frame import get_folded_structure, trim_3p, trim_5p

from rna_lib_design import logger, settings, util
from rna_lib_design.design import (
    DesignOptions,
    str_to_range,
    HelixRandomizer,
    get_best_designs_in_dataframe,
)
from rna_lib_design.structure import rna_structure, common_structures, Structure
from rna_lib_design.structure_set import (
    StructureSet,
    AddType,
    get_single_struct_set,
)
from rna_lib_design.barcode import (
    SingleBarcoder,
    DoubleBarcode,
    TripleBarcode,
    CustomSingleBarcoder,
    CustomDoubleBarcode,
    output_results,
)

log = logger.setup_applevel_logger()

common_sequence_options = option_group(
    "Common sequence options",
    "options for controling what is added to every sequence for priming and RT",
    option(
        "-p5seq",
        "--p5-sequence",
        help="p5 sequence. The sequence that will appear at the 5' of every "
        "designed sequence. [default=GGAAGAUCGAGUAGAUCAAA]",
    ),
    option(
        "-p5ss",
        "--p5-structure",
        default=None,
        help="the p5 secondary structure, if not supplied will fold p5 sequence "
        "with vienna fold. [default=None]",
    ),
    option(
        "-p5n",
        "--p5-name",
        default=None,
        help="supply a named p5 sequence/structure",
    ),
    option(
        "-p3seq",
        "--p3-sequence",
        help="p3 sequence. The sequence that will appear at the 3' of every "
        "designed sequence. [default=AAAGAAACAACAACAACAAC]",
    ),
    option(
        "-p3ss",
        "--p3-structure",
        default=None,
        help="the p3 secondary structure, if not supplied will fold p3 sequence "
        "with vienna fold. [default=None]",
    ),
    option(
        "-p3n",
        "--p3-name",
        default=None,
        help="supply a named p3 sequence/structure",
    ),
    option("--no-p5", is_flag=True, help="Do not add a common p5 sequence"),
    option("--no-p3", is_flag=True, help="Do not add a common p3 sequence"),
    option(
        "-p5bseq",
        "--p5b-sequence",
        help="p5 buffer sequence. A sequence between p5 and each construct. "
        "[default=None]",
    ),
    option(
        "-p5bss",
        "--p5b-structure",
        default=None,
        help="the p5 buffer secondary structure, if not supplied will fold p5 "
        "buffer sequence with vienna fold. [default=None]",
    ),
    option(
        "-p3bseq",
        "--p3b-sequence",
        help="p3 buffer sequence. A sequence between p3 and each construct. "
        "[default=None]",
    ),
    option(
        "-p3bss",
        "--p3b-structure",
        default=None,
        help="the p3 buffer secondary structure, if not supplied will fold p3 "
        "buffer sequence with vienna fold. [default=None]",
    ),
    option(
        "-lseq",
        "--loop-sequence",
        type=str,
        help="sequence of loop to have in the hairpin (default=CUUCGG)",
    ),
    option(
        "-lss",
        "--loop-structure",
        type=str,
        help="structure of loop to have in the hairpin (default=None)",
    ),
    option(
        "-ln",
        "--loop-name",
        type=str,
        help="specify existing loop sequence/structure from common sequences"
        "(default=None)",
    ),
)

csv_options = option_group(
    "CSV options",
    "options to apply to the supplied csv file",
    option("--trim_5p", default=0, type=int),
    option("--trim_3p", default=0, type=int),
)

design_options = option_group(
    "Design options",
    "options to control how generate each construct",
    option(
        "--max_ens_defect",
        default=2,
        type=int,
        help="the max increase of ensemble defect by adding additional sequences. "
        "[default=2]",
    ),
    option(
        "--max_design_attempts",
        default=100,
        type=int,
        help="how many times to try generating a design before giving up. "
        "[default=100]",
    ),
    option(
        "--max_design_solutions",
        default=10,
        type=int,
        help="how many solutions to generate before returning the best. "
        "[default=10]",
    ),
)

output_options = option_group(
    "Output options",
    "options to control output",
    option("-ll", "--log-level", default="info"),
    option(
        "-n",
        "--name",
        default="results",
        help="what to call the results folder",
    ),
    option("-on", "--opool-name", default="opool"),
    option("-m", "--max", default=1000000, type=str),
)


class CLIParser(object):
    @staticmethod
    def set_log_level(args):
        log.info(f"barcode type is {args['btype']}")
        if args["log_level"] == "debug":
            log.setLevel(logger.logging.DEBUG)
        elif args["log_level"] == "info":
            log.setLevel(logger.logging.INFO)

    @staticmethod
    def get_dataframe(args) -> pd.DataFrame:
        df = pd.read_csv(args["csv"])
        if "name" not in df:
            df["name"] = [f"seq_{x}" for x in range(len(df))]
        if "sequence" not in df:
            log.error(f"csv must contain a `sequence` column")
            exit()
        if "structure" not in df:
            get_folded_structure(df)
        if args["trim_5p"] > 0:
            trim_5p(df, args["trim_5p"])
        if args["trim_3p"] > 0:
            trim_3p(df, args["trim_3p"])
        return df

    @staticmethod
    def get_design_opts(args) -> DesignOptions:
        opts = DesignOptions(
            args["max_ens_defect"],
            args["max_design_attempts"],
            args["max_design_solutions"],
        )
        log.debug(f"design options are {opts}")
        return opts

    @staticmethod
    def get_p5(args) -> Optional[StructureSet]:
        if args["no_p5"]:
            return None
        p5_seq, p5_ss = args["p5_sequence"], args["p5_structure"]
        p5_name = args["p5_name"]
        p5_sequences = pd.read_csv(
            settings.RESOURCES_PATH + "/p5_sequences.csv"
        )
        if p5_seq is None and p5_name is None:
            row = p5_sequences.loc[0]
            p5 = rna_structure(row["sequence"], row["structure"])
            log.info(f"no p5 sequence supplied using: {p5.sequence}")
        elif p5_seq is not None and p5_name is not None:
            raise ValueError("cannot supply both p5_sequence and p5_name")
        elif p5_name in p5_sequences["name"].unique():
            row = p5_sequences[p5_sequences["name"] == p5_name].iloc[0]
            log.info(f"p5 supplied: {p5_name}")
            log.info(
                f"p5 sequence: {row['sequence']}, structure: {row['structure']}, "
                f"code: {row['code']}"
            )
            p5 = rna_structure(row["sequence"], row["structure"])
        elif p5_seq is not None and p5_ss is not None:
            log.info(f"p5 sequence supplied: {p5_seq}")
            log.info(f"p5 structure supplied: {p5_ss}")
            p5 = rna_structure(p5_seq, p5_ss)
        else:
            r = vienna.fold(p5_seq)
            log.info(f"p5 sequence supplied: {p5_seq}")
            log.info(
                f"p5 sequence has a folded structure of {r.dot_bracket} with "
                f"ens_defect of {r.ens_defect}"
            )
            p5 = rna_structure(p5_seq, r.dot_bracket)
        return get_single_struct_set(p5, AddType.LEFT)

    @staticmethod
    def get_p3(args) -> Optional[StructureSet]:
        if args["no_p3"]:
            return None
        p3_seq, p3_ss = args["p3_sequence"], args["p3_structure"]
        p3_name = args["p3_name"]
        p3_sequences = pd.read_csv(
            settings.RESOURCES_PATH + "/p3_sequences.csv"
        )
        if p3_seq is None and p3_name is None:
            row = p3_sequences.loc[0]
            p3 = rna_structure(row["sequence"], row["structure"])
            log.info(f"no p3 sequence supplied using: {p3.sequence}")
        elif p3_seq is not None and p3_name is not None:
            raise ValueError("cannot supply both p3_sequence and p3_name")
        elif p3_name in p3_sequences["name"].unique():
            row = p3_sequences[p3_sequences["name"] == p3_name].iloc[0]
            log.info(f"p3 supplied: {p3_name}")
            log.info(
                f"p3 sequence: {row['sequence']}, structure: {row['structure']}, "
                f"code: {row['code']}"
            )
            p3 = rna_structure(row["sequence"], row["structure"])
        elif p3_seq is not None and p3_ss is not None:
            log.info(f"p3 sequence supplied: {p3_seq}")
            log.info(f"p3 structure supplied: {p3_ss}")
            p5 = rna_structure(p3_seq, p3_ss)
        else:
            r = vienna.fold(p3_seq)
            log.info(f"p3 sequence supplied: {p3_seq}")
            log.info(
                f"p3 sequence has a folded structure of {r.dot_bracket} with "
                f"ens_defect of {r.ens_defect}"
            )
            p3 = rna_structure(p3_seq, r.dot_bracket)
        return get_single_struct_set(p3, AddType.RIGHT)

    @staticmethod
    def get_loop(args) -> Structure:
        loop_seq, loop_ss = args["loop_sequence"], args["loop_structure"]
        loop_name = args["loop_name"]
        common_structs = common_structures()
        if loop_seq is not None and loop_name is not None:
            log.error("cannot supply both loop_sequence and loop_name!")
            exit()
        elif loop_seq is None and loop_name is None:
            loop_struct = common_structs["uucg_loop"]
            log.info("no loop info was suplied will be using CUUCGG if needed")
        elif loop_name in common_structs:
            loop_struct = common_structs[loop_name]
            log.info(f"loop name supplied: {loop_name}")
            log.info(
                f"seq: {loop_struct.sequence} ss: {loop_struct.dot_bracket}"
            )
        elif loop_seq is not None and loop_ss is not None:
            loop_struct = rna_structure(loop_seq, loop_ss)
        else:
            r = vienna.fold(loop_seq)
            log.info(f"loop sequence supplied: {loop_seq}")
            log.info(
                f"loop sequence has a folded structure of {r.dot_bracket} with "
                f"ens_defect of {r.ens_defect}"
            )
            loop_struct = rna_structure(loop_seq, r.dot_bracket)
        return loop_struct

    @staticmethod
    def get_p5_buffer(args) -> Optional[StructureSet]:
        p5b_seq, p5b_ss = args["p5b_sequence"], args["p5b_structure"]
        if p5b_seq is None:
            return None
        log.info(f"p5 buffer sequence supplied {p5b_seq}")
        if p5b_ss is not None:
            log.info(f"p5 buffer sequence supplied {p5b_ss}")
            p5_buffer_struct = rna_structure(p5b_seq, p5b_ss)
        else:
            p5_buffer_struct = rna_structure(
                p5b_seq, vienna.fold(p5b_seq).dot_bracket
            )
        p5_buffer = get_single_struct_set(p5_buffer_struct, AddType.LEFT)
        return p5_buffer

    @staticmethod
    def get_p3_buffer(args) -> Optional[StructureSet]:
        p3b_seq, p3b_ss = args["p3b_sequence"], args["p3b_structure"]
        if p3b_seq is None:
            return None
        log.info(f"p3 buffer sequence supplied {p3b_seq}")
        if p3b_ss is not None:
            log.info(f"p3 buffer sequence supplied {p3b_ss}")
            p3_buffer_struct = rna_structure(p3b_seq, p3b_ss)
        else:
            p3_buffer_struct = rna_structure(
                p3b_seq, vienna.fold(p3b_seq).dot_bracket
            )
        p3_buffer = get_single_struct_set(p3_buffer_struct, AddType.RIGHT)
        return p3_buffer


@dataclass(frozen=True, order=True)
class CLIData:
    df: pd.DataFrame
    design_opts: DesignOptions
    p5: Optional[StructureSet]
    p3: Optional[StructureSet]
    loop: Optional[Structure] = None
    p5_buffer: Optional[StructureSet] = None
    p3_buffer: Optional[StructureSet] = None


def parse_cli(args: Dict):
    CLIParser.set_log_level(args)
    df = CLIParser.get_dataframe(args)
    opts = CLIParser.get_design_opts(args)
    p5 = CLIParser.get_p5(args)
    p3 = CLIParser.get_p3(args)
    loop = CLIParser.get_loop(args)
    p5_buffer = CLIParser.get_p5_buffer(args)
    p3_buffer = CLIParser.get_p3_buffer(args)
    return CLIData(df, opts, p5, p3, loop, p5_buffer, p3_buffer)


@cloup.group()
def cli():
    pass


@cli.command(help="Add a single barcode to an rna library")
@cloup.argument("csv", type=cloup.Path(exists=True))
@cloup.option_group(
    "main options",
    cloup.option(
        "-t",
        "--btype",
        type=cloup.Choice(["helix", "hairpin"]),
        default="helix",
        help="what type of barcoding? (default=`helix`)",
    ),
    cloup.option(
        "-l",
        "--length",
        default=6,
        type=int,
        help="how long should the barcode be? (default=6)",
    ),
    cloup.option(
        "--add_5p",
        is_flag=True,
        help="if barcode is a hairpin add it on the 5' of each sequence instead"
        " of the 3'",
    ),
)
@csv_options
@output_options
@common_sequence_options
@design_options
def barcode(**kwargs):
    os.makedirs(kwargs["name"], exist_ok=True)
    global log
    log = logger.setup_applevel_logger(
        file_name=kwargs["name"] + "/barcode.log"
    )
    cli_data = parse_cli(kwargs)
    add_type = AddType.RIGHT
    if kwargs["add_5p"]:
        add_type = AddType.LEFT
    btype = kwargs["btype"].lower()
    bcoder = SingleBarcoder(btype, kwargs["length"], add_type)
    bcoder.set_p5_and_p3(cli_data.p5, cli_data.p3)
    bcoder.set_loop(cli_data.loop)
    bcoder.set_p5_buffer(cli_data.p5_buffer)
    bcoder.set_p3_buffer(cli_data.p3_buffer)
    df_result = bcoder.barcode(cli_data.df, cli_data.design_opts)
    output_results(df_result, kwargs)


@cli.command(help="Add double barcodes to an rna library")
@cloup.argument("csv", type=cloup.Path(exists=True))
@cloup.option_group(
    "main options",
    cloup.option(
        "-t",
        "--btype",
        type=cloup.Choice(
            ["helix_hairpin", "hairpin_helix", "hairpin_hairpin"]
        ),
        default="helix",
        help="what type of barcoding? (default=`helix_hairpin`)",
    ),
    cloup.option(
        "-l",
        "--lengths",
        nargs=2,
        default=(6, 6),
        help="how long should the barcode be? (default=(6,6))",
    ),
)
@csv_options
@output_options
@common_sequence_options
@design_options
def barcode2(**kwargs):
    os.makedirs(kwargs["name"], exist_ok=True)
    global log
    log = logger.setup_applevel_logger(
        file_name=kwargs["name"] + "/barcode.log"
    )
    cli_data = parse_cli(kwargs)
    btype = kwargs["btype"].lower()
    bcoder = DoubleBarcode(btype, kwargs["lengths"], cli_data.loop)
    bcoder.set_p5_and_p3(cli_data.p5, cli_data.p3)
    bcoder.set_p5_buffer(cli_data.p5_buffer)
    bcoder.set_p3_buffer(cli_data.p3_buffer)
    df_result = bcoder.barcode(cli_data.df, cli_data.design_opts)
    output_results(df_result, kwargs)


@cli.command(help="Add triple barcodes to an rna library")
@cloup.argument("csv", type=cloup.Path(exists=True))
@cloup.option_group(
    "main options",
    cloup.option(
        "-t",
        "--btype",
        type=cloup.Choice(["hairpin_helix_hairpin"]),
        default="hairpin_helix_hairpin",
        help="what type of barcoding? (default=`hairpin_helix_hairpin`)",
    ),
    cloup.option(
        "-l",
        "--lengths",
        nargs=3,
        default=(6, 6, 6),
        help="how long should the barcode be? (default=(6,6,6))",
    ),
)
@csv_options
@output_options
@common_sequence_options
@design_options
def barcode3(**kwargs):
    os.makedirs(kwargs["name"], exist_ok=True)
    global log
    log = logger.setup_applevel_logger(
        file_name=kwargs["name"] + "/barcode.log"
    )
    cli_data = parse_cli(kwargs)
    btype = kwargs["btype"].lower()
    bcoder = TripleBarcode(btype, kwargs["lengths"], cli_data.loop)
    bcoder.set_p5_and_p3(cli_data.p5, cli_data.p3)
    bcoder.set_p5_buffer(cli_data.p5_buffer)
    bcoder.set_p3_buffer(cli_data.p3_buffer)
    df_result = bcoder.barcode(cli_data.df, cli_data.design_opts)
    output_results(df_result, kwargs)


@cli.command(help="Add a single barcode to an rna library")
@cloup.argument("csv", type=cloup.Path(exists=True))
@cloup.option_group(
    "main options",
    cloup.option(
        "-t",
        "--btype",
        type=cloup.Choice(["helix", "hairpin"]),
        default="helix",
        help="what type of barcoding? (default=`helix`)",
    ),
    cloup.option(
        "-bf", "--barcode-file", type=str, help="specify barcode file"
    ),
    cloup.option(
        "--add_5p",
        is_flag=True,
        help="if barcode is a hairpin add it on the 5' of each sequence instead"
        " of the 3'",
    ),
    cloup.option("-p", "--pooled", type=int, default=-1),
)
@csv_options
@output_options
@common_sequence_options
@design_options
def cbarcode(**kwargs):
    os.makedirs(kwargs["name"], exist_ok=True)
    global log
    log = logger.setup_applevel_logger(
        file_name=kwargs["name"] + "/barcode.log"
    )
    cli_data = parse_cli(kwargs)
    add_type = AddType.RIGHT
    if kwargs["add_5p"]:
        add_type = AddType.LEFT
    btype = kwargs["btype"].lower()
    bcoder = CustomSingleBarcoder(btype, kwargs["barcode_file"], add_type)
    bcoder.set_p5_and_p3(cli_data.p5, cli_data.p3)
    bcoder.set_loop(cli_data.loop)
    bcoder.set_p5_buffer(cli_data.p5_buffer)
    bcoder.set_p3_buffer(cli_data.p3_buffer)
    bcoder.multiprocess = kwargs["pooled"]
    df_result = bcoder.barcode(cli_data.df, cli_data.design_opts)
    output_results(df_result, kwargs)


@cli.command(help="Add two barcodes to an rna library")
@cloup.argument("csv", type=cloup.Path(exists=True))
@cloup.option_group(
    "main options",
    cloup.option(
        "-t",
        "--btype",
        type=cloup.Choice(
            ["helix_hairpin", "hairpin_helix", "hairpin_hairpin"]
        ),
        default="helix_hairpin",
        help="what type of barcoding? (default=`helix_hairpin`)",
    ),
    cloup.option(
        "-bf", "--barcode-files", type=str, help="specify barcode files"
    ),
    cloup.option(
        "--add_5p",
        is_flag=True,
        help="if barcode is a hairpin add it on the 5' of each sequence instead"
        " of the 3'",
    ),
    cloup.option("-p", "--pooled", type=int, default=-1),
)
@csv_options
@output_options
@common_sequence_options
@design_options
def cbarcode2(**kwargs):
    os.makedirs(kwargs["name"], exist_ok=True)
    global log
    log = logger.setup_applevel_logger(
        file_name=kwargs["name"] + "/barcode.log"
    )
    cli_data = parse_cli(kwargs)
    btype = kwargs["btype"].lower()
    bcoder = CustomDoubleBarcode(btype, kwargs["barcode_files"], cli_data.loop)
    bcoder.set_p5_and_p3(cli_data.p5, cli_data.p3)
    bcoder.set_p5_buffer(cli_data.p5_buffer)
    bcoder.set_p3_buffer(cli_data.p3_buffer)
    bcoder.multiprocess = kwargs["pooled"]
    df_result = bcoder.barcode(cli_data.df, cli_data.design_opts)
    output_results(df_result, kwargs)


@cli.command(
    help="randomize the helices of rnas will perseving their secondary structure"
)
@cloup.argument("csv", type=cloup.Path(exists=True))
@cloup.option("-e", "--exclude", default=None)
@cloup.option("-es", "--exclude-seq", default=None)
@csv_options
@output_options
def randhelix(**kwargs):
    df = CLIParser.get_dataframe(kwargs)
    if kwargs["log_level"] == "debug":
        log.setLevel(logger.logging.DEBUG)
    print(util.compute_edit_distance(df))
    exclude = []
    if kwargs["exclude"] is not None:
        exclude = str_to_range(kwargs["exclude"])
    else:
        exclude = []
    if kwargs["exclude_seq"] is not None:
        exclude_seqs = kwargs["exclude_seq"].split(",")
    else:
        exclude_seqs = None
    df["ens_defect"] = np.nan
    for i, row in df.iterrows():
        hr = HelixRandomizer()
        sol = hr.run(row["sequence"], row["structure"], exclude, exclude_seqs)
        df.at[i, ["sequence", "structure", "ens_defect"]] = [
            str(sol.design.sequence),
            str(sol.design.dot_bracket),
            sol.ens_defect,
        ]
        log.info(sol)
    print(util.compute_edit_distance(df))
    df.to_csv("output.csv", index=False)


@cli.command()
@csv_options
@output_options
@common_sequence_options
@design_options
@cloup.argument("csv", type=cloup.Path(exists=True))
def addcommon(**kwargs):
    kwargs["btype"] = "na"
    os.makedirs(kwargs["name"], exist_ok=True)
    global log
    log = logger.setup_applevel_logger(
        file_name=kwargs["name"] + "/barcode.log"
    )
    cli_data = parse_cli(kwargs)
    add_ons = [cli_data.p5_buffer, cli_data.p3_buffer, cli_data.p5, cli_data.p3]
    sets = []
    for ao in add_ons:
        if ao is not None:
            sets.append(ao)
    df_result = get_best_designs_in_dataframe(
        cli_data.df, sets, cli_data.design_opts
    )
    output_results(df_result, kwargs)


if __name__ == "__main__":
    cli()
