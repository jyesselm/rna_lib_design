import multiprocessing
from collections import Counter
import os
from tabulate import tabulate
import yaml
import json

import pandas as pd
import numpy as np
from pathlib import Path

from dataclasses import dataclass
from vienna import fold
from vienna.vienna import FoldResults

from seq_tools import (
    fold as fold_seqs_in_df,
    to_dna,
    to_dna_template,
    to_fasta,
    calc_edit_distance,
)
from rna_lib_design.structure_set import (
    SequenceStructure,
    SequenceStructureSetParser,
)
from rna_lib_design.logger import get_logger
from rna_lib_design.util import get_seq_fwd_primer

log = get_logger("DESIGN")


def parse_build_str(seq_str):
    """
    Parses a string of the form
    "P5-HPBARCODE-HBARCODE6A-SOI-HBARCODE6B-AC-P3" to return a dictionary
    containing the order and direction in which the other segments appear
    in reference to SOI.

    :param seq_str: The sequence string to parse.
    :type seq_str: str
    :return: A dictionary containing the order and direction in
    which the other segments appear in reference to SOI.
    :rtype: dict[str, int]
    """
    segments = seq_str.split("-")
    soi_index = segments.index("SOI")
    if soi_index == -1:
        raise ValueError("no SOI in sequence string")
    order = {}
    for i, segment in enumerate(segments):
        if i == soi_index:
            continue
        distance = i - soi_index
        key = segment.rstrip("AB")
        if segment.endswith("A"):
            if key in order:
                order[key] = min(distance, order[key])
            else:
                order[key] = distance
        elif segment.endswith("B"):
            if key in order:
                order[key] = max(distance, order[key])
            else:
                order[key] = distance
        else:
            if key in order:
                order[key] = min(distance, order[key])
            else:
                order[key] = distance
    return order


@dataclass(frozen=True, order=True)
class DesignOpts:
    increase_ens_defect: float = 2.0
    max_ens_defect: float = 5.0
    max_attempts: int = 10
    max_solutions: int = 10
    score_method: str = "increase"
    allowed_ss_mismatch: int = 2
    allowed_ss_mismatch_barcodes: int = 2


@dataclass(frozen=True, order=True)
class DesignerResults:
    df_results: pd.DataFrame
    failures: dict


class Designer:
    def __init__(self):
        self.symbols = [
            "1",
            "2",
            "3",
            "4",
            "5",
            "6",
            "7",
            "8",
            "9",
            "0",
            "*",
            "#",
            "$",
            "%",
            "^",
            "*",
            "B",
            "D",
            "F",
        ]
        self.opts = DesignOpts()
        self.failures = {
            "high_ens_defect": 0,
            "ss_mismatches": 0,
            "ss_mismatches_barcodes": 0,
        }

    def setup(self, opts: DesignOpts):
        self.opts = opts

    def design(self, df_sequences, build_str, params):
        df_results = self.__setup_dataframe(df_sequences)
        build_up = self.__get_build_up(len(df_sequences), build_str, params)
        count = -1
        for i, row in df_results.iterrows():
            count += 1
            if count % 100 == 0 and count > 0:
                log.info(f"processed {count} sequences")
                print(f"processed {count} sequences")
            try:
                soi_seq_struct = SequenceStructure(
                    row["org_sequence"], row["org_structure"]
                )
            except:
                log.error(
                    f"failed process {row['name'] }{row['org_sequence']} - {row['org_structure']} skipping"
                )
                continue
            seq_struct, iterating_sets = self.__setup_seq_struct(
                soi_seq_struct, build_up
            )
            df_results.at[i, "design_sequence"] = seq_struct.sequence
            df_results.at[i, "design_structure"] = seq_struct.structure
            results = self.__get_designed_seq_struct(
                seq_struct, row["name"], row["org_ens_defect"], iterating_sets
            )
            # no design found
            if results[0] == "":
                continue
            df_results.at[i, "sequence"] = results[0]
            df_results.at[i, "structure"] = results[1]
            df_results.at[i, "ens_defect"] = results[2]
            df_results.at[i, "mfe"] = results[3]
        df_results = df_results[df_results["sequence"] != ""]
        return DesignerResults(df_results, self.failures)

    def __setup_dataframe(self, df):
        df = df.copy()
        if "sequence" not in df.columns:
            raise ValueError("no sequence column in dataframe")
        if "structure" not in df.columns:
            log.info("no 'structure' column folding it now")
            if "ens_defect" in df.columns:
                log.warn(
                    "ens_defect column found but not structure weird behavior will happen"
                )
            df = fold_seqs_in_df(df)
        if "ens_defect" not in df.columns:
            log.info("no 'ens_defect' column - adding one")
            if "structure" in df.columns:
                log.info("structure column will be overwritten with folded structure")
            df = fold_seqs_in_df(df)
        df.rename(
            columns={
                "sequence": "org_sequence",
                "structure": "org_structure",
                "ens_defect": "org_ens_defect",
            },
            inplace=True,
        )
        df["sequence"] = ""
        df["structure"] = ""
        df["design_sequence"] = ""
        df["design_structure"] = ""
        df["ens_defect"] = -999
        df["mfe"] = 999
        col_order = ["name", "sequence", "structure", "ens_defect", "mfe"]
        df = df.reindex(col_order + list(df.columns.difference(col_order)), axis=1)
        return df

    def __get_designed_seq_struct(
        self, seq_struct, name, org_ens_defect, iterating_sets
    ):
        best = []
        best_seq_struct = SequenceStructure("", "")
        best_r = FoldResults("", 999, 999, [])
        num_solutions = 0
        no_solution = True
        fails = []
        for i in range(0, self.opts.max_attempts):
            sequence = seq_struct.sequence
            structure = seq_struct.structure
            used = []
            for replace_str, cur_set in iterating_sets:
                sol = cur_set.get_random()
                used.append(sol)
                strands = sol.split_strands()
                for rs, strand in zip(replace_str, strands):
                    sequence = sequence.replace(rs, strand.sequence, 1)
                    structure = structure.replace(rs, strand.structure, 1)
            r = fold(sequence)
            result = self.__score_design(structure, seq_struct.structure, r)
            if result != "SUCCESS":
                fails.append(result)
                continue
            no_solution = False
            if r.ens_defect < best_r.ens_defect:
                best_r = r
                best_seq_struct = SequenceStructure(sequence, r.dot_bracket)
                best = used
            num_solutions += 1
            if num_solutions >= self.opts.max_solutions:
                break
        if no_solution:
            count = Counter(fails)
            self.failures[count.most_common(1)[0][0]] += 1
            log.debug("no design found for sequence: " + seq_struct.sequence)
        if self.opts.score_method == "increase":
            diff = best_r.ens_defect - org_ens_defect
            if diff > self.opts.increase_ens_defect:
                self.failures["high_ens_defect"] += 1
                log.debug(
                    f"design ens_defect increase too large: {str(diff)} for seq {name}"
                )
                return ["", "", -999, 999]
        elif self.opts.score_method == "max":
            if best_r.ens_defect > self.opts.max_ens_defect:
                self.failures["high_ens_defect"] += 1
                log.debug(f"design ens_defect too large: " + str(best_r.ens_defect))
                return ["", "", -999, 999]
        else:
            raise ValueError("unknown score method: " + self.opts.score_method)

        if len(best) != 0:
            for sol, sss in zip(best, iterating_sets):
                sss[1].set_used(sol)
        return [
            best_seq_struct.sequence,
            best_seq_struct.structure,
            best_r.ens_defect,
            best_r.mfe,
        ]

    def __score_design(self, structure, design_structure, r) -> str:
        if r.dot_bracket == structure:
            return "SUCCESS"
        total_score = 0
        barcode_score = 0
        for _, (s1, s2, ds) in enumerate(
            zip(structure, r.dot_bracket, design_structure)
        ):
            if s1 != s2:
                total_score += 1
            if s1 != s2 and ds not in ["(", ")", "."]:
                barcode_score += 1
        if total_score > self.opts.allowed_ss_mismatch:
            return "ss_mismatches"
        if barcode_score > self.opts.allowed_ss_mismatch_barcodes:
            return "ss_mismatch_barcodes"
        return "SUCCESS"

    def __get_build_up(self, num_seqs, build_str, params):
        parser = SequenceStructureSetParser()
        sets = parser.parse(num_seqs, params)
        segments = parse_build_str(build_str)
        pos = 1
        seen = []
        build_up = []
        while True:
            found = False
            for seg_name, seg_pos in segments.items():
                if pos != abs(seg_pos):
                    continue
                if seg_name in seen:
                    continue
                if seg_name not in sets:
                    raise ValueError(f"no set for {seg_name}")
                found = True
                seen.append(seg_name)
                direction = "3PRIME"
                if seg_pos < 0:
                    direction = "5PRIME"
                cur_set = sets[seg_name]
                # hacky way to check if helix
                if cur_set.get_random().sequence.count("&") > 0:
                    direction = "HELIX"
                build_up.append((direction, seg_name, cur_set))
            if not found:
                pos += 1
            if len(seen) == len(segments):
                break

        return build_up

    def __setup_seq_struct(self, seq_struct, build_up):
        iterating_sets = []
        pos = -1
        for direction, seg_name, cur_set in build_up:
            pos += 1
            symbol = self.symbols[pos]
            if direction == "HELIX":
                if len(cur_set) > 1:
                    length = len(cur_set.get_random().split_strands()[0])
                    seq = symbol * length + "&" + symbol * length
                    h_seq_struct = SequenceStructure(seq, seq).split_strands()
                    seq_struct = h_seq_struct[0] + seq_struct + h_seq_struct[1]
                    iterating_sets.append([[symbol * length, symbol * length], cur_set])
                else:
                    h_seq_struct = cur_set.get_random().split_strands()
                    seq_struct = h_seq_struct[0] + seq_struct + h_seq_struct[1]
            elif direction == "5PRIME":
                if len(cur_set) > 1:
                    seq = symbol * len(cur_set.get_random())
                    seq_struct = SequenceStructure(seq, seq) + seq_struct
                    iterating_sets.append([[seq], cur_set])
                else:
                    seq_struct = cur_set.get_random() + seq_struct
            elif direction == "3PRIME":
                if len(cur_set) > 1:
                    seq = symbol * len(cur_set.get_random())
                    seq_struct = seq_struct + SequenceStructure(seq, seq)
                    iterating_sets.append([[seq], cur_set])
                else:
                    seq_struct = seq_struct + cur_set.get_random()
            else:
                raise ValueError(f"unknown direction {direction}")
        return seq_struct, iterating_sets


# design interface to be used with single core or multicore
def design(n_processes, df_sequences, build_str, params, design_opts) -> pd.DataFrame:
    """
    design interface to be used with single core or multicore
    :param n_processes: number of processes to use
    :param df_sequences: dataframe of sequences to design
    :param build_str: build string
    :param params: params
    :param design_opts: design options
    :return: dataframe of designed sequences
    """
    log.info("starting design")
    # need to fix this here and not in the design object as it wont work with
    # multiprocessing
    if "name" not in df_sequences.columns:
        log.info("no 'name' column was in dataframe - adding one")
        df_sequences["name"] = [f"seq_{i}" for i in range(0, len(df_sequences))]
    designer = Designer()
    designer.setup(design_opts)
    # single core run
    if n_processes == 1:
        log.info("running on single core")
        return designer.design(df_sequences, build_str, params)
    # multicore runs
    log.info(f"running on {n_processes} cores with mutliprocessing")
    with multiprocessing.Pool(n_processes) as pool:
        results = pool.starmap(
            designer.design,
            [
                (df_s, build_str, params)
                for df_s in np.array_split(df_sequences, n_processes)
            ],
        )
        dfs = []
        failures = {}
        for r in results:
            dfs.append(r.df_results)
            for key, value in r.failures.items():
                if key in failures:
                    failures[key] += value
                else:
                    failures[key] = value
        return DesignerResults(pd.concat(dfs), failures)


def design_and_save_output(df, output_dir, params):
    os.makedirs(output_dir, exist_ok=True)
    design_opts = DesignOpts(**params["design_opts"])
    yaml.dump(params, open(f"{output_dir}/params.yml", "w"))
    log.info(f"Using parameters:\n{json.dumps(params, indent=4)}")
    results = design(
        params["num_of_processes"],
        df,
        params["build_str"],
        params["segments"],
        design_opts,
    )
    log_failed_design_sequences(results)
    df_results = results.df_results
    write_output_dir(df_results, Path(output_dir))
    if not params["postprocess"]["skip_edit_distance"]:
        edit_dist = calc_edit_distance(df_results)
        log.info(f"the edit distance of lib is: {edit_dist}")
    else:
        log.info("skipping edit distance calculation")
    return df_results


def log_failed_design_sequences(results) -> None:
    """
    add failed sequences to log
    """
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


def write_output_dir(df: pd.DataFrame, output_dir) -> None:
    """
    writes out of the results of a design run to a directory
    :param df: dataframe of results
    """
    if not Path(output_dir).exists():
        raise ValueError(f"output path {output_dir} does not exist")
    log.info(
        f"{output_dir}/results-all.csv contains all information generated from run"
    )
    df.to_csv(f"{output_dir}/results-all.csv", index=False)
    log.info(
        f"{output_dir}/results-rna.csv contains only information related to the RNA sequence"
    )
    df = df[["name", "sequence", "structure", "ens_defect", "mfe"]]
    df.to_csv(f"{output_dir}/results-rna.csv", index=False)
    df_sub = df[["name", "sequence"]].copy()
    # get primer for sequencing
    p5_seq = get_seq_fwd_primer(df_sub)
    if p5_seq is None:
        log.warning("no p5 sequence found")
    else:
        log.info("p5 seq -> " + str(p5_seq))
    df_sub = to_dna(df_sub)
    to_fasta(df_sub, f"{output_dir}/results.fasta")
    df_sub = to_dna_template(df_sub)
    df_sub.to_csv(f"{output_dir}/results-dna.csv", index=False)
    df_sub = df_sub.rename(columns={"name": "Pool name", "sequence": "Sequence"})
    df_sub["Pool name"] = Path(output_dir).stem
    df_sub.to_excel(f"{output_dir}/results-opool.xlsx", index=False)
    df_sub.to_csv(f"{output_dir}/results-opool.csv", index=False)
