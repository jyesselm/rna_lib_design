import logging

import re
import click
import pandas as pd
import numpy as np
from dataclasses import dataclass
from typing import List, Dict
import rna_secstruct as rl

import vienna
from seq_tools import to_dna
from rna_lib_design import structure_set, logger, structure, util
from rna_lib_design.structure import rna_structure

log = logger.setup_applevel_logger()


################################################################################
# classes and functions                                                        #
################################################################################


@dataclass(frozen=True, order=True)
class DesignResults:
    ens_defect: float
    design: structure.Structure


@dataclass(frozen=True, order=True)
class DesignOptions:
    max_ens_defect: float = 2.0
    max_design_attempts: int = 100
    max_design_solutions: int = 10
    use_structure: bool = True


def get_best_design(
    sd_dicts: List[structure_set.StructureSet],
    struct: structure.Structure,
    opts: DesignOptions,
    org_ens_defect: float = 999.0,
) -> DesignResults:
    count = 0
    best_score = 1000
    best = None
    best_uses = []
    if org_ens_defect > 100:
        if len(struct) > 5:
            org_ens_defect = vienna.fold(struct).ens_defect
        else:
            org_ens_defect = 1
    while 1:
        count += 1
        final_struct = structure_set.apply(sd_dicts, struct)
        vr = vienna.fold(final_struct.sequence)
        if count > opts.max_design_attempts:
            return DesignResults(999, structure.rna_structure("A", "."))
        # score = final_struct.dot_bracket_difference(structure.DotBracket(vr.dot_bracket))
        # print(final_struct.dot_bracket, vr.ensemble_diversity, score)
        if final_struct.dot_bracket != vr.dot_bracket:
            continue
        if opts.max_ens_defect < vr.ens_defect - org_ens_defect:
            continue
        if vr.ens_defect < best_score:
            best_score = vr.ens_defect
            best = final_struct
            best_uses = [x.get_current_pos() for x in sd_dicts]
        if count > opts.max_design_solutions:
            break
    if opts.use_structure:
        for i, bu in enumerate(best_uses):
            sd_dicts[i].set_used(bu)
    return DesignResults(best_score, best)


def get_best_designs_in_dataframe(
    df: pd.DataFrame,
    sets: List[structure_set.StructureSet],
    opts: DesignOptions,
    max: int = 1000000,
) -> pd.DataFrame:
    df["org_sequence"] = df["sequence"]
    df["org_structure"] = df["structure"]
    if "ens_defect" in df:
        df["org_ens_defect"] = df["ens_defect"]
    else:
        df["ens_defect"] = np.nan
    count = 0
    for i, row in df.iterrows():
        if i > max:
            break
        count += 1
        if count % 100 == 0:
            log.info(f"barcoded {count} constructs!")
        struct = structure.rna_structure(row["sequence"], row["structure"])
        if "org_ens_defect" in df:
            sol = get_best_design(sets, struct, opts, row["org_ens_defect"])
        else:
            sol = get_best_design(sets, struct, opts)
        if sol.ens_defect > 100:
            log.error(
                f"unable to find a design for {row['name']} below the cutoff "
                f"{opts.max_ens_defect}"
            )
            df.at[i, ["sequence", "structure", "ens_defect"]] = [
                np.nan,
                np.nan,
                np.nan,
            ]
        else:
            log.debug(f"solution found: name={row['name']}, sol={sol},")
            df.at[i, ["sequence", "structure", "ens_defect"]] = [
                str(sol.design.sequence),
                str(sol.design.dot_bracket),
                sol.ens_defect,
            ]
    if max < len(df):
        df = df[:max]
    return df


def write_results_to_file(
    df: pd.DataFrame, fname="results", opool_name="opool"
) -> None:
    log.info(f"{fname}-all.csv contains all information generated from run")
    df.to_csv(f"{fname}-all.csv", index=False)
    log.info(
        f"{fname}-rna.csv contains only information related to the RNA sequence"
    )
    df.to_csv(f"{fname}-rna.csv", index=False)
    df_sub = df[["name", "sequence"]].copy()
    df_sub = to_dna(df_sub)
    f = open(f"{fname}.fasta", "w")
    for i, row in df_sub.iterrows():
        f.write(f">{row['name']}\n{row['sequence']}\n")
    f.close()
    edit_dist = util.compute_edit_distance(df_sub)
    log.info(f"the edit distance of lib is: {edit_dist}")
    df_sub["sequence"] = [
        "TTCTAATACGACTCACTATA" + seq for seq in df_sub["sequence"]
    ]
    p5_seq = util.indentify_p5_sequence(df_sub["sequence"])
    fwd_primer = util.indentify_fwd_primer(df_sub["sequence"])
    rev_primer = util.indentify_rev_primer(df_sub["sequence"])
    log.info("p5 seq -> " + str(p5_seq))
    log.info("fwd primer -> " + str(fwd_primer))
    log.info("rev primer -> " + str(rev_primer))
    df_sub.to_csv(f"{fname}-dna.csv", index=False)
    df_sub = df_sub.rename(
        columns={"name": "Pool name", "sequence": "Sequence"}
    )
    df_sub["Pool name"] = opool_name
    df_sub.to_excel(f"{fname}-opool.xlsx", index=False)
    df_sub.to_csv(f"{fname}-opool.csv", index=False)


def str_to_range(x):
    return sum(
        (
            i if len(i) == 1 else list(range(i[0], i[1] + 1))
            for i in (
                [int(j) for j in i if j]
                for i in re.findall("(\d+),?(?:-(\d+))?", x)
            )
        ),
        [],
    )


class HelixRandomizer(object):
    def __init__(self):
        pass

    def __randomize_helix(self, h, exclude):
        for i in range(100):
            org_seq = h.sequence.split("&")
            org_s = util.compute_stretches(org_seq[0], org_seq[1])
            seq1, seq2 = self.__generate_helix_sequence(h, exclude)
            new_s = util.compute_stretches(seq1, seq2)
            if (
                org_s.max_gc_stretch > 3
                and new_s.max_gc_stretch > org_s.max_gc_stretch
            ):
                continue
            elif new_s.max_gc_stretch > 3:
                continue
            if (
                org_s.max_stretch_1 > 3
                and new_s.max_stretch_1 > org_s.max_stretch_1
            ):
                continue
            elif new_s.max_stretch_1 > 3:
                continue
            if (
                org_s.max_stretch_2 > 3
                and new_s.max_stretch_2 > org_s.max_stretch_2
            ):
                continue
            elif new_s.max_stretch_2 > 3:
                continue
            return seq1 + "&" + seq2

    def __generate_helix_sequence(self, h, exclude):
        if exclude is None:
            exclude = []
        strand1, strand2 = h.strands
        seq1, seq2 = "", ""
        for i, (s1, s2) in enumerate(zip(strand1, strand2[::-1])):
            # dont change end pairs
            if i == 0 or i == len(strand1) - 1:
                seq1 += self.sequence[s1]
                seq2 += self.sequence[s2]
            elif s1 in exclude or s2 in exclude:
                seq1 += self.sequence[s1]
                seq2 += self.sequence[s2]
            else:
                bp = util.random_weighted_basepair()
                seq1 += bp[0]
                seq2 += bp[1]
        return seq1, seq2[::-1]

    def run(self, sequence, structure, exclude=None, exclude_seqs=None):
        sequence = sequence.replace("T", "U")
        # log.debug("excluded: " + str(exclude))
        if exclude is None:
            exclude = []
        if exclude_seqs is not None:
            for es in exclude_seqs:
                pos = sequence.find(es)
                if pos != -1:
                    exclude.extend(list(range(pos, pos + len(es) + 1)))
        # log.debug("final excluded: " + str(exclude))
        self.sequence, self.structure = sequence, structure
        s = rl.SecStruct(sequence, structure)
        best = 1000
        best_seq = ""
        for _ in range(100):
            for h in s:
                if not h.is_helix():
                    continue
                new_seq = self.__randomize_helix(h, exclude)
                s.change_motif(h.m_id, new_seq, h.structure)
            r = vienna.fold(s.sequence)
            if r.dot_bracket != structure:
                continue
            if r.ens_defect < best:
                best = r.ens_defect
                best_seq = s.sequence
        if len(best_seq) == 0:
            return DesignResults(best, rna_structure("A", "."))
        else:
            return DesignResults(best, rna_structure(best_seq, structure))


def replace_junction(df, seq, ss):
    pass
