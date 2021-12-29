import vienna
import pandas as pd
from rna_lib_design import structure_set, logger, structure
from dataclasses import dataclass
from typing import List

log = logger.setup_applevel_logger()


@dataclass(frozen=True, order=True)
class DesignResults:
    ens_defect: float
    design: structure.Structure


def get_best_design(
    sd_dicts: List[structure_set.StructureSet],
    struct: structure.Structure,
    cutoff=5.0,
    max_count=100,
    keep_best=10,
) -> DesignResults:
    count = 0
    best_score = 1000
    best = None
    best_uses = []
    while 1:
        count += 1
        final_struct = structure_set.apply(sd_dicts, struct)
        vr = vienna.fold(final_struct.sequence)
        if count > max_count:
            return None
        # score = final_struct.dot_bracket_difference(structure.DotBracket(vr.dot_bracket))
        # print(final_struct.dot_bracket, vr.ensemble_diversity, score)
        if final_struct.dot_bracket != vr.dot_bracket:
            continue
        if cutoff < vr.ensemble_diversity:
            continue
        if vr.ensemble_diversity < best_score:
            best_score = vr.ensemble_diversity
            best = final_struct
            best_uses = [x.get_current_pos() for x in sd_dicts]
        if count > keep_best:
            break
    for i, bu in enumerate(best_uses):
        sd_dicts[i].set_used(bu)
    return DesignResults(best_score, best)


def write_results_to_csv(constructs, fname="final"):
    f_rna = open(fname + "_rna.csv", "w")
    f_rna.write("name,sequence,structure,ens_defect\n")
    f_dna = open(fname + "_dna.csv", "w")
    f_dna.write("name,sequence\n")
    for c in constructs:
        f_rna.write(
            c[0]
            + ","
            + c[1].sequence.str()
            + ","
            + str(c[1].dot_bracket)
            + ","
            + str(c[2])
            + "\n"
        )
        f_dna.write(
            c[0] + "," + "TTCTAATACGACTCACTATA" + c[1].sequence.to_dna().str() + "\n"
        )
    f_rna.close()
    f_dna.close()


def helix_barcode(
    df: pd.DataFrame, h_length: int, p5: structure.Structure, p3: structure.Structure,
) -> pd.DataFrame:
    h_set = structure_set.get_optimal_helix_set(h_length, len(df) * 1.1)


def cli():
    pass


if __name__ == "__main__":
    cli()
