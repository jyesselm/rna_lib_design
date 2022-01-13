import pandas as pd
import vienna
from rna_lib_design import logger, structure_set, structure, settings

log = logger.setup_applevel_logger()


def get_p5_from_str(p5_common) -> structure_set.StructureSet:
    common_structs = structure.common_structures()
    p5_sequences = pd.read_csv(settings.RESOURCES_PATH + "/p5_sequences.csv")
    if p5_common is None:
        p5 = common_structs["ref_hairpin_5prime"]
        log.info(f"no p5 sequence supplied using: {p5.sequence}")
    elif p5_common in p5_sequences["name"].unique():
        row = p5_sequences[p5_sequences["name"] == p5_common].iloc[0]
        log.info(f"p5 supplied: {p5_common}")
        log.info(
            f"p5 sequence: {row['sequence']}, structure: {row['structure']}, "
            f"code: {row['code']}"
        )
        p5 = structure.rna_structure(row["sequence"], row["structure"])
    else:
        r = vienna.fold(p5_common)
        log.info(f"p5 sequence supplied: {p5_common}")
        log.info(
            f"p5 sequence has a folded structure of {r.dot_bracket} with "
            f"ens_defect of {r.ensemble_diversity}"
        )
        p5 = structure.rna_structure(p5_common, r.dot_bracket)
    return structure_set.get_single_struct_set(p5, structure_set.AddType.LEFT)


def get_p3_from_str(p3_common) -> structure_set.StructureSet:
    common_structs = structure.common_structures()
    if p3_common is None:
        p3 = common_structs["rt_tail"]
        log.info(f"no p3 sequence supplied using: {p3.sequence}")
    else:
        r = vienna.fold(p3_common)
        log.info(f"p5 sequence supplied: {p3_common}")
        log.info(
            f"p5 sequence has a folded structure of {r.dot_bracket} with "
            f"ens_defect of {r.ensemble_diversity}"
        )
        p3 = structure.rna_structure(p3_common, r.dot_bracket)
    return structure_set.get_single_struct_set(p3, structure_set.AddType.RIGHT)


def get_loop_from_str(loop) -> structure.Structure:
    if loop is None:
        loop_struct = structure.get_common_struct("uucg_loop")
    else:
        r = vienna.fold(loop)
        log.info(f"loop sequence supplied: {loop}")
        log.info(
            f"loop sequence has a folded structure of {r.dot_bracket} with "
            f"ens_defect of {r.ensemble_diversity}"
        )
        loop_struct = structure.rna_structure(loop, r.dot_bracket)
    return loop_struct
