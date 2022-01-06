import vienna
from rna_lib_design.design import get_best_design
from rna_lib_design import logger, structure_set, structure

log = logger.setup_applevel_logger()


def get_p5_from_str(p5_common) -> structure_set.StructureSet:
    common_structs = structure.common_structures()
    if p5_common is None:
        p5 = common_structs["ref_hairpin_5prime"]
        log.info(f"no p5 sequence supplied using: {p5.sequence}")
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
