import pandas as pd
from rna_lib_design import structure_set, settings, structure, testing


def test_sstrand_structure_set():
    struct_set = testing.get_test_hairpin()
    struct_1 = struct_set.get_random()[0]
    struct_set.set_used()
    struct_2 = struct_set.get_random()[0]
    # should not recieve the same structure twice
    assert struct_1 != struct_2
    rna_struct = structure.rna_structure("GGGAAAACCC", "(((....)))")
    rna_struct_new = struct_set.apply_random(rna_struct)
    pos = struct_set.get_current_pos()
    struct_3 = struct_set.get(pos)[0]
    assert rna_struct_new == rna_struct + struct_3
    assert rna_struct_new == struct_set.apply(rna_struct, pos)


def test_helix_structure_set():
    struct_set = testing.get_test_helices()
    structs = struct_set.get_random()
    struct_1 = structs[0] + structs[1]
    struct_set.set_used()
    structs = struct_set.get_random()
    struct_2 = structs[0] + structs[1]
    assert struct_1 != struct_2
    rna_struct = structure.rna_structure("GGGAAAACCC", "(((....)))")
    rna_struct_new = struct_set.apply_random(rna_struct)
    pos = struct_set.get_current_pos()
    structs = struct_set.get(pos)
    assert rna_struct_new == structs[0] + rna_struct + structs[1]
    assert rna_struct_new == struct_set.apply(rna_struct, pos)


def test_hairpin_structure_set():
    struct_set = testing.get_test_hairpin()
    struct_1 = struct_set.get_random()[0]
    struct_set.set_used()
    struct_2 = struct_set.get_random()[0]
    assert struct_1 != struct_2
    rna_struct = structure.rna_structure("GGGAAAACCC", "(((....)))")
    rna_struct_new = struct_set.apply_random(rna_struct)
    pos = struct_set.get_current_pos()
    struct_3 = struct_set.get(pos)[0]
    assert rna_struct_new == rna_struct + struct_3
    assert rna_struct_new == struct_set.apply(rna_struct, pos)


def test_hairpin_structure_set_no_buffer():
    hp_set = testing.get_test_hairpin()
    hp_set.set_buffer(None)
    rna_struct = structure.rna_structure("GGGAAAACCC", "(((....)))")
    rna_struct_new = hp_set.apply(rna_struct, 0)
    assert rna_struct_new.dot_bracket == "(((....)))(((((((....)))))))"


def test_single_structure_set():
    rna_struct = structure.rna_structure("GGGAAAACCC", "(((....)))")
    struct_set = structure_set.get_single_struct_set(
        rna_struct, structure_set.AddType.RIGHT
    )
    struct_1 = struct_set.get_random()[0]
    struct_set.set_used()
    struct_2 = struct_set.get_random()[0]
    assert struct_1 == struct_2
    pos = struct_set.get_current_pos()
    struct_3 = struct_set.get(pos)[0]
    rna_struct_new = struct_set.apply_random(rna_struct)
    assert rna_struct_new == rna_struct + struct_3
    assert rna_struct_new == struct_set.apply(rna_struct, pos)


def test_apply():
    rna_struct = structure.rna_structure("GGGAAAACCC", "(((....)))")
    struct_set = testing.get_test_helices()
    new_struct = structure_set.apply([struct_set], rna_struct)
    assert rna_struct != new_struct


def test_apply_2():
    rna_struct = structure.rna_structure("GGGAAAACCC", "(((....)))")
    sets = [testing.get_test_helices(), testing.get_test_helices()]
    new_struct = structure_set.apply(sets, rna_struct)
    assert len(new_struct) == 34


def test_apply_blank():
    rna_struct = structure.rna_structure("", "")
    struct_set = testing.get_test_helices()
    new_struct = structure_set.apply([struct_set], rna_struct)


def test_get_helices():
    struct_set = structure_set.get_optimal_helix_set(5, 10)
    assert len(struct_set) > 10
    structs = struct_set.get(0)
    assert len(structs[0]) == 5


def test_get_hairpins():
    loop = structure.rna_structure("GAAAAC", "(....)")
    hp_set = structure_set.get_optimal_hairpin_set(
        5, loop, 10, structure_set.AddType.RIGHT
    )
    assert len(hp_set) > 10
    struct = hp_set.get(0)[0]
    assert len(struct) == 19
    hp_set.set_buffer(None)
    struct = hp_set.get(0)[0]
    assert len(struct) == 16


def test_get_sstrands():
    struct_set = structure_set.get_optimal_sstrand_set(5, 10)
    assert len(struct_set) > 10
    struct = struct_set.get(0)[0]
    assert len(struct) == 5
    assert struct.dot_bracket == "....."


def get_common_seq_structure_set():
    struct_set = structure_set.get_common_seq_structure_set(
        "ref_hairpin_5prime"
    )
    assert len(struct_set) == 1
