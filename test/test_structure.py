from rna_lib_design.structure import (
    Sequence,
    rna_sequence,
    dna_sequence,
    DotBracket,
)
from rna_lib_design.structure import (
    rna_structure,
    rna_sequence,
    common_structures,
)
import pytest


#######################################################################################
# sequence tests                                                                      #
#######################################################################################


def test_sequence():
    seq = Sequence(list("GGGGUUUUCCCC"))
    assert seq == "GGGGUUUUCCCC"
    seq = Sequence(list("GGGGUUUUCCCC"), seq_type="DNA")
    assert seq == "GGGGTTTTCCCC"


def test_rna_sequence():
    seq = rna_sequence("GGGGUUUUCCCC")
    assert seq.is_rna()
    assert seq.is_dna() == False
    assert seq.get_complement() == "CCCCAAAAGGGG"
    assert seq.get_dna_complement() == "CCCCAAAAGGGG"
    assert seq.get_rna_complement() == "CCCCAAAAGGGG"


def test_dna_sequence():
    seq = dna_sequence("GGGGUUUUCCCC")
    assert seq.is_dna()
    assert str(seq) == "GGGGTTTTCCCC"
    assert seq.get_complement() == "CCCCAAAAGGGG"


def test_insert_sequence():
    seq = rna_sequence("GGGGUUUUCCCC")
    seq2 = seq.insert("A", 0)
    assert seq2 == "AGGGGUUUUCCCC"
    seq2 = seq.insert("A", len(seq) - 1)
    assert seq2 == "GGGGUUUUCCCCA"
    seq2 = seq.insert("A", 1)
    assert seq2 == "GAGGGUUUUCCCC"
    seq2 = seq.insert("AA", 1)
    assert seq2 == "GAAGGGUUUUCCCC"


def test_insert_segments_sequence():
    seq = rna_sequence("GGGGUUUUCCCC")
    seq2 = seq.insert_segments("A", "U", 0, len(seq) - 1)
    assert seq2 == "AGGGGUUUUCCCCU"


def test_str_concat_sequence():
    seq = rna_sequence("GGGGUUUUCCCC")
    seq_str = seq + "A"
    assert seq_str == "GGGGUUUUCCCCA"
    seq_str = "A" + seq
    assert seq_str == "AGGGGUUUUCCCC"


def test_remove_sequence():
    seq = rna_sequence("GGGGUUUUCCCC")
    seq_1 = seq[1:-1]
    assert type(seq_1) == Sequence
    assert seq_1 == "GGGUUUUCCC"


def test_remove_segment_sequence():
    seq = rna_sequence("GGGGUUUUCCCC")
    seq_1 = seq.remove_segment(4, 8)
    assert seq_1 == "GGGGCCCC"


def test_remove_segments():
    seq = rna_sequence("GGGGUUUUCCCC")
    seq_1 = seq.remove_segments([0, 1], [11, 12])
    assert seq_1 == "GGGUUUUCCC"


#######################################################################################
# dot brackests                                                                       #
#######################################################################################


def test_dot_bracket():
    db = DotBracket("((()))")
    assert db.is_paired(0)
    assert db.get_pair_partner(0) == 5
    assert db == "((()))"


def test_add_db():
    db1 = DotBracket("((((")
    db2 = DotBracket("))))")
    assert db1.is_valid() == False
    db = db1 + db2
    assert db.is_valid()
    assert db.is_paired(0)
    assert db.get_pair_partner(1) == 6


def test_remove_db():
    db = DotBracket("((((....))))")
    db_1 = db.remove_segments([0, 1], [11, 12])
    assert db_1 == "(((....)))"


#######################################################################################
# structure tests                                                                     #
#######################################################################################


def test_rna():
    rna_struct = rna_structure("GGGGAAAACCCC", "((((....))))")
    assert rna_struct.is_paired(1)
    assert not rna_struct.is_paired(5)


def test_add():
    rna_1 = rna_structure("GGGGAAAACCCC", "((((....))))")
    rna_2 = rna_structure("GGGGAAAACCCC", "((((....))))")
    rna = rna_1 + rna_2
    assert rna.sequence == "GGGGAAAACCCCGGGGAAAACCCC"
    assert rna.dot_bracket == "((((....))))((((....))))"


def test_insert():
    rna_1 = rna_structure("GGGGAAAACCCC", "((((....))))")
    rna_2 = rna_1.insert("A", ".", 0)
    assert rna_2.sequence == "AGGGGAAAACCCC"
    assert rna_2.dot_bracket == ".((((....))))"


def test_insert_segment():
    rna_1 = rna_structure("GGGGAAAACCCC", "((((....))))")
    rna_2 = rna_1.insert_segments("A", ".", "A", ".", 0, len(rna_1) - 1)
    assert rna_2.sequence == "AGGGGAAAACCCCA"
    assert rna_2.dot_bracket == ".((((....))))."
    rna_2 = rna_1.insert_segments("AA", "..", "AA", "..", 0, len(rna_1) - 1)
    assert rna_2.sequence == "AAGGGGAAAACCCCAA"
    assert rna_2.dot_bracket == "..((((....)))).."


def test_insert_bp():
    rna_1 = rna_structure("GGGGAAAACCCC", "((((....))))")
    rna_2 = rna_1.insert_bps("A&U", 0)
    assert rna_2.sequence == "AGGGGAAAACCCCU"
    assert rna_2.dot_bracket == "(((((....)))))"
    rna_2 = rna_1.insert_bps("GA&UC", 1)


def test_remove():
    rna_1 = rna_structure("GGGGAAAACCCC", "((((....))))")
    rna_2 = rna_1[1:-1]
    assert rna_2.sequence == "GGGAAAACCC"
    assert rna_2.dot_bracket == "(((....)))"

    rna_2 = rna_1.remove_segments([1, 2], [11, 12])
    assert rna_2.sequence == "GGGAAAACCC"
    assert rna_2.dot_bracket == "(((....)))"


def test_remove_bp():
    rna_1 = rna_structure("GGGGAAAACCCC", "((((....))))")
    rna_2 = rna_1.remove_bp(0)
    assert rna_2.sequence == "GGGAAAACCC"
    assert rna_2.dot_bracket == "(((....)))"

    rna_1 = rna_structure("GAGGAAAACCUC", "((((....))))")
    rna_2 = rna_1.remove_bp(1)
    assert rna_2.sequence == "GGGAAAACCC"
    assert rna_2.dot_bracket == "(((....)))"


def test_blank_structure():
    rna_1 = rna_structure("", "")
    assert len(rna_1) == 0
    rna_2 = rna_structure("GGGGAAAACCCC", "((((....))))")
    assert len(rna_2) == 12
    assert len(rna_1 + rna_2) == 12


#######################################################################################
# function tests                                                                      #
#######################################################################################


def test_common_structures():
    common_structs = common_structures()
    struct = common_structs["ref_hairpin_5prime"]
    assert struct.sequence == "GGAAGAUCGAGUAGAUCAAA"
    assert struct.dot_bracket == "....((((.....))))..."
