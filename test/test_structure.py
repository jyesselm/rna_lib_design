from rna_lib_design.structure import Sequence, rna_sequence, dna_sequence, DotBracket
from rna_lib_design.structure import Structure, rna_structure, rna_sequence
import pytest

#######################################################################################
# sequence tests                                                                      #
#######################################################################################


def test_sequence():
    seq = Sequence(list("GGGGUUUUCCCC"))
    assert seq == "GGGGUUUUCCCC"
    seq = Sequence(list("GGGGUUUUCCCC"), seq_type="DNA")
    assert seq == "GGGGTTTTCCCC"


def test_rna():
    seq = rna_sequence("GGGGUUUUCCCC")
    assert seq.is_rna()
    assert seq.is_dna() == False
    assert seq.get_complement() == "CCCCAAAAGGGG"
    assert seq.get_dna_complement() == "CCCCAAAAGGGG"
    assert seq.get_rna_complement() == "CCCCAAAAGGGG"


def test_dna():
    seq = dna_sequence("GGGGUUUUCCCC")
    assert seq.is_dna()
    assert str(seq) == "GGGGTTTTCCCC"
    assert seq.get_complement() == "CCCCAAAAGGGG"


def test_insert():
    seq = rna_sequence("GGGGUUUUCCCC")
    seq2 = seq.insert("A", 0)
    assert seq2 == "AGGGGUUUUCCCC"
    seq2 = seq.insert("A", len(seq) - 1)
    assert seq2 == "GGGGUUUUCCCCA"
    seq2 = seq.insert("A", 1)
    assert seq2 == "GAGGGUUUUCCCC"
    seq2 = seq.insert("AA", 1)
    assert seq2 == "GAAGGGUUUUCCCC"


def test_insert_segments():
    seq = rna_sequence("GGGGUUUUCCCC")
    seq2 = seq.insert_segments("A", "U", 0, len(seq) - 1)
    assert seq2 == "AGGGGUUUUCCCCU"


"""
class SequenceUnittests(unittest.TestCase):
     

    

    def test_remove(self):
        seq = rna_sequence("GGGGUUUUCCCC")
        seq_1 = seq[1:-1]
        self.assertTrue(type(seq_1) == Sequence)

    def test_remove_segment(self):
        seq = rna_sequence("GGGGUUUUCCCC")
        seq_1 = seq.remove_segment(4, 8)
        self.assertTrue(seq_1 == "GGGGCCCC")

    def test_remove_segments(self):
        seq = rna_sequence("GGGGUUUUCCCC")
        seq_1 = seq.remove_segments([0, 1], [11, 12])
        self.assertTrue(seq_1 == "GGGUUUUCCC")


class DotBracketUnittests(unittest.TestCase):
    def test(self):
        db = DotBracket("((()))")
        self.assertTrue(db.is_paired(0))
        self.assertTrue(db.get_pair_partner(0) == 5)
        self.assertTrue(db == "((()))")

    def test_add(self):
        db1 = DotBracket("((((")
        db2 = DotBracket("))))")
        db = db1 + db2
        self.assertTrue(db.is_valid())
        self.assertTrue(db.is_paired(0))
        self.assertTrue(db.get_pair_partner(1) == 6)

    def test_remove(self):
        db = DotBracket("((((....))))")
        db_1 = db.remove_segments([0, 1], [11, 12])
        self.assertTrue(db_1 == "(((....)))")

    def test_get_helical_segments(self):
        db_1 = DotBracket("((((....))))")
        print(db_1.get_helical_segments())

        db_2 = DotBracket("(((...))((...)))")
        # print(db_2.get_helical_segments())


class StructureUnittests(unittest.TestCase):
    def test_rna(self):
        rna_struct = rna_structure("GGGGAAAACCCC", "((((....))))")
        self.assertTrue(rna_struct.is_paired(1))
        self.assertFalse(rna_struct.is_paired(5))

    def test_add(self):
        rna_1 = rna_structure("GGGGAAAACCCC", "((((....))))")
        rna_2 = rna_structure("GGGGAAAACCCC", "((((....))))")
        rna = rna_1 + rna_2
        self.assertTrue(rna.sequence == "GGGGAAAACCCCGGGGAAAACCCC")
        self.assertTrue(rna.dot_bracket == "((((....))))((((....))))")

    def test_insert(self):
        rna_1 = rna_structure("GGGGAAAACCCC", "((((....))))")
        rna_2 = rna_1.insert("A", ".", 0)
        self.assertTrue(rna_2.sequence == "AGGGGAAAACCCC")
        self.assertTrue(rna_2.dot_bracket == ".((((....))))")

    def test_insert_segment(self):
        rna_1 = rna_structure("GGGGAAAACCCC", "((((....))))")
        rna_2 = rna_1.insert_segments("A", ".", "A", ".", 0, len(rna_1) - 1)
        self.assertTrue(rna_2.sequence == "AGGGGAAAACCCCA")
        self.assertTrue(rna_2.dot_bracket == ".((((....)))).")
        rna_2 = rna_1.insert_segments("AA", "..", "AA", "..", 0, len(rna_1) - 1)
        self.assertTrue(rna_2.sequence == "AAGGGGAAAACCCCAA")
        self.assertTrue(rna_2.dot_bracket == "..((((....))))..")

    def test_insert_bp(self):
        rna_1 = rna_structure("GGGGAAAACCCC", "((((....))))")
        rna_2 = rna_1.insert_bps("A&U", 0)
        self.assertTrue(rna_2.sequence == "AGGGGAAAACCCCU")
        self.assertTrue(rna_2.dot_bracket == "(((((....)))))")
        rna_2 = rna_1.insert_bps("GA&UC", 1)

    def test_remove(self):
        rna_1 = rna_structure("GGGGAAAACCCC", "((((....))))")
        rna_2 = rna_1[1:-1]
        self.assertTrue(rna_2.sequence == "GGGAAAACCC")
        self.assertTrue(rna_2.dot_bracket == "(((....)))")

        rna_2 = rna_1.remove_segments([1, 2], [11, 12])
        self.assertTrue(rna_2.sequence == "GGGAAAACCC")
        self.assertTrue(rna_2.dot_bracket == "(((....)))")

    def test_remove_bp(self):
        rna_1 = rna_structure("GGGGAAAACCCC", "((((....))))")
        rna_2 = rna_1.remove_bp(0)
        self.assertTrue(rna_2.sequence == "GGGAAAACCC")
        self.assertTrue(rna_2.dot_bracket == "(((....)))")

        rna_1 = rna_structure("GAGGAAAACCUC", "((((....))))")
        rna_2 = rna_1.remove_bp(1)
        self.assertTrue(rna_2.sequence == "GGGAAAACCC")
        self.assertTrue(rna_2.dot_bracket == "(((....)))")
"""
