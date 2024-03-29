import pandas as pd
import pytest

from rna_lib_design.settings import get_resources_path, get_test_path
from rna_lib_design.structure_set import (
    SequenceStructure,
    SequenceStructureSet,
    SequenceStructureSetParser,
    get_named_seq_structs,
    get_named_seq_struct,
    get_optimal_sstrand_set,
    get_optimal_helix_set,
    get_optimal_hairpin_set,
)

TEST_RESOURCES = get_test_path() / "resources"


class TestResources:
    @staticmethod
    def get_test_sstrand():
        csv_path = TEST_RESOURCES / "barcodes" / "sstrand.csv"
        return SequenceStructureSet.from_csv(csv_path)


class TestSequenceStructureSet:
    def test_init(self):
        ss1 = SequenceStructure("ATCG", "((((")
        ss2 = SequenceStructure("CGAT", "))))")
        sss = SequenceStructureSet([ss1, ss2])
        assert len(sss.seqstructs) == 2
        assert not any(sss.used)
        assert not sss.allow_duplicates

    def test_from_single(self):
        ss = SequenceStructure("ATCG", "((((")
        sss = SequenceStructureSet.from_single(ss)
        assert len(sss.seqstructs) == 1
        assert not any(sss.used)
        assert sss.allow_duplicates

    def test_get_random(self):
        ss1 = SequenceStructure("ATCG", "((((")
        ss2 = SequenceStructure("CGAT", "))))")
        sss = SequenceStructureSet([ss1, ss2])
        sec_struct = sss.get_random()
        assert sec_struct in [ss1, ss2]
        sss.set_used(ss1)
        assert any(sss.used)

    def test_set_used(self):
        ss1 = SequenceStructure("ATCG", "((((")
        ss2 = SequenceStructure("CGAT", "))))")
        sss = SequenceStructureSet([ss1, ss2])
        sss.set_used(ss1)
        assert sss.used[0]

    def test_set_last_used(self):
        ss1 = SequenceStructure("ATCG", "((((")
        ss2 = SequenceStructure("CGAT", "))))")
        sss = SequenceStructureSet([ss1, ss2])
        sss.get_random()
        sss.set_last_used()
        assert sss.used[sss.last]

    def test_from_csv(self):
        csv_path = get_resources_path() / "barcodes/helices/len_1/md_0_gu_0_0.csv"
        sss = SequenceStructureSet.from_csv(csv_path)
        assert len(sss.seqstructs) == 4
        assert not any(sss.used)
        assert not sss.allow_duplicates

    # old tests from version 1.0
    def test_sstrand_structure_set(self):
        struct_set = TestResources.get_test_sstrand()
        struct_1 = struct_set.get_random()
        struct_set.set_last_used()
        struct_2 = struct_set.get_random()
        # should not recieve the same structure twice
        assert struct_1 != struct_2
        seq_struct = SequenceStructure("GGGAAAACCC", "(((....)))")
        seq_struct_applied = struct_1 + seq_struct
        assert seq_struct_applied.structure == "......(((....)))"


class TestSequenceStructureSetFromParams:
    @classmethod
    def setup_class(cls):
        cls.parser = SequenceStructureSetParser()

    def test_helix(self):
        params = {"H1": {"m_type": "HELIX", "gu": True, "length": "5-6"}}
        set_dict = self.parser.parse(10, params)
        assert len(set_dict["H1"]) == 27

    def test_sstrand(self):
        params = {"SS1": {"m_type": "SSTRAND", "length": "5"}}
        set_dict = self.parser.parse(10, params)
        assert len(set_dict["SS1"]) == 32

    def test_hairpin(self):
        params = {
            "HP1": {
                "m_type": "HAIRPIN",
                "length": "5",
                "loop_seq": "CAAAG",
                "loop_ss": "(...)",
            }
        }
        set_dict = self.parser.parse(10, params)
        assert len(set_dict["HP1"]) == 11
        assert len(set_dict["HP1"].get_random()) == 18

    def test_single(self):
        params = {
            "SS1": {
                "sequence": "CAAAG",
                "structure": "(...)",
            }
        }
        set_dict = self.parser.parse(10, params)
        assert len(set_dict["SS1"]) == 1
        assert len(set_dict["SS1"].get_random()) == 5

    def test_complex_parse(self):
        params = {
            "H1": {"m_type": "HELIX", "gu": True, "length": "5-6"},
            "SS1": {"m_type": "SSTRAND", "length": "5"},
            "HP1": {
                "m_type": "HAIRPIN",
                "length": "5",
                "loop_seq": "CAAAG",
                "loop_ss": "(...)",
            },
        }
        set_dict = self.parser.parse(10, params)
        assert len(set_dict) == 3
        assert len(set_dict["H1"].seqstructs) == 27
        assert len(set_dict["SS1"].seqstructs) == 32
        assert len(set_dict["HP1"].seqstructs) == 11


class TestNamedSequenceStructure:
    def test_get_all(self):
        df = get_named_seq_structs()
        df_sub = df[df["name"] == "uucg_p5_rev_primer"]
        assert len(df_sub) == 1

    def test_get_one(self):
        ss = get_named_seq_struct("uucg_p5_rev_primer")
        assert ss.sequence == "GGAACAGCACUUCGGUGCAAA"
        assert ss.structure == "......((((....))))..."

    def test_get_one_not_found(self):
        with pytest.raises(ValueError):
            get_named_seq_struct("not_a_real_name")


def test_get_optimal_sstrand_set():
    sset = get_optimal_sstrand_set(5, 10)
    assert len(sset) == 32


def test_split_set():
    csv_path = get_resources_path() / "barcodes/helices/len_1/md_0_gu_0_0.csv"
    sss = SequenceStructureSet.from_csv(csv_path)
    sets = sss.split(2)
    assert len(sets) == 2
    assert len(sets[0]) == 2


def test_split_single_set():
    seq_struct = SequenceStructure("GGGAAAACCC", "(((....)))")
    sss = SequenceStructureSet.from_single(seq_struct)
    sets = sss.split(2)
    sets[0].set_used(sets[0].get_random())
    assert sets[0].num_available() == 1
    assert sets[0].num_used() == 0
