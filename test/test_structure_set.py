import pandas as pd

from rna_lib_design.settings import get_resources_path, get_test_path
from rna_lib_design.structure_set import (
    SequenceStructure,
    SequenceStructureSet,
    SequenceStructureSetParser,
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
        csv_path = (
            get_resources_path() / "barcodes/helices/len_1/md_0_gu_0_0.csv"
        )
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
        assert len(set_dict['H1']) == 27

    def test_sstrand(self):
        params = {"SS1": {"m_type": "SSTRAND", "length": "5"}}
        set_dict = self.parser.parse(10, params)
        assert len(set_dict["SS1"]) == 32

    def test_hairpin(self):
        params = {
            "HP1": {
                "m_type": "HAIRPIN",
                "length": "5",
                "sequence": "CAAAG",
                "structure": "(...)",
            }
        }
        set_dict = self.parser.parse(10, params)
        assert len(set_dict['HP1']) == 11
        assert len(set_dict['HP1'].get_random()) == 18

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
            "H1" : {"m_type": "HELIX", "gu": True, "length": "5-6"},
            "SS1": {"m_type": "SSTRAND", "length": "5"},
            "HP1": {
                "m_type"   : "HAIRPIN",
                "length"   : "5",
                "sequence" : "CAAAG",
                "structure": "(...)",
            },
        }
        set_dict = self.parser.parse(10, params)
        assert len(set_dict) == 3
        assert len(set_dict["H1"].seqstructs) == 27
        assert len(set_dict["SS1"].seqstructs) == 32
        assert len(set_dict["HP1"].seqstructs) == 11
