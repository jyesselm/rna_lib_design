from seq_tools import SequenceStructure
from rna_lib_design.mockup import (
    parse_sequence_string,
    sequence_structure_set_from_params,
    SequenceStructureSet,
)
from rna_lib_design.settings import get_resources_path


class TestParseSequenceString:
    def test_basic_sequence_string(self):
        seq_str = "P5-HPBARCODE-HBARCODE6A-SOI-HBARCODE6B-AC-P3"
        expected_order = {
            "HPBARCODE": -2,
            "HBARCODE6": 1,
            "AC": 2,
            "P5": -3,
            "P3": 3,
        }
        assert parse_sequence_string(seq_str) == expected_order

    def test_single_helix_sequence_string(self):
        seq_str = "H1A-SOI-H1B"
        expected_order = {"H1": 1}
        assert parse_sequence_string(seq_str) == expected_order

    def test_multi_helix_sequence_string(self):
        seq_str = "H1A-SOI-AC-H1B"
        expected_order = {"H1": 2, "AC": 1}
        assert parse_sequence_string(seq_str) == expected_order

    def test_multi_hairpin_sequence_string(self):
        seq_str = "HP1-SOI-HP2-AC"
        expected_order = {"HP1": -1, "HP2": 1, "AC": 2}
        assert parse_sequence_string(seq_str) == expected_order

    def test_helix_with_a_and_b(self):
        seq_str = "H1A-SOI-H1B"
        expected_order = {"H1": 1}
        assert parse_sequence_string(seq_str) == expected_order

    def test_multi_helix_with_a_and_b(self):
        seq_str = "H1A-H2A-SOI-H2B-AC-H1B"
        expected_order = {"H1": 3, "H2": 1, "AC": 2}
        assert parse_sequence_string(seq_str) == expected_order


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


class TestSequenceStructureSetFromParams:
    def test_helix(self):
        params = {"H1": {"m_type": "HELIX", "gu": True, "length": "5-6"}}
        sss = sequence_structure_set_from_params(10, params['H1'])
        assert len(sss.seqstructs) == 27

    def test_sstrand(self):
        params = {"SS1": {"m_type": "SSTRAND", "length" : "5"}}
        sss = sequence_structure_set_from_params(10, params['SS1'])

