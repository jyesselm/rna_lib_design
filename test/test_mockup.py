from seq_tools import SequenceStructure
from rna_lib_design.mockup import (
    get_design_setup,
    design,
    parse_sequence_string,
    parse_sequence_structure_sets,
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


def test_parse_sequence_structure_sets():
    params = {
        "H1": {"m_type": "HELIX", "gu": True, "length": "5-6"},
        "SS1": {"m_type": "SSTRAND", "length": "5"},
        "HP1": {
            "m_type": "HAIRPIN",
            "length": "5",
            "sequence": "CAAAG",
            "structure": "(...)",
        },
    }
    sss_dict = parse_sequence_structure_sets(10, params)
    assert len(sss_dict) == 3
    assert len(sss_dict["H1"].seqstructs) == 27
    assert len(sss_dict["SS1"].seqstructs) == 32
    assert len(sss_dict["HP1"].seqstructs) == 11


def test_get_design_setup():
    seq_str = "P5-HPBARCODE-HBARCODE6A-SOI-HBARCODE6B-AC-P3"
    params = {
        "P5": {"name": "ref_hairpin_5prime"},
        "P3": {"name": "rt_tail"},
        "HPBARCODE": {
            "m_type": "HAIRPIN",
            "length": "5",
            "sequence": "CAAAG",
            "structure": "(...)",
        },
        "HBARCODE6": {"m_type": "HELIX", "length": "5-6"},
        "AC": {"sequence": "AC", "structure": ".."},
    }
    build_up = get_design_setup(10, seq_str, params)
    assert len(build_up) == 5
    assert build_up[0][1] == "HBARCODE6"
    assert build_up[1][1] == "HPBARCODE"
    assert build_up[2][1] == "AC"
    assert build_up[3][1] == "P5"
    assert build_up[4][1] == "P3"


def test_design():
    seq_str = "P5-HPBARCODE-HBARCODE6A-SOI-HBARCODE6B-AC-P3"
    params = {
        "P5": {"name": "ref_hairpin_5prime"},
        "P3": {"name": "rt_tail"},
        "HPBARCODE": {
            "m_type": "HAIRPIN",
            "length": "5",
            "sequence": "CAAAG",
            "structure": "(...)",
        },
        "HBARCODE6": {"m_type": "HELIX", "length": "5-6"},
        "AC": {"sequence": "AC", "structure": ".."},
    }
    build_up = get_design_setup(10, seq_str, params)
    seq_struct = SequenceStructure("GGGGAAAACCCC", "((((....))))")
    output = design(seq_struct, build_up)
