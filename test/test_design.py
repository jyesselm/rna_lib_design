import pandas as pd
from seq_tools import SequenceStructure
from rna_lib_design.design import (
    parse_build_str,
    Designer,
)
from rna_lib_design.settings import get_resources_path


class TestResources:
    @staticmethod
    def get_complex_params():
        params = {
            "P5": {"name": "ref_hairpin_5prime"},
            "P3": {"name": "rt_tail"},
            "HPBARCODE": {
                "m_type": "HAIRPIN",
                "length": "5",
                "loop_sequence": "CAAAG",
                "loop_structure": "(...)",
            },
            "HBARCODE6": {"m_type": "HELIX", "length": "5-6"},
            "AC": {"sequence": "AC", "structure": ".."},
        }
        return params

    @staticmethod
    def get_simple_sequence_df():
        return pd.DataFrame(
            {
                "sequence": ["GGGAAAACCC"],
                "structure": ["(((....)))"],
            }
        )


class TestParseBuildStr:
    def test_basic_sequence_string(self):
        seq_str = "P5-HPBARCODE-HBARCODE6A-SOI-HBARCODE6B-AC-P3"
        expected_order = {
            "HPBARCODE": -2,
            "HBARCODE6": 1,
            "AC": 2,
            "P5": -3,
            "P3": 3,
        }
        assert parse_build_str(seq_str) == expected_order

    def test_single_helix_sequence_string(self):
        seq_str = "H1A-SOI-H1B"
        expected_order = {"H1": 1}
        assert parse_build_str(seq_str) == expected_order

    def test_multi_helix_sequence_string(self):
        seq_str = "H1A-SOI-AC-H1B"
        expected_order = {"H1": 2, "AC": 1}
        assert parse_build_str(seq_str) == expected_order

    def test_multi_hairpin_sequence_string(self):
        seq_str = "HP1-SOI-HP2-AC"
        expected_order = {"HP1": -1, "HP2": 1, "AC": 2}
        assert parse_build_str(seq_str) == expected_order

    def test_helix_with_a_and_b(self):
        seq_str = "H1A-SOI-H1B"
        expected_order = {"H1": 1}
        assert parse_build_str(seq_str) == expected_order

    def test_multi_helix_with_a_and_b(self):
        seq_str = "H1A-H2A-SOI-H2B-AC-H1B"
        expected_order = {"H1": 3, "H2": 1, "AC": 2}
        assert parse_build_str(seq_str) == expected_order


def test_design():
    build_str = "P5-HPBARCODE-HBARCODE6A-SOI-HBARCODE6B-AC-P3"
    params = TestResources.get_complex_params()
    df_sequences = TestResources.get_simple_sequence_df()
    designer = Designer()
    df_results = designer.design(df_sequences, build_str, params)
    #row = df_results.iloc[0]
    #print()
    #print(row)

