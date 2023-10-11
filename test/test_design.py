import pandas as pd
from seq_tools import SequenceStructure
from rna_lib_design.design import (
    parse_build_str,
    get_seq_struct_designer,
    design,
    Designer,
    DesignOpts,
)
from rna_lib_design.settings import get_resources_path, get_test_path


class TestResources:
    @staticmethod
    def get_complex_params():
        params = {
            "P5": {"name": "org_minittr_pool_rev_seq_primer"},
            "P3": {"name": "rt_tail"},
            "HPBARCODE": {
                "m_type": "HAIRPIN",
                "length": "5",
                "loop_seq": "CAAAG",
                "loop_ss": "(...)",
            },
            "HBARCODE6": {"m_type": "HELIX", "length": "6"},
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


class TestSequencerDesigner:
    def test_get_sequence_designer(self):
        build_str = "P5-HPBARCODE-HBARCODE6A-SOI-HBARCODE6B-AC-P3"
        params = TestResources.get_complex_params()
        df_sequences = TestResources.get_simple_sequence_df()
        sd = get_seq_struct_designer(1, build_str, params)
        seq_struct = df_sequences.iloc[0]
        d_seq_struct = sd.get_designable_seq_struct(seq_struct)
        assert (
            d_seq_struct.sequence
            == "GGAAGAUCGAGUAGAUCAAA222222222222222222111111GGGAAAACCC111111ACAAAGAAACAACAACAACAAC"
        )
        assert (
            d_seq_struct.structure
            == "....((((.....))))...222222222222222222111111(((....)))111111......................"
        )
        final_seq_struct = sd.apply(d_seq_struct)
        assert sd.steps[0].set.num_used() == 0
        sd.accept_design()
        assert sd.steps[0].set.num_used() == 1

    def test_split(self):
        build_str = "P5-HPBARCODE-HBARCODE6A-SOI-HBARCODE6B-AC-P3"
        params = TestResources.get_complex_params()
        sd = get_seq_struct_designer(10, build_str, params)
        sds = sd.split(2)
        seq_struct = SequenceStructure("GGGAAAACCC", "(((....)))")
        d_seq_struct = sds[0].get_designable_seq_struct(seq_struct)
        assert (
            d_seq_struct.sequence
            == "GGAAGAUCGAGUAGAUCAAA222222222222222222111111GGGAAAACCC111111ACAAAGAAACAACAACAACAAC"
        )
        assert (
            d_seq_struct.structure
            == "....((((.....))))...222222222222222222111111(((....)))111111......................"
        )

    def test_get_solution(self):
        build_str = "P5-HPBARCODE-HBARCODE6A-SOI-HBARCODE6B-AC-P3"
        params = TestResources.get_complex_params()
        sd = get_seq_struct_designer(10, build_str, params)
        seq_struct = SequenceStructure("GGGAAAACCC", "(((....)))")
        d_seq_struct = sd.get_designable_seq_struct(seq_struct)
        sd.apply(d_seq_struct)
        solution = sd.get_solution()
        sd.accept_previous_solution(solution)
        assert sd.steps[0].set.num_used() == 1


def _test_designer():
    build_str = "P5-HPBARCODE-HBARCODE6A-SOI-HBARCODE6B-AC-P3"
    params = TestResources.get_complex_params()
    df_sequences = TestResources.get_simple_sequence_df()
    sd = get_seq_struct_designer(len(df_sequences), build_str, params)
    designer = Designer()
    results = designer.design(df_sequences, sd)
    df_results = results.df_results
    row = df_results.iloc[0]
    assert (
        row["design_sequence"]
        == "GGAAGAUCGAGUAGAUCAAA222222222222222222111111GGGAAAACCC111111ACAAAGAAACAACAACAACAAC"
    )


def _test_design():
    build_str = "P5-HPBARCODE-HBARCODE6A-SOI-HBARCODE6B-AC-P3"
    params = TestResources.get_complex_params()
    df_sequences = TestResources.get_simple_sequence_df()
    results = design(1, df_sequences, build_str, params, DesignOpts())
    assert len(results.df_results) == 1


def test_design_w_multithreading():
    build_str = "P5-HPBARCODE-HBARCODE6A-SOI-HBARCODE6B-AC-P3"
    params = TestResources.get_complex_params()
    df_sequences = pd.read_csv(get_test_path() / "resources/libs/C0098.csv")
    results = design(2, df_sequences, build_str, params, DesignOpts())
