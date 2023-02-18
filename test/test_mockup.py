import pytest
from rna_lib_design.mockup import parse_sequence_string


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
