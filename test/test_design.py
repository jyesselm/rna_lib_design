import pandas as pd
import os
import pytest
from rna_lib_design import design, structure, structure_set, settings, testing
from click.testing import CliRunner

TEST_RESOURCES = settings.RESOURCES_PATH + "/testing/"


def test_design():
    rna = structure.rna_structure("GGGGAAAACCCC", "((((....))))")
    set_1 = testing.get_test_helices()
    sol = design.get_best_design([set_1], rna, design.DesignOptions())
    assert sol.design.dot_bracket == "((((((((((....))))))))))"
    assert sol.ens_defect < 1.0


def test_design_2():
    rna = structure.rna_structure("GGGGAAAACCCC", "((((....))))")
    set_1 = testing.get_test_helices()
    set_2 = structure_set.get_common_seq_structure_set("ref_hairpin_5prime")
    set_3 = structure_set.get_tail_structure_set()
    sol = design.get_best_design(
        [set_1, set_2, set_3], rna, design.DesignOptions()
    )
    assert (
        sol.design.dot_bracket
        == "....((((.....))))...((((((((((....))))))))))...................."
    )
    assert sol.ens_defect < 1.0
