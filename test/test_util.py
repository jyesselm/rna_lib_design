from rna_lib_design import util


def test_max_stretch():
    assert util.max_stretch("GGGGC") == 4
    assert util.max_stretch("CGGGG") == 4
    assert util.max_stretch("CGGGC") == 3
