#from rna_lib_design import util, settings

"""
def test_max_stretch():
    assert util.max_stretch("GGGGC") == 4
    assert util.max_stretch("CGGGG") == 4
    assert util.max_stretch("CGGGC") == 3


def test_indentify_p5():
    seqs = ["GGAACAGCACUUCGGUGCAAAGGGCCCGAGUAGGGUCCAAAGCCUCCAAGGGUUGCUUCGGCA"]
    p5 = util.indentify_p5_sequence(seqs)
    assert p5.name == "uucg_p5_rev_primer"


def test_indentify_p5_2():
    # should not match to anything
    seqs = ["UAUGGAGGCAAAGAAACAACAACAACAAC"]
    p5 = util.indentify_p5_sequence(seqs)
    assert p5.name == ""


def test_indentify_p5_3():
    # should not match to anything no overlap between sequences
    seqs = [
        "GGAACAGCACUUCGGUGCAAAGGGCCCGAGUAGGGUCCAAAGCCUCCAAGGGUUGCUUCGGCA",
        "UAUGGAGGCAAAGAAACAACAACAACAAC",
    ]
    p5 = util.indentify_p5_sequence(seqs)
    assert p5.name == ""


def test_indentify_p5_4():
    # has common sequence
    seqs = [
        "TTCTAATACGACTCACTATAGGAACAGCACUUCGGUGCAAAGGG",
        "TTCTAATACGACTCACTATAGGAACAGCACUUCGGUGCAAACCC",
    ]
    p5 = util.indentify_p5_sequence(seqs)
    assert p5.name == "uucg_p5_rev_primer"


def test_indentify_fwd_primer():
    seqs = ["TTCTAATACGACTCACTATAGGUAUGGAGGCAAAGAAACAACAACAACAAC"]
    fwd_p = util.indentify_fwd_primer(seqs)
    assert fwd_p.name == "t7_fwd_p_pur"


def test_indentify_fwd_primer_2():
    # close but should match to t7_fwd_p_pur
    seqs = [
        "TTCTAATACGACTCACTATAGGAUCGGAGGCAAAGAAACAACAACAACAAC",
        "TTCTAATACGACTCACTATAGGAACGG",
    ]
    fwd_p = util.indentify_fwd_primer(seqs)
    assert fwd_p.name == "t7_fwd_p_pur"


def test_indentify_rev_primer():
    # removed first T should not match now
    seqs = ["TCTAATACGACTCACTATAGGUAUGGAGACAAAGAAACAACAACAACAAC"]
    rev_p = util.indentify_rev_primer(seqs)
    assert rev_p.name == "tail_rev_p_pur"


def test_find_valid_subsequences():
    seqs = [
        "TTCTAATACGACTCACTATAGGAACC",
        "TTCTAATACGACTCACTATAGGAACA",
    ]
    df = util.get_primer_dataframe(settings.RESOURCES_PATH + "fwd_primers.csv")
    df_sub = util.find_valid_subsequences(df, seqs)
    assert len(df_sub) == 2
"""