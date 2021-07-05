from rna_lib_design import structure_dict, settings, structure


def test_sstrand():
    path = settings.RESOURCES_PATH + "/barcodes/sstrand/5_no_gg.csv"
    sd = structure_dict.SStrandDict(path)
    assert len(sd) > 0


def test_sstrand_apply_next():
    sd = structure_dict.SingleDict(structure.rna_structure("AAAAA", "....."), "RIGHT")
    struct1 = structure.rna_structure("GGGAAAACCC", "(((....)))")
    sstrand = sd.get_next()
    struct2 = sd.apply_next(struct1)
    target = struct1 + sstrand
    assert target == struct2


def test_helix_apply_next():
    path = settings.RESOURCES_PATH + "/barcodes/helices/helix_barcode_length_5_min_dist_2.csv"
    sd = structure_dict.HelixDict(path)
    struct1 = structure.rna_structure("GGGAAAACCC", "(((....)))")
    helix = sd.get_next()
    struct2 = helix[0] + struct1 + helix[1]


def test_hairpin_apply_next():
    path = settings.RESOURCES_PATH + "/barcodes/helices/helix_barcode_length_5_min_dist_2.csv"
    helices = structure_dict.HelixDict(path)
    loop = structure.rna_structure('CGAGUAG', '(.....)')
    sd = structure_dict.HairpinDict(helices, loop)
    struct1 = structure.rna_structure("GGGAAAACCC", "(((....)))")
    hp = sd.get_next()
    final = hp + struct1
    assert len(final) == 30



def test_apply():
    path = settings.RESOURCES_PATH + "/barcodes/sstrand/5_no_gg.csv"
    sd1 = structure_dict.SingleDict(structure.rna_structure("GGGGG", "....."))
    sd2 = structure_dict.SingleDict(structure.rna_structure("AAAAA", "....."), "RIGHT")
    sds = [sd1, sd2]
    struct1 = structure.rna_structure("GGGAAAACCC", "(((....)))")
    struct_final = structure_dict.apply(sds, struct1)
    target = structure.rna_structure("GGGGGGGGAAAACCCAAAAA", ".....(((....))).....")
    assert target == struct_final


def test_apply_helix():
    sd1 = structure_dict.SingleDict(structure.rna_structure("GGGGG", "....."))
    left_s = structure.rna_structure("GAGA", "((((")
    right_s = structure.rna_structure("UCUC", "))))")
    sd2 = structure_dict.SingleHelixDict(left_s, right_s)
    sds = [sd1, sd2]
    struct1 = structure.rna_structure("GGGAAAACCC", "(((....)))")
    struct_final = structure_dict.apply(sds, struct1)
    target = structure.rna_structure(
            "GGGGGGAGAGGGAAAACCCUCUC", ".....(((((((....)))))))")
    assert target == struct_final


def test_trim_next():
    left_s = structure.rna_structure("GAGA", "((((")
    right_s = structure.rna_structure("UCUC", "))))")
    sd1 = structure_dict.SingleHelixDict(left_s, right_s)
    sd1.set_trim(2)
    next = sd1.get_next()
    assert next[0].sequence == "GA"
    assert next[1].sequence == "UC"

def test_get_helices():
    helices = structure_dict.get_helices(10)
    assert helices is not None

def test_get_sstrands():
    sstrands = structure_dict.get_sstrands(10)
    assert sstrands is not None
