from seq_tools import structure
import collections as col

dna_di = {
    "AA": 27400,
    "AC": 21200,
    "AG": 25000,
    "AT": 22800,
    "CA": 21200,
    "CC": 14600,
    "CG": 18000,
    "CT": 15200,
    "GA": 25200,
    "GC": 17600,
    "GG": 21600,
    "GT": 20000,
    "TA": 23400,
    "TC": 16200,
    "TG": 19000,
    "TT": 16800,
}

dna_mono = {"A": 15400, "C": 7400, "G": 11500, "T": 8700}

rna_di = {
    "AA": 27400,
    "AC": 21200,
    "AG": 25000,
    "AU": 24000,
    "CA": 21200,
    "CC": 14600,
    "CG": 18000,
    "CU": 16200,
    "GA": 25200,
    "GC": 17600,
    "GG": 21600,
    "GU": 21200,
    "UA": 24600,
    "UC": 17200,
    "UG": 20000,
    "UU": 19600,
}

rna_mono = {"A": 15400, "C": 7400, "G": 11500, "U": 9900}

bracket_left = "([{<ABCDEFGHIJKLMNOPQRSTUVWXYZ"
bracket_right = ")]}>abcdefghijklmnopqrstuvwxyz"


def inverse_brackets(bracket):
    res = col.defaultdict(int)
    for i, a in enumerate(bracket):
        res[a] = i
    return res


def dotbracket_to_pairtable(struct):
    """
    Converts arbitrary structure in dot bracket format to pair table
    (ViennaRNA format).
    """
    if len(struct) == 0:
        raise ValueError("Cannot convert empty structure to pairtable")
    pt = [0] * ((len(struct)) - struct.count("&"))
    # pt[0] = len(struct) - struct.count("&")

    stack = col.defaultdict(list)
    inverse_bracket_left = inverse_brackets(bracket_left)
    inverse_bracket_right = inverse_brackets(bracket_right)

    i = 0
    for a in struct:
        if a == "&":
            continue
        i += 1
        # print i,a, pt
        if a == ".":
            pt[i - 1] = -1
        else:
            if a in inverse_bracket_left:
                stack[inverse_bracket_left[a]].append(i)
            else:
                assert a in inverse_bracket_right
                if len(stack[inverse_bracket_right[a]]) == 0:
                    raise ValueError("Too many closing brackets!")
                j = stack[inverse_bracket_right[a]].pop()
                pt[i - 1] = j - 1
                pt[j - 1] = i - 1

    if len(stack[inverse_bracket_left[a]]) != 0:
        raise ValueError("Too many opening brackets!")

    return pt


def reverse_complement(seq, type):
    complement = ""
    for e in seq:
        if e == "A" and type == "RNA":
            complement += "U"
        elif e == "C":
            complement += "G"
        elif e == "G":
            complement += "C"
        elif e == "U":
            complement += "A"
        elif e == "T":
            complement += "A"
    return complement[::-1]


def get_mono_contribution_rna(seq):
    total = 0
    for e in seq[1:-1]:
        total += rna_mono[e]
    return total


def get_mono_contribution_dna(seq):
    total = 0
    for e in seq[1:-1]:
        total += dna_mono[e]
    return total


def get_di_contribution_rna(seq):
    total = 0
    for i in range(0, len(seq) - 1):
        distep = seq[i] + seq[i + 1]
        total += rna_di[distep]
    return total


def get_di_contribution_dna(seq):
    total = 0
    for i in range(0, len(seq) - 1):
        distep = seq[i] + seq[i + 1]
        total += dna_di[distep]
    return total


def get_hypochromicity_dna(seq):
    frac_at = 0
    for e in seq:
        if e == "A" or e == "T":
            frac_at += 1
    frac_at /= len(seq)
    return frac_at * 0.287 + (1 - frac_at) * 0.059


def get_hypochromicity_rna(seq, ss):
    pairtable = dotbracket_to_pairtable(ss)
    frac_au = 0
    frac_gc = 0
    for i in range(0, len(pairtable)):
        if pairtable[i] == -1:
            continue
        pos = pairtable[i]
        pos2 = pairtable[pos]
        name = seq[pos] + seq[pos2]
        if name == "AU" or name == "UA":
            frac_au += 1
        if name == "GC" or name == "CG":
            frac_gc += 1
    frac_au /= len(seq)
    frac_gc /= len(seq)
    return frac_au * 0.26 + frac_gc * 0.059


def get_coefficient_dna(seq, ds):
    mono = get_mono_contribution_dna(seq)
    di = get_di_contribution_dna(seq)
    strand1 = di - mono
    if not ds:
        return strand1
    rc = reverse_complement(seq, 'DNA')
    strand2 = get_di_contribution_dna(rc) - get_mono_contribution_dna(rc)
    hc = get_hypochromicity_dna(seq)
    final = round((1 - hc) * (strand1 + strand2))
    return final


def get_coefficient_rna(seq, ss):
    mono = get_mono_contribution_rna(seq)
    di = get_di_contribution_rna(seq)
    if ss is not None:
        hc = get_hypochromicity_rna(seq, ss)
        return round((1 - hc) * (di - mono))
    else:
        return di - mono
