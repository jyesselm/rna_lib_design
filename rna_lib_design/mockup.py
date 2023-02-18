from seq_tools import SequenceStructure


def generate_helix(length, symbol):
    seq = symbol * length
    return SequenceStructure(seq + "&" + seq, "1" * length + "&" + "1" * length)


def parse_sequence_string(seq_str):
    """
    Parses a string of the form
    "P5-HPBARCODE-HBARCODE6A-SOI-HBARCODE6B-AC-P3" to return a dictionary
    containing the order and direction in which the other segments appear
    in reference to SOI.

    :param seq_str: The sequence string to parse.
    :type seq_str: str
    :return: A dictionary containing the order and direction in
    which the other segments appear in reference to SOI.
    :rtype: dict[str, int]
    """
    segments = seq_str.split("-")
    soi_index = segments.index("SOI")
    order = {}
    for i, segment in enumerate(segments):
        if i != soi_index:
            distance = i - soi_index
            key = segment.rstrip("AB")
            if segment.endswith("A"):
                if key in order:
                    order[key] = min(distance, order[key])
                else:
                    order[key] = distance
            elif segment.endswith("B"):
                if key in order:
                    order[key] = max(distance, order[key])
                else:
                    order[key] = distance
            else:
                if key in order:
                    order[key] = min(distance, order[key])
                else:
                    order[key] = distance
    return order


""" 
def unsued():
    symbols = [
        "1",
        "2",
        "3",
        "4",
        "5",
        "6",
        "7",
        "8",
        "9",
        "0",
        "*",
        "#",
        "$",
        "%",
        "^",
        "*",
        "B",
        "D",
        "F",
    ]
    replaces = []
    pos = 0
    return seqstruct
    for info in params:
        if info.startswith("HELIX-"):
            helix_len = int(info.split("-")[1])
            helix_strands = generate_helix(
                helix_len, symbols[pos]
            ).split_strands()
            seqstruct = helix_strands[0] + seqstruct + helix_strands[1]
    print(seqstruct)

    return replaces, seqstruct
"""
