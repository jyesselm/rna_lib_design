import re
import math
import pickle
import pandas as pd
from enum import Enum
from abc import ABC, abstractmethod
from multipledispatch import dispatch

ALLOWED_PAIRS = {"AU", "UA", "GC", "CG", "UG", "GU"}
HAS_RNA_DRAW = False


class MotifType(Enum):
    UNASSIGNED = 0
    SINGLESTRAND = 1
    HELIX = 2
    HAIRPIN = 3
    JUNCTION = 4


TYPE_MAPPER = {
    MotifType.SINGLESTRAND: "SingleStrand",
    MotifType.HELIX: "Helix",
    MotifType.HAIRPIN: "Hairpin",
    MotifType.JUNCTION: "Junction",
}


def connectivity_list(structure):
    connections, pairs = [-1] * len(structure), []

    for index, db in enumerate(structure):
        if db == "(":
            pairs.append(index)
        elif db == ")":
            complement = pairs.pop()
            connections[complement] = index
            connections[index] = complement

    assert len(pairs) == 0

    return connections


def is_circular(start, connections):
    it = start + 1
    while True:
        while it < len(connections) and connections[it] == -1:
            it += 1
        if it == len(connections):
            return False

        it = connections[it] + 1
        if it == start or it < start:
            return True


class Motif(ABC):
    def __init__(self, **kwargs):
        self.type_ = MotifType.UNASSIGNED
        self.parent_ = None
        self.sequence_ = str()
        self.children_ = []
        self.strands_ = []
        self.depth_ = int()
        self.structure_ = str()
        self.token_ = str()
        self.start_pos_ = math.inf
        self.end_pos_ = -1
        self.positions_ = set()
        self.id_ = None

        if "sequence" in kwargs:
            self.sequence_ = kwargs["sequence"]

        if "strands" in kwargs:
            self.strands_ = kwargs["strands"]

        if len(self.sequence_) == 0:
            return

        strand_seqs = []
        for strand in self.strands_:
            self.start_pos_ = min(min(strand), self.start_pos_)
            self.end_pos_ = max(max(strand), self.end_pos_)

            for pos in strand:
                self.positions_.add(pos)

            strand_seqs.append("".join([self.sequence_[index] for index in strand]))

        self.sequence_ = "&".join(strand_seqs)

    def link_children_(self, depth):

        if self.type() == MotifType.SINGLESTRAND:
            depth = 0

        self.depth(depth)

        # there might be empty lists in here... remove them if that is the case
        self.children_ = [
            child for child in self.children() if not isinstance(child, list)
        ]

        for child in self.children():
            child.parent(self)

        for child in self.children():
            child.link_children_(depth + 1)

    def str(self):
        if self.id_ is not None:
            identification = f"ID: {self.id_}, "
        else:
            identification = ""

        depth = "\t" * self.depth()
        if not self.has_children():
            return f"{depth}{identification}{self.token_} {self.structure_} {self.sequence_}"
        else:
            contents = [""]
            for index, child in enumerate(self.children()):
                contents.append(child.str())

                if (
                    index < len(self.children()) - 1
                    and self.children()[index + 1].is_singlestrand()
                ):
                    tks = contents[-1].splitlines()[0]
                    first_line = contents[-1].splitlines()[0]
                    length = len(first_line) - len(first_line.lstrip())
                    # contents.append('<-' + '-'*len(contents[-1].splitlines()[0]) )
            # children = '\n'.join([""] + [child.str() for child in self.children()])
            children = "\n".join(contents)
            return f"{depth}{identification}{self.token_} {self.structure_} {self.sequence_}{children}"

    def __eq__(self, other) -> bool:
        return (
            self.type_ == other.type_
            and self.sequence_ == other.sequence_
            and self.token_ == other.token_
        )

    def __str__(self) -> str:
        return TYPE_MAPPER[self.type_] + "," + self.sequence_ + "," + self.structure_

    # def __str__(self):
    #    return self.str()

    # def __repr__(self):
    #    return self.str()

    def is_helix(self):
        return False

    def is_singlestrand(self):
        return False

    def is_hairpin(self):
        return False

    def is_junction(self):
        return False

    def type(self):
        return self.type_

    def children(self):
        return self.children_

    def set_children(self, other):
        self.children_ = other

    def parent(self, other=None):
        if other is not None:
            self.parent_ = other
        else:
            return self.parent_

    def token(self, tk=None):
        if tk is None:
            return self.token_
        else:
            self.token_ = tk

    def structure(self, secstruct=None):
        if secstruct is None:
            return self.structure_
        else:
            self.structure_ = secstruct

    def strands(self):
        return self.strands_

    def sequence(self, seq=None):
        if seq is None:
            return self.sequence_
        else:
            self.sequence_ = seq

    @dispatch()
    def id(self):
        return self.id_

    @dispatch(int)
    def id(self, new_id):
        self.id_ = new_id

    @dispatch()
    def depth(self):
        return self.depth_

    @dispatch(int)
    def depth(self, value):
        self.depth_ = value

    @abstractmethod
    def buffer(self):
        pass

    def has_children(self):
        return len(self.children_) > 0

    def has_parent(self):
        return self.parent_ is not None

    @abstractmethod
    def recursive_sequence(self):
        pass

    @abstractmethod
    def recursive_structure(self):
        pass

    @abstractmethod
    def has_non_canonical(self):
        pass

    def start_pos(self):
        return self.start_pos_

    def end_pos(self):
        return self.end_pos_

    def contains(self, pos):
        return pos in self.positions_


class SingleStrand(Motif):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.type_ = MotifType.SINGLESTRAND

        self.structure_ = "." * len(self.sequence_)
        self.token_ = f"SingleStrand{len(self.sequence_)}"

    def buffer(self):
        return -1

    def is_singlestrand(self):
        return True

    def recursive_structure(self):
        result = self.structure_
        for child in self.children():
            result += child.recursive_structure()
        return result

    def recursive_sequence(self):
        result = self.sequence_
        for child in self.children():
            result += child.recursive_sequence()
        return result

    def has_non_canonical(self):
        return False


class Helix(Motif):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.type_ = MotifType.HELIX
        self.size_ = (len(self.sequence_) - 1) // 2
        self.structure_ = f'{"("*self.size_}&{")"*self.size_}'
        self.token_ = f"Helix{self.size_}"
        self.pairs_ = []

        for index in range(self.size_):
            self.pairs_.append(self.sequence_[index] + self.sequence_[-index - 1])

    def size(self, val=None):
        if val is None:
            return self.size_
        else:
            self.size_ = val

    def buffer(self):
        return self.size_

    def pairs(self):
        return self.pairs_

    def is_helix(self):
        return True

    def recursive_structure(self):
        result = self.structure_.split("&")

        result = result[0] + self.children()[0].recursive_structure() + result[1]
        if len(self.children()) == 2:
            result += self.children()[1].recursive_structure()

        return result

    def recursive_sequence(self):
        result = self.sequence_.split("&")

        result = result[0] + self.children()[0].recursive_sequence() + result[1]
        if len(self.children()) == 2:
            result += self.children()[1].recursive_sequence()

        return result

    def has_non_canonical(self):
        for pair in self.pairs():
            if pair not in ALLOWED_PAIRS:
                return True
        return False


class Hairpin(Motif):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.type_ = MotifType.HAIRPIN

        size = len(self.sequence_) - 2
        self.structure_ = "(" + "." * size + ")"
        self.token_ = f"Hairpin{size}"

    def buffer(self):
        return self.parent().buffer()

    def is_hairpin(self):
        return True

    def recursive_structure(self):
        return self.structure_[1:-1]

    def recursive_sequence(self):
        return self.sequence_[1:-1]

    def has_non_canonical(self):
        pair = self.sequence_[0] + self.sequence_[-1]
        return pair not in ALLOWED_PAIRS


class Junction(Motif):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.type_ = MotifType.JUNCTION
        self.closing_pairs_ = []
        self.gaps_ = [len(strand) - 2 for strand in self.strands_]

        self.structure_ = ["."] * len(self.sequence_)
        self.closing_pairs_.append(self.sequence_[0] + self.sequence_[-1])

        for index, nt in enumerate(self.sequence_):
            if nt == "&":
                self.structure_[index - 1] = "("
                self.structure_[index] = "&"
                self.structure_[index + 1] = ")"
                self.closing_pairs_.append(
                    self.sequence_[index - 1] + self.sequence_[index + 1]
                )

        self.structure_[0] = "("
        self.structure_[-1] = ")"

        self.structure_ = "".join(self.structure_)
        self.token_ = f"Junction{len(self.strands())}_" + "|".join(
            [str(len(strand) - 2) for strand in self.strands_]
        )
        self.num_branches_ = len(self.strands())
        self.symmetric_ = len(set([len(strand) - 2 for strand in self.strands_])) == 1

    def buffer(self):
        buffers = [self.parent.buffer()]
        for child in self.children():
            buffers.append(child.buffer())
        return buffers

    def gaps(self):
        return self.gaps_

    def is_junction(self):
        return True

    def recursive_structure(self):
        result = str()
        secstruct_chunks = [subseq[1:-1] for subseq in self.structure_.split("&")]
        for secstruct_chunk, child in zip(secstruct_chunks, self.children()):
            result += secstruct_chunk
            result += child.recursive_structure()
        result += secstruct_chunks[-1]
        return result

    def recursive_sequence(self):
        result = str()
        seq_chunks = [subseq[1:-1] for subseq in self.sequence_.split("&")]
        for seq_chunk, child in zip(seq_chunks, self.children()):
            result += seq_chunk
            result += child.recursive_sequence()
        result += seq_chunks[-1]
        return result

    def closing_pairs(self):
        return self.closing_pairs_

    def has_non_canonical(self):
        for pair in self.closing_pairs_:
            if pair not in ALLOWED_PAIRS:
                return True
        return False

    def number_branches(self):
        return self.num_branches_

    def symmetric(self):
        return self.symmetric_


# wrap these into a class or make private with __func_name
def helix_length_(connections, start):
    complement = connections[start]
    length = 0
    while (
        connections[start + length] == complement - length
        and connections[complement - length] == start + length
    ):
        length += 1
    return length


def get_junction_or_hairpin_(sequence, connections, start):

    strands = []
    offset = 1
    it = start

    while True:
        next = [it]
        it += 1

        while connections[it] == -1:
            next.append(it)
            it += 1

        next.append(it)
        strands.append(next)
        it = connections[it]

        if it == start:
            break
    if len(strands) > 1:
        # JUNCTION
        junction = Junction(strands=strands, sequence=sequence)
        for strand in strands[:-1]:
            junction.children_.append(get_motifs_(sequence, connections, strand[-1]))
        return junction
    else:
        hp = Hairpin(sequence=sequence, strands=strands)
        return hp


def get_helix_(sequence, connections, start):
    # figure out how big this current helix i
    helix_len = helix_length_(connections, start)
    lhs, rhs = [], []
    for index in range(start, start + helix_len):
        lhs.append(index)
        rhs.append(connections[index])
    rhs.reverse()
    helix = Helix(sequence=sequence, strands=[lhs, rhs])
    # figure out what is at the other end of the junction
    if connections[start + helix_len - 1] > start:
        helix.children_.append(
            get_junction_or_hairpin_(sequence, connections, start + helix_len - 1)
        )

    if not is_circular(rhs[-1], connections):
        helix.children_.append(get_motifs_(sequence, connections, rhs[-1] + 1))
    return helix


def get_singlestrand_(sequence, connections, start):
    offset = 0
    while start + offset < len(connections) and connections[start + offset] == -1:
        offset += 1

    strand = list(range(start, start + offset))
    singlestrand = SingleStrand(sequence=sequence, strands=[strand])
    # checking here if there is an outer nt
    if start + offset < len(connections):
        singlestrand.children_.append(
            get_motifs_(sequence, connections, start + offset)
        )

    return singlestrand


def get_motifs_(sequence, connections, start):
    if start >= len(connections):
        return []
    elif connections[start] < 0:
        return get_singlestrand_(sequence, connections, start)
    return get_helix_(sequence, connections, start)


def parse_to_motifs(structure, sequence):
    # basic sanity checks
    assert len(structure) == len(sequence)
    assert len(re.sub("[\(\.\)]", "", structure)) == 0
    assert len(re.sub("[ACGUT]", "", sequence)) == 0
    assert structure.count("(") == structure.count(")")
    for ii in range(3):
        invalid = "(" + "." * ii + ")"
        assert structure.find(invalid) == -1
    # generate the connections
    connections = connectivity_list(structure)
    # get the motifs
    root = get_motifs_(sequence, connections, 0)
    # reverse link them to the parents
    root.link_children_(0)
    return root


def traverse(motif):
    assert not motif.has_non_canonical()
    for c in motif.children():
        # print('chilren')
        # print('\t',c)
        traverse(c)


class SecStruct:
    # what do we want this to do?
    # 1. serve as a place to hold the motifs
    # 2. allow for insertion and changing of motifs

    # dont use named variables when there was only 2 variables
    # switch sequence and secstruct
    def __init__(self, secstruct: str, sequence: str):
        self.root_ = None
        self.structure_ = ""
        self.sequence_ = ""
        self.id_mapping_ = {}
        self.total_structures_ = 1

        assert len(secstruct) == len(sequence) and len(secstruct) > 0
        assert len(re.sub("[\(\.\)]", "", secstruct)) == 0
        assert len(re.sub("[ACGUT]", "", sequence)) == 0
        assert secstruct.count("(") == secstruct.count(")")
        self.has_rna_draw_ = HAS_RNA_DRAW

        self.structure_ = secstruct
        self.sequence_ = sequence

        self.root_ = parse_to_motifs(secstruct, sequence)
        # helpe for recusion
        self.counter_ = 0
        self.set_ids_(self.root_)

    def set_ids_(self, m: Motif):
        m.id(self.counter_)
        self.id_mapping_[self.counter_] = m

        self.counter_ += 1

        for child in m.children():
            self.set_ids_(child)

    def display(self):
        if HAS_RNA_DRAW:
            rna_draw.rna_draw(
                ss=self.structure_,
                seq=self.sequence_,
                out=f"structure{self.total_structures_}.png",
            )
        print(self.root_.str())

    @property
    def sequence(self):
        return self.sequence_

    @property
    def structure(self):
        return self.structure_

    def helix_replace_(self, id, secstruct, sequence):
        helix = self.id_mapping_[id]
        helix.structure(secstruct)
        helix.sequence(sequence)

    def motif_replace_(self, id, new_secstruct, new_sequence):
        motif = self.id_mapping_[id]
        new_motif = parse_to_motifs(new_secstruct, new_sequence)

        if motif.has_parent():
            # reset the parent of the new motif
            new_motif.parent(motif.parent())
            motif.parent(None)

            cleaned_children = []
            for child in motif.parent().children():
                if child.id() == id:
                    cleaned_children.append(new_motif)
                else:
                    cleaned_children.append(child)

            new_motif.parent().set_children(cleaned_children)

        else:
            # this is the case where you are actually just replacing the root? seems like a dumb thing to do
            # but who knows what people will try these days
            pass

    def change_motif(self, id, new_secstruct, new_sequence):
        # performing some basic checks
        assert len(new_sequence) == len(new_secstruct)
        assert len(re.sub("[AUCGT&]", "", new_sequence)) == 0
        assert len(re.sub("[.)(&]", "", new_secstruct)) == 0
        assert (
            new_secstruct.count("&") == new_sequence.count("&")
            and new_sequence.count("&") <= 1
        )
        assert id in self.id_mapping_

        helix_replace = new_secstruct.find("&") != -1

        if helix_replace:
            self.helix_replace_(id, new_secstruct, new_sequence)
        else:
            self.motif_replace_(id, new_secstruct, new_sequence)

        full_secstruct, full_sequence = (
            self.root_.recursive_structure(),
            self.root_.recursive_sequence(),
        )

        self.structure_ = full_secstruct
        self.sequence_ = full_sequence

        self.root_ = parse_to_motifs(full_secstruct, full_sequence)
        # helpe for recusion
        self.counter_ = 0
        self.id_mapping_ = dict()
        self.set_ids_(self.root_)
        self.total_structures_ += 1

    def get_sequence_structure(self):
        return self.sequence_, self.structure_

    def _get_ids_internal(self, m, ids, mtype):
        if m.type() == mtype:
            ids.append((m.id(), m.token(), m.sequence()))

        for c in m.children():
            self._get_ids_internal(c, ids, mtype)

    def get_ids(self, motif_type):
        ids = []
        for mtype, string in TYPE_MAPPER.items():
            if string == motif_type:
                self._get_ids_internal(self.root_, ids, mtype)
                return ids
        else:
            print(
                f"{motif_type} is an invalid motif type. Allowed values are:\n{', '.join(list(TYPE_MAPPER.values()))}"
            )

    def get_motif(self, id) -> Motif:
        motif = self.id_mapping_[id]
        return motif

    def get_substructure(self, id1, id2=None):
        # get substructure between id1 and id2
        # or get all children of id1
        # make
        pass


if __name__ == "__main__":
    ss, seq = (
        ".....((((..((((...........))))..((((...........))))..))))...........(((((((....))))))).....................",
        "GGAAAGUCCAAGAGCGAACAAACAAAGCUCAAGUACGAACCAACAAAGUACAAGGACAAAGAAAGAAAGCUGUCCUUCGGGGCAGCAAAAGAAACAACAACAACAAC",
    )
    SS = SecStruct(secstruct=ss, sequence=seq)
    SS.display()
    SS.change_motif(3, "(((((((&)))))))", "GGGGGGG&CCCCCCC")
    SS.display()
    SS.change_motif(3, "(...)", "GAAAC")
    SS.display()
    SS.change_motif(
        3, "(((..(((....)))..(((....)))..)))", "GGGAAGGGAAAACCCAAGGGAAAACCCAACCC"
    )
    SS.display()
    print(SS.get_ids("Hairpin"))
