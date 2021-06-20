import os
import subprocess

from seq_tools import structure


class Results(object):
    def __init__(self, struct, energy, ensemble_prob, ensemble_diversity):
        self.dot_bracket = structure.DotBracket(struct)
        self.mfe = energy
        self.ens_prob = ensemble_prob
        self.ensemble_diversity = ensemble_diversity


class InverseResults(object):
    def __init__(self, seqs, scores):
        self.seqs, self.scores = seqs, scores

    def __len__(self):
        return len(self.seqs)


def fold(seq):
    if len(seq) == 0:
        raise ValueError("must supply a sequence longer then 0")
    output = subprocess.check_output(
            'echo "' + str(seq) + '" | ' + "RNAfold -p --noLP -d2",
            shell=True,
    )
    lines = output.decode("utf-8").split("\n")
    spl1 = lines[1].split()
    spl2 = lines[-2].split()
    ensemble_prob = float(spl2[6][:-1])
    ensemble_diversity = float(spl2[-1])
    structure = spl1[0]
    energy = float(lines[1].split("(")[-1][:-1])
    results = Results(structure, energy, ensemble_prob, ensemble_diversity)
    return results

def folded_structure(seq):
    r = fold(seq)
    return str(r.dot_bracket)

def cofold(seq):
    if len(seq) == 0:
        raise ValueError("must supply a sequence longer then 0")
    os.system(
            'echo "' + seq + '" | ' + self.bin_path + "RNAcofold -p > rnafold_dump"
    )
    f = open("rnafold_dump")
    lines = f.readlines()
    f.close()
    # print lines
    try:
        last_line = lines.pop()
    except:
        return None
    spl = last_line.split()
    ensemble_prob = float(spl[6][:-1])
    ensemble_diversity = float(spl[-1].split("=")[-1])
    spl = lines[1].split()
    spl2 = lines[1].split("(")
    structure = spl[0]
    energy = float(spl2[-1][:-2].rstrip())
    results = FoldingParams(structure, energy, ensemble_prob, ensemble_diversity)
    os.remove("rnafold_dump")
    os.remove("rna.ps")
    os.remove("dot.ps")
    return results


def inverse(ss, constraint, n_sol=100, discard_misfolds=True):
    f = open("seqsecstruct.txt", "w")
    f.writelines([ss + "\n", constraint])
    f.close()

    try:
        output = subprocess.check_output(
                "RNAinverse -Fmp -f 0.5 -d2 -R{} < seqsecstruct.txt".format(n_sol),
                shell=True,
        )
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
                "command '{}' return with error (code {}): {}".format(
                        e.cmd, e.returncode, e.output
                )
        )
    lines = output.decode("utf-8").split("\n")
    seqs = []
    scores = []
    for e in lines:
        print(e)
        spl = e.split()
        if len(spl) != 2:
            continue
        seqs.append(spl[0])
        scores.append(int(spl[1]))
    return ViennaInverseFoldingResults(seqs, scores)
