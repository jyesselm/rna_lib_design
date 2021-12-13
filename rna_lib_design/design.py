import vienna
from rna_lib_design import structure_dict

def get_solution(sd_dicts, struct, cutoff=5.0, max_count=1000, keep_best=10):
    count = 0
    best_score = 1000
    best = None
    best_uses = []
    while 1:
        count += 1
        final_struct = structure_dict.apply(sd_dicts, struct)
        vr = vienna.fold(final_struct.sequence)
        if count > max_count:
            return None, None
        #score = final_struct.dot_bracket_difference(structure.DotBracket(vr.dot_bracket))
        #print(final_struct.dot_bracket, vr.ensemble_diversity, score)
        if final_struct.dot_bracket != vr.dot_bracket:
            continue
        if cutoff < vr.ensemble_diversity:
            continue
        if vr.ensemble_diversity < best_score:
            best_score = vr.ensemble_diversity
            best = final_struct
            best_uses = [x.last for x in sd_dicts]
        if count > keep_best:
            break
    for i, bu in enumerate(best_uses):
        sd_dicts[i].set_used(bu)
    return best, best_score


def write_results_to_csv(constructs, fname='final'):
    f_rna = open(fname + '_rna.csv', "w")
    f_rna.write("name,sequence,structure,ens_defect\n")
    f_dna = open(fname + "_dna.csv", "w")
    f_dna.write("name,sequence\n")
    for c in constructs:
        f_rna.write(
                c[0] + "," + c[1].sequence.str() + "," + str(c[1].dot_bracket) + "," +
                str(c[2]) + "\n")
        f_dna.write(c[0] + "," + 'TTCTAATACGACTCACTATA'+c[1].sequence.to_dna().str() + "\n")
    f_rna.close()
    f_dna.close()