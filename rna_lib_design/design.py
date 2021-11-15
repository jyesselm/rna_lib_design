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