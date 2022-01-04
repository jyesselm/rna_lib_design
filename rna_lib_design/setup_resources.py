import click
import os
import itertools

import pandas as pd
import vienna
from rna_lib_design.util import (
    max_gc_stretch,
    max_stretch,
    hamming,
    random_helix,
)
from rna_lib_design import logger, structure_set, params, structure

log = logger.setup_applevel_logger()


def generate_helix_barcodes(length, min_distance, gus):
    barcodes = []
    count = 0
    while True:
        count += 1
        if count > 1000:
            break
        seq_1, seq_2 = random_helix(length, gu=gus)
        barcode = seq_1 + "&" + seq_2
        if barcode in barcodes:
            continue
        if max_stretch(seq_1) > 3 or max_stretch(seq_2) > 3:
            continue
        if max_gc_stretch(seq_1, seq_2) > 3:
            continue
        fail = 0
        for b in barcodes:
            dist = hamming(b, barcode)
            if dist < min_distance:
                fail = 1
                break
        if fail:
            continue
        if len(barcodes) % 1000 == 0 and len(barcodes) != 0:
            log.info(f"{len(barcodes)} barcodes found so far")
        barcodes.append(barcode)
        count = 0
        if len(barcodes) > 200000:
            log.warn("reached max num of barcodes: 200000")
            break

    return barcodes


def write_barcodes_to_file(fname, barcodes):
    f = open(fname, "w")
    f.write("seq_1,seq_2,ss_1,ss_2,dg\n")
    for b in barcodes:
        spl = b.split("&")
        seq = spl[0] + "CUCUUCGGAG" + spl[1]
        r = vienna.fold(seq)
        spl2 = str(r.dot_bracket).split("(((....)))")
        if len(spl2) != 2:
            continue
        f.write(
            spl[0]
            + ","
            + spl[1]
            + ","
            + spl2[0]
            + ","
            + spl2[1]
            + ","
            + str(r.mfe)
            + "\n"
        )
    f.close()


@click.group()
def cli():
    pass


@cli.command()
@click.option("-l", "--length", type=int, required=True)
@click.option("-md", "--min-dist", type=int, required=True)
@click.option("-gu", "--gus", default=0, type=int)
@click.option("-o", "--output", default="helices.csv")
def hcodes(length, min_dist, gus, output):
    barcodes = generate_helix_barcodes(length, min_dist, gus)
    log.info(f"{len(barcodes)} barcodes found!")
    write_barcodes_to_file(output, barcodes)


@cli.command()
@click.option("-lmin", "--length-min", type=int, required=True)
@click.option("-lmax", "--length-max", type=int, required=True)
def hcodesweep(length_min, length_max):
    os.makedirs("helices", exist_ok=True)
    df = pd.DataFrame(columns="length diff gu size path".split())
    df_pos = 0
    for pos in range(length_min, length_max + 1):
        for md in range(pos - 1, pos * 2, 1):
            for gu in range(0, int(pos / 2), 1):
                for i in range(0, 1):
                    barcodes = generate_helix_barcodes(pos, md, gu)
                    if len(barcodes) < 10:
                        continue
                    log.info(
                        f"hlen={pos}\tmin_dist={md}\tgu={gu}\t{len(barcodes)} barcodes found!"
                    )
                    if not os.path.isdir(f"helices/len_{pos}"):
                        os.mkdir(f"helices/len_{pos}")
                    write_barcodes_to_file(
                        f"helices/len_{pos}/md_{md}_gu_{gu}_{i}.csv", barcodes
                    )
                    df.loc[df_pos] = [
                        pos,
                        md,
                        gu,
                        len(barcodes),
                        f"helices/len_{pos}/md_{md}_gu_{gu}_{i}.csv",
                    ]
                    df_pos += 1
    df.to_csv("helices.csv", index=False)


def __check_folds(sx_struct, sy_struct, h_set, hp_set, n=100, just_list=False):
    total = 0
    for i in range(n):
        h1, h2 = h_set.get(i)
        hp = hp_set.get(i)[0]
        struct = h1 + sx_struct + hp + sy_struct + h2
        vr = vienna.fold(struct.sequence)
        if struct.dot_bracket == vr.dot_bracket:
            total += 1
    return total


@cli.command()
@click.argument("sx", type=int)
@click.argument("sy", type=int)
@click.option("--list", "just_list")
def twoways(sx, sy, just_list):
    h_set = structure_set.get_optimal_helix_set(6, 100)
    loop_struct = structure.get_common_struct("uucg_loop")
    hp_set = structure_set.get_optimal_hairpin_set(6, loop_struct, 100)
    hp_set.remove_dots()
    hp_set.set_buffer(None)
    bases = "ACGU"
    basepairs = params.basepairs
    data = []
    for bp1 in basepairs:
        for bp2 in basepairs:
            print(bp1, bp2)
            if sx > 0:
                sx_combos = itertools.product(*[bases for _ in range(sx)])
            else:
                sx_combos = [""]
            for c1 in sx_combos:
                sx_seq_str = bp1[0] + "".join(c1) + bp2[0]
                sx_ss_str = "(" + "." * len(c1) + "("
                sx_struct = structure.rna_structure(sx_seq_str, sx_ss_str)
                if sy > 0:
                    sy_combos = itertools.product(*[bases for _ in range(sy)])
                else:
                    sy_combos = [""]
                for c2 in sy_combos:
                    sy_seq_str = bp2[1] + "".join(c2) + bp1[1]
                    sy_ss_str = ")" + "." * len(c2) + ")"
                    sy_struct = structure.rna_structure(sy_seq_str, sy_ss_str)
                    score = __check_folds(
                        sx_struct, sy_struct, h_set, hp_set, 10, just_list
                    )
                    if score == 0:
                        continue
                    data.append([sx_seq_str, sy_seq_str, sx_ss_str, sy_ss_str, score])
    df = pd.DataFrame(data, columns="seq_1,seq_2,ss_1,ss_2,score".split(","))
    print(df)
    df.to_csv("out.csv")


if __name__ == "__main__":
    cli()
