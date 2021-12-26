import click
import os

import pandas as pd
import vienna
from rna_lib_design.util import max_gc_stretch, max_stretch, hamming, random_helix
from rna_lib_design import logger

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
                spl[0] + "," + spl[1] + "," + spl2[0] + "," + spl2[1] + "," + str(r.mfe) + "\n"
        )
    f.close()


@click.group()
def cli():
    pass


@cli.command()
@click.option('-l', '--length', type=int, required=True)
@click.option('-md', '--min-dist', type=int, required=True)
@click.option('-gu', '--gus', default=0, type=int)
@click.option('-o', '--output', default="helices.csv")
def hcodes(length, min_dist, gus, output):
    barcodes = generate_helix_barcodes(length, min_dist, gus)
    log.info(f"{len(barcodes)} barcodes found!")
    write_barcodes_to_file(output, barcodes)


@cli.command()
@click.option('-lmin', '--length-min', type=int, required=True)
@click.option('-lmax', '--length-max', type=int, required=True)
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
                    log.info(f"hlen={pos}\tmin_dist={md}\tgu={gu}\t{len(barcodes)} barcodes found!")
                    if not os.path.isdir(f"helices/len_{pos}"):
                        os.mkdir(f"helices/len_{pos}")
                    write_barcodes_to_file(f"helices/len_{pos}/md_{md}_gu_{gu}_{i}.csv", barcodes)
                    df.loc[df_pos] = [
                        pos, md, gu, len(barcodes), f"helices/len_{pos}/md_{md}_gu_{gu}_{i}.csv"]
                    df_pos += 1
    df.to_csv('helices.csv', index=False)


if __name__ == "__main__":
    cli()
