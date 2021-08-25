import random

import click
import vienna
from rna_lib_design.util import max_gc_stretch, max_stretch, hamming, random_helix


def generate_barcodes(length, min_distance):
    barcodes = []
    count = 0
    while True:
        count += 1
        if count > 1000:
            break
        seq_1, seq_2 = random_helix(length, gu=1)
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
        barcodes.append(barcode)
        count = 0
        if len(barcodes) > 500000:
            print("reached max num of barcodes: 500000")
            break

    return barcodes


@click.command()
@click.option("--length", required=True, type=int, help="helix length")
@click.option("--min_diff", required=True, type=int, help="min sequence diff")
def main(length, min_diff):
    barcodes = generate_barcodes(length, min_diff)
    f = open("org_helices.csv", "w")
    f.write("seq_1,seq_2\n")
    for b in barcodes:
        spl = b.split("&")
        f.write(spl[0] + "," + spl[1] + "\n")
    f.close()
    f = open("helices.csv", "w")
    f.write("seq_1,seq_2,ss_1,ss_2,dg\n")
    for b in barcodes:
        spl = b.split("&")
        seq = spl[0] + "CUUCGG" + spl[1]
        r = vienna.fold(seq)
        spl2 = str(r.dot_bracket).split("(....)")
        if len(spl2) != 2:
            continue
        f.write(
                spl[0] + "," + spl[1] + "," + spl2[0] + "," + spl2[1] + "," + str(r.mfe) + "\n"
        )
    f.close()


if __name__ == "__main__":
    main()
