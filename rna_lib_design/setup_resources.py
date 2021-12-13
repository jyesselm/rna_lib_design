import click
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
        if len(barcodes) % 1000 == 0:
            log.info(f"{len(barcodes)} barcodes found so far")
        barcodes.append(barcode)
        count = 0
        if len(barcodes) > 200000:
            log.warn("reached max num of barcodes: 200000")
            break

    return barcodes


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
    f = open(output, "w")
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


if __name__ == "__main__":
    cli()
