import random

from seq_tools import vienna

def generate_barcodes(length, min_distance):
    barcodes = []
    count = 0
    while True:
        count += 1
        if count > 1000:
            break
        seq_1, seq_2 = random_helix(length)
        barcode = seq_1 + "&" + seq_2
        if barcode in barcodes:
            continue
        if max_stretch(seq_1) > 2 or max_stretch(seq_2) > 2:
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
        if len(barcodes) > 100000:
            print("reached max num of barcodes: 100000")
            break

    return barcodes


def main():
    for i in range(9, 15):
        for j in range(6, 12):
            if i * 2 < j:
                continue
            for k in range(1):
                barcodes = generate_barcodes(i, j)
                f = open(
                    "helix_barcode_length_{}_min_dist_{}_{}.csv".format(i, j, k),
                    "w",
                )
                f.write("seq_1,seq_2,ss_1,ss_2\n")
                print(i, j, k, len(barcodes))
                for b in barcodes:
                    spl = b.split("&")
                    seq = spl[0] + "CUUCGG" + spl[1]
                    db = vienna.fold(seq).dot_bracket
                    spl2 = str(db).split("(....)")
                    if len(spl2) != 2:
                        continue
                    f.write(
                        spl[0] + "," + spl[1] + "," + spl2[0] + "," + spl2[1] + "\n"
                    )
                f.close()


if __name__ == "__main__":
    main()
