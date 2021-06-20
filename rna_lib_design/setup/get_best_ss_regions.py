import pandas as pd
import glob

from seq_tools import settings


def main():
    path = "/Users/josephyesselman/projects/sequence_tools/data/common_seqs/singlestrand/"
    r_path = settings.RESOURCES_PATH + "/barcodes/sstrand/"
    files = glob.glob(path + "/*.csv")
    for fname in files:
        df = pd.read_csv(fname)
        df = df[df.ensemble_prob > 0.98]
        df = df[~df.seq.str.contains("AAAA")]
        df = df.sort_values("ensemble_prob", ascending=False)
        spl = fname.split("/")
        f = open(r_path + spl[-1], "w")
        f.write("seq,ss,ensemble_prob\n")
        for i, row in df.iterrows():
            f.write(row["seq"] + "," + "." * len(row["seq"]) + "," + str(row["ensemble_prob"]) + "\n")
        f.close()


if __name__ == "__main__":
    main()
