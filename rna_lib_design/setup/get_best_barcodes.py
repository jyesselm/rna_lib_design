import pandas as pd
import glob

from seq_tools import settings


def get_best(path):
    files = glob.glob(path)
    best_df = None
    best_size = 0
    best_name = ""
    total_dots = 0
    for f in files:
        df = pd.read_csv(f)
        df["left_dots"] = [x.count(".") for x in df["ss_1"]]
        df["right_dots"] = [x.count(".") for x in df["ss_2"]]
        df["total_dots"] = df["left_dots"] + df["right_dots"]
        df = df[(df["left_dots"] < 3) & (df["right_dots"] < 3)]
        df.drop(columns=["left_dots", "right_dots"], inplace=True)
        if len(df) > best_size:
            best_size = len(df)
            best_df = df
            best_name = f
    if best_df is None:
        return [None, None]
    best_df = best_df.sort_values(["total_dots"])
    return [best_df, best_name]


def main():
    path = settings.RESOURCES_PATH + "barcodes/"
    summary_df = pd.DataFrame(columns="length min_dist count path gu".split())
    pos = 0
    for i in range(3, 20):
        for j in range(2, 12):
            s = "data/with_gu/helix_barcode_length_{}_min_dist_{}_*".format(i, j)
            s_spl = s.split("/")
            df, df_name = get_best(s)
            if df is not None:
                name = "_".join(s_spl[2].split("_")[:-1]) + ".csv"
                print(name)
                summary_df.loc[pos] = [i, j, len(df), path + "/helices/" + name, 1]
                pos += 1
                df.to_csv(path + "/helices/" + name, index=False)

            s = "data/no_gu/helix_barcode_no_gu_length_{}_min_dist_{}_*".format(i, j)
            s_spl = s.split("/")
            df, df_name = get_best(s)
            if df is not None:
                name = "_".join(s_spl[2].split("_")[:-1]) + ".csv"
                summary_df.loc[pos] = [i, j, len(df), path + "/helices/" + name, 0]
                pos += 1
                df.to_csv(path + "/helices/" + name, index=False)

    summary_df.to_csv(path + "helices.csv", index=False)


if __name__ == "__main__":
    main()
