import unittest
import pandas as pd

from seq_tools import barcoder, settings, structure, vienna


def get_test_dataframe(path, prime5, prime3):
    df = pd.read_csv(path)
    if "structure" not in df:
        structures = []
        for i, row in df.iterrows():
            struct = vienna.fold("CG" + row["sequence"] + "CG").dot_bracket
            structures.append(str(struct)[2:-2])
            print(structures[0])
            exit()
        df["structure"] = structures

    df["prime5_sequence"] = str(prime5.sequence)
    df["prime5_structure"] = str(prime5.dot_bracket)
    df["prime3_sequence"] = str(prime3.sequence)
    df["prime3_structure"] = str(prime3.dot_bracket)
    return df


class BarcoderUnittest(unittest.TestCase):
    def test(self):
        structures = structure.common_structures()
        path = settings.UNITTEST_PATH + "/resources/opool_final_2.csv"
        df = get_test_dataframe(
            path, structures["ref_hairpin_5prime"], structures["rt_tail"]
        )
        bcoder = barcoder.HelixBarcoder()
        df_new = bcoder.barcode(df, 3)
        print(df_new)


def main():
    unittest.main()


if __name__ == "__main__":
    main()
