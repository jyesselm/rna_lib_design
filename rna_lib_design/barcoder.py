import pandas as pd

from seq_tools import logger, settings, structure, vienna

log = logger.init_logger("barcoder.py")


def get_barcodes(fname, type):
    df = pd.read_csv(fname)
    barcodes = []
    for i, row in df.iterrows():
        if type.lower() == "helix":
            left_struct = structure.rna_structure(row["seq_1"], row["ss_1"])
            right_struct = structure.rna_structure(row["seq_2"], row["ss_2"])
            barcodes.append(HelixBarcode(left_struct, right_struct))
    return barcodes


class Barcode(object):
    def __init__(self):
        self._used = False
        self._successful = False

    def is_used(self):
        return self._used

    def set_used(self):
        self._used = True

    def get_barcode_sequence(self):
        pass


class HelixBarcode(Barcode):
    def __init__(self, left_struct, right_struct):
        super().__init__()
        self.__left_struct = left_struct
        self.__right_struct = right_struct

    def get_barcode_sequence(self):
        return (
                str(self.__left_struct.sequence) + "&" + str(self.__right_struct.sequence)
        )

    @property
    def left_struct(self):
        return self.__left_struct

    @property
    def right_struct(self):
        return self.__right_struct


class Barcoder(object):
    def __init__(self, fname):
        self.barcodes = get_barcodes(fname)

    def get_barcode(self):
        pass


class HelixBarcoder(object):
    def __init__(self):
        pass

    def __get_barcodes(self, length, required_barcodes, gu):
        path = settings.RESOURCES_PATH + "/barcodes/helices.csv"
        df_barcodes = pd.read_csv(path)
        if length is not None:
            df_barcodes = df_barcodes[df_barcodes["length"] == length]
        if len(df_barcodes) == 0:
            log.error("no barcode options of length: {}".format(length))
            exit(1)
        safe_num_barcodes = required_barcodes * 1.05
        df_barcodes = df_barcodes[df_barcodes["count"] > safe_num_barcodes]
        df_barcodes = df_barcodes[df_barcodes["gu"] == gu]
        if len(df_barcodes) == 0:
            log.error(
                    "no barcode options that can safely barcode {} constructs".format(
                            safe_num_barcodes
                    )
            )
            exit(1)
        df_barcodes = df_barcodes.sort_values(["min_dist"], ascending=False)
        df = pd.read_csv(df_barcodes.iloc[0]["path"])
        log.info(
                "barcoding info: length -> {}, min_dist -> {}, count -> {}, path -> {}, gu -> {}".format(
                        *df_barcodes.iloc[0]
                )
        )
        required_col = "seq_1 seq_2 ss_1 ss_2".split()
        for col in required_col:
            if col not in df:
                log.error(
                        "{} is not in {} which is a required column for barcoding dataframe".format(
                                col, (df_barcodes.iloc[0]["path"])
                        )
                )

        barcodes = []
        for i, row in df.iterrows():
            b1 = structure.rna_structure(row["seq_1"], row["ss_1"])
            b2 = structure.rna_structure(row["seq_2"], row["ss_2"])
            barcodes.append(HelixBarcode(b1, b2))
        return barcodes

    def barcode(self, df: pd.DataFrame, length=None, gu=0) -> pd.DataFrame:
        barcodes = self.__get_barcodes(length, len(df), gu)
        final_sequences, final_structures = [], []
        order_sequences = []
        used_barcodes = []
        for i, row in df.iterrows():
            prime5_struct = structure.rna_structure(
                    row["prime5_sequence"], row["prime5_structure"]
            )
            prime3_struct = structure.rna_structure(
                    row["prime3_sequence"], row["prime3_structure"]
            )
            construct = structure.rna_structure(row["sequence"], row["structure"])
            found = 0
            for b in barcodes:
                if b.is_used():
                    continue
                rna_struct = b.apply(prime5_struct, prime3_struct, construct)
                if not b.is_sucessful():
                    continue
                b.set_used()
                used_barcodes.append(b.get_barcode_sequence())
                final_sequences.append(str(rna_struct.sequence))
                final_structures.append(str(rna_struct.dot_bracket))
                order_sequences.append(
                        "TTCTAATACGACTCACTATA" + str(rna_struct.sequence.to_dna())
                )
                found = 1
                break
            if found:
                continue
            log.error("could not find barcode for {} ".format(construct.sequence))

        df["final_sequence"] = final_sequences
        df["final_structure"] = final_structures
        df["order_sequence"] = order_sequences
        df["barcode"] = used_barcodes
        return df
