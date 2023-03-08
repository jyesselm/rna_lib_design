import os
import pandas as pd
from rna_lib_design.settings import get_resources_path


def main():
    """
    main function for script
    """
    input_dir = get_resources_path() / "barcodes"
    output_dir = "barcodes"
    os.makedirs(output_dir, exist_ok=True)
    for root, dirs, files in os.walk(get_resources_path()):
        for file in files:
            if not file.endswith(".csv"):
                continue
            input_path = os.path.join(root, file)
            rel_path = os.path.relpath(input_path, input_dir)
            output_path = os.path.join(output_dir, rel_path)
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            df = pd.read_csv(input_path)
            if "seq_1" in df:
                df["sequence"] = df["seq_1"] + "&" + df["seq_2"]
                df["structure"] = df["ss_1"] + "&" + df["ss_2"]
                df.drop(
                    columns=["seq_1", "seq_2", "ss_1", "ss_2"], inplace=True
                )
                if "dg" in df:
                    df = df[["sequence", "structure", "dg"]]
            elif "seq" in df:
                df.rename(
                    columns={"seq": "sequence", "ss": "structure"}, inplace=True
                )
            df.to_csv(output_path, index=False)


# pylint: disable=no-value-for-parameter
if __name__ == "__main__":
    main()
