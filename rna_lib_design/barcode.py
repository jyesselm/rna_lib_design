from rna_lib_design.logger import get_logger
from rna_lib_design.parameters import parse_parameters_from_file
from rna_lib_design.settings import get_resources_path

log = get_logger("BARCODE")


def get_preset_parameters(btype: str, barcode_name: str):
    """
    get the preset parameters for a given barcode type and name
    """
    full_name = f"{barcode_name}_{btype}.yml"
    if not (get_resources_path() / "presets" / full_name).exists():
        presets = (get_resources_path() / "presets").glob(f"{barcode_name}*.yml")
        log.error(presets)
        raise ValueError(f"Invalid barcode type: {btype}")
    schema_file = get_resources_path() / "schemas" / f"{barcode_name}.json"
    params_file = get_resources_path() / "presets" / full_name
    params = parse_parameters_from_file(params_file, schema_file)
    return params


def single_barcode(df, params):
    pass
