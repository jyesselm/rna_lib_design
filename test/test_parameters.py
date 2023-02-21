
import json
import yaml
from rna_lib_design.settings import get_resources_path
from rna_lib_design.parameters import validate_parameters

def test_schema():
    schema_path = get_resources_path() / "schemas" / "single_barcode.json"
    params_path = get_resources_path() / "presets" / "single_barcode_helix.yml"
    schema = json.load(schema_path.open())
    params = yaml.safe_load(params_path.open())
    validate_parameters(params, schema)
    assert params["design_opts"]["max_attempts"] == 10