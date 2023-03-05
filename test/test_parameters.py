import json
import yaml
import pytest
from rna_lib_design.settings import get_resources_path
from rna_lib_design.parameters import (
    validate_parameters,
    validate_segment_parameters,
)


def test_schema():
    schema_path = get_resources_path() / "schemas" / "single_barcode.json"
    params_path = get_resources_path() / "presets" / "single_barcode_helix.yml"
    schema = json.load(schema_path.open())
    params = yaml.safe_load(params_path.open())
    validate_parameters(params, schema)
    assert params["design_opts"]["max_attempts"] == 10


class TestValidateSegmentParameters:
    def test_mtype(self):
        params = {"m_type": "HELIX", "length": "5-6"}
        try:
            validate_segment_parameters(params)
        except ValueError:
            pytest.fail("Unexpected ValueError")

    def test_mtype_name(self):
        params = {"m_type": "HELIX", "length": "5-6", "name": "test"}
        with pytest.raises(ValueError):
            validate_segment_parameters(params)

    def test_mtype_sequence(self):
        params = {"m_type": "HELIX", "length": "5-6", "sequence": "ACGU"}
        with pytest.raises(ValueError):
            validate_segment_parameters(params)

    def test_mtype_structure(self):
        params = {"m_type": "HELIX", "length": "5-6", "structure": "...."}
        with pytest.raises(ValueError):
            validate_segment_parameters(params)

    def test_sequence(self):
        params = {"sequence": "ACGU", "structure": "...."}
        try:
            validate_segment_parameters(params)
        except ValueError:
            pytest.fail("Unexpected ValueError")
