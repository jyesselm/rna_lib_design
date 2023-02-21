import yaml
import json
import jsonschema
from jsonschema import Draft4Validator, validators
from rna_lib_design.settings import get_py_path

# TODO validate possible combination of key for segments parameters

def extend_with_default(validator_class):
    validate_properties = validator_class.VALIDATORS["properties"]

    def set_defaults(validator, properties, instance, schema):
        for property_, subschema in properties.items():
            if "default" in subschema and not isinstance(instance, list):
                instance.setdefault(property_, subschema["default"])

        for error in validate_properties(
            validator,
            properties,
            instance,
            schema,
        ):
            yield error

    return validators.extend(
        validator_class,
        {"properties": set_defaults},
    )


def validate_parameters(params, schema):
    """
    Validate the parameters against a schema
    :params: parameters in dictionary form
    """
    # Validate the params against the schema
    FillDefaultValidatingDraft4Validator = extend_with_default(Draft4Validator)
    try:
        FillDefaultValidatingDraft4Validator(schema).validate(params)
    except jsonschema.exceptions.ValidationError as e:
        raise ValueError(e.message)


def parse_parameters_from_file(param_file, schema_file):
    """
    Parse a YAML file and validate from a schema file loaded from json
    """
    # load param_file and validate against schema
    with open(param_file) as f:
        params = yaml.safe_load(f)
    if params is None:
        params = {}
    with open(schema_file) as f:
        schema = json.load(f)
    validate_parameters(params, schema)
    return params
