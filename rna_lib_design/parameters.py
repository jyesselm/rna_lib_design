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

    # check to make sure each segment dict is valid
    if "segments" not in params:
        return
    # for segment in params["segments"]:
    #    validate_segment_parameters(segment)


def validate_segment_parameters(params):
    if "m_type" in params:
        if "sequence" in params or "structure" in params:
            raise ValueError(
                "Cannot specify sequence or structure for a segment with m_type"
            )
        elif "name" in params:
            raise ValueError("Cannot specify name for a segment with m_type")
    elif "sequence" in params:
        if "structure" not in params:
            raise ValueError("Must specify structure for a segment with sequence")
        elif "name" in params:
            raise ValueError("Cannot specify name for a segment with sequence")
    elif "name" in params:
        if "structure" in params or "sequence" in params:
            raise ValueError(
                "Cannot specify sequence or structure for a segment with name"
            )
        elif "m_type" in params:
            raise ValueError("Cannot specify m_type for a segment with name")
    else:
        raise ValueError("Must specify either m_type, sequence, or name for a segment")
    if "m_type" in params and params["m_type"] == "HAIRPIN":
        print("made it")


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
