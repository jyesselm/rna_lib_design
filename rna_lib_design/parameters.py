import yaml
import json
import jsonschema
from jsonschema import Draft4Validator, validators
from rna_lib_design.settings import get_py_path, get_resources_path
from rna_lib_design.logger import get_logger

log = get_logger("PARAMETERS")


def extend_with_default(validator_class):
    """
    a wrapper to add default values to the schema
    """
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
    :schema: schema in dictionary form
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
    used_segments = params["build_str"].split("-")
    segments = params["segments"].copy()
    for name, segment in segments.items():
        found = 0
        # make sure to catch BARCODE1A and BARCODE1B
        for used_segment in used_segments:
            if used_segment.startswith(name):
                found = 1
                break
        if not found:
            del params["segments"][name]
            continue
        log.debug(name)
        validate_segment_parameters(segment)


def validate_segment_parameters(params):
    log.debug(params)
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
    # if "m_type" in params and params["m_type"] == "HAIRPIN":
    #    print("made it")


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


def combine_params(original_dict, new_dict):
    for key, value in new_dict.items():
        if isinstance(value, dict):
            original_dict[key] = combine_params(original_dict.get(key, {}), value)
        else:
            # safe guard since this might really miss stuff up
            if key == "m_type":
                if original_dict[key] != value:
                    raise ValueError(
                        f"Cannot change m_type from {original_dict[key]} to {value}"
                    )
            original_dict[key] = value
    return original_dict


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
