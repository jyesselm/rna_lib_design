import os
from pathlib import Path


def get_lib_path():
    file_path = os.path.realpath(__file__)
    spl = file_path.split("/")
    base_dir = "/".join(spl[:-2])
    return Path(base_dir)


def get_py_path():
    return get_lib_path() / "rna_lib_design"


def get_resources_path():
    return get_py_path() / "resources"


def get_test_path():
    return get_lib_path() / "test"
