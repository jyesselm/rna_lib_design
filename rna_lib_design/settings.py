import os

def get_lib_path():
    file_path = os.path.realpath(__file__)
    spl = file_path.split("/")
    base_dir = "/".join(spl[:-2])
    return base_dir

LIB_PATH = get_lib_path()
RESOURCES_PATH = LIB_PATH + "/rna_lib_design/resources/"
TEST_PATH = LIB_PATH + "/test/"
DATA_PATH = LIB_PATH + "/data/"