from setuptools import setup, find_packages
from shutil import which


def is_tool(name):
    """Check whether `name` is on PATH and marked as executable."""
    return which(name) is not None


if not is_tool("RNAfold"):
    print("RNAfold from vienna must be installed for this package to work!")
    exit()

with open("requirements.txt") as f:
    requirements = f.read().splitlines()

setup(
    name="rna_lib_design",
    version="0.0.2",
    author="Joe Yesselman",
    author_email="jyesselm@unl.edu",
    packages=["rna_lib_design"],
    py_modules=[
        "rna_lib_design/assemble",
        "rna_lib_design/barcode",
        "rna_lib_design/cli"
        "rna_lib_design/design",
        "rna_lib_design/logger",
        "rna_lib_design/params",
        "rna_lib_design/query_resources",
        "rna_lib_design/settings",
        "rna_lib_design/setup_resources",
        "rna_lib_design/structure",
        "rna_lib_design/structure_set",
        "rna_lib_design/testing",
        "rna_lib_design/util",
    ],
    include_package_data=True,
    install_requires=requirements,
    entry_points={
        "console_scripts": [
            "rld = rna_lib_design.cli:cli"
        ]
    },
)
