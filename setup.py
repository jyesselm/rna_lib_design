from setuptools import setup, find_packages


with open("requirements.txt") as f:
    requirements = f.read().splitlines()

setup(
    name="rna_lib_design",
    version="0.0.3",
    author="Joe Yesselman",
    author_email="jyesselm@unl.edu",
    packages=["rna_lib_design"],
    py_modules=[
        "rna_lib_design/cli" "rna_lib_design/design",
        "rna_lib_design/logger",
        "rna_lib_design/params",
        "rna_lib_design/settings",
        "rna_lib_design/setup_resources",
        "rna_lib_design/structure_set",
        "rna_lib_design/util",
    ],
    include_package_data=True,
    #install_requires=requirements,
    entry_points={"console_scripts": ["rld = rna_lib_design.cli:cli"]},
)
