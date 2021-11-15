from setuptools import setup, find_packages

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
        name='rna_lib_design',
        version='0.0.2',
        author='Joe Yesselman',
        author_email='jyesselm@unl.edu',
        packages=["rna_lib_design"],
        py_modules=[
            "rna_lib_design/structure",
            "rna_lib_design/logger",
            "rna_lib_design/util",
            "rna_lib_design/structure_dict"
            "rna_lib_design/design"
        ],
        include_package_data=True,
        install_requires=requirements,
        entry_points={
        }

)
