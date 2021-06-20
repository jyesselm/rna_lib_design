from setuptools import setup, find_packages

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name='seq_tools',
    version='0.0.2',
    author='Joe Yesselman',
    author_email='jyesselm@unl.edu',
    packages=["seq_tools"],
    py_modules=["seq_tools/motif", "seq_tools/structure", "seq_tools/barcoder",
                "seq_tools/logger", "seq_tools/util", "seq_tools/vienna", "seq_tools/clt"],
    include_package_data=True,
    install_requires=requirements,
    entry_points={
        'console_scripts': [
            'seq_tools = seq_tools.clt:main',
            'gen_opool = seq_tools.clt:gen_opool'
        ]
    }

)