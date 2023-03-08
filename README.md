# rna_lib_design

[![formatting: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)


This is package is the library design tool for the Yesselman Lab (https://yesselmanlab.com/) 


## how to install 

It is highly recommended to install with requirements in a conda environment

```shell
conda create --name py3 python=3.8
conda activate py3 
```

This package is not meant to be on pypi atm so the option is local install
```shell
git clone https://github.com/jyesselm/rna_lib_design.git
cd rna_lib_design
pip install .
```
This will create a executable in your path called `rld` 

## how to use 

```shell
rld --help


```

## things todo 
get docker working for automated testing
    allow vienna rna to be installed from source so it works on all operating systems
new types of barcodes?
    triple barcoding?
    double barcodes with single strands?