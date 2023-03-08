# rna_lib_design

[![formatting: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![install with conda](https://github.com/jyesselm/rna_lib_design/actions/workflows/conda.yml/badge.svg)](https://github.com/jyesselm/rna_lib_design/actions/workflows/conda.yml)

This is package is the library design tool for the Yesselman Lab (https://yesselmanlab.com/). This tool takes a set of RNA sequences and finalizes them so they can be ordered from IDT, agilent or Twist or others. There are options to add unique barcodes to increase the diversity of libraries. There are also several checks to stop ordering libraries that will not work.

All libraries are assumed to be RNA and a T7 promoter (TTCTAATACGACTCACTA) will be added to the 5' end of each sequence. So that the library will be able to be transcribed into RNA 


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

the rld command line tool has 5 sub commands. `list` `barcode` `barcode2` `add-common` and `edit-distance`


```shell
rld --help

Usage: rld [OPTIONS] COMMAND [ARGS]...

  This is a simple tool to generate barcodes for a library of RNAs to be ordered
  as an DNA oligo library.

Options:
  --help  Show this message and exit.

Commands:
  add-common     add common p5/p3 sequences
  barcode        adds a single barcode
  barcode2       adds two barcodes
  edit-distance  compute edit distance of library
  list           lists resources available
```

`list` displays the available barcodes and common sequences that can be used in the other sub commands
`edit-distance` will compute the edit distance of the library

The other sub commands are different methods of finalizing RNA libraries so they can be ordered
as DNA oligo pools.

### understanding input 
All libraries should be supplied as csvs. Lets start with a simple example such as 

```shell
$ cat test/resources/libs/simple.csv

sequence
GGGGGGAAAACCCCCC
CCAAAACCCCUUUUGG
```

### understanding the parameter file

All parameters to be used by each sub command is are stored in parameter files. You can see all of these files in the `rna_lib_design/resources/presets` directory. The parameters are stored in a yaml file. These parameters are validated by the `jsonschema` package. The schema for the parameters is stored in `rna_lib_design/resources/schemas/`.

Here is a breakdown of the parameters.

#### build_str 

The build string defines what segments of RNA are we going to add to each sequence. For example.

```yaml
build_str: P5-P5EXT-SOI-P3EXT-P3
```

Each build string MUST contain `SOI` or sequence of interest. All other are optional but if they appear in the build string they must also be included in `segments` which we will discuss in a second. The position of each segment in the build string is important. The order of the segments will be the order they are added to the sequence of interest. So for this build string the segment of P5EXT will be added immediately to the 5' end of SOI and P3EXT will be added immediately to the 3' end of SOI. P5 will be added to the 5' end after P5EXT and etc.

#### segments
Each named segment in the build_str must be defined in another section called `segments` here each name must have its own parameters defined. There are many possibilites. using the `name` option defines a set sequence that will be added to each sequence of interest in the library. You can see the corresponding values foreahc name in the rna_lib_design/resources/named_eqs directory. 


```yaml
segments:
  P5:
    name: "uucg_p5_rev_primer"
  P3:
    name: "rt_tail"
  P3EXT:
    sequence: "AC"
    structure: ".."
  P5EXT:
    sequence: ""
    structure: ""
```


We will use this file as an example for the rest of the documentation.

### adding commmon sequences 
The add-common sub command will add 5' and 3' sequences to each sequence in a csv file.
```
rld add-common --help

#### command line interface

Usage: rld add-common [OPTIONS] CSV

  add common p5/p3 sequences

Main options:
  These are the main options for generating a library
  -t, --btype TEXT             what type of barcode to use see full list in
                               resources/presets
  --param-file PATH            supply a new param file to override specific
                               present or to manuallydetermine each option
  -o, --output TEXT            the path to save results to
  -p, --num-processes INTEGER  number of processes to run simultaneously
  --debug                      turn on debug logging for the application
  --skip-edit-dist             skip the edit distance calculation
  --trim-p5 INTEGER            trim sequence at 5' end by this length
  --trim-p3 INTEGER            trim sequence at 3' end by this length

Other options:
  --help                       Show this message and exit.
```

#### examples 
Lets add the standard p5/p3 sequences to our simple library. If no additional options are supplied the standard preset will be used. The standard preset is defined in `rna_lib_design/resources/presets/add_common_standard.yml` and can be overridden with the `--param-file` option.

All the parameters are displayed in the log output. Every single one can be changed. The parameters are stored in a yaml file and can be edited by hand. The parameters are also stored in the results directory for future reference.

```shell
$ rld add-common test/resources/libs/simple.csv 

INFO     RLD.CLI      Using csv: test/resources/libs/simple.csv
INFO     RLD.CLI      Using output dir: results
INFO     RLD.CLI      Copying test/resources/libs/simple.csv to results/input.csv
INFO     RLD.CLI      csv has 2 sequences
INFO     RLD.CLI      Writing parameters to results/params.yml
INFO     RLD.CLI      No preset or param file supplied, using standard preset
INFO     RLD.CLI      Using parameters:
{
    "build_str": "P5-P5EXT-SOI-P3EXT-P3",
    "segments": {
        "P5": {
            "name": "uucg_p5_rev_primer"
        },
        "P3": {
            "name": "rt_tail"
        },
        "P3EXT": {
            "sequence": "AC",
            "structure": ".."
        },
        "P5EXT": {
            "sequence": "",
            "structure": ""
        }
    },
    "design_opts": {
        "increase_ens_defect": 2.0,
        "max_ens_defect": 5.0,
        "max_attempts": 10,
        "max_solutions": 10,
        "score_method": "increase",
        "allowed_ss_mismatch": 2,
        "allowed_ss_mismatch_barcodes": 2
    }
}
INFO     RLD.DESIGN   starting design
INFO     RLD.DESIGN   no 'name' column was in dataframe - adding one
INFO     RLD.DESIGN   running on single core
INFO     RLD.DESIGN   no 'structure' column folding it now
INFO     RLD.SSET     P5 is using a named sequence/structure: uucg_p5_rev_primer
INFO     RLD.SSET     P5 -> SequenceStructure(sequence='GGAACAGCACUUCGGUGCAAA', structure='......((((....))))...')
INFO     RLD.SSET     P3 is using a named sequence/structure: rt_tail
INFO     RLD.SSET     P3 -> SequenceStructure(sequence='AAAGAAACAACAACAACAAC', structure='....................')
INFO     RLD.CLI      no sequences discarded
INFO     RLD.DESIGN   results/results-all.csv contains all information generated from run
INFO     RLD.DESIGN   results/results-rna.csv contains only information related to the RNA sequence
INFO     RLD.DESIGN   p5 seq -> SequenceInfo(name='uucg_p5_rev_primer', sequence='GGAACAGCACUUCGGUGCAAA', code='P0058')
INFO     RLD.CLI      the edit distance of lib is: 12.0

```



## things todo 
get docker working for automated testing
    allow vienna rna to be installed from source so it works on all operating systems
new types of barcodes?
    triple barcoding?
    double barcodes with single strands?