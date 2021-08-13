# Hierarchical Analysis of Component linKs (HACK)
This repository offers a Multi-level decomposition approach to Stabilised Independent Component Analysis for genomic data

The algorithm works in 3 steps:
* Pre-processing of given data
* Deconvolution over multiple dimensions using stabilised ICA (https://github.com/ncaptier/stabilized-ica)
* Creation of a hierarchical graph

## Installation

```
$ pip install git+https://github.com/NicolasSompairac/HACK#egg=hack
```

## Example

A Jupyter notebook [Hierarchical ICA](Hierarchical_ICA.ipynb) can be found, describing the different steps needed to be performed as well as the required packages to run.
The notebook can be used with the example data files using the following [Zenodo repository](https://zenodo.org/record/4720408).

For an easier and faster use, a deconvolution run on the example data has already been performed from 2 to 100 components.
Simply download the data and decompress them in the same repository as the Jupyter Notebook.

## Data

The data used in the example is coming from an RNA-Seq Breast Cancer dataset taken from TCGA (https://portal.gdc.cancer.gov)

## Usage

The list of required packages and their version are listed in the Jupyter notebook. The required Stabilised ICA functions ["sica"](https://github.com/ncaptier/stabilized-ica/tree/master/sica) must be installed separately.

For a direct use of the notebook, it is recommended to install the Stabilised ICA package directly via "pip".
It is also possible to download the ["sica"](https://github.com/ncaptier/stabilized-ica/tree/master/sica) folder manually and put it in the same directory as the notebook.
