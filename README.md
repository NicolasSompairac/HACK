# Hierarchical_ICA
This repository offers a Hierarchical approach to Stabilised Independent Component Analysis for genomic data

The algorithm works in 3 steps:
* Pre-processing of given data
* Deconvolution over multiple dimensions using stabilised ICA (https://github.com/ncaptier/Stabilized_ICA)
* Creation of a hierarchical graph

## Installation

```
$ pip install git+https://github.com/NicolasSompairac/Hierarchical_ICA#egg=hica
```

## Example

A Jupyter notebook [Hierarchical ICA](Hierarchical_ICA.ipynb) can be found, describing the different steps needed to be performed as well as the required packages to run.
The notebook can be used with the example data found as .zip files in the folder [example_data](example_data). For faster use, a deconvolution run on the example data has already been performed from 2 to 100 components and can be found as .zip files in the folder [example_deconvolutions](example_deconvolutions).

## Data

The data used in the example is an RNA-Seq Breast Cancer dataset taken from TCGA (https://portal.gdc.cancer.gov)

## Usage

The list of required packages and their version are listed in the Jupyter notebook. The required Stabilised ICA folder ["sica"](https://github.com/ncaptier/Stabilized_ICA/tree/master/sica) must be downloaded separately and put in the same folder as the notebook.
For a direct use of the notebook, it is recommended to put the example data in the same directory as the notebook as well.
