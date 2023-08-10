**Content of the source code documentation**  

1. Name and short description of the software, authors, date of initial development
1. Main features
1. Main requirements
1. Further information:
    1. [Input examples and explanations, step-by-step tutorial](doc/input.md)
    1. [More detailed description of scientific approach and input variables reference](doc/method.md)
    1. [Validity range of the parameters](doc/parameters.md)
    1. [License, bug tracker, references, citations](doc/further.md)
    1. [Source code description](doc/sphinxdoc.md) - functions and classes, modules, variables

## MEGS: Morphological Evaluation of Galactic Structure 

This project, developed by Ufuk Çakir, introduces the MEGS software, designed to execute a PCA-based model for the morphological evaluation of galactic structures.
*Interdisciplinary Center for Scientific Computing (IWR), Heidelberg University, 06/2023*

The MEGS software includes preprocessing scripts for galaxy images, PCA computation routines, and modules for the low-dimensional projection of these galaxy images. For usage details, see [input](doc/input.md). The methods are exhaustively described in [method](doc/method.md). More detailed information on the input parameters can be found in [parameters](doc/parameters.md). A detailed source code description is given through the sphinx documentation hosted on [Read the docs](https://megs.readthedocs.io/en/latest/).

The software requires a Python environment with scientific computing libraries such as `numpy` and `scikit-learn` installed.

For installation, run  
`source setup.sh`
