[![Build Status](https://travis-ci.org/SUNCAT-Center/catmap.svg)](https://travis-ci.org/SUNCAT-Center/catmap)
[![Software License](https://img.shields.io/badge/license-GPLv3-brightgreen.svg?style=flat-square)](COPYING.txt)
[![Documentation Status](https://readthedocs.org/projects/catmap/badge/?version=latest)](http://catmap.readthedocs.org/en/latest/?badge=latest)

# CatMAP

## INSTALLATION

### pip/setup.py

CatMap can be installed directly via pip: 

    pip install --upgrade https://github.com/SUNCAT-Center/catmap/zipball/master
    
or download/clone the repository and run 

    python setup.py install

as of the repository root folder.

### via add-to Path

To use the package add this directory to the PYTHONPATH, e.g. in bash
shell:

    export PYTHONPATH=$HOME/THIS_FOLDER_PATH:$PYTHONPATH

or in cshell:

    setenv PYTHONPATH $HOME/THIS_FOLDER_PATH:$PYTHONPATH

You will need to ensure that you are running python version 2.5 or
greater, and that the mpmath, numpy, scipy, ase, and matplotlib python
libraries are installed. A significant speed increase can also be
obtained by including the gmpy package.

The installation can be verified by starting python and typing the
following commands:

    import numpy
    import mpmath
    import matplotlib
    import catmap

See the [documentation](http://catmap.readthedocs.org) for more details
and tutorials.


## Cite
If you find CatMAP useful to your research, please cite:
> Medford, A. J., Shi, C., Hoffmann, M. J., Lausche, A. C., Fitzgibbon, S. R., Bligaard, T., & Nørskov, J. K. (2015). CatMAP: a software package for descriptor-based microkinetic mapping of catalytic trends. Catalysis Letters, 145, 794-807.

If you are using the current version of CatMAP, please also cite:
> Vijay, S., H. Heenen, H., Singh, A. R., Chan, K., & Voss, J. (2024). Number of sites‐based solver for determining coverages from steady‐state mean‐field micro‐kinetic models. Journal of Computational Chemistry, 45(9), 546-551.

which details the implementation of the numbers solver, the current default solver used by CatMAP for improved numerical stability. The behavior of previous catmap versions ( <=v0.3.2 ) can be reproduced via `use_numbers_solver = False`