[![Build Status](https://travis-ci.org/SUNCAT-Center/catmap.svg)](https://travis-ci.org/SUNCAT-Center/catmap)
[![Software License](https://img.shields.io/badge/license-GPLv3-brightgreen.svg?style=flat-square)](COPYING.txt)
[![Documentation Status](https://readthedocs.org/projects/catmap/badge/?version=latest)](http://catmap.readthedocs.org/en/latest/?badge=latest)

# CatMAP

## INSTALLATION

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

[![Analytics](https://suncat-analytics.appspot.com/UA-75027967-2/catmap/githubreadme)](https://github.com/igrigorik/ga-beacon)
