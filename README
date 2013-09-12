############################################################################
                               INSTALLATION
############################################################################

To use the package add this directory to the PYTHONPATH, e.g. in bash shell:

export PYTHONPATH=$HOME/THIS_FOLDER_PATH:$PYTHONPATH

or in cshell:

setenv PYTHONPATH $HOME/THIS_FOLDER_PATH:$PYTHONPATH

You will need to ensure that you are running python version 2.5 or greater,
and that the mpmath, numpy, scipy, ase, and matplotlib python libraries are 
installed. A significant speed increase can also be obtained by including the 
gmpy package. Note that if you are running this from the SUNCAT cluster all 
dependencies are already installed, but you will need to pick up the proper
environment for python2.5 by:

setenv NUMPYDIR /opt/TWWfsw/numpy12/lib/python25
setenv SCIPYDIR /opt/TWWfsw/scipy06/lib/python25
setenv MPMATH /nfs/slac/g/suncatfs/sw/external/mpmath-0.17/install/lib/python
setenv GMPY /nfs/slac/g/suncatfs/sw/external/gmpy-1.16/install/lib/python
setenv ASE /nfs/slac/g/suncatfs/sw/external/ase/3.6.1.2694
setenv PYTHONPATH ${NUMPYDIR}:${SCIPYDIR}:${MPMATH}:{GMPY}:${ASE}:${PYTHONPATH}
alias python /usr/local/bin/python2.5 #optional 

Changing the environment to python2.5 may cause other scripts to break, so use
cautiously.

The installation can be verified by starting python and typing the following
commands:

import numpy
import mpmath
import matplotlib
import kinetics
import simplekinetics
import mkm

############################################################################
                             CONTENTS
############################################################################

The directory contains three modules:

kinetics: old (deprecated) version included for compatibility for people using
the previous versions.
simplekinetics: a separate and more script-like micro-kinetic modelling module
wich also inclues adsorbate-adsorbate interactions. Written by Simon
Brodersen.
mkm: the current micro-kinetic modeling package.

In addition, the following folders are included:


demos: some example scripts, separated by which package
they correspond to. These should provide a starting point for understanding
how the code works.

doc: documentation, including a brief overview of the code structure and 
philosophy.

dev: contains dev.txt which is an intro and to-do list for anyone who wants to
help develop the code. Also home to testing scripts/utilities or anything
which is used strictly for code development.
