Installation
============

CatMAP is currently in alpha testing and thus only can only be installed from
source via GitHub. The code runs directly from source, so it can be "installed"
by cloning the GitHub repository. Before installing make sure that the required
dependencies are installed:

- python 2.5 or greater
- numpy 
- scipy
- matplotlib 
- mpmath 
- ase
- gmpy (optional - gives 2-3x speedup)

You can check that the dependencies are installed by starting a
python session and importing them:

.. code:: sh

    bash $ python Python >2.5 (....) 
    Type "help", "copyright", "credits" or "license" for more information.  
    >>>import mpmath
    >>>import matplotlib 
    >>>import numpy 
    >>>import scipy 
    >>>import ase
    >>>import gmpy

After ensuring that you have the dependencies, change to the
directory where you want the CatMAP source to be installed (we will
call it $CATMAP): 

.. code:: bash
    
    $ cd $CATMAP 
    $ git clone https://github.com/ajmedford/catmap.git 

This will clone the
repository into a directory called "catmap". Next, you need to add
the location of the CatMAP source code to the $PYTHONPATH
environment variable so that python knows where to find it:

BASH: 

.. code:: bash

    bash $ export PYTHONPATH=$CATMAP/catmap:$PYTHONPATH

CSHELL: 

.. code:: csh

    $ setenv PYTHONPATH=$CATMAP/catmap:$PYTHONPATH

You can verify that everything went smoothly by importing the
CatMAP module:

.. code::

    $ python >>>import catmap

Documentation (this wiki) is located in the catmap/docs folder, and
is available online at http://catmap.readthedocs.org/.
The best place to start learning how to use the code is the
:doc:`tutorials/index`.
