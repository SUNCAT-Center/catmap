#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
Catalysis Micro-kinetic Analysis Package (CatMAP)
"""

import os
import sys

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

# from catmap import __version__ as version
__version__ = "0.4.1"
__python_version__ = sys.version


setup(
    author="Andrew J. Medford",
    author_email="ajmedfor@slac.stanford.edu",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Environment :: X11 Applications :: GTK",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Developers",
        "OSI Approved :: GNU General Public License (GPL)",
        "Natural Language :: English",
        "Operating System :: POSIX :: Linux",
        "Operating System :: POSIX :: Windows",
        "Programming Language :: Fortran",
        "Programming Language :: Python",
        "Topic :: Education",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Visualization",
    ],
    description=__doc__,
    license="COPYING.txt",
    long_description=open("README.md", "r", encoding="utf-8").read(),
    maintainer="Andrew J. Medford",
    maintainer_email="ajmedfor@slac.stanford.edu",
    name="python-catmap",
    package_data={"catmap": []},
    package_dir={"catmap": "catmap"},
    packages=[
        "catmap",
        "catmap.analyze",
        "catmap.api",
        "catmap.data",
        "catmap.mappers",
        "catmap.parsers",
        "catmap.scalers",
        "catmap.solvers",
        "catmap.thermodynamics",
    ],
    platforms=["linux", "windows"],
    scripts=[] if os.name == "nt" else ["tools/catmap"],
    url="https://github.com/ajmedford/catmap",
    version=__version__,
    install_requires=[
        "matplotlib",
        "mpmath",
        "numpy",
        "graphviz",
        "ase" if sys.version_info >= (3, 5) else "ase==3.17",
    ],
    tests=[
        "pgtest~=1.3",
        "pytest~=6.0",
        "pytest-regressions~=1.0",
    ],
    extras_require={
        "tests": [
            "pgtest~=1.3",
            "pytest~=6.0",
            "pytest-regressions~=1.0",
        ],
    },
)
