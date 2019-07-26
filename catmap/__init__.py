# Standard dependencies
import os
import sys
import inspect
import time
try:
    import cPickle as pickle
except (ImportError, ModuleNotFoundError):
    import _pickle as pickle

import re
from copy import copy
from string import Template

# Non-standard dependencies
import numpy as np
try:
    from scipy.interpolate import InterpolatedUnivariateSpline as spline
except ImportError:
    # input kwarg k is intentionally ignored.
    def spline_wrapper(x_data, y_data, k=3):
        # behaves like scipy.interpolate.InterpolatedUnivariateSpline for k=1
        def spline_func(x):
            # loss of precision here
            return np.interp(x, map(float,x_data), map(float,y_data))
        return spline_func
    spline = spline_wrapper

import matplotlib as mpl
mpl.use('Agg')
import pylab as plt
import matplotlib.transforms as mtransforms
from scipy.interpolate import griddata as sp_griddata

def griddata(*args, **kwargs):
    """Wrapper function to avoid annoying griddata errors"""
    try:
        return sp_griddata(*args, **kwargs)
    except RuntimeError:
        kwargs['interp'] = 'linear'
        return sp_griddata(*args, **kwargs)

import mpmath as mp
try:
    from ase.symbols import string2symbols
except:
    from ase.atoms import string2symbols
from ase.thermochemistry import IdealGasThermo, HarmonicThermo
try:
    from ase.build import molecule
except ImportError:
    from ase.structure import molecule
from ase.thermochemistry import IdealGasThermo, HarmonicThermo, HinderedThermo
from catmap.model import ReactionModel
from . import data

__version__ = "0.3.1"

def load(setup_file):
    rxm = ReactionModel(setup_file = setup_file)
    return rxm

modified = []
class ReactionModelWrapper:
    #def __getattribute__(self,attr):
        #"Force use of custom getattr"
        #return object.__getattr__(self,attr)

    def __getattr__(self,attr):
        "Return the value of the reaction model instance if its there. Otherwise return the instances own value (or none if the instance does not have the attribute defined and the attribute is not private)"
        if attr == '_rxm':
            return object.__getattribute__(self,attr)

        elif hasattr(self._rxm,attr):
            return getattr(self._rxm,attr)
        else:
            if attr in self.__dict__:
                val =  object.__getattribute__(self,attr)
                del self.__dict__[attr]
                #this makes sure that the attr is read from _rxm
                setattr(self._rxm,attr,val)
                return val
            elif attr.startswith('_'):
                raise AttributeError("Attribute {attr} in invalid".format(**locals()))
            else:
                return None

    def __setattr__(self,attr,val):
        "Set attribute for the instance as well as the reaction_model instance"
        accumulate = ['_required','_log_strings','_function_strings']
        if attr == '_rxm':
            self.__dict__[attr] = val
        elif attr in accumulate:
            self._rxm.__dict__[attr].update(val)
        else:
            setattr(self._rxm,attr,val)

