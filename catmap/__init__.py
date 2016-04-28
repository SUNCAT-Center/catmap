#Standard dependencies
import os
import sys
import inspect
import time
import cPickle as pickle
import re
from copy import copy
from string import Template

#Non-standard dependencies
import numpy as np
try:
    from scipy.interpolate import InterpolatedUnivariateSpline as spline
except ImportError:
    def spline_wrapper(x_data, y_data, k=3):  # input kwarg k is intentionally ignored
        # behaves like scipy.interpolate.InterpolatedUnivariateSpline for k=1
        def spline_func(x):
            return np.interp(x, map(float,x_data), map(float,y_data))  # loss of precision here
        return spline_func
    spline = spline_wrapper

import matplotlib as mpl
mpl.use('Agg')
import pylab as plt
import matplotlib.transforms as mtransforms
from matplotlib.mlab import griddata as mlab_griddata

def griddata(*args, **kwargs):
    """Wrapper function to avoid annoying griddata errors"""
    try:
        return mlab_griddata(*args, **kwargs)
    except RuntimeError:
        kwargs['interp'] = 'linear'
        return mlab_griddata(*args, **kwargs)

import mpmath as mp
from ase.atoms import string2symbols
from ase.thermochemistry import IdealGasThermo, HarmonicThermo
from ase.structure import molecule
from catmap.model import ReactionModel
import data

__version__ = "0.2.270"

def griddata(*args, **kwargs):
    """Wrapper function to avoid annoying griddata errors"""
    try:
        return mlab_griddata(*args, **kwargs)
    except RuntimeError:
        kwargs['interp'] = 'linear'
        return mlab_griddata(*args, **kwargs)

def load(setup_file):
    rxm = ReactionModel(setup_file = setup_file)
    return rxm

modified = []
class ReactionModelWrapper:
    def __getattribute__(self,attr):
        "Force use of custom getattr"
        return self.__getattr__(self,attr)

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

