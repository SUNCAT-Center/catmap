import catmap
from catmap import ReactionModelWrapper
from catmap.model import ReactionModel
from catmap.functions import smooth_piecewise_linear
from catmap.functions import offset_smooth_piecewise_linear
from catmap.functions import parse_constraint
from catmap.thermodynamics import FirstOrderInteractions
import pylab as plt
import numpy as np
try:
    from scipy import integrate
except ImportError:
    integrate = None

class SecondOrderInteractions(FirstOrderInteractions,ReactionModelWrapper):
    """Class for implementing 'first-order adsorbate interaction model. 
    Should be sub-classed by scaler."""

    def __init__(self,reaction_model=ReactionModel()):
        FirstOrderInteractions.__init__(self,reaction_model)

    @staticmethod
    def smooth_piecewise_linear_response(*args,**kwargs):
        #Note these need to override first-order functions
        #since second-derivatives are needed
        return smooth_piecewise_linear(*args,**kwargs)

    @staticmethod
    def offset_smooth_piecewise_linear_response(*args,**kwargs):
        return offset_smooth_piecewise_linear(*args,**kwargs)

    @staticmethod
    def piecewise_linear_response(*args,**kwargs):
        kwargs['smoothing'] = 0
        return smooth_piecewise_linear(*args,**kwargs)

    @staticmethod
    def linear_response(*args,**kwargs):
        kwargs['smoothing'] = 0
        kwargs['cutoff'] = 0
        return smooth_piecewise_linear(*args,**kwargs)
