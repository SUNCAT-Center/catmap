"""Class for `mapping' equilibrium coverages and rates through
descriptor space. This class acts as a base class to be inherited
by other mapper classes, but is not functional on its own.

get_rxn_parameter_map(descriptor_ranges,resolution): Uses a
    scaler object to determine the reaction parameters as a function of
    descriptor space. May be useful for debugging or providing
    intuition about rate determining steps. Should return a list of
    the form
[[descriptor_1,descriptor_2,...],[rxn_parameter1, rxn_parameter2, ...]]

save_map(map,map_file): creates a pickle of the "map" list and dumps it
    to the map_file

load_map(map_file):  loads a "map" list by loading a pickle from
    the map_file

A functional derived mapper class must also contain the methods:

get_coverage_map(descriptor_ranges,resolution): a function which
    returns a list of the form
    [[descriptor_1,descriptor_2,...], [cvg_ads1,cvg_ads2,...]]

get_rates_map(descriptor_ranges,resolution): a function which returns
    a list of the form
    [[descriptor_1,descriptor_2,...], [rate_rxn1,rate_rxn2,...]]

"""
from matplotlib.mlab import griddata
import numpy as np
import mpmath as mp
import cPickle as pickle
import os
from copy import copy
from catmap.model import ReactionModel
from catmap import ReactionModelWrapper
from catmap import plt

class MapperBase(ReactionModelWrapper):
    # XXX : Having an instantiated object as default parameter
    # may have side-effects since every instance of MapperBase will have
    # the identical instance of ReactionModel as its attribute
    # Unless this is deliberately so, one should better use e.g. None
    # as the default value and then instantiate ReactionModel in the
    # function body of __init__ .
    def __init__(self,reaction_model=ReactionModel()):
        self._rxm = reaction_model
        self._solver_output = ['coverage','rate', #outputs requiring solver
                'turnover_frequency','selectivity','rate_control',
                'noninteracting_coverages']

    def get_point_output(self,descriptors,*args,**kwargs):
        self.solver.compile()
        self._output_variables = [v for v in self.output_variables]

        self._descriptors = descriptors
        params = self.scaler.get_rxn_parameters(descriptors)
        self._params = params

        if True in [v in self._solver_output for v in self.output_variables]:
            if 'coverage' not in self._output_variables:
                self._output_variables = ['coverage'] + self._output_variables
            elif self._output_variables[0] != 'coverage':
                self._output_variables.remove('coverage')
                self._output_variables = ['coverage'] + self._output_variables
            self.output_labels['coverage'] = self.adsorbate_names
            self.output_labels['rate'] = self.elementary_rxns
                # Need coverages for solver vars

        for out in self._output_variables:
            if getattr(self,'get_point_'+out):
                val = getattr(self,'get_point_'+out)(descriptors,*args,**kwargs)
                setattr(self,'_'+out,val)

        self.solver.set_output_attrs(params)
        self.scaler.set_output_attrs(descriptors)

        for out in self.output_variables:
            mapp = getattr(self,'_'+out+'_temp',{})
            mapp[repr(descriptors)] = getattr(self,'_'+out)
            setattr(self,'_'+out+'_temp',mapp)

    def get_output_map(self,descriptor_ranges,resolution,*args,**kwargs):
        self.solver.compile()
        self._output_variables = [v for v in self.output_variables]
        if True in [v in self._solver_output for v in self.output_variables]:
            #determines whether or not solver is needed
            if 'coverage' not in self._output_variables:
                self._output_variables = ['coverage'] + self._output_variables
                self._coverage_map = None
            elif self._output_variables[0] != 'coverage':
                self._output_variables.remove('coverage')
                self._output_variables = ['coverage'] + self._output_variables

        # Need coverages for solver vars
        ismapped = False
        for out in self._output_variables:
            if getattr(self,'get_'+out+'_map'):
                val = getattr(self,'get_'+out+'_map')(
                        descriptor_ranges,resolution,*args,**kwargs)
                setattr(self,out+'_map',val)
                ismapped = True

        if ismapped == False:
            d1Vals, d2Vals = self.process_resolution()
            for d1V in d1Vals:
                for d2V in d2Vals:
                    self._descriptors = [d1V,d2V]
                    self.get_point_output(self._descriptors)

        for out in self.output_variables:
#            if getattr(self,out+'_map'):
#                mapp = getattr(self,out+'_map')
#            else:
            map_dict = getattr(self,'_'+out+'_temp',[])
            mapp = []
            for key in map_dict:
                mapp.append([eval(key),map_dict[key]])
            setattr(self,out+'_map',mapp)

            if getattr(self,out+'_map_file'):
                outfile = getattr(self,out+'_map_file')
                self.save_map(mapp,outfile)
    
    def process_resolution(self, descriptor_ranges = None, resolution = None):
        if not descriptor_ranges:
            descriptor_ranges = self.descriptor_ranges
        if resolution is None:
            resolution = self.resolution
        resolution = np.array(resolution)
        if resolution.size == 1:
            resx = resy = float(resolution)
        elif resolution.size ==2:
            resx = resolution[0]
            resy = resolution[1]
        else:
            raise ValueError('Resolution is not the correct shape')

        d1min, d1max = descriptor_ranges[0]
        d2min, d2max = descriptor_ranges[1]
        d1Vals = np.linspace(d1min, d1max, resx)
        d2Vals = np.linspace(d2min, d2max, resy)
        return d1Vals, d2Vals

