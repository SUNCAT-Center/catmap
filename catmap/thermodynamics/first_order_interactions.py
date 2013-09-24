import catmap
from catmap import ReactionModelWrapper
from catmap.model import ReactionModel
from catmap.functions import smooth_piecewise_linear
from catmap.functions import parse_constraint
import pylab as plt
import numpy as np
from scipy import integrate

class FirstOrderInteractions(ReactionModelWrapper):
    """Class for implementing 'first-order adsorbate interaction model. 
    Should be sub-classed by scaler."""

    def __init__(self,reaction_model=ReactionModel()):
        self._rxm = reaction_model
        defaults = dict(
                cross_interaction_mode='geometric_mean',
                transition_state_cross_interaction_mode='intermediate_state',
                max_self_interaction = 'Pd',
                interaction_response_function = 'linear',
                interaction_response_parameters = {'max_coverage':1,'cutoff':0.25,
                    'smoothing':0.05},
                interaction_fitting_mode=None,
                default_interaction_constraints = [None]*(len(self.descriptor_names)+1),
                )
        self._rxm.update(defaults)
        self._required = {'cross_interaction_mode':str,
                'transition_state_cross_interaction_mode':str,
                'interaction_fitting_mode':None
                }

    def parameterize_interactions(self):
        if self.interaction_fitting_mode:
            if self.interaction_fitting_mode == 'average_self':
                def linearizer(theta,params):
                    theta += 1e-10 # avoid zero-division errors
                    Fint = lambda x: self.interaction_response_function(x,**params)[0]*x
                    integrated =  integrate.quad(Fint,0,theta)
                    return integrated[0]/theta
            elif self.interaction_fitting_mode == 'differential_self':
                def linearizer(theta,params):
                    return self.interaction_response_function(theta,**params)[0]*theta
            else:
                raise AttributeError('Invalid interaction fitting mode.')
            self.self_interaction_parameter_dict = {}
            if self.max_self_interaction:
                if self.max_self_interaction in self.surface_names:
                    max_idx = self.surface_names.index(self.max_self_interaction)
                else:
                    raise TypeError('max_self_interaction must be the name of a '+\
                        'surface or a dictionary with species names as keys and maxs as vals.')
            else:
                max_idx = None

            for sp in self.species_definitions:
                site = self.species_definitions[sp]['site']
                cov_dep_E = self.species_definitions[sp].get('coverage_dependent_energy',None)
                if cov_dep_E:
                    params = self.species_definitions[site].get('interaction_response_parameters',
                            self.interaction_response_parameters)
                    int_params = []
                    for theta_E in cov_dep_E:
                        if theta_E and len(theta_E) >= 2:
                            theta,E = zip(*theta_E)
                            theta_l = [linearizer(t_i,params) for t_i in theta]
                            m,b = plt.polyfit(theta_l,E,1)
                            int_params.append(m)
                        else:
                            int_params.append(None)
                    if max_idx:
                        max_p = int_params[max_idx]
                        if max_p is None:
                            print 'Warning: No interaction parameter for '+sp+' on '+\
                                    self.max_self_interaction+'. No maximum will be used'
                        else:
                            self.species_definitions[sp]['max_self_interaction'] = max_p

                    self.species_definitions[sp]['self_interaction_parameter'] = int_params
                    self.self_interaction_parameter_dict[sp] = int_params
                    all_ads = self.adsorbate_names + self.transition_state_names
                    self.parameter_names = list(all_ads)
                    for pi in all_ads:
                        for pj in all_ads:
                            self.parameter_names.append(pi + '&' + pj)

    def get_interaction_scaling_matrix(self):
        interaction_dict = {}
        n_ads = len(self.adsorbate_names)
        cross_term_names = []
        for a in self.adsorbate_names:
            interaction_dict[a] = self.species_definitions[a]['self_interaction_parameter']
            cross_params = self.species_definitions[a].get('cross_interaction_parameters',{})
            for cp in cross_params:
                if cp not in self.adsorbate_names+self.transition_state_names:
                    raise ValueError('Cross parameter name must be in adsorbate names. The '+\
                            'name ' + cp + ' is not in ' + str(self.adsorbate_names+self.transition_state_names))
                params = cross_params[cp]
                if len(params) != len(self.surface_names):
                    raise ValueError('Cross parameters must be specified for each surface. '+\
                            str(len(params)) + ' parameters were specified, but there are '+\
                            str(len(self.surface_names)) + 'surfaces for '+ a +','+cp+'.' )
                ads_a, ads_b = sorted([a,cp])
                name = '&'.join([ads_a,ads_b])
                if name not in cross_term_names:
                    cross_term_names.append(name)
                    interaction_dict[name] = params
        
        constraint_dict = {}
        for ads in self.scaling_constraint_dict:
            if '-' not in ads:
                constr = self.scaling_constraint_dict[ads]
                new_constr = []
                #preserve only 0 constraints
                for ci in constr:
                    if ci != 0:
                        new_constr.append(None)
                    else:
                        new_constr.append(0)
                constraint_dict[ads] = new_constr

        for ads in self.scaling_constraint_dict:
            if '&' in ads:
                a,b = ads.split('&')
                a,b = sorted([a,b])
                new_ads = '&'.join([a,b])
                constraint_dict[new_ads] = self.scaling_constraint_dict[ads]
            else:
                constraint_dict[ads] = self.scaling_constraint_dict[ads]

        cross_term_names = tuple(cross_term_names) 
        param_names = self.adsorbate_names + cross_term_names

        if cross_term_names:
            self.interaction_cross_term_names = cross_term_names

        #get mins/maxs
        interaction_mins = []
        interaction_maxs = []
        for p in param_names:
            if p not in constraint_dict:
                constr = self.default_interaction_constraints
            else:
                constr = constraint_dict[p]
            minvs,maxvs = parse_constraint(constr,p)
            interaction_mins.append(minvs)
            interaction_maxs.append(maxvs)

        #get coefficient matrix
        C = catmap.functions.scaling_coefficient_matrix(
                interaction_dict, self.descriptor_dict, 
                self.surface_names, 
                param_names,
                interaction_mins,
                interaction_maxs)
        self.interaction_coefficient_matrix = C.T

        #get TS scaling matrix equivalent
        if self.transition_state_cross_interaction_mode == 'transition_state_scaling':
            TS_weight_matrix = []
            for params in list(self.transition_state_scaling_matrix):
                TS_weight_matrix.append(params[:-1])

        else:
            #take TS interaction parameters as some linear combination of initial/final state
            if self.transition_state_cross_interaction_mode == 'initial_state':
                weight = 0
            elif self.transition_state_cross_interaction_mode == 'final_state':
                weight = 1
            elif self.transition_state_cross_interaction_mode.startswith('intermediate_state'):
                if self.transition_state_cross_interaction_mode == 'intermediate_state':
                    weight = 0.5
                else:
                    mname = self.transition_state_cross_interaction_mode
                    junk,w = mname.rsplit('(',1)
                    w,junk = w.rsplit(')',1)
                    try:
                        weight = float(w)
                    except:
                        raise AttributeError('Intermediate state weight is not floatable:'+w)
            else:
                raise AttributeError('Undefined transition_state_cross_interaction_mode: ' +\
                        self.transition_state_cross_interaction_mode)
            TS_weight_matrix = self.get_TS_weight_matrix(weight)

        self.interaction_transition_state_scaling_matrix = TS_weight_matrix

        return self.interaction_coefficient_matrix

    def get_interaction_matrix(self,descriptors):
        full_descriptors = list(descriptors) + [1.]
        param_vector = np.dot(self.coefficient_matrix,full_descriptors)

        n_ads = len(self.adsorbate_names)
        n_TS = len(self.transition_state_names)
        all_names = self.adsorbate_names + self.transition_state_names
        n_tot = len(all_names)
        epsilon_matrix = np.zeros((n_tot,n_tot))
        self_interactions = param_vector[n_tot:n_tot+n_ads]
        cross_interactions = param_vector[n_tot+n_ads:]
        for ads,param in zip(self.adsorbate_names,self_interactions):
            max_val = self.species_definitions[ads].get('max_self_interaction',None)

            if max_val is not None:
                param = min(param,max_val)
                self_interactions[self.adsorbate_names.index(ads)] = param

        for i,e_ii in enumerate(self_interactions):
            epsilon_matrix[i,i] = e_ii

        for i, e_i in enumerate(self_interactions):
            for j, e_j in enumerate(self_interactions):
                if self.cross_interaction_mode == 'geometric_mean':
                    epsilon_matrix[i,j] = np.sqrt(abs(e_i)*abs(e_j))
                elif self.cross_interaction_mode == 'arithmetic_mean':
                    epsilon_matrix[i,j] = (e_i+e_j)/2.
                elif self.cross_interaction_mode == 'neglect':
                    epsilon_matrix[i,j] = 0

        if self.non_interacting_site_pairs:
            for site_1, site_2 in self.non_interacting_site_pairs:
                idxs1 = [self.adsorbate_names.index(sp) for 
                        sp in self.adsorbate_names if
                        self.species_definitions[sp]['site']==site_1]
                idxs2 = [self.adsorbate_names.index(sp) for 
                        sp in self.adsorbate_names if
                        self.species_definitions[sp]['site']==site_2]
                for i in idxs1:
                    for j in idxs2:
                        epsilon_matrix[i,j] = epsilon_matrix[j,i] = 0
        
        for i,TS_params in enumerate(list(self.interaction_transition_state_scaling_matrix)):
            i += n_ads
            TS_params = list(TS_params)
            for j, epsilon_params in enumerate(list(epsilon_matrix[0:n_ads,0:n_ads])):
                e_TS = np.dot(TS_params,epsilon_params)
                epsilon_matrix[i,j] = epsilon_matrix[j,i] = e_TS

        if len(cross_interactions)>0:
            for name, param in zip(self.interaction_cross_term_names,cross_interactions):
                a,b = name.split('&')
                i = all_names.index(a)
                j = all_names.index(b)
                epsilon_matrix[i,j] = epsilon_matrix[j,i] = param

        epsilon_matrix *= self.interaction_strength
        self._interaction_matrix = epsilon_matrix.tolist()
        return epsilon_matrix

    def get_TS_weight_matrix(self,weight):
        weight_matrix = []
        for TS in self.transition_state_names:
            weight_TS = [0]*len(self.adsorbate_names)
            IS = None
            FS = None
            for rxn in self.elementary_rxns:
                if TS in rxn[1]:
                    if IS or FS:
                        print('Warning: Ambiguous initial/final state for'+TS)
                    IS = rxn[0]
                    FS = rxn[-1]
            for ads in IS:
                if ads in self.adsorbate_names:
                    idx = self.adsorbate_names.index(ads)
                    weight_TS[idx] += (1-weight)
            for ads in FS:
                if ads in self.adsorbate_names:
                    idx = self.adsorbate_names.index(ads)
                    weight_TS[idx] += weight
            weight_matrix.append(weight_TS)
        return weight_matrix

    @staticmethod
    def smooth_piecewise_linear_response(*args,**kwargs):
        return smooth_piecewise_linear(*args,**kwargs)

    @staticmethod
    def piecewise_linear_response(*args,**kwargs):
        kwargs['smoothing'] = 0
        return smooth_piecewise_linear(*args,**kwargs)

    @staticmethod
    def linear_response(*args,**kwargs):
        kwargs['smoothing'] = 0
        kwargs['cutoff'] = 0
        return smooth_piecewise_linear(*args,**kwargs)
