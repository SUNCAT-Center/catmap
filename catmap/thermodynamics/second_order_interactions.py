import catmap
from catmap import ReactionModelWrapper
from catmap.model import ReactionModel
from catmap.functions import smooth_piecewise_linear
from catmap.functions import offset_smooth_piecewise_linear
from catmap.functions import parse_constraint
from catmap.thermodynamics import FirstOrderInteractions
import pylab as plt
import numpy as np
from scipy import integrate

class SecondOrderInteractions(FirstOrderInteractions,ReactionModelWrapper):
    """Class for implementing 'first-order adsorbate interaction model. 
    Should be sub-classed by scaler."""

    def __init__(self,reaction_model=ReactionModel()):
        self._rxm = reaction_model
        defaults = dict(
                cross_interaction_mode='geometric_mean',
                transition_state_cross_interaction_mode='intermediate_state',
                interaction_response_function = 'linear',
                interaction_response_parameters = {'max_coverage':1,'cutoff':0.25,
                    'smoothing':0.05},
                interaction_fitting_mode=None,
                interaction_strength = 1,
                #weight interaction parameters by this
                )
        self._rxm.update(defaults)
        self._required = {'cross_interaction_mode':str,
                'transition_state_cross_interaction_mode':str,
                'interaction_fitting_mode':None
                }

    def get_linearizer(self):
        if self.interaction_fitting_mode:
            if self.interaction_fitting_mode == 'average_self':
                def linearizer(theta,params):
                    fi = self.interaction_response_function(theta,**params)[0]
                    return 0.5*(fi**2)*theta
                self._linearizer = linearizer
            elif self.interaction_fitting_mode == 'differential_self':
                def linearizer(theta,params):
                    fi,dfi,d2fi = self.interaction_response_function(theta,**params)
                    return fi*dfi*theta**2 + (fi**2)*theta
                self._linearizer = linearizer
            else:
                raise AttributeError('Invalid interaction fitting mode.')

    def params_to_matrix(self,param_vector):
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
                if not epsilon_matrix[i,j]:
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
        
        f,df,d2f = self.interaction_response_function(1.0)
        beta = df/f #Need to make this site-specific
        for i,TS_params in enumerate(list(self.interaction_transition_state_scaling_matrix)):
            i += n_ads
            TS_params = list(TS_params)
            for j, epsilon_params in enumerate(list(epsilon_matrix[0:n_ads,0:n_ads])):
                e_TS = 0
                Nn = 2 #Need a better way of getting this number...
                for k, ep_k in enumerate(epsilon_params):
                    e_TS += ep_k*TS_params[k]
                e_TS = (1+(beta/Nn))*e_TS
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

    @staticmethod
    def smooth_piecewise_linear_response(*args,**kwargs):
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
