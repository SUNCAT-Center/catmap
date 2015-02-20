import catmap
from catmap import ReactionModelWrapper
from catmap.model import ReactionModel
from catmap.functions import smooth_piecewise_linear
from catmap.functions import parse_constraint
import pylab as plt
import numpy as np
from scipy import integrate
from scipy.optimize import fmin

class FirstOrderInteractions(ReactionModelWrapper):
    """Class for implementing 'first-order adsorbate interaction model. 
    Should be sub-classed by scaler."""

    def __init__(self,reaction_model=ReactionModel()):
        self._rxm = reaction_model
        defaults = dict(
                cross_interaction_mode='geometric_mean',
                transition_state_cross_interaction_mode='intermediate_state',
                default_self_interaction_parameter = 0,
                interaction_response_function = 'linear',
                interaction_response_parameters = {'max_coverage':1,'cutoff':0.25,
                    'smoothing':0.05},
                interaction_fitting_mode=None,
                input_overrides_fit = True, #user inputs override fitted values
                interaction_strength = 1,
                #weight interaction parameters by this
                )
        self._rxm.update(defaults)
        self._required = {'cross_interaction_mode':str,
                'transition_state_cross_interaction_mode':str,
                'interaction_fitting_mode':None
                }

    def parameterize_interactions(self):
        self._parameterized = True
        self.get_interaction_transition_state_scaling_matrix()
        
        if self.interaction_fitting_mode is not None:
            self.fit()

    def get_energy_error(self, epsilon_ij, theta, Ediff, Eint, 
            parameter_name, surface_name):
        i_surf = self.surface_names.index(surface_name)
        all_ads = self.adsorbate_names + self.transition_state_names
        E0_list = [self.species_definitions[sp]['formation_energy'][i_surf] 
                    for sp in all_ads]

        if parameter_name in self.interaction_response_parameters:

            if self.interaction_cross_term_names:
                param_names = self.adsorbate_names + self.interaction_cross_term_names
            else:
                param_names = self.adsorbate_names

            info = self.thermodynamics.adsorbate_interactions.get_interaction_info()
            params = [info[pi][i_surf] for pi in param_names]
            int_strength = self.interaction_strength
            self.interaction_strength = 1.0
            eps_matrix = self.params_to_matrix(E0_list+params)
            self.interaction_strength = int_strength
            eps_list = list(eps_matrix.ravel())

            old_params = self.interaction_response_parameters.copy()
            self.interaction_response_parameters[parameter_name] = epsilon_ij
            self.solver._compiled = False
            self.solver.compile()

            if 'differential' in self.interaction_fitting_mode:
                diff_idx = theta.index(min([t for t in theta if t != 0]))
                Ediff_model = self.interaction_function(theta,E0_list,eps_list,
                        self.interaction_response_function,False)[0][diff_idx]
                diff_err = (Ediff-Ediff_model)
            else:
                diff_err = None

            if 'integral' in self.interaction_fitting_mode:
                Eint_model = self.integral_interaction_function(theta,E0_list,
                        eps_list,self.interaction_response_function,False)
                int_err = (Eint-Eint_model)
            else:
                int_err = None

            self.interaction_response_parameters = old_params
            self.solver._compiled = False
            self.solver.compile()


        else:
            if '&' not in parameter_name:
                i = j = all_ads.index(parameter_name)
            else:
                ads_a, ads_b = parameter_name.split('&')
                i = all_ads.index(ads_a)
                j = all_ads.index(ads_b)

            if self.interaction_cross_term_names:
                param_names = self.adsorbate_names + self.interaction_cross_term_names
            else:
                param_names = self.adsorbate_names

            info = self.thermodynamics.adsorbate_interactions.get_interaction_info()
            params = [info[pi][i_surf] for pi in param_names]
            int_strength = self.interaction_strength
            self.interaction_strength = 1.0
            eps_matrix = self.params_to_matrix(E0_list+params)
            self.interaction_strength = int_strength
            eps_matrix[i,j] = eps_matrix[j,i] = epsilon_ij
            eps_list = list(eps_matrix.ravel())

            if 'differential' in self.interaction_fitting_mode:
                Ediff_model = self.interaction_function(theta,E0_list,eps_list,
                        self.interaction_response_function,False)[0][i]
                diff_err = (Ediff-Ediff_model)
            else:
                diff_err = None

            if 'integral' in self.interaction_fitting_mode:
                Eint_model = self.integral_interaction_function(theta,E0_list,
                        eps_list,self.interaction_response_function,False)
                int_err = (Eint-Eint_model)
            else:
                int_err = None

        return diff_err, int_err

    def error_norm(self,diff_err, int_err):
        if diff_err is None:
            return abs(int_err)
        elif int_err is None:
            return abs(diff_err)
        elif self.integral_error_fitting_weight:
            w = self.integral_error_fitting_weight
            return np.sqrt(sum([diff_err**2 + w*int_err**2])) 
        else:
            return abs(diff_err) + abs(int_err)
#            return np.sqrt(sum([diff_err**2 + int_err**2])) 

    def fit_interaction_parameter(self,theta_list,E_diffs,E_ints,param_name,surf_name):

       def target_function(eps_ij, theta_list, E_diffs,E_ints,param_name,surf_name):
           error = 0
           for theta,Ed,Ei in zip(theta_list,E_diffs,E_ints):
               d_err, i_err = self.get_energy_error(eps_ij, theta,Ed,Ei,param_name,surf_name)
               err_norm = self.error_norm(d_err, i_err)
               error += err_norm
#           print param_name, eps_ij, error/len(theta_list), len(theta_list)
           return error/len(theta_list)
        
       minimize = lambda x: target_function(x,theta_list,E_diffs,E_ints,param_name,surf_name)

       if param_name in self.interaction_response_parameters:
           x0 = self.interaction_response_parameters[param_name]
       else:
           x0 = 0

       print param_name
       eps_ij = fmin(minimize,[x0],disp=True)[0]
       print eps_ij
       return eps_ij


    def fit(self):
        all_adsnames = self.adsorbate_names+self.transition_state_names
        f_list = []
        E0_list = []

        #Remove
        for ads in all_adsnames:
            site = ads.split('_')[-1]
            total_cvg = self.species_definitions[site].get('total',1)
            F_params = self.species_definitions[site].get('interaction_response_parameters',
                    self.interaction_response_parameters)
            if 'max_coverage' not in F_params:
                F_params['max_coverage'] = total_cvg
            else:
                F_params['max_coverage'] *= total_cvg
            fi = lambda x: self.interaction_response_function(x,**F_params)
            f_list.append(fi)

        for i, surf in enumerate(self.surface_names):
            fitting_info = {}
            for sp in self.species_definitions:
                cvg_dep = self.species_definitions[sp].get('coverage_dependent_energy',
                        [None]*len(self.surface_names))[i]
                if cvg_dep:
                    for theta_E in cvg_dep:
                        theta,Ed,Ei = theta_E
                        n_theta = sum(ti != 0 for ti in theta)
                        ads_idx = all_adsnames.index(sp)
                        if n_theta == 1:
                            if sp in fitting_info:
                                if theta_E not in fitting_info[sp]:
                                    fitting_info[sp].append(theta_E)
                            else:
                                fitting_info[sp] = [theta_E]
                        elif n_theta == 2:
                            coads_idx = [j for j,cvg in enumerate(theta) if cvg and j != ads_idx][0]
                            coads = all_adsnames[coads_idx]
                            ads_a, ads_b = sorted([sp,coads])
                            name = '&'.join([ads_a,ads_b])
                            if name in fitting_info:
                                if theta_E not in fitting_info[name]:
                                    fitting_info[name].append(theta_E)
                            else:
                                fitting_info[name] = [theta_E]
                        else:
                            print('Warning: Ignoring coverage dependent entry for '+sp)
            
            for key in fitting_info.keys():
                if '&' not in key and len(fitting_info[key]) == 1:
                    #move into cross-parameter info if possible
                    cross_keys = [k for k in fitting_info.keys() if ('&' in k and key in k.split('&'))]
                    for k in cross_keys:
                        fitting_info[k] = fitting_info[key] + fitting_info[k]

                    del fitting_info[key] #a single entry is the same one used for Ef
                    
                elif self.input_overrides_fit: #check to make sure the user didn't input 
                    if '&' not in key:
                        if self.species_definitions[key].get('self_interaction_parameter',
                                [None]*len(self.surface_names))[i] is not None:
                            del fitting_info[key]
                    else:
                        ads_a, ads_b = key.split('&')
                        #have to check for cross interaction specified either way
                        ab_crossint = self.species_definitions[ads_a].get(
                               'cross_interaction_parameters',{})
                        if ads_b in ab_crossint:
                            del fitting_info[key]
                        ba_crossint = self.species_definitions[ads_b].get(
                               'cross_interaction_parameters',{})
                        if ads_a in ba_crossint:
                            del fitting_info[key]

            self_params = [k for k in fitting_info if '&' not in k]
            cross_params = [k for k in fitting_info if '&' in k]
            E0_list = [self.species_definitions[sp]['formation_energy'][i] for sp in all_adsnames]
            fitted_params = {}
            for key in self_params:
                thetas,Ediffs,Eints = zip(*fitting_info[key])
                eps_ii = self.fit_interaction_parameter(thetas,Ediffs,Eints,key,surf)
                eps = self.species_definitions[key].get('self_interaction_parameter',
                        [None]*len(self.surface_names))
                eps[i] = eps_ii
                fitted_params[key] = [eps_ii, thetas,Ediffs,Eints]
                self.species_definitions[key]['self_interaction_parameter'] = eps
            
            all_thetas = []
            all_Ediffs = []
            all_Eints = []
            for key in cross_params:
                thetas,Ediffs,Eints = zip(*fitting_info[key])
                all_thetas +=thetas
                all_Ediffs +=Ediffs
                all_Eints +=Eints

                ads_a, ads_b = key.split('&')
                
                ##MUST BE MORE GENERAL##
                if ads_a == 'CO_s':
                    ads_a, ads_b = ads_b, ads_a

                fit_param = '&'.join([ads_a,ads_b])
                eps_ij = self.fit_interaction_parameter(thetas,Ediffs,Eints,
                        fit_param,surf)
                eps_dict_b = self.species_definitions[ads_b].get(
                        'cross_interaction_parameters',
                        {})
                if ads_a not in eps_dict_b:
                    eps_dict_b[ads_a] = [None]*len(self.surface_names)
                eps_dict_b[ads_a][i] = eps_ij
                fitted_params[key] = [eps_ij, thetas,Ediffs,Eints]
                self.species_definitions[ads_b]['cross_interaction_parameters'] = eps_dict_b
            
#            for key in ['offset']:#,'smoothing','cutoff']:
            for key in []:#,'smoothing','cutoff']:
                param = self.fit_interaction_parameter(all_thetas,all_Ediffs,all_Eints,key,surf)
                print key, param
                self.interaction_response_parameters[key] = param

            self.solver._compiled=False
            self.solver.compile()

            ###TEMPORARY FITTING CHECK## 
            n = len(fitting_info)
            fig_diff = plt.figure()
            fig_int = plt.figure()
            n_x = int(np.sqrt(n))
            n_y = n_x+1
            info = self.thermodynamics.adsorbate_interactions.get_interaction_info()
            cross_names = self.interaction_cross_term_names
            if cross_names:
                param_names = self.adsorbate_names + cross_names
            else:
                param_names = self.adsorbate_names
            params = [info[pi][0] for pi in param_names]
            int_strength = self.interaction_strength
            self.interaction_strength = 1.0
            eps_matrix = self.params_to_matrix(E0_list+params)
            self.interaction_strength = int_strength
            eps_list = list(eps_matrix.ravel())
            all_ads = self.adsorbate_names+self.transition_state_names
            CO_idx = all_ads.index('CO_s')
            
            ##add in third point for cross-interaction fitting 
            for k,sp in enumerate(fitting_info):
                eps, thetas, Ediffs, Eints = fitted_params[sp]
                
                if '&' in sp:
                    species = sp.split('&')                
                else:
                    species = [sp]
                if 'CO_s' in species:
                    species.remove('CO_s')
                
                E_diff_model = []
                E_int_model = []
                theta_tots = []
                for theta_CO in np.linspace(0,0.7,20):
                    theta_i = [0]*len(all_ads)
                    theta_i[CO_idx] = theta_CO
                    if species:
                        if len(species) > 1:
                            raise ValueError(species)
                        else:
                            sp_idx = all_ads.index(species[0])
                            theta_i[sp_idx] = 0.11 #very not general
                    else:
                        sp_idx = CO_idx
                    
                    theta_tot = sum(theta_i)
                    Ediff_i = self.interaction_function(theta_i,E0_list,eps_list,
                            self.interaction_response_function,False)[0][sp_idx]
                    Eint_i = self.integral_interaction_function(theta_i,E0_list,
                            eps_list,self.interaction_response_function,False)

                    E_diff_model.append(Ediff_i)
                    E_int_model.append(Eint_i)
                    theta_tots.append(theta_tot)

                ax_d = fig_diff.add_subplot(n_x,n_y,k)
                ax_i = fig_int.add_subplot(n_x,n_y,k)
                ax_d.set_title(sp)
                ax_i.set_title(sp)
                ax_d.plot(theta_tots,E_diff_model,'-k')
                ax_i.plot(theta_tots,E_int_model,'-k')
                for theta_j, Edi, Eii in zip(thetas,Ediffs,Eints):
                    ax_i.plot(sum(theta_j),Eii,'or')
                    ax_d.plot(sum(theta_j),Edi,'or')
                
            plt.show()    



    def get_interaction_info(self):
        interaction_dict = {}
        n_ads = len(self.adsorbate_names)
        cross_term_names = []
        for a in self.adsorbate_names:
            interaction_dict[a] = self.species_definitions[a].get('self_interaction_parameter',None)
            if not interaction_dict[a]:
                eps_ii0 = self.default_self_interaction_parameter
                if not getattr(self,'_self_interaction_warned',None):
                    print('Warning: No self-interaction parameter specified for '+a+'. ' + \
                            'Assuming default ('+str(eps_ii0)+')')
                    self._self_interaction_warned = True
                interaction_dict[a] = [eps_ii0]*len(self.surface_names)
#                self.species_definitions[a]['self_interaction_parameter'] = interaction_dict[a]
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
                ads_a, ads_b = sorted([a,cp]) #sort to avoid duplicates 
                name = '&'.join([ads_a,ads_b])
                if name not in cross_term_names:
                    cross_term_names.append(name)
                    interaction_dict[name] = params

        if cross_term_names:
            cross_term_names = tuple(cross_term_names) 
            self.interaction_cross_term_names = cross_term_names
    
        return interaction_dict

    def get_interaction_scaling_matrix(self):
        interaction_dict = self.get_interaction_info()

        cross_names = self.interaction_cross_term_names
        if cross_names:
            param_names = self.adsorbate_names + cross_names
        else:
            param_names = self.adsorbate_names
        constraint_dict = {}
        if not self.interaction_scaling_constraint_dict:
            self.interaction_scaling_constraint_dict = self.scaling_constraint_dict

        for ads in self.interaction_scaling_constraint_dict:
            if '-' not in ads:
                constr = self.interaction_scaling_constraint_dict[ads]
                new_constr = []
                #preserve only 0 constraints
                for ci in constr:
                    if ci != 0:
                        new_constr.append(None)
                    else:
                        new_constr.append(0)
                constraint_dict[ads] = new_constr

        for ads in self.interaction_scaling_constraint_dict:
            if '&' in ads:
                a,b = ads.split('&')
                a,b = sorted([a,b])
                new_ads = '&'.join([a,b])
                constraint_dict[new_ads] = self.interaction_scaling_constraint_dict[ads]
            else:
                constraint_dict[ads] = self.interaction_scaling_constraint_dict[ads]

        #get mins/maxs
        interaction_mins = []
        interaction_maxs = []
        if self.default_interaction_constraints is None:
            self.default_interaction_constraints = [None]*(len(self.descriptor_names)+1)
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

        return self.interaction_coefficient_matrix

    def get_interaction_transition_state_scaling_matrix(self):
        #get TS scaling matrix equivalent
        if self.transition_state_cross_interaction_mode == 'transition_state_scaling':
            if self.transition_state_scaling_matrix is None:
                raise AttributeError('Transition state scaling can only be used '+\
                ' for interactions if the transition_state_scaling_matrix is defined.')
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

    def get_interaction_matrix(self,descriptors):
        full_descriptors = list(descriptors) + [1.]
        param_vector = np.dot(self.coefficient_matrix,full_descriptors)
        eps_vector = self.params_to_matrix(param_vector)
        return eps_vector

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
        """Helper function to get `weights' of how
        to distribute TS-cross interactions between IS/FS.
        Should not be called externally."""

        weight_matrix = []
        for TS in self.transition_state_names:
            weight_TS = [0]*len(self.adsorbate_names)
            IS = None
            FS = None
            for rxn in self.elementary_rxns:
                if TS in rxn[1]:
                    if IS or FS:
                        print('Warning: Ambiguous initial/final state for '+TS)
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
        return smooth_piecewise_linear(*args,**kwargs)[:2]

    @staticmethod
    def piecewise_linear_response(*args,**kwargs):
        kwargs['smoothing'] = 0
        return smooth_piecewise_linear(*args,**kwargs)[:2]

    @staticmethod
    def linear_response(*args,**kwargs):
        kwargs['smoothing'] = 0
        kwargs['cutoff'] = 0
        return smooth_piecewise_linear(*args,**kwargs)[:2]
