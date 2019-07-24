from .scaler_base import *
import numpy as np

class ConcentrationScaler(ScalerBase):
    """Scaler which uses concentrations (fractions) as descriptors and 
    treates energetics as a constant"""
    def __init__(self,reaction_model=None):
        ScalerBase.__init__(self,reaction_model)       

    def get_electronic_energies(self,descriptors):
        if len(self.surface_names) > 1: 
            raise IndexError('Concentration scaler works only with a \
                    single surface.')

        if self.adsorbate_interaction_model not in [None,'ideal']:
            if not getattr(self.thermodynamics.adsorbate_interactions,'_parameterized',None):
                self.thermodynamics.adsorbate_interactions.parameterize_interactions() 

        energy_dict = {}
        for species in self.adsorbate_names+self.transition_state_names:
            energy_dict[species] = self.species_definitions[species]['formation_energy'][0]
        for species in self.gas_names+self.site_names:
            energy_dict[species] = self.species_definitions[species]['formation_energy']
        return energy_dict

    def get_thermodynamic_energies(self,descriptors,**kwargs):
        thermo_state = {}
        #synchronize all thermodynamic varibles
        for var,val in zip(self.descriptor_names,descriptors):
            thermo_state[var] = val
            setattr(self,var,val)
        
        if self.pressure_mode == 'static':
            #static pressure doesn't make sense if
            #pressure is a descriptor
            self.pressure_mode = 'concentration'

        assert(hasattr(self,'pressure'))

        thermo_dict =  self.thermodynamics.get_thermodynamic_corrections(
                **kwargs)

        for key in self.site_names:
            if key not in thermo_dict:
                thermo_dict[key] = 0
        return thermo_dict

    def get_rxn_parameters(self, descriptors, *args, **kwargs):
        self.force_recalculation = False
        self.thermodynamics.concentration_pressure()
        assert(all([_ in self.gas_names for _ in self.descriptor_names[0].split(':')]))
        self.balance_sps = set(self.gas_names)^set(self.descriptor_names[0].split(':'))
        self.scale_sps = self.descriptor_names[0].split(':')
        self._gas_pressures = [0. for _ in range(len(self.gas_names))]
        if len(self.scale_sps) != 2 and \
            self.descriptor_names[1] not in self.gas_names:
                raise NotImplementedError('Either ther first parameter must be a collon-divided \
                                          ratio or both the first and second must be gas species.')
        elif len(self.scale_sps) != 2 and \
            self.descriptor_names[1] in self.gas_names:
            if sum([max(_) for _ in self.parameter_names]):
                raise OverflowError('Maximum sum of concentrations must be one.')
            else:
                self.scale_sps += [self.descriptor_names[1]]
                self.balance_sps ^= {self.descriptor_names[1]}
                total_balance_concentration = sum([self.species_definitions[_]['concentration'] for _ in self.balance_sps])
                for _1, _2 in enumerate(self.gas_names):
                    if _2 in self.balance_sps:
                        self._gas_pressures[_1] = (self.species_definitions[_2]['concentration']\
                                           /total_balance_concentration)*self.pressure*(1.-sum(descriptors))
                    else:
                        if _2 == self.scale_sps[0]:
                            self._gas_pressures[_1] = descriptors[0]*self.pressure
                        else:
                            self._gas_pressures[_1] = descriptors[1]*self.pressure
        else:
            balance_frac = sum([self.species_definitions[_]['concentration'] for _ in self.balance_sps])
            for _1, _2 in enumerate(self.gas_names):
                if _2 in self.balance_sps:
                    self._gas_pressures[_1] = self.species_definitions[_2]['concentration']*self.pressure*descriptors[1]/balance_frac
                else:
                    if _2 == self.scale_sps[0]:
                        self._gas_pressures[_1] = (1.-descriptors[1])*self.pressure* \
                                            (descriptors[0]/(1.+descriptors[0]))                          
                    else:
                        self._gas_pressures[_1] = (1.-descriptors[1])*self.pressure* \
                                            (1./(1.+descriptors[0]))

        if self.adsorbate_interaction_model not in ['ideal',None]:
            params =  self.get_formation_energy_interaction_parameters()
            return params
        else:
            params = self.get_formation_energy_parameters(descriptors)
            return params

    def get_formation_energy_parameters(self,descriptors):
        self.parameter_names = self.adsorbate_names + self.transition_state_names
        free_energy_dict = self.get_free_energies(descriptors)
        params =  [free_energy_dict[sp] for sp in self.adsorbate_names+self.transition_state_names]
        return params

    def get_formation_energy_interaction_parameters(self):
        E_f = self.get_formation_energy_parameters()
        if self.interaction_cross_term_names:
            param_names = self.adsorbate_names + self.interaction_cross_term_names
        else:
            param_names = self.adsorbate_names
        
        if not self.interaction_parameters:
            info = self.thermodynamics.adsorbate_interactions.get_interaction_info()
            params = [info[pi][0] for pi in param_names]
            params_valid = []
            for p,pname in zip(params,param_names):
                if p is not None:
                    params_valid.append(p)
                else:
                    raise ValueError('No interaction parameter specified for '+pname)
            self.interaction_parameters = params_valid

        epsilon = self.thermodynamics.adsorbate_interactions.params_to_matrix(E_f+self.interaction_parameters)
        epsilon = list(epsilon.ravel())
        return E_f + epsilon
