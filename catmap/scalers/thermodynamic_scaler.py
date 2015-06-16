from scaler_base import *
import numpy as np

class ThermodynamicScaler(ScalerBase):
    """Scaler which uses temperature/pressure/potential as descriptors and 
    treates energetics as a constant"""
    def __init__(self,reaction_model):
        ScalerBase.__init__(self,reaction_model)

    def get_electronic_energies(self,descriptors):
        if len(self.surface_names) > 1: 
            raise IndexError('Thermodynamic scaler works only with a \
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
        if 'pressure' in self.descriptor_names:
            P = thermo_state['pressure']
        elif 'logPressure' in self.descriptor_names:
            P = 10**thermo_state['logPressure']
        else:
            P = 1

        if 'pressure' in self.descriptor_names or 'logPressure' in self.descriptor_names:
            if self.pressure_mode == 'static':
                #static pressure doesn't make sense if
                #pressure is a descriptor
                self.pressure_mode = 'concentration'

        self.pressure = P

        thermo_dict =  self.thermodynamics.get_thermodynamic_corrections(
                **kwargs)

        for key in self.site_names:
            if key not in thermo_dict:
                thermo_dict[key] = 0
        return thermo_dict

    def get_rxn_parameters(self,descriptors, *args, **kwargs):
        if self.adsorbate_interaction_model not in ['ideal',None]:
            params =  self.get_formation_energy_interaction_parameters(descriptors)
            return params
        else:
            params = self.get_formation_energy_parameters(descriptors)
            return params

    def get_formation_energy_parameters(self,descriptors):
        self.parameter_names = self.adsorbate_names + self.transition_state_names
        free_energy_dict = self.get_free_energies(descriptors)
        params =  [free_energy_dict[sp] for sp in self.adsorbate_names+self.transition_state_names]
        return params

    def get_formation_energy_interaction_parameters(self,descriptors):
        E_f = self.get_formation_energy_parameters(descriptors)
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
