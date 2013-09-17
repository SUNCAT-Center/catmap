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
            self.pressure_mode = 'concentration'
        elif 'logPressure' in self.descriptor_names:
            P = 10**thermo_state['logPressure']
            self.pressure_mode = 'concentration'
        else:
            P = 1

        self.pressure = P

        thermo_dict =  self.thermodynamics.get_thermodynamic_corrections(
                **kwargs)

        for key in self.site_names:
            if key not in thermo_dict:
                thermo_dict[key] = 0
        return thermo_dict

    def get_rxn_parameters(self,descriptors):
        self.parameter_names = self.adsorbate_names + self.transition_state_names
        free_energy_dict = self.get_free_energies(descriptors)
        params =  [free_energy_dict[sp] for sp in self.adsorbate_names+self.transition_state_names]
        return params

