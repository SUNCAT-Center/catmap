from .scaler_base import *

class NullScaler(ScalerBase):
    """Scaler which passes descriptor values directly to solver"""

    def get_electronic_energies(self,descriptors):
        E_dict = {}
        for g in self.gas_names:
            E_dict[g] = self.species_definitions[g]['formation_energy']

        for i,di in enumerate(self.descriptor_names):
            E_dict[di] = descriptors[i]

        return E_dict

    def get_rxn_parameters(self,descriptors):
        Gs = self.get_free_energies(descriptors)
        return [Gs[d] for d in self.descriptor_names]
#        return descriptors
