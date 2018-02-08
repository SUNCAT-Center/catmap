from .solver_base import *
from copy import copy

class IntegratedRateControlSolver(SolverBase):
    """Class for estimating rates based on the degree of rate control
    screening method {citation after published}"""

    def __init__(self,reaction_model=None):
        SolverBase.__init__(self,reaction_model)
        valid_DRCs = [sp for sp in self.species_definitions if self.species_definitions[sp].get('static_rate_control',0)]
        defaults = dict(
                include_DRCs = valid_DRCs,
                )
        self._rxm.update(defaults)
        self._required = {
                'include_DRCs':list,
                'reference_surface':str,
                'reference_rate':float,
                }


    def set_output_attrs(self,params):
        SolverBase.set_output_attrs(self,params)
        if 'integrated_DRC_rate' in self.output_variables:
            DRC_rate = self.get_integrated_DRC_rate(params)
            self._integrated_DRC_rate = DRC_rate
            self.output_labels['integrated_DRC_rate'] = ['rate']

    def compile(self):
        return

    def get_integrated_DRC_rate(self,params,*args,**kwargs):
        if not hasattr(self,'_initialized'):
            DRCs = {}
            ref_idx = self.surface_names.index(self.reference_surface)
            for dname in self.include_DRCs:
                Ef = self.species_definitions[dname]['formation_energy'][ref_idx]
                X = self.species_definitions[dname].get('static_rate_control',None)
                if X is None:
                    raise ValueError('static_rate_control not specified for: '+dname)
                else:
                    DRCs[dname] = X

            desc_energies = []
            for dname in self.descriptor_names:
                Ef = self.species_definitions[dname]['formation_energy'][ref_idx]
                desc_energies.append(Ef)
            
            ref_dict = self.scaler.get_free_energies(desc_energies)
            self.reference_free_energies = ref_dict
            self._static_rate_controls = DRCs
            self._initialized = True
        
        Gs = {}
        for key,G in zip(self.adsorbate_names + self.transition_state_names, params):
            Gs[key] = G

        dGs = {}
        for key in self._static_rate_controls:
            dGs[key] = Gs[key] - self.reference_free_energies[key]

        exp_term = 0
        for key in dGs:
            Xi = self._static_rate_controls[key]
            exp_term += Xi*(-dGs[key]/(self._kB*self.temperature))

        multiplier = np.exp(float(exp_term))
        DRC_rate = [multiplier*self.reference_rate]
        self._integrated_DRC_rate = DRC_rate
        return DRC_rate
