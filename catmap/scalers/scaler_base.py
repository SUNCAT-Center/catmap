import catmap
from catmap import ReactionModelWrapper
from catmap.model import ReactionModel
np = catmap.np
re = catmap.re
copy = catmap.copy
string2symbols = catmap.string2symbols

class ScalerBase(ReactionModelWrapper):
    def __init__(self,reaction_model = ReactionModel()):
        """Class for `scaling' descriptors to free energies of reaction and 
        activation (or other parameters). This class acts as a base class 
        to be inherited by other scaler classes, but is not 
        functional on its own. 

        This class contains the description of the microkinetic model 
        (adsorbate_names, gas_names, etc.) along with the temperature and 
        gas_pressures. In most cases these will automatically be populated 
        by the parent reaction_model class. 
        The scaler-specific attributes are:

        gas_energies: defines the energies of the gas-phase species. 
            This sets the references for the system.
        gas_thermo_mode: the mode for obtaining thermal contributions in 
            the gas phase. Default is to use the ideal gas approxmation.
        adsorbate_thermo_mode: the mode for obtaining thermal contributions 
            from adsorbed species. Default is to use the harmonic 
            adsorbate approximation.
        frequency_dict: a dictionary of vibrational frequencies (in eV) for 
            each gas/adsorbate. Should be of the form 
            frequency_dict[ads] = [freq1, freq2,...]. Needed for ideal gas or 
            harmonic adsorbate approximations.

        A functional derived scaler class must also contain the methods:

        get_electronic_energies(descriptors): a function to `scale' the 
            descriptors to electronic energies. Returns a dictionary of 
            the electronic energies of each species in the model.
        get_energetics(descriptors): a function to obtain the reaction 
            energetics from the descriptors. Should return a list of 
            length N (number of elementary reactions): 
            [[E_rxn1,E_a1],[E_rxn2,E_a2],...[E_rxnN,E_aN]]
        get_rxn_parameters(descriptors): a function to obtain all necessary 
            reaction parameters from the descriptors. Should return a list of 
            length N (number of elementary reactions): 
            [[param1_rxn1,param2_rxn1...]...[param1_rxnN,param2_rxnN...]]. 
            For a simple model this could be the same as get_energetics, 
            but models accounting for interactions may require more 
            parameters which can be scaled.
        """
        self._rxm = reaction_model
        defaults = dict(
		        parameter_mode = 'formation_energy',
                )

        self._rxm.update(defaults,override=False)

    def set_output_attrs(self,descriptors):
        "Function to set output information."
        ads = self.gas_names+self.adsorbate_names+self.transition_state_names
        if 'rxn_parameter' in self.output_variables:
            params = self.get_rxn_parameters(descriptors)
            self._rxn_parameter = params
            self.output_labels['rxn_parameter'] = self.parameter_names

        if 'frequency' in self.output_variables:
            self._frequency = [self.frequency_dict.get(a,[]) for a in ads]
            max_len = 0
            for f in self._frequency:
                if len(f) > max_len:
                    max_len = len(f)
            nu_labels = ['nu_'+str(i+1) for i in range(max_len)]
            self.output_labels['frequency'] = [ads,nu_labels]

        if 'electronic_energy' in self.output_variables:
            electronic_energy_dict = self.get_electronic_energies(descriptors)
            self._electronic_energy = [electronic_energy_dict[a] 
                    for a in ads]
            self.output_labels['electronic_energy'] = ads

        if 'free_energy' in self.output_variables:
            free_energy_dict = self.get_free_energies(descriptors)
            self._free_energy = [free_energy_dict[a] for a in ads]
            self.output_labels['free_energy'] = ads

        if 'gas_pressure' in self.output_variables:
            self._gas_pressure = [p for p in self.gas_pressures]
            self.output_labels['gas_pressure'] = self.gas_names

        if 'boltzmann_coverage' in self.output_variables:
            cvgs = self.thermodynamics.boltzmann_coverages(
                    self.get_free_energies(descriptors))
            self._boltzmann_coverage = cvgs
            self.output_labels['boltzmann_coverage'] = self.adsorbate_names

        if set(['enthalpy','entropy','zero_point_energy']).issubset(
                set(self.output_variables)):
            self.thermodynamics.get_thermodynamic_corrections()
            self._zero_point_energy = [self._zpe_dict[a] for a in ads]
            self._enthalpy = [self._enthalpy_dict[a] for a in ads]
            self._entropy = [self._entropy_dict[a] for a in ads]
            self.output_labels['enthalpy'] = ads
            self.output_labels['entropy'] = ads
            self.output_labels['zero_point_energy'] = ads

        if 'interaction_matrix' in self.output_variables:
            self._interaction_matrix = getattr(
                self.thermodynamics.adsorbate_interactions,
                '_interaction_matrix',
                None)
            all_names = self.adsorbate_names + self.transition_state_names
            self.output_labels['interaction_matrix'] = [all_names,all_names]

    def get_energetics(self,descriptors):
        raise AttributeError('Scaler class does not contain this method. \
                Please ensure that it is defined in order for the scaler \
                to be useful.')

    def get_rxn_parameters(self,descriptors):
        raise AttributeError('Scaler class does not contain this method. \
                Please ensure that it is defined in order for the scaler \
                to be useful.')

    def get_electronic_energies(self,descriptors):
        raise AttributeError('Scaler class does not contain this method. \
                Please ensure that it is defined in order for the scaler \
                to be useful.')

    def get_thermodynamic_energies(self,**kwargs):
        thermo_dict =  self.thermodynamics.get_thermodynamic_corrections(
                **kwargs)
        return thermo_dict

    def get_free_energies(self,descriptors,**kwargs):
        electronic_energy_dict = self.get_electronic_energies(descriptors)
        self._electronic_energy_dict = electronic_energy_dict
        thermodynamic_energy_dict = self.get_thermodynamic_energies(
                descriptors=descriptors,**kwargs)
        free_energy_dict = {}
        for key in electronic_energy_dict:
            if key in thermodynamic_energy_dict:
                E_DFT = electronic_energy_dict[key]
                G = thermodynamic_energy_dict[key]
                if G is None:
                    raise ValueError('No free energy for '+key)
                elif E_DFT is None:
                    raise ValueError('No formation energy for '+key)
                free_energy_dict[key] =  E_DFT + G
        self._gas_energies = [free_energy_dict[g] for g in self.gas_names] 
        self._site_energies = [free_energy_dict.get(s,0) for s in self.site_names] 
        return free_energy_dict

    def summary_text(self):
        return ''

