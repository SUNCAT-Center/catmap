import catmap
from catmap import ReactionModelWrapper
from catmap.model import ReactionModel
from catmap.functions import get_composition, add_dict_in_place
try:
    from scipy.optimize import fmin_powell
except ImportError:
    fmin_powell = None

import warnings
from mpmath import mpf
from math import exp, log
IdealGasThermo = catmap.IdealGasThermo
HarmonicThermo = catmap.HarmonicThermo
molecule = catmap.molecule
np = catmap.np
copy = catmap.copy

class ThermoCorrections(ReactionModelWrapper):
    """Class for including thermodynamic corrections.

    The function "get_thermodynamic_corrections" automatically does all the work
    assuming the correct functions are in place.

    thermodynamic_corrections: List of fundamentally different types of 
        corrections which could be included. Defaults are gas and adsorbate
        but other possibilities might be interface, electrochemical, etc.

    thermodynamic_variables: List of variables which define a thermodynamic
        state. If these attributes of the underlying reaction model do not
        change then the thermodynamic corrections will not be recalculated
        in order to save time.

    To add a new correction type (called custom_correction):
        1) Define the function which performs the correction as an attribute.
            Assume the function is called "simple_custom_correction".
        2) Place the "custom_correction" in the "thermodynamic_corrections" list
        3) Place any variables which the custom correction depends on in
            the thermodynamic_variables list
        4) Set the "custom_correction_thermo_mode" attribute of the 
            underlying reaction model to "simple_custom_correction"

    If these steps are followed then the correction should automatically be
    included in all calculations.

    """

    _kJmol2eV = 0.01036427
    _bar2Pa = 1e5

    def __init__(self,reaction_model=ReactionModel()):
        self._rxm = reaction_model
        self._log_strings = {
        'harmonic_transition_state_warning':
        'averaging initial/final state thermal contributions for ${TS}',
        'shomate_warning':
        'temperature below shomate minumum for ${gas};'+
        ' Cp(${T}) and S(${T}) are used below ${T}.'
        }

        #set defaults
        defaults = dict(
                gas_thermo_mode = 'ideal_gas',
                adsorbate_thermo_mode = 'harmonic_adsorbate',
                electrochemical_thermo_mode = 'simple_electrochemical',
                pressure_mode = 'static',
                thermodynamic_corrections = ['gas','adsorbate'],
                thermodynamic_variables = ['temperature','gas_pressures','voltage','beta','pH'],
                ideal_gas_params = catmap.data.ideal_gas_params,
                fixed_entropy_dict = catmap.data.fixed_entropy_dict,
                shomate_params = catmap.data.shomate_params,
                hbond_dict = catmap.data.hbond_dict,
                atoms_dict = {},
                frequency_dict = {},
                force_recalculation = False,
                )
        self._required = {'thermodynamic_corrections':list,
                'thermodynamic_variables':list,
                }

        self._zpe_dict = {}
        self._enthalpy_dict = {}
        self._entropy_dict = {}
        self._rxm.update(defaults)
        for corr in self.thermodynamic_corrections:
            self._required[corr+'_thermo_mode'] = str
            self.thermodynamic_variables.append(corr+'_thermo_mode')

    def get_thermodynamic_corrections(self,**kwargs):
        """
        Calculate all ``thermodynamic'' corrections beyond the energies
        in the input file. This master function will call sub-functions
        depending on the ``thermo mode'' of each class of species
        """
        l = self.thermodynamic_corrections
        if 'electrochemical' in l:
            self.force_recalculation = True
        state_dict = {}
        for v in self.thermodynamic_variables:
            state_dict[v] = getattr(self,v)
        for key in kwargs:
            if key in state_dict:
                state_dict[key] = kwargs[key]
        current_state = [repr(state_dict[v]) 
                for v in self.thermodynamic_variables]

        for sp in self.species_definitions:
            self.frequency_dict[sp] = \
                    self.species_definitions[sp].get('frequencies',[])
        frequency_dict = self.frequency_dict.copy()

        if (
                getattr(self,'_current_state',None) == current_state and 
                getattr(self,'_frequency_dict',None) == frequency_dict and
                not self.force_recalculation
                ): #if the thermodynamic state (and frequencies) 
            #has not changed then don't bother re-calculating the corrections.
            return self._correction_dict

        correction_dict = {}
        self._correction_dict = correction_dict
        self._current_state = current_state
        self._frequency_dict = frequency_dict

        # apply corrections in self.thermodynamic_corrections on top of each other
        for correction in self.thermodynamic_corrections:
            mode = getattr(self,correction+'_thermo_mode')
            thermo_dict = getattr(self,mode)()
            add_dict_in_place(correction_dict, thermo_dict)

        if self.pressure_mode:
            getattr(self,self.pressure_mode+'_pressure')()

        if 'electrochemical' in l:
            correction_dict = self._get_echem_corrections(correction_dict)
        return correction_dict

    def _get_echem_corrections(self, correction_dict):
        """
        Perform the thermodynamic corrections relevant to electrochemistry but
        are not specific to any particular mode.
        """
        # pH corrections to proton and hydroxide species
        if any(ads in ['ele_g', 'H_g', 'OH_g'] for ads in self.species_definitions.keys()):
            G_H2 = self._electronic_energy_dict['H2_g'] + self._correction_dict['H2_g']
            G_H = 0.5*G_H2 - .0592*self.pH
            G_H2O = self._electronic_energy_dict['H2O_g'] + self._correction_dict['H2O_g']
            H2O_index = self.gas_names.index('H2O_g')
            G_OH = G_H2O - G_H + self._kB*self.temperature*np.log(self.gas_pressures[H2O_index]*1.e14)
            correction_dict['H_g'] = G_H
            correction_dict['OH_g'] = G_OH

        # pressure corrections to species in the echem double layer based on kH
        if 'dl' in self.species_definitions.keys():
            dl_species = [spec for spec in self.species_definitions.keys()
                            if '_dl' in spec and '*' not in spec]
            for spec in dl_species:
                tempname = spec.split('_')[0]
                gas_spec = tempname+'_g'
                C_H2O = 55.
                KH_gas = self.species_definitions[spec].get('kH', 55.)  # defaults to no correction
                P_gas = C_H2O / KH_gas
                P_corr = np.log(P_gas) * self._kB * self.temperature
                correction_dict[spec] = correction_dict[gas_spec] + P_corr

        # Generate energy for fake echem transition states after all other corrections
        if len(self.echem_transition_state_names) > 0:
            echem_thermo_dict = self.generate_echem_TS_energies()
            add_dict_in_place(correction_dict, echem_thermo_dict)

        return correction_dict

    def ideal_gas(self):
        """Calculate the thermal correction to the free energy of 
        an ideal gas using the IdealGasThermo class in ase.thermochemistry 
        along with the molecular structures in ase.data.molecules.

        gas_names = the chemical formulas of the gasses of interest (usually 
            ending in _g to denote that they are in the gas phase).
        freq_dict = dictionary of vibrational frequencies for each gas 
            of interest. Vibrational frequencies should be in eV. 
            The dictionary should be of the form 
            freq_dict[gas_name] = [freq1, freq2, ...]
        ideal_gas_params = dictionary of the symetry number, 
            geometry keyword, and spin of the gas. If no dictionary 
            is specified then the function will attempt to look the 
            gas up in the hard-coded gas_params dictionary. 
            The dictionary should be of the form 
            ideal_gas_params[gas_name] = [symmetry_number,geometry, spin]
        atoms_dict = dictionary of ase atoms objects to use for 
            calculating rotational contributions. If none is specified 
            then the function will look in ase.data.molecules.

        """

        freq_dict = self.frequency_dict
        gas_param_dict =self.ideal_gas_params
        temperature= float(self.temperature)
        gas_names = self.gas_names

        thermo_dict = {}
        if temperature == 0: temperature = 1e-99

        gas_renames = {'CH2O_g':'H2CO_g'}

        ase_atoms_dict = {}
        for gas in self.gas_names:
            if gas in gas_renames:
                atom_name = gas_renames[gas].replace('_g','')
            else:
                atom_name = gas.replace('_g','')
            try:
                ase_atoms_dict[gas] = molecule(atom_name)
            except NotImplementedError:
                pass

        ase_atoms_dict.update(self.atoms_dict)
        self.atoms_dict = ase_atoms_dict
        atoms_dict = self.atoms_dict

        for gas in gas_names:
            # Hard coding corrections for fictitious gas molecules used in echem
            if gas in ['pe_g', 'ele_g', 'H_g', 'OH_g',
                        self.proton_species, self.hydroxide_species]:
                thermo_dict[gas] = 0.
                self._zpe_dict[gas] = 0.
                self._enthalpy_dict[gas] = 0.
                self._entropy_dict[gas] = 0.
                continue
            gpars = gas_param_dict[gas]
            symmetry,geometry,spin = gpars[:3]
            fugacity = self._bar2Pa
            #Pressures should not be used in microkinetic 
            #modelling; they are implicitly included in the 
            #rate expressions via the thermodynamic derivations.

            atoms = atoms_dict[gas]
            therm = IdealGasThermo(
                    freq_dict[gas], geometry, 
                    atoms=atoms, symmetrynumber=symmetry, 
                    spin=spin)

            ZPE = 0.5*(sum(freq_dict[gas]))

            H = therm.get_enthalpy(temperature, verbose=False)
            S = therm.get_entropy(
                    temperature, fugacity, verbose=False)

            free_energy = H-temperature*S
            thermo_dict[gas] = free_energy #use thermodynamic state 
                    #from ase.thermochemistry to calculate thermal corrections.
            self._zpe_dict[gas] = ZPE
            self._enthalpy_dict[gas] = H
            self._entropy_dict[gas] = S

        return thermo_dict

    def shomate_gas(self):
        """
        Calculate free energy corrections using shomate equation
        """
        gas_names = self.gas_names
        temperature = float(self.temperature)
        temperature_ref = 298.15

        shomate_params = self.shomate_params

        def H(T,params):
            A,B,C,D,E,F,G,H = params
            t = T/1000.0
            H = A*t + (B/2.0)*t**2 + (C/3.0)*t**3 + (D/4.0)*t**4 - E/t + F - H 
            #kJ/mol
            return H

        def S(T,params):
            A,B,C,D,E,F,G,H = params
            t = T/1000.0
            S = A*np.log(t) + B*t + (C/2.0)*t**2 + (D/3.0)*t**3 - E/(2.0*t**2) \
                    + G #J/mol*K
            return S

        def Cp(T,params):
            A,B,C,D,E,F,G,H = params
            t = T/1000.0
            Cp = A + B*t + C*t**2 + D*t**3 +E/(t**2)
            return Cp
        
        thermo_dict = {}
        for gas in gas_names:
            for key in shomate_params.keys():
                gas_key,T_range = key.split(':')
                T_min,T_max = [float(t) for t in T_range.split('-')]
                if (gas == gas_key 
                        and temperature >= T_min 
                        and temperature <= T_max
                        ):
                    params = shomate_params[key]
                    Cp_ref = Cp(temperature_ref,params)
                    dH = H(temperature,params) - H(temperature_ref,params)
                    #deltaH(298-T) = shomate(T) - shomate(298)
                    dS = S(temperature,params)
                    dH = (temperature_ref*Cp_ref/1000.0 + dH)*(self._kJmol2eV) #eV
                    #dH = 298*Cp(298) + dH(298-T)
                    dS = dS*(self._kJmol2eV/1e3) #eV/K
                    ZPE = sum(self.frequency_dict[gas])/2.0 
                    free_energy = ZPE +  dH - temperature*dS
                    self._zpe_dict[gas] = ZPE
                    self._enthalpy_dict[gas] = dH
                    self._entropy_dict[gas] = dS
                    thermo_dict[gas] = free_energy
                elif temperature < T_min and T_min < 300:
                    params = shomate_params[key]
                    Cp_ref = Cp(T_min,params)
                    dS = S(T_min,params)
                    dH = (temperature*Cp_ref/1000.0)*(self._kJmol2eV) #eV
                    dS = dS*(self._kJmol2eV/1e3) #eV/K
                    ZPE = sum(self.frequency_dict[gas])/2.0 
                    free_energy = ZPE +  dH - temperature*dS
                    self._zpe_dict[gas] = ZPE
                    self._enthalpy_dict[gas] = dH
                    self._entropy_dict[gas] = dS
                    thermo_dict[gas] = free_energy
                    self.log('shomate_warning',gas=gas,T=T_min)
        for key in gas_names:
            not_there = []
            if key not in thermo_dict:
                not_there.append(key)
            if not_there:
                raise ValueError('No Shomate parameters specified for '+' '.join(not_there))

        return thermo_dict

    def fixed_entropy_gas(self,include_ZPE=True):
        """
        Add entropy based on fixed_entropy_dict (entropy contribution to free energy assumed linear with temperature) and ZPE 
        """
        thermo_dict = {}
        gas_names = self.gas_names
        temperature = self.temperature
        entropy_dict = self.fixed_entropy_dict
        if temperature == 0: temperature = 1e-99

        freq_dict = self.frequency_dict

        for gas in gas_names:
            if include_ZPE == True:
                ZPE = 0.5*sum(freq_dict[gas])
            else:
                ZPE = 0
            if gas in entropy_dict.keys():
                S = entropy_dict[gas]
            else:
                S = entropy_dict['other']
            free_energy = ZPE-temperature*S
            thermo_dict[gas] = free_energy
            self._zpe_dict[gas] = ZPE
            self._enthalpy_dict[gas] = 0
            self._entropy_dict[gas] = S
        return thermo_dict

    def frozen_fixed_entropy_gas(self):
        """
        Do not add ZPE, calculate fixed entropy correction.
        """
        return self.fixed_entropy_gas(False)

    def zero_point_gas(self):
        """
        Add zero point energy correction to gasses.
        """
        gas_names = self.gas_names
        freq_dict = self.frequency_dict
        thermo_dict = {}
        for gas in gas_names:
            ZPE = 0.5*sum(freq_dict[gas])
            self._zpe_dict[gas] = ZPE
            self._enthalpy_dict[gas] = 0
            self._entropy_dict[gas] = 0
            thermo_dict[gas] = ZPE
        return thermo_dict

    def frozen_gas(self):
        """
        Neglect all thermal contributions, including the zero point energy.
        """
        gas_names = self.gas_names
        thermo_dict = {}
        for gas in gas_names:
            self._zpe_dict[gas] = 0
            self._enthalpy_dict[gas] = 0
            self._entropy_dict[gas] = 0
            thermo_dict[gas] = 0
        return thermo_dict

    def fixed_enthalpy_entropy_gas(self,gas_names=None):
        """
        Calculate free energy corrections based on input enthalpy, entropy, ZPE
        """
        thermo_dict = {}
        if not gas_names:
            gas_names = self.gas_names
        for gas in gas_names:
            G = 0
            species_def = self.species_definitions[gas]
            for key in ['zero_point_energy','enthalpy','entropy']:
                if key not in species_def:
                    print('Warning: No '+key+' found for '+gas+'. Using 0')
            ZPE = species_def.get('zero_point_energy',0)
            enthalpy = species_def.get('enthalpy',0)
            entropy = species_def.get('entropy',0)
            self._zpe_dict[gas] = ZPE
            self._enthalpy_dict[gas] = enthalpy
            self._entropy_dict[gas] = entropy
            thermo_dict[gas] = ZPE + enthalpy - self.temperature*entropy
        return thermo_dict

    def harmonic_adsorbate(self):
        """Calculate the thermal correction to the free energy of 
        an adsorbate in the harmonic approximation using the HarmonicThermo 
        class in ase.thermochemistry.

        adsorbate_names = the chemical formulas of the gasses of interest 
            (usually ending in _g to denote that they are in the gas phase).
            freq_dict = dictionary of vibrational frequencies for each gas of 
            interest. Vibrational frequencies should be in eV. The dictionary 
            should be of the form freq_dict[gas_name] = [freq1, freq2, ...]
        """
        adsorbate_names = self.adsorbate_names+self.transition_state_names
        temperature = float(self.temperature)
        freq_dict = self.frequency_dict

        thermo_dict = {}
        if temperature == 0: temperature = 1e-99

        avg_TS = []
        self._freq_cutoffs = {}

        for ads in adsorbate_names:
            if ads in freq_dict:
                if '-' in ads and freq_dict[ads] in [None,[],()]:
                    avg_TS.append(ads)
                frequencies = freq_dict[ads]
                if self.max_entropy_per_mode:
                    if temperature in self._freq_cutoffs:
                        nu_min = self._freq_cutoffs[temperature]
                    else:
                        kB_multiplier = float(self.max_entropy_per_mode/self._kB)
                        nu_min = self.get_frequency_cutoff(kB_multiplier,float(temperature))
                        nu_min /= 1000.
                        self._freq_cutoffs[temperature] = nu_min

                    frequencies = [max(nu,nu_min) for nu in frequencies]
                therm = HarmonicThermo(frequencies)
                try:
                    free_energy = therm.get_helmholtz_energy(
                            temperature,verbose=False)
                except AttributeError:
                    warnings.warn('HarmonicThermo.get_free_energy is deprecated.'
                                   'Update your ASE version.')
                    free_energy = therm.get_free_energy(
                            temperature,verbose=False)
                ZPE = sum(frequencies)/2.0 
                dS = therm.get_entropy(temperature,verbose=False)
                dH = therm.get_internal_energy(temperature,verbose=False) - ZPE
                self._zpe_dict[ads] = ZPE
                self._enthalpy_dict[ads] = dH
                self._entropy_dict[ads] = dS
                thermo_dict[ads] = free_energy #use thermodynamic state from 
                #ase.thermochemistry to calculate thermal corrections.
            elif '-' in ads:
                avg_TS.append(ads)
            else:
                raise IndexError('Missing vibrational frequencies for '+ads)

        ts_thermo = self.average_transition_state(thermo_dict,avg_TS)
        thermo_dict.update(ts_thermo)

        return thermo_dict
    
    def zero_point_adsorbate(self):
        """
        Add zero point energy correction to adsorbate energy.
        """
        adsorbate_names = self.adsorbate_names+self.transition_state_names
        freq_dict = self.frequency_dict
        thermo_dict = {}
        avg_TS = []
        for ads in adsorbate_names:
            if freq_dict.get(ads,None):
                ZPE = 0.5*sum(freq_dict[ads])
                self._zpe_dict[ads] = ZPE
                self._enthalpy_dict[ads] = 0
                self._entropy_dict[ads] = 0
                thermo_dict[ads] = ZPE
            elif '-' in ads:
                avg_TS.append(ads)
            else:
                raise IndexError('Missing vibrational frequencies for '+ads)

        ts_thermo = self.average_transition_state(thermo_dict,avg_TS)
        thermo_dict.update(ts_thermo)

        return thermo_dict

    def frozen_adsorbate(self):
        """
        Neglect all zero point, enthalpy, entropy corrections to adsorbate energy.
        """
        thermo_dict = {}
        for ads in self.adsorbate_names+self.transition_state_names:
            self._zpe_dict[ads] = 0
            self._enthalpy_dict[ads] = 0
            self._entropy_dict[ads] = 0
            thermo_dict[ads] = 0
        return thermo_dict

    def fixed_enthalpy_entropy_adsorbate(self):
        """
        Return free energy corrections based on input enthalpy, entropy, ZPE
        """
        return self.fixed_enthalpy_entropy_gas(self.adsorbate_names+self.transition_state_names)

    def average_transition_state(self,thermo_dict,transition_state_list = []):
        """
        Return transition state thermochemical corrections as average of IS and FS corrections 
        """
        if transition_state_list is None:
            transition_state_list = self.transition_state_names

        def state_thermo(therm_dict,rx,site_defs,rx_id):
            return sum([therm_dict[s] for s in rx[rx_id] if (
                            s not in site_defs and not 
                            s.endswith('_g'))])

        for ads in transition_state_list:
            self.log('harmonic_transition_state_warning',TS=ads)
            rx = [rx for rx in self.elementary_rxns if ads in rx[1]][0]
            for therm_dict in [thermo_dict,self._zpe_dict,
                    self._enthalpy_dict,self._entropy_dict]:
                IS = state_thermo(therm_dict,rx,self.site_names,0)
                FS = state_thermo(therm_dict,rx,self.site_names,-1)
                therm_dict[ads] = (IS+FS)/2.0
        return thermo_dict

    def generate_echem_TS_energies(self):
        """ 
        Give real energies to the fake echem transition states
        """
        echem_TS_names = self.echem_transition_state_names
        voltage = self.voltage
        beta = self.beta
        thermo_dict = {}

        def get_E_to_G(state, E_to_G_dict):
            E_to_G = 0.
            for ads in state:
                if ads in E_to_G_dict:
                    E_to_G += E_to_G_dict[ads]
            return E_to_G

        for echem_TS in echem_TS_names:
            preamble, site = echem_TS.split('_')
            echem, rxn_index, barrier = preamble.split('-')
            rxn_index = int(rxn_index)
            rxn = self.elementary_rxns[rxn_index]
            if rxn_index in self.rxn_options_dict['beta'].keys():
                beta = float(self.rxn_options_dict['beta'][rxn_index])
            IS = rxn[0]
            FS = rxn[-1]
            E_IS = self.get_state_energy(IS, self._electronic_energy_dict)
            E_FS = self.get_state_energy(FS, self._electronic_energy_dict)
            G_IS = E_IS + get_E_to_G(IS, self._correction_dict)
            G_FS = E_FS + get_E_to_G(FS, self._correction_dict)
            dG = G_FS - G_IS
            G_TS = G_FS + float(barrier) + (1 - beta) *  -dG  # only tested for reductions
            # make sure we're "correcting" the right value
            assert(self._electronic_energy_dict[echem_TS]) == 0.
            self._correction_dict[echem_TS] = 0.

            thermo_dict[echem_TS] = G_TS
        return thermo_dict

    def get_rxn_index_from_TS(self, TS):
        """ 
        Take in the name of a transition state and return the reaction index of
        the elementary rxn from which it belongs
        """
        for rxn_index, eq in enumerate(self.elementary_rxns):
            if len(eq) == 3 and TS in eq[1]:
                return rxn_index

    def simple_electrochemical(self):
        """
        Calculate electrochemical (potential) corrections to free energy. Transition state energies are corrected by a beta*voltage term.  
        """
        thermo_dict = {}
        gas_names = [gas for gas in self.gas_names if gas.split('_')[0] in ['pe', 'ele']]
        TS_names = [TS for TS in self.transition_state_names if
            'pe' in TS.split('_')[0] or 'ele' in TS.split('_')[0]]
        voltage = self.voltage
        beta = self.beta

        # scale pe thermo correction by voltage (vs RHE)
        for gas in gas_names:
            thermo_dict[gas] = -voltage

        # no hbond correction for simple_electrochemical

        # correct TS energies with beta*voltage (and hbonding?)
        for TS in TS_names:
            rxn_index = self.get_rxn_index_from_TS(TS)
            if rxn_index in self.rxn_options_dict['beta'].keys():
                beta = float(self.rxn_options_dict['beta'][rxn_index])
            thermo_dict[TS] = -voltage * (1 - beta)

        return thermo_dict

    def hbond_electrochemical(self):
        """
        Update simple_electrochemical with hbonding corrections as if they were on Pt(111)
        """
        thermo_dict = self.simple_electrochemical()
        TS_names = [TS for TS in self.transition_state_names if
            'pe' in TS.split('_')[0] or 'ele' in TS.split('_')[0]]
        hbond_dict = self.hbond_dict
        for ads in list(self.adsorbate_names) + TS_names:
            if ads in hbond_dict:
                if ads in thermo_dict:
                    thermo_dict[ads] += hbond_dict[ads]
                else:
                    thermo_dict[ads] = hbond_dict[ads]
            elif ads.split('_')[0] in hbond_dict:
                if ads in thermo_dict:
                    thermo_dict[ads] += hbond_dict[ads.split('_')[0]]
                else:
                    thermo_dict[ads] = hbond_dict[ads.split('_')[0]]

        return thermo_dict

    def estimate_hbond_corr(formula):
        """
        Generate hydrogen bonding corrections given a formula and estimations
        for various functional groups used in Peterson(2010) - valid mostly for Pt(111)
        This is a very simplistic function.  If you need more advanced descriptions of
        hydrogen bonding, consider setting your own hbond_dict.
        """
        num_OH = formula.count('OH')
        num_O = get_composition(formula.split('_s')[0]).setdefault('O',0)
        num_ketone = num_O - num_OH
        if num_ketone < 0:
            print "(number of O) - (number of OH) is negative??"
            assert(False)
        if formula in ["OH_s", "OH"]:
            corr = -0.50
        else:
            corr = -0.25 * num_OH + -0.10 * num_ketone
        print "estimated hbond correction for", formula, "is", corr
        return corr

    def hbond_with_estimates_electrochemical(self):
        """
        Add hbond corrections to transition states involving pe and ele (coupled proton-electron transfers and electron transfers)
        """
        thermo_dict = self.hbond_electrochemical()
        TS_names = [TS for TS in self.transition_state_names if
            'pe' in TS.split('_')[0] or 'ele' in TS.split('_')[0]]

        for ads in list(self.adsorbate_names) + TS_names:
            if ads not in hbond_dict:
                if ads in thermo_dict:
                    thermo_dict[ads] += estimate_hbond_corr(ads)
                else:
                    thermo_dict[ads] = estimate_hbond_corr(ads)

        return thermo_dict

    def boltzmann_coverages(self,energy_dict):
        """
        Return coverages based on Boltzmann distribution
        """
        #change the reference
        reservoirs = getattr(self,'atomic_reservoir_dict',None)
        if reservoirs:
            comp_dict = {}
            for sp in energy_dict.keys():
                comp_dict[sp] = self.species_definitions[sp]['composition']
            energy_dict = self.convert_formation_energies(
                    energy_dict,reservoirs,comp_dict)[0]

        #calculate coverages
        cvgs = [0]*len(self.adsorbate_names)
        for site in self.site_names:
            if site not in energy_dict:
                energy_dict[site] = 0
            relevant_ads = [a for a in self.adsorbate_names if 
                    self.species_definitions[a]['site'] == site]
            free_energies = [energy_dict[a] for a in relevant_ads]+[energy_dict[site]]
            boltz_sum = sum([self._math.exp(-G/(self._kB*self.temperature)) 
                for G in free_energies])
            for ads in relevant_ads:
                if ads in self.adsorbate_names:
                    i_overall = self.adsorbate_names.index(ads)
                    i_rel = relevant_ads.index(ads)
                    if self.species_definitions[site]['type'] not in ['gas']:
                        cvgs[i_overall] = self._math.exp(-free_energies[i_rel]/(
                            self._kB*self.temperature))/boltz_sum
        return cvgs

    def static_pressure(self):
        self.gas_pressures = [self.species_definitions[g]['pressure'] for g in self.gas_names]

    def concentration_pressure(self):
        if 'pressure' not in self.thermodynamic_variables:
            self.thermodynamic_variables += ['pressure']
        for g in self.gas_names:
            if 'concentration' not in self.species_definitions[g]:
                raise UserWarning('Concentration not specified for '+g+'.')
        self.gas_pressures = [self.species_definitions[g]['concentration']*self.pressure for g in self.gas_names]

    def approach_to_equilibrium_pressure(self):
        """Set product pressures based on approach to equilibrium. Requires the following attributes
        to be set:
        global_reactions - a list of global reactions in the same syntax as elementary expressions,
            with each one followed by its respective approach to equilibrium.
        pressure_mode - must be set to 'approach_to_equilibrium'
        Note that this function is not well-tested and should be used with caution.
        """

           
        if 'pressure' not in self.thermodynamic_variables:
            self.thermodynamic_variables += ['pressure']

        self.pressure_mode = None #avoid recursion...
        G_dict = self.scaler.get_free_energies(self._descriptors)
        self.pressure_mode = 'approach_to_equilibrium'

        def set_product_pressures(rxn,G_dict,gamma,gas_pressures,product_pressures={}):

            dG = mpf(self.get_rxn_energy(rxn,G_dict)[0])
            K_eq = exp(-dG/(self._kB*self.temperature))
            K_eq *= mpf(gamma)
            PK = K_eq

            for g in gas_pressures:
                n_g = mpf(rxn[0].count(g))
                p_g = mpf(gas_pressures[g])
                PK *= pow(p_g,n_g)

            #remove any products that have pressures specified
            products = rxn[-1]
            for sp in products:
                if sp in product_pressures:
                    gas_pressures[sp] = product_pressures[sp]
                    n_prod = products.count(sp)
                    PK /= pow(mpf(product_pressures[sp]),n_prod)
                    products = [p for p in products if p != sp]

            N_products = len(products)
            P_i = pow(PK,(mpf(1.)/mpf(N_products)))
            for gi in set(products):
                gas_pressures[gi] += P_i

        gas_pressures = {}
        global_rxns = []
        gammas = []
        reactants = []
        products = []
        for gi,gamma_i in self.global_reactions:
            rxn = self.expression_string_to_list(gi)
            global_rxns.append(rxn)
            reactants += rxn[0]
            products += rxn[1]
            gammas.append(gamma_i)
        
        for sp in set(reactants):
            c = self.species_definitions[sp].get('concentration',None)
            if c is None:
                raise UserWarning('Concentrations must be specified for all reactants. No concentration '
                        'was found for '+str(sp)+'.')
            gas_pressures[sp] = c*self.pressure

        product_pressures = {}
        for sp in set(products):
            gas_pressures[sp] = 0
            c = self.species_definitions[sp].get('concentration',None)
            if c is not None:
                product_pressures[sp] = c*self.pressure
                gas_pressures[sp] = c*self.pressure

        #Ensure that all gasses are either products or reactants
        if sorted(list(set(products))+list(set(reactants))) != sorted(self.gas_names):
            raise UserWarning('All gasses must appear as products or reactants in self.global_reactions.')

        for gi,gamma_i in zip(global_rxns,gammas):
            set_product_pressures(gi,G_dict,gamma_i,gas_pressures,product_pressures)

        self.gas_pressures = [gas_pressures[gi] for gi in self.gas_names] 
        print [float(i) for i in self.gas_pressures]

    def get_frequency_cutoff(self,kB_multiplier,temperature=None):
        kB = float(self._kB)
        if temperature is None:
            T = self.temperature
        else:
            T = temperature

        def get_entropy(nu,T):
            nu_eV = nu*1e-3 #input in meV
            HT = HarmonicThermo([float(nu_eV)])
            return HT.get_entropy(float(T),verbose=False)

        def target(nu,kB_multiplier=kB_multiplier,T=T):
            nu = max(1e-8,nu)
            SbykB = (get_entropy(nu,T)/kB)
            return (kB_multiplier - SbykB)**2

        nu_cutoff = fmin_powell(target,10,disp=0)
        return nu_cutoff

    def summary_text(self):
        return ''


def fit_shomate(Ts, Cps, Hs, Ss, params0,plot_file = None):
    try:
        from scipy.optimize import leastsq
    except ImportError:
        print("Scipy not installed, returning initial guess")
        leastsq = lambda resid, initial, **kwargs: [initial, True]

    def H(t,A,B,C,D,E,F,H_c):
        H = A*t + (B/2.0)*t**2 + (C/3.0)*t**3 + (D/4.0)*t**4 - E/t + F - H_c 
        #kJ/mol
        return H
    def H_resid(params,H_act,t):
        A,B,C,D,E,F,H_c = params
        return H_act - H(t,A,B,C,D,E,F,H_c)

    def S(t,A,B,C,D,E,G):
        S = A*np.log(t) + B*t + (C/2.0)*t**2 + (D/3.0)*t**3 - E/(2.0*t**2) + G 
        #J/mol*K
        return S
    def S_resid(params,S_act,t):
        A,B,C,D,E,G = params
        return S_act - S(t,A,B,C,D,E,G)

    def Cp(t,A,B,C,D,E):
        Cp = A + B*t + C*t**2 + D*t**3 +E/(t**2)
        return Cp

    def Cp_resid(params,Cp_act,t):
        A,B,C,D,E = params
        return Cp_act - Cp(t,A,B,C,D,E)

    A0,B0,C0,D0,E0,F0,G0,H0 = params0
    Cps = np.array(Cps)
    Hs = np.array(Hs)
    Ss = np.array(Ss)
    ts = np.array(Ts)
    ts = ts/1000.0

    [A,B,C,D,E,F,H_c],flag = leastsq(
            H_resid,[A0,B0,C0,D0,E0,F0,H0],args=(Hs,ts))
    [A,B,C,D,E,G],flag = leastsq(S_resid,[A,B,C,D,E,G0],args=(Ss,ts))


    if plot_file:
        import pylab as plt
        fig = plt.figure()
        ax1 = fig.add_subplot(131)
        ax2 = fig.add_subplot(132)
        ax3 = fig.add_subplot(133)
        ax1.plot(ts,H(ts,A,B,C,D,E,F,H_c),'-k')
        ax1.plot(ts,Hs,'ok')
        ax1.set_title('H')
        ax2.plot(ts,S(ts,A,B,C,D,E,G),'-r')
        ax2.plot(ts,Ss,'or')
        ax2.set_title('S')
        ax3.plot(ts,Cp(ts,A,B,C,D,E),'-b')
        ax3.plot(ts,Cps,'ob')
        ax2.set_title('Cps')
        fig.savefig(plot_file)


    return [A,B,C,D,E,F,G,H_c]
