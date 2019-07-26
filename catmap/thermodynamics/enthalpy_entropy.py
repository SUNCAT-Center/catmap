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
HinderedThermo = catmap.HinderedThermo
molecule = catmap.molecule
np = catmap.np
# copy = catmap.copy

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

    def __init__(self,reaction_model=None):
        if reaction_model is None:
            reaction_model = ReactionModel()
        self._rxm = reaction_model
        self._log_strings = {
        'harmonic_transition_state_warning':
        'averaging initial/final state thermal contributions for ${TS}',
        'shomate_warning':
        'temperature below shomate minimum for ${gas};'+
        ' Cp(${T}) and S(${T}) are used below ${T}.',
        'force_recompilation':
        'Enabling model.force_recompilation = True.  Necessary for field corrections',
        }

        #set defaults
        defaults = dict(
                gas_thermo_mode = 'ideal_gas',
                adsorbate_thermo_mode = 'harmonic_adsorbate',
                electrochemical_thermo_mode = 'simple_electrochemical',
                pressure_mode = 'static',
                thermodynamic_corrections = ['gas','adsorbate'],
                thermodynamic_variables = ['temperature','gas_pressures','voltage','beta','pH','Upzc'],
                ideal_gas_params = catmap.data.ideal_gas_params,
                hindered_ads_params = {},
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

        self.kBT = self._kB*self.temperature
        
        # apply corrections in self.thermodynamic_corrections on top of each other
        for correction in self.thermodynamic_corrections:
            mode = getattr(self,correction+'_thermo_mode')
            thermo_dict = getattr(self,mode)()
            add_dict_in_place(correction_dict, thermo_dict)

        if self.pressure_mode:
            getattr(self,self.pressure_mode+'_pressure')()
            
        #if hasattr(self,'equilibrated') and self.equilibrated:
        #    self.set_equilibrated()

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
            G_OH = G_H2O - G_H # Do not need Kw, just need to make sure equilibria are satisfied
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
        ideal_gas_params = dictionary of the symmetry number, 
            geometry keyword, and spin of the gas. If no dictionary 
            is specified then the function will attempt to look the 
            gas up in the hard-coded gas_params dictionary. 
            The dictionary should be of the form 
            ideal_gas_params[gas_name] = [symmetry_number, geometry, spin]
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
            except (NotImplementedError, KeyError):
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
            S = therm.get_entropy(temperature, fugacity, verbose=False)

            free_energy = H-temperature*S
            thermo_dict[gas] = free_energy #use thermodynamic state 
                    #from ase.thermochemistry to calculate thermal corrections.
            self._zpe_dict[gas] = ZPE
            self._enthalpy_dict[gas] = H
            self._entropy_dict[gas] = S

        return thermo_dict

    def shomate_gas(self):
        """
        Calculate free energy corrections using Shomate equation
        """
        gas_names = self.gas_names
        temperature = float(self.temperature)

        # added to avoid unnecessary inner loop ahead
        shomate_params = {}        
        for _ in self.shomate_params.keys():
            _n = _.split(':')[0].replace('*','')
            T_min, T_max = sorted([float(__) for __ in _.split(':')[1].split('-')])
            if _n in shomate_params.keys():
                shomate_params[_n] += [{'params':self.shomate_params[_],'T_min':T_min,'T_max':T_max}]
            else:
                shomate_params[_n] = [{'params':self.shomate_params[_],'T_min':T_min,'T_max':T_max}]
                
        thermo_dict = {}
        not_there = []
        for gas in gas_names:
            if gas in shomate_params.keys():
                params = shomate_params[gas]
                loc = [_ for _ in [_ if ((temperature>=params[_]['T_min'] and temperature<=params[_]['T_max']) or\
                             (temperature<params[_]['T_min'] and params[_]['T_min']<300)) else None for _ in range(len(params))] if _ is not None]
                if len(loc) > 0:
                    params = params[loc[0]]
                    if (temperature >= params['T_min'] and temperature <= params['T_max']):
                        dH, dS, Cp_ref = self._shomate_eq(params['params'])
                    elif temperature < params['T_min'] and params['T_min'] < 300:
                        dH, dS, Cp_ref = self._shomate_eq(params['params'],temperature=params['T_min'])
                        self.log('shomate_warning_gas',gas=gas,T=params['T_min'])
                    else:
                        raise Exception('Logical error.')
                    ZPE = sum(self.frequency_dict[gas])/2.0     
                    free_energy = ZPE +  dH - temperature*dS 
                    self._zpe_dict[gas] = ZPE
                    self._enthalpy_dict[gas] = dH
                    self._entropy_dict[gas] = dS
                    thermo_dict[gas] = free_energy
                else:
                    print(shomate_params.keys)
                    print(shomate_params[gas])
                    raise ValueError('No Shomate parameters available for T = {}  specified for {}'.format(temperature,gas))
            else:
                not_there.append(gas)
        if not_there:
            raise ValueError('No Shomate parameters specified for '+' '.join(not_there))

        return thermo_dict

    def fixed_entropy_gas(self,include_ZPE=True):
        """
        Add entropy based on fixed_entropy_dict (entropy contribution to free 
        energy assumed linear with temperature) and ZPE
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

        adsorbate_names = the chemical formulas of the adsorbates of interest.
        freq_dict = dictionary of vibrational frequencies for each adsorbate of 
            interest. Vibrational frequencies should be in eV. The dictionary 
            should be of the form freq_dict[ads_name] = [freq1, freq2, ...]
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
                        kB_multiplier = float(
                                self.max_entropy_per_mode/self._kB)
                        nu_min = self.get_frequency_cutoff(
                                        kB_multiplier,float(temperature))
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
    
    def _shomate_eq(self,params,temperature=[]):
        if not temperature:
            temperature = self.temperature
        temperature_ref = 298.15
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
        
        dH = H(temperature,params) - H(temperature_ref,params)
        dS = S(temperature,params)
        Cp_ref = Cp(temperature_ref,params)
        dH = (temperature_ref*Cp_ref/1000.0 + dH)*(self._kJmol2eV) #eV
        dS = dS*(self._kJmol2eV/1e3) #eV/K
        return dH, dS, Cp_ref
    
    def shomate_adsorbate(self):
        """Calculate the thermal correction to the free energy of an
        adsorbate using pre-fitted shomate parameters.        
        """
        
        adsorbate_names = self.adsorbate_names+self.transition_state_names
        temperature = float(self.temperature)
        
        thermo_dict = {}
        if temperature == 0: temperature = 1e-99

        avg_TS = []
        self._freq_cutoffs = {}
        
        # added to avoid unnecessary inner loop ahead
        shomate_params = {}        
        for _ in self.shomate_params.keys():
            _n = _.split(':')[0].replace('*','')
            T_min, T_max = sorted([float(__) for __ in _.split(':')[1].split('-')])
            if _n in shomate_params.keys():
                shomate_params[_n] += [{'params':self.shomate_params[_],'T_min':T_min,'T_max':T_max}]
            else:
                shomate_params[_n] = [{'params':self.shomate_params[_],'T_min':T_min,'T_max':T_max}]
        
        def _thermo(params):
            if (temperature >= params['T_min'] and temperature <= params['T_max']):
                dH, dS, Cp_ref = self._shomate_eq(params['params'])
            elif temperature < params['T_min'] and params['T_min'] < 300:
                dH, dS, Cp_ref = self._shomate_eq(params['params'],temperature=params['T_min'])
                self.log('shomate_warning_ads',adsorbate=ads,T=params['T_min'])
            else:
                raise Exception('Logical error.')
            return dH, dS, Cp_ref
        
        for ads in adsorbate_names:
            # this will also work for TS if shomate parameters are available.
            if ads in shomate_params:
                params = shomate_params[ads]
                loc = list(filter(None.__ne__,[_ if ((temperature>=params[_]['T_min'] and temperature<=params[_]['T_max']) or\
                             (temperature<params[_]['T_min'] and params[_]['T_min']<300)) else None for _ in range(len(params))]))
                if len(loc) > 0:
                    if self.frequency_dict[ads]:
                        ZPE = sum(self.frequency_dict[ads])/2.0 
                        self._zpe_dict[ads] = ZPE
                    else:
                        self.average_transition_state(thermo_dict,transition_state_list = [ads], thermo_vars = [self._zpe_dict])
                        ZPE = self._zpe_dict[ads]
                    params = params[loc[0]]
                    dH, dS, Cp_ref = _thermo(params)
                    self._enthalpy_dict[ads] = dH
                    self._entropy_dict[ads] = dS
                    free_energy = ZPE +  dH - temperature*dS 
                    thermo_dict[ads] = free_energy
                else:
                    raise ValueError('No Shomate parameters available for T = {}  specified for {}'.format(temperature,ads))
            elif '-' in ads:
                avg_TS.append(ads)
            else:
                raise ValueError('No Shomate parameters specified for '+' '.join(ads))
        ts_thermo = self.average_transition_state(thermo_dict,avg_TS)
        thermo_dict.update(ts_thermo)
        return thermo_dict
    
    def hindered_adsorbate(self):
        """Calculate the thermal correction to the free energy of an 
        adsorbate in the hindered translator and hindered rotor model using 
        the HinderedThermo class in ase.thermochemistry along with the 
        molecular structures in ase.data.molecules. Requires ase version 
        3.12.0 or greater.

        adsorbate_names = the chemical formulas of the adsorbates of interest. 
        freq_dict = dictionary of vibrational frequencies for each adsorbate of 
            interest. Vibrational frequencies should be in eV. The dictionary 
            should be of the form freq_dict[ads_name] = [freq1, freq2, ...]
        hindered_ads_params = dictionary containing for each adsorbate 
            [0] = translational energy barrier in eV (barrier for the 
                  adsorbate to diffuse on the surface)
            [1] = rotational energy barrier in eV (barrier for the adsorbate 
                  to rotate about an axis perpendicular to the surface)
            [2] = surface site density in cm^-2 
            [3] = number of equivalent minima in full adsorbate rotation 
            [4] = mass of the adsorbate in amu (can be unspecified by putting 
                  None, in which case mass will attempt to be calculated from 
                  the ase atoms class)
            [5] = reduced moment of inertia of the adsorbate in amu*Ang^-2 
                  (can be unspecified by putting None, in which case inertia 
                  will attempt to be calculated from the ase atoms class)
            [6] = symmetry number of the adsorbate (number of symmetric arms 
                  of the adsorbate which depends upon how it is bound to the 
                  surface. For example, propane bound through its end carbon 
                  has a symmetry number of 1 but propane bound through its 
                  middle carbon has a symmetry number of 2. For single atom 
                  adsorbates such as O* the symmetry number is 1.)
            The dictionary should be of the form 
            hindered_ads_params[ads_name] = [barrierT, barrierR, site_density, 
            rotational_minima, mass, inertia, symmetry_number]
        atoms_dict = dictionary of ase atoms objects to use for calculating 
            mass and rotational inertia. If none is specified then the function 
            will look in ase.data.molecules. Can be omitted if both mass and 
            rotational inertia are specified in hindered_ads_params.

        """
        adsorbate_names = self.adsorbate_names+self.transition_state_names
        temperature = float(self.temperature)
        freq_dict = self.frequency_dict
        ads_param_dict = self.hindered_ads_params

        thermo_dict = {}
        if temperature == 0: temperature = 1e-99

        ase_atoms_dict = {}
        for ads in self.adsorbate_names:
            atom_name = ads.rsplit('_',1)[0]
            try:
                ase_atoms_dict[ads] = molecule(atom_name)
            except(NotImplementedError,KeyError):
                pass
        ase_atoms_dict.update(self.atoms_dict)
        self.atoms_dict = ase_atoms_dict
        atoms_dict = self.atoms_dict

        avg_TS = []
        self._freq_cutoffs = {}

        for ads in adsorbate_names:
            if '-' in ads and (freq_dict[ads] in [None,[],()] or 
                    ads not in ads_param_dict):
                avg_TS.append(ads)
                break

            #get frequencies or throw error
            if freq_dict[ads] not in [None,[],()]:
                frequencies = freq_dict[ads]
            else:
                raise IndexError('Missing vibrational frequencies for '+ads)
            if self.max_entropy_per_mode:
                if temperature in self._freq_cutoffs:
                    nu_min = self._freq_cutoffs[temperature]
                else:
                    kB_multiplier = float(
                            self.max_entropy_per_mode/self._kB)
                    nu_min = self.get_frequency_cutoff(
                            kB_multiplier,float(temperature))
                    nu_min /= 1000.
                    self._freq_cutoffs[temperature] = nu_min
                frequencies = [max(nu,nu_min) for nu in frequencies]

            #get all other parameters or throw error
            if ads in ads_param_dict:
                apars = ads_param_dict[ads]
            else:
                raise IndexError('Missing hindered_ads_params for '+ads)
            barrierT = apars[0]
            barrierR = apars[1]
            sitedensity = apars[2]
            rotationalminima = apars[3]
            mass = apars[4]
            inertia = apars[5]
            symmetrynumber = apars[6]
            try:
                atoms = atoms_dict[ads]
            except:
                atoms = {}
            if not ((mass and inertia) or atoms):
                if '-' in ads:
                    avg_TS.append(ads)
                    break
                else:
                    raise IndexError('Missing either mass and inertia of '+ads+
                                     ' or atoms object for '+ads)
            therm = HinderedThermo(
                    frequencies, barrierT, barrierR, sitedensity, 
                    rotationalminima, mass=mass, inertia=inertia, 
                    atoms=atoms, symmetrynumber=symmetrynumber)

            free_energy = therm.get_helmholtz_energy(
                    temperature,verbose=False)
            ZPE = therm.get_zero_point_energy(verbose=False)
            dS = therm.get_entropy(temperature,verbose=False)
            dH = therm.get_internal_energy(temperature,verbose=False) - ZPE
            self._zpe_dict[ads] = ZPE
            self._enthalpy_dict[ads] = dH
            self._entropy_dict[ads] = dS
            thermo_dict[ads] = free_energy #use thermodynamic state from 
            #ase.thermochemistry to calculate thermal corrections.

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

    def average_transition_state(self,thermo_dict,transition_state_list = [], thermo_vars = []):
        """
        Return transition state thermochemical corrections as average of IS and FS corrections 
        """
        if transition_state_list is None:
            transition_state_list = self.transition_state_names

        def state_thermo(therm_dict,rx,site_defs,rx_id):
            return sum([therm_dict[s] for s in rx[rx_id] if (
                            s not in site_defs and not 
                            s.endswith('_g'))])
        
        if thermo_vars:
            if thermo_dict not in thermo_vars:
                thermo_vars = [thermo_dict] + thermo_vars
            else:
                pass
        else:
            thermo_vars = [thermo_dict,self._zpe_dict,
                    self._enthalpy_dict,self._entropy_dict]

        for ads in transition_state_list:
            self.log('harmonic_transition_state_warning',TS=ads)
            rx = [rx for rx in self.elementary_rxns if ads in rx[1]][0]
            for therm_dict in thermo_vars:
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

    def homogeneous_field(self):
        """
        Update simple_electrochemical with field corrections for adsorbates that respond to a field
        """
        thermo_dict = self.simple_electrochemical()
        voltage = self.voltage
        Upzc = self.Upzc
        #distance = self.distance
        #field = (voltage - Upzc)/distance
        field = self.field
        for ads in self.adsorbate_names + self.transition_state_names:
            if 'field_params' in self.species_definitions[ads]:
                mu,alpha = self.species_definitions[ads]['field_params']
                if ads in thermo_dict:
                    thermo_dict[ads] += mu*field + alpha*(field)**2
                else:
                    thermo_dict[ads] = mu*field + alpha*(field)**2

        return thermo_dict

    def local_field_electrochemical(self):
        """
        Obtains corrections to thermo_dict in the presence of ions
        hey you need to specify these things:
        model.Upzc (float)
        model.CH (float)
        model.field_site_name
        model.unfield_site_name
        and DO NOT specify beta
        """
        if not self.force_recompilation:
                self.force_recompilation = True
                self.log('force_recompilation')
        thermo_dict = self.simple_electrochemical()
        voltage = self.voltage
        Upzc = self.Upzc
        CH = self.CH
        field = self.field      #Typically, the local field exerted by a cation is on the order of -1.0 V/Angstrom
        theta_ion = -CH * (voltage - Upzc)
        if theta_ion < 0:
            theta_ion = 0.0
        if theta_ion > 1:
            theta_ion = 1.0
        #self.rxn_options_dict['beta'] = {}
        self.species_definitions['b']['total'] = theta_ion
        self.species_definitions['a']['total'] = 1 - theta_ion
        for ads in self.adsorbate_names + self.transition_state_names:
            if 'field_params' in self.species_definitions[ads]:
                mu,alpha = self.species_definitions[ads]['field_params']
                if ads in thermo_dict:
                    thermo_dict[ads] += mu*field - (alpha/2.0)*(field)**2
                else:
                    thermo_dict[ads] = mu*field - (alpha/2.0)*(field)**2

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
            print("(number of O) - (number of OH) is negative??")
            assert(False)
        if formula in ["OH_s", "OH"]:
            corr = -0.50
        else:
            corr = -0.25 * num_OH + -0.10 * num_ketone
        print("estimated hbond correction for {formula:s} is {corr:s}".format(formula=formula, corr=corr))
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
        if hasattr(self,'_gas_pressures'):
            self.gas_pressures = self._gas_pressures
        else:
            self.gas_pressures = [self.species_definitions[g]['concentration']*self.pressure for g in self.gas_names]

    def approach_to_equilibrium_pressure(self):
        """Set product pressures based on approach to equilibrium. Requires the following attributes
        to be set:
        global_reactions - a list of global reactions in the same syntax as elementary expressions,
            with each one followed by its respective approach to equilibrium.
        pressure_mode - must be set to 'approach_to_equilibrium'
        Note that this function is not well-tested and should be used with caution.
        """

        print('APPROACH TO EQULIBRIUM PRESSURE')   
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
        print([float(i) for i in self.gas_pressures])

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
   
    ##############################################################
    ### _equilibrium functions intend to serve as utility      ### 
    ###  for initial guess estimation routines. (under dev)    ###  
    ##############################################################

    def set_affine_pressure_equilibrium(self,alpha,x0=[]):
        """ defines gas_pressure as an affine combination between the actual pressure
        and that of equilibrium by an alpha factor
        """
        eq_dict = self.get_pressure_equilibrium()
        xeq = eq_dict['xeq']; 
        if not x0:
            x0 = eq_dict['x0']
        else:
            pass
        xalpha = alpha*np.array(x0)+(1.-alpha)*xeq
        self.gas_pressures = self.pressure*xalpha
        return xalpha, xeq, x0        
    
    def get_pressure_equilibrium(self,xguess=None,ftol=1e-5):
        # Calaculate Equilirium Conversion
        gas_species = self.gas_names
        sps_dict = self.species_definitions
        n_atoms = max([len(self.species_definitions[_]['composition'].keys()) \
                                for _ in gas_species])
        # Create atomic constraints 
        Aeq = np.zeros([n_atoms,len(gas_species)], dtype=np.float128)
        beq = np.zeros(n_atoms,dtype=np.float128)
        x0 = np.zeros(len(gas_species),dtype=np.float128)
        atoms = {}
        self.pressure = float(self.pressure)
        self.gas_pressures = np.array([self.pressure*_/sum(self.gas_pressures) for _ in self.gas_pressures],dtype=np.float128)
        for _ in range(len(gas_species)):
            sps = gas_species[_]
            comp = self.species_definitions[sps]['composition']
            for atom in comp.keys():
                if atom in atoms.keys():
                    pass
                else:
                    atoms[atom] = len(atoms)
                Aeq[atoms[atom],_] = comp[atom]
                if hasattr(self,'gas_pressures'):
                        beq[atoms[atom]] += comp[atom]*(self.gas_pressures[_]/self.pressure)#*self.pressure/self.kBT
                        x0[_] = self.gas_pressures[_]/self.pressure
                elif 'concentration' in sps_dict[sps].keys() and\
                    sps_dict[sps]['concentration'] > 0:
                        beq[atoms[atom]] += comp[atom]*sps_dict[sps]['concentration']#*self.pressure/self.kBT
                        x0[_] = sps_dict[sps]['concentration']
                elif 'pressure' in sps_dict[sps].keys() and\
                    sps_dict[sps]['pressure'] > 0:
                        beq[atoms[atom]] += comp[atom]*(sps_dict[sps]['pressure']/self.pressure)#*sps_dict[sps]['pressure']/self.kBT
                        x0[_] = sps_dict[sps]['pressure']/self.pressure
       
        # Optimization section
        from scipy.optimize import minimize, Bounds
    
        # First minimization problem, x should be a molar fraction assuming static pressure
        # Calculate the systems' ultimate equilibrium composition
        
        thermo_corrections = self._correction_dict
        energies = np.array([sps_dict[_]['formation_energy']+thermo_corrections[_] for _ in gas_species],dtype=np.float128)
        
        log0 = lambda x : np.nan_to_num(np.log(x)*(np.isfinite(np.log(x))))
        
        # Objective function
        def total_gibbs(x):
            G0 = np.dot(energies,x)/float(self.kBT)
            Sconf = np.dot(log0(x),x)
            return (G0+Sconf+np.log(self.pressure*sum(x)))
        
        # Jacobian
        
        jac_total_gibbs = lambda x : np.array(energies/float(self.kBT)+np.log(x+1e-300)+np.ones(len(x)),dtype=np.float128)
        
        # Bounds
        bounds = Bounds(*[[0. for _ in range(len(gas_species))],[np.inf for _ in range(len(gas_species))]])
    
        # Equality constraints
        #Aeq[-1,:] = np.array(1.,dtype=np.float128)
        #beq[-1] = np.array(1.,dtype=np.float128)
        eq_cons = {'type': 'eq',
                   'fun' : lambda x: np.dot(Aeq,x)-beq,
                   'jac' : lambda x: Aeq}
        
        ineq_cons = {'type': 'ineq',
                   'fun' : lambda x: x,
                   'jac' : lambda x: np.eye(*np.shape(x))}
        
        xeq = minimize(total_gibbs, xguess if xguess else x0, method='SLSQP', jac=jac_total_gibbs,
                       constraints=[eq_cons,ineq_cons], options={'ftol': ftol, 'disp': False, 'maxiter':1000,'eps': 1.4901161193847656e-08},
                       bounds=bounds)

        if not xeq.success and not xeq.status == 4:
            print(xeq)
            print(Aeq.dot(xeq.x))
            print(Aeq)
            print(beq)
            print(Aeq.dot(xeq.x)-beq)
            print(sum(xeq.x))
            raise Exception()
        
        return dict(zip(('xeq','x0','Aeq','beq','feq','G0'),(xeq.x, x0, Aeq, beq, xeq.fun, total_gibbs(x0))))
    
    def set_equilibrated(self):
        """ Set reactants/products as their equilibrium composition a priori
            such that close-to-equilibrium species TOF do not overshadow thos of
            species that are far from it.
            It will look for `.equilibrated` parameter.
            This function assume 'static' pressure, such the extent of reactions
            toward equilibrium would not lead to significant change in the systems'
            total pressure.
        """
        print('SET EQUILIBRATED.')
        
        gas_species = self.gas_names
        
        xeq, x0, Aeq, beq, Geq, G0 = self.get_pressure_equilibrium()
        
        # Total Gibbs
        def total_gibbs(x):
            G0 = np.dot(energies,x)/float(self.kBT)
            Sconf = np.dot(log0(x),x)
            return G0+Sconf+np.log(self.pressure)
        
        # Second optimization problem
        # Fix specified components at final equilibrium composition and find
        # the initial compositions with minimum change from the initial one 
        # that still satisfied the imposed constraints
        fix_pos = list(filter(lambda x : isinstance(x,int), [_ if gas_species[_] in self.equilibrated else None \
                        for _ in range(len(gas_species))]))
        Aeq2 = np.zeros([len(self.equilibrated),len(gas_species)],dtype=np.float128)
        Aeq2[range(len(self.equilibrated)),fix_pos] = 1.
        Aeq2 = np.concatenate((Aeq,Aeq2),0)
        beq2 = np.concatenate((beq.flatten(),xeq[fix_pos]),0)
        
        # second objective function
        diff_x_xeq = lambda x : np.dot(x-x0,x-x0)
        
        # second jacobian
        jacobian_diff_x_xeq = lambda x : 2.*(x-x0)
        
        # Bounds
        bounds = Bounds(*[[0. for _ in range(len(gas_species))],[np.inf for _ in range(len(gas_species))]])
    
        # Equality constraints 2
        eq_cons2 = {'type': 'eq',
                   'fun' : lambda x: np.dot(Aeq2,x)-beq2,
                   'jac' : lambda x: Aeq2}       
    
        xf = minimize(diff_x_xeq, x0, method='SLSQP', jac=jacobian_diff_x_xeq,
                       constraints=eq_cons2, options={'ftol': 1e-8, 'disp': False, 'maxiter':1000},
                       bounds=bounds)        
        
        # set pressures
        self.gas_pressures = (xeq.x/sum(xeq.x))*self.pressure
        
        """
        print('x0: {}'.format(x0))
        print('xeq: {}'.format(xeq.x))
        print('xf: {}'.format(xf.x))
        """
        return dict(zip(('xeq','xf','x0','Aeq','beq','Aeq2','beq2','Geq','Gf','G0'),(xeq, xf.x, x0, Aeq, beq, Aeq2, beq2, Geq, total_gibbs(xf.x),G0)))

    ##############################################################
    ##############################################################

def fit_shomate(Ts, Cps, Hs, Ss, params0=[], plot_file = None):
    """This regression functionality has been updated from a non-linear version
        to a linearized one which does not need initial guesses (params0).
        params0 was kept as parameter to the function for backwards compatibiliy.
        It should be following deprecated.
    """
    try:
        from scipy.linalg import solve
    except ImportError:
        print("Scipy not installed, returning initial guess")
        return None
        #leastsq = lambda resid, initial, **kwargs: [initial, True]
    
    Ts, Cps, Hs, Ss = [np.array(_) for _ in [Ts,Cps,Hs,Ss]]
    ts = Ts/1000.
     
    # H-collocation matrix
    H_matrix = np.zeros([len(ts),7])
    H_matrix[:,0] = ts
    H_matrix[:,1] = ts**2/2.
    H_matrix[:,2] = ts**3/3.
    H_matrix[:,3] = ts**4/4.
    H_matrix[:,4] = -1./ts
    H_matrix[:,5] = 1.
    H_matrix[:,6] = 0.
    
    # S-collocation matrix
    S_matrix = np.zeros([len(ts),7])
    S_matrix[:,0] = np.log(ts)
    S_matrix[:,1] = ts
    S_matrix[:,2] = ts**2/2.
    S_matrix[:,3] = ts**3/3.
    S_matrix[:,4] = -1./(2.0*ts**2)
    S_matrix[:,5] = 0.
    S_matrix[:,6] = 1.
    
    # Cp-collocation matrix
    Cp_matrix = np.zeros([len(ts),7])
    Cp_matrix[:,0] = 1.
    Cp_matrix[:,1] = ts
    Cp_matrix[:,2] = ts**2.
    Cp_matrix[:,3] = ts**3.
    Cp_matrix[:,4] = 1./(ts**2)
    Cp_matrix[:,5] = 0.
    Cp_matrix[:,6] = 0.
    
    M = np.concatenate((H_matrix,-np.diag(ts).dot(S_matrix),np.diag(ts).dot(Cp_matrix)),0)
    b = np.concatenate((Hs,-np.diag(ts).dot(Ss),np.diag(ts).dot(Cps)),0)
    
    A, B, C, D, E, F, G = solve(M.T.dot(M),M.T.dot(b))
    
    if plot_file:
        import pylab as plt
        print('R2: {}, sumE: {}'.format(np.corrcoef(M.dot([A, B, C, D, E, F, G]),b),(M.dot([A, B, C, D, E, F, G])-b).T.dot(M.dot([A, B, C, D, E, F, G])-b)))
        fig = plt.figure()
        ax1 = fig.add_subplot(131)
        ax2 = fig.add_subplot(132)
        ax3 = fig.add_subplot(133)
        ax1.plot(ts,H_matrix.dot(np.array([A,B,C,D,E,F,G])),'-k')
        ax1.plot(ts,Hs,'ok')
        ax1.set_title('H')
        ax2.plot(ts,S_matrix.dot(np.array([A,B,C,D,E,F,G])),'-r')
        ax2.plot(ts,Ss,'or')
        ax2.set_title('S')
        ax3.plot(ts,Cp_matrix.dot(np.array([A,B,C,D,E,F,G])),'-b')
        ax3.plot(ts,Cps,'ob')
        ax2.set_title('Cps')
        fig.savefig(plot_file)

    H_c = 0
    return [A,B,C,D,E,F,G,H_c]

def harmonic_to_shomate(frequencies,Tmin,Tmax,resolution):
    """Generate Shomate parameters as of frequency data by using `fit_shomate`.
    """
    from scipy.interpolate import pchip_interpolate
    frequency_unit_conversion = 1.239842e-4;
    therm = HarmonicThermo([_*frequency_unit_conversion for _ in frequencies])
    _kJmol2eV = 0.01036427
    
    Ts = np.linspace(Tmin,Tmax,resolution)
    Ss = []
    Hs = []
    ZPE = frequency_unit_conversion*sum(frequencies)/2.0 
    for t in Ts:
        Ss += [therm.get_entropy(t,verbose=False)/_kJmol2eV*1e3]
        Hs += [(therm.get_internal_energy(t,verbose=False) - ZPE)/_kJmol2eV]
    Cps = pchip_interpolate(Ts,Hs,Ts,der=1)*1000.
    
    
    def _shomate_eq(params,temperature=[]):
        _kJmol2eV=0.01036427
        temperature_ref = 298.15
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
        
        dH = H(temperature,params) - H(temperature_ref,params)
        dS = S(temperature,params)
        Cp_ref = Cp(temperature_ref,params)
        dH = (temperature_ref*Cp_ref/1000.0 + dH)*(_kJmol2eV) #eV
        dS = dS*(_kJmol2eV/1e3) #eV/K
        return dH, dS, Cp_ref
    
    return fit_shomate(Ts,Cps,Hs,Ss)

   
