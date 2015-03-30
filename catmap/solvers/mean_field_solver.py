from solver_base import *
from catmap.data import templates
from copy import copy
import mpmath as mp
from catmap.functions import numerical_jacobian

class MeanFieldSolver(SolverBase):
    """Class for handling mean-field type kinetic models. Can be sub-classed to
    get functionality for steady-state solutions, sabatier solutions, etc."""

    def __init__(self,reaction_model=ReactionModel()):
        SolverBase.__init__(self,reaction_model)
        defaults = dict(
                tolerance = 1e-35,
                perturbation_size = 1e-14,
                )
        self._rxm.update(defaults)
        self._log_strings = {
                            'jacobian_fail':
                            "stagnated or diverging (residual = ${resid})."+\
                            " Assuming Jacobian is 0.",
                            }
    
    def get_rxn_rates(self,coverages,rate_constants):
        rates = self.elementary_rates(
                rate_constants,
                coverages,
                self.gas_pressures,
                self._mpfloat,
                self._matrix
                )
        return [ri for ri in rates]

    def get_rate(self, rxn_parameters, 
            coverages=None, verify_coverages=True,
            **coverage_kwargs):
        if not coverages: 
            coverages = self._coverage
        if not coverages: 
            raise ValueError('Input coverages to use as an initial guess '+\
                    'for the solver.')
        if verify_coverages == True:
            coverages = self.get_coverage(
                    rxn_parameters,coverages,**coverage_kwargs)
        self._coverage = coverages
        rate_constants = self.get_rate_constants(rxn_parameters,coverages)
        rates =  self.get_rxn_rates(coverages,rate_constants)
        return rates

    def get_turnover_frequency(self,rxn_parameters,rates=None,verify_coverages=True):
        rxn_parameters = list(rxn_parameters)
        if rates is None:
            rates = self.get_rate(rxn_parameters,verify_coverages=verify_coverages)

        def gas_to_idxs(gas):
            idxs = []
            for i,rxn in enumerate(self.elementary_rxns):
                if gas in rxn[0]:
                    idxs.append(-(i+1))
                elif gas in rxn[-1]:
                    idxs.append(i+1)
            return idxs

        def idx_to_rate(idx,rates):
            if idx < 0:
                mult = -1.0
            elif idx > 0:
                mult = 1.0
            i = abs(idx) -1
            return mult*rates[i]

        turnover_freq = []
        for g in self.gas_names:
            idxs = gas_to_idxs(g)
            tof = sum([idx_to_rate(idx,rates) for idx in idxs])
            turnover_freq.append(tof)
        
        self._turnover_frequency = turnover_freq
        return turnover_freq

    def get_selectivity(self,rxn_parameters):
        tofs = self.get_turnover_frequency(rxn_parameters)
        if self.products is None:
            self.products = [g for g,r in zip(self.gas_names,tofs) if r >0]
        if self.reactants is None:
            self.reactants = [g for g,r in zip(self.gas_names,tofs) if r <=0]
        prod_rate = sum([max(r,0)
            for g,r in zip(self.gas_names,tofs) if g in self.products])
        reac_rate = sum([max(-r,0)
            for g,r in zip(self.gas_names,tofs) if g in self.reactants])

        selectivities = []
        for g,r in zip(self.gas_names,tofs):
            if g in self.products and prod_rate:
                sel = max(r,0)/prod_rate

            elif g in self.reactants and reac_rate:
                sel = max(-r,0)/reac_rate
            else:
                sel = 0
            selectivities.append(sel)

        self._selectivities = selectivities
        return selectivities

    def get_rate_control(self,rxn_parameters):
        kT = self._kB*self.temperature
        eps = self._mpfloat(self.perturbation_size)
        try:
            dRdG = numerical_jacobian(self.get_turnover_frequency,rxn_parameters,self._matrix,eps)
        except ValueError, strerror:
            resid = str(strerror).rsplit('=',1)[1]
            resid = resid.replace(')','')
            resid.strip()
            self.log('jacobian_fail',resid=resid)
            dRdG = np.zeros((len(self.gas_names),len(self.adsorbate_names+self.transition_state_names)))

        t0 = self.get_turnover_frequency(rxn_parameters)
        dRdG *= -kT
        dRdG = dRdG.tolist()
        DRC = []
        for ti, Ji in zip(t0,dRdG):
           if ti == 0:
                DRC.append([0.0]*len(Ji))
           else:
                DRC.append([float(Jj/ti) for Jj in Ji])
	return DRC

    def get_interacting_energies(self,rxn_parameters):
        all_ads = self.adsorbate_names + self.transition_state_names
        N_ads = len(all_ads)
        energies = rxn_parameters[:N_ads]
        eps_vector = rxn_parameters[N_ads:]
        cvg = self._coverage + [0]*len(self.transition_state_names)
        E_int = self.interaction_function(cvg,energies,eps_vector,self.thermodynamics.adsorbate_interactions.interaction_response_function,False)[0]
        return E_int

    def get_selectivity_control(self,rxn_parameters):
        kT = self._kB*self.temperature
        eps = self._mpfloat(self.perturbation_size)
        try:
            dSdG = numerical_jacobian(self.get_selectivity,rxn_parameters,self._matrix,eps)
        except ValueError,strerror:
            resid = str(strerror).rsplit('=',1)[1]
            resid = resid.replace(')','')
            resid.strip()
            self.log('jacobian_fail',resid=resid)
            dRdG = np.zeros((len(self.gas_names),len(self.adsorbate_names+self.transition_state_names)))

        s0 = self.get_selectivity(rxn_parameters)
        dSdG *= -kT
        dSdG = dSdG.tolist()
        DSC = []
        for si, Ji in zip(s0,dSdG):
            if si == 0:
                DSC.append([0.0]*len(Ji))
            else:
                DSC.append([float(Jj/si) for Jj in Ji])
        return DSC
        
    def get_rxn_order(self,rxn_parameters,epsilon=1e-10):
        current_tofs = self.get_turnover_frequency(rxn_parameters)
        current_Ps = [p for p in self.gas_pressures]
        DRC = []

        for i,p in enumerate(current_Ps):
            new_p = copy(current_Ps)
            new_p[i] = current_Ps[i]*(1+epsilon)
            self._rxm.gas_pressures = new_p ##HACK
            #setting self.gas_pressures = new_p inexplicably breaks the solver.
            new_tofs = self.get_turnover_frequency(rxn_parameters)
            DRC_i = []
            for j,old_tof,new_tof,gas in zip(
                    range(0,len(current_tofs)),current_tofs,new_tofs,self.gas_names):
                if old_tof == 0:
                    if new_tof == 0:
                        dTOF = 0
                    else:
                        dTOF = (new_tof-old_tof)/new_tof
                else:
                    dTOF = (new_tof-old_tof)/old_tof
                dP = (new_p[i] - current_Ps[i])/current_Ps[i]
                DRC_i.append(float(dTOF/dP))
            DRC.append(DRC_i)
        self._rxm.gas_pressures = current_Ps ##HACK
        #setting self.gas_pressures = current_Ps inexplicably breaks the solver.
        self._rxn_order = DRC
        return DRC

    def get_apparent_activation_energy(self,rxn_parameters,epsilon=1e-10):
        #returns apparent Arrhenius activation energies (in units of R)
        #for production/consumption of each gas phase species.
        #Calculated as
        #E_app = T^2(dlnr_+/dT)=(T^2/r_+)(dr_+/dT), where r+ is the TOF
        current_tofs = self.get_turnover_frequency(rxn_parameters)
        current_T = self.temperature
        new_T = current_T*(1+epsilon)
        dT = new_T-current_T
        self.temperature = new_T
        descriptors = list(self._rxm.mapper._descriptors) #don't overwrite them, if temperature is a descriptor
        if 'temperature' in self._rxm.descriptor_names:
                index = self._rxm.descriptor_names.index('temperature')
                descriptors[index] = new_T
        rxn_parameters_newT = self._rxm.scaler.get_rxn_parameters(descriptors)
        new_tofs = self.get_turnover_frequency(rxn_parameters_newT)
        E_apps = []
        R = 8.31447e-3/96.485307#units of eV

        for i,gas in enumerate(self.gas_names):
            barriers_i = []
            dlnTOF = mp.log(new_tofs[i])-mp.log(current_tofs[i]) #this will fail if any of the TOFs are 0.
            E_app = R*float(dlnTOF.real)/dT*(current_T**2)
            E_apps.append(E_app)

        self.temperature = current_T
        self._apparent_activation_energy = E_apps
        #self.get_turnover_frequency(rxn_parameters)
	print E_apps
        return E_apps

    def summary_text(self):
        return ''

    def rate_equation_term(self,species_list,rate_constant_string,d_wrt=None):
        """Function to compose a term in the rate equation - e.g. kf[1]*theta[0]*p[0]"""

        #This clause allows for multiple site types. 
        site_indices={}

        gas_idxs = [self.gas_names.index(gas) 
                for gas in species_list if gas in self.gas_names]
        ads_idxs = [self.adsorbate_names.index(ads) 
                for ads in species_list if ads in self.adsorbate_names]
        sites = [s for s in species_list 
                if s in self.site_names] #allows for multiple site types
        if len(gas_idxs+ads_idxs+sites) != len(species_list):
            raise ValueError('Undefined species in '+','.join(species_list))

        rate_string = rate_constant_string

        if not d_wrt:
            for id in gas_idxs:
                rate_string += '*p['+str(id)+']'
            for id in ads_idxs:
                rate_string += '*theta['+str(id)+']'
            for s in sites:
                rate_string += '*s['+str(self.site_names.index(s))+']'
            return rate_string
        else:
            d_idx = self.adsorbate_names.index(d_wrt)
            d_site = self.species_definitions[d_wrt]['site']
            if (d_idx in ads_idxs #expression is a function of d_wrt
                    and d_site not in sites): #empty site not there -> easy
                multiplier = ads_idxs.count(d_idx) #get order
                ads_idxs.remove(d_idx) #reduce order by 1   
                if multiplier != 1:
                    rate_string = str(multiplier)+'*'+rate_string
            elif (d_site in sites #empty site appears, 
                    and d_idx not in ads_idxs): #but not the adsorbate
                multiplier = sites.count(d_site)
                multiplier = -1*multiplier #this accounds for the 
                #fact that d_site/d_ads = -1
                sites.remove(d_site) #reduce the order of site by 1
                rate_string = str(multiplier)+'*'+rate_string
            elif (d_site in sites #function of adsorbate
                    and d_idx in ads_idxs): #and function of site (1-theta_i)
                #need to use chain rule...
                ads_mult = ads_idxs.count(d_idx)
                site_mult = sites.count(d_site)*-1 #negative 1 to account for 
                #d_site/d_ads
                ads_str = 'theta['+str(d_idx)+']'
                site_str = 's['+str(self.site_names.index(d_site))+']'
                sites = [s for s in sites if s != d_site] #remove d_site 
                #from the site list
                ads_idxs = [a for a in ads_idxs if a != d_idx] #remove d_idx 
                #for adsorbate list
                if (-site_mult-1):
                    op = '*'
                else:
                    op = ''
                mult_rule = '('+str(site_mult)+'*'+ads_str+op+'*'.join(
                        [site_str]*(-site_mult-1)) # cvg*d_site
                mult_rule += ' + '+str(ads_mult)+'*'+site_str+'*'.join(
                        [ads_str]*(ads_mult-1))+ ')' 
                rate_string += '*'+mult_rule
            else:
                return '0'
            for id in gas_idxs:
                rate_string += '*p['+str(id)+']'
            for id in ads_idxs:
                rate_string += '*theta['+str(id)+']'
            for s in sites:
                rate_string += '*s['+str(self.site_names.index(s))+']'
            return rate_string

    def site_string_list(self):
        """Function to compose an alalytic expression for the coverage of empty sites"""
        site_strings=[]
        site_totals={}
        for site in self.site_names:
            site_totals[site] = self._mpfloat(self.species_definitions[site]['total'])
            site_idxs = [[idx,self.species_definitions[ads]['n_sites']] for idx,ads in enumerate(self.adsorbate_names)
                    if self.species_definitions[ads]['site'] == site]
            site_str = repr(site_totals[site])
            for idx_i,nsites in site_idxs:
                site_str += ' - theta['+str(idx_i)+']'
            site_strings.append('('+site_str+')')
        return site_strings

    def substitutions_dict(self):
        """Dictionary of substitutions needed for static compiled functions"""
        subdict = {}
        subdict['temperature'] = 'T = '+repr(self.temperature)
        subdict['kB'] = 'kB = '+repr(self._kB)
        subdict['h'] = 'h = '+repr(self._h)
        subdict['prefactor_list'] =  'prefactor_list = [' + ', '.join(self.prefactor_list) + ']'
        subdict['n_adsorbates'] = 'n_adsorbates = '+str(len(self.adsorbate_names))
        subdict['n_transition_states'] = 'n_transition_states = '+str(len(self.transition_state_names))
        
        max_cvg_list = []
        for ads in self.adsorbate_names:
            site = self.species_definitions[ads]['site']
            default_max = self.species_definitions[site]['total']
            max_cvg_list.append(self.species_definitions[ads].get('max_coverage',default_max))
        subdict['max_coverage_list'] = 'max_coverage_list = ' + repr(max_cvg_list)

        idx_dict = {}
        surf_species = self.adsorbate_names+self.transition_state_names
        for s in self.site_names:
            idxs = [surf_species.index(a) for a in surf_species if
                    self.species_definitions[a]['site'] == s]
            if idxs:
                if self.adsorbate_interaction_model not in ['ideal',None]:
                    default_params = getattr(
                            self.thermodynamics.adsorbate_interactions,
                            'interaction_response_parameters',{})
                else:
                    default_params = {}
                F_params = self.species_definitions[s].get('interaction_response_parameters',default_params)
                idx_dict[s] = [idxs,self.species_definitions[s]['total'],F_params]
        subdict['site_info_dict'] =  'site_info_dict = ' + repr(idx_dict)
        return subdict

    def rate_equations(self):
        """Compose analytical expressions for the reaction rates and 
        change of surface species wrt time (dc/dt).
        Assumes:
        
        kf is defined as a list of forward rate-constants
        kr is defined as a list of reverse rate-constants
        theta is defined as a list of coverages
        p is defined as a list of pressures
        """

        site_strings = self.site_string_list()

        rate_strings = []
        rate_strings.append('s = [0]*'+str(len(self.site_names)))
        for i,s in enumerate(site_strings):
            rate_strings.append('s['+str(i)+'] = '+site_strings[i])
        for i,rxn in enumerate(self.elementary_rxns):
            fRate_string = self.rate_equation_term(rxn[0],'kf['+str(i)+']')
            rRate_string = self.rate_equation_term(rxn[-1],'kr['+str(i)+']')
            rateString = 'r['+str(i)+'] = '+fRate_string + ' - ' + rRate_string
            rate_strings.append(rateString)

        dcdt_strings = []
        for i,ads in enumerate(self.adsorbate_names):
            dcdt_str = 'dtheta_dt['+str(i)+'] = '
            for j,rxn in enumerate(self.elementary_rxns):
                rxnCounts = [-1.0*rxn[0].count(ads), 1.0*rxn[-1].count(ads)]
                rxnOrder = [o for o in rxnCounts if o]
                if rxnOrder: 
                    rxnOrder = rxnOrder[0]
                    dcdt_str += ' + ' + str(int(rxnOrder))+'*r['+str(j)+']'
            if dcdt_str.endswith('= '):
                dcdt_str += '0'
            dcdt_strings.append(dcdt_str)

        all_strings = rate_strings + dcdt_strings
        return all_strings

    def jacobian_equations(self,adsorbate_interactions=True):
        """Composes analytical expressions for the Jacobian matrix.
        Assumes:

        kf is defined as a list of forward rate-constants
        kr is defined as a list of reverse rate-constants
        theta is defined as a list of coverages
        p is defined as a list of pressures

        If the rate constants depend on coverage, use 
        adsorbate_interactions = True. 
        Assumes:

        kB is defined as Boltzmann's constant
        T is defined as the temperature
        dEf is defined as a list of lists where dEf[i][j] is the
            derivative of forward activation free energy i wrt coverage j
        dEr is defined as a list of lists where dEr[i][j] is the
            derivative of reverse activation free energy i wrt coverage j
        """
        
        site_strings = self.site_string_list()

        J_strings = []
        J_strings.append('s = [0]*'+str(len(self.site_names)))
        for i,s in enumerate(site_strings):
            J_strings.append('s['+str(i)+'] = '+site_strings[i])

        if adsorbate_interactions == True:
            J_strings.append('kfkBT = [0]*'+str(len(self.elementary_rxns)))
            J_strings.append('krkBT = [0]*'+str(len(self.elementary_rxns)))
            for i in range(len(self.elementary_rxns)):
                J_strings.append('kfkBT['+str(i)+'] = kf['+str(i)+']/kBT')
                J_strings.append('krkBT['+str(i)+'] = kr['+str(i)+']/kBT')

        for i,ads_i in enumerate(self.adsorbate_names):
            for j,ads_j in enumerate(self.adsorbate_names):
                J_str = 'J['+str(i)+']['+str(j)+'] = 0'
                for k,rxn in enumerate(self.elementary_rxns):
                    rxnCounts = [-1.0*rxn[0].count(ads_i), 
                            1.0*rxn[-1].count(ads_i)]
                    rxnOrder = [o for o in rxnCounts if o]
                    if rxnOrder: 
                        rxnOrder = rxnOrder[0]
                        fRate_string = self.rate_equation_term(rxn[0],'kf['+str(k)+']',ads_j)
                        rRate_string = self.rate_equation_term(rxn[-1],'kr['+str(k)+']',ads_j)
                        if adsorbate_interactions == True:
                            dfRate_string = self.rate_equation_term(rxn[0],'(kfkBT['+str(k)+'])*dEf['+str(k)+']['+str(j)+']')
                            drRate_string = self.rate_equation_term(rxn[-1],'(krkBT['+str(k)+'])*dEr['+str(k)+']['+str(j)+']')
                            fRate_string += ' + ' + dfRate_string
                            rRate_string += ' - ' + drRate_string
                        if fRate_string != '0' and rRate_string != '0':
                            dr_dx = '('+fRate_string + ' - ' + rRate_string+')'
                        elif fRate_string != '0':
                            dr_dx = fRate_string
                        elif rRate_string != '0':
                            dr_dx = '-1*'+rRate_string
                        elif fRate_string == rRate_string == '0':
                            dr_dx = None
                        if dr_dx:
                            J_str += ' + ' + str(int(rxnOrder))+'*'+dr_dx
                J_strings.append(J_str)
        return J_strings

    def reaction_energy_equations(self,adsorbate_interactions=True):
        """Composes a list of analytical expressions which give the reaction 
        and activation energies for elementary steps. Note that while this 
        is useful primarily for models with adsorbate-interactions 
        (otherwise these energetics can easily be obtained by the reaction
        model itself), they are technically valid for all mean-field models.
        Assumes:
       
        Gf is a list of formation energies ordered as 
            adsorbate_names+transition_state_names

        If model includes adsorbate interactions then use 
        adsorbate_interactions = True to include dEa/dtheta in the output.
        Assumes:

        dGs is a matrix/array of derivatives of free energies wrt coverages
            such that dGs[:,i] is a vector of derivatives of the free energy
            of species i wrt each coverage ordered as adsorbate_names
            
        """

        idx_dict = {}
        for i,ads in enumerate(self.adsorbate_names):
            idx_dict[ads] = str(i)
        for i,TS in enumerate(self.transition_state_names):
            idx_dict[TS] = str(i + len(self.adsorbate_names))

        expressions = []
        n_rxns = len(self.elementary_rxns)
        expressions.append('G_IS = [0]*'+str(n_rxns))
        expressions.append('G_TS = [0]*'+str(n_rxns))
        expressions.append('G_FS = [0]*'+str(n_rxns))
        expressions.append('G_af = [0]*'+str(n_rxns))
        expressions.append('G_ar = [0]*'+str(n_rxns))
        if adsorbate_interactions == True:
            expressions.append('dG_IS = [0]*'+str(n_rxns))
            expressions.append('dG_TS = [0]*'+str(n_rxns))
            expressions.append('dG_FS = [0]*'+str(n_rxns))

        def species_strings(state_list,list_name,include_constants=True,type='list'):
            species_strs = []
            for species in state_list:
                if species in idx_dict:
                    if type == 'list':
                        species_strs.append(list_name+'['+idx_dict[species]+']')
                    elif type == 'matrix_row':
                        species_strs.append(list_name+'['+idx_dict[species]+',:]')
                    elif type == 'matrix_col':
                        species_strs.append(list_name+'[:,'+idx_dict[species]+']')
                elif include_constants == True:
                    if species in self.gas_names:
                        idx = self.gas_names.index(species)
                        species_strs.append('gas_energies['+str(idx)+']')
                    elif species in self.site_names:
                        idx = self.site_names.index(species)
                        species_strs.append('site_energies['+str(idx)+']')
                    else:
                        raise ValueError('Undefined species '+species)
            if not species_strs:
                species_strs = ['[0]*'+str(len(self.adsorbate_names))]
            return species_strs

        for i,rx in enumerate(self.elementary_rxns):
            IS = rx[0]
            FS = rx[-1]
            if len(rx) > 2:
                TS = rx[1]
            else:
                TS = None

            txt = 'G_IS['+str(i)+'] = ' + ' + '.join(species_strings(IS,'Gf')) + '\n    '
            txt += 'G_FS['+str(i)+'] = ' + ' + '.join(species_strings(FS,'Gf')) + '\n    '
            if adsorbate_interactions == True:
                txt += 'dG_IS['+str(i)+'] = element_wise_addition([' + ' , '.join(species_strings(IS,'dGs',False,'list')) + '])\n    '
                txt += 'dG_FS['+str(i)+'] = element_wise_addition([' + ' , '.join(species_strings(FS,'dGs',False,'list')) + '])\n    '
            if TS:
                TS_strings = species_strings(TS,'Gf')
                txt += 'G_TS['+str(i)+'] = ' + ' + '.join(species_strings(TS,'Gf')) + '\n    '
                if adsorbate_interactions == True:
                    txt += 'dG_TS['+str(i)+'] = element_wise_addition([' + ' , '.join(species_strings(TS,'dGs',False,'list')) + '])\n    '
            else:
                txt += 'G_TS['+str(i)+'] = max([G_IS['+str(i)+'],G_FS['+str(i)+']])' + '\n    '
                if adsorbate_interactions == True:
                    txt += 'dG_TS['+str(i)+'] = None #determined later\n    '

            txt = txt.replace(' + 0', '')
            expressions.append(txt)
        return expressions

