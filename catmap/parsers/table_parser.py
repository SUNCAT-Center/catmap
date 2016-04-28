import numpy as np
import catmap
from parser_base import *
string2symbols = catmap.string2symbols
Template = catmap.Template

class TableParser(ParserBase):
    """Parses attributes based on column headers and filters.

    Additional functionality may be added by inheriting and defining
        the parse_{header_name} function where header_name is the 
        column header for the additional variable to be parsed.
    """
    def __init__(self,reaction_model,**kwargs):
        ParserBase.__init__(self,reaction_model)
        defaults = dict(
                estimate_frequencies = 1, #Use frequencies from different sites 
                #if available (set variable to 1 or True). 
                #Use dissociated state frequencies for TS (set to 2)
                #If no frequencies are available from other sites then 
                #concatenate frequencies from 
                #individual atoms (set to 3).
                #If no frequencies can be found, use empty frequency set
                #(set to >3)
                frequency_surface_names = [], #Use frequencies from a specific 
                #surface_name only. If "None" or empty then an average of 
                #the frequencies from all available surfaces will be used.
                required_headers = ['species_name','surface_name','site_name'
                                    ,'formation_energy','frequencies',
                                    'reference'],
                parse_headers = ['formation_energy','frequencies'],
                frequency_unit_conversion = 1.239842e-4, # conversion factor to 
                coverage_headers = ['coverage','coadsorbate_coverage'],
                #go from input units to eV
                standard_coverage = 'min',
                standard_coadsorbate_coverage = 'min',
                #coverage to use as the "base" in coverage-dependent input file
                #use "min" to take the minimum or specify explicitly
                interaction_surface_names = None,
                #use a different set of (more) surfaces to form interaction matrix.
                #If none then only the surfaces in the model will be used.
                )

        self._linebreak = '\n'
        self._separator = '\t'
        self._rxm.update(kwargs,override=True)
        self._rxm.update(defaults,override=False)
        self._required = {'input_file':str,'estimate_frequencies':bool,
                               'required_headers':list,
                               'parse_headers':list,
                               'frequency_unit_conversion':float,
                               'frequency_surface_names':None}

    def parse(self,**kwargs):
        f = open(self.input_file)
        lines = f.read().split(self._linebreak)
        lines = [L for L in lines if L]
        f.close()
        
        self._baseparse()


        headers = lines.pop(0).split(self._separator)
        headers = [h.strip() for h in headers]
        if not set(self.required_headers).issubset(set(headers)):
            raise ValueError('Required headers are missing! '+\
                             'Please be sure that all headers '+\
                             'are specified: '+' '.join(self.required_headers))
        linedicts = []
        for L in lines:
            linedict = {}
            for k, v in zip(headers, 
                    L.split(self._separator, len(headers))):
                linedict[k] = v
            if len(linedict) != len(headers):
                print("Input line " + str(linedict) + " does not have all required fields.  Ignoring.")
                continue
            sites = [s for s in self.species_definitions if
                    self.species_definitions[s].get('type',None) == 'site' and 
                    linedict['site_name'] in 
                    self.species_definitions[s]['site_names']
                    and '*' not in s]
            if not sites:
                sites = ['?']
            adskey = [linedict['species_name']+'_'+site_i for site_i in sites]
            linedict['species_keys'] = adskey
            linedicts.append(linedict)


        self._line_dicts = linedicts
        self._headers = headers

        for p in self.parse_headers:
            if callable(getattr(self,'parse_'+p)):
                getattr(self,'parse_'+p)()
            else:
                raise AttributeError('No parsing function defined for '+p)

    def parse_formation_energy(self,**kwargs):
        "Parse in basic info for reaction model"

        self.__dict__.update(kwargs)
       
        all_ads = [k for k in self.species_definitions.keys()
                   if self.species_definitions[k].get('type',None) != 'site']

        for adsdef in all_ads:
            if adsdef in ['ele_g', 'H_g', 'OH_g']:
                continue  # Ignore electrochemical species that should not be manually defined.
            ads = self.species_definitions[adsdef].get('name',None)
            if ads is None:
                del self.species_definitions[adsdef]
                print('Warning: Species with undefined "name" was encountered ('+adsdef+'). '+\
                     'Ensure that all species which are explicitly set in "species_definitions" '+\
                     'are also defined in the reaction network ("rxn_expressions"). This definition '+\
                     'will be ignored.')
            else:
                site = self.species_definitions[adsdef]['site']
                alternative_names = self.species_definitions[adsdef].get(
                        'alternative_names',[])
                adsnames = [ads]+alternative_names

                sites = self.species_definitions[site]['site_names']
                infodict = {}
                for linedict in self._line_dicts:
                    if (
                            linedict['species_name'] in adsnames and 
                            linedict['site_name'] in sites and 
                            linedict['surface_name'] in list(self.surface_names)+['None']
                            ):
                        
                        #The following clause ensures that the low-coverage limit
                        #is used unless otherwise specified. 
                        #It should probably be abstracted out into something cleaner.
                        pass_dict = {}
                        surf = linedict['surface_name']
                        for cvg_key in ['coverage','coadsorbate_coverage']:
                            pass_dict[cvg_key] = True
                            if cvg_key in linedict:
                                standard_cvg = getattr(self,'standard_'+cvg_key, None)
                                if standard_cvg in ['min','minimum',None]:
                                    if surf in infodict:
                                        if linedict[cvg_key] > infodict[surf][cvg_key]:
                                            pass_dict[cvg_key] = False
                                else:
                                    if linedict[cvg_key] != standard_cvg:
                                        pass_dict[cvg_key] = False

                        if False not in pass_dict.values():
                            infodict[surf] = linedict
                
                paramlist = []
                sources = []
                if self.species_definitions[adsdef]['type'] not in ['gas']:
                    for surf in self.surface_names:
                        if surf in infodict:
                            E = float(infodict[surf]['formation_energy'])
                            paramlist.append(E)
                            sources.append(infodict[surf]['reference'].strip())
                        else:
                            paramlist.append(None)
                    self.species_definitions[adsdef]['formation_energy'] = paramlist
                    self.species_definitions[adsdef]['formation_energy_source'] = sources
                else:
                    if 'None' in infodict:
                        E = float(infodict['None']['formation_energy'])
                        self.species_definitions[adsdef]['formation_energy'] = E
                        self.species_definitions[adsdef]['formation_energy_source'] = \
                                infodict['None']['reference'].strip()
                    else:
                        raise ValueError('No formation energy found for '+str(adsdef)+'. Check input file.')

    def parse_frequencies(self,**kwargs):
        self.__dict__.update(kwargs)
        allfreqdict = {}
        frequency_dict = {}

        #Parse in all available frequencies
        for linedict in self._line_dicts:
            if eval(linedict['frequencies']):
                freqs = eval(linedict['frequencies'])
                freqs = [self.frequency_unit_conversion*f for f in freqs]
                if linedict['species_name'] not in allfreqdict:
                    allfreqdict[linedict['species_name']] = \
                        [[linedict['surface_name'], 
                         linedict['site_name'],
                         freqs]] #Store frequency info for parsing later
                else:
                    frq = [linedict['surface_name'], 
                         linedict['site_name'],
                         freqs]
                    if frq not in allfreqdict[linedict['species_name']]:
                        allfreqdict[linedict['species_name']].append(frq)

        def freq_handler(freqdict_entry,site,ads):
            """
            Returns a single list of frequencies from a freqdict_entry, which is a list
            of all frequency data for a given species.  Entries matching both site
            and surface (if specified in self.frequency_surface_names) are preferred over
            those that only match surface.  If more than match of the highest validity
            is found, the mean of those frequencies is returned.
            """
            perfect_matches = []
            partial_matches = []
            if self.frequency_surface_names is None:
                self.frequency_surface_names = []
            for entry in freqdict_entry:
                masked = [entry[0] in self.frequency_surface_names,
                        entry[1] in self.species_definitions.get(site,{'site_names':[]})['site_names'],
                          entry[2]]
                if not self.frequency_surface_names:
                    if site in self._gas_sites and entry[0] == 'None':
                        masked[0] = True
                    elif site not in self._gas_sites:
                        masked[0] = True
                else:
                    if site in self._gas_sites and entry[0] == 'None':
                        masked[0] = True

                if all(masked):
                    perfect_matches.append(masked[-1])
                elif masked[0] and site not in self._gas_sites: #Surface matches but site might not...
                    if entry[1] != 'gas': #HACK... this whole function needs to be cleaned up.
                        partial_matches.append(masked[-1])
            
            def match_handler(perfect_matches):
                if len(perfect_matches) == 1:
                    return perfect_matches[0]
                elif len(perfect_matches) > 1:
                    if len(set([len(pm) for pm in perfect_matches]))>1:
                        raise ValueError('Frequency vectors have different '+\
                                'lengths for '+ str(ads))
                    matcharray = np.array(perfect_matches)
                    freqout = matcharray.mean(0) #average valid frequencies
                    return list(freqout)
                else: #No valid frequencies are found...
                    return []
            
            if len(perfect_matches) > 0:
                return match_handler(perfect_matches)
            elif self.estimate_frequencies:
                return match_handler(partial_matches)
            else:
                return []

        all_ads = [k for k in self.species_definitions.keys()
                   if self.species_definitions[k]['type'] != 'site']

        for adsdef in all_ads+allfreqdict.keys(): #format all freqs
            if '_' in adsdef:
                adsname,site = adsdef.split('_')
            elif adsdef in allfreqdict.keys():
                adsname = adsdef
                site = self._default_site

            if adsname in allfreqdict:
                frequency_dict[adsdef] = freq_handler(allfreqdict[adsname],site
                        ,adsname)
            elif self.estimate_frequencies > 3:
                frequency_dict[adsdef] = []

        for adsdef in all_ads:
            adsname,site = [self.species_definitions[adsdef][k] 
                    for k in ['name','site']]
            #Use single-atom frequencies...
            if (
                    not frequency_dict.get(adsdef,None) and 
                    self.estimate_frequencies > 2 and 
                    '-' not in adsname #Don't include TS's
                    ):
                symbols = string2symbols(adsname)
                freqs = []
                if set(symbols).issubset(set(frequency_dict.keys())):
                    for s in symbols:
                        freqs += frequency_dict[s]
                frequency_dict[adsdef] = freqs

        for adsdef in all_ads:
            #Use dissosciated TS frequencies
            adsname,site = [self.species_definitions[adsdef][k] 
                    for k in ['name','site']]
            if (
                    not frequency_dict.get(adsdef,None) and
                    self.estimate_frequencies > 1 and
                    '-' in adsname
                    ):
                A,B = adsname.split('-')
                frequency_dict[adsdef] = frequency_dict[A] + frequency_dict[B]

        for key in self.species_definitions.keys():
            self.species_definitions[key]['frequencies'] = frequency_dict.get(key,[])

    def parse_coverage(self,**kwargs):
        self.__dict__.update(kwargs)
        
        n = len(self.adsorbate_names)
        surfaces = self.surface_names

        info_dict = {}
        ads_names = self.adsorbate_names+self.transition_state_names

        for surf in surfaces:
            cvg_dict = {}
            for linedict in self._line_dicts:
                for skey in linedict['species_keys']:
                    if (skey in self.adsorbate_names+self.transition_state_names
                            and linedict['surface_name'] == surf):

                        ads = skey
                        if 'delta_theta' in linedict:
                            self.species_definitions[ads]['delta_theta'] = float(
                                    linedict['delta_theta'])

                        theta_vec = [0]*len(ads_names)
                        idx_i = ads_names.index(ads)
                        theta_i = float(linedict['coverage'])
                        theta_vec[idx_i] += theta_i
                        for coads_name in ['coadsorbate','coadsorbate2']: 
                            #could add coadsorbate3, coadsorbate4,... as needed
                            if coads_name+'_name' in linedict:
                                if linedict[coads_name+'_name'] != 'None':
                                    coads = linedict[coads_name+'_name']
                                    site = ads.split('_')[-1]
                                    site = linedict.get(coads_name+'_site',site)
                                    coads += '_'+site #assume coads on same site as ads if not specified
                                    theta_j = float(linedict[coads_name+'_coverage'])
                                    if coads in ads_names:
                                        idx_j = ads_names.index(coads)
                                        theta_vec[idx_j] += theta_j
                                    else:
                                        names_only = [n.split('_')[0] for n in ads_names]
                                        coads_name = coads.split('_')[0]
                                        if coads_name not in names_only:
                                            print 'Warning: Could not find co-adsorbed species '\
                                            +coads+' (adsorbate '+ads+'). Ignoring this entry.'
                                        else:
                                            idx_j = names_only.index(coads_name)
                                            actual_ads = ads_names[idx_j]
                                            print 'Warning: Could not find co-adsorbed species '\
                                            +coads+' (adsorbate '+ads+'). Using '+actual_ads+'.'
                                            theta_vec[idx_j] += theta_j
                        E_diff = float(linedict['formation_energy'])
                        E_int = linedict.get('integral_formation_energy',None)
                        if E_int:
                            E_int = float(E_int)

                        theta_E = [theta_vec,
                                E_diff,E_int]

                        if ads in cvg_dict:
                            cvg_dict[ads].append(theta_E)
                        else:
                            cvg_dict[ads] = [theta_E]
            info_dict[surf] = cvg_dict
        
        for i_ads,ads in enumerate(self.adsorbate_names+self.transition_state_names):
            cvg_dep_E = [None]*len(surfaces)
            for surf in surfaces:
                cvgs = info_dict[surf].get(ads,None)
                if cvgs is None:
                    pass
                else:
                    cvg_dep_E[self.surface_names.index(surf)] = cvgs
            self.species_definitions[ads]['coverage_dependent_energy'] = cvg_dep_E
