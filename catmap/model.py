import os
import sys
import inspect
from copy import copy
import numpy as np
import catmap
from string import Template
import functions
import re
from data import regular_expressions
string2symbols = catmap.string2symbols
pickle = catmap.pickle
plt = catmap.plt
griddata = catmap.griddata

class ReactionModel:
    """
The central object that defines a microkinetic model consisting of:

- active sites
- species
- possible reaction steps
- rate constant expressions
- descriptors and descriptor ranges
- data files for energies
- external parameters (temperature, pressures)
- other more technical settings related to the solver and mapper

    """
    def __init__(self,**kwargs): #
        """Class for managing microkinetic models.

           :param setup_file: Specify <mkm-file> from which to load model.
           :type setup_file: str
            """

        #Set static utility functions
        for f in dir(functions):
            if not f.startswith('_') and callable(getattr(functions,f)):
                setattr(self,f,getattr(functions,f))

        self._kB = 8.617332478e-5 #eV/K (from NIST)
        self._h = 4.135667516e-15 #eV s (from NIST)
        self._default_site = 's'
        self._gas_sites = ['g']
        self._required = {'adsorbate_names':tuple,
                'transition_state_names':tuple,
                'gas_names':tuple,
                'descriptor_names':tuple,
                'surface_names':tuple,
                'species_definitions':dict,
                'elementary_rxns':tuple,
                'thermodynamics':None,
                'numerical_representation':None,
                'scaler':None,
                'solver':None,
                'mapper':None,
                } #These are requirements of any reaction model.

        self._classes = ['parser','scaler', 'thermodynamics',
                'solver','mapper']
        self._solved = None
        #keeps track of whether or not the model was solved

        #attributes for logging
        self._log_lines = []
        self._log_dict = {}

        self._log_strings = {'input_success':
                             'loaded all data from input file',
                             'header_fail':
                             'could not save ${save_txt}'}
        #modules to import in the log file to allow opening with python -i
        self._log_imports = "from numpy import array\n\n"+\
                            "import cPickle as pickle\n\n"
        #attrs taking up more space than this will be dumpted to self.data_file
        self._max_log_line_length = 8000 #100 lines at 80 cols
        #Character used for loop depth in std out
        self._pickle_attrs = [] #attributes to pickle rather than ASCII
        #place to store warnings
        self._warned = []

        #attributes for function generation
        self._function_templates = {}
        #dictionary of strings to be substituted into function templates
        self._function_substitutions = {}
        #dictionary of strings corresponding to all functions generated
        self._function_strings = {}

        #Possible outputs:

        ##Solver level
        #coverage
        #rate
        #production_rate
        #consumption_rate
        #turnover_frequency
        #rxn_direction
        #selectivity
        #rate_control
        #rate_constant
        #equilibrium_constant

        ##Scaler level
        #rxn_parameter
        #frequency
        #electronic_energy
        #free_energy
        #zero_point_energy
        #enthalpy
        #entropy

        self.output_variables = ['coverage','rate']
        #dictionary to save labels for plotting
        self.output_labels = {}

        for key in kwargs:
            setattr(self,key,kwargs[key])

        if hasattr(self,'setup_file'): #parse in setup file.
            #Note that the syntax is simply python variable definitions.
            #This is NOT idiot proof.
            self.model_name = self.setup_file.rsplit('.',1)[0]
            self.load(self.setup_file)


    # Functions for executing the kinetic model

    def run(self,**kwargs):
        """Run the microkinetic model. If recalculate is True then
        data which is re-loaded will be used as an initial guess; otherwise it
        will be assumed to be correct.

        :param recalculate: If True solve model again using previous results as initial guess
        :type recalculate: bool

        """

        for key in kwargs:
            setattr(self,key,kwargs[key])

        #ensure resolution has the proper dimensions
        if not hasattr(self.resolution,'__iter__'):
            self.resolution = [self.resolution]*len(self.descriptor_names)

        #set numerical representation
        if self.numerical_representation == 'mpmath':
            import mpmath as mp
            mp.mp.dps = self.decimal_precision + 25 #pad decimal precision
            self._math = mp
            self._log_imports += "\nfrom mpmath import mpf \n\n"
            self._kB = mp.mpf('8.617332478e-5') #eV/K (from NIST)
            self._h = mp.mpf('4.135667516e-15') #eV s (from NIST)
            self._mpfloat = mp.mpf
            self._matrix = mp.matrix
            self._Axb_solver = mp.lu_solve
            self._math.infnorm = lambda x: mp.norm(x,'inf')
        elif self.numerical_representation in ['numpy','python']:
            self._log_imports += "\nfrom numpy import matrix \n\n"
            self._math = np
            self._math.infnorm = lambda x: np.linalg.norm(x,np.inf)
            self._mpfloat = float
            def matrixT(*args,**kwargs):
                array = np.array(args[0])
                while 1 in array.shape:
                    sum_idx = array.shape.index(1)
                    array = array.sum(sum_idx)
                return array.T

            self._matrix = matrixT
            def Axb_solver(A,b):
                try:
                    return np.linalg.solve(A,b.T)
                except np.linalg.linalg.LinAlgError:
                    raise ZeroDivisionError
            self._Axb_solver = Axb_solver
            if self.decimal_precision > 15:
                print(
                'Warning: Max precision with numpy/python is 16 digits')
                self.decimal_precision = 15
        else:
            raise AttributeError(
            'Numerical representation must be mpmath, numpy, or python.')

        #set up interaction model
        if self.adsorbate_interaction_model == 'first_order':
            interaction_model = catmap.thermodynamics.FirstOrderInteractions(self)
            interaction_model.get_interaction_info()
            response_func = interaction_model.interaction_response_function
            if not callable(response_func):
                int_function = getattr(interaction_model,
                        response_func+'_response')
                interaction_model.interaction_response_function = int_function
            self.thermodynamics.__dict__['adsorbate_interactions'] = interaction_model

        elif self.adsorbate_interaction_model == 'second_order':
            interaction_model = catmap.thermodynamics.SecondOrderInteractions(self)
            interaction_model.get_interaction_info()
            response_func = interaction_model.interaction_response_function
            if not callable(response_func):
                int_function = getattr(interaction_model,
                        response_func+'_response')
                interaction_model.interaction_response_function = int_function
            self.thermodynamics.__dict__['adsorbate_interactions'] = interaction_model

        elif self.adsorbate_interaction_model in ['ideal',None]:
            self.thermodynamics.adsorbate_interactions = None
        else:
            raise AttributeError(
                    'Invalid adsorbate_interaction_model specified.')

        self.compatibility_check()

        #determine whether or not to (re-)solve model
        has_all = True
        for v in self.output_variables:
            if not hasattr(self,v+'_map'):
                has_all = False

        if not hasattr(self,'stdout'):
            #any re-loaded model will have stdout
            has_all = False

        if self._solved == self._token() and not getattr(self,'recalculate',None):
            #Do not solve the same model twice
            has_all = True

        elif has_all and not getattr(self,'recalculate',None):
            #All data has been loaded and no verification => solved.
            self._solved = self._token()
            self.log('input_success')

        else:
            ran_dsa = False #When no map exists, run descriptor space analysis first
            if not getattr(self,'coverage_map',None):
                #Make "volcano plot"
                if getattr(self,'descriptor_ranges',None) and getattr(self,'resolution',None):
                    self.descriptor_space_analysis()
                    ran_dsa = True

            #Get rates at single points
            if getattr(self,'descriptors',None):
                self.single_point_analysis(self.descriptors)
            #Get rates at multiple points
            if getattr(self,'descriptor_values',None):
                self.multi_point_analysis()

            #If a map exists, run descriptor space analysis last (so that single-point guesses are
            #not discarded)
            if not ran_dsa and getattr(self,'descriptor_ranges',None) and getattr(self,'resolution',None):
                    self.descriptor_space_analysis()


            #Save long attrs in data_file
            for attr in dir(self):
                if (not attr.startswith('_') and
                        not callable(getattr(self,attr)) and
                        attr not in self._classes):
                    if (len(repr(getattr(self,attr))) >
                            self._max_log_line_length):
                        #line is too long for logfile -> put into pickle
                        self._pickle_attrs.append(attr)
            pickled_data = {}
            for attr in self._pickle_attrs:
                pickled_data[attr] = getattr(self,attr)
            pickle.dump(pickled_data,open(self.data_file,'w'))

            #Make logfile
            log_txt = self._log_imports
            log_txt += self._header(exclude_outputs=self._pickle_attrs)
            self._stdout = '\n'.join(self._log_lines)
            log_txt += 'stdout = ' + '"'*3+self._stdout+'"'*3
            #this construction means that self.stdout will only be set
            #for models which have been re-loaded.
            self._solved = self._token()
            if hasattr(self,'log_file'):
                logfile = self.log_file
            else:
                name,suffix = self.setup_file.rsplit('.',1)
                if suffix != 'log':
                    suffix = 'log'
                else:
                    suffix = 'out'
                logfile = '.'.join([name,suffix])

            f = open(logfile,'w')
            f.write(log_txt)
            f.close()

            if getattr(self,'create_standalone',None):
                self.make_standalone()

    def descriptor_space_analysis(self):
        """
        Use mapper to create map/volcano-plot of rates/coverages.
        """
        self.mapper.get_output_map(self.descriptor_ranges,self.resolution)

    def single_point_analysis(self,pt):
        """
        Find rates/coverages at a single point.

        .. todo:: Explain pt argument
        """
        self.mapper.get_point_output(pt)
        for out in self.output_variables:
            mapp = getattr(self,out+'_map',[])
            if mapp is None:
                mapp = []
            mapp.append([pt,getattr(self,'_'+out)])
            setattr(self,out+'_map',mapp)

    def multi_point_analysis(self):
        """
        .. todo:: __doc__
        """
        for pt in self.descriptor_values:
            self.single_point_analysis(pt)

    def generate_static_functions(self):
        "Dynamically compile static functions"
        for func_name in self._function_templates:
            if not callable(getattr(self,func_name,None)):
                template = Template(self._function_templates[func_name])
                func_string = template.substitute(self._function_substitutions)
                self._function_strings[func_name] = func_string
                locs = {}
                exec func_string in globals(), locs
                setattr(self,func_name,locs[func_name])

   #File IO functions
    def load(self,setup_file): #
        """Load a 'setup file' by importing it and assigning all local
        variables as attributes of the kinetic model. Special attributes
        mapper, parser, scaler, solver will attempt to convert strings
        to modules."""
        defaults = dict(mapper='MinResidMapper',
                parser='TableParser',
                scaler='GeneralizedLinearScaler',
                solver='SteadyStateSolver',
                thermodynamics='ThermoCorrections',
                numerical_representation = 'mpmath',
                adsorbate_interaction_model = 'ideal',
                prefactor_list=None,
                interaction_fitting_mode=None,
                decimal_precision = 75,
                verbose = 1,
                data_file = 'data.pkl')
        globs = {}
        locs = defaults

        execfile(setup_file,globs,locs)
        for var in locs.keys():
            if var in self._classes:
                #black magic to auto-import classes
                try:
                    if locs[var]:
                        if not var.endswith('s'):
                            pyfile = var + 's'
                        else:
                            pyfile = var
                        basepath=os.path.dirname(
                                inspect.getfile(inspect.currentframe()))
                        if basepath not in sys.path:
                            sys.path.append(basepath)
                        sublocs = {}
                        _temp = __import__(pyfile,globals(),sublocs, [locs[var]])
                        class_instance = getattr(_temp,locs[var])(self)
                        setattr(self,var,class_instance)
                    else:
                        setattr(self,var,None)
                except ImportError:
                    raise AttributeError(var.capitalize()+' '+locs[var]+ \
                            ' could not be imported. Ensure that the class ' +\
                            'exists and is spelled properly.')
            elif var == 'rxn_expressions':
                self.parse_elementary_rxns(locs[var])
                setattr(self,var,locs[var])
            else:
                setattr(self,var,locs[var])


        if self.parser:
            if self.input_file:
                if getattr(self,'interaction_fitting_mode',None):
                    #automatically parse in "coverage" if fitting interaction params
                    self.parse_headers += ['coverage']
                self.parse() #Automatically parse in data.

        self.load_data_file()

        self.generate_echem_TS()

        self.verify()

    def load_data_file(self,overwrite=False):
        """
        Load in output data from external files.
        """
        if os.path.exists(self.data_file):
            pickled_data = pickle.load(open(self.data_file,'r'))
            for attr in pickled_data:
                if not overwrite:
                    if not getattr(self,attr,None): #don't over-write
                        setattr(self,attr,pickled_data[attr])
                else:
                    setattr(self,attr,pickled_data[attr])

    def parse(self,*args, **kwargs): #
        """
        .. todo:: __doc__
        """
        self.parser.parse(*args, **kwargs)

    def log(self,event,**kwargs):
        """
        .. todo:: __doc__
        """
        message = self._log_strings[event]
        loop, status = event.rsplit('_',1)
        kwargs['loop'] = loop
        kwargs['status'] = status

        if not hasattr(self,'_descriptors'):
            self._descriptors = [0,0]

        if 'pt' not in kwargs:
            kwargs['pt'] = self.print_point(self._descriptors,
                    self.descriptor_decimal_precision)
        if 'priority' not in kwargs:
            kwargs['priority'] = 0
        if 'n_iter' in kwargs:
            log_string = \
                    '${loop}_iteration_${n_iter}: ${status} - '+message
        else:
            log_string = '${loop}_evaluation: ${status} - '+message
        log_string = Template(log_string).substitute(kwargs)
        if log_string not in self._warned:
            lg = self._log_dict.get(kwargs['pt'],[])
            lg.append(log_string)
            self._log_dict[kwargs['pt']] = lg

            if self.verbose > kwargs['priority']:
                self._log_lines.append(log_string)
                print(log_string)
            if status == 'warning':
                self._warned.append(log_string)

    #Parsing and formatting functions

    def parse_elementary_rxns(self, equations): #
        """
        .. todo:: __doc__
        """
        elementary_rxns = []
        gas_names = []
        adsorbate_names = []
        transition_state_names = []
        site_names = []
        echem_transition_state_names = []
        for rxn_index, eq in enumerate(equations):
            #Replace separators with ' '
            regex = re.compile(regular_expressions['species_separator'][0])
            eq = regex.sub(' ',eq)
            state_dict = functions.match_regex(eq,
                    *regular_expressions['initial_transition_final_states'])
            rxn_list = []
            for key in ['initial_state','transition_state','final_state']:
                state_list = []
                state_str = state_dict[key]
                if state_str:
                    state_strings = [si for si in state_str.split() if si]
                    # echem transition state syntax parsing. generate echem transition state from TS like "^0.6eV_a"
                    if key == 'transition_state' and len(state_strings) == 1 and state_strings[0].startswith('^'):
                        try:
                            preamble, site = state_strings[0].split('_')
                            barrier = preamble.lstrip('^').rstrip('eV')
                            echem_TS_name = '-'.join(["echemTS", str(rxn_index), barrier]) + "_" + site
                            rxn_list.append([echem_TS_name])
                            echem_transition_state_names.append(echem_TS_name)
                            continue
                        except:
                            raise ValueError('improper specification of electrochemical transition state.  should be\
                                of the form "^0.6eV_a" for an 0.6 eV barrier on site a')
                    elif key == 'transition_state' and len(state_strings) == 1 and state_strings[0].startswith('echemTS'):
                        try:
                            echem_TS = state_strings[0]
                            preamble, site = echem_TS.split('_')
                            echem, i, barrier = preamble.split('-')
                            i = int(i)
                            assert(rxn_index == i)
                            barrier = float(barrier)
                            echem_transition_state_names.append(echem_TS)
                            continue
                        except:
                            raise ValueError('improper specification of electrochemical transition state.  should be\
                                of the form "echemTS-10-0.6_a" for an 0.6 eV barrier on site a for rxn 10')
                    for st in state_strings:
                        species_dict = functions.match_regex(st,
                                *regular_expressions['species_definition'])
                        if species_dict['stoichiometry'] == '':
                            species_dict['stoichiometry'] = 1
                        else:
                            try:
                                species_dict['stoichiometry'] = int(
                                        species_dict['stoichiometry'])
                            except:
                                raise ValueError('Non-integer stoichomtry: '+st)
                        if species_dict['site'] is None:
                            species_dict['site'] = self._default_site
                        site_names.append(species_dict['site'])
                        if species_dict['name'] != '*':
                            species_key = species_dict['name']+'_'+species_dict['site']
                        else:
                            species_key = species_dict['site']
                        if key in ['initial_state','final_state']:
                            if species_dict['name'] != '*':
                                if species_dict['site'] not in  self._gas_sites:
                                    adsorbate_names.append(species_key)
                                else:
                                    gas_names.append(species_key)
                        else:
                            if species_dict['name'] != '*':
                                if species_key not in adsorbate_names+gas_names:
                                    transition_state_names.append(species_key)
                        state_list += ([species_key]*species_dict['stoichiometry'])
                    rxn_list.append(state_list)
                elif key == 'transition_state':
                    pass # No transition-state
                else:
                    raise ValueError('Initial or final state is undefined: '+eq)
            elementary_rxns.append(rxn_list)
        gas_names = list(set(gas_names))
        adsorbate_names = list(set(adsorbate_names))
        transition_state_names = list(set(transition_state_names))
        site_names = tuple(set(site_names))
        for ts in transition_state_names:
            if ts in gas_names + adsorbate_names:
                transition_state_names.remove(ts)
        def sort_list(species_list):
            """
            .. todo:: __doc__
            """
            if not species_list:
                return ()
            new_list = []
            for sp in species_list:
                if '_' not in sp:
                    site = self._default_site
                else:
                    name,site = sp.rsplit('_',1)
                new_list.append([site,sp])
            new_list.sort()
            return zip(*new_list)[-1]

        self.gas_names = sort_list(gas_names)
        self.adsorbate_names = sort_list(adsorbate_names)
        self.transition_state_names = sort_list(transition_state_names)
        self.elementary_rxns = elementary_rxns
        self.site_names = site_names
        self.echem_transition_state_names = echem_transition_state_names

    def texify(self,ads): #
        """Generate LaTeX representation of an adsorbate.

        :param ads: Adsorbate short-hand.
        :type ads: str

        """
        sub_nums = [str(n) for n in range(2,15)]
        ads_def = ads
        ads = ads.replace('-','\mathrm{-}')
        if self.species_definitions[ads_def]['type'] == 'site':
            return '*_'+ads
        elif '_' in ads:
            adsN,site = ads.split('_')
        else:
            adsN = ads
            site = 's'
        for num in sub_nums:
            adsN = adsN.replace(num,'_{'+num+'}')
        if site == 'g':
            site = '(g)'
            tex_ads = adsN + '_{'+site+'}'
        else:
            site = '*_'+site
            tex_ads = adsN + '^{'+site+'}'
        tex_ads = tex_ads.replace('}_{','')
        return tex_ads

    def print_rxn(self,rxn,mode='latex',include_TS=True,print_out = False): #
        """
        .. todo:: __doc__
        """
        if mode == 'latex':
            def texify(ads):
                return self.texify(ads)
            def leftrightarrow():
                return ' \leftrightarrow '
            def rightarrow():
                return r' \rightarrow '

        elif mode == 'text':
            def texify(ads):
                return ads
            def leftrightarrow():
                return '\t<->\t'
            def rightarrow():
                return '\t->\t'

        IS = '+'.join([texify(a) for a in rxn[0]])
        if len(rxn) > 2:
            TS = '+'.join([texify(a) for a in rxn[1]])
        else:
            TS = None
        FS = '+'.join([texify(a) for a in rxn[-1]])
        if TS and include_TS == True:
            rxn_str = IS + leftrightarrow() + TS + rightarrow() + FS
        else:
            rxn_str = IS + leftrightarrow() + FS
        if print_out == True:
            print rxn_str
        return rxn_str

    @staticmethod
    def print_point(descriptors,n = 2):
        """
        .. todo:: __doc__
        """
        string = '['
        for d in descriptors:
            d = float(d)
            string += str('%5.'+str(n)+'f') % d
            string += ','
        string = string[:-1] #remove trailing comma
        string += ']'
        return string

    def model_summary(self,summary_file='summary.tex'):
        """
            Write summary of model into TeX file.

            :param summary_file: Filename where TeX summary of model is written.
            :type summary_file: str

        """
        if not hasattr(self,'summary_file'):
            self.summary_file = summary_file

        latex_template = catmap.data.templates['latex_summary']
        longtable_template = catmap.data.templates['latex_longtable']

        subs_dict = {}
        #reactions
        out_txt = ''
        out_txt += r'\section{Elementary Reactions}'
        for n,rx in enumerate(self.elementary_rxns):
            out_txt+= '\n'+r'\begin{equation}'
            out_txt+= '\n'+r'\label{reaction_'+str(n)+'}'
            out_txt+= '\n'+self.print_rxn(rx)
            out_txt+= '\n'+r'\end{equation}'

        #scaling relations
        out_txt += r'\section{Scaling Summary}'+'\n'
        out_txt += self.scaler.summary_text()

        #thermodynamics
        out_txt += r'\section{Thermodynamics Summary}'+'\n'
        out_txt += self.thermodynamics.summary_text()

        #solver
        out_txt += r'\section{Solver Summary}' + '\n'
        out_txt += self.solver.summary_text()

        #inputs
        out_txt += r'\section{Input Summary}' + '\n'
        longtable_txt = ''

        max_freqs = 3
        #This could be cleaned up a lot using templates...
        for spec in self.gas_names:
            energy = self.species_definitions[spec]['formation_energy']
            freqs = [str(round(v*1e3,1)) for v in self.species_definitions[spec].get('frequencies',[])]
            frequencies = '\parbox[t]{3cm}{'
            while freqs:
                freq_subset = [freqs.pop(0) for i in range(0,max_freqs) if freqs]
                frequencies += ', '.join(freq_subset)+r'\\'
            frequencies = frequencies[:-2] + '}'
            def cleanstring(string):
                return string.replace('\r','').replace('"','').replace("'",'')
            ads,site = spec.rsplit('_',1)
            facet = 'gas'
            surf = 'None'
            ref_tag = '_'.join([surf,facet,ads])
            reference = '\parbox[t]{4cm}{'+self.species_definitions[spec]['formation_energy_source']+'}'
            spec_tex = '$'+self.texify(spec)+'$'
            tabrow = '&'.join(
                    [spec_tex,'gas',str(
                        round(energy,2)),frequencies,reference]) + r'\\'
            longtable_txt += tabrow + '\n'

        for spec in self.adsorbate_names + self.transition_state_names:
            energy = self.parameter_dict[spec]
            for e,surf in zip(energy,self.surface_names):
                if e and e != '-':
                    e = str(round(e,2))
                    if self.species_definitions[spec].get('frequencies',[]):
                        freqs = [str(round(v*1e3,1))
                                for v in self.species_definitions[spec].get('frequencies',[])]
                        frequencies = '\parbox[t]{3cm}{'
                        while freqs:
                            freq_subset = [freqs.pop(0)
                                    for i in range(0,max_freqs) if freqs]
                            frequencies += ', '.join(freq_subset)+r'\\'
                        frequencies = frequencies[:-2] + '}'
                    else:
                        frequencies = ''
                    if '_' in spec:
                        ads,site = spec.rsplit('_',1)
                    else:
                        ads = spec
                        site = 's'
                    facet = self.species_definitions[site]['site_names']
                    if str(facet) != facet:
                        facet = ' or '.join(facet)
                    ref_tag = '_'.join([surf,facet,ads])
#                    if ref_tag in self.reference_dict:
#                        reference = '\parbox[t]{4cm}{'+ \
#                               r';\\'.join(
#                               [ref.replace('\r','').replace('"','').replace("'",'')
#                               for ref in [self.reference_dict[ref_tag]]])+'}'
#                        spec_tex = '$'+self.texify(spec)+'$'
#                        tabrow = '& '.join(
#                                [spec_tex,surf,e,frequencies,reference]) + r'\\'
#                        longtable_txt += tabrow + '\n'

        subs_dict['longtable_txt'] = longtable_txt

        longtable = Template(longtable_template).substitute(subs_dict)

        out_txt += longtable

        subs_dict['summary_txt'] = out_txt

        summary = Template(latex_template).substitute(subs_dict)

        f = open(self.summary_file,'w')
        f.write(summary)
        f.close()

    #Self checks, debugging, and code-structure related functions.

    def update(self,dict,override = False):
        """
        .. todo:: __doc__
        """
        if override == False:
            dict.update(self.__dict__)
            self.__dict__ = dict
        else:
            self.__dict__.update(dict)

    def compatibility_check(self):
        """
        Check that the reaction model has all required attributes.
        """
        total_dict = self._required
        for var in total_dict:
            val = getattr(self,var,None)
            if total_dict[var]:
                try:
                    right_type = total_dict[var](val)
                    setattr(self,var,right_type)
                except:
                    raise TypeError('The attribute '+str(var)+' must be '+
                            'compatible with the function ' +
                            str(total_dict[var]))
            if not hasattr(self,var):
                raise AttributeError('Reaction model must contain the '+
                        'attribute '+str(var))

    def verify(self):
        """
Run several consistency check on the model, such as :

- all gas ratios add to 1.
- all mass and site balances are fulfilled.
- prefactors are set in the correct format.
- a mapping resolution is set (the default 15 data points per descriptor axis).

        """


        #Check gas_ratios
        if hasattr(self,'gas_ratios') and self.gas_ratios:
            if sum(self.gas_ratios) != 1.0:
                print('Warning: gas_ratios do not sum to 1.'+
                     ' The total pressure will be '+str(sum(self.gas_ratios))
                     +' times the pressure indicated in output')

        #Check for mass/site balances on elementary equations
        def composition(state_list,type='atoms'):
            """
            .. todo:: __doc__
            """
            total_comp = {}
            for sp in state_list:
                if type == 'atoms':
                    comp = self.species_definitions[sp]['composition']
                    for atom in comp:
                        if atom not in total_comp:
                            total_comp[atom] = comp[atom]
                        else:
                            total_comp[atom] += comp[atom]
                elif type == 'sites':
                    n_sites = self.species_definitions[sp].get('n_sites',1)
                    if n_sites:
                        site = self.species_definitions[sp]['site']
                        if site not in total_comp:
                            total_comp[site] = n_sites
                        else:
                            total_comp[site] += n_sites
            for key in total_comp.keys():
                if total_comp[key] == 0:
                    del total_comp[key]
            return total_comp

        for rxn in self.elementary_rxns:
            if len(rxn) == 2:
                IS, FS = rxn
                TS = IS
            else:
                IS,TS,FS = rxn
            # ignore composition checking for echemTS stuff
            if 'echemTS' in TS[0]:
                TS = IS
            if composition(IS) == composition(TS) == composition(FS):
                pass
            else:
                raise ValueError('Mass balance is not satisfied for ' + \
                        self.print_rxn(rxn,mode='text'))
            if composition(IS,'sites') == composition(TS,'sites') == \
                    composition(FS,'sites'):
                pass
            else:
                raise ValueError('Site balance is not satisfied for ' + \
                        self.print_rxn(rxn,mode='text')+'. Make sure to '+\
                        'set n_sites in species_definitions for species '+\
                        'occupying more than 1 site')
        #Check that frequencies are defined if necessary

        #Check prefactor_list is in the right format
        default_prefactor = 'kB*T/h'
        default_prefactor_list = [default_prefactor] * len(self.elementary_rxns)

        if not self.prefactor_list:
            self.prefactor_list = default_prefactor_list
        elif isinstance(self.prefactor_list, list) and len(self.prefactor_list) == len(self.elementary_rxns):
            prefactor_list = []
            for prefactor in self.prefactor_list:
                if prefactor == None:
                    prefactor_list.append(default_prefactor)
                else:
                    prefactor_list.append(str(prefactor))
            self.prefactor_list = prefactor_list
        else:
            raise ValueError('prefactor_list must be None or a list ' + \
                'containing a prefactor for each elementary rxn.  The ' + \
                'elements of this list may contain None if you wish to use ' + \
                'the default prefactor of kB*T/h for that rxn')

        if not hasattr(self, 'resolution') or self.resolution is None:
            self.resolution = 15
            print("Info: set resolution to {self.resolution} as default.".format(**locals()))

    #Data manipulation and conversion

    def _header(self,exclude_outputs=[],re_parse=False):
        """Create a string which acts as a header for the log file."""
        header = ''
        for attr in self._classes:
            inst = getattr(self,attr)
            if inst:
                name = str(inst.__class__)
                junk,name = name.rsplit('.',1)
                name = name.replace("'",'').replace('>','')
                if re_parse:
                    header += attr + ' = ' + "'"+name+"'\n\n"
                elif attr != 'parser': #don't include a parser -> no reparsing
                    header += attr + ' = ' + "'"+name+"'\n\n"
                elif attr == 'parser': #don't include a parser -> no reparsing
                    header += attr + ' = ' + 'None\n\n'
            else:
                header += attr + ' = ' + 'None\n\n'
        exec(self._log_imports) #needed for testing evaluation of log file
        #load in pickled data at the beginning
        header += 'binary_data = ' + 'pickle.load(open("' + \
                                                self.data_file +'"))\n\n'
        header += 'locals().update(binary_data)\n\n'
        for attr in dir(self):
            if (not attr.startswith('_') and
                    not callable(getattr(self,attr)) and
                    attr not in self._classes):
                val = repr(getattr(self,attr))
                new_line = ''
                if attr not in self._pickle_attrs:
                    new_line = attr + ' = '+ val+ '\n\n'
                else:
                    new_line += attr + ' = binary_data[' + attr + ']'
                if new_line:
                    try:
                        if 'binary_data[' not in new_line:
                            #pickled data fails since it is saved after logfile
                            exec(new_line)
                            header+= new_line
                    except:
                        self.log('header_fail',save_txt = new_line)
        return header

    def _token(self):
        """Create a 'token' which uniquely identifies the model
        based on the user-input parameters. Two models with
        identical tokens should have identical solutions, although
        this is not guaranteed if an expert user changes some private
        attributes."""
        token = self._header()
        token = token.replace(' ','').replace('\n','')
        return token

    def retrieve_data(self,mapp,point,precision=2):
        """
        .. todo:: __doc__
        """
        if not mapp:
            return None
        n = precision
        if not hasattr(self,'_dict_maps'):
            self._dict_maps = {}
        if id(mapp) not in self._dict_maps:
            self._dict_maps[id(mapp)] = {}
            for pt,cvg in mapp:
                pt = tuple([round(v,n) for v in pt])
                self._dict_maps[id(mapp)][pt] = cvg
        newpt = tuple([round(v,n) for v in point])
        try:
            return self._dict_maps[id(mapp)][newpt]
        except KeyError:
            for pt,cvg in mapp:
                pt = tuple([round(v,n) for v in pt])
                self._dict_maps[id(mapp)][pt] = cvg
        return self._dict_maps[id(mapp)].get(newpt,None)

    @staticmethod
    def map_to_array(mapp,descriptor_ranges,resolution,
            log_interpolate=False,minval=None,maxval=None):
        """
        .. todo:: __doc__
        """
        desc_rngs = copy(descriptor_ranges)
        pts,datas = zip(*mapp)
        cols = zip(*datas)
        if len(pts[0]) == 1:
            xData = np.array(zip(*pts)[0])
            sorted_order = xData.argsort()
            maparray = np.zeros((resolution[0],len(datas[0])))  # resolution assumed to be [x, 1]
            x_range = desc_rngs[0]
            xi = np.linspace(x_range[0],x_range[1],resolution[0])
            datas = zip(*datas)
            for i,yData in enumerate(datas):
                yData = np.array(yData)
                y_sp = catmap.spline(xData[sorted_order],yData[sorted_order],k=1)
                yi = y_sp(xi)
                maparray[:,i] = yi

        elif len(pts[0]) == 2:
            xData,yData = zip(*pts)
            maparray = np.zeros((resolution[1],resolution[0],len(datas[0])))
            datas = zip(*datas)
            x_range,y_range = desc_rngs
            xi = np.linspace(*x_range+[resolution[0]])
            yi = np.linspace(*y_range+[resolution[1]])
            for i,Zdata in enumerate(datas):
                if minval:
                    Zdata = np.array([max(zn,minval) for zn in Zdata])
                if maxval:
                    Zdata = np.array([min(zn,maxval) for zn in Zdata])
                if log_interpolate == True:
                    Zdata_log = np.array(
                            [np.log(abs(float(zn))) for zn in Zdata])
                    z_sign = np.sign(griddata(xData,yData,Zdata,xi,yi,interp='linear'))
                    z_num = griddata(xData,yData,Zdata_log,xi,yi,interp='linear')
                    zi = np.exp(z_num)*z_sign
                else:
                    zi = griddata(xData,yData,Zdata,xi,yi,interp='linear')
                maparray[:,:,i] = zi
        return maparray

    @staticmethod
    def array_to_map(array,descriptor_ranges,resolution):
        """
        .. todo:: __doc__
        """
        dim = len(array.shape)
        xy = []
        ij = []

        def make_ptlist(descriptor_ranges,resolution,pt_list=[],ij_list=[]):
            if descriptor_ranges:
                rng = descriptor_ranges.pop(0)
                vals = list(np.linspace(*rng+[resolution]))
                if not pt_list:
                    new_pts = [[v] for v in vals]
                    new_ij = [[i] for i in range(0,resolution)]
                else:
                    new_pts = []
                    new_ij = []
                    for ij,pt in zip(ij_list,pt_list):
                        for i,v in enumerate(vals):
                            new_pts.append(pt+[v])
                            new_ij.append(ij+[i])
                ij_list, pt_list = \
                        make_ptlist(descriptor_ranges,resolution,new_pts,new_ij)
                return ij_list, pt_list
            else:
                return ij_list, pt_list


        rngs = copy(descriptor_ranges)
        pt_idxs, pts = make_ptlist(rngs,resolution)

        def get_next_dim(array,xy_vec):
            if xy_vec:
                new_array = array[xy_vec.pop(-1)]
                new_array = get_next_dim(new_array,xy_vec)
                return new_array
            else:
                return list(array)

        datas = []
        for ij in pt_idxs:
            data = get_next_dim(array,list(ij))
            datas.append(data)

        mapp = zip(*[pts,datas])
        return mapp

    #Commonly used convenience functions

    @staticmethod
    def same_rxn(rxn1,rxn2):
        """Determine if two reactions *rxn1* and *rxn2* are identical.

        """

        def same_state(state1, state2):
            state1 = [s.replace('_s','') for s in state1]
            state2 = [s.replace('_s','') for s in state2]
            if sorted(state1) == sorted(state2):
                return True
            else:
                return False

        def get_IS_TS_FS(rxn):
            IS = rxn[0]
            FS = rxn[-1]
            if len(rxn) == 2:
                TS = []
            elif len(rxn) > 2:
                TS = rxn[2]
                if same_state(TS,IS) or same_state(TS,FS):
                    TS = []
            return IS, TS, FS
        IS1, TS1, FS1 = get_IS_TS_FS(rxn1)
        IS2, TS2, FS2 = get_IS_TS_FS(rxn2)
        if (
                same_state(IS1,IS2) and
                same_state(TS1, TS2) and
                same_state(FS1, FS2)
                    ):
            return True
        else:
            return False

    @staticmethod
    def reverse_rxn(rxn):
        """
        .. todo:: __doc__
        """
        return [rxn[-1],rxn[1],rxn[0]]

    def get_rxn_energy(self,rxn,energy_dict):
        """
        .. todo:: __doc__
        """

        IS = rxn[0]
        FS = rxn[-1]
        if len(rxn) <= 2:
            TS = None
        else:
            TS = rxn[1]
        E_IS = self.get_state_energy(IS,energy_dict)
        E_FS = self.get_state_energy(FS,energy_dict)
        if TS:
            E_TS = self.get_state_energy(TS,energy_dict)
        else:
            E_TS = max(E_IS,E_FS)

        dE = E_FS - E_IS
        E_a = max(0,(E_TS - E_IS),(E_FS - E_IS))
        return dE, E_a

    def get_state_energy(self,rxn_state,energy_dict):
        """
        .. todo:: __doc__
        """
        energy = 0
        for species in rxn_state:
            if species in energy_dict:
                energy += energy_dict[species]
            elif self.species_definitions[species]['type'] == 'site':
                energy += self._site_energies[self.site_names.index(species)]
        return energy

    def adsorption_to_reaction_energies(self,free_energy_dict):
        "Convert adsorption energies to reaction energies/barriers."
        Grxn_Ga = []
        for rxn in self.elementary_rxns:
            dG, G_a = self.get_rxn_energy(rxn,free_energy_dict)
            Grxn_Ga.append([dG,G_a])
        return Grxn_Ga

    def make_standalone(self, standalone_script='stand_alone.py'):
        """Create a stand alone script containing the current model.

           :param standalone_script: The name of the file where the standalone script is created [stand_alone.py].
           :type standalone_script: str


        """
        txt = ''
        for func in self._function_strings:
            txt+= self._function_strings[func]
            txt += '\n'*3
        f = open(standalone_script, 'w')
        f.write(txt)
        f.close()

    def generate_echem_TS(self):
        """generates fake transition state species from self.echem_transition_state_names
        and populates self.species_definitions with them.
        """
        for echem_TS in self.echem_transition_state_names:
            preamble, site = echem_TS.split('_')
            echem, rxn_index, barrier = preamble.split('-')
            rxn_index = int(rxn_index)
            barrier = float(barrier)
            self.species_definitions[echem_TS] = {}
            self.species_definitions[echem_TS]['type'] = 'transition_state'
            self.species_definitions[echem_TS]['site'] = site
            self.species_definitions[echem_TS]['formation_energy_source'] = ["Generated by CatMAP"] * len(self.surface_names)
            self.species_definitions[echem_TS]['formation_energy'] = [0.] * len(self.surface_names)
            self.species_definitions[echem_TS]['frequencies'] = []
            self.species_definitions[echem_TS]['name'] = preamble
            self.species_definitions[echem_TS]['n_sites'] = 1  # Someone may want to change this to be user-specified at some point

            # look up what IS and FS of rxn_index are
            regex = re.compile(regular_expressions['species_separator'][0])
            eq = regex.sub(' ',self.rxn_expressions[rxn_index])
            state_dict = functions.match_regex(eq,
                *regular_expressions['initial_transition_final_states'])
            IS_species = [si for si in state_dict['initial_state'].split(' ') if si]
            FS_species = [si for si in state_dict['final_state'].split(' ') if si]
            # assume composition is already balanced - set TS composition to IS composition
            total_composition = {}
            for species in IS_species:
                functions.add_dict_in_place(total_composition, self.species_definitions[species]['composition'])
            self.species_definitions[echem_TS]['composition'] = total_composition

        # add echem TSs to regular TSes - this might be more trouble than it's worth
        self.transition_state_names += tuple(self.echem_transition_state_names)
