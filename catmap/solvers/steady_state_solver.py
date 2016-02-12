from solver_base import *
from mean_field_solver import *
from catmap import string2symbols
try:
    from scipy.optimize import fmin_powell as fmin
except ImportError:
    fmin = None

from catmap.functions import numerical_jacobian
import math
from string import Template
import random
import re

class SteadyStateSolver(MeanFieldSolver):

    def __init__(self,reaction_model=ReactionModel()):
        MeanFieldSolver.__init__(self,reaction_model)
        defaults = dict(
                max_rootfinding_iterations = 50,
                internally_constrain_coverages = True,
                residual_threshold = 0.9,
                analytical_jacobian = True,
                optimize_analytical_expressions = False,
                )
        self._rxm.update(defaults)
        self._rate_constant_memoize = {}
        self._steady_state_memoize = {}
        self._required = {'max_rootfinding_iterations':int,
                          'internally_constrain_coverages':None,
                          'residual_threshold':float,
                          'analytical_jacobian':bool,
                          }
        self._log_strings = {'rootfinding_fail':
                            "stagnated or diverging (residual = ${resid})",
                            'rootfinding_maxiter':
                            "exceeded maximum iterations (residual = ${resid})",
                            'rootfinding_cancel':
                  "stationary point or singular jacobian (residual = ${resid})",
                            'rootfinding_success':
                            "found solution at point ${pt}",
                            'rootfinding_status':
                            "converging (residual = ${resid})"}

    def get_rate_constants(self,rxn_parameters,coverages):
        """Return rate constants for given sequence of reaction parameters and coverages.

        :param rxn_parameters: Sequence of reaction parameters.
        :type rxn_parameters: [float]
        :param coverages: Sequence of coverages.
        :type coverages: [float]

        """

        if self.adsorbate_interaction_model not in [None,'ideal']:
            memo = tuple(rxn_parameters) + tuple(coverages) + tuple(self._gas_energies)
        else:
            memo = tuple(rxn_parameters) + tuple(self._gas_energies)
        if memo in self._rate_constant_memoize:
            kf, kr = self._rate_constant_memoize[memo]
            self._kf = kf
            self._kr = kr
            return kf+kr
        kfs, krs, dkfs, dkrs = self.rate_constants(rxn_parameters,coverages,
            self._gas_energies,self._site_energies,
            self.temperature,self.interaction_response_function,
            self._mpfloat,self._matrix,self._math.exp)
        self._kf = kfs
        self._kr = krs
        self._rate_constant_memoize[memo] = [kfs,krs]
        return kfs + krs

    def get_coverage(self,rxn_parameters,c0=None,findrootArgs={}):
        """Return coverages for given reaction parameters and coverage constraints.

            :param rxn_parameters: Sequence of rxn_parameters
            :type rxn_parameters: [float]
            :param c0: Coverage constraints.
            :type c0: **TODO**
            :param findrootArgs: *deprecated*
        """
        if self.adsorbate_interaction_model in [None,'ideal'] or self.interaction_strength == 0:
            return self.get_ideal_coverages(rxn_parameters,c0,True,findrootArgs)
        else:
            return self.get_interacting_coverages(rxn_parameters,c0,1.0,findrootArgs)

    def get_steady_state_coverage(self,rxn_parameters,steady_state_fn, jacobian_fn,
            c0=None,findrootArgs={}):
        """Return steady-state coverages using catmap.solvers.solver_base.NewtonRoot .

        :param rxn_parameters: Sequence of reaction parameters.
        :type rxn_parameters: [float]
        :param steady_state_fn: **TODO**
        :type steady_state_fn: **TODO**
        :param jacobian_fn: **TODO**
        :type jacobian_fn: **TODO**
        :param c0: Coverage constraints
        :type c0: **TODO**
        :param findrootArgs: *deprecated*

        """

        n_tot = len(self.adsorbate_names) + len(self.transition_state_names)
        if c0 is None:
            raise ValueError("No initial coverage supplied. Mapper must supply initial guess")

        c0 = self.constrain_coverages(c0)
        self.steady_state_function = steady_state_fn
        self.steady_state_jacobian = jacobian_fn
        self._coverage = [self._mpfloat(ci) for ci in c0]
        self._rxn_parameters = rxn_parameters

        #Enter root finding algorithm
        f = steady_state_fn
        f_resid = lambda x: self.get_residual(x,True,False)
        norm = self._math.infnorm

        if self.internally_constrain_coverages == True:
            constraint = self.constrain_coverages
        else:
            constraint = lambda x: x
        solver = NewtonRoot

        if f_resid(c0) <= self.tolerance:
            self._coverage = c0
            return c0

        solver_kwargs = dict(
                norm = norm,
                verbose = self.verbose,
                constraint = constraint,
                )

        if self.analytical_jacobian == True:
            solver_kwargs['J'] = jacobian_fn
        else:
            def J(x):
                return numerical_jacobian(f,x,self._matrix)
            solver_kwargs['J'] = J


        iterations = solver(f,c0, self._matrix, self._mpfloat,
                            self._Axb_solver, **solver_kwargs)
        old_error = 1e99
        coverages = None
        maxiter = self.max_rootfinding_iterations
        iterations.maxiter = maxiter
        i = 0
        x = c0
        for x,error in iterations:
            self.log('rootfinding_status',
                    n_iter=i,
                    resid=float(error),
                    priority=1)
            i+=1
            if error < self.tolerance:
                if f_resid(x) < self.tolerance:
                    coverages = self.constrain_coverages(x)
                    self.log('rootfinding_success',
                            n_iter = i,
                            priority = 1)
                    break
                else:
                    x = self.constrain_coverages(x)
                    error = f_resid(x)

            elif error >= self.residual_threshold*old_error:
                self.log('rootfinding_fail',
                        n_iter=i,
                        resid = float(f_resid(x)),
                        )
                raise ValueError('Stagnated or diverging residual (resid='+\
                        str(float(f_resid(x)))+')')

            old_error = error

            if i >= maxiter:
                self.log('rootfinding_maxiter',
                        n_iter=i,
                        resid = float(f_resid(x)))
                raise ValueError('Out of iterations (resid='+\
                        str(float(f_resid(x)))+')')
            self._coverage = x
            self._error = error

        if coverages:
            self._coverage = [c for c in coverages]
            return [c for c in coverages]
        else:
            if f_resid(x) < self.tolerance:
                coverages = self.constrain_coverages(x)
            else:
                self.log('rootfinding_cancel',
                        n_iter=i,
                        resid=float(f_resid(x)))
                raise ValueError('Solver cancellation. (resid='+\
                        str(float(f_resid(x)))+')')

    def get_ideal_coverages(self,rxn_parameters,c0=None,
            refresh_rate_constants=True,findrootArgs={}):
        """Return

            :TODO:

        """
        if refresh_rate_constants:
            self.get_rate_constants(rxn_parameters,[0]*len(self.adsorbate_names))
        return self.get_steady_state_coverage(rxn_parameters,self.ideal_steady_state_function,
                self.ideal_steady_state_jacobian,c0,findrootArgs)

    def get_interacting_coverages(self,rxn_parameters,c0=None,
            interaction_strength=1.0,findrootArgs={}):
        """:TODO:

        """

        #weight interaction parameters - useful for using non-interacting systems as guess
        n_tot = len(self.adsorbate_names +self.transition_state_names)
        rxn_parameters = list(rxn_parameters)
        if interaction_strength is not None:
            if len(rxn_parameters) == n_tot+n_tot**2:
                rxn_parameters = rxn_parameters[:n_tot] + [pi*interaction_strength for pi in rxn_parameters[-n_tot**2:]]
            elif len(rxn_parameters) == n_tot and interaction_strength == 0:
                pass
            else:
                raise ValueError('System does not have enough parameters for interactions')
        cvgs =  self.get_steady_state_coverage(rxn_parameters,self.interacting_steady_state_function,
            self.interacting_steady_state_jacobian,c0,findrootArgs)
        return cvgs

    def bisect_interaction_strength(self,rxn_parameters,valid_strength,valid_coverages,target_strength,max_bisections,findrootArgs={}):
        """
            :TODO:
        """
        n_tot = len(self.adsorbate_names+self.transition_state_names)
        bisect_iter = 0
        n_bisects = 0
        while n_bisects <= max_bisections:
            new_strength = valid_strength + (target_strength-valid_strength)/(2**n_bisects)
            new_params = rxn_parameters[:n_tot] + [pi*new_strength for pi in rxn_parameters[-n_tot**2:]]
            try:
                valid_coverages =  self.get_steady_state_coverage(new_params,self.interacting_steady_state_function,
                    self.interacting_steady_state_jacobian,valid_coverages,findrootArgs)
                print 'Successfully bisected with strength: ',new_strength
                valid_strength = new_strength
                n_bisects = 0
                if valid_strength > 0.95*target_strength:
                    return valid_coverages
            except ValueError:
                print 'Failed to bisect with strength: ',new_strength
                n_bisects += 1
        return valid_coverages


    def get_initial_coverage(self,rxn_parameters):
        """Return coverages based on probabilties according to the Boltzmann distribution
        and the adsorption energies for a given sequence of rxn_parameters.

        :param rxn_parameters: Sequence of reaction parameters
        :type rxn_parameter: [float]

        """

        energy_dict = {}
        for ads,E in zip(self.adsorbate_names,rxn_parameters):
            energy_dict[ads] = E
        for gas,E in zip(self.gas_names,self._gas_energies):
            energy_dict[gas] = E
        if not self.atomic_reservoir_dict:
            #check all possibilities, return min residual
            min_resid = 1e99
            boltz_cvgs = [[0]*len(self.adsorbate_names)] #include empty coverage as possible guess
            for ref_dict in self.atomic_reservoir_list:
                self.atomic_reservoir_dict = ref_dict
                cvgs = self.thermodynamics.boltzmann_coverages(energy_dict)
                boltz_cvgs.append(cvgs)

        else:
            boltz_cvgs = [self.thermodynamics.boltzmann_coverages(energy_dict)]

        return boltz_cvgs

    def get_residual(self, coverages,
            validate_coverages = True, refresh_rate_constants = True):
        """
            :TODO:
        """

        if validate_coverages == True:
            coverages = self.constrain_coverages(coverages)
        self._coverage = coverages
        if refresh_rate_constants == True:
            self._rxn_parameters = self.scaler.get_rxn_parameters(
                    self._descriptors)
            self.get_rate_constants(self._rxn_parameters,coverages)
#        cvg_rates = self.steady_state_function(None)
        cvg_rates = self.steady_state_function(coverages)
        residual = max([abs(r) for r in cvg_rates])
        return residual


    def interacting_steady_state_function(self, coverages):
        """
            :TODO:
        """

        memo = tuple(self._rxn_parameters) + tuple(self._gas_energies) + \
                tuple(self._site_energies) + tuple(coverages) + tuple(self.gas_pressures+[self.temperature])
        if memo in self._steady_state_memoize:
            return self._steady_state_memoize[memo]
        else:
            c = self.interacting_mean_field_steady_state(
                    self._rxn_parameters,coverages,self.gas_pressures,
                    self._gas_energies, self._site_energies,
                    self.temperature,self.interaction_response_function,
                    self._mpfloat, self._matrix,self._math.exp)
            self._steady_state_memoize[memo] = c
            return c

    def ideal_steady_state_function(self,coverages):
        """
            :TODO:
        """

        memo = tuple(self._kf) + tuple(self._kr) + tuple(coverages) + tuple(self.gas_pressures+[self.temperature])
        if memo in self._steady_state_memoize:
            return self._steady_state_memoize[memo]
        else:
            c = self.ideal_mean_field_steady_state(
                    self._kf,self._kr,coverages,self.gas_pressures,
                    self._mpfloat, self._matrix)
            self._steady_state_memoize[memo] = c
            return c

    def interacting_steady_state_jacobian(self,coverages):
        """
            :TODO:
        """

        J = self.interacting_mean_field_jacobian(
                self._rxn_parameters,coverages,self.gas_pressures,
                self._gas_energies,self._site_energies,
                self.temperature,self.interaction_response_function,
                self._mpfloat, self._matrix,self._math.exp)
        return J

    def ideal_steady_state_jacobian(self,coverages):
        """
            :TODO:
        """
        J = self.ideal_mean_field_jacobian(
                self._kf,self._kr,coverages,self.gas_pressures,
                self._mpfloat, self._matrix)
        return J

    def constrain_coverages(self,cvgs):
        """
            :TODO:
        """

        min_cvg = self._mpfloat(10**(-(self.decimal_precision)))
        cvgs = self.constrain_coverage_function(list(cvgs),self._mpfloat,min_cvg)
        return cvgs

    def compile(self):
        """
            :TODO:
        """
        if not self._compiled:
            intermediate_subs = {}
            self._function_substitutions.update(
                    self.substitutions_dict())

            #make 2 versions of rate-constants function
            energy_expressions_noderivs = '\n    '.join(self.reaction_energy_equations(adsorbate_interactions = False))
            energy_expressions_derivs = '\n    '.join(self.reaction_energy_equations(adsorbate_interactions = True))

            templates['rate_constants_no_derivatives'] = Template(templates['rate_constants']).safe_substitute({'elementary_step_energetics':energy_expressions_noderivs})
            templates['rate_constants_with_derivatives'] = Template(templates['rate_constants']).safe_substitute({'elementary_step_energetics':energy_expressions_derivs})

            #make steady-state expressions
            ss_eqs = self.rate_equations()
            self._function_substitutions['steady_state_expressions'] = '\n    '.join(ss_eqs)

            #make jacobian expressions
            jac_eqs = self.jacobian_equations(adsorbate_interactions=True)
            self._function_substitutions['jacobian_expressions'] = '\n    '.join(jac_eqs)
            jac_eqs_nd = self.jacobian_equations(adsorbate_interactions=False)
            self._function_substitutions['jacobian_expressions_no_derivatives'] = '\n    '.join(jac_eqs_nd)

            def indent_string(string,levels):
                lines = string.split('\n')
                indention = '\n'+'    '*levels
                return indention.join(lines)

            #pre-substitute the interaction function into rate_constants (needed because its nested 2 levels)
            indented = indent_string(templates[self.adsorbate_interaction_model+'_interaction_function'],1)
            indented = Template(indented).substitute(self._function_substitutions)
            self._function_substitutions['interaction_function'] = indented

            for f in ['rate_constants_no_derivatives','rate_constants_with_derivatives']:
                templates[f] = Template(templates[f]).substitute(self._function_substitutions)

            #indent rate_constant functions because they are nested 1 level
            indented_funcs = [
                    ['rate_constants_no_derivatives',''],
                    ['rate_constants_with_derivatives',''],
                    ]

            for func,tempname in indented_funcs:
                if not tempname:
                    tempname = func
                self._function_substitutions[func] = indent_string(templates[tempname],1)
            compiled_funcs = [
                    ['rate_constants','rate_constants_with_derivatives'],
                    ['interaction_function',self.adsorbate_interaction_model+'_interaction_function'],
                    ['interacting_mean_field_steady_state',''],
                    ['ideal_mean_field_steady_state',''],
                    ['interacting_mean_field_jacobian',''],
                    ['ideal_mean_field_jacobian',''],
                    ['constrain_coverage_function','constrain_coverages'],
                    ['elementary_rates',''],
                    ]

            for func,tempname in compiled_funcs:
                if not tempname:
                    tempname = func
                self._function_templates[func] = templates[tempname]
            self.generate_static_functions()

            if self.optimize_analytical_expressions:
                test_theta = [self._mpfloat(random.random()) for a in self.adsorbate_names]
                test_params = [self._mpfloat(random.random()) for a in self.adsorbate_names+self.transition_state_names]
                test_kfs = [self._mpfloat(random.random()) for a in self.elementary_rxns]
                test_krs = [self._mpfloat(random.random()) for a in self.elementary_rxns]
                test_p = [self._mpfloat(random.random()) for a in self.gas_names]
                test_gas_E = [self._mpfloat(random.random()) for a in self.gas_names]
                test_site_E = [self._mpfloat(random.random()) for a in self.site_names]
                test_T = 500
                test_smearing = 0.02
                arg_dict = {
                        'interacting_mean_field_steady_state':[
                            test_params,test_theta,test_p,test_gas_E,test_site_E,
                            test_T,self.interaction_response_function,
                            self._mpfloat,self._matrix,self._math.exp],
                        'ideal_mean_field_steady_state':[
                            test_kfs, test_krs, test_theta, test_p,
                            self._mpfloat, self._matrix],
                        'interacting_mean_field_jacobian':[
                            test_params,test_theta,test_p,test_gas_E,test_site_E,
                            test_T,self.interaction_response_function,
                            self._mpfloat, self._matrix, self._math.exp],
                        'ideal_mean_field_jacobian':[
                            test_kfs, test_krs, test_theta, test_p,
                            self._mpfloat, self._matrix]
                        }

                for func in arg_dict:
                    args= arg_dict[func]
                    old_str = self._function_strings[func]
                    insertion_line = None

                    for i,li in enumerate(old_str.split('\n')):
                        if '    s['+str(len(self.site_names)-1)+']' in li:
                            insertion_line = i+1

                    if insertion_line is None:
                        raise ValueError('Could not find line for inserting optimizations')

                    #optimize string
                    func_string = self.optimize_analytical_function(func,old_str,
                            insertion_line,1,*args)

                    #re-compile optimized function
                    self._function_strings[func] = func_string
                    locs = {}
                    exec func_string in globals(), locs
                    setattr(self,func,locs[func])

            self._compiled = True

            if self.adsorbate_interaction_model not in [None,'ideal']:
                self.steady_state_function = self.interacting_steady_state_function
            else:
                self.steady_state_function = self.ideal_steady_state_function

    def optimize_analytical_function(self,func_name,func_string,insertion_line,indention_level,*test_args):
        """Replace some common multiplication terms to speed up functions.

        """

        locs = {}
        exec func_string in globals(), locs
        unoptimized = locs[func_name]

        #replace common multiplications with substitution
        mult_regex = '((?:(?:k(?:f|r)\[\d\]|theta\[\d\]|p\[\d\]|s\[\d\])\*?)+ +)'
        opt_strs   = []
        match = re.findall(mult_regex,func_string)
        matches = []
        for m in match:
            if m.count('*') and func_string.count(m) > 1:
                matches.append((-1*m.count('*'),func_string.count(m),m))
        matches = list(set(matches))
        matches = sorted(matches)
        for i,m in enumerate(matches):
            opt_strs.append(m[-1].strip())
            func_string = func_string.replace(m[-1].strip(),'subs['+str(i)+'] ')

        sub_dict = {
                    '':r'( +1\*|\-1\*\-1\*| +0 +\+| +1\.0\*)+',
                    ' + ':r'( +\- *\-1\*| +\+ *1\*| +\+ *\+ +)+',
                    ' - ':r'( +\+ *\-1\*| +\- *1\*| +\+ *\- +| +0 +\-|\-1\*|\-1\*1\*)',
                    }

        opt_str = '\n'+'    '*indention_level+'subs = [0]*'+str(len(opt_strs))
        for i,op in enumerate(opt_strs):
            opt_str += '\n'+'    '*indention_level+ 'subs['+str(i)+'] = '+op

        func_lines = func_string.split('\n')
        func_lines.insert(insertion_line, opt_str)
        func_string = '\n'.join(func_lines)

        mtch = 1
        while mtch:
            new_fn = func_string
            mtch = None
            for sub in sub_dict:
                rgx = sub_dict[sub]
                mtch = re.findall(rgx,func_string)
                for m in list(set(mtch)):
                    func_string = func_string.replace(m,sub)

        locs = {}
        exec func_string in globals(), locs

        optimized = locs[func_name]
        delta = np.array((optimized(*test_args) - unoptimized(*test_args)).tolist()).max()

        if delta > 10**(-(self.decimal_precision-1)):
            print('Warning: Function optimization failed for '+\
                    func_name+'. Using unoptimized functions')
            func_string = '\n'.join(func_lines)
            return func_string
        else:
            return func_string
