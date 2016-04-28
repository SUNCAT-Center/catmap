import catmap
from catmap.model import ReactionModel
from catmap import ReactionModelWrapper
import numpy as np
import mpmath as mp
from ase.atoms import string2symbols

class SolverBase(ReactionModelWrapper):
    def __init__(self,reaction_model=ReactionModel()):
        """
        Class for `solving' for equilibrium coverages and rates as a 
        function of reaction parameters. This class acts as a base class 
        to be inherited by other solver classes, but is not 
        functional on its own. 

        rxn_parameters: list of necessary parameters to solve the kinetic 
        system. This will usually be populated by the scaler.

        A functional derived solver class must also contain the methods:

        get_coverage(): a function which returns coverages for each 
            adsorbate as a list [cvg_ads1,cvg_ads2,...]

        get_rate(): a function which returns steady-state reaction 
            rates for each elementary step as a list [rate_rxn1,rate_rxn2,...]

        get_residual(): a function for computing the norm of the residual. This
            is the condition which will be minimized to reach steady-state.

        compile(): a function to set-up/compile the solver.

        """

        self._rxm = reaction_model
        self._compiled = False

    def set_output_attrs(self,rxn_parameters):
        """
        :param rxn_parameters: Reaction parameters.
        :type rxn_parameters: list
        """
        if True in [v in self.mapper._solver_output 
                for v in self.output_variables]:
            cvgs = self._coverage
            self._coverage = list(self.solver.get_coverage( 
                    rxn_parameters,c0=cvgs)) #verify coverage

            self._rate = list(self.solver.get_rate(rxn_parameters,
                    coverages=self._coverage))

        if (
                'turnover_frequency' in self.output_variables or
                'production_rate' in self.output_variables or
                'consumption_rate' in self.output_variables or
                'rxn_direction' in self.output_variables):

            self._turnover_frequency = self.get_turnover_frequency(
                    rxn_parameters)
        
        if 'selectivity' in self.output_variables:
            self._selectivity = self.get_selectivity(rxn_parameters)
            self.output_labels['selectivity'] = self.gas_names

        if 'carbon_selectivity' in self.output_variables:
            weights = []
            for g in self.gas_names:
                name,site = g.split('_')
                weight = string2symbols(name).count('C')
                weights.append(weight)
            self._carbon_selectivity = self.get_selectivity(rxn_parameters,weights=weights)
            self.output_labels['carbon_selectivity'] = self.gas_names


        if 'rate_control' in self.output_variables:
            self._rate_control = self.get_rate_control(rxn_parameters)
            self.output_labels['rate_control'] = [self.gas_names,self.parameter_names]

        if 'selectivity_control' in self.output_variables:
            self._selectivity_control = self.get_selectivity_control(rxn_parameters)
            self.output_labels['selectivity_control'] = [self.gas_names,self.parameter_names]

        if 'rxn_order' in self.output_variables:
            self._rxn_order = self.get_rxn_order(rxn_parameters)
            self.output_labels['rxn_order'] = [self.gas_names,self.gas_names]

        if 'apparent_activation_energy' in self.output_variables:
                self._apparent_activation_energy = self.get_apparent_activation_energy(rxn_parameters)
                self.output_labels['apparent_activation_energy'] = self.gas_names

        if 'interacting_energy' in self.output_variables:
            if self.adsorbate_interaction_model in [None,'ideal']:
                self._interacting_energy = rxn_parameters
            else:
                self._interacting_energy = self.get_interacting_energies(rxn_parameters)
            self.output_labels['interacting_energy'] = self.adsorbate_names+self.transition_state_names

        if 'directional_rates' in self.output_variables:
            self._directional_rates = self.get_directional_rates(rxn_parameters)
            self.output_labels['directional_rates'] = [str(rxn) + ' forward' for rxn in self.elementary_rxns] + \
                [str(rxn) + ' reverse' for rxn in self.elementary_rxns]

        if 'turnover_frequency' in self.output_variables:
            self.output_labels['turnover_frequency'] = self.gas_names

        for out in self.output_variables:
            if out == 'production_rate':
                self._production_rate = [max(0,r) 
                        for r in self._turnover_frequency]
                self.output_labels['production_rate'] = self.gas_names
            if out == 'consumption_rate':
                self._consumption_rate = [max(0,-r) 
                        for r in self._turnover_frequency]
                self.output_labels['consumption_rate'] = self.gas_names
            if out == 'forward_rate':
                self._forward_rate = [max(0,r) 
                        for r in self._rate]
                self.output_labels['forward_rate'] = self.elementary_rxns
            if out == 'reverse_rate':
                self._reverse_rate = [max(0,-r) 
                        for r in self._rate]
                self.output_labels['reverse_rate'] = self.elementary_rxns
            if out == 'rxn_direction':
                self._rxn_direction = [np.sign(r) 
                        for r in self._turnover_frequency]
                self.output_labels['rxn_direction'] = self.elementary_rxns
            if out == 'rate_constant':
                self._rate_constant = list(self._kf)+list(self._kr)
                self.output_labels['rate_constant'] = self.elementary_rxns + self.elementary_rxns
            if out == 'forward_rate_constant':
                self._forward_rate_constant = list(self._kf)
                self.output_labels['forward_rate_constant'] = self.elementary_rxns
            if out == 'reverse_rate_constant':
                self._reverse_rate_constant = list(self._kr)
                self.output_labels['reverse_rate_constant'] = self.elementary_rxns
            if out == 'equilibrium_constant':
                self._equilibrium_constant = [kf/kr 
                        for kf,kr in zip(self._kf,self._kr)]
                self.output_labels['equilibrium_constant'] = self.elementary_rxns


class NewtonRoot:
    """
    Hacked from MDNewton in mpmath/calculus/optimization.py in order
    to allow for constraints on the solution.

    Find the root of a vector function numerically using Newton's method.

    f is a vector function representing a nonlinear equation system.

    x0 is the starting point close to the root.

    J is a function returning the Jacobian matrix for a point.

    Supports overdetermined systems.

    Use the 'norm' keyword to specify which norm to use. Defaults to max-norm.
    The function to calculate the Jacobian matrix can be given using the
    keyword 'J'. Otherwise it will be calculated numerically.

    Please note that this method converges only locally. Especially for high-
    dimensional systems it is not trivial to find a good starting point being
    close enough to the root.

    It is recommended to use a faster, low-precision solver from SciPy [1] or
    OpenOpt [2] to get an initial guess. Afterwards you can use this method for
    root-polishing to any precision.

    [1] http://scipy.org

    [2] http://openopt.org
    """

    maxsteps = 10

    def __init__(self, f, x0, matrix, mpfloat, Axb_solver, **kwargs):
        self._matrix = matrix
        self._mpfloat = mpfloat
        self._Axb = Axb_solver
        self.f = f
        self.x0 = x0
        if 'J' in kwargs:
            self.J = kwargs['J']

            #the following is useful for debugging/benchmarking
            #analytical derivatives, and should be commented out
            #for any production code.
#            import time
#            def J(x): #Use this to confirm the analytical jacobian is correct
#                t0 = time.time()
#                analytical = kwargs['J'](x)
#                t_a = time.time() - t0
#                t0 = time.time()
#                numerical = catmap.functions.numerical_jacobian(f,x,matrix,1e-300)
#                t_n = time.time() - t0
#                error = analytical - numerical
#                error = error.tolist()
#                max_error = -1
#                max_pos = None
#                for i,ei in enumerate(error):
#                    for j,ej in enumerate(ei):
#                        if abs(ej) > 1e-10:
#                            print 'big error', ej, [i,j]
#                            pass
#                        if abs(ej) > max_error:
#                            max_error = abs(ej)
#                            max_pos = [i,j]
#                print 'max_error', max_error, max_pos
#                print 't_analytic/t_numerical', t_a/t_n
#                return numerical
#            self.J = J

#            def J_numerical(x): #Use this to confirm the analytical jacobian is correct
#                numerical = catmap.functions.numerical_jacobian(f,x,matrix,1e-50)
#                return numerical
#            self.J = J_numerical


        else:
            raise ValueError('No method for estimating Jacobian.')
        if 'constraint' in kwargs:
            self.constraint = kwargs['constraint']
        else:
            def constraint(x):
                return x
            self.constraint = constraint
        self.norm = kwargs['norm']
        self.verbose = kwargs['verbose']
        self.max_damping = 10

    def __iter__(self):
        f = self.f
        x0 = self.constraint(self.x0)
        norm = self.norm
        J = self.J
        fx = self._matrix(f(x0))
        fxnorm = norm(fx)
        cancel = False
        x0 = self._matrix(x0)
        while not cancel:
            # get direction of descent
            fxn = -fx
            Jx = J(x0)
            try:
                s = self._Axb(Jx, fxn)
            except ZeroDivisionError:
                cancel = True
                break
            # damping step size TODO: better strategy (hard task)
            l = self._mpfloat('1.0')
            x1 = x0 + l*s
            damp_iter = 0
            while True:
                damp_iter += 1
                if x1.tolist() == x0.tolist() or damp_iter > self.max_damping:
                    if self.verbose > 1:
                        print("Solver: Found stationary point.")
                    cancel = True
                    break
                x1 = self._matrix(self.constraint(x1))
                fx = self._matrix(f(list(x1)))
                newnorm = norm(fx)
                if newnorm <= fxnorm:
                    # new x accepted
                    fxnorm = newnorm
                    x0 = x1
                    break
                l /= 2.0
                x1 = x0 + l*s
            yield (x0, fxnorm)

