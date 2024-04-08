import catmap
from catmap.model import ReactionModel
from catmap import ReactionModelWrapper
import numpy as np
import mpmath as mp
from catmap import string2symbols
import logging

class SolverBase(ReactionModelWrapper):
    def __init__(self,reaction_model=None):
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
        if reaction_model is None:
            reaction_model = ReactionModel()
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
            # If the numbers solver is being used
            # make sure to pass along the numbers and not
            # the coverage as an initial guess
            if self.use_numbers_solver:
                numbers = self._numbers
                self._coverage = list(self.solver.get_coverage(
                        rxn_parameters,c0=numbers)) #verify coverage
            else:
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

            # the following is useful for debugging/benchmarking
            # analytical derivatives, and should be commented out
            # for any production code.
            # import time
            # def J(x): #Use this to confirm the analytical jacobian is correct
            #     t0 = time.time()
            #     analytical = kwargs['J'](x)
            #     t_a = time.time() - t0
            #     t0 = time.time()
            #     numerical = catmap.functions.numerical_jacobian(f,x,matrix,1e-300)
            #     t_n = time.time() - t0
            #     error = analytical - numerical
            #     error = error.tolist()
            #     max_error = -1
            #     max_pos = None
            #     for i,ei in enumerate(error):
            #         for j,ej in enumerate(ei):
            #             if abs(ej) > 1e-10:
            #                 print('big error', ej, [i,j])
            #                 pass
            #             if abs(ej) > max_error:
            #                 max_error = abs(ej)
            #                 max_pos = [i,j]
            #     print('max_error', max_error, max_pos)
            #     print('t_analytic/t_numerical', t_a/t_n)
            #     return numerical
            # self.J = J

            # def J_numerical(x): #Use this to confirm the analytical jacobian is correct
            #     numerical = catmap.functions.numerical_jacobian(f,x,matrix,1e-50)
            #     return numerical
            # self.J = J_numerical


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
            except:
                try:
                    s = mp.qr_solve(Jx, fxn)[0]
                except ZeroDivisionError:
                    cancel = True
                    break
                except TypeError:
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

class ODESolver:
    def __init__(self, f, x0, matrix, mpfloat, **kwargs):
        self._matrix = matrix
        self._mpfloat = mpfloat
        self.f = f
        self.x0 = x0
        tol = kwargs.get('tol', self._mpfloat('1e-1'))
        self.dt = kwargs.get('dt', self._mpfloat('1e-15'))
        self.t0 = kwargs.get('t0', self._mpfloat('0.0'))

        self.dtheta_dt = lambda t, x: self.f(x)
        self.theta_t = mp.odefun(self.dtheta_dt, self.t0, self.x0, tol=tol, degree=2)
        self.t = self.t0

    def __iter__(self):
        cancel = False
        while not cancel:
            self.t += self.dt
            theta_t = self.theta_t(self.t)
            dtheta_dt = self.f(theta_t)
            yield (self.t, theta_t, dtheta_dt)

class NewtonRootNumbers:
    """
    Hacked from MDNewton in mpmath/calculus/optimization.py to
    implement the numbers solver method.

    Find the root of a vector function numerically using Newton's method.

    f is a vector function representing a nonlinear equation system.

    x0 is the starting point close to the root.

    J is a function returning the Jacobian matrix for a point.

    """

    def __init__(self, f, x0, math, matrix, mpfloat, Axb_solver, **kwargs):
        # Store an inner log file and get what to store
        self.verbose = kwargs['verbose']

        # Get the information about how the math should be handled
        self._math = math
        self._matrix = matrix
        self._mpfloat = mpfloat

        # Solver for the inner optimization problem
        self._Axb = Axb_solver

        # This f is the steady state function taking in coverages
        self.f = f

        # x0 is the starting coverage
        self.x0 = x0

        # Decide if x_star is fixed or free
        self.fix_x_star = kwargs['fix_x_star']

        if "rcond" in kwargs:
            self.rcond = kwargs['rcond']
        else:
            self.rcond = None

        # Decide if the Jacobian should be checked
        if kwargs['DEBUG']:
            # The Jacobian is checked based on the numerical Jacobian
            # determined by finite differences at the same point
            self.check_jacobian = True
            self.DEBUG = True
        else:
            # If not debugging, then don't check the Jacobian
            self.check_jacobian = False
            self.DEBUG = False

        # The Jacobian matrix which is a function of coverages
        if self.check_jacobian:
            self.J = self.J_confirm
        else:
            self.J = kwargs['J']

        # Store also the analytical Jacobian
        self.J_analytical = kwargs['J']

        # Convergence criteria
        self.norm = kwargs['norm']
        self.max_damping = kwargs['max_damping']

        # conversion from coverages to numbers and conversion
        # dtheta/dx matrix to get the Jacobian matrix
        self.conversion_function = kwargs['conversion_function']
        self.dtheta_dx_function = kwargs['dtheta_dx_function']

        # Store the precision
        self.precision = mp.power(self._mpfloat('10'), -1 * kwargs['precision'])

        # Decide on the number of empty sites
        self.esites = kwargs['esites']
        # The total coverage will be the number of empty sites * 1
        self.total_coverage = self.esites * self._mpfloat('1.0')


    def J_numerical(self, theta):
        """Pass the coverages with the f to get the numerical Jacobian."""
        numerical = catmap.functions.numerical_jacobian(self.f, theta, self._matrix, 1e-50)
        return numerical

    def J_confirm(self, theta):
        """The following is useful for debugging/benchmarking
        analytical derivatives, and should be commented out
        for any production code. It is a modified form
        of the J_confirm used in the solver above."""

        # Get the analytical Jacobian
        analytical = self.J_analytical(theta)
        # Get the numerical Jacobian
        numerical = self.J_numerical(theta)
        # Compare the analytical and numerical Jacobian
        error = analytical - numerical
        if self.fix_x_star:
            error = error[:-self.esites,:-self.esites]
        error = error.tolist()
        # Set bounds on the max and min errors
        max_error = -1
        max_pos = None
        for i,ei in enumerate(error):
            for j,ej in enumerate(ei):
                if abs(ej) > 1e-10:
                    print('big error in Jacobian:', ej, 'at position', [i,j])
                    pass
                if abs(ej) > max_error:
                    max_error = abs(ej)
                    max_pos = [i,j]
        if self.fix_x_star:
            return analytical[:-self.esites,:-self.esites], numerical[:-self.esites,:-self.esites]
        else:
            return analytical, numerical

    def moore_penrose_inverse(self, M, rcond=None):
        """Determine the moore-penrose inverse of a matrix."""
        # Take the SVD of M
        U, S, V = self._math.svd_r(M)
        if rcond:
            cutoff = rcond * mp.fabs(S[0])
        else:
            cutoff = self.precision

        # Take the inverse of S
        S_inv = []
        for i in range(len(S)):
            if mp.fabs(S[i]) < cutoff:
                S_inv.append(self._mpfloat('0.0'))
            else:
                S_inv.append(self._mpfloat('1.0') / S[i])

        # S is a diagonal matrix
        S_plus = self._math.diag(S_inv)

        # Take the pseudo-inverse of M
        M_plus = V.T * S_plus * U.T

        return M_plus

    def __iter__(self):
        # There are two possible routines for finding the root
        # through the multi-dimensional Newton method.
        # The first is to consider all x as independent variables
        # The second is to consider x_star = 1 and change the other
        if not self.fix_x_star:
            # Note that this is the objective function
            # That still depends on theta
            f = self.f
            # The Jacobian which depends on coverage
            J = self.J
            # matrix to get_dtheta_dx
            dtheta_dx_function = lambda x: self.dtheta_dx_function(list(x))
        else:
            # In case a fixed x_star is required, truncate the steady
            # state function, the jacobian in both theta and x spaces
            f = lambda theta: self.f(theta)[:-self.esites]
            if not self.check_jacobian:
                J = lambda theta: self.J(theta)[:-self.esites,:-self.esites]
            else:
                J = self.J # Checked in the function itself
            dtheta_dx_function = lambda x: self.dtheta_dx_function(x)[:-self.esites,:-self.esites]

        # Define the norm
        norm = self.norm

        # Initial guess in x
        x0 = self._matrix(self.x0)

        # conversion function that changes x within the iteration to theta
        conversion_function = lambda x: self.conversion_function(list(x))

        # get the initial coverages
        theta = conversion_function(x0)
        theta = self._matrix(theta)

        if self.DEBUG:
            # Check to make sure that the sum of coverages is not greater than 1
            assert mp.fabs(mp.fsum(theta) - self.total_coverage) < self.precision, "The sum of coverages is greater than 1.0"

        # Cancel if we have reached the right answer
        cancel = False

        while not cancel:
            # Begin the outer loop iteration
            fx = self._matrix(f(theta))

            if not self.fix_x_star:
                if self.DEBUG:
                    # Ensure that the sum of fx is 0
                    assert mp.fabs(mp.fsum(fx)) < self.precision, "fx is not zero"

                # Check if the objective function is the
                # same length as theta
                assert fx.rows == len(theta)

            if self.fix_x_star:
                assert fx.rows == len(theta) - self.esites

            # Get the norm of the objective function
            # this number is our error
            fxnorm = norm(fx)

            # get direction of descent
            fxn = -fx

            # Get the Jacobian matrix as a function of theta
            if self.check_jacobian:
                Jtheta, Jtheta_numerical = J(theta)
            else:
                Jtheta = J(theta)
                Jtheta_numerical = None

            # dtheta_dx matrix to convert Jacobian to x space
            dtheta_dx = dtheta_dx_function(x0)

            # Perform checks on the Jacobian(theta) matrix
            assert dtheta_dx.rows == Jtheta.rows, "Jacobian and dtheta_dx have different number of rows"
            assert dtheta_dx.cols == Jtheta.cols, "Jacobian and dtheta_dx have different number of columns"

            # Find Jx by taking the dot product of J(theta)
            # and dtheta/dx
            Jx =  Jtheta * dtheta_dx

            if self.check_jacobian:
                # Find the numerical equivalent of Jx
                Jx_numerical = Jtheta_numerical * dtheta_dx
            else:
                Jx_numerical = None

            if self.DEBUG:
                if not self.fix_x_star:
                    # The sum over all columns of the Jacobian
                    # matrix must be 0. This is because the term would be
                    # d/dtheta_i (sum_i f_i)
                    # where sum_i f_i = 0
                    for i in range(Jtheta.cols):
                        assert mp.fabs(mp.fsum(Jtheta[:, i])) < self.precision, "Jacobian column is not zero"
                    for i in range(dtheta_dx.cols):
                        assert mp.fabs(mp.fsum(dtheta_dx[:, i])) < self.precision, "dtheta_dx column is not zero"

                    # The Jacobian in x-space must also have the same
                    # checks as Jacobian in theta-space
                    assert Jx.rows == Jtheta.rows, "Jacobian and dtheta_dx have different number of rows"
                    assert Jx.cols == Jtheta.cols, "Jacobian and dtheta_dx have different number of columns"
                    for i in range(Jx.cols):
                        assert mp.fabs(mp.fsum(Jx[:, i])) < self.precision, "Jacobian column is not zero"


            if self.fix_x_star:
                try:
                    s = self._Axb(Jx, fxn)[0]
                except Exception as e:
                    try:
                        # Is the conventional solver fails used
                        # the Moore-Penrose pseudo-inverse to approx
                        # the inverse of the Jacobian
                        Jx_plus = self.moore_penrose_inverse(M=Jx)
                        s = Jx_plus * fxn
                    except Exception as e:
                        raise ValueError("Jx cannot be inverted because it is numerically singular")

                # Add 0 to the s matrix such that the empty site
                # x value for empty site is never updated
                s = list(s)
                for i in range(self.esites):
                    s.append(self._mpfloat('0.0'))
                s = self._matrix(s)

            else:
                # Since we are not fixing x_star, we are
                # solving an unconstrained problem, with more
                # variables than equations, so we need to use
                # a least-squares method, which starts with
                # taking the Moore-Penrose inverse of Jx
                # The conditioning number rcond is 10^5 times
                # the machine precision
                prec = -1*int(mp.log10(self.precision)/2)
                rcond = self._mpfloat(f'1e{prec}') * self.precision
                Jx_plus = self.moore_penrose_inverse(M=Jx, rcond=self.rcond)
                s = Jx_plus * fxn

            # damping step size
            l = self._mpfloat('1.0')
            x1 = x0 + l*s
            damp_iter = 0

            while True:
                if x1.tolist() == x0.tolist():
                    if self.verbose > 1:
                        print(f"x1 == x0, stationary point.")
                    cancel = True
                    break
                if damp_iter > self.max_damping:
                    if self.verbose > 1:
                        print(f"Damping step size exceeded max_damping.")
                    cancel = True
                    break

                # Convert to theta space for new iteration
                theta = conversion_function(x1)
                fx = self._matrix(f(theta))
                newnorm = norm(fx)

                if newnorm <= fxnorm:
                    # new x accepted
                    fxnorm = newnorm
                    x0 = x1
                    break
                l /= self._mpfloat('2.0')
                x1 = x0 + l*s
                damp_iter += 1

            if self.DEBUG:
                # Return the norm of the Jacobian as well
                Jxnorm = norm(Jx)
                Jx_numerical_norm = norm(Jx_numerical)
                yield (x0, fxnorm, Jxnorm, Jx_numerical_norm)
            else:
                yield (x0, fxnorm)
