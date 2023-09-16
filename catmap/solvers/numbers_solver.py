import abc
import csv
from typing import List, Dict, Any, Sequence

from .solver_base import NewtonRootNumbers


class BaseNumbersSolver(abc.ABC):
    def __init__(self, *args, **kwargs):
        pass

    @abc.abstractmethod
    def change_x_to_theta(self, x_including_surface, *args, **kwargs):
        pass

    @abc.abstractmethod
    def get_conversion_matrix(self, x_including_surface, *args, **kwargs):
        pass


class SquaredNumbersToTheta(BaseNumbersSolver):
    def __init__(self):
        pass

    @classmethod
    def change_x_to_theta(
        cls,
        x_including_surface,
        return_sum_sq_numbers=False,
        adsorbate_names: Sequence[str] = None,
        site_names: Sequence[str] = None,
        species_definitions: Dict[str, Any] = None,
        math_obj=None,
    ):
        # Each site has its own coverage which will
        # Split the adsorbate names into the different sites
        # The subtraction of 1 is to remove the 'g' site
        num_sites = len(site_names) - 1
        site_list = [[] for i in range(num_sites)]
        site_names = list(site_names)
        site_names.remove("g")
        species_list = list(adsorbate_names) + site_names
        species_index = []
        for i, ads in enumerate(species_list):
            site = species_definitions[ads]["site"]
            site_index = site_names.index(site)
            site_list[site_index].append(x_including_surface[i])
            species_index.append(site_index)
        sum_sq_numbers = [[] for i in range(num_sites)]
        for i, site_x in enumerate(site_list):
            sum_sq_numbers[i] = math_obj.fsum([math_obj.power(n, 2) for n in site_x])
        theta = []
        for i, x in enumerate(x_including_surface):
            coverage = math_obj.power(x, 2) / sum_sq_numbers[species_index[i]]
            theta.append(coverage)
        if return_sum_sq_numbers:
            return theta, sum_sq_numbers, species_index
        else:
            return theta

    @classmethod
    def get_conversion_matrix(
        cls,
        x_including_surface,
        adsorbate_names: Sequence[str] = None,
        site_names: Sequence[str] = None,
        species_definitions: Dict[str, Any] = None,
        math_obj=None,
        mpfloat=None,
    ):
        coverages, sum_sq_numbers, species_index = cls.change_x_to_theta(
            x_including_surface,
            adsorbate_names=adsorbate_names,
            site_names=site_names,
            species_definitions=species_definitions,
            math_obj=math_obj,
            return_sum_sq_numbers=True,
        )
        # Populate the diagonal of M with the dtheta/dx terms
        dtheta_dx_matrix = math_obj.matrix(len(coverages), len(coverages))
        for i in range(dtheta_dx_matrix.rows):
            for k in range(dtheta_dx_matrix.cols):
                if i == k:
                    # This is a diagonal element and so has the additional delta_{ik}
                    # term.
                    dtheta_dx_matrix[i, k] += mpfloat("2.0") * x_including_surface[i]
                    dtheta_dx_matrix[i, k] /= sum_sq_numbers[species_index[i]]
                # This extra term is to to account for the influence of
                # all other x's on a particular theta value, this
                # term is present only if the species index of both
                # the x's i and k are the same
                if species_index[i] == species_index[k]:
                    dtheta_dx_matrix[i, k] -= (
                        math_obj.power(x_including_surface[i], 2)
                        / math_obj.power(sum_sq_numbers[species_index[k]], 2)
                        * mpfloat("2.")
                        * x_including_surface[k]
                    )
        # The dtheta/dx matrix has all the terms for all the species + empty sites
        return dtheta_dx_matrix


def debug_writer(filename, writeout):
    with open(filename, "a") as csvfile:
        writer = csv.writer(
            csvfile, delimiter=",", quotechar="|", quoting=csv.QUOTE_MINIMAL
        )
        writer.writerow(writeout)


def get_theta_converter(numbers_type):
    if numbers_type == "squared" or numbers_type == None:
        return SquaredNumbersToTheta()
    else:
        raise ValueError(
            """\
Unknown type of numbers_type, use either squared or exponential."""
        )


class SteadyStateNumbersSolver:
    def __init__(self):
        pass

    def get_steady_state_numbers(
        self, rxn_parameters, steady_state_fn, jacobian_fn, c0=None
    ):
        if c0 is None:
            raise ValueError(
                "No initial numbers supplied. Mapper must supply initial guess"
            )
        totheta = get_theta_converter(self.numbers_type)
        self._rxn_parameters = rxn_parameters
        # Get the number of empty sites, remove 1 because we do
        # not want to include 'g' in the sites
        esites = len(self.site_names) - 1
        # The length of the initial numbers guess must include the
        # adsorbates as well as the slab value of exp(0).
        assert len(c0) == len(self.adsorbate_names) + esites
        # The steady state function is a function of the coverages
        # f(theta); Objective function which is a function of the coverage
        self.steady_state_function = steady_state_fn
        f = steady_state_fn
        solver = NewtonRootNumbers
        norm = self._math.lsqnorm
        change_x_to_theta = lambda x: totheta.change_x_to_theta(
            x,
            adsorbate_names=self.adsorbate_names,
            site_names=self.site_names,
            species_definitions=self.species_definitions,
            math_obj=self._math,
            return_sum_sq_numbers=False,
        )
        get_conversion_matrix = lambda x: totheta.get_conversion_matrix(
            x,
            adsorbate_names=self.adsorbate_names,
            site_names=self.site_names,
            species_definitions=self.species_definitions,
            math_obj=self._math,
            mpfloat=self._mpfloat,
        )
        # Before starting the calculation check if the job is a return job
        # That is, some part of this code has already used this function
        # and is returning the coverages instead of the numbers
        # The simple condition to check this is to see if the norm is lower
        # than the tolerance
        theta0 = change_x_to_theta(c0)
        if self.DEBUG:
            coverages_0 = change_x_to_theta(c0)
            writeout = self._descriptors + [0] + [float(_c) for _c in coverages_0]
            debug_writer("solution.csv", writeout)
            writeout = self._descriptors + [0] + [float(_x) for _x in c0]
            debug_writer("solution_numbers.csv", writeout)
            _norm_error = norm(f(theta0))
            writeout = self._descriptors + [0, _norm_error]
            debug_writer("error_log.csv", writeout)
        initial_norm_error = norm(f(theta0))
        if initial_norm_error <= self.tolerance:
            self._coverage = theta0
            self._numbers = c0
            return theta0
        solver_kwargs = dict(norm=norm, verbose=self.verbose)
        # The Jacobian function which is dependent on theta
        solver_kwargs["J"] = jacobian_fn
        # Conversion functions for converting x to theta
        solver_kwargs["conversion_function"] = change_x_to_theta
        solver_kwargs["dtheta_dx_function"] = get_conversion_matrix
        solver_kwargs["fix_x_star"] = self.fix_x_star
        solver_kwargs["esites"] = esites
        solver_kwargs["max_damping"] = self.max_damping_iterations
        solver_kwargs["DEBUG"] = self.DEBUG
        # Writes inner loop details to separate file
        solver_kwargs["verbose"] = 2
        # Store the decimal precision
        solver_kwargs["precision"] = self.decimal_precision
        # Run the solver; iterations is a generator
        iterations = solver(
            f,
            c0,
            self._math,
            self._matrix,
            self._mpfloat,
            self._math.qr_solve,
            **solver_kwargs
        )
        coverages = None
        maxiter = self.max_rootfinding_iterations
        iterations.maxiter = maxiter
        i = 1
        for data_iteration in iterations:
            if self.DEBUG:
                x, error, Jxnorm, Jxnorm_numer = data_iteration
                coverages = change_x_to_theta(list(x))
                debug_writer("error_log.csv", self._descriptors + [i, error])
                writeout = self._descriptors + [i, Jxnorm, Jxnorm_numer]
                debug_writer("jacobian_norm.csv", writeout)
                writeout = self._descriptors + [i] + [float(_c) for _c in coverages]
                debug_writer("solution.csv", writeout)
                writeout = self._descriptors + [i] + [float(_x) for _x in x]
                debug_writer("solution_numbers.csv", writeout)
            else:
                x, error = data_iteration
                coverages = change_x_to_theta(list(x))
            self.log("rootfinding_status", n_iter=i, resid=float(error), priority=1)
            i += 1
            if error < self.tolerance:
                self.log("rootfinding_success", n_iter=i, priority=1)
                self._error = error
                break
            if i >= maxiter:
                self.log("rootfinding_maxiter", n_iter=i, resid=float(error))
                raise ValueError("Out of iterations (resid=" + str(float(error)) + ")")
        if coverages:
            self._coverage = [c for c in coverages]
            self._numbers = [n for n in x]
            return [c for c in coverages]
        else:
            self.log("rootfinding_cancel", n_iter=i, resid=float(error))
            raise ValueError("Solver cancellation. (resid=" + str(float(error)) + ")")
