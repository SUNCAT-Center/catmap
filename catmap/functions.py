import numpy as np
import catmap
import re
from copy import copy
from ase.atoms import string2symbols
import warnings

def get_composition(species_string):
    """
    Convert string of species into a dictionary of species and the number of each species.

    :param species_string: A string of the reaction species. Should be a chemical formula string
                           that may also contain '-','&',or,'pe'. 'pe' is a special case corresponding
                           to a proton-electron pair and has the compositon of H, while ele corresponds to an electron and has no associated atoms.
    :type species: str

    """
    composition = {}
    # clean up transition states and electrochem
    species_string = species_string.replace('-','')
    species_string = species_string.replace('pe','H')
    species_string = species_string.replace('&','')
    species_string = species_string.replace('ele','')
    try:
        symbs = string2symbols(species_string)
        for a in set(symbs):
            composition[a] = symbs.count(a)
    except ValueError:
        composition = None
    return composition

def cartesian_product(*args, **kwds):
    """
    Take the Cartesian product

    .. todo:: Explain what the args and kwds are
    """
    # product('ABCD', 'xy') --> Ax Ay Bx By Cx Cy Dx Dy
    # product(range(2), repeat=3) --> 000 001 010 011 100 101 110 111
    pools = map(tuple, args) * kwds.get('repeat', 1)
    result = [[]]
    for pool in pools:
                result = [x+[y] for x in result for y in pool]
    for prod in result:
                yield tuple(prod)

def convert_formation_energies(energy_dict,atomic_references,composition_dict):
    """
    Convert dictionary of energies, atomic references and compositions into a dictionary of formation energies

    :param energy_dict: Dictionary of energies for all species.
                        Keys should be species names and values
                        should be energies.
                        
    :type energy_dict: dict

    :param atomic_references: Dictionary of atomic reference energies (?)

    :type atomic_references: dict

    :param composition_dict: Dictionary of compositions

    :type composition_dict: dict

    .. todo:: Explain the keys and values for energy_dict, atomic_references, and composition_dict
    """
    n = len(atomic_references)
    R = np.zeros((n,n))
    e = []
    ref_offsets = {}
    atoms = sorted(atomic_references)
    for i,a in enumerate(atoms):
        composition = composition_dict[atomic_references[a]]
        e.append(energy_dict[atomic_references[a]])
        for j,a in enumerate(atoms):
            n_a = composition.get(a,0)
            R[i,j] = n_a
    if not np.prod([R[i,i] for i in range(0,n)]):
        raise ValueError('Reference set is not valid.')
    e = np.array(e)
    try:
        R_inv = np.linalg.solve(R,np.eye(n))
    except np.linalg.linalg.LinAlgError:
        raise ValueError('Reference set is not valid.')
    x = list(np.dot(R_inv,e))
    for a,v in zip(atoms,x):
        ref_offsets[a] = v
    new_data = {}
    for key in energy_dict:
        composition = composition_dict[key]
        E = energy_dict[key]
        for symb in composition:
            E -= ref_offsets[symb]*composition[symb]
        new_data[key] = round(E,5)
    return new_data,ref_offsets

def parse_constraint(minmaxlist,name):
    """
    Parse constraints for the relation. Returns two lists of minimum and maximum constraints

    :param minmaxlist: List of minimum and maximum constraints.

    :type minmaxlist: list

    :param name: Name for the list of constraints.

    :type name: str

    .. todo:: Explain minmaxlist and name
    """
    minlist = []
    maxlist = []
    for mm in minmaxlist:
        try:
            if mm is None:
                minv = -1e99
                maxv = 1e99
            elif str(mm).count(':') == 1:
                minv,maxv = [float(v) for v in mm.split(':')]
            elif mm == '+':
                minv = 0
                maxv = 1e99
            elif mm == '-':
                minv = -1e99
                maxv = 0
            else:
                minv = float(mm)
                maxv = float(mm)
        except (ValueError,TypeError,AttributeError):
            raise ValueError('Could not parse constraint for '+\
                    name+': '+str(minmaxlist))
        minlist.append(minv)
        maxlist.append(maxv)
    return minlist,maxlist

def constrained_relaxation(
        A,b,x0,x_min,x_max,max_iter = 100000,tolerance = 1e-10):
    """
    Solve Ax=b subject to the constraints that 
    x_i > x_min_i and x_i < x_max_i. Algorithm is from Axelson 1996.
    
    Note that x_min/Max are both lists/arrays of length equal to x

    :param A: A matrix.

    :type A: numpy.array

    :param b: b vector.

    :type b: numpy.array

    :param x0: x vector

    :type x0: numpy.array

    :param x_min: Minimum constraints.

    :type x_min: array_like

    :param x_max: Maximum constraints.

    :type x_max: array_like

    :param max_iter: Maximum number of iterations.

    :type max_iter: int, optional

    :param tolerance: Tolerance.

    :type tolerance: float, optional

    .. todo:: Check to make sure docstring is correct.
    """

    #define functional corresponding to Ax=b
    def J(x,A,b):
        """
        Functional of x which corresponds to Ax=b for 
        the relaxation method used.
        
        :param x: x vector.

        :type x: array_like

        :param A: A matrix.

        :type A: numpy.array

        :param b: b vector.

        :param b: array_like

        .. todo:: Check that docstring is correct
        """
        answer =  np.dot(
                np.dot(np.dot(x.T,A.T),A),x) - 2*np.dot(np.dot(b.T,A),x)
        return answer

    ai = x_min
    bi = x_max
    N = len(x0)

    def find_min(q):
        """
        Find minimum

        .. todo:: Explain what this does in the context of constrained_relaxation.
        """
        u[q] = 0
        v = np.dot(A,u)
        num1 = 0
        num2 = 0
        denom = 0
        Aq = A[:,q]
        for k in range(0,N):
            num2 += np.dot(v,Aq)
            num1 += np.dot(b,Aq)
            denom += np.dot(Aq,Aq)
        zeta = (num1-num2)/denom
        if zeta > bi[q]: zeta = bi[q]
        if zeta < ai[q]: zeta = ai[q]
        return zeta

    u = catmap.copy(x0)
    nIter =0
    converged = False
    while nIter < max_iter and converged == False:
        nIter += 1
        fOld = J(u,A,b)
        for j in range(0,N):
            u[j] = find_min(j)
        fNew = J(u,A,b)
        fDiff = fOld - fNew
        if np.linalg.norm(fDiff) < tolerance:
            converged = True

    if converged == True:
        return u
    else:
        raise ValueError('Constrained relaxation did not converge.'+
                'Residual was '+str(np.linalg.norm(fDiff)))

def scaling_coefficient_matrix(
        parameter_dict, descriptor_dict, surface_names, 
        parameter_names=None, coeff_mins = 0, coeff_maxs = 1e99, 
        return_error_dict = False):
    """Class for determining adsorption and transition-state energies
    as a linear function of descriptors. 

    :param parameter_dict: Dictionary where the key is adsorbate name 
                           and the value is a list of adsorption energies for each surface. 
                           If some surfaces do not have an adsorption energy use None
                           as a placeholder.

    :type parameter_dict: dict
    
    :param descriptor_dict: Dictionary where the key is surface name and the 
                            value is a list of descriptor values for each surface.
    
    :type descriptor_dict: dict

    :param surface_names: List of surfaces which defines the order of 
                          surface adsorption energies in parameter_dict.

    :type surface_names: list
    
    :param parameter_names: List of adsorbates which defines the order 
                            of adsorption coefficients in the output. Default is the 
                            order of parameter_dict.keys().

    :type parameter_names: list, optional

    :param coeff_mins: Defines the minimum value of the coefficient 
                       for each descriptor. Should be a matrix/array/list of lists
                       which matches the shape of the expected output.

    :type coeff_mins: float, optional

    :param coeff_maxs: Same as coeff_mins but for the maximum value of the coefficient.

    :type coeff_maxs: float, optional

    :param return_error_dict: Specify whether or not to return a dictionary of the errors.

    :type return_error_dict: bool, optional
    """

    #Define adsorbate order if it isn't
    if not parameter_names: parameter_names = parameter_dict.keys()

    #Check that input dictionaries are valid.
    def check_lengths(dictionary,force_numeric = False):
        """
        Check that the input dictionaries are valid.

        :param dictionary: Input dictionary.

        :type dictionary: dict

        :param force_numeric: Ensure that all values in the dictionary are numeric.

        :type force_numeric: bool, optional
        """
        for val in dictionary.values():
            if len(val) != len(dictionary.values()[0]):
                key_len = '\n'.join([key+':'+str(len(dictionary[key])) 
                    for key in dictionary])
                raise ValueError('All values must be lists of same length.'+
                         'Use None as placeholder. \nKey:Length\n'+key_len)
            if force_numeric:
                try:
                    [float(num) for num in val]
                except:
                    raise ValueError('All values must be numeric. '+
                            'Error when parsing ' + str(val))

    check_lengths(parameter_dict,False) 
    check_lengths(descriptor_dict,True) 

    #initialize coefficient matrix that will be returned.
    C = np.zeros(
            (len(descriptor_dict.values()[0])+1,len(parameter_dict.keys())))

    #initialize error dictionary that will be returned 
    #(if return_error_dict=True)
    error_dict = {}
    for key in parameter_dict:
        error_dict[key] = [None]*len(surface_names)

    #initialize descriptor matrix that will be used if 
    #all surfaces are present for a given adsorbate
    Dtotal = np.zeros((len(surface_names),len(descriptor_dict.values()[0])+1))
    for i,Dsurf in enumerate(surface_names):
        Dtotal[i,-1] = 1 #constant term.
        for j,DE in enumerate(descriptor_dict[Dsurf]):
            Dtotal[i,j] = float(DE)

    for Nads,ads in enumerate(parameter_names):
        #construct vectors for relaxation method.

        A = []

        #if mins and maxs are equal then the system is fully constrained and
        #there is no reason to solve for the parameters. However, in order to 
        #preserve ordering we use the known coeffs to put in "scaled" parameter
        #energy values and solve for the coefficients which, by definition, will
        #come out to be the same as the constraints.
        if coeff_mins[Nads] == coeff_maxs[Nads]:
            coeffs = coeff_mins[Nads]
            surfs = surface_names
            for surf in surfs:
                descriptors_i = descriptor_dict[surf]
                ads_i = sum([ci*di 
                    for ci,di in zip(coeffs,descriptors_i)]) + coeffs[-1]
                A.append(ads_i)

        #construct parameter vector A from parameters
        #if coefficients are not totally constrained
        if not A:
            A = [val for val in parameter_dict[ads] if val is not None]
            surfs = [surface_names[i] 
                    for i,val in enumerate(parameter_dict[ads]) 
                    if val is not None] #determine the surfaces which have 
                          #parameter energy values for this adsorbate

        try:
            A = np.array([float(val) for val in A])
        except ValueError:
            raise ValueError('Non-numeric value for the '
                    'parameters of '+ads+'.')

        #construct "descriptor" matrix (note that this is done inside the 
        #for loop to allow different parameters to have different number 
        #of surfaces)

        if len(surfs) <= len(descriptor_dict.values()[0])+1:
            warnings.warn('Number of energies specified is less than the number'
                          'of free parameters for '+ads+'. Scaling is not reliable'
                          'unless parameters are explicitly specified in '
                          ' constraints_dict.')

        if len(surfs) == len(surface_names):
            D = Dtotal
        else:
            D = np.zeros((len(surfs),len(descriptor_dict.values()[0])+1))
            for i,Dsurf in enumerate(surfs):
                D[i,-1] = 1 #constant term.
                for j,DE in enumerate(descriptor_dict[Dsurf]):
                    D[i,j] = float(DE)

        #find initial guess for the "coefficient" matrix by solving the 
        #unconstrained least-squares problem Dc=A using psuedo-inverse 
        #of D (note that this is not efficient, but it doesn't matter 
        #for such small matrices)
        D = np.array(D)
        A = np.array(A)
        if len(A) > 1:
            Dinv = np.linalg.pinv(D)
            c0 = np.dot(Dinv,A)

            #use relaxation method to solve the problem subject to the 
            #constraints specified by coeff_mins/Maxs.

            cMin = coeff_mins[Nads]
            cMax = coeff_maxs[Nads]

            c = constrained_relaxation(
                    D,A,c0,cMin,cMax)

        elif coeff_mins[Nads] == coeff_maxs[Nads]:
            c = coeff_mins[Nads]

        elif A:
            #If there is only one data point, assume constant.
            warnings.warn('Assuming constant value for: '+ads)
            c = [0]*len(Dtotal[i,:])
            c[-1] = A[0]
        else:
            warnings.warn('No data found for : '+ads+', assuming scaling parameters are 0.')
            c = [0]*len(Dtotal[i,:])

        for Ndesc,coeff in enumerate(c):
            C[Ndesc,Nads] = np.round(coeff,5)

        if return_error_dict == True:
            err = np.dot(D,c) - A
            for surf,errVal in zip(surfs,err):
                index = surface_names.index(surf)
                error_dict[ads][index] = np.round(errVal,6)
    if return_error_dict == True:
        return C, error_dict
    else:
        return C

def linear_regression(x,y,constrain_slope=None):
    """
    Perform linear regression on x and y and return the slope and intercept.

    :param x: x-coordinates.

    :type x: array_like

    :param y: y-coordinates.

    :type y: array_like

    :param constrain_slope: Slope constraint

    :type constrain_slope: float, optional
    """
    if constrain_slope is None:
        m,b = catmap.plt.polyfit(x,y,1)
    else:
        m = float(constrain_slope)
        b = sum([yi-m*xi for xi,yi in zip(x,y)])/len(x)
    return m,b

def match_regex(string,regex,group_names):
    """
    Find matching regular expression in string and return a dictionary of the matched expressions.

    :param string: String.

    :type string: str

    :param regex: Regular expression.

    :type regex: str

    :param group_names: Corresponding names for each matched group.

    :type group_names: list

    .. todo:: Check that this docstring is correct.
    """
    match_dict = {}
    match = re.match(regex,string)
    if match:
        for name,val in zip(group_names,match.groups()):
            match_dict[name] = val
        return match_dict
    else:
        return None
            
def numerical_jacobian(f, x, matrix, h = 1e-10,diff_idxs=None):
    """
    Calculate the Jacobian matrix of a function at the point x0.

    This is the first derivative of a vectorial function:

        f : R^m -> R^n with m >= n

    Hacked from mpmath.calculus.optimize

    :param f: Function.

    :type f: callable

    :param x:

    :type x:

    :param matrix:

    :type matrix: 

    :param h:

    :type h: float, optional

    .. todo:: Fill in the descriptions for f, x, matrix, and h
    """
    x = matrix(x)
    fx = matrix(f(x))
    m = len(fx)
    n = len(x)
    J = matrix(m, n)
    if not diff_idxs:
        diff_idxs = xrange(n)
    for j in diff_idxs:
        xj = x.copy()
        delta = abs(h*xj[j])
        delta = max(delta,h)
        #using delta proportional to xj is more stable
        #for very small numbers.
        xj[j] += delta
        Jj = (matrix(f(xj)) - fx)/(delta)
        for i in xrange(m):
            J[i,j] = Jj[i]
    return J

def smooth_piecewise_linear(theta_tot,slope=1,cutoff=0.25,smoothing=0.05):
    """
    Smooth piecewise linear function.

    :param theta_tot:

    :type theta_tot:

    :param slope: slope of line

    :type slope: float, optional

    :param cutoff: Cutoff.

    :type cutoff: float, optional

    :param smoothing: Amount of smoothing.

    :type smoothing: float, optional

    .. todo:: Fill in descriptions and types for theta_tot
    """
    x1 = cutoff + smoothing
    x0 = cutoff - smoothing
    if theta_tot <= x0:
        c_0 = 0
        dC = 0
        d2C = 0
    elif theta_tot <= x1:
        alpha = slope/(2*(x1-x0))
        c_0 = (alpha*(theta_tot-x0)**2)/theta_tot
        dC = alpha*(1-(x0/theta_tot)**2)
        d2C = (2*alpha*x0**2)/(theta_tot**3)
    else:
        c_0 = slope*(theta_tot - cutoff)/theta_tot
        dC = slope*(cutoff/(theta_tot**2))
        d2C = (-2*slope*cutoff)/(theta_tot**3)
    return c_0, dC, d2C

def offset_smooth_piecewise_linear(theta_tot,slope=1,cutoff=0.25, smoothing=0.05, offset=0.1):
    """
    Piecewise linear function with an offset. Not equivalent to piecewise linear
    for second-order interactions

    :param theta_tot:

    :type theta_tot:

    :param max_coverage: Maximum coverage.

    :type max_coverage: int, optional

    :param cutoff: Cutoff.

    :type cutoff: float, optional

    :param smoothing: Smoothing.

    :type smoothing: smoothing, optional

    :param offset: Offset.

    :type offset: float, optional

    .. todo:: Fill in description for theta_tot
    """
    c_0, dC, d2C = smooth_piecewise_linear(theta_tot,slope,cutoff,smoothing)
    c_0 += offset
    return c_0, dC, d2C

def add_dict_in_place(dict1, dict2):
    """
    Updates dict1 with elements in dict2 if they do not exist.  otherwise,
    add the value for key in dict2 to the value for that key in dict1

    :param dict1: Dictionary.

    :type dict1: dict

    :param dict2: Dictionary.

    :type dict2: dict
    """
    for k, v in dict2.iteritems():
        if k in dict1:
            dict1[k] += dict2[k]
        else:
            dict1[k] = dict2[k]


def fetch_all_output_variables():
    """Use code-inspection to extract all processed output variables from

            catmap.scalers.scaler_base.ScalerBase.set_output_attrs,
            catmap.solvers.solver_base.SolverBase.set_output_attrs

    New keywords should work out of the box if they are added
    in one of those functions and using one of the kind of if-statements
    that are already in place.

    """

    import ast
    import inspect

    import catmap
    import catmap.scalers
    import catmap.solvers

    keywords = []
    for infunc in [
        catmap.scalers.scaler_base.ScalerBase.set_output_attrs,
        catmap.solvers.solver_base.SolverBase.set_output_attrs
    ]:
        source = inspect.getsource(infunc)
        unsource = ''
        for line in source.split('\n'):
            unsource += line[4:] + '\n'
        node = ast.parse(unsource)
        f = node.body[0]

        class IfLister(ast.NodeVisitor):
            """ast.NodeVisitor subclass that extracts string  keywords
            from a function that parses keyword variables using
            code inspection.
            """
            def visit_If(self, node):
                """Override default behavior for visiting an if-statement.
                """
                if hasattr(node.test, 'left'):
                    if hasattr(node.test.left, 's'):
                        keywords.append(node.test.left.s)  # Done
                    elif hasattr(node.test, 'comparators'):
                        for comparator in node.test.comparators:
                            if hasattr(comparator, 's'):
                                keywords.append(comparator.s)
                            else:
                                pass  # dead-end
                    else:
                        pass  # dead-end

                elif hasattr(node.test, 'func'):
                    if hasattr(node.test.func, 'value'):
                        if hasattr(node.test.func.value, 'args'):
                            if hasattr(node.test.func.value.args[0], 'elts'):
                                for elt in node.test.func.value.args[0].elts:
                                    keywords.append(elt.s)
                elif hasattr(node.test, 'values'):
                    for elt in node.test.values:
                        if hasattr(elt, 'left') and hasattr(elt.left, 's'):
                            keywords.append(elt.left.s)
                else:
                    pass  # dead-end

        IfLister().visit(f)

    return keywords
