templates = {}

#LaTeX templates

templates['latex_longtable'] = """
\centering
\begin{longtable}{lllll}
\toprule
Species & Surface & Formation Energy [eV] & Vibrational Frequencies [meV] & Reference \\
\midrule
\endfirsthead
\toprule
Species & Surface & Formation Energy [eV] & Vibrational Frequencies [meV] & Reference \\
\midrule
\endhead
\hline
\endfoot
\hline
\endlastfoot

${longtable_txt}

\bottomrule
\end{longtable}
"""

templates['latex_summary'] = """
\documentclass[a4paper,8pt]{report}
\usepackage{booktabs}
\usepackage{longtable}
\usepackage[hmargin=1cm,vmargin=3cm]{geometry}
\begin{document}

${summary_txt}

\end{document}
"""

#Dynamically compiled function templates

templates['constrain_coverages'] = """
def constrain_coverage_function(cvgs,mpf,c_min):
    cvgs = [max(ci,c_min) for ci in cvgs]

    ${n_adsorbates}
    ${site_info_dict}
    ${max_coverage_list}

    #enforce explicit maxima
    cvgs = [min(ci,maxi) for ci,maxi in zip(cvgs,max_coverage_list)]

    #enforce site conservation
    for s in site_info_dict:
        idxs = [idx for idx in site_info_dict[s][0] if idx < n_adsorbates]
        tot = site_info_dict[s][1]
        sum_a = sum([cvgs[id] for id in idxs])
        if sum_a > tot:
            for id in idxs:
                cvgs[id] = cvgs[id]/sum_a
                cvgs[id] = cvgs[id]*tot
    return cvgs
"""

templates['rate_constants'] = """
def rate_constants(rxn_parameters,theta,gas_energies,site_energies,T,smearing,mpf,matrix,mpexp,exp,log,include_derivatives=True):
    ${interaction_function}

    kfs = []
    krs = []
    dEfs = []
    dErs = []

    ${n_adsorbates}
    ${n_transition_states}
    n_tot = n_adsorbates+n_transition_states
    if len(theta) == n_adsorbates:
        theta = list(theta) + [0]*n_transition_states #put transition-state coverages to 0
    elif len(theta) != n_adsorbates+n_transition_states:
        raise ValueError('Coverage vector was not correct length')
    energies = rxn_parameters[:n_tot]
    if len(rxn_parameters) == n_tot + n_tot**2:
        interaction_vector = rxn_parameters[-n_tot**2:]
    elif len(rxn_parameters) == n_tot:
        interaction_vector = [0]*n_tot**2
    else:
        raise ValueError('Length of reaction parameters is not correct. '+ str(rxn_parameters))

    Gf,dGs =  interaction_function(theta,energies,interaction_vector,smearing,exp,log,include_derivatives=include_derivatives)
#    if dGs is not None:
#        dGs = matrix(dGs)
    ${kB}
    ${h}
    ${prefactor_list}

    def element_wise_addition(lists):
        return [sum(L) for L in zip(*lists)]

    ${elementary_step_energetics}

    n_rxns = len(G_IS)
    for i in range(n_rxns):
        G_list = [G_IS[i],G_FS[i],G_TS[i]]
        G_TS_i = max(G_list)
        max_idx = G_list.index(G_TS_i)
        G_af[i] = G_TS_i - G_IS[i]
        G_ar[i] = G_TS_i - G_FS[i]
        kf = prefactor_list[i]*mpexp(-G_af[i]/(kB*T))
        kr = prefactor_list[i]*mpexp(-G_ar[i]/(kB*T))
        kfs.append(kf)
        krs.append(kr)
        if include_derivatives:
            dG_TS_i = [dG_IS[i],dG_FS[i],dG_TS[i]][max_idx]

            dG_TS_i = list(dG_TS_i[:n_adsorbates])
            dG_IS_i = list(dG_IS[i][:n_adsorbates])
            dG_FS_i = list(dG_FS[i][:n_adsorbates])

            dEf = [(-dG_TS_j + dG_IS_j) for dG_TS_j,dG_IS_j in zip(dG_TS_i,dG_IS_i)]
            dEr = [(-dG_TS_j + dG_FS_j) for dG_TS_j,dG_FS_j in zip(dG_TS_i,dG_FS_i)]
            dEfs.append(dEf)
            dErs.append(dEr)
        else:
            dEfs = dErs = None
    return kfs, krs, dEfs, dErs
"""


templates['elementary_rates'] = """
def elementary_rates(rate_constants,theta,p,mpf,matrix):

    kf = rate_constants[0:len(rate_constants)/2]
    kr = rate_constants[len(rate_constants)/2:]

    r = matrix([0]*len(kf))
    dtheta_dt = matrix([0]*len(theta))
    
    ${steady_state_expressions}
    
    return r

"""

templates['interacting_mean_field_steady_state'] = """
def interacting_mean_field_steady_state(rxn_parameters,theta,p,gas_energies,site_energies,T,smearing,mpf,matrix,mpexp,exp,log):

    ${rate_constants_no_derivatives}
    
    kf, kr, junk, junk= rate_constants(rxn_parameters,theta,gas_energies,site_energies,T,smearing,
            mpf,matrix,mpexp,exp,log,include_derivatives=False)
   
    r = [0]*len(kf)
    dtheta_dt = [0]*len(theta)
    
    ${steady_state_expressions}
 
    r = matrix(r)
    dtheta_dt = matrix(dtheta_dt)

    return dtheta_dt
"""

templates['ideal_mean_field_steady_state'] = """
def ideal_mean_field_steady_state(kf,kr,theta,p,mpf,matrix):

    r = [0]*len(kf)
    dtheta_dt = [0]*len(theta)
    
    ${steady_state_expressions}
    
    r = matrix(r)
    dtheta_dt = matrix(dtheta_dt)
    
    return dtheta_dt
"""

templates['interacting_mean_field_jacobian'] = """
def interacting_mean_field_jacobian(rxn_parameters,theta,p,gas_energies,site_energies,T,smearing,mpf,matrix,mpexp,exp,log):

    ${kB}
    ${n_adsorbates}
    kBT = kB*T


    J = [[0 for i in range(n_adsorbates)] for j in range(n_adsorbates)]

    ${rate_constants_with_derivatives}

    kf, kr, dEf, dEr = rate_constants(
                          rxn_parameters,theta,gas_energies,site_energies,T,smearing,
                          mpf,matrix,mpexp,exp,log,
                          include_derivatives=True)
    ${jacobian_expressions}
    
    J = matrix(J)
    return J
"""

templates['ideal_mean_field_jacobian'] = """
def ideal_mean_field_jacobian(kf,kr,theta,p,mpf,matrix):
    ${n_adsorbates}
    J = [[0 for i in range(n_adsorbates)] for j in range(n_adsorbates)]

    ${jacobian_expressions_no_derivatives}

    J = matrix(J)
    return J
"""

templates['linear_interaction_function'] = """
def interaction_function(coverages,energies,interaction_vector,smearing,exp,log,include_derivatives=True): 
#    Function for evaluating coverage-dependent intearction energies. 
#
#    -coverages is a vector of n_coverages coverages (thetas).
#    -energies is a vector of energies (length n_coverages).
#    -interaction_vector is a raveled matrix of self and cross interaction parameters
#    (length n_coverages^2).
#    -smearing is a constant for smoothing the piecewise linear function - corresponds
#    to smearing of a fermi distribution which is integrated to get piecewise linear.
#    -include_derivatives: True or False
#
#    Function evaluates:
#
#    E_i = E_0 + sum_j(epsilon_ij*theta_j*L(theta_tot)/theta_tot)
#
#    dE_i/dtheta_n = sum_j(epsilon_ij*((theta_tot*delta_nj - theta_j)/(theta_tot^2))*L(theta_tot)) +
#                    sum_j(epsilon_ij*theta_j*(dL/dtheta_tot)/theta_tot)
#
#    where sum_j is summation over index j, epsilon_ij is the interaction parameter between 
#    adsorbate i and adsorbat j, theta_tot is a sum on theta, delta_nj is the kronecker delta,
#    theta_j is the coverage of adsorbate j, and L is the smoothed approximation of piecewise
#    linearity obtained by integrateing the Fermi distribution.


    sigma_0 = 1./smearing
    Es = []
    dEs = []
    n_cvg = len(coverages)

    if interaction_vector is None or sum([abs(pi) for pi in interaction_vector]) == 0:
        return energies, [[0]*n_cvg]*n_cvg
#    Dictionary with site names as keys and values a list of:
#    [indices_corresponding_to_site, max_coverage_of_site, piecewise_linearity_threshold]
    
    ${site_info_dict}

#    def piecewise_linear(c_tot,sigma,cutoff,max_coverage):
#        expC = exp(sigma*c_tot)
#        expCutoff = exp(sigma*cutoff)
#        c_0 = (1./(sigma*max_coverage))*(log(1+expC/expCutoff))
#        if include_derivatives:
#            dC = (1./max_coverage)*(expC)/(expC + expCutoff)
#        else:
#            dC = None
#        return c_0, dC

    def piecewise_linear(c_tot,sigma,cutoff,max_coverage):
        sigma = 1./sigma
        x1 = cutoff + sigma
        x0 = cutoff - sigma
        slope = (1./max_coverage)
        if c_tot <= x0:
            c_0 = 0
            dC = 0
        elif c_tot <= x1:
            alpha = slope/(2*(x1-x0))
            c_0 = alpha*(c_tot-x0)**2
            dC = 2*alpha*(c_tot-x0)
        else:
            c_0 = slope*(c_tot - cutoff)
            dC = slope
        return c_0, dC

    L_vector = [0]*len(coverages)
    dL_vector = [0]*len(coverages)
    sites = [0]*len(coverages)
    normed_coverages = [0]*len(coverages)
    c_tots = [0]*len(coverages)

    for s in site_info_dict:
        idxs, max_cvg, cutoff = site_info_dict[s]
        cvgs = [coverages[j] for j in idxs]
        c_tot_i = sum(cvgs)
        if c_tot_i < 1e-10:
            #function is numerically unstable below c_tot < 1e-10
            if c_tot_i == 0:
                for idx in idxs:
                    coverages[idx] = 1e-10/float(len(idxs))
            else:
                cmax_idx = cvgs.index(max(cvgs))
                coverages[idxs[cmax_idx]] = 1e-11
            cvgs = [coverages[k] for k in idxs]
            c_tot_i = sum(cvgs)
        L, dL = piecewise_linear(c_tot_i,sigma_0,cutoff,max_cvg)
        for idx in idxs:
            L_vector[idx] = L
            dL_vector[idx] = dL
            sites[idx] = s
            c_tots[idx] = c_tot_i
            normed_coverages[idx] = coverages[idx]/c_tot_i

    for i,theta_i in enumerate(coverages):
        E_0 = energies[i]
        epsilon_vector = interaction_vector[i*n_cvg:(i+1)*n_cvg]
        L_epsilon_dot_theta = 0
        for ep,nc,L in zip(epsilon_vector,normed_coverages,L_vector):
            L_epsilon_dot_theta += L*ep*nc

        E_i = E_0 + L_epsilon_dot_theta
        Es.append(E_i) #HACK

        if include_derivatives:
            diff_vector = []
            for n, theta_n in enumerate(coverages):
                quotient_rule = []
                for j,theta_j in enumerate(coverages):
                    c_tot = c_tots[j]
                    if sites[n] == sites[j]:
                        if j == n:
                            quotient_rule.append((c_tot-theta_j)/(c_tot*c_tot))
                        else:
                            quotient_rule.append((-theta_j)/(c_tot*c_tot))
                    else:
                        if j == n:
                            quotient_rule.append((1.)/(c_tot))
                        else:
                            quotient_rule.append(0)
                            
                L_epsilon_dot_dQ = 0
                for ep,nc,L in zip(epsilon_vector,quotient_rule,L_vector):
                    L_epsilon_dot_dQ += ep*nc*L

                dL_epsilon_dot_theta = 0
                for j,ep,nc,dL in zip(range(len(epsilon_vector)),epsilon_vector,normed_coverages,dL_vector):
                    if sites[j] == sites[n]:
                        dL_epsilon_dot_theta += dL*ep*nc

                diff_vector.append(L_epsilon_dot_dQ+dL_epsilon_dot_theta)
            dEs.append(diff_vector)
    if include_derivatives:
        pass
#        dEs = [list(dE) for dE in zip(*dEs)]
    else:
        dEs = None
    return Es, dEs
    """
