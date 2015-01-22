templates = {}

#LaTeX templates

templates['latex_longtable'] = r"""
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

templates['latex_summary'] = r"""
\documentclass[a4paper,8pt]{report}
\usepackage{booktabs}
\usepackage{longtable}
\usepackage[hmargin=1cm,vmargin=3cm]{geometry}
\begin{document}

${summary_txt}

\end{document}
"""

#Dynamically compiled function templates

templates['constrain_coverages'] = r"""
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

templates['rate_constants'] = r"""
def rate_constants(rxn_parameters,theta,gas_energies,site_energies,T,F,mpf,matrix,mpexp,include_derivatives=True):
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

    Gf,dGs =  interaction_function(theta,energies,interaction_vector,F,include_derivatives=include_derivatives)

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


templates['elementary_rates'] = r"""
def elementary_rates(rate_constants,theta,p,mpf,matrix):

    kf = rate_constants[0:len(rate_constants)/2]
    kr = rate_constants[len(rate_constants)/2:]

    r = matrix([0]*len(kf))
    dtheta_dt = matrix([0]*len(theta))
    
    ${steady_state_expressions}
    
    return r
"""

templates['interacting_mean_field_steady_state'] = r"""
def interacting_mean_field_steady_state(rxn_parameters,theta,p,gas_energies,site_energies,T,F,mpf,matrix,mpexp):

    ${rate_constants_no_derivatives}
    
    kf, kr, junk, junk= rate_constants(rxn_parameters,theta,gas_energies,site_energies,T,F,
            mpf,matrix,mpexp,include_derivatives=False)
   
    r = [0]*len(kf)
    dtheta_dt = [0]*len(theta)
    
    ${steady_state_expressions}
 
    r = matrix(r)
    dtheta_dt = matrix(dtheta_dt)

    return dtheta_dt
"""

templates['ideal_mean_field_steady_state'] = r"""
def ideal_mean_field_steady_state(kf,kr,theta,p,mpf,matrix):

    r = [0]*len(kf)
    dtheta_dt = [0]*len(theta)
    
    ${steady_state_expressions}

    r = matrix(r)
    dtheta_dt = matrix(dtheta_dt)
    
    return dtheta_dt
"""

templates['interacting_mean_field_jacobian'] = r"""
def interacting_mean_field_jacobian(rxn_parameters,theta,p,gas_energies,site_energies,T,F,mpf,matrix,mpexp):

#    print 'rxn_parameters = ', rxn_parameters
#    print 'theta = ', theta
#    print 'p = ', p
#    print 'gas_energies = ', gas_energies
#    print 'site_energies = ', site_energies
#    print 'T = ', T
    ${kB}
    ${n_adsorbates}
    kBT = kB*T


    J = [[0 for i in range(n_adsorbates)] for j in range(n_adsorbates)]

    ${rate_constants_with_derivatives}

    kf, kr, dEf, dEr = rate_constants(
                          rxn_parameters,theta,gas_energies,site_energies,T,F,
                          mpf,matrix,mpexp,include_derivatives=True)
    ${jacobian_expressions}
    
    J = matrix(J)
    return J
"""

templates['ideal_mean_field_jacobian'] = r"""
def ideal_mean_field_jacobian(kf,kr,theta,p,mpf,matrix):
    ${n_adsorbates}
    J = [[0 for i in range(n_adsorbates)] for j in range(n_adsorbates)]

    ${jacobian_expressions_no_derivatives}

    J = matrix(J)
    return J
"""

templates['first_order_interaction_function'] = r"""
def interaction_function(coverages,energies,interaction_vector,F,include_derivatives=True): 

#    Function for evaluating coverage-dependent intearction energies. 
#
#    -coverages is a vector of n_coverages coverages (thetas).
#    -energies is a vector of energies (length n_coverages).
#    -interaction_vector is a raveled matrix of self and cross interaction parameters
#    (length n_coverages^2).
#    -include_derivatives: True or False
#
#    Function evaluates:
#
#    E_i = E_0 + F(theta_tot)*sum_j(epsilon_ij*theta_j)
#
#    dE_i/dtheta_n = sum_j(dF/dtheta_tot*dtheta_tot/dtheta_n*epsilon_ij*theta_j + F*epsilon_ij*delta_jn)
#
#    where sum_j is summation over index j, epsilon_ij is the interaction parameter between 
#    adsorbate i and adsorbat j, theta_tot is a sum on theta, delta_nj is the kronecker delta,
#    theta_j is the coverage of adsorbate j, and F is the "response function".

    Es = []
    dEs = []
    n_cvg = len(coverages)

#    Dictionary with site names as keys and values a list of:
#    [indices_corresponding_to_site, max_coverage_of_site, piecewise_linearity_threshold]

    ${site_info_dict}

    f_vector = [0]*len(coverages)
    df_vector = [0]*len(coverages)
    sites = [0]*len(coverages)
    c_tots = [0]*len(coverages)
    for s in site_info_dict:
        idxs, max_cvg, F_params = site_info_dict[s]
        F_params['max_coverage'] = max_cvg
        cvgs = [coverages[j] for j in idxs]
        c_tot_i = sum(cvgs)
        f, df = F(c_tot_i,**F_params)
        for idx in idxs:
            f_vector[idx] = f
            df_vector[idx] = df
            sites[idx] = s
            c_tots[idx] = c_tot_i
    Es = []
    dEs = []

    for i,theta_i in enumerate(coverages):
        E_0 = energies[i]
        epsilon_vector = interaction_vector[i*n_cvg:(i+1)*n_cvg]
        f_epsilon_dot_theta = 0
        for ep,nc,f in zip(epsilon_vector,coverages,f_vector):
            f_epsilon_dot_theta += f*ep*nc

        E_i = E_0 + f_epsilon_dot_theta
        Es.append(E_i)

        if include_derivatives:
            diff_vector = []
            for n,theta_n in enumerate(coverages):
                df_sum = 0
                f_sum = 0
                for j,theta_j in enumerate(coverages):
                    df= df_vector[j]
                    f = f_vector[j]
                    ep = epsilon_vector[j]
                    if j==n:
                        delta = 1
                    else:
                        delta = 0
                    if sites[j] == sites[n]:
                        dsum = 1
                    else:
                        dsum = 0
                    df_sum += df*dsum*ep*theta_j
                    f_sum += f*ep*delta
                diff_vector.append(df_sum+f_sum)
            dEs.append(diff_vector)
        else:
            dEs = None
    return Es, dEs
    """

templates['second_order_interaction_function'] = r"""
def interaction_function(coverages,energies,epsilon,F,include_derivatives=True):
    ${site_info_dict}

    N_ads = len(coverages)
    N_sites = len(site_info_dict)
    idx_lists = []
    f = []
    df = []
    d2f = []
    for s in site_info_dict:
        idxs,max_cvg,F_params = site_info_dict[s]
        F_params['max_coverage'] = max_cvg
        idx_lists.append(site_info_dict[s][0])
        theta_tot = sum([coverages[i] for i in idxs])
        fs,dfs,d2fs = F(theta_tot,**F_params)
        f.append(fs)
        df.append(dfs)
        d2f.append(d2fs)

    term_1 = [0]*N_ads
    term_1_sum = [[0]*N_ads for i in range(N_ads)]
    term_2 = [0]*N_ads
    term_2_sum = [[0]*N_sites for i in range(N_sites)]

    for s in range(N_sites):
        for q in range(N_sites):
            for j in idx_lists[s]:
                for k in idx_lists[q]:
                    term_2_sum[q][s] += epsilon[j*N_ads+k]*coverages[j]*coverages[k]

    for q in range(N_sites):
        for n in range(N_ads):
            for j in idx_lists[q]:
                term_1_sum[n][q] += epsilon[j*N_ads+n]*coverages[j]

    #Double-check:
        #Should be possible to pull fk out of term_1
        #Make sure s,q index trick takes care of 1/2 term of sum2 properly
    for s in range(N_sites):
        for q in range(N_sites):
            for n in idx_lists[s]:
                term_2[n] += df[s]*f[q]*term_2_sum[q][s]
                term_1[n] += f[s]*f[q]*term_1_sum[n][q]
                    
    E_diff = [a+b+c for a,b,c in zip(energies,term_1,term_2)]

    if include_derivatives:
        E_jacob = [[0]*N_ads for i in range(N_ads)]
        for s in range(0,N_sites):
            for q in range(0,N_sites):
                for n in idx_lists[s]:
                    for m in idx_lists[s]:
                        prod_nmsq = df[s]*f[q]*term_1_sum[n][q]
                        E_jacob[n][m] += prod_nmsq 
                        E_jacob[m][n] += prod_nmsq
                        E_jacob[n][m] += d2f[s]*f[q]*term_2_sum[q][s]
                    for m in idx_lists[q]:
                        prod_nmsq2 = df[q]*f[s]*term_1_sum[n][q]
                        E_jacob[n][m] += prod_nmsq2
                        E_jacob[m][n] += prod_nmsq2
                        E_jacob[n][m] += df[s]*df[q]*term_2_sum[q][s]
                        E_jacob[n][m] += epsilon[m*N_ads+n]*f[s]*f[q]
    else:
        E_jacob = None

    return E_diff, E_jacob
"""

templates['ideal_interaction_function'] = r"""
def interaction_function(coverages,energies,interaction_vector,F,include_derivatives=True): 
    #Dummy function for non-interacting
    derivs = [[0]*len(coverages)]
    derivs = derivs*len(coverages)
    return energies, derivs
    """
