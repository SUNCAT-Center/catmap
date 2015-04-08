
def constrain_coverage_function(cvgs,mpf,c_min):
    cvgs = [max(ci,c_min) for ci in cvgs]

    n_adsorbates = 8
    site_info_dict = {'h': [[0], 0.33, {'cutoff': 0.62, 'max_coverage': 0.4, 'offset': 0.7, 'smoothing': 0.05}], 's': [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14], 1, {'cutoff': 0.62, 'max_coverage': 0.4, 'offset': 0.7, 'smoothing': 0.05}]}
    max_coverage_list = [0.33, 1, 1, 1, 1, 1, 1, 1]

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




def elementary_rates(rate_constants,theta,p,mpf,matrix):

    kf = rate_constants[0:len(rate_constants)/2]
    kr = rate_constants[len(rate_constants)/2:]

    r = matrix([0]*len(kf))
    dtheta_dt = matrix([0]*len(theta))
    
    s = [0]*3
    s[0] = (mpf('0.330000000000000015543122344752191565930843353271484375') - theta[0])
    s[1] = (mpf('1.0') - theta[1] - theta[2] - theta[3] - theta[4] - theta[5] - theta[6] - theta[7])
    s[2] = (mpf('0.0'))
    r[0] = kf[0]*p[3]*s[0]*s[0] - kr[0]*theta[0]*theta[0]
    r[1] = kf[1]*p[1]*s[1] - kr[1]*theta[6]
    r[2] = kf[2]*theta[6]*theta[0] - kr[2]*theta[4]*s[0]
    r[3] = kf[3]*theta[4]*theta[0] - kr[3]*theta[3]*s[0]
    r[4] = kf[4]*theta[3]*s[1] - kr[4]*theta[5]*theta[7]
    r[5] = kf[5]*theta[5]*theta[0] - kr[5]*theta[1]*s[0]
    r[6] = kf[6]*theta[1]*theta[0] - kr[6]*theta[2]*s[0]
    r[7] = kf[7]*theta[2]*theta[0] - kr[7]*p[0]*s[1]*s[0]
    r[8] = kf[8]*theta[7]*theta[0] - kr[8]*p[2]*s[1]*s[0]
    dtheta_dt[0] =  + 2*r[0] + -1*r[2] + -1*r[3] + -1*r[5] + -1*r[6] + -1*r[7] + -1*r[8]
    dtheta_dt[1] =  + 1*r[5] + -1*r[6]
    dtheta_dt[2] =  + 1*r[6] + -1*r[7]
    dtheta_dt[3] =  + 1*r[3] + -1*r[4]
    dtheta_dt[4] =  + 1*r[2] + -1*r[3]
    dtheta_dt[5] =  + 1*r[4] + -1*r[5]
    dtheta_dt[6] =  + 1*r[1] + -1*r[2]
    dtheta_dt[7] =  + 1*r[4] + -1*r[8]
    
    return r




def interacting_mean_field_jacobian(rxn_parameters,theta,p,gas_energies,site_energies,T,F,mpf,matrix,mpexp):

#    print 'rxn_parameters = ', rxn_parameters
#    print 'theta = ', theta
#    print 'p = ', p
#    print 'gas_energies = ', gas_energies
#    print 'site_energies = ', site_energies
#    print 'T = ', T
    kB = mpf('0.000086173324780000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000012')
    n_adsorbates = 8
    kBT = kB*T


    J = [[0 for i in range(n_adsorbates)] for j in range(n_adsorbates)]

    
    def rate_constants(rxn_parameters,theta,gas_energies,site_energies,T,F,mpf,matrix,mpexp,include_derivatives=True):
        
        def interaction_function(coverages,energies,epsilon,F,include_derivatives=True):
            site_info_dict = {'h': [[0], 0.33, {'cutoff': 0.62, 'max_coverage': 0.4, 'offset': 0.7, 'smoothing': 0.05}], 's': [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14], 1, {'cutoff': 0.62, 'max_coverage': 0.4, 'offset': 0.7, 'smoothing': 0.05}]}
        
            N_ads = len(coverages)
            N_sites = len(site_info_dict)
        
            idx_lists = []
            f = []
            df = []
            d2f = []
            for s in site_info_dict:
                idxs,max_cvg,F_params = site_info_dict[s]
                #Things might get strange when max_coverage != 1...
                if 'max_coverage' not in F_params:
                    F_params['max_coverage'] = max_cvg
                else:
                    F_params['max_coverage'] *= max_cvg
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
        
    
        kfs = []
        krs = []
        dEfs = []
        dErs = []
    
        n_adsorbates = 8
        n_transition_states = 7
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
    
        kB = mpf('0.000086173324780000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000012')
        h = mpf('0.0000000000000041356675160000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000024')
        prefactor_list = [kB*T/h, kB*T/h, kB*T/h, kB*T/h, kB*T/h, kB*T/h, kB*T/h, kB*T/h, kB*T/h]
    
        def element_wise_addition(lists):
            return [sum(L) for L in zip(*lists)]
    
        G_IS = [0]*9
        G_TS = [0]*9
        G_FS = [0]*9
        G_af = [0]*9
        G_ar = [0]*9
        dG_IS = [0]*9
        dG_TS = [0]*9
        dG_FS = [0]*9
        G_IS[0] = gas_energies[3] + site_energies[0] + site_energies[0]
        G_FS[0] = Gf[0] + Gf[0]
        dG_IS[0] = element_wise_addition([[0]*8])
        dG_FS[0] = element_wise_addition([dGs[0] , dGs[0]])
        G_TS[0] = max([G_IS[0],G_FS[0]])
        dG_TS[0] = None #determined later
        
        G_IS[1] = site_energies[1] + gas_energies[1]
        G_FS[1] = Gf[6]
        dG_IS[1] = element_wise_addition([[0]*8])
        dG_FS[1] = element_wise_addition([dGs[6]])
        G_TS[1] = max([G_IS[1],G_FS[1]])
        dG_TS[1] = None #determined later
        
        G_IS[2] = Gf[6] + Gf[0]
        G_FS[2] = Gf[4] + site_energies[0]
        dG_IS[2] = element_wise_addition([dGs[6] , dGs[0]])
        dG_FS[2] = element_wise_addition([dGs[4]])
        G_TS[2] = Gf[12] + site_energies[0]
        dG_TS[2] = element_wise_addition([dGs[12]])
        
        G_IS[3] = Gf[4] + Gf[0]
        G_FS[3] = Gf[3] + site_energies[0]
        dG_IS[3] = element_wise_addition([dGs[4] , dGs[0]])
        dG_FS[3] = element_wise_addition([dGs[3]])
        G_TS[3] = Gf[14] + site_energies[0]
        dG_TS[3] = element_wise_addition([dGs[14]])
        
        G_IS[4] = Gf[3] + site_energies[1]
        G_FS[4] = Gf[5] + Gf[7]
        dG_IS[4] = element_wise_addition([dGs[3]])
        dG_FS[4] = element_wise_addition([dGs[5] , dGs[7]])
        G_TS[4] = Gf[9] + site_energies[1]
        dG_TS[4] = element_wise_addition([dGs[9]])
        
        G_IS[5] = Gf[5] + Gf[0]
        G_FS[5] = Gf[1] + site_energies[0]
        dG_IS[5] = element_wise_addition([dGs[5] , dGs[0]])
        dG_FS[5] = element_wise_addition([dGs[1]])
        G_TS[5] = Gf[8] + site_energies[0]
        dG_TS[5] = element_wise_addition([dGs[8]])
        
        G_IS[6] = Gf[1] + Gf[0]
        G_FS[6] = Gf[2] + site_energies[0]
        dG_IS[6] = element_wise_addition([dGs[1] , dGs[0]])
        dG_FS[6] = element_wise_addition([dGs[2]])
        G_TS[6] = Gf[10] + site_energies[0]
        dG_TS[6] = element_wise_addition([dGs[10]])
        
        G_IS[7] = Gf[2] + Gf[0]
        G_FS[7] = gas_energies[0] + site_energies[1] + site_energies[0]
        dG_IS[7] = element_wise_addition([dGs[2] , dGs[0]])
        dG_FS[7] = element_wise_addition([[0]*8])
        G_TS[7] = Gf[11] + site_energies[0]
        dG_TS[7] = element_wise_addition([dGs[11]])
        
        G_IS[8] = Gf[7] + Gf[0]
        G_FS[8] = gas_energies[2] + site_energies[1] + site_energies[0]
        dG_IS[8] = element_wise_addition([dGs[7] , dGs[0]])
        dG_FS[8] = element_wise_addition([[0]*8])
        G_TS[8] = Gf[13] + site_energies[0]
        dG_TS[8] = element_wise_addition([dGs[13]])
        
    
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
    

    kf, kr, dEf, dEr = rate_constants(
                          rxn_parameters,theta,gas_energies,site_energies,T,F,
                          mpf,matrix,mpexp,include_derivatives=True)
    s = [0]*3
    s[0] = (mpf('0.330000000000000015543122344752191565930843353271484375') - theta[0])
    s[1] = (mpf('1.0') - theta[1] - theta[2] - theta[3] - theta[4] - theta[5] - theta[6] - theta[7])
    s[2] = (mpf('0.0'))
    kfkBT = [0]*9
    krkBT = [0]*9
    kfkBT[0] = kf[0]/kBT
    krkBT[0] = kr[0]/kBT
    kfkBT[1] = kf[1]/kBT
    krkBT[1] = kr[1]/kBT
    kfkBT[2] = kf[2]/kBT
    krkBT[2] = kr[2]/kBT
    kfkBT[3] = kf[3]/kBT
    krkBT[3] = kr[3]/kBT
    kfkBT[4] = kf[4]/kBT
    krkBT[4] = kr[4]/kBT
    kfkBT[5] = kf[5]/kBT
    krkBT[5] = kr[5]/kBT
    kfkBT[6] = kf[6]/kBT
    krkBT[6] = kr[6]/kBT
    kfkBT[7] = kf[7]/kBT
    krkBT[7] = kr[7]/kBT
    kfkBT[8] = kf[8]/kBT
    krkBT[8] = kr[8]/kBT
    J[0][0] = 0 + 2*(-2*kf[0]*p[3]*s[0] + (kfkBT[0])*dEf[0][0]*p[3]*s[0]*s[0] - 2*kr[0]*theta[0] - (krkBT[0])*dEr[0][0]*theta[0]*theta[0]) + -1*(kf[2]*theta[6] + (kfkBT[2])*dEf[2][0]*theta[6]*theta[0] - -1*kr[2]*theta[4] - (krkBT[2])*dEr[2][0]*theta[4]*s[0]) + -1*(kf[3]*theta[4] + (kfkBT[3])*dEf[3][0]*theta[4]*theta[0] - -1*kr[3]*theta[3] - (krkBT[3])*dEr[3][0]*theta[3]*s[0]) + -1*(kf[5]*theta[5] + (kfkBT[5])*dEf[5][0]*theta[5]*theta[0] - -1*kr[5]*theta[1] - (krkBT[5])*dEr[5][0]*theta[1]*s[0]) + -1*(kf[6]*theta[1] + (kfkBT[6])*dEf[6][0]*theta[1]*theta[0] - -1*kr[6]*theta[2] - (krkBT[6])*dEr[6][0]*theta[2]*s[0]) + -1*(kf[7]*theta[2] + (kfkBT[7])*dEf[7][0]*theta[2]*theta[0] - -1*kr[7]*p[0]*s[1] - (krkBT[7])*dEr[7][0]*p[0]*s[1]*s[0]) + -1*(kf[8]*theta[7] + (kfkBT[8])*dEf[8][0]*theta[7]*theta[0] - -1*kr[8]*p[2]*s[1] - (krkBT[8])*dEr[8][0]*p[2]*s[1]*s[0])
    J[0][1] = 0 + 2*(0 + (kfkBT[0])*dEf[0][1]*p[3]*s[0]*s[0] - 0 - (krkBT[0])*dEr[0][1]*theta[0]*theta[0]) + -1*(0 + (kfkBT[2])*dEf[2][1]*theta[6]*theta[0] - 0 - (krkBT[2])*dEr[2][1]*theta[4]*s[0]) + -1*(0 + (kfkBT[3])*dEf[3][1]*theta[4]*theta[0] - 0 - (krkBT[3])*dEr[3][1]*theta[3]*s[0]) + -1*(0 + (kfkBT[5])*dEf[5][1]*theta[5]*theta[0] - kr[5]*s[0] - (krkBT[5])*dEr[5][1]*theta[1]*s[0]) + -1*(kf[6]*theta[0] + (kfkBT[6])*dEf[6][1]*theta[1]*theta[0] - 0 - (krkBT[6])*dEr[6][1]*theta[2]*s[0]) + -1*(0 + (kfkBT[7])*dEf[7][1]*theta[2]*theta[0] - -1*kr[7]*p[0]*s[0] - (krkBT[7])*dEr[7][1]*p[0]*s[1]*s[0]) + -1*(0 + (kfkBT[8])*dEf[8][1]*theta[7]*theta[0] - -1*kr[8]*p[2]*s[0] - (krkBT[8])*dEr[8][1]*p[2]*s[1]*s[0])
    J[0][2] = 0 + 2*(0 + (kfkBT[0])*dEf[0][2]*p[3]*s[0]*s[0] - 0 - (krkBT[0])*dEr[0][2]*theta[0]*theta[0]) + -1*(0 + (kfkBT[2])*dEf[2][2]*theta[6]*theta[0] - 0 - (krkBT[2])*dEr[2][2]*theta[4]*s[0]) + -1*(0 + (kfkBT[3])*dEf[3][2]*theta[4]*theta[0] - 0 - (krkBT[3])*dEr[3][2]*theta[3]*s[0]) + -1*(0 + (kfkBT[5])*dEf[5][2]*theta[5]*theta[0] - 0 - (krkBT[5])*dEr[5][2]*theta[1]*s[0]) + -1*(0 + (kfkBT[6])*dEf[6][2]*theta[1]*theta[0] - kr[6]*s[0] - (krkBT[6])*dEr[6][2]*theta[2]*s[0]) + -1*(kf[7]*theta[0] + (kfkBT[7])*dEf[7][2]*theta[2]*theta[0] - -1*kr[7]*p[0]*s[0] - (krkBT[7])*dEr[7][2]*p[0]*s[1]*s[0]) + -1*(0 + (kfkBT[8])*dEf[8][2]*theta[7]*theta[0] - -1*kr[8]*p[2]*s[0] - (krkBT[8])*dEr[8][2]*p[2]*s[1]*s[0])
    J[0][3] = 0 + 2*(0 + (kfkBT[0])*dEf[0][3]*p[3]*s[0]*s[0] - 0 - (krkBT[0])*dEr[0][3]*theta[0]*theta[0]) + -1*(0 + (kfkBT[2])*dEf[2][3]*theta[6]*theta[0] - 0 - (krkBT[2])*dEr[2][3]*theta[4]*s[0]) + -1*(0 + (kfkBT[3])*dEf[3][3]*theta[4]*theta[0] - kr[3]*s[0] - (krkBT[3])*dEr[3][3]*theta[3]*s[0]) + -1*(0 + (kfkBT[5])*dEf[5][3]*theta[5]*theta[0] - 0 - (krkBT[5])*dEr[5][3]*theta[1]*s[0]) + -1*(0 + (kfkBT[6])*dEf[6][3]*theta[1]*theta[0] - 0 - (krkBT[6])*dEr[6][3]*theta[2]*s[0]) + -1*(0 + (kfkBT[7])*dEf[7][3]*theta[2]*theta[0] - -1*kr[7]*p[0]*s[0] - (krkBT[7])*dEr[7][3]*p[0]*s[1]*s[0]) + -1*(0 + (kfkBT[8])*dEf[8][3]*theta[7]*theta[0] - -1*kr[8]*p[2]*s[0] - (krkBT[8])*dEr[8][3]*p[2]*s[1]*s[0])
    J[0][4] = 0 + 2*(0 + (kfkBT[0])*dEf[0][4]*p[3]*s[0]*s[0] - 0 - (krkBT[0])*dEr[0][4]*theta[0]*theta[0]) + -1*(0 + (kfkBT[2])*dEf[2][4]*theta[6]*theta[0] - kr[2]*s[0] - (krkBT[2])*dEr[2][4]*theta[4]*s[0]) + -1*(kf[3]*theta[0] + (kfkBT[3])*dEf[3][4]*theta[4]*theta[0] - 0 - (krkBT[3])*dEr[3][4]*theta[3]*s[0]) + -1*(0 + (kfkBT[5])*dEf[5][4]*theta[5]*theta[0] - 0 - (krkBT[5])*dEr[5][4]*theta[1]*s[0]) + -1*(0 + (kfkBT[6])*dEf[6][4]*theta[1]*theta[0] - 0 - (krkBT[6])*dEr[6][4]*theta[2]*s[0]) + -1*(0 + (kfkBT[7])*dEf[7][4]*theta[2]*theta[0] - -1*kr[7]*p[0]*s[0] - (krkBT[7])*dEr[7][4]*p[0]*s[1]*s[0]) + -1*(0 + (kfkBT[8])*dEf[8][4]*theta[7]*theta[0] - -1*kr[8]*p[2]*s[0] - (krkBT[8])*dEr[8][4]*p[2]*s[1]*s[0])
    J[0][5] = 0 + 2*(0 + (kfkBT[0])*dEf[0][5]*p[3]*s[0]*s[0] - 0 - (krkBT[0])*dEr[0][5]*theta[0]*theta[0]) + -1*(0 + (kfkBT[2])*dEf[2][5]*theta[6]*theta[0] - 0 - (krkBT[2])*dEr[2][5]*theta[4]*s[0]) + -1*(0 + (kfkBT[3])*dEf[3][5]*theta[4]*theta[0] - 0 - (krkBT[3])*dEr[3][5]*theta[3]*s[0]) + -1*(kf[5]*theta[0] + (kfkBT[5])*dEf[5][5]*theta[5]*theta[0] - 0 - (krkBT[5])*dEr[5][5]*theta[1]*s[0]) + -1*(0 + (kfkBT[6])*dEf[6][5]*theta[1]*theta[0] - 0 - (krkBT[6])*dEr[6][5]*theta[2]*s[0]) + -1*(0 + (kfkBT[7])*dEf[7][5]*theta[2]*theta[0] - -1*kr[7]*p[0]*s[0] - (krkBT[7])*dEr[7][5]*p[0]*s[1]*s[0]) + -1*(0 + (kfkBT[8])*dEf[8][5]*theta[7]*theta[0] - -1*kr[8]*p[2]*s[0] - (krkBT[8])*dEr[8][5]*p[2]*s[1]*s[0])
    J[0][6] = 0 + 2*(0 + (kfkBT[0])*dEf[0][6]*p[3]*s[0]*s[0] - 0 - (krkBT[0])*dEr[0][6]*theta[0]*theta[0]) + -1*(kf[2]*theta[0] + (kfkBT[2])*dEf[2][6]*theta[6]*theta[0] - 0 - (krkBT[2])*dEr[2][6]*theta[4]*s[0]) + -1*(0 + (kfkBT[3])*dEf[3][6]*theta[4]*theta[0] - 0 - (krkBT[3])*dEr[3][6]*theta[3]*s[0]) + -1*(0 + (kfkBT[5])*dEf[5][6]*theta[5]*theta[0] - 0 - (krkBT[5])*dEr[5][6]*theta[1]*s[0]) + -1*(0 + (kfkBT[6])*dEf[6][6]*theta[1]*theta[0] - 0 - (krkBT[6])*dEr[6][6]*theta[2]*s[0]) + -1*(0 + (kfkBT[7])*dEf[7][6]*theta[2]*theta[0] - -1*kr[7]*p[0]*s[0] - (krkBT[7])*dEr[7][6]*p[0]*s[1]*s[0]) + -1*(0 + (kfkBT[8])*dEf[8][6]*theta[7]*theta[0] - -1*kr[8]*p[2]*s[0] - (krkBT[8])*dEr[8][6]*p[2]*s[1]*s[0])
    J[0][7] = 0 + 2*(0 + (kfkBT[0])*dEf[0][7]*p[3]*s[0]*s[0] - 0 - (krkBT[0])*dEr[0][7]*theta[0]*theta[0]) + -1*(0 + (kfkBT[2])*dEf[2][7]*theta[6]*theta[0] - 0 - (krkBT[2])*dEr[2][7]*theta[4]*s[0]) + -1*(0 + (kfkBT[3])*dEf[3][7]*theta[4]*theta[0] - 0 - (krkBT[3])*dEr[3][7]*theta[3]*s[0]) + -1*(0 + (kfkBT[5])*dEf[5][7]*theta[5]*theta[0] - 0 - (krkBT[5])*dEr[5][7]*theta[1]*s[0]) + -1*(0 + (kfkBT[6])*dEf[6][7]*theta[1]*theta[0] - 0 - (krkBT[6])*dEr[6][7]*theta[2]*s[0]) + -1*(0 + (kfkBT[7])*dEf[7][7]*theta[2]*theta[0] - -1*kr[7]*p[0]*s[0] - (krkBT[7])*dEr[7][7]*p[0]*s[1]*s[0]) + -1*(kf[8]*theta[0] + (kfkBT[8])*dEf[8][7]*theta[7]*theta[0] - -1*kr[8]*p[2]*s[0] - (krkBT[8])*dEr[8][7]*p[2]*s[1]*s[0])
    J[1][0] = 0 + 1*(kf[5]*theta[5] + (kfkBT[5])*dEf[5][0]*theta[5]*theta[0] - -1*kr[5]*theta[1] - (krkBT[5])*dEr[5][0]*theta[1]*s[0]) + -1*(kf[6]*theta[1] + (kfkBT[6])*dEf[6][0]*theta[1]*theta[0] - -1*kr[6]*theta[2] - (krkBT[6])*dEr[6][0]*theta[2]*s[0])
    J[1][1] = 0 + 1*(0 + (kfkBT[5])*dEf[5][1]*theta[5]*theta[0] - kr[5]*s[0] - (krkBT[5])*dEr[5][1]*theta[1]*s[0]) + -1*(kf[6]*theta[0] + (kfkBT[6])*dEf[6][1]*theta[1]*theta[0] - 0 - (krkBT[6])*dEr[6][1]*theta[2]*s[0])
    J[1][2] = 0 + 1*(0 + (kfkBT[5])*dEf[5][2]*theta[5]*theta[0] - 0 - (krkBT[5])*dEr[5][2]*theta[1]*s[0]) + -1*(0 + (kfkBT[6])*dEf[6][2]*theta[1]*theta[0] - kr[6]*s[0] - (krkBT[6])*dEr[6][2]*theta[2]*s[0])
    J[1][3] = 0 + 1*(0 + (kfkBT[5])*dEf[5][3]*theta[5]*theta[0] - 0 - (krkBT[5])*dEr[5][3]*theta[1]*s[0]) + -1*(0 + (kfkBT[6])*dEf[6][3]*theta[1]*theta[0] - 0 - (krkBT[6])*dEr[6][3]*theta[2]*s[0])
    J[1][4] = 0 + 1*(0 + (kfkBT[5])*dEf[5][4]*theta[5]*theta[0] - 0 - (krkBT[5])*dEr[5][4]*theta[1]*s[0]) + -1*(0 + (kfkBT[6])*dEf[6][4]*theta[1]*theta[0] - 0 - (krkBT[6])*dEr[6][4]*theta[2]*s[0])
    J[1][5] = 0 + 1*(kf[5]*theta[0] + (kfkBT[5])*dEf[5][5]*theta[5]*theta[0] - 0 - (krkBT[5])*dEr[5][5]*theta[1]*s[0]) + -1*(0 + (kfkBT[6])*dEf[6][5]*theta[1]*theta[0] - 0 - (krkBT[6])*dEr[6][5]*theta[2]*s[0])
    J[1][6] = 0 + 1*(0 + (kfkBT[5])*dEf[5][6]*theta[5]*theta[0] - 0 - (krkBT[5])*dEr[5][6]*theta[1]*s[0]) + -1*(0 + (kfkBT[6])*dEf[6][6]*theta[1]*theta[0] - 0 - (krkBT[6])*dEr[6][6]*theta[2]*s[0])
    J[1][7] = 0 + 1*(0 + (kfkBT[5])*dEf[5][7]*theta[5]*theta[0] - 0 - (krkBT[5])*dEr[5][7]*theta[1]*s[0]) + -1*(0 + (kfkBT[6])*dEf[6][7]*theta[1]*theta[0] - 0 - (krkBT[6])*dEr[6][7]*theta[2]*s[0])
    J[2][0] = 0 + 1*(kf[6]*theta[1] + (kfkBT[6])*dEf[6][0]*theta[1]*theta[0] - -1*kr[6]*theta[2] - (krkBT[6])*dEr[6][0]*theta[2]*s[0]) + -1*(kf[7]*theta[2] + (kfkBT[7])*dEf[7][0]*theta[2]*theta[0] - -1*kr[7]*p[0]*s[1] - (krkBT[7])*dEr[7][0]*p[0]*s[1]*s[0])
    J[2][1] = 0 + 1*(kf[6]*theta[0] + (kfkBT[6])*dEf[6][1]*theta[1]*theta[0] - 0 - (krkBT[6])*dEr[6][1]*theta[2]*s[0]) + -1*(0 + (kfkBT[7])*dEf[7][1]*theta[2]*theta[0] - -1*kr[7]*p[0]*s[0] - (krkBT[7])*dEr[7][1]*p[0]*s[1]*s[0])
    J[2][2] = 0 + 1*(0 + (kfkBT[6])*dEf[6][2]*theta[1]*theta[0] - kr[6]*s[0] - (krkBT[6])*dEr[6][2]*theta[2]*s[0]) + -1*(kf[7]*theta[0] + (kfkBT[7])*dEf[7][2]*theta[2]*theta[0] - -1*kr[7]*p[0]*s[0] - (krkBT[7])*dEr[7][2]*p[0]*s[1]*s[0])
    J[2][3] = 0 + 1*(0 + (kfkBT[6])*dEf[6][3]*theta[1]*theta[0] - 0 - (krkBT[6])*dEr[6][3]*theta[2]*s[0]) + -1*(0 + (kfkBT[7])*dEf[7][3]*theta[2]*theta[0] - -1*kr[7]*p[0]*s[0] - (krkBT[7])*dEr[7][3]*p[0]*s[1]*s[0])
    J[2][4] = 0 + 1*(0 + (kfkBT[6])*dEf[6][4]*theta[1]*theta[0] - 0 - (krkBT[6])*dEr[6][4]*theta[2]*s[0]) + -1*(0 + (kfkBT[7])*dEf[7][4]*theta[2]*theta[0] - -1*kr[7]*p[0]*s[0] - (krkBT[7])*dEr[7][4]*p[0]*s[1]*s[0])
    J[2][5] = 0 + 1*(0 + (kfkBT[6])*dEf[6][5]*theta[1]*theta[0] - 0 - (krkBT[6])*dEr[6][5]*theta[2]*s[0]) + -1*(0 + (kfkBT[7])*dEf[7][5]*theta[2]*theta[0] - -1*kr[7]*p[0]*s[0] - (krkBT[7])*dEr[7][5]*p[0]*s[1]*s[0])
    J[2][6] = 0 + 1*(0 + (kfkBT[6])*dEf[6][6]*theta[1]*theta[0] - 0 - (krkBT[6])*dEr[6][6]*theta[2]*s[0]) + -1*(0 + (kfkBT[7])*dEf[7][6]*theta[2]*theta[0] - -1*kr[7]*p[0]*s[0] - (krkBT[7])*dEr[7][6]*p[0]*s[1]*s[0])
    J[2][7] = 0 + 1*(0 + (kfkBT[6])*dEf[6][7]*theta[1]*theta[0] - 0 - (krkBT[6])*dEr[6][7]*theta[2]*s[0]) + -1*(0 + (kfkBT[7])*dEf[7][7]*theta[2]*theta[0] - -1*kr[7]*p[0]*s[0] - (krkBT[7])*dEr[7][7]*p[0]*s[1]*s[0])
    J[3][0] = 0 + 1*(kf[3]*theta[4] + (kfkBT[3])*dEf[3][0]*theta[4]*theta[0] - -1*kr[3]*theta[3] - (krkBT[3])*dEr[3][0]*theta[3]*s[0]) + -1*(0 + (kfkBT[4])*dEf[4][0]*theta[3]*s[1] - 0 - (krkBT[4])*dEr[4][0]*theta[5]*theta[7])
    J[3][1] = 0 + 1*(0 + (kfkBT[3])*dEf[3][1]*theta[4]*theta[0] - 0 - (krkBT[3])*dEr[3][1]*theta[3]*s[0]) + -1*(-1*kf[4]*theta[3] + (kfkBT[4])*dEf[4][1]*theta[3]*s[1] - 0 - (krkBT[4])*dEr[4][1]*theta[5]*theta[7])
    J[3][2] = 0 + 1*(0 + (kfkBT[3])*dEf[3][2]*theta[4]*theta[0] - 0 - (krkBT[3])*dEr[3][2]*theta[3]*s[0]) + -1*(-1*kf[4]*theta[3] + (kfkBT[4])*dEf[4][2]*theta[3]*s[1] - 0 - (krkBT[4])*dEr[4][2]*theta[5]*theta[7])
    J[3][3] = 0 + 1*(0 + (kfkBT[3])*dEf[3][3]*theta[4]*theta[0] - kr[3]*s[0] - (krkBT[3])*dEr[3][3]*theta[3]*s[0]) + -1*(kf[4]*(-1*theta[3] + 1*s[1]) + (kfkBT[4])*dEf[4][3]*theta[3]*s[1] - 0 - (krkBT[4])*dEr[4][3]*theta[5]*theta[7])
    J[3][4] = 0 + 1*(kf[3]*theta[0] + (kfkBT[3])*dEf[3][4]*theta[4]*theta[0] - 0 - (krkBT[3])*dEr[3][4]*theta[3]*s[0]) + -1*(-1*kf[4]*theta[3] + (kfkBT[4])*dEf[4][4]*theta[3]*s[1] - 0 - (krkBT[4])*dEr[4][4]*theta[5]*theta[7])
    J[3][5] = 0 + 1*(0 + (kfkBT[3])*dEf[3][5]*theta[4]*theta[0] - 0 - (krkBT[3])*dEr[3][5]*theta[3]*s[0]) + -1*(-1*kf[4]*theta[3] + (kfkBT[4])*dEf[4][5]*theta[3]*s[1] - kr[4]*theta[7] - (krkBT[4])*dEr[4][5]*theta[5]*theta[7])
    J[3][6] = 0 + 1*(0 + (kfkBT[3])*dEf[3][6]*theta[4]*theta[0] - 0 - (krkBT[3])*dEr[3][6]*theta[3]*s[0]) + -1*(-1*kf[4]*theta[3] + (kfkBT[4])*dEf[4][6]*theta[3]*s[1] - 0 - (krkBT[4])*dEr[4][6]*theta[5]*theta[7])
    J[3][7] = 0 + 1*(0 + (kfkBT[3])*dEf[3][7]*theta[4]*theta[0] - 0 - (krkBT[3])*dEr[3][7]*theta[3]*s[0]) + -1*(-1*kf[4]*theta[3] + (kfkBT[4])*dEf[4][7]*theta[3]*s[1] - kr[4]*theta[5] - (krkBT[4])*dEr[4][7]*theta[5]*theta[7])
    J[4][0] = 0 + 1*(kf[2]*theta[6] + (kfkBT[2])*dEf[2][0]*theta[6]*theta[0] - -1*kr[2]*theta[4] - (krkBT[2])*dEr[2][0]*theta[4]*s[0]) + -1*(kf[3]*theta[4] + (kfkBT[3])*dEf[3][0]*theta[4]*theta[0] - -1*kr[3]*theta[3] - (krkBT[3])*dEr[3][0]*theta[3]*s[0])
    J[4][1] = 0 + 1*(0 + (kfkBT[2])*dEf[2][1]*theta[6]*theta[0] - 0 - (krkBT[2])*dEr[2][1]*theta[4]*s[0]) + -1*(0 + (kfkBT[3])*dEf[3][1]*theta[4]*theta[0] - 0 - (krkBT[3])*dEr[3][1]*theta[3]*s[0])
    J[4][2] = 0 + 1*(0 + (kfkBT[2])*dEf[2][2]*theta[6]*theta[0] - 0 - (krkBT[2])*dEr[2][2]*theta[4]*s[0]) + -1*(0 + (kfkBT[3])*dEf[3][2]*theta[4]*theta[0] - 0 - (krkBT[3])*dEr[3][2]*theta[3]*s[0])
    J[4][3] = 0 + 1*(0 + (kfkBT[2])*dEf[2][3]*theta[6]*theta[0] - 0 - (krkBT[2])*dEr[2][3]*theta[4]*s[0]) + -1*(0 + (kfkBT[3])*dEf[3][3]*theta[4]*theta[0] - kr[3]*s[0] - (krkBT[3])*dEr[3][3]*theta[3]*s[0])
    J[4][4] = 0 + 1*(0 + (kfkBT[2])*dEf[2][4]*theta[6]*theta[0] - kr[2]*s[0] - (krkBT[2])*dEr[2][4]*theta[4]*s[0]) + -1*(kf[3]*theta[0] + (kfkBT[3])*dEf[3][4]*theta[4]*theta[0] - 0 - (krkBT[3])*dEr[3][4]*theta[3]*s[0])
    J[4][5] = 0 + 1*(0 + (kfkBT[2])*dEf[2][5]*theta[6]*theta[0] - 0 - (krkBT[2])*dEr[2][5]*theta[4]*s[0]) + -1*(0 + (kfkBT[3])*dEf[3][5]*theta[4]*theta[0] - 0 - (krkBT[3])*dEr[3][5]*theta[3]*s[0])
    J[4][6] = 0 + 1*(kf[2]*theta[0] + (kfkBT[2])*dEf[2][6]*theta[6]*theta[0] - 0 - (krkBT[2])*dEr[2][6]*theta[4]*s[0]) + -1*(0 + (kfkBT[3])*dEf[3][6]*theta[4]*theta[0] - 0 - (krkBT[3])*dEr[3][6]*theta[3]*s[0])
    J[4][7] = 0 + 1*(0 + (kfkBT[2])*dEf[2][7]*theta[6]*theta[0] - 0 - (krkBT[2])*dEr[2][7]*theta[4]*s[0]) + -1*(0 + (kfkBT[3])*dEf[3][7]*theta[4]*theta[0] - 0 - (krkBT[3])*dEr[3][7]*theta[3]*s[0])
    J[5][0] = 0 + 1*(0 + (kfkBT[4])*dEf[4][0]*theta[3]*s[1] - 0 - (krkBT[4])*dEr[4][0]*theta[5]*theta[7]) + -1*(kf[5]*theta[5] + (kfkBT[5])*dEf[5][0]*theta[5]*theta[0] - -1*kr[5]*theta[1] - (krkBT[5])*dEr[5][0]*theta[1]*s[0])
    J[5][1] = 0 + 1*(-1*kf[4]*theta[3] + (kfkBT[4])*dEf[4][1]*theta[3]*s[1] - 0 - (krkBT[4])*dEr[4][1]*theta[5]*theta[7]) + -1*(0 + (kfkBT[5])*dEf[5][1]*theta[5]*theta[0] - kr[5]*s[0] - (krkBT[5])*dEr[5][1]*theta[1]*s[0])
    J[5][2] = 0 + 1*(-1*kf[4]*theta[3] + (kfkBT[4])*dEf[4][2]*theta[3]*s[1] - 0 - (krkBT[4])*dEr[4][2]*theta[5]*theta[7]) + -1*(0 + (kfkBT[5])*dEf[5][2]*theta[5]*theta[0] - 0 - (krkBT[5])*dEr[5][2]*theta[1]*s[0])
    J[5][3] = 0 + 1*(kf[4]*(-1*theta[3] + 1*s[1]) + (kfkBT[4])*dEf[4][3]*theta[3]*s[1] - 0 - (krkBT[4])*dEr[4][3]*theta[5]*theta[7]) + -1*(0 + (kfkBT[5])*dEf[5][3]*theta[5]*theta[0] - 0 - (krkBT[5])*dEr[5][3]*theta[1]*s[0])
    J[5][4] = 0 + 1*(-1*kf[4]*theta[3] + (kfkBT[4])*dEf[4][4]*theta[3]*s[1] - 0 - (krkBT[4])*dEr[4][4]*theta[5]*theta[7]) + -1*(0 + (kfkBT[5])*dEf[5][4]*theta[5]*theta[0] - 0 - (krkBT[5])*dEr[5][4]*theta[1]*s[0])
    J[5][5] = 0 + 1*(-1*kf[4]*theta[3] + (kfkBT[4])*dEf[4][5]*theta[3]*s[1] - kr[4]*theta[7] - (krkBT[4])*dEr[4][5]*theta[5]*theta[7]) + -1*(kf[5]*theta[0] + (kfkBT[5])*dEf[5][5]*theta[5]*theta[0] - 0 - (krkBT[5])*dEr[5][5]*theta[1]*s[0])
    J[5][6] = 0 + 1*(-1*kf[4]*theta[3] + (kfkBT[4])*dEf[4][6]*theta[3]*s[1] - 0 - (krkBT[4])*dEr[4][6]*theta[5]*theta[7]) + -1*(0 + (kfkBT[5])*dEf[5][6]*theta[5]*theta[0] - 0 - (krkBT[5])*dEr[5][6]*theta[1]*s[0])
    J[5][7] = 0 + 1*(-1*kf[4]*theta[3] + (kfkBT[4])*dEf[4][7]*theta[3]*s[1] - kr[4]*theta[5] - (krkBT[4])*dEr[4][7]*theta[5]*theta[7]) + -1*(0 + (kfkBT[5])*dEf[5][7]*theta[5]*theta[0] - 0 - (krkBT[5])*dEr[5][7]*theta[1]*s[0])
    J[6][0] = 0 + 1*(0 + (kfkBT[1])*dEf[1][0]*p[1]*s[1] - 0 - (krkBT[1])*dEr[1][0]*theta[6]) + -1*(kf[2]*theta[6] + (kfkBT[2])*dEf[2][0]*theta[6]*theta[0] - -1*kr[2]*theta[4] - (krkBT[2])*dEr[2][0]*theta[4]*s[0])
    J[6][1] = 0 + 1*(-1*kf[1]*p[1] + (kfkBT[1])*dEf[1][1]*p[1]*s[1] - 0 - (krkBT[1])*dEr[1][1]*theta[6]) + -1*(0 + (kfkBT[2])*dEf[2][1]*theta[6]*theta[0] - 0 - (krkBT[2])*dEr[2][1]*theta[4]*s[0])
    J[6][2] = 0 + 1*(-1*kf[1]*p[1] + (kfkBT[1])*dEf[1][2]*p[1]*s[1] - 0 - (krkBT[1])*dEr[1][2]*theta[6]) + -1*(0 + (kfkBT[2])*dEf[2][2]*theta[6]*theta[0] - 0 - (krkBT[2])*dEr[2][2]*theta[4]*s[0])
    J[6][3] = 0 + 1*(-1*kf[1]*p[1] + (kfkBT[1])*dEf[1][3]*p[1]*s[1] - 0 - (krkBT[1])*dEr[1][3]*theta[6]) + -1*(0 + (kfkBT[2])*dEf[2][3]*theta[6]*theta[0] - 0 - (krkBT[2])*dEr[2][3]*theta[4]*s[0])
    J[6][4] = 0 + 1*(-1*kf[1]*p[1] + (kfkBT[1])*dEf[1][4]*p[1]*s[1] - 0 - (krkBT[1])*dEr[1][4]*theta[6]) + -1*(0 + (kfkBT[2])*dEf[2][4]*theta[6]*theta[0] - kr[2]*s[0] - (krkBT[2])*dEr[2][4]*theta[4]*s[0])
    J[6][5] = 0 + 1*(-1*kf[1]*p[1] + (kfkBT[1])*dEf[1][5]*p[1]*s[1] - 0 - (krkBT[1])*dEr[1][5]*theta[6]) + -1*(0 + (kfkBT[2])*dEf[2][5]*theta[6]*theta[0] - 0 - (krkBT[2])*dEr[2][5]*theta[4]*s[0])
    J[6][6] = 0 + 1*(-1*kf[1]*p[1] + (kfkBT[1])*dEf[1][6]*p[1]*s[1] - kr[1] - (krkBT[1])*dEr[1][6]*theta[6]) + -1*(kf[2]*theta[0] + (kfkBT[2])*dEf[2][6]*theta[6]*theta[0] - 0 - (krkBT[2])*dEr[2][6]*theta[4]*s[0])
    J[6][7] = 0 + 1*(-1*kf[1]*p[1] + (kfkBT[1])*dEf[1][7]*p[1]*s[1] - 0 - (krkBT[1])*dEr[1][7]*theta[6]) + -1*(0 + (kfkBT[2])*dEf[2][7]*theta[6]*theta[0] - 0 - (krkBT[2])*dEr[2][7]*theta[4]*s[0])
    J[7][0] = 0 + 1*(0 + (kfkBT[4])*dEf[4][0]*theta[3]*s[1] - 0 - (krkBT[4])*dEr[4][0]*theta[5]*theta[7]) + -1*(kf[8]*theta[7] + (kfkBT[8])*dEf[8][0]*theta[7]*theta[0] - -1*kr[8]*p[2]*s[1] - (krkBT[8])*dEr[8][0]*p[2]*s[1]*s[0])
    J[7][1] = 0 + 1*(-1*kf[4]*theta[3] + (kfkBT[4])*dEf[4][1]*theta[3]*s[1] - 0 - (krkBT[4])*dEr[4][1]*theta[5]*theta[7]) + -1*(0 + (kfkBT[8])*dEf[8][1]*theta[7]*theta[0] - -1*kr[8]*p[2]*s[0] - (krkBT[8])*dEr[8][1]*p[2]*s[1]*s[0])
    J[7][2] = 0 + 1*(-1*kf[4]*theta[3] + (kfkBT[4])*dEf[4][2]*theta[3]*s[1] - 0 - (krkBT[4])*dEr[4][2]*theta[5]*theta[7]) + -1*(0 + (kfkBT[8])*dEf[8][2]*theta[7]*theta[0] - -1*kr[8]*p[2]*s[0] - (krkBT[8])*dEr[8][2]*p[2]*s[1]*s[0])
    J[7][3] = 0 + 1*(kf[4]*(-1*theta[3] + 1*s[1]) + (kfkBT[4])*dEf[4][3]*theta[3]*s[1] - 0 - (krkBT[4])*dEr[4][3]*theta[5]*theta[7]) + -1*(0 + (kfkBT[8])*dEf[8][3]*theta[7]*theta[0] - -1*kr[8]*p[2]*s[0] - (krkBT[8])*dEr[8][3]*p[2]*s[1]*s[0])
    J[7][4] = 0 + 1*(-1*kf[4]*theta[3] + (kfkBT[4])*dEf[4][4]*theta[3]*s[1] - 0 - (krkBT[4])*dEr[4][4]*theta[5]*theta[7]) + -1*(0 + (kfkBT[8])*dEf[8][4]*theta[7]*theta[0] - -1*kr[8]*p[2]*s[0] - (krkBT[8])*dEr[8][4]*p[2]*s[1]*s[0])
    J[7][5] = 0 + 1*(-1*kf[4]*theta[3] + (kfkBT[4])*dEf[4][5]*theta[3]*s[1] - kr[4]*theta[7] - (krkBT[4])*dEr[4][5]*theta[5]*theta[7]) + -1*(0 + (kfkBT[8])*dEf[8][5]*theta[7]*theta[0] - -1*kr[8]*p[2]*s[0] - (krkBT[8])*dEr[8][5]*p[2]*s[1]*s[0])
    J[7][6] = 0 + 1*(-1*kf[4]*theta[3] + (kfkBT[4])*dEf[4][6]*theta[3]*s[1] - 0 - (krkBT[4])*dEr[4][6]*theta[5]*theta[7]) + -1*(0 + (kfkBT[8])*dEf[8][6]*theta[7]*theta[0] - -1*kr[8]*p[2]*s[0] - (krkBT[8])*dEr[8][6]*p[2]*s[1]*s[0])
    J[7][7] = 0 + 1*(-1*kf[4]*theta[3] + (kfkBT[4])*dEf[4][7]*theta[3]*s[1] - kr[4]*theta[5] - (krkBT[4])*dEr[4][7]*theta[5]*theta[7]) + -1*(kf[8]*theta[0] + (kfkBT[8])*dEf[8][7]*theta[7]*theta[0] - -1*kr[8]*p[2]*s[0] - (krkBT[8])*dEr[8][7]*p[2]*s[1]*s[0])
    
    J = matrix(J)
    return J




def ideal_mean_field_steady_state(kf,kr,theta,p,mpf,matrix):

    r = [0]*len(kf)
    dtheta_dt = [0]*len(theta)
    
    s = [0]*3
    s[0] = (mpf('0.330000000000000015543122344752191565930843353271484375') - theta[0])
    s[1] = (mpf('1.0') - theta[1] - theta[2] - theta[3] - theta[4] - theta[5] - theta[6] - theta[7])
    s[2] = (mpf('0.0'))
    r[0] = kf[0]*p[3]*s[0]*s[0] - kr[0]*theta[0]*theta[0]
    r[1] = kf[1]*p[1]*s[1] - kr[1]*theta[6]
    r[2] = kf[2]*theta[6]*theta[0] - kr[2]*theta[4]*s[0]
    r[3] = kf[3]*theta[4]*theta[0] - kr[3]*theta[3]*s[0]
    r[4] = kf[4]*theta[3]*s[1] - kr[4]*theta[5]*theta[7]
    r[5] = kf[5]*theta[5]*theta[0] - kr[5]*theta[1]*s[0]
    r[6] = kf[6]*theta[1]*theta[0] - kr[6]*theta[2]*s[0]
    r[7] = kf[7]*theta[2]*theta[0] - kr[7]*p[0]*s[1]*s[0]
    r[8] = kf[8]*theta[7]*theta[0] - kr[8]*p[2]*s[1]*s[0]
    dtheta_dt[0] =  + 2*r[0] + -1*r[2] + -1*r[3] + -1*r[5] + -1*r[6] + -1*r[7] + -1*r[8]
    dtheta_dt[1] =  + 1*r[5] + -1*r[6]
    dtheta_dt[2] =  + 1*r[6] + -1*r[7]
    dtheta_dt[3] =  + 1*r[3] + -1*r[4]
    dtheta_dt[4] =  + 1*r[2] + -1*r[3]
    dtheta_dt[5] =  + 1*r[4] + -1*r[5]
    dtheta_dt[6] =  + 1*r[1] + -1*r[2]
    dtheta_dt[7] =  + 1*r[4] + -1*r[8]

    r = matrix(r)
    dtheta_dt = matrix(dtheta_dt)
    
    return dtheta_dt




def ideal_mean_field_jacobian(kf,kr,theta,p,mpf,matrix):
    n_adsorbates = 8
    J = [[0 for i in range(n_adsorbates)] for j in range(n_adsorbates)]

    s = [0]*3
    s[0] = (mpf('0.330000000000000015543122344752191565930843353271484375') - theta[0])
    s[1] = (mpf('1.0') - theta[1] - theta[2] - theta[3] - theta[4] - theta[5] - theta[6] - theta[7])
    s[2] = (mpf('0.0'))
    J[0][0] = 0 + 2*(-2*kf[0]*p[3]*s[0] - 2*kr[0]*theta[0]) + -1*(kf[2]*theta[6] - -1*kr[2]*theta[4]) + -1*(kf[3]*theta[4] - -1*kr[3]*theta[3]) + -1*(kf[5]*theta[5] - -1*kr[5]*theta[1]) + -1*(kf[6]*theta[1] - -1*kr[6]*theta[2]) + -1*(kf[7]*theta[2] - -1*kr[7]*p[0]*s[1]) + -1*(kf[8]*theta[7] - -1*kr[8]*p[2]*s[1])
    J[0][1] = 0 + -1*-1*kr[5]*s[0] + -1*kf[6]*theta[0] + -1*-1*-1*kr[7]*p[0]*s[0] + -1*-1*-1*kr[8]*p[2]*s[0]
    J[0][2] = 0 + -1*-1*kr[6]*s[0] + -1*(kf[7]*theta[0] - -1*kr[7]*p[0]*s[0]) + -1*-1*-1*kr[8]*p[2]*s[0]
    J[0][3] = 0 + -1*-1*kr[3]*s[0] + -1*-1*-1*kr[7]*p[0]*s[0] + -1*-1*-1*kr[8]*p[2]*s[0]
    J[0][4] = 0 + -1*-1*kr[2]*s[0] + -1*kf[3]*theta[0] + -1*-1*-1*kr[7]*p[0]*s[0] + -1*-1*-1*kr[8]*p[2]*s[0]
    J[0][5] = 0 + -1*kf[5]*theta[0] + -1*-1*-1*kr[7]*p[0]*s[0] + -1*-1*-1*kr[8]*p[2]*s[0]
    J[0][6] = 0 + -1*kf[2]*theta[0] + -1*-1*-1*kr[7]*p[0]*s[0] + -1*-1*-1*kr[8]*p[2]*s[0]
    J[0][7] = 0 + -1*-1*-1*kr[7]*p[0]*s[0] + -1*(kf[8]*theta[0] - -1*kr[8]*p[2]*s[0])
    J[1][0] = 0 + 1*(kf[5]*theta[5] - -1*kr[5]*theta[1]) + -1*(kf[6]*theta[1] - -1*kr[6]*theta[2])
    J[1][1] = 0 + 1*-1*kr[5]*s[0] + -1*kf[6]*theta[0]
    J[1][2] = 0 + -1*-1*kr[6]*s[0]
    J[1][3] = 0
    J[1][4] = 0
    J[1][5] = 0 + 1*kf[5]*theta[0]
    J[1][6] = 0
    J[1][7] = 0
    J[2][0] = 0 + 1*(kf[6]*theta[1] - -1*kr[6]*theta[2]) + -1*(kf[7]*theta[2] - -1*kr[7]*p[0]*s[1])
    J[2][1] = 0 + 1*kf[6]*theta[0] + -1*-1*-1*kr[7]*p[0]*s[0]
    J[2][2] = 0 + 1*-1*kr[6]*s[0] + -1*(kf[7]*theta[0] - -1*kr[7]*p[0]*s[0])
    J[2][3] = 0 + -1*-1*-1*kr[7]*p[0]*s[0]
    J[2][4] = 0 + -1*-1*-1*kr[7]*p[0]*s[0]
    J[2][5] = 0 + -1*-1*-1*kr[7]*p[0]*s[0]
    J[2][6] = 0 + -1*-1*-1*kr[7]*p[0]*s[0]
    J[2][7] = 0 + -1*-1*-1*kr[7]*p[0]*s[0]
    J[3][0] = 0 + 1*(kf[3]*theta[4] - -1*kr[3]*theta[3])
    J[3][1] = 0 + -1*-1*kf[4]*theta[3]
    J[3][2] = 0 + -1*-1*kf[4]*theta[3]
    J[3][3] = 0 + 1*-1*kr[3]*s[0] + -1*kf[4]*(-1*theta[3] + 1*s[1])
    J[3][4] = 0 + 1*kf[3]*theta[0] + -1*-1*kf[4]*theta[3]
    J[3][5] = 0 + -1*(-1*kf[4]*theta[3] - kr[4]*theta[7])
    J[3][6] = 0 + -1*-1*kf[4]*theta[3]
    J[3][7] = 0 + -1*(-1*kf[4]*theta[3] - kr[4]*theta[5])
    J[4][0] = 0 + 1*(kf[2]*theta[6] - -1*kr[2]*theta[4]) + -1*(kf[3]*theta[4] - -1*kr[3]*theta[3])
    J[4][1] = 0
    J[4][2] = 0
    J[4][3] = 0 + -1*-1*kr[3]*s[0]
    J[4][4] = 0 + 1*-1*kr[2]*s[0] + -1*kf[3]*theta[0]
    J[4][5] = 0
    J[4][6] = 0 + 1*kf[2]*theta[0]
    J[4][7] = 0
    J[5][0] = 0 + -1*(kf[5]*theta[5] - -1*kr[5]*theta[1])
    J[5][1] = 0 + 1*-1*kf[4]*theta[3] + -1*-1*kr[5]*s[0]
    J[5][2] = 0 + 1*-1*kf[4]*theta[3]
    J[5][3] = 0 + 1*kf[4]*(-1*theta[3] + 1*s[1])
    J[5][4] = 0 + 1*-1*kf[4]*theta[3]
    J[5][5] = 0 + 1*(-1*kf[4]*theta[3] - kr[4]*theta[7]) + -1*kf[5]*theta[0]
    J[5][6] = 0 + 1*-1*kf[4]*theta[3]
    J[5][7] = 0 + 1*(-1*kf[4]*theta[3] - kr[4]*theta[5])
    J[6][0] = 0 + -1*(kf[2]*theta[6] - -1*kr[2]*theta[4])
    J[6][1] = 0 + 1*-1*kf[1]*p[1]
    J[6][2] = 0 + 1*-1*kf[1]*p[1]
    J[6][3] = 0 + 1*-1*kf[1]*p[1]
    J[6][4] = 0 + 1*-1*kf[1]*p[1] + -1*-1*kr[2]*s[0]
    J[6][5] = 0 + 1*-1*kf[1]*p[1]
    J[6][6] = 0 + 1*(-1*kf[1]*p[1] - kr[1]) + -1*kf[2]*theta[0]
    J[6][7] = 0 + 1*-1*kf[1]*p[1]
    J[7][0] = 0 + -1*(kf[8]*theta[7] - -1*kr[8]*p[2]*s[1])
    J[7][1] = 0 + 1*-1*kf[4]*theta[3] + -1*-1*-1*kr[8]*p[2]*s[0]
    J[7][2] = 0 + 1*-1*kf[4]*theta[3] + -1*-1*-1*kr[8]*p[2]*s[0]
    J[7][3] = 0 + 1*kf[4]*(-1*theta[3] + 1*s[1]) + -1*-1*-1*kr[8]*p[2]*s[0]
    J[7][4] = 0 + 1*-1*kf[4]*theta[3] + -1*-1*-1*kr[8]*p[2]*s[0]
    J[7][5] = 0 + 1*(-1*kf[4]*theta[3] - kr[4]*theta[7]) + -1*-1*-1*kr[8]*p[2]*s[0]
    J[7][6] = 0 + 1*-1*kf[4]*theta[3] + -1*-1*-1*kr[8]*p[2]*s[0]
    J[7][7] = 0 + 1*(-1*kf[4]*theta[3] - kr[4]*theta[5]) + -1*(kf[8]*theta[0] - -1*kr[8]*p[2]*s[0])

    J = matrix(J)
    return J




def interaction_function(coverages,energies,epsilon,F,include_derivatives=True):
    site_info_dict = {'h': [[0], 0.33, {'cutoff': 0.62, 'max_coverage': 0.4, 'offset': 0.7, 'smoothing': 0.05}], 's': [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14], 1, {'cutoff': 0.62, 'max_coverage': 0.4, 'offset': 0.7, 'smoothing': 0.05}]}

    N_ads = len(coverages)
    N_sites = len(site_info_dict)

    idx_lists = []
    f = []
    df = []
    d2f = []
    for s in site_info_dict:
        idxs,max_cvg,F_params = site_info_dict[s]
        #Things might get strange when max_coverage != 1...
        if 'max_coverage' not in F_params:
            F_params['max_coverage'] = max_cvg
        else:
            F_params['max_coverage'] *= max_cvg
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




def rate_constants(rxn_parameters,theta,gas_energies,site_energies,T,F,mpf,matrix,mpexp,include_derivatives=True):
    
    def interaction_function(coverages,energies,epsilon,F,include_derivatives=True):
        site_info_dict = {'h': [[0], 0.33, {'cutoff': 0.62, 'max_coverage': 0.4, 'offset': 0.7, 'smoothing': 0.05}], 's': [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14], 1, {'cutoff': 0.62, 'max_coverage': 0.4, 'offset': 0.7, 'smoothing': 0.05}]}
    
        N_ads = len(coverages)
        N_sites = len(site_info_dict)
    
        idx_lists = []
        f = []
        df = []
        d2f = []
        for s in site_info_dict:
            idxs,max_cvg,F_params = site_info_dict[s]
            #Things might get strange when max_coverage != 1...
            if 'max_coverage' not in F_params:
                F_params['max_coverage'] = max_cvg
            else:
                F_params['max_coverage'] *= max_cvg
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
    

    kfs = []
    krs = []
    dEfs = []
    dErs = []

    n_adsorbates = 8
    n_transition_states = 7
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

    kB = mpf('0.000086173324780000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000012')
    h = mpf('0.0000000000000041356675160000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000024')
    prefactor_list = [kB*T/h, kB*T/h, kB*T/h, kB*T/h, kB*T/h, kB*T/h, kB*T/h, kB*T/h, kB*T/h]

    def element_wise_addition(lists):
        return [sum(L) for L in zip(*lists)]

    G_IS = [0]*9
    G_TS = [0]*9
    G_FS = [0]*9
    G_af = [0]*9
    G_ar = [0]*9
    dG_IS = [0]*9
    dG_TS = [0]*9
    dG_FS = [0]*9
    G_IS[0] = gas_energies[3] + site_energies[0] + site_energies[0]
    G_FS[0] = Gf[0] + Gf[0]
    dG_IS[0] = element_wise_addition([[0]*8])
    dG_FS[0] = element_wise_addition([dGs[0] , dGs[0]])
    G_TS[0] = max([G_IS[0],G_FS[0]])
    dG_TS[0] = None #determined later
    
    G_IS[1] = site_energies[1] + gas_energies[1]
    G_FS[1] = Gf[6]
    dG_IS[1] = element_wise_addition([[0]*8])
    dG_FS[1] = element_wise_addition([dGs[6]])
    G_TS[1] = max([G_IS[1],G_FS[1]])
    dG_TS[1] = None #determined later
    
    G_IS[2] = Gf[6] + Gf[0]
    G_FS[2] = Gf[4] + site_energies[0]
    dG_IS[2] = element_wise_addition([dGs[6] , dGs[0]])
    dG_FS[2] = element_wise_addition([dGs[4]])
    G_TS[2] = Gf[12] + site_energies[0]
    dG_TS[2] = element_wise_addition([dGs[12]])
    
    G_IS[3] = Gf[4] + Gf[0]
    G_FS[3] = Gf[3] + site_energies[0]
    dG_IS[3] = element_wise_addition([dGs[4] , dGs[0]])
    dG_FS[3] = element_wise_addition([dGs[3]])
    G_TS[3] = Gf[14] + site_energies[0]
    dG_TS[3] = element_wise_addition([dGs[14]])
    
    G_IS[4] = Gf[3] + site_energies[1]
    G_FS[4] = Gf[5] + Gf[7]
    dG_IS[4] = element_wise_addition([dGs[3]])
    dG_FS[4] = element_wise_addition([dGs[5] , dGs[7]])
    G_TS[4] = Gf[9] + site_energies[1]
    dG_TS[4] = element_wise_addition([dGs[9]])
    
    G_IS[5] = Gf[5] + Gf[0]
    G_FS[5] = Gf[1] + site_energies[0]
    dG_IS[5] = element_wise_addition([dGs[5] , dGs[0]])
    dG_FS[5] = element_wise_addition([dGs[1]])
    G_TS[5] = Gf[8] + site_energies[0]
    dG_TS[5] = element_wise_addition([dGs[8]])
    
    G_IS[6] = Gf[1] + Gf[0]
    G_FS[6] = Gf[2] + site_energies[0]
    dG_IS[6] = element_wise_addition([dGs[1] , dGs[0]])
    dG_FS[6] = element_wise_addition([dGs[2]])
    G_TS[6] = Gf[10] + site_energies[0]
    dG_TS[6] = element_wise_addition([dGs[10]])
    
    G_IS[7] = Gf[2] + Gf[0]
    G_FS[7] = gas_energies[0] + site_energies[1] + site_energies[0]
    dG_IS[7] = element_wise_addition([dGs[2] , dGs[0]])
    dG_FS[7] = element_wise_addition([[0]*8])
    G_TS[7] = Gf[11] + site_energies[0]
    dG_TS[7] = element_wise_addition([dGs[11]])
    
    G_IS[8] = Gf[7] + Gf[0]
    G_FS[8] = gas_energies[2] + site_energies[1] + site_energies[0]
    dG_IS[8] = element_wise_addition([dGs[7] , dGs[0]])
    dG_FS[8] = element_wise_addition([[0]*8])
    G_TS[8] = Gf[13] + site_energies[0]
    dG_TS[8] = element_wise_addition([dGs[13]])
    

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




def interacting_mean_field_steady_state(rxn_parameters,theta,p,gas_energies,site_energies,T,F,mpf,matrix,mpexp):

    
    def rate_constants(rxn_parameters,theta,gas_energies,site_energies,T,F,mpf,matrix,mpexp,include_derivatives=True):
        
        def interaction_function(coverages,energies,epsilon,F,include_derivatives=True):
            site_info_dict = {'h': [[0], 0.33, {'cutoff': 0.62, 'max_coverage': 0.4, 'offset': 0.7, 'smoothing': 0.05}], 's': [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14], 1, {'cutoff': 0.62, 'max_coverage': 0.4, 'offset': 0.7, 'smoothing': 0.05}]}
        
            N_ads = len(coverages)
            N_sites = len(site_info_dict)
        
            idx_lists = []
            f = []
            df = []
            d2f = []
            for s in site_info_dict:
                idxs,max_cvg,F_params = site_info_dict[s]
                #Things might get strange when max_coverage != 1...
                if 'max_coverage' not in F_params:
                    F_params['max_coverage'] = max_cvg
                else:
                    F_params['max_coverage'] *= max_cvg
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
        
    
        kfs = []
        krs = []
        dEfs = []
        dErs = []
    
        n_adsorbates = 8
        n_transition_states = 7
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
    
        kB = mpf('0.000086173324780000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000012')
        h = mpf('0.0000000000000041356675160000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000024')
        prefactor_list = [kB*T/h, kB*T/h, kB*T/h, kB*T/h, kB*T/h, kB*T/h, kB*T/h, kB*T/h, kB*T/h]
    
        def element_wise_addition(lists):
            return [sum(L) for L in zip(*lists)]
    
        G_IS = [0]*9
        G_TS = [0]*9
        G_FS = [0]*9
        G_af = [0]*9
        G_ar = [0]*9
        G_IS[0] = gas_energies[3] + site_energies[0] + site_energies[0]
        G_FS[0] = Gf[0] + Gf[0]
        G_TS[0] = max([G_IS[0],G_FS[0]])
        
        G_IS[1] = site_energies[1] + gas_energies[1]
        G_FS[1] = Gf[6]
        G_TS[1] = max([G_IS[1],G_FS[1]])
        
        G_IS[2] = Gf[6] + Gf[0]
        G_FS[2] = Gf[4] + site_energies[0]
        G_TS[2] = Gf[12] + site_energies[0]
        
        G_IS[3] = Gf[4] + Gf[0]
        G_FS[3] = Gf[3] + site_energies[0]
        G_TS[3] = Gf[14] + site_energies[0]
        
        G_IS[4] = Gf[3] + site_energies[1]
        G_FS[4] = Gf[5] + Gf[7]
        G_TS[4] = Gf[9] + site_energies[1]
        
        G_IS[5] = Gf[5] + Gf[0]
        G_FS[5] = Gf[1] + site_energies[0]
        G_TS[5] = Gf[8] + site_energies[0]
        
        G_IS[6] = Gf[1] + Gf[0]
        G_FS[6] = Gf[2] + site_energies[0]
        G_TS[6] = Gf[10] + site_energies[0]
        
        G_IS[7] = Gf[2] + Gf[0]
        G_FS[7] = gas_energies[0] + site_energies[1] + site_energies[0]
        G_TS[7] = Gf[11] + site_energies[0]
        
        G_IS[8] = Gf[7] + Gf[0]
        G_FS[8] = gas_energies[2] + site_energies[1] + site_energies[0]
        G_TS[8] = Gf[13] + site_energies[0]
        
    
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
    
    
    kf, kr, junk, junk= rate_constants(rxn_parameters,theta,gas_energies,site_energies,T,F,
            mpf,matrix,mpexp,include_derivatives=False)
   
    r = [0]*len(kf)
    dtheta_dt = [0]*len(theta)
    
    s = [0]*3
    s[0] = (mpf('0.330000000000000015543122344752191565930843353271484375') - theta[0])
    s[1] = (mpf('1.0') - theta[1] - theta[2] - theta[3] - theta[4] - theta[5] - theta[6] - theta[7])
    s[2] = (mpf('0.0'))
    r[0] = kf[0]*p[3]*s[0]*s[0] - kr[0]*theta[0]*theta[0]
    r[1] = kf[1]*p[1]*s[1] - kr[1]*theta[6]
    r[2] = kf[2]*theta[6]*theta[0] - kr[2]*theta[4]*s[0]
    r[3] = kf[3]*theta[4]*theta[0] - kr[3]*theta[3]*s[0]
    r[4] = kf[4]*theta[3]*s[1] - kr[4]*theta[5]*theta[7]
    r[5] = kf[5]*theta[5]*theta[0] - kr[5]*theta[1]*s[0]
    r[6] = kf[6]*theta[1]*theta[0] - kr[6]*theta[2]*s[0]
    r[7] = kf[7]*theta[2]*theta[0] - kr[7]*p[0]*s[1]*s[0]
    r[8] = kf[8]*theta[7]*theta[0] - kr[8]*p[2]*s[1]*s[0]
    dtheta_dt[0] =  + 2*r[0] + -1*r[2] + -1*r[3] + -1*r[5] + -1*r[6] + -1*r[7] + -1*r[8]
    dtheta_dt[1] =  + 1*r[5] + -1*r[6]
    dtheta_dt[2] =  + 1*r[6] + -1*r[7]
    dtheta_dt[3] =  + 1*r[3] + -1*r[4]
    dtheta_dt[4] =  + 1*r[2] + -1*r[3]
    dtheta_dt[5] =  + 1*r[4] + -1*r[5]
    dtheta_dt[6] =  + 1*r[1] + -1*r[2]
    dtheta_dt[7] =  + 1*r[4] + -1*r[8]
 
    r = matrix(r)
    dtheta_dt = matrix(dtheta_dt)

    return dtheta_dt



