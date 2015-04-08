
def interaction_function(coverages,energies,epsilon,F,include_derivatives=True):

    ##this dictionary is passed in via the "template" so that it can be compiled
    site_info_dict = {
                      'h': [ #each site has an entry with the static information: 
                          [0], #the indexes of all adsorbates on this site
                          0.33,#the maximum allowed coverage
                          {'cutoff': 0.62, 'max_coverage': 0.4, 'offset': 0.7, 'smoothing': 0.05} # the parameters of F for this site
                          ], 
                      's': [
                          [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14], 
                          1, 
                          {'cutoff': 0.62, 'max_coverage': 0.4, 'offset': 0.7, 'smoothing': 0.05}]
                      }

    N_ads = len(coverages)
    N_sites = len(site_info_dict)

    idx_lists = []
    f = []
    df = []
    d2f = []

    ##for the "multisite" model this should yield f,df,d2f matrices rather than vectors
    for s in site_info_dict:
        idxs,max_cvg,F_params = site_info_dict[s]
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
   
    ## This is just an implementation of the indicial notation expression
    ## for the differential energy
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

    for s in range(N_sites):
        for q in range(N_sites):
            for n in idx_lists[s]:
                term_2[n] += df[s]*f[q]*term_2_sum[q][s]
                term_1[n] += f[s]*f[q]*term_1_sum[n][q]
                    
    E_diff = [a+b+c for a,b,c in zip(energies,term_1,term_2)]

    ## When calculating an analytical jacobian the derivatives of the
    ## differential energy are needed.
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
