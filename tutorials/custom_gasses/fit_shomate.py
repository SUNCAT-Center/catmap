import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from catmap.thermodynamics import fit_shomate, ThermoCorrections
import catmap
import numpy as np

def H2O_shomate(output_file=None):
    Ts = np.array([100,200,298.15,300,400,500,600])
    Cps = np.array([33.299,33.349,33.59,33.596,34.262,35.226,36.325])
    Hs = np.array([-6.615,-3.282,0,0.062,3.452,6.925,10.501])
    Ss = np.array([152.388,175.485,188.834,189.042,198.788,206.534,213.052])
    params0 = [30.09200,6.832514,6.793435,
            -2.534480,0.082139,-250.8810,223.3967,-241.8264]
    params = fit_shomate(Ts,Cps,Hs,Ss,params0,output_file)
    params[-3] -= params[-1]
    params[-1] -= params[-1]
    return params

def CH3OH_shomate(output_file=None):
    """ Raw data from CRC handbook, 91st edition
    H is constrained to 0 since it can be lumped with F
    """
    Ts =  [298.15, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0, 1300.0, 1400.0, 1500.0]
    Ss =  [239.865, 240.139, 253.845, 266.257, 277.835, 288.719, 298.987, 308.696, 317.896, 326.629, 334.93, 342.833, 350.367, 357.558]
    Cps =  [44.101, 44.219, 51.713, 59.8, 67.294, 73.958, 79.838, 85.025, 89.597, 93.624, 97.165, 100.277, 103.014, 105.422]
    Hs =  [0.0, 0.082, 4.864, 10.442, 16.803, 23.873, 31.569, 39.817, 48.553, 57.718, 67.262, 77.137, 87.304, 97.729]
    params0 = [30.09200,6.832514,6.793435,
            -2.534480,0.082139,-250.8810,223.3967,-241.8264]
    params = fit_shomate(Ts,Cps,Hs,Ss,params0,output_file)
    params[-3] -= params[-1]
    params[-1] -= params[-1]
    return params

def CH3CH2OH_shomate(output_file=None):
    """Raw data from CRC handbook, 91st edition
    H is constrained to 0 since it can be lumped with F
    """
    Ts =  [298.15, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0]
    Ss =  [281.622, 282.029, 303.076, 322.75, 341.257, 358.659, 375.038, 390.482, 405.075, 418.892, 431.997]
    Cps =  [65.652, 65.926, 81.169, 95.4, 107.656, 118.129, 127.171, 135.049, 141.934, 147.958, 153.232]
    Hs =  [0.0, 0.122, 7.474, 16.318, 26.487, 37.79, 50.065, 63.185, 77.042, 91.543, 106.609]
    params0 = [30.09200,6.832514,6.793435,
            -2.534480,0.082139,-250.8810,223.3967,-241.8264]
    params = fit_shomate(Ts,Cps,Hs,Ss,params0,output_file)
    params[-3] -= params[-1]
    params[-1] -= params[-1]
    return params

def CH3CHO_shomate(output_file=None):
    """Raw data from CRC handbook, 91st edition
    H is constrained to 0 since it can be lumped with F
    """
    Ts =  [298.15, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0, 1300.0, 1400.0, 1500.0]
    Ss =  [263.84, 264.18, 281.62, 297.54, 312.36, 326.23, 339.26, 351.52, 363.1, 374.04, 384.4, 394.23, 403.57, 412.46]
    Cps =  [55.318, 55.51, 66.282, 76.675, 85.942, 94.035, 101.07, 107.19, 112.49, 117.08, 121.06, 124.5, 127.49, 130.09]
    Hs =  [0.0, 0.103, 6.189, 13.345, 21.486, 30.494, 40.258, 50.698, 61.669, 73.153, 85.065, 97.344, 109.954, 122.834]
    params0 = [30.09200,6.832514,6.793435,
            -2.534480,0.082139,-250.8810,223.3967,-241.8264]
    params = fit_shomate(Ts,Cps,Hs,Ss,params0,output_file)
    params[-3] -= params[-1]
    params[-1] -= params[-1]
    return params

def ideal_shomate_comparison(): 

    """ Compare ideal gas and shomate corrections
    """
    thermo = ThermoCorrections()
    shomate_gasses = [g.split(':')[0] for g in catmap.data.shomate_params.keys()]
    thermo.frequency_dict = catmap.data.experimental_gas_frequencies
    ideal_gasses = set(thermo.ideal_gas_params.keys()) & set(thermo.frequency_dict.keys())
    thermo.gas_names = list(set(shomate_gasses) & ideal_gasses)
    thermo.gas_pressures = [1]*len(thermo.gas_names)

    T_range = np.linspace(300,1000,100)
    err_dict = {}
    labels = thermo.gas_names
    for l in labels:
        err_dict[l] = []
    for T in T_range:
        thermo.temperature = T
        ideal_dict = thermo.ideal_gas()
        shomate_dict = thermo.shomate_gas()
        for key in thermo.gas_names:
            err = ideal_dict[key] - shomate_dict[key]
            err_dict[key].append(err)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel('Error (eV)')
    for label in labels:
        ax.plot(T_range,err_dict[label],label=label)
    plt.legend()
    fig.savefig('shomate_ideal_comparison.pdf')

print 'methanol',CH3OH_shomate('methanol_shomate.pdf')
print 'ethanol',CH3CH2OH_shomate('ethanol_shomate.pdf')
print 'acetaldehyde',CH3CHO_shomate('acetaldehyde_shomate.pdf')

ideal_shomate_comparison()
