# J. Montoya - making input file for catmap, CO2 reduction to CO for Cu, Ag, Au
import pickle
import os,sys
from ase.atoms import string2symbols
import numpy as np

edir = '/Volumes/suncatnfs/montoyjh/usr/pylib/hori/data/electronic-energies/'
vdir = '/Volumes/suncatnfs/montoyjh/usr/pylib/hori/data/generic-vibrations/'

freqconv = 1.239842e-4 # Frequency conversion factor eV/wavenumber

metals = ['Au', 'Al', 'Cd', 'Cu','Rh','Zn', 'Ag','Pt','Pd','Ni']
gases = ['CO', 'CO2', 'H2O','C2H4', 'H2']
adsorbates = ['COOH', 'CO','CHO','OCCHO','OH','O','OCHCHO', 'H']
              #TSes = ['OC-CHO', 'H-CO']

vibdofdict = {'CO':1,
              'CO2':4,
              'H2O':3,
              'C2H4':12,
              'H2':1,}

pcorr = {'CO2':+0.45,
         'CO' :+0.0,
         'H2O':+0.0,
         'H2':+0.0,
         'C2H4':+0.0} # Peterson correction for gas phase molecules
f = open('CO2_red_input.txt','w')
import string
f.write(string.join(['surface_name','site_name','species_name','formation_energy','bulk_structure','frequencies','other_parameters','reference\n'],'\t'))
refdict,tempdict = {},{}

refdict['H'] = 0.5*float(pickle.load(open(edir+'H2'))['electronic energy'])
refdict['O'] = float(pickle.load(open(edir+'H2O'))['electronic energy'])-2*refdict['H']
refdict['C'] = float(pickle.load(open(edir+'CH4'))['electronic energy'])-4*refdict['H']

beta = 0.5 # Symmetry factor for transition states

# Store H as proton electron pair with voltage correction,note that the free energy correction here is going to be weird... might have to do it manually
# V = float(sys.argv[1])
# f.write(string.join(['None','gas','pe',"0.",'None','['+str(0.5473/freqconv)+']','[]','gas phase calcs\n'],'\t'))
f.write(string.join(['None','gas','pe',"0.",'None','[]','[]','gas phase calcs\n'],'\t'))

# change voltage
mkm = open('CO2_reduction_template.mkm','rt')
text = mkm.read()
mkm.close()
mkm = open('CO2_reduction.mkm','wt')
new_text = text.replace('voltage = -0.5', 'voltage = '+str(sys.argv[1]))
mkm.write(new_text)
mkm.close()

for gas in gases:
    gasdict = pickle.load(open(edir+gas))
    e = float(gasdict['electronic energy']+pcorr[gas])
    formula = gas.replace('-','')
    composition = string2symbols(formula)
    for atom in composition:
        e -= refdict[atom]
    e = round(e,3)
    # Used to include frequencies for gas phase 
    vibstring = ','.join([str(np.real(gasdict['vibrations'][-i-1])/freqconv) for i in range(vibdofdict[gas])])
    freqstring = '['+vibstring+']'
    f.write(string.join(['None','gas',gas,str(e),'None',freqstring,'[]','gas phase calcs\n'],'\t'))
    tempdict[gas+'_g'] = e

for metal in metals:
    surf = float(pickle.load(open(edir+metal+'-fcc211_4x3'))['electronic energy'])
    for ads in adsorbates+['OC-CHO','OHC-CHO']:
        vibdict = pickle.load(open(vdir+ads))
        e = float(pickle.load(open(edir+ads+'_'+metal+'-fcc211_4x3'))['electronic energy'])
        e-=surf
        formula = ads.replace('-','')
        composition = string2symbols(formula)
        for atom in composition:
            e -= refdict[atom]
        e = round(e,3)
        vibdof = 3*len(composition) # 3N vibrational degrees of freedom
        vibstringlist = [str(np.real(vibdict['vibrations'][-i-1])/freqconv) for i in range(vibdof)]
        vibstringlist = filter(lambda a: a != '0.0',vibstringlist)
        vibstring = ','.join(vibstringlist)
        freqstring = '['+vibstring+']'
        f.write(string.join([metal,'211',ads,str(e),'fcc',freqstring,'[]','my ads calcs\n'],'\t'))
        tempdict[ads] = e
        print string.join([metal,'211',ads,str(e),'fcc',freqstring,'[]','my ads calcs\n'],'\t')
        print metal+' '+ads

    # Transition states for electrochemical steps
    barrier_CHO = 0.5
    barrier_O = 0.25
    barrier_OH = 0.25
    U0_cho_co = -(tempdict['CHO'] - tempdict['CO'] - 0.)
    e_cho_co = tempdict['CO'] + barrier_CHO # + beta*(V - U0_cho_co)
    #print metal+':'+str(U0_co2_cooh)
    f.write(string.join([metal,'211','pe-CO',str(e_cho_co),'None','[]','[]','my ads calcs\n'],'\t'))
    
    U0_oh_o = -(tempdict['OH'] - tempdict['O'] - 0.)
    e_oh_o = tempdict['O'] + barrier_O # + beta*(V - U0_oh_o) 
    f.write(string.join([metal,'211','O-pe',str(e_oh_o),'None','[]','[]','my ads calcs\n'],'\t'))

    U0_h2o_oh = -(tempdict['H2O_g'] - tempdict['OH'] - 0.)
    e_h2o_oh = tempdict['OH'] + barrier_OH # + beta*(V - U0_h2o_oh) 
    f.write(string.join([metal,'211','pe-OH',str(e_h2o_oh),'None',freqstring,'[]','my ads calcs\n'],'\t'))
    U0_co2_cooh = -(tempdict['COOH'] - tempdict['CO2_g'] - 0.)
    e_oco_h = tempdict['COOH'] + barrier_O # + beta*(V - U0_co2_cooh)
    #print metal+':'+str(U0_co2_cooh)
    f.write(string.join([metal,'211','OCO-pe',str(e_oco_h),'None','[]','[]','my ads calcs\n'],'\t'))
    U0_cooh_co = -(tempdict['CO'] + tempdict['H2O_g'] - tempdict['COOH'] - 0.)
    e_oc_hoh = tempdict['COOH'] + barrier_OH # + beta*(V - U0_cooh_co) 
    f.write(string.join([metal,'211','CO-peOH',str(e_oc_hoh),'None','[]','[]','my ads calcs\n'],'\t'))
    # Transition states for desorption
    e_co_d = tempdict['CO_g']
    vibdict = pickle.load(open(vdir+'CO'))
    vibstring = ','.join([str(np.real(vibdict['vibrations'][-i-1])/freqconv) for i in range(6)])
    freqstring = '['+vibstring+']'
    f.write(string.join([metal,'211','CO-',str(e_co_d),'None',freqstring,'[]','my ads calcs\n'],'\t'))
f.close()
