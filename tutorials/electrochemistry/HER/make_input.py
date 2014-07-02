from chuan_common import get_binding_energy_from_hori, get_gas_energy_from_hori, get_vibrations_from_hori
import sys
import numpy as np
from ase.all import string2symbols

# change voltage
mkm = open('CO2_reduction_template.mkm','rt')
text = mkm.read()
mkm.close()
mkm = open('CO2_reduction.mkm','wt')
new_text = text.replace('voltage = -0.5', 'voltage = '+str(sys.argv[1]))
mkm.write(new_text)
mkm.close()

freqconv = 1.239842e-4 # Frequency conversion factor eV/wavenumber

metals = ['Au', 'Ag', 'Cu','Rh','Ni','Pt','Pd']
gases = ['CO', 'CO2', 'H2O', 'H2', 'HCOOH']
adsorbates = ['COOH', 'CO','OH','O', 'H', 'OCHO']
TSes = ['COO-pe', 'H-H', 'pe-H', 'OC-peO', 'HCOO-pe', 'O-pe', 'pe-OH', 'COOH-pe']
facet = 111

vibdofdict = {'CO':1,
              'CO2':4,
              'H2O':3,
              'C2H4':12,
              'H2':1,
              'HCOOH':9}

pcorr = {'CO2':+0.45,
         'CO' :+0.0,
         'H2O':+0.0,
         'H2':+0.0,
         'C2H4':+0.0,
         'HCOOH':+0.0} # Peterson correction for gas phase molecules

# H2 dissociation barriers on fcc(111) metals from CatApp
H2_dissoc_dict = {
    'Rh':0.25,
    'Ag':1.0,
    'Pd':0.12,
    'Pt':0.19,
    'Cu':0.63,
    'Au':1.1,
}

tempdict = {}

f = open('CO2_red_input.txt','w')
lines = []

header = ['surface_name','site_name','species_name','formation_energy','bulk_structure','frequencies','other_parameters','reference']
lines.append(header)

pe = ['None','gas','pe',"0.",'None','[]','[]','gas phase calcs']
lines.append(pe)

for gas in gases:
    energy = get_gas_energy_from_hori(gas) + pcorr[gas]
    vibrations = get_vibrations_from_hori(gas)

    vibstring = ','.join([str(np.real(vibrations[-i-1])/freqconv) for i in range(vibdofdict[gas])])
    freqstring = '['+vibstring+']'
    lines.append(['None','gas',gas,str(energy),'None',freqstring,'[]','gas phase calcs'])
    tempdict[gas+'_g'] = energy

for metal in metals:
    for ads in adsorbates:
        energy = get_binding_energy_from_hori(ads, metal, facet)
        if np.isnan(energy):
            continue
        vibrations = get_vibrations_from_hori(ads, metal, facet)

        formula = ads.replace('-','')
        composition = string2symbols(formula)
        energy = round(energy,3)
        vibdof = 3*len(composition) # 3N vibrational degrees of freedom
        vibstringlist = [str(np.real(vibrations[-i-1])/freqconv) for i in range(vibdof)]
        vibstringlist = filter(lambda a: a != '0.0',vibstringlist)
        vibstring = ','.join(vibstringlist)
        freqstring = '['+vibstring+']'
        lines.append([metal,str(facet),ads,str(energy),'fcc',freqstring,'[]','my ads calcs'])
        tempdict[ads] = energy

    # Transition states for electrochemical steps
    barrier_COO_COOpe = 0.2
    barrier_H_pe = 0.7
    barrier_COO_OCpeO = 0.8
    barrier_OCHO_HCOOpe = 0.2
    barrier_O_Ope = 0.2
    barrier_OH_peOH = 0.3
    barrier_COOH_COOHpe = 0.5

    # COO-pe
    higher = max(tempdict['CO2_g'] + 0.5*tempdict['H2_g'], tempdict['COOH'])  
    assert(tempdict['H2_g'] == 0.)
    e_ts = higher + barrier_COO_COOpe
    lines.append([metal,str(facet),'COO-pe',str(e_ts),'fcc','[]','[]','my ads calcs'])
    
    # pe-H
    higher = max(tempdict['H'] + 0.5*tempdict['H2_g'], tempdict['H2_g'])
    e_ts = higher + barrier_H_pe
    lines.append([metal,str(facet),'pe-H',str(e_ts),'fcc','[]','[]','my ads calcs'])

    # OC-peO
    higher = max(tempdict['CO2_g'] + 0.5*tempdict['H2_g'], tempdict['OCHO'])
    e_ts = higher + barrier_COO_OCpeO
    lines.append([metal,str(facet),'OC-peO',str(e_ts),'fcc','[]','[]','my ads calcs'])

    # HCOO-pe
    higher = max(tempdict['OCHO'] + 0.5*tempdict['H2_g'], tempdict['HCOOH_g'])
    e_ts = higher + barrier_OCHO_HCOOpe
    lines.append([metal,str(facet),'HCOO-pe',str(e_ts),'fcc','[]','[]','my ads calcs'])

    # O-pe
    higher = max(tempdict['O'] + 0.5*tempdict['H2_g'], tempdict['OH'])
    e_ts = higher + barrier_O_Ope
    lines.append([metal,str(facet),'O-pe',str(e_ts),'fcc','[]','[]','my ads calcs'])

    # pe-OH
    higher = max(tempdict['OH'] + 0.5*tempdict['H2_g'], tempdict['H2O_g'])
    e_ts = higher + barrier_OH_peOH
    lines.append([metal,str(facet),'pe-OH',str(e_ts),'fcc','[]','[]','my ads calcs'])

    # COOH-pe
    higher = max(tempdict['COOH'] + 0.5*tempdict['H2_g'], tempdict['H2O_g'] + tempdict['CO'])
    e_ts = higher + barrier_COOH_COOHpe
    lines.append([metal,str(facet),'COOH-pe',str(e_ts),'fcc','[]','[]','my ads calcs'])

    # H-H
    if metal in H2_dissoc_dict:
        e_ts = H2_dissoc_dict[metal] + tempdict['H2_g']
        lines.append([metal,str(facet),'H-H',str(e_ts),'fcc','[]','[]','CatApp'])        

for line in lines:
    to_write = "\t".join(line) + "\n"
    f.write(to_write)

f.close()
