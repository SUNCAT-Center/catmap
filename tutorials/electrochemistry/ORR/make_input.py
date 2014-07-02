from chuan_common import get_binding_energy_from_hori, get_gas_energy_from_hori, get_vibrations_from_hori
import sys
import numpy as np
from ase.all import string2symbols
import string

freqconv = 1.239842e-4 # Frequency conversion factor eV/wavenumber

gases = ['H2O', 'O2', 'H2', 'H2O2']
gas_energy_dict = {
    'H2O':0.0,
    'H2':0.0,
    'O2':5.19,  # could be 1.59 + 0.9*4 or 4.92
    'H2O2':1.86 + 0.9*2,
}
ads_energy_dict = {
    'O2_a':1.39 + 0.9*4,
    'OOH_a':1.21 + 0.9*3,
    'O_a':-0.10 + 0.9*2,
    'OH_a':-0.15 + 0.9,
    'H2O2_a':1.53 + 2*0.9,
    'O_b':0.53 + 0.9*2,
    'OH_b':0.20 + 0.9,
    '*_a':0.,
    '*_b':0.,
    'H_a':-0.3,
}
ads_squiggly_dict = {
    'O2_a':1.2,
    'OOH_a':1,
    'O_a':2,
    'OH_a':1,
    'H2O2_a':0.004,
    'O_b':2,
    'OH_b':1,
    'H_a':0.001,
}
OH_BEs = [-0.15
-0.1,
-0.05,
0.0,
0.05,
0.10,
0.15,
0.20,
0.25,
0.30,
]
labels = string.ascii_lowercase[:len(OH_BEs)]


step_to_TS_name = {
    '2':'O2-_a',
    '3':'OO-pe_a',
    '4':'OOH-pe_a',
    '5':'O-pe_a',
    '6':'OH-pe_a',
    '7':'O-O_a',
    '8':'O-pe_b',
    '9':'OH-pe_b',
    '10':'OH-O_a',
    '11':'O-peOH_a',
    '12':'OH-OH_a',
    '13':'HOOH-_a',
}
chemical_TS_energy_dict = {
    '2':gas_energy_dict['O2'] + 0.4076,  # from the 8 * 10^5 prefactor of O2(aq) -> O2(dl)
    '7':ads_energy_dict['O2_a'] + 0.48,
    '10':ads_energy_dict['OOH_a'] + 0.37,
    '12':ads_energy_dict['H2O2_a'] + 0.46,
    '13':0.32 + 0.284,  # dG_Pt + 1 * 10^8 prefactor
}
chemical_TS_gamma_dict = {
    '2':0.,
    '7':0.69,
    '10':0.39,
    '12':0.19,
    '13':0.,
}

echem_TS_effective_barrier = 0.26 + 0.22  # tripkovic barrier + water reorganization 'barrier'
echem_TS_ISs = {  # assumes pe is also in the initial state
    '3':'O2_a',
    '4':'OOH_a',
    '5':'O_a',
    '6':'OH_a',
    '8':'O_b',
    '9':'OH_b',
    '11':'OOH_a',
}

echem_TS_FSs = {
    '3':'OOH_a',
    '4':'O_a',
    '5':'OH_a',
    '6':'*_a',
    '8':'OH_b',
    '9':'*_b',
    '11':'H2O2_a',
}

tempdict = {}

f = open('ORR_input.txt','w')
lines = []

header = ['surface_name','site_name','species_name','formation_energy','bulk_structure','frequencies','other_parameters','reference']
lines.append(header)

pe = ['None','gas','pe',"0.",'None','[]','[]','gas phase calcs']
lines.append(pe)

for gas in gases:
    energy = gas_energy_dict[gas]
    vibrations = []
    lines.append(['None','gas',gas,str(energy),'None',str(vibrations),'[]','Hansen 2014'])
    tempdict[gas+'_g'] = energy

assert(len(labels) == len(OH_BEs))
for i, label in enumerate(labels):
    surface = 'Hansenite_' + label
    for ads in ads_energy_dict.keys():
        ads_species, facet = ads.split('_')
        if ads_species != '*':
            ads_energy = ads_energy_dict[ads] + ads_squiggly_dict[ads]*OH_BEs[i]
            lines.append([surface, facet, ads_species, str(ads_energy), 'fcc', '[]', '[]', 'Hansen 2014'])
        else:
            ads_energy = 0.
        tempdict[ads] = ads_energy

    for chem_TS in chemical_TS_energy_dict.keys():
        energy = chemical_TS_energy_dict[chem_TS] + chemical_TS_gamma_dict[chem_TS]*OH_BEs[i]
        split_up_identity = step_to_TS_name[chem_TS].split('_')
        facet = split_up_identity[-1]
        species = split_up_identity[0]
        lines.append([surface, str(facet), species, str(energy), 'fcc', '[]', '[]', 'Hansen 2014'])

    for echem_TS in echem_TS_ISs.keys():
        IS_energy = tempdict[echem_TS_ISs[echem_TS]]
        FS_energy = tempdict[echem_TS_FSs[echem_TS]]
        higher = max(IS_energy, FS_energy)
        e_TS = higher + echem_TS_effective_barrier
        split_up_identity = step_to_TS_name[echem_TS].split('_')
        facet = split_up_identity[-1]
        species = split_up_identity[0]
        lines.append([surface, str(facet), species, str(e_TS), 'fcc', '[]', '[]', 'Hansen 2014'])


for line in lines:
    to_write = "\t".join(line) + "\n"
    f.write(to_write)

f.close()
