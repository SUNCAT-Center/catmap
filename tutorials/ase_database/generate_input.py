# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 10:47:15 2016

@author: mhangaard

Title: db2catmap_input

Description:
    Creates CATMAP input table from ase db files

Input:
    1) string: <name>


Dependencies:
    mol.db: db file which must include the key-value fields:
        "enrgy" with float value of DFT calculated electronic energy.
        "vacuum" with float value specifying amount of vacuum.
        "pw" float with plane wave cutoff.

    <name>.db: db file, which must include the key-value fields:
        "name" with str value that identifies the catalyst stochiometry.
        "phase" with str value that identifies the catalyst phase.
        "facet" with str value that identifies the site.
        "series" with string value of 'slab' or adsorbate chemical formula.
        "enrgy" with float value of DFT calculated electronic energy.
        "pw" float with plane wave cutoff.
        "kpts" str with k-points in the plane dimensions separated by 'x'


Output:
    CATMAP input file: <name>.txt

Example:
    python db2catmap_input therm_dft
    will produce "therm_dft.txt"
"""
from sys import argv
from catmap.asedb2catmap import db2catmap

# Initialize thermodynamics module.
dehydro = db2catmap()

# Define search filters.

# Search strings for molecules.
mol_select = ['vacuum=8', 'fmaxout<0.05']

# Import molecules from ase database.
dehydro.get_molecules('molecules.db', selection=mol_select)

# Search strings for all slabs, adsorbates.
fixed_p = ['ads!=FBL', 'ads!=NEB', 'layers=5', 'pw=500', 'psp=gbrv1.5pbe',
           'site!=off', 'C<3']

# Search strings for subsets of slabs, adsorbates.
surfaces1 = ['facet=1x1x1', 'phase=fcc', 'kpts=4x6', 'supercell=3x2'] + fixed_p
surfaces2 = ['facet=0x0x1', 'phase=hcp', 'kpts=4x6', 'supercell=3x2',
             'name!=Co'] + fixed_p
surfaces3 = ['surf_lattice=hexagonal', 'kpts=6x6', 'supercell=1x1'] + fixed_p

# Import three different subsets of slabs, adsorbates.
dehydro.get_surfaces('surfaces.db', selection=surfaces1)
dehydro.get_surfaces('surfaces.db', selection=surfaces2)
dehydro.get_surfaces('surfaces.db', selection=surfaces3)

# Get transition states.
# dehydro.get_transition_states('fbl.db', selection=['Re=0'])
# dehydro.get_transition_states('fbl.db', selection=['Re', 'species!=CH2CH2-H'])

# Calculate formation energies with custom elemental references.
dehydro.calc_formation_energies(references={'H': 'H2_gas',
                                            'O': 'H2O_gas',
                                            'C': 'CH4_gas'})

# Save catmap input file.
file_name = argv[1] + '.txt'
dehydro.make_input_file(file_name)
