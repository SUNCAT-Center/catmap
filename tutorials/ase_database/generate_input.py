# -*- coding: utf-8 -*-
"""
Title: generate_input.py

Description
-----------
    Creates CATMAP input table from ase databases.

Input
----
    file_name : string
        Output filename, e.g. my_input.txt


Dependencies
------------
        db file which must include the key-value fields:
        "epot" with float value of DFT calculated electronic energy.

    <name>.db: db file
        file for molecules.
        Must include the key-value-pairs:
            "name" : str
                identifies the catalyst.
            "energy" or "epot" : float
                calculated potential energy.

    <name>.db: db file
        file for surfaces and adsorbate structures.
        Must include the key-value-pairs:
            "name" : str
                identifies the catalyst.
            "energy" or "epot" : float
                calculated potential energy.
            "species" : str
                adsorbate chemical formula.
        Recommended key-value-pairs:
            "phase" : str
                identifies the catalyst crystal phase.
            "facet" : str 
                identifies the site.
            "supercell" : str
                formatted XxY
                identifies the supercell.
            "layers" : int
                number of layers in the slab.

Output
------
    <name>.txt : text file
        CATMAP input file.

Example
-------
    python generate_input.py my_input.txt
"""
from sys import argv
from catmap.api.ase_data import energy_landscape

# Initialize thermodynamics module.
project = energy_landscape()

# Step 1: Define search filters. These are needed to select comparable data.
# They are typically calculation parameters, supercell sizes,
# number of layers in slabs
# Search strings for molecules.
mol_select = ['vacuum=8', 'fmaxout<0.05', 'pw=500']

# Step 2: Import from an ase database.
# Import molecules.
print('Importing molecules.')
project.get_molecules('molecules.db', selection=mol_select, frequency_db='frequencies.db')

# Search strings for all slabs, adsorbates.
fixed_p = ['ads!=FBL', 'ads!=NEB', 'layers=5', 'pw=500', 'psp=gbrv1.5pbe',
           'site!=off', 'C<3']

# Search strings for subsets of slabs, adsorbates.
surfaces1 = ['facet=(111)', 'phase=fcc', 'kpts=4x6', 'supercell=3x2'] + fixed_p
surfaces2 = ['facet=(001)', 'phase=hcp', 'kpts=4x6', 'supercell=3x2',
             'name!=Co'] + fixed_p
surfaces3 = ['surf_lattice=hexagonal', 'kpts=6x6', 'supercell=1x1'] + fixed_p

# Import three different subsets of slabs, adsorbates.
print('Importing surfaces.')
project.get_surfaces('surfaces.db', selection=surfaces1, frequency_db='frequencies.db')
project.get_surfaces('surfaces.db', selection=surfaces2)
project.get_surfaces('surfaces.db', selection=surfaces3)

# Get transition states.
# dehydro.get_transition_states('neb.db')

# Step 3: Calculate formation energies of adsorbates and transition states.
# refences is an optional parameter,
# which defines the energy references for each element.
project.calc_formation_energies(references=(('H', 'H2_gas'),
                                            ('O', 'H2O_gas'),
                                            ('C', 'CH4_gas'),))
# The defaults are as hydrogen, water and methane.

# Step 4: Save catmap input file.
file_name = 'my_input.txt'
project.make_input_file(file_name, site_specific='facet')
