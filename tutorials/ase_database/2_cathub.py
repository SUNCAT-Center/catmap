
# coding: utf-8

# # Import data from Catalysishub.org for CatMAP
# 
# In this tutorial we will
#     -   Download a set of formation energies from a publication and export them in a CatMAP EnergyLandscape object.
#     -   Create an ASE-db sqlite3 file containing the corresponding atomic structures.

# In[ ]:


# Import modules.
import os
import ase.db
from ase.visualize import view

from catmap.api.cathub import CatalysisHub


# ### Download formation energies.
# 
# First we need to write a query for a publication.

# In[ ]:


# Instantiate cathub interface object.
cathub = CatalysisHub()

# GraphQL search string to the publications table.
publication = "WO3 surface doping for OER to be published"


# Formation energies as used by CatMAP are simply reaction energies with respect to a fixed reference.
# Therefore, you only need to query for reactions from the relevant gas phase references, in order to download the relevant set of formation energies.

# In[ ]:


# Choose your references.
references = ['H2gas', 'H2Ogas']

# Fetch energies and create an EnergyLandscape object.
energy_landscape = cathub.publication_energy_landscape(publication, references, site_specific=True, limit=10)


# We have now retrieved a list of dictionaries, `reactions`. The reaction energies can be attached to a `catmap.api.energy_landscape` object as formation energies.

# Finally, as usual, we export a CatMAP input file.

# In[ ]:


fname = 'my_energies.txt'
energy_landscape.make_input_file(fname)


# In[ ]:


# Take a peak at the file.
with open(fname) as fp:
    for line in fp.readlines()[:5]:
        print(line)


# Notice the `reference` column contains catalysis-hub ids corresponding to the atomic structure.

# ### Atomic structures.

# Next, we will retrieve atomic structures from the publication.

# In[ ]:


# Return a list of atoms objects.
images = cathub.get_publication_atoms(publication, limit=10)


# Finally, we can save them to an ase database, keeping the catalysis-hub ids, to connect them with the energy data file.

# In[ ]:


# Save them to an ASE-db file.
os.remove('my_asedb.db')
c = ase.db.connect('my_asedb.db')

for atoms in images:
    c.write(atoms, key_value_pairs=atoms.info['key_value_pairs'])

