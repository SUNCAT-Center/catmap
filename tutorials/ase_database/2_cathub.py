
# coding: utf-8

# # Import data from Catalysishub.org for CatMAP
# 
# In this tutorial we will
#     -   Download a set of formation energies from a publication and export them as a CatMAP data file.
#     -   Create an ASE-db sqlite3 file containing the corresponding atomic structures.

# In[ ]:


# Import modules
from os import environ

import ase.db
from ase.visualize import view

from catmap.api.cathub import CatalysisHub


# ### Download formation energies.
# 
# First we need to query the database and get the `pubId` for the publication we want to obtain data from.

# In[ ]:


# Instantiate cathub interface object.
cathub = CatalysisHub()

# GraphQL search string.
publication = """ title: \"WO3 surface doping for OER to be published\" """

# Return the pubId.
pubid = cathub.get_publication_id(publication)


# Once we have `pubId`, it is easy to query reactions limiting them to a publication.
# 
# Formation energies as used by CatMAP are simply reaction energies with respect to a fixed reference.
# Therefore, you only need to query for reactions from the relevant gas phase references, in order to download the relevant set of formation energies.

# In[ ]:


# GraphQL search string.
query_string = """ reactants: \"~H2Ogas + ~H2gas\" pubId: \"""" + pubid + """\" first: 10 """

reactions = cathub.get_reactions(query_string)
print(len(reactions))


# We have now retrieved a list of dictionaries, `reactions`. The reaction energies can be attached to a `catmap.api.energy_landscape` object as formation energies.

# In[ ]:


# Download and attach reaction energies.
energy_landscape = cathub.attach_reaction_energies(reactions, site_specific=True)


# Finally, as usual, we export a CatMAP input file.

# In[ ]:


energy_landscape.make_input_file(pubid + '.txt')


# ### Atomic structures.
# 
# Retrieving structures directly from the Catalysishub backend requires a login.
# 
# You can store the login details as environment variables by typing
# 
#     export CATHUB_USER=<username>
#     export CATHUB_PASSWORD=<password>
#     
# on the command line, or by adding it to your .bash_profile or .bashrc, or just type them in below.

# In[ ]:


# Retrieve the login details, assuming you have stored them in an environment variable.
username = environ['CATHUB_USER']
password = environ['CATHUB_PASSWORD']

# Instantiate the catmap.api.cathub.CatalysisHub class and store the login details in the class.
cathub = CatalysisHub(username=username, password=password)


# Next, we will search for a set of structures from a publication.

# In[ ]:


# Return a list of atoms objects.
images = cathub.get_publication_atoms(""" pubId: \"""" + pubid + """\" """, limit=10)


# In[ ]:


# Save them to an ASE-db file.
c = ase.db.connect('my_asedb.db')

for atoms in images:
    c.write(atoms, key_value_pairs=atoms.info['key_value_pairs'])

