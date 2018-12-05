
# coding: utf-8

# ## How to import the data
# 
#     1. Define search filters. This is needed if some data has to be filtered out.
#     2. Import data from ase databases.
#     3. Store references and calculate formation energies.
#     4. Export catmap input file.
# 
# It is important to pay attention to the search filters. If you get garbage results, it is likely because the search filters are not sufficient for your dataset. Make sure you filter calculator parameters such as XC-functional, basis sets cutoffs, k-point sampling, ect., when necessary.
# 
# Importing data from correctly formatted .db files:

# In[ ]:


# Import and instantiate energy_landscape object.
from catmap.api.ase_data import EnergyLandscape
energy_landscape = EnergyLandscape()

# Import all gas phase species from db.
search_filter_gas = []
energy_landscape.get_molecules('molecules.db', selection=search_filter_gas)

# Import all adsorbates and slabs from db.
search_filter_slab = []
energy_landscape.get_surfaces('surfaces.db', selection=search_filter_slab, site_specific=False)


# The `site_specific` option accepts `True`, `False` or a string. In the latter case, the `site` key is recognized only if the value matches the string, while all other sites are treated as identical.
# 
# Your data is now stored in dictionaries that are attached to your `EnergyLandscape` object.
# 
# ## Get formation energies and export to catmap format.
# 
# Formation energies are calculated:

# In[ ]:


references = (('H', 'H2_gas'), ('O', 'H2O_gas'), ('C', 'CH4_gas'),)
energy_landscape.calc_formation_energies(references)


# `references` is a required parameter, that should contain gas phase references. If a gas phase reference is dependent on another, order the dependent one after the latter.

# In[ ]:


file_name = 'my_input.txt'
energy_landscape.make_input_file(file_name)


# finally saves the input file, compatible with catmap.
# 
# ## How to import frequencies.
# 
# The field `data` always contains a dictionary. Use it with the key `frequencies` to make them accessible to the asedb2catmap module.
# 
# Importing frequencies is handled by the methods `get_surfaces` and `get_molecules`, which we have already used. It is necessary to pass the parameter `frequency_db` to it to import frequencies along with atomic structures like so:

# In[ ]:


energy_landscape.get_molecules('molecules.db', frequency_db='frequencies.db', selection=search_filter_gas)


# ## How to import transition states and pathways.
# 
# Transition states and paths have the mandatory key value pairs:
# 
#  - `path_id`
#  - `step` or `image`
# 
# `step` or `image` is used to order the images.
# There is one additional recommended key value pair:
# 
#  - `distance`
# 
# which is useful for making plots of the energy versus a reaction coordinate.
# To add formation energies of transition states to your catmap input, you can use the method:

# In[ ]:


energy_landscape.get_transition_states('neb.db')

