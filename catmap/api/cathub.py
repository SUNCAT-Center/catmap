"""
API for catalysis-hub.org

Requires a development branch of ASE: https://gitlab.com/ase/ase/tree/database
"""
import requests
import ase.db
from catmap.api.ase_data import EnergyLandscape
try:
    from tqdm import tqdm
except ImportError:
    def tqdm(interable):
        return interable
import ast


class CatalysisHub(object):
    """ API for importing energetic data from Catalysis-hub.org """
    def __init__(self, username=None, password=None, limit=100):
        self.root = "https://api.catalysis-hub.org/graphql"

        self.protocol_prefix = 'postgresql://'
        self.port = '5432'
        self.path = '/catalysishub'
        self.database_address = 'catalysishub.c8gwuc8jwb7l' + \
                                '.us-west-2.rds.' + \
                                'amazonaws.com' + ':' + self.port + self.path
        self.username = username
        self.password = password

        self.limit = limit
        self.c = None

    def get_reactions(self, reaction):
        """Return a list of reaction dictionaries from a catalysishub query.

        Parameters
        ----------
        reaction : str
            GraphQL search string. This search string for the reactions table.
            The string must be enclosed in three double hyphens and additional
            double hyphens must be preceeded by a backslash.
        """
        if 'first' not in reaction:
            reaction += """, first: """ + str(self.limit)

        query = \
            """
            query{reactions
            """ + \
            """(""" + reaction + """) {
                totalCount
                edges {
                  node {
                      id
                      reactants
                      products
                      facet
                      sites
                      reactionEnergy
                      activationEnergy
                      surfaceComposition
                    }
                  }
                }
              }
            """
        data = requests.post(self.root, {'query': query}).json()
        if 'errors' in data:
            raise KeyError(data['errors'])

        reactions = data['data']['reactions']['edges']

        return reactions

    def get_publication_id(self, publication):
        """Return a dictionary from a catalysishub query."""
        if 'first' not in publication:
            publication += """, first: """ + str(self.limit)

        query = \
            """{
              publications(""" + publication + """) {
                edges {
                  node {
                    pubId
                  }
                }
              }
            }
            """
        data = requests.post(self.root, {'query': query}).json()
        if 'errors' in data:
            raise KeyError(data['errors'])

        pub_id = data['data']['publications']['edges'][0]['node']['pubId']
        return pub_id

    def attach_reaction_energies(self, reactions,
                                 previous=None, site_specific=False):
        """Return CatMAP energy_landscape object with formation energies.

        Parameters
        ----------
        fname : str
            Path and filename of candidate ase database file.
        database_ids : list
            Database ids.
        prediction : list
            Predicted means in the same order as database_ids.
        uncertainty : list
            Predicted uncertainties in the same order as database_ids.
        catmap : object
            CatMAP energy_landscape object.
        site_specific : bool
            If True: Dinstinguish sites using the site key value pair, and
            stores a the potential energy of adsorbates on each site.
            Else: Use the minimum ab initio energy, disregarding the site.
        """
        if previous is None:
            energy_landscape = EnergyLandscape()
        else:
            energy_landscape = previous

        for rxn in reactions:
            node = rxn['node']
            rxn_id = node['id']
            name = node['surfaceComposition']
            facet = node['facet']
            fstate = ast.literal_eval(node['products'])
            sites = ast.literal_eval(node['sites'])
            for product in fstate.keys():
                if product != 'star' and 'star' in product:
                    species = product.rstrip('star')
                    n = fstate[product]
                    site = sites[species]
                    break

            key = '_'.join(['1', species, name, 'phase', 'surf',
                            facet, 'cell', site])
            if key not in energy_landscape.formation_energies:
                energy_landscape.formation_energies[key] = \
                    node['reactionEnergy'] / n
                energy_landscape.dbid[key] = rxn_id
            elif (energy_landscape.formation_energies[key] >
                  node['reactionEnergy'] / n):
                energy_landscape.formation_energies[key] = \
                    node['reactionEnergy'] / n
                energy_landscape.dbid[key] = rxn_id

        return energy_landscape

    def get_atoms(self, unique_ids, limit=100):
        """Return a list of atoms objects.

        Parameters
        ----------
        unique_ids : list
            Unique id's of atoms objects in the catalysishub database.
        """
        if self.c is None:
            self.c = ase.db.connect(self.protocol_prefix +
                                    self.username + ':' +
                                    self.password + '@' +
                                    self.database_address)

        if limit is not None and len(unique_ids) > limit:
            unique_ids = unique_ids[:limit]

        images = []
        for uid in tqdm(unique_ids):
            atoms = self.c.get_atoms(['unique_id=' + uid],
                                     add_additional_information=True)
            images.append(atoms)
        return images

    def get_publication_atoms(self, publication, limit=10):
        """Return a list of ASE-db unique_id's from a catalysishub query.

        Parameters
        ----------
        reaction : str
            GraphQL search string. This search string for the publications
            table. The string must be enclosed in three double hyphens and
            additional double hyphens must be preceeded by a backslash.
        """
        query = \
            """{
              publications(""" + publication + """) {
                edges {
                  node {
                    systems {
                      uniqueId
                    }
                  }
                }
              }
            }
            """

        data = requests.post(self.root, {'query': query}).json()
        if 'errors' in data:
            raise KeyError(data['errors'])

        unique_ids = []
        for dat in data['data']['publications']['edges'][0]['node']['systems']:
            unique_ids.append(str(dat['uniqueId']))

        return self.get_atoms(unique_ids, limit=limit)

    def import_energies(self, unique_ids, previous=None, site_specific=False):
        """Return CatMAP energy_landscape object with potential energies.

        Parameters
        ----------
        fname : str
            Path and filename of candidate ase database file.
        database_ids : list
            Database ids.
        prediction : list
            Predicted means in the same order as database_ids.
        uncertainty : list
            Predicted uncertainties in the same order as database_ids.
        catmap : object
            CatMAP energy_landscape object.
        site_specific : bool
            If True: Dinstinguish sites using the site key value pair, and
            stores a the potential energy of adsorbates on each site.
            Else: Use the minimum ab initio energy, disregarding the site.
        """
        if previous is None:
            energy_landscape = EnergyLandscape()
        else:
            energy_landscape = previous

        if self.c is None:
            self.c = ase.db.connect(self.protocol_prefix +
                                    self.username + ':' +
                                    self.password + '@' +
                                    self.database_address)

        for uid in tqdm(unique_ids):
            d = self.c.get(['unique_id=' + uid])
            if str(d.state) == 'gas':
                key = str(d.formula) + '_gas'
            else:
                try:
                    n, species, name, phase, surf_lattice, facet, cell = \
                        energy_landscape._get_adsorbate_fields(d)
                except AttributeError:
                    for key in d:
                        print(key, d[key])
                    raise
                # layers are inconsistent in catalysishub.
                x, y = cell[0], cell[2]
                cell = 'x'.join([x, y])
                if str(d.state) == 'star':
                    site = 'slab'
                elif site_specific and 'site' in d:
                    site = str(d.site)
                else:
                    site = 'site'
                key = '_'.join([str(n), species, name, phase, surf_lattice,
                                facet, cell, site])
            epot = float(d.energy)
            if key not in energy_landscape.epot:
                energy_landscape.epot[key] = epot
                energy_landscape.dbid[key] = uid
            elif energy_landscape.epot[key] > epot:
                energy_landscape.epot[key] = epot
                energy_landscape.dbid[key] = uid

        return energy_landscape
