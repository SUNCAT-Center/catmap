"""
API for catalysis-hub.org.

For building of micro-kinetic models using catalysis-hub data.

This module serves to build an EnergyLandscape object, importing the necessary
data from catalysis-hub.org.
"""
import warnings
import io
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
    def __init__(self):
        # Address of catalysis-hub.
        self.root = "https://api.catalysis-hub.org/graphql"
        self.protocol_prefix = 'postgresql://'
        self.port = '5432'
        self.path = '/catalysishub'
        self.database_address = 'catalysishub.c8gwuc8jwb7l' + \
                                '.us-west-2.rds.' + \
                                'amazonaws.com' + ':' + self.port + self.path

    def publication_energy_landscape(self, publication, references,
                                     site_specific=False, limit=100):
        """Return CatMAP EnergyLandscape object with potential energies
        from a publication query.

        Parameters
        ----------
        publication : str
            Title of publication.
        references : list
            List of gas phase references, e.g. ['H2', 'H2O', 'CH4']
        site_specific : bool
            TRUE distinguishes between adsorbate sites.
        limit : int
            Maximum number of formation energies to return.

        Returns
        ----------
            EnergyLandscape : object
                CatMAP EnergyLandscape object.
        """
        pub_string = """ title: \"""" + publication + """\" """ + \
            """ first: """ + str(limit) + """ """
        pubid = self.find_pubid(pub_string)

        ref_q = '~' + ' + ~'.join(references)

        query_string = """ reactants: \" """ + ref_q + """ \" pubId: \"""" + \
            pubid + """\" first: """ + str(limit) + """ """

        reactions = self.get_reactions(query_string)
        energy_landscape = self.attach_reaction_energies(reactions)

        return energy_landscape

    def get_publication_atoms(self, publication, limit=100):
        """
        Parameters
        ----------
        publication : str
            Title of publication.
        limit : int
            Maximum number of atoms objects to return.

        Returns
        ----------
        images : list
            List of ASE Atoms objects.
        """
        pub_string = """ title: \"""" + publication + """\" """ + \
            """ first: """ + str(limit) + """ """

        unique_ids = self.get_publication_uids(pub_string)
        images = self.get_atoms_from_uids(unique_ids, limit=limit)

        return images

    def attach_reaction_energies(self, reactions, energy_landscape=None,
                                 site_specific=False):
        """Return CatMAP EnergyLandscape object with formation energies.

        Parameters
        ----------
        reactions : str
            reactions GraphQL string.
        energy_landscape : object
            CatMAP EnergyLandscape object.
        site_specific : bool
            If True: Dinstinguish sites using the site key value pair, and
            stores a the potential energy of adsorbates on each site.
            Else: Use the minimum ab initio energy, disregarding the site.

        Returns
        ----------
            CatMAP EnergyLandscape object.
        """
        if energy_landscape is None:
            energy_landscape = EnergyLandscape()

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
            reaction += """, first: 100 """
            warnings.warn("Appending 'first: 100' to query." +
                          "A limit is mandatory")

        # GraphQL query to reactions table.
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
                      reactionSystems {
                              name
                              aseId
                              }
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

    def find_pubid(self, publication):
        """Return a publication id (pubid) from a catalysishub
        publication query.

        Parameters
        ----------
        publication : str
            GraphQL search string to the publications table.
        """
        if 'first' not in publication:
            publication += """, first: 100 """
            warnings.warn("Appending 'first: 100' to query." +
                          "A limit is mandatory")

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

    def get_atoms_from_uids(self, unique_ids, limit=100):
        """Return a list of atoms objects.

        Parameters
        ----------
        unique_ids : list
            Unique id's of atoms objects in the catalysishub database.
        limit : int
            Maximum number of atoms objects to return.
        """

        images = []
        for uid in tqdm(unique_ids[:limit]):
            query = \
                """
                query{
              systems(uniqueId: \"""" + uid + """\") {
                      totalCount
                      edges {
                              node {
                                      id
                                      uniqueId
                                      energy
                                      Trajdata
                                      }
                              }
                              }
                              }
                """
            systems = requests.post(self.root, {'query': query}).json()

            for j, _ in enumerate(systems['data']):
                with io.StringIO() as tmp_file:
                    system = systems['data']['systems']['edges'][j].pop('node')
                    atoms_dict = system.pop('Trajdata')
                    tmp_file.write(atoms_dict)
                    tmp_file.seek(0)
                    atoms = ase.io.read(tmp_file, format='json')
                    atoms.info['key_value_pairs'] = {}
                    atoms.info['key_value_pairs']['cathub_id'] = \
                        system.pop('id')
            images.append(atoms)
        return images

    def get_publication_uids(self, publication):
        """Return a list of ASE-db unique_id's from a catalysishub query.

        Parameters
        ----------
        reaction : str
            GraphQL search string to the publications table.
        """
        if 'first' not in publication:
            publication += """, first: 100 """
            warnings.warn("Appending 'first: 100' to query." +
                          "A limit is mandatory")

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

        return unique_ids
