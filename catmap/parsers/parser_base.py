#Standard imports
import os
import inspect

#Non-standard imports
import catmap
from catmap import ReactionModelWrapper
from catmap.model import ReactionModel

from ase.atoms import string2symbols

class ParserBase(ReactionModelWrapper):
    def __init__(self,reaction_model=ReactionModel()):
        """Class for `parsing' information from raw data 
        (databases, spreadsheets, text files, trajectories, etc.) into a 
        structure which is useful to the microkinetic model. This class acts 
        as a base class to be inherited by other parser classes, but it is 
        not functional on its own.

        input_file: defines the file path or object to get data from

        A functional derived parser class must also contain the methods:

        parse(input_file): a function to parse the input_file file/object and 
        return properly formatted data. The parse function should save all 
        necessary attributes to the Parser class. After parsing the parent 
        microkinetic model class will update itself from the Parser attributes.

        """
        self._rxm = reaction_model
        self._required = {} #No user-defined attributes are required.

    def _baseparse(self):
        #Make dictionary of useful information about species in model
        if not self.species_definitions:
            self.species_definitions = {}
        for species in (self.gas_names+self.adsorbate_names+
                self.transition_state_names):

            ads_info = {}
            if '_' in species:
                name,site = species.rsplit('_',1)
            else:
                name =  species
                site = self._default_site
            ads_info['name'] = name
            ads_info['site'] = site
            if species in self.gas_names:
                ads_info['type'] = 'gas'
                ads_info['n_sites'] = 0
            elif species in self.adsorbate_names:
                ads_info['type'] = 'adsorbate'
                ads_info['n_sites'] = 1
            elif species in self.transition_state_names:
                ads_info['type'] = 'transition_state'
                ads_info['n_sites'] = 1
            else:
                ads_info['type'] = 'unknown'

            composition = {}
            try:
                symbs = string2symbols(name.replace('-',''))
                for a in set(symbs):
                    composition[a] = symbs.count(a)
            except ValueError:
                pass
            ads_info['composition'] = composition
            if species in self.species_definitions:
                ads_info.update(self.species_definitions[species])
            if not ads_info['composition']:
                raise ValueError('Could not determine composition for '+species)
            self.species_definitions[species] = ads_info

        for species in self.species_definitions.keys(): #set site definitions
            site = self.species_definitions[species].get('site',None)
            if site:
                ads_info = {}
                ads_info['type'] = 'site'
                ads_info['site'] = site
                ads_info['formation_energy'] = 0
                if site not in self._gas_sites:
                    ads_info['n_sites'] = 1
                else:
                    ads_info['n_sites'] = 0
                    ads_info['site_names'] = ['gas']
                    ads_info['total'] = 0
                ads_info['composition'] = {}
                if self.site_definitions: #Deprecate later...
                    site_names = self.site_definitions[site]
                    if isinstance(site_names,basestring):
                        site_names = [site_names]
                    ads_info['site_names'] = site_names
                if self.site_totals: #Deprecate later...
                    ads_info['total'] = self.site_totals[site]
                if not self.adsorbate_interactions:
                    ads_info['interaction_threshold'] = 0
                if site in self.species_definitions:
                    ads_info.update(self.species_definitions[site])
                self.species_definitions[site] = self.species_definitions['*_'+site] \
                = ads_info

        if not self.atomic_reservoir_list:
            #Make list of valid reference sets for e.g. boltzmann coverages
            cart_product = []
            all_atoms = []
            composition_dict = {}
            dummy_dict = {}
            for sp in self.gas_names:
                composition_dict[sp] = self.species_definitions[sp]['composition']
                dummy_dict[sp] = 0 #dummy dict of energies
                for key in composition_dict[sp].keys():
                    if key not in all_atoms:
                        all_atoms.append(key)
            for key in all_atoms:
                possibles = []
                for sp in self.gas_names:
                    if composition_dict[sp].get(key,None):
                        possibles.append(sp)
                cart_product.append(possibles)
            
            ref_sets = []
            for prod in catmap.functions.cartesian_product(*cart_product):
                refdict = {}
                for ai,pi in zip(all_atoms,prod):
                    refdict[ai] = pi

                if (sorted(list(refdict.values())) == 
                   sorted(list(set(refdict.values()))) and 
                   sorted(list(refdict.values())) not in 
                   [sorted(list(rs.values())) for rs in ref_sets]):
                    if refdict and dummy_dict and composition_dict:
                        try:
                            self.convert_formation_energies(dummy_dict,
                                    refdict,composition_dict)
                            ref_sets.append(refdict)
                        except ValueError:
                            pass
            self.atomic_reservoir_list = ref_sets
