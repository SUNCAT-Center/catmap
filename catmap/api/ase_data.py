"""
Description:
    Module to convert ase.db files into a catmap input file

    <fname>.db: db file
        Contains all adsorbate states on
        surfaces and clean slabs for reference.
        Mandatory key value pairs:
        --------------------------
            "name" : str
                Value that identifies the catalyst composition.
            "species" : str
                adsorbate chemical formula.
                '' should be the value for clean slabs.
                '-' should be inserted between seperate fragments.
            "energy" or "epot" : float
                potential energy from DFT.

        Recommended key value pairs:
        ----------------------------
            "site" : str
                name of adsorption site.
            "phase" : str
                Value that identifies the catalyst phase.
            "facet" : str
                Name of the facet, preferably in hkl notation, e.g. '(100)'.
                Note that integer strings like '100' are not accepted.
            "surf_lattice" : str
                Name of the surface lattice geometry.
                E.g. HCP(001) and FCC(111) has "hexagonal" surface lattices.
            "layers" : int
                Number of atomic layers in slab.
            "supercell" : str
                Supercell size separated by 'x', e.g. '2x2'
            "n": int
                number of identical adsorbates.

        Recommended keys in "data":
        ---------------------------
            "BEEFvdW_contribs" : list
                32 non-selfconsistent BEEF-vdW energies.
            "frequencies" : list
                vibrational frequencies.

    Optional file dependencies:
    ---------------------------
    <fname>.db : db file
        Stores vibrational frequencies along with atomic structures
        and energies.
"""
import warnings
import os
from uuid import uuid4
import numpy as np
import ase.db
from catmap import string2symbols
from ase.data import covalent_radii, atomic_numbers
from ase.calculators.singlepoint import SinglePointDFTCalculator
from catmap.api.bee import BEEFEnsemble as bee


try:
    from tqdm import tqdm
except (ImportError, ModuleNotFoundError):
    def tqdm(iterable):
        return iterable


class EnergyLandscape(object):
    """Class for converting raw data from ASE db to an energy txt output.
    The class is made for treating atomic structures in the db as points on
    or over a global potential energy surface.
    """
    def __init__(self, beef_size=2000, beef_seed=0):
        """Initialize class."""
        self.bee = bee(size=beef_size, seed=beef_seed)
        self.epot = {}
        self.freq = {}
        self.ens = {}
        self.de_dict = {}
        self.dbid = {}
        self.reference_epot = {}
        self.reference_ens = {}
        self.rxn_paths = {}
        self.formation_energies = {}
        self.std = {}

    def get_molecules(self, fname, selection=[], frequency_db=None):
        """ Method for importing molecules.

        Parameters
        ----------
        fname : str
            path and filename of an ase database file containing molecules.
        selection : list
            ASE database filter strings.
        frequency_db : str
            path and filename of an ASE-db file containing atomic structures
            and vibrational frequencies.
        """
        [mol_epot,
         mol_freq,
         mol_ens,
         mol_dbid] = self._db2mol(fname, selection=selection,
                                  freq_path=frequency_db)
        self.epot.update(mol_epot)
        self.freq.update(mol_freq)
        self.ens.update(mol_ens)
        self.dbid.update(mol_dbid)

    def get_surfaces(self, fname, selection=[], frequency_db=None,
                     site_specific=False):
        """ Method for importing clean slabs and slabs with adsorbates.

        Parameters
        ----------
        fname : str
            path and filename of an ase database file containing slabs.
        selection : list
            ASE database filter strings.
        frequency_db : str
            path and filename of an ASE-db file containing atomic structures
            and vibrational frequencies.
        site_specific : bool
            flag for distinguishing sites or not, when importing adsorbates.
        """
        # Select and import from database file. Return most stable states.
        [surf_epot,
         surf_freq,
         surf_ens,
         surf_dbid] = self._db2surf(fname, selection=selection,
                                    freq_path=frequency_db,
                                    site_specific=site_specific)
        # Store data in dictionaries.
        self.epot.update(surf_epot)
        self.freq.update(surf_freq)
        self.ens.update(surf_ens)
        self.dbid.update(surf_dbid)

    def get_transition_states(self, fname, selection=[], frequency_db=None,
                              site_specific=False):
        """ Method for importing surface transition states.

        Parameters
        ----------
        fname : str
            path and filename of an ase database file
            containing reaction images.
        selection : list
            ASE database filter strings.
        frequency_db : str
            path and filename of an ASE-db file containing atomic structures
            and vibrational frequencies.
        site_specific : bool
            flag for distinguishing sites or not, when importing adsorbates.
        """
        # Select and import images from a database file.
        rxn_paths = self._db2pes(fname, selection=selection,
                                 site_specific=site_specific)
        # Return lowest saddle points.
        self.rxn_paths.update(rxn_paths)
        [surf_epot,
         surf_freq,
         surf_ens,
         surf_dbid] = self.pes2ts(freq_path=frequency_db)
        # Store data in dictionaries.
        self.epot.update(surf_epot)
        self.freq.update(surf_freq)
        self.ens.update(surf_ens)
        self.dbid.update(surf_dbid)

    def calc_formation_energies(self, references, beef=True):
        """ Method for generating formation energies.

        Parameters
        ----------
        references : list of tuples of strings.
            The first item in each tuple must be an atomic symbol, and the
            second item in each tuple must be the self.epot dictionary key
            of a reference gas phase species, <species name>_gas.
        """
        # Get atomic references.
        [self.atomic_e,
         self.atomic_ens] = self._mol2ref(references=references)
        # Get dictionaries with slab references and atomic references.
        [self.reference_epot,
         self.reference_ens] = self._get_refs()
        #
        self.formation_energies = self._get_formation_energies()
        if beef:
            self.de_dict, self.std = self._get_BEEstd()

    def correction(self, key, correction):
        """Apply energy correction to a formation energy.

        Parameters
        ----------
        key : str
            Key from self.formation_energies
        correction : float
            Energy correction to be added.
        """
        self.formation_energies[key] += correction

    def db_attach_reference_id(self, slab_db, ads_db, overwrite=True):
        slab_dict = self._slabs()
        c_ads = ase.db.connect(ads_db)
        for key in slab_dict:
            slab_id = int(slab_dict[key]['id'])
            for ads_id in slab_dict[key]['ads_ids']:
                if overwrite:
                    c_ads.update(ads_id, slab_id=slab_id)
                elif 'slab_id' not in c_ads.get(ads_id):
                    c_ads.update(ads_id, slab_id=slab_id)

    def _slabs(self):
        """Return a dictionary constaining keys of slabs and dictionaries with
        associated adsorbate keys.

        Parameters
        ----------

        """
        ads_slab_pairs = {}
        missing_slab = []
        for key in list(self.epot):
            if 'gas' in key:
                continue
            n, species, cat, pha, lattice, fac, cell, site = key.split('_')
            site_name = '_'.join(['0_', cat, pha, lattice, fac, cell,
                                  'slab'])
            if 'slab' not in key:
                ads = species + '_' + site
                dbid = self.dbid[key]
                try:
                    slab_id = self.dbid[site_name]
                except KeyError:
                    missing_slab.append(str(self.dbid[key]))
                    continue
                if site_name not in ads_slab_pairs:
                    ads_slab_pairs.update({site_name: {'id': slab_id,
                                                       'species': [ads],
                                                       'ads_ids': [dbid]}})
                else:
                    ads_slab_pairs[site_name]['species'].append(ads)
                    ads_slab_pairs[site_name]['ads_ids'].append(dbid)
        if len(missing_slab) > 0:
            print('Missing slabs: ' + ','.join(missing_slab))
        return ads_slab_pairs

    def _db2mol(self, fname, selection=[], freq_path=None):
        """ Returns four dictionaries containing:
            ab initio energies,
            frequecies,
            non-selfconsistent BEEF perturbations,
            database ids.

        Parameters
        ----------
        fname : str
            path and filename of ase db file.
        selection : list
            ASE database filter strings.
        freq_path : str
            path and filename of an ASE-db file containing atomic structures
            and vibrational frequencies.

        File dependencies:
        ------------------
        fname : ase-db file
            Contains molecular reference states.

            Mandatory key value pairs:
            --------------------------
                    "energy" or "epot" : float
                        DFT calculated potential energy.
            Optional key value pairs:
            --------------------------
                    "data.BEEFens" : list
                        32 non-selfconsistent BEEF-vdW energies.
        """
        # Connect to a database.
        cmol = ase.db.connect(fname)
        # Select data using search filters.
        smol = cmol.select(selection)
        # Connect to a database with frequencies.
        if freq_path is not None:
            c_freq = ase.db.connect(freq_path)
        abinitio_energies = {}
        freq_dict = {}
        dbids = {}
        ens_dict = {}
        # Iterate over molecules.
        for d in smol:
            if 'energy' in d:
                abinitio_energy = float(d.energy)
            else:
                abinitio_energy = float(d.epot)
            species_name = str(d.formula)
            # Attempt to retrieve the 32 BEEF perturbations.
            try:
                contribs = d.data.BEEFvdW_contribs
                ens = self.bee.get_ensemble_perturbations(contribs)
            except AttributeError:
                ens = 0
            # Store the most stable state of each molecule.
            if species_name + '_gas' not in abinitio_energies:
                abinitio_energies[species_name+'_gas'] = abinitio_energy
                dbids[species_name + '_gas'] = int(d.id)
                ens_dict[species_name + '_gas'] = ens
                if freq_path is not None:
                    try:
                        d_freq = c_freq.get(['formula=' + species_name])
                        frequencies = d_freq.data.frequencies
                        freq_dict.update({species_name + '_gas': frequencies})
                    except KeyError:
                        continue
            elif abinitio_energies[species_name+'_gas'] > abinitio_energy:
                abinitio_energies[species_name+'_gas'] = abinitio_energy
                dbids[species_name + '_gas'] = int(d.id)
                ens_dict[species_name + '_gas'] = ens
        return abinitio_energies, freq_dict, ens_dict, dbids

    def _db2surf(self, fname, selection=[], freq_path=None,
                 site_specific=False):
        """ Returns four dictionaries containing:
            ab initio energies,
            frequecies,
            non-selfconsistent BEEF perturbations,
            database ids.

            File dependencies
            -----------------
            fname : ase.db file.

            Parameters
            ----------
            fname : string
                path/filname.
            selection : list
                Optional ASE-db filter strings.
            site_specific : boolean
                If True: Dinstinguish sites using the site key value pair, and
                stores a the potential energy of adsorbates on each site.
                Else: Use the minimum ab initio energy, disregarding the site.
        """
        # Connect to a database with clean surfaces and/or adsorbates.
        csurf = ase.db.connect(fname)
        # Connect to a database with frequencies.
        if freq_path is not None:
            c_freq = ase.db.connect(freq_path)
        ssurf = csurf.select(selection)
        abinitio_energies = {}
        freq_dict = {}
        dbids = {}
        ens_dict = {}
        # Loop over states.
        for d in ssurf:
            [n, species, name, phase, surf_lattice, facet,
             cell] = self._get_adsorbate_fields(d)
            # Skip any transition states.
            if '-' in species:
                continue
            if 'ads' in d:
                ads = str(d.ads)
            else:
                ads = species
            if 'energy' in d:
                abinitio_energy = float(d.energy)
            else:
                abinitio_energy = float(d.epot)
            if species == '' or ('ads' in d and
                                 (ads == 'slab' or ads == 'clean')):
                species = ''
                ads = 'slab'
                site = 'slab'
                n = 0
            elif 'site' in d and site_specific:
                site = str(d.site)
            else:
                site = 'site'
            # Make key unique to a physical state of a site.
            key = '_'.join([str(n), species, name, phase, surf_lattice,
                            facet, cell, site])
            freq_key = key
            # Attempt to import the 32 BEEF perturbations.
            try:
                contribs = d.data.BEEFvdW_contribs
                ens = self.bee.get_ensemble_perturbations(contribs)
            except AttributeError:
                ens = 0
            if key not in abinitio_energies:
                abinitio_energies[key] = abinitio_energy
                dbids[key] = int(d.id)
                ens_dict[key] = ens
                if species != '' and ads != 'slab' and freq_path is not None:
                    try:
                        freqsearch = ['species=' + species, 'name=' + name]
                        if site_specific is True:
                            freqsearch.append('site=' + site)
                        d_freq = c_freq.get(freqsearch)
                        frequencies = d_freq.data.frequencies
                        freq_dict.update({freq_key: frequencies})
                    except KeyError:
                        continue
            elif abinitio_energies[key] > abinitio_energy:
                abinitio_energies[key] = abinitio_energy
                dbids[key] = int(d.id)
                ens_dict[key] = ens
        return abinitio_energies, freq_dict, ens_dict, dbids

    def _get_adsorbate_fields(self, d):
        """Return a set of fields characterizing an adsorbate state.

        Parameters
        ----------
        d : dictionary
            ASE database row.

        Returns
        ----------
        n : str
            Number of adsorbates.
        species : str
            Species name. Must be chemical formula.
        name : str
            Name of catalyst.
        phase : str
            Crystal structure.
        surf_lattice : str
            Name of surface structure.
        facet : str
            Facet <hkl>
        cell : str
            Slab size <XxYxL>, where L denotes layers.
        """
        if 'species' in d:
            species = str(d.species)
        else:
            species = ''
        name = str(d.name)
        if 'supercell' in d:
            cell = str(d.supercell)
        else:
            cell = 'XxY'
        if 'layers' in d:
            cell += 'x' + str(d.layers)
        if 'crystal' in d:
            phase = str(d.crystal)
        elif 'phase' in d:
            phase = str(d.phase)
        elif 'spacegroup' in d:
            phase = str(d.spacegroup)
        else:
            phase = ''
        if 'surf_lattice' in d:
            surf_lattice = str(d.surf_lattice)
        else:
            surf_lattice = ''
        if 'facet' in d:
            facet = str(d.facet)
        else:
            facet = 'facet'
        if 'n' in d:
            n = int(d.n)
        elif species == '':
            n = 0
        else:
            n = 1
        return n, species, name, phase, surf_lattice, facet, cell

    def _db2pes(self, fname, selection=[], site_specific=False):
        """Returns a dictionary containing potential energy surfaces and
        meta data.

        Dependencies
        ----------
        fname : ase.db file.

        Parameters
        ----------
        fname : str
            path/filname.
        selection : list of strings.
            Optional ase.db selection strings.
        site_specific : boolean
            If True: Dinstinguish sites using the site key value pair, and
            stores a the potential energy of adsorbates on each site.
            Else: Use the minimum ab initio energy, disregarding the site.
        """
        c = ase.db.connect(fname)
        s = c.select(selection)
        rxn_paths = {}
        # Loop over states from ase .db
        for d in s:
            # Store variables identifying the atomic states.
            species = str(d.species)
            # - identifies transition states.
            if '-' not in species:
                continue
            # Most fiels are optional.
            if 'supercell' in d:
                cell = str(d.supercell)
            else:
                cell = 'XxY'
            if 'layers' in d:
                cell += 'x' + str(d.layers)
            if 'phase' in d:
                phase = str(d.phase)
            else:
                phase = ''
            if 'surf_lattice' in d:
                surf_lattice = str(d.surf_lattice)
            else:
                surf_lattice = ''
            if 'facet' in d:
                facet = str(d.facet)
            else:
                str(d.facet)
            surf = '_'.join([str(d.name), phase, surf_lattice, facet, cell])
            # 'energy' only exist if a calculator is attached to the atoms.
            if 'energy' in d:
                abinitio_energy = float(d.energy)
            else:
                abinitio_energy = float(d.epot)
            dbid = int(d.id)
            # Try to import non-selfconsistent BEEF contributions.
            try:
                BEEFvdW_contribs = d.data.BEEFvdW_contribs
                ens = self.bee.get_ensemble_perturbations(BEEFvdW_contribs)
            except AttributeError:
                ens = 0  # np.zeros(self.bee.size)
            if 'path_id' in d:
                rxn_id = str(d.path_id)
            else:
                rxn_id = uuid4().hex
            if 'distance' in d:
                distance = float(d.distance)
            else:
                distance = np.nan
            if 'step' in d:
                step = int(d.step)
            else:
                step = np.nan
            if rxn_id in rxn_paths:
                rxn_paths[rxn_id]['pes'].append(abinitio_energy)
                rxn_paths[rxn_id]['dbids'].append(dbid)
                rxn_paths[rxn_id]['ens'].append(ens)
                rxn_paths[rxn_id]['distance'].append(distance)
                # except AttributeError:
                #    d0 = c.get(rxn_paths[rxn_id]['dbids'][0])
                #    atoms0 = c.get_atoms(rxn_paths[rxn_id]['dbids'][0])
                #    species0 = str(d0.species)
                #    g1, g2 = species0.split('-')
                #    cons = atoms0.constraints
                #    fbl = cons[-1].todict()['kwargs']['pairs'][0]
                #    atom1 = int(fbl[0])
                #    atom2 = int(fbl[1])
                #    assert (atoms0[atom2].symbol == g1[0] and
                #            atoms0[atom1].symbol == g2[0]) or \
                #        (atoms0[atom1].symbol == g1[0] and
                #         atoms0[atom2].symbol == g2[0])
                #    atoms = c.get_atoms(dbid)
                #    dist = atoms.get_distance(atom1, atom2, mic=True)
                #    print(dbid, dist, atom1, atom2)
                #    c.update(dbid, distance=dist)
                if 'image' in d:
                    rxn_paths[rxn_id]['images'].append(int(d.image))
                else:
                    rxn_paths[rxn_id]['images'].append(int(d.step))
            else:
                if 'site' in d and site_specific:
                    site = str(d.site)
                else:
                    site = 'site'
                rxn_paths[rxn_id] = {'surface_name': surf,
                                     'species': species,
                                     'pes': [abinitio_energy],
                                     'dbids': [dbid],
                                     'ens': [ens],
                                     'site': site,
                                     'images': [step],
                                     'distance': [distance]}
        return rxn_paths

    def _mol2ref(self, references):
        """ Returns two dictionaries containing:
            abinitio energy references for atoms
            ensemble non-selfconsistent perturbations for atoms.
        """
        atomic_e = {}
        atomic_ens = {}
        for t in references:
            key = t[0]
            species = t[1]
            atomic_e[key] = self.epot[species]
            if species in self.ens:
                atomic_ens[key] = np.array(self.ens[species])
            composition = string2symbols(species.split('_')[0])
            unique_element, count = np.unique(composition, return_counts=True)
            n = None
            for i, symbol in enumerate(unique_element):
                if symbol == key:
                    n = count[i]
                else:
                    atomic_e[key] -= atomic_e[symbol] * float(count[i])
                    if key in atomic_ens and symbol in atomic_ens:
                        atomic_ens[key] = atomic_ens[key] - \
                            np.array(atomic_ens[symbol]) * float(count[i])
            atomic_e[key] = atomic_e[key] / float(n)
            if key in atomic_ens:
                atomic_ens[key] = atomic_ens[key] / float(n)
        return atomic_e, atomic_ens

    def _get_refs(self):
        """ Returns dictionaries with referece energies of slabs and
        single atoms.

        Parameters
        ----------
        references : list of tuples of strings.
            The first item in each tuple must be an atomic symbol, and the
            second item in each tuple must be the self.epot dictionary key
            of a reference gas phase species, <species name>_gas.
        """
        ref_dict = self.atomic_e
        ref_ens = self.atomic_ens
        for key in self.epot.keys():
            if 'slab' in key:
                ref_dict[key] = self.epot[key]
                if key in self.ens:
                    ref_ens[key] = self.ens[key]
        return ref_dict, ref_ens

    def _get_formation_energies(self):
        """ Returns a dictionary with formation energies of adsorbates.

        Dependencies
        ----------
        self : db2catmap object
            self.epot : dictionary
                Each key is named in the format:
                n_species_name_phase_facet_supercell{x}layers_site,
                and contains the potential energy of a slab an adsorbate or
                a molecule.
            self.reference_epot : dictionary
                Each key is either an atomic symbol and contains the reference
                potential energy of that atom,
                or the key is named in the format: _name_phase_facet_slab and
                it contains the reference potential energy of the slab.
        """
        formation_energies = {}
        missing_slab = []
        for key in list(self.epot):
            E0 = 0
            if 'gas' in key:
                species, site_name = key.split('_')  # Split key into name/site
            else:
                n, species, cat, pha, lattice, fac, cell, site = key.split('_')
                site_name = '_'.join(['0_', cat, pha, lattice, fac, cell,
                                      'slab'])
                try:
                    E0 -= self.reference_epot[site_name]
                except KeyError:
                    missing_slab.append(str(self.dbid[key]))
                    continue
            if 'slab' not in key:
                composition = string2symbols(species.replace('-', ''))
                E0 += self.epot[key]
                for atom in composition:
                    E0 -= self.reference_epot[atom]
                formation_energies[key] = E0
                if abs(E0) / len(composition) > 5.:
                    warnings.warn('Large formation energy: ' +
                                  str(E0 / len(composition)) +
                                  ' eV per atom. ' + str(self.dbid[key]) +
                                  ': ' + key)
        if len(missing_slab) > 0:
            print('Missing slabs: ' + ','.join(missing_slab))
        return formation_energies

    def _get_BEEstd(self):
        """ Returns a dictionary with BEEF ensembles and one with
        BEEF standard deviations on formation energies.

        """
        de_dict = {}
        std_dict = {}
        for key in self.ens:
            de = np.zeros(self.bee.size)
            if 'gas' in key:
                species, site_name = key.split('_')  # Split key into name/site
            else:
                n, species, cat, pha, lattice, fac, cell, site = key.split('_')
                site_name = '_'.join(['0_', cat, pha, lattice, fac, cell,
                                      'slab'])
                try:
                    de -= self.reference_ens[site_name]
                except KeyError:
                    continue
            if 'slab' not in key:
                composition = string2symbols(species.replace('-', ''))
                de += self.ens[key]
                for atom in composition:
                    de -= self.reference_ens[atom]
                de_dict[key] = de
                sigma = np.std(de)
                std_dict[key] = sigma
                if sigma / len(composition) > 0.5:
                    msg = "Large BEEF 1 sigma: " + \
                        str(sigma / len(composition)) + " eV per atom. " + \
                        str(self.dbid[key]) + ": " + key
                    warnings.warn(msg)
        return de_dict, std_dict

    def get_ellipses(self, ads_x, ads_y,
                     site_x=None, site_y=None):
        """ Returns three dictionaries, width, height and angle with the
        parameters for plotting the covariance ellipses showing the
        +/- 1 sigma confidence intervals.

        Parameters
        ----------
        de_dict : dict
            contains beef perturbations of adsorbates on slabs,
            where keys are named: adsorbate_name_phase_facet
        ref_de : dict
            contains beef perturbations of references,
            where keys are refernce elements, e.g: 'C','H',
            and also slabs references.
        ads_x : str
            adsorbate first dimension
        ads_y : str
            adsorbate second dimension
        """
        widths = {}
        heights = {}
        angles = {}
        # Loop over reference surfaces in ref_de.
        for slab in self.reference_epot.keys():
            # Ignore gas species.
            if 'slab' not in slab:
                continue
            key_x = self._create_state(slab, ads_x, site_x)
            if isinstance(key_x, list):
                de_x = np.zeros(self.bee.size)
                continue_outer = False
                for k_x in key_x:
                    if k_x in self.de_dict:
                        de_x += self.de_dict[k_x]
                    else:
                        continue_outer = True
                        break
                if continue_outer is True:
                    continue
            elif key_x in self.de_dict:
                de_x = self.de_dict[key_x]
            else:
                continue
            key_y = self._create_state(slab, ads_y, site_y)
            if key_y in self.de_dict:
                de_y = self.de_dict[key_y]
            else:
                continue
            if np.isclose(de_x, 0).all() or np.isclose(de_y, 0).all():
                continue
            width, height, angle = self.bee.get_ellipse(de_x, de_y)
            widths[slab] = width
            heights[slab] = height
            angles[slab] = angle
        self.width = widths
        self.height = heights
        self.angle = angles
        return widths, heights, angles

    def pes2ts(self, freq_path=None, rtol=1.1):
        """ Returns dictionaries containing transition state energies.

            Parameters
            ----------
            freq_path : str
                path/folder where frequency database is located.
            rtol : float
                relative tolerance of the threshold distance, where fixed bond
                lenght calculations are considered complete.
        """
        if freq_path is not None:
            c_freq = ase.db.connect(freq_path)
        abinitio_energies = {}
        freq_dict = {}
        dbids = {}
        ens_dict = {}
        calculate = []
        for rxn_id in self.rxn_paths:
            warn = False
            species = self.rxn_paths[rxn_id]['species']
            m = self.rxn_paths[rxn_id]['surface_name']
            site = self.rxn_paths[rxn_id]['site']
            key = '1_' + species + '_' + m + '_' + site
            images = self.rxn_paths[rxn_id]['images']
            pes = np.array(self.rxn_paths[rxn_id]['pes'])
            if len(images) > 1:
                s = np.argsort(images)
                # Look for local minima and maxima.
                localmins = np.where(np.r_[True, pes[s][1:] < pes[s][:-1]] &
                                     np.r_[pes[s][:-1] < pes[s][1:], True])[0]
                localmaxs = np.where(np.r_[True, pes[s][1:] > pes[s][:-1]] &
                                     np.r_[pes[s][:-1] > pes[s][1:], True])[0]
                # Measure path roughness
                differences = np.diff(pes[s])
                roughness = np.std(differences)
                # For fixed bond length (drag) calculations
                g1, g2 = species.split('-')
                dbond = covalent_radii[atomic_numbers[g1[0]]]
                if len(g2) > 0:
                    dbond += covalent_radii[atomic_numbers[g2[0]]]
                else:
                    # Assume the bond is with the surface.
                    try:
                        dbond += covalent_radii[atomic_numbers[
                            m.split('_')[0]]]
                    except KeyError:
                        print("Bond not defined.")
                if len(np.unique(images)) != len(images):
                    warn = True
                    print('non unique image number!')
                    print('Warning!', species, m, roughness,
                          len(localmaxs), len(localmins), images)
                    continue
                if (len(localmaxs) > 1 or
                    len(localmins) > 2 or
                   len(localmins) == 1):
                    warn = True
                try:
                    shortest = np.nanmin(self.rxn_paths[rxn_id]['distance'])
                    if shortest > dbond * rtol:
                        warn = True
                        s_last = np.argmax(self.rxn_paths[rxn_id]['images'])
                        calculate.append(
                            self.rxn_paths[rxn_id]['dbids'][s_last])
                        continue
                except KeyError:
                    print('Distances missing in reaction path.')
                if warn:
                    warnings.warn("Warning! " + species + "*" + m + " " +
                                  str(round(dbond * rtol, 3)) + " AA. " +
                                  str(round(shortest, 3)) + " AA. " +
                                  str(roughness) + " eV. " +
                                  str(len(localmaxs)) + " local maxima. " +
                                  str(len(localmins)) + " local minima. " +
                                  str(len(images)) + " images.")
            tst = np.argmax(pes)
            if key not in abinitio_energies:
                abinitio_energies[key] = pes[tst]
                dbids[key] = self.rxn_paths[rxn_id]['dbids'][tst]
                ens_dict[key] = self.rxn_paths[rxn_id]['ens'][tst]
                if freq_path is not None:
                    try:
                        d_freq = c_freq.get('path_id=' + rxn_id)
                        frequencies = d_freq.data.frequencies
                        freq_dict.update({key: frequencies})
                    except KeyError:
                        continue
            elif abinitio_energies[key] > pes[tst]:
                abinitio_energies[key] = pes[tst]
                dbids[key] = self.rxn_paths[rxn_id]['dbids'][tst]
                ens_dict[key] = self.rxn_paths[rxn_id]['ens'][tst]
        if len(calculate) > 0:
            incomplete = ','.join([str(int(a)) for a in np.unique(calculate)])
            print('Incomplete:', incomplete)
        return abinitio_energies, freq_dict, ens_dict, dbids

    def scaling_analysis(self, x, y, lattice=None, site_x=None, site_y=None):
        """Returns the scaling relation information between the species
        x and y on the surface geometry 'lattice' and site.

        Parameters
        ----------
        x : str or list
            species x
        y : str
            species y
        lattice : str or None
            surface lattice of y
        site_x='site' : str
            Site of x. If _db2surf was run with site_specific=False, use the
            default value, None or 'site'.
        site_y='site' : str
            Site of y. If _db2surf was run with site_specific=False, use the
            default value, None or 'site'.

        Returns
        ----------
        slope : float
            Scaling relation slope.
        intercept : float
            Scaling relation intercept.
        X : list
            List of formation energies.
        Y : list
            List of formation energies.
        labels : list
            List of keys referring to slabs.
        """
        X = []
        Y = []
        labels = []
        # Loop over slabs.
        for slab in self.reference_epot.keys():
            if 'slab' not in slab:
                continue
            # Filter by surface lattice.
            if lattice is not None and lattice not in slab:
                continue
            # Get formation energy if it is stored.
            key_y = self._create_state(slab, y, site_y)
            if key_y not in self.formation_energies:
                continue
            key_x = self._create_state(slab, x, site_x)
            # If the x species is a list, sum them.
            if isinstance(key_x, list):
                DeltaE_x = 0
                continue_outer = False
                for k_x in key_x:
                    if k_x in self.formation_energies:
                        DeltaE_x += self.formation_energies[k_x]
                    else:
                        continue_outer = True
                        break
                if continue_outer:
                    continue
                else:
                    X.append(DeltaE_x)
                    Y.append(self.formation_energies[key_y])
            else:
                if key_x in self.formation_energies:
                    X.append(self.formation_energies[key_x])
                    Y.append(self.formation_energies[key_y])
                else:
                    continue
            # Store the list of slabs.
            labels.append(slab)
        # Get the scaling relation.
        slope, intercept = np.polyfit(X, Y, deg=1)
        return slope, intercept, X, Y, labels

    def _create_state(self, slab, species, site=None):
        """Return a key for a hypothetical adsorbate state. This is useful
        for reading or filling in formation energies of states that are related
        by scaling relations.

        Parameters
        ----------
        slab : str
            key of the slab
        species : str
            adsorbate chemical formula
        site : str
            optional site.
        """
        if isinstance(species, list):
            key = []
            for x in range(len(species)):
                fields = slab.split('_')
                fields[0] = '1'
                fields[1] = species[x]
                if site is None:
                    fields[-1] = 'site'
                else:
                    fields[-1] = site[x]
                key.append('_'.join(fields))
        else:
            fields = slab.split('_')
            fields[0] = '1'
            fields[1] = species
            if site is None:
                fields[-1] = 'site'
            else:
                fields[-1] = site
            key = '_'.join(fields)
        return key

    def insert_interpolated_states(self, x, y, lattice=None, site_y=None,
                                   slope=None, intercept=None):
        """ Update the formation_energy dictionary with interpolated values.
        This is intended for use with thermodynamic scalers only.

        Parameters
        ----------
        x : list of strings
            keys from self.formation_energy
        y : str
            species
        site_y : str
            site of species y
        slope : float or None
        intercept : float or None
        """
        if slope is None or intercept is None:
            raise NotImplementedError("Call scaling_analysis.")
            [slope,
             intercept,
             X,
             Y,
             labels] = self.scaling_analysis(x, y, lattice=lattice,
                                             site_x=None, site_y=None)
        for key_x in x:
            X = self.formation_energies[key_x]
            Y = slope * X + intercept
            fields = key_x.split('_')
            fields[1] = y
            if site_y is not None:
                fields[-1] = site_y
            key_y = '_'.join(fields)
            self.formation_energies.update({key_y: Y})

    def _insert_state(self, key, descriptor, slopes, intercept):
        """ Update the formation_energy dictionary with interpolated values.
        This is intended for use with thermodynamic scalers only.

        Parameters
        ----------
        key : str
            exact key of the species to be inserted self.formation_energy
        descriptor : str or list
            species names of descriptors.
        slope : float or None
        intercept : float or None
        """
        interpolated = intercept
        fields = key.split('_')
        for a in range(len(descriptor)):
            fields[1] = descriptor[a]
            fields[-1] = 'site'
            energy_key = '_'.join(fields)
            interpolated += slopes[a] * self.formation_energies[energy_key]
        state = {key: interpolated}
        self.formation_energies.update(state)

    def _insert_frequencies(self, key, frequencies):
        """ Update the formation_energy dictionary with interpolated values.
        This is intended for use with thermodynamic scalers only.

        Parameters
        ----------
        key : str
            exact key of the species to be inserted self.freq
        frequencies : list
            list of frequencies.
        slope : list
        """
        state = {key: frequencies}
        self.freq.update(state)

    def insert_rscaled_states(self, x, y, site_y=None,
                              slope=None, intercept=None):
        raise NotImplementedError("Coming in 2018.")

    def make_ensemble_input_files(self, prefix, suffix, site_specific=False):
        """ Save catmap input files for ensemble models.
        It is advisable to use a smaller than default beef_size.
        """
        if self.beef_size >= 2000:
            warnings.warn("It is advisable to use a smaller than default " +
                          "beef_size for BEEF ensemble propagation through" +
                          "the micro-kinetic model.")
        raise NotImplementedError("Coming in 2018.")
        # Create a header.
        # headerlist = ['surface_name', 'phase', 'site_name',
        #              'species_name', 'formation_energy',
        #              'frequencies', 'reference', 'coverage', 'std']
        # header = '\t'.join(headerlist)

    def db_attach_formation_energy(self, fname, key_name, overwrite=True):
        """ Update a database file to append formation energies.

        Parameters
        ----------
            fname : str
                path and filename of ase database.
        """
        c = ase.db.connect(fname)
        for key in tqdm(list(self.formation_energies)):
            if 'gas' not in str(key):
                if overwrite:
                    kvp = {key_name: float(self.formation_energies[key])}
                    c.update(int(self.dbid[key]), **kvp)
                elif key_name not in c.get(int(self.dbid[key])):
                    kvp = {key_name: float(self.formation_energies[key])}
                    c.update(int(self.dbid[key]), **kvp)

    def make_input_file(self, file_name, site_specific='facet',
                        catalyst_specific=False, covariance=None,
                        reference=None):
        """ Saves the catmap input file.

        Parameters
        ----------
        file_name : string
            path and name of the output file.
        site_specific : boolean or string
            Decides what to export to the site key.
            True exports the site field from the db.
            False exports the lattice field from the db.
            'facet' exports the facet field from the db.
            str : another string is treated as False, except if that string is
            found in the site field.
        covariance : tuple
            Must contains two strings, which are species names, between which
            BEEF covariance ellipses will be stored in the
            width, heigh and angle columns of the catmap data file.
        """
        # Create a header.
        headerlist = ['surface_name', 'phase', 'site_name',
                      'species_name', 'formation_energy',
                      'frequencies', 'reference', 'coverage', 'std']
        if covariance is not None:
            headerlist += ['width', 'height', 'angle', 'covariance']
        header = '\t'.join(headerlist)
        # List of lines in the output.
        lines = []
        for key in self.formation_energies.keys():  # Iterate through keys
            if key in self.dbid and reference is None:
                ref = self.dbid[key]
            else:
                ref = reference
            E = round(self.formation_energies[key], 4)
            if 'gas' in key:
                name, site = key.split('_')  # Split key into name/site
                try:
                    frequency = self.freq[key]
                except KeyError:
                    frequency = []
                try:
                    std = round(self.std[key], 4)
                except KeyError:
                    std = np.NaN
            else:
                n, name, cat, pha, lattice, facet, cell, site = key.split('_')
                if catalyst_specific and not catalyst_specific == cat:
                    continue
            if 'slab' not in key:  # Do not include empty site energy (0)
                try:
                    freq_key = key
                    frequency = self.freq[freq_key]
                except KeyError:
                    frequency = []
                try:
                    std = round(self.std[key], 4)
                except KeyError:
                    std = np.NaN
                if site == 'gas':
                    surface = None
                    phase = ''
                    lattice = ''
                    site_name = 'gas'
                    coverage = 0.
                else:
                    surface = cat
                    phase = pha
                    try:
                        coverage = round(float(n) /
                                         (float(cell[0]) * float(cell[2])), 3)
                    except ValueError:
                        coverage = 0.
                    if site_specific is True:
                        site_name = lattice + '_' + site
                    elif site_specific is False:
                        site_name = lattice
                    elif site_specific == 'facet':
                        site_name = facet
                    else:
                        if site == site_specific:
                            site_name = site
                        else:
                            site_name = lattice
                E = float(E)
                outline = [surface, phase, site_name, name,
                           E, list(frequency), ref, coverage, std]
                if covariance is not None:
                    if covariance[1] == name:
                        slab = '_'.join(['0_', cat, pha, lattice, facet, cell,
                                         'slab'])
                        if slab in self.width:
                            width = round(self.width[slab], 4)
                            height = round(self.height[slab], 4)
                            angle = round(self.angle[slab], 4)
                            ellipse_ref = covariance[0]
                        else:
                            width = ''
                            height = ''
                            angle = ''
                            ellipse_ref = ''
                    else:
                        width = ''
                        height = ''
                        angle = ''
                        ellipse_ref = ''
                    outline += [width, height, angle, ellipse_ref]
                line = '\t'.join([str(w) for w in outline])
                lines.append(line)
        # The file is easier to read if sorted (optional).
        lines.sort()
        # Add header to top.
        lines = [header] + lines
        # Join the lines with a line break.
        input_file = '\n'.join(lines)
        # Open the file name in write mode.
        input = open(file_name, 'w')
        # Write the text.
        input.write(input_file)
        # Close the file.
        input.close()
        print("Formation energies exported to " + file_name)

    def make_nested_folders(self, project, reactions, surfaces=None,
                            site='site', mol_db=None,
                            slab_db=None, ads_db=None, ts_db=None,
                            publication='', url='', xc='xc', code='code'):
        """Saves a nested directory structure.
        The folder structure for catalysis-hub.org should be
        <project>
            <code>
                <xc>
                    <catalyst>
                        <facet>
                            <reaction>@<site>
                                <species>.traj or
                                <species>_<slab>.traj or
                                'TS'.traj
                        <catalyst>_<phase>_<bulk>
                    <gas>
                        <species>.traj>

        Parameters
        ----------
        project : str
            parent folder name.
        reactions : list
            catmap's rxn_expressions. A list of strings.
        surfaces : list
            List of catalyst names.
        site : str
            Site name.
        mol_db : str
            Path and filename of database containing gas species
        slab_db : str
            Path and filename of database containing slabs
        ads_db : str
            Path and filename of database containing adsorbate/slab structures.
        ts_db : str
            Path and filename of database containing reaction paths.
        publication : str
            Author or publication reference.
        url : str
            url to publication.
        """
        data_folder = project + '/' + code + '/' + xc

        # Create a header
        # spreadsheet = [['chemical_composition', 'facet', 'reactants',
        #                'products', 'reaction_energy',
        #                'beef_standard_deviation',
        #                'activation_energy', 'DFT_code', 'DFT_functional',
        #                'reference', 'url']]

        # Connect to databases.
        if surfaces is None:
            surfaces = [s for s in self.reference_epot.keys() if 'slab' in s]
        if mol_db is None:
            mol_db = self.mol_db
        if ads_db is None:
            ads_db = self.ads_db
        if slab_db is None:
            slab_db = ads_db
        if ts_db is None:
            ts_db = self.ts_db
        c_mol = ase.db.connect(mol_db)
        c_ads = ase.db.connect(ads_db)
        if slab_db == ads_db:
            c_slab = c_ads
        else:
            c_slab = ase.db.connect(slab_db)
        if ts_db is not None:
            c_ts = ase.db.connect(ts_db)

        # Iterate over surfaces
        Nsurf = 0
        Nrxn = 0
        for slabkey in surfaces:
            [n, species, name, phase,
             lattice, facet, cell, slab] = slabkey.split('_')
            catalyst_name = name.replace('/', '-') + '_' + phase
            path_surface = data_folder + '/' + catalyst_name
            facet_name = facet.replace('(', '').replace(')', '') + '_' + cell
            path_facet = path_surface + '/' + facet_name

            # Loop over reaction expressions.
            for i, rxn in enumerate(reactions):
                # Separate left and right side of reaction, and ts.
                states = rxn.replace(' ', '').split('<->')
                if len(states) == 1:
                    states = states[0].split('->')
                    if len(states) == 1:
                        states = states[0].split('<-')
                elif len(states) < 3:
                    states = [states[0]] + states[-1].split('->')
                    if len(states) < 3:
                        states = states[0].split('<-') + states[1:]

                # List individual species.
                rname, reactants = self._state2species(states[0])
                pname, products = self._state2species(states[-1])

                reaction_name = '__'.join([rname, pname])
                path_reaction = path_facet + '/' + reaction_name
                DeltaE = 0.
                de = np.zeros(self.bee.size)
                ea = 0.
                intermediates_exist = True
                totraj = {}
                # Find reactant structures and energies
                for reactant in reactants:
                    species, sitesymbol = reactant.split('_')
                    n, species = self._coefficient_species(species)
                    if sitesymbol == 'g':
                        rkey = species + '_gas'
                        fname = data_folder + '/gas/' + rkey
                    elif species == '*':
                        rkey = slabkey
                        fname = path_facet + '/empty_slab'
                    else:
                        rkey = '_'.join(['1', species.replace('*', ''), name,
                                         phase, lattice, facet, cell, site])
                        fname = path_reaction + '/' + species
                    if rkey not in self.dbid:
                        intermediates_exist = False
                        break
                    totraj.update({rkey:
                                   {'dbid': self.dbid[rkey],
                                    'fname': fname}})
                    if species != '*':
                        DeltaE -= n * self.formation_energies[rkey]
                        de -= n * self.de_dict[rkey]
                if not intermediates_exist:
                    continue
                # Find transition state structures and energies.
                if ts_db is not None:
                    tstates = states[1].split('+')
                    for ts in tstates:
                        if '-' not in ts:
                            continue
                        species, sitesymbol = ts.split('_')
                        tskey = '_'.join(['1', species, name, phase, lattice,
                                          facet, cell, site])
                        if tskey not in self.dbid:
                            continue
                        totraj.update({tskey:
                                       {'dbid': self.dbid[tskey],
                                        'fname': path_reaction + '/TS'}})
                        ea += self.formation_energies[tskey] - DeltaE
                # Find product structures and energies.
                for product in products:
                    species, sitesymbol = product.split('_')
                    if species[0].isdigit():
                        n = int(species[0])
                        species = species[1:]
                    else:
                        n = 1
                    if sitesymbol == 'g':
                        pkey = species + '_gas'
                        fname = data_folder + '/gas/' + pkey
                    elif species == '*':
                        pkey = slabkey
                        fname = path_facet + '/empty_slab'
                    else:
                        pkey = '_'.join(['1', species.replace('*', ''), name,
                                         phase, lattice, facet, cell, site])
                        fname = path_reaction + '/' + species
                    if pkey not in self.dbid:
                        intermediates_exist = False
                        break
                    totraj.update({pkey:
                                   {'dbid': self.dbid[pkey],
                                    'fname': fname}})
                    if species != '*':
                        DeltaE += n * self.formation_energies[pkey]
                        de += n * self.de_dict[pkey]
                # If all states are found for this surface, write.
                if intermediates_exist:
                    # Loop over states to export.
                    for trajkey in totraj.keys():
                        fname = totraj[trajkey]['fname'] + '.traj'
                        # Load the atomic structure from appropriate db.
                        if 'gas' in trajkey:
                            atoms = c_mol.get_atoms(self.dbid[trajkey])
                            d = c_mol.get(self.dbid[trajkey])
                        elif '-' in trajkey.split('_')[1]:
                            atoms = c_ts.get_atoms(self.dbid[trajkey])
                            d = c_ts.get(self.dbid[trajkey])
                        elif 'slab' in trajkey.split('_')[-1]:
                            atoms = c_slab.get_atoms(self.dbid[trajkey])
                            d = c_slab.get(self.dbid[trajkey])
                        else:
                            atoms = c_ads.get_atoms(self.dbid[trajkey])
                            d = c_ads.get(self.dbid[trajkey])
                        if atoms.calc is None:
                            # Require a calculator.
                            calc = SinglePointDFTCalculator(atoms)
                            calc.results['energy'] = float(d.epot)
                            atoms.set_calculator(calc)
                        if 'data' not in atoms.info:
                            atoms.info['data'] = {}
                        if trajkey in self.freq:
                            # Attach vibrational frequencies.
                            atoms.info['data'].update(
                                {'frequencies': self.freq[trajkey]})
                        # Save trajectory file.
                        folder_structure = fname.split('/')
                        for depth in range(1, len(folder_structure)-1):
                            directory = '/'.join(folder_structure[:depth+1])
                            if not os.path.isdir(directory):
                                os.mkdir(directory)
                        atoms.write(fname)
                    # Store rows for spreadsheet.
                    # std = np.std(de)
                    # spreadsheet.append([name, facet, rname, pname,
                    #                     DeltaE, std, ea,
                    #                     code, xc,
                    #                     publication, url])
                    # width, height, angle, covariance
                    Nrxn += 1
                else:
                    continue
            Nsurf += 1
        # with open(project + '/data.csv', 'wb') as f:
        #     writer = csv.writer(f)
        #     writer.writerows(spreadsheet)
        print(Nrxn, 'reactions imported.')
        print(Nsurf, 'surfaces saved.')

    def _state2species(self, state):
        """Parse one side of a CatMAP rxn expression, i.e. a chemical state.

        Parameters
        ----------
        state : str
            Left or right side of CatMAP rxn expression, excluding arrows.

        Returns
        ----------
        state_name : str
            Name of chemical state formatted for folder naming.
        slist : list
            List of species. <coefficient><structure formula or Hill formula>.
        """
        slist = state.split('+')
        species = []
        for specie in slist:
            if '_g' in specie:
                # Gas species.
                species.append(specie.split('_')[0] + 'gas')
            elif '*_' in specie:
                # Empty sites.
                species.append(specie.split('_')[0].replace('*', 'star'))
            else:
                # Adsorbates.
                species.append(specie.split('_')[0] + 'star')
        return '_'.join(species), slist

    def _coefficient_species(self, species):
        """Return the stochiometric coefficient and the species type.

        Parameters
        ----------
        species : str
            <coefficient><structure formula or Hill formula>.

        Returns
        ----------
        n : int
            Stochiometric coefficient.
        species : str
            Species name.
        """
        i = 0
        n = 1
        if species[0] == '-':
            # Negative coefficient allowed.
            i += 1
            n = -1
        while species[i].isdigit():
            i += 1
        if i > 1:
            n = int(species[:i])
        elif i == 1 and n != -1:
            n = int(species[0])

        return n, species[i:]
