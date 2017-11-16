# -*- coding: utf-8 -*-
"""
Description:
    Module to convert ase.db files into a catmap input file

    <surface_name>.db: db file
        Contains all adsorbate states on
        surfaces and clean slabs for reference.
        Mandatory key value pairs:
        --------------------------
            "name" : str
                Value that identifies the catalyst composition.
            "species" : str
                value of the adsorbate chemical formula.
                '' should be the value for clean slabs.
                '-' should be inserted between seperate fragments.
            "epot" : float
                value with potential energy from DFT.

        Recommended key value pairs:
        ----------------------------
            "site" : str
                name of adsorption site.
            "phase" : str
                Value that identifies the catalyst phase.
            "facet" : str
                Name of the facet,
                preferably in hkl notation and separated by 'x', e.g. '0x0x0'.
            "surf_lattice" : str
                Name of the surface lattice geometry.
                E.g. hcp001 and fcc111 has "hexagonal" surface lattices.
            "layers" : int
                Number of atomic layers in slab.
            "supercell" : str
                Supercell size separated by x, e.g. 2x2
            "n": int
                number of identical adsorbates.

        Recommended keys in "data":
        ---------------------------
            "BEEFvdW" : list
                32 non-selfconsistent BEEF-vdW energies.

    Optional file dependencies:
    ---------------------------
        json files in the working directory containing a dictionary with
        frequencies. The file name and the key name should be
            species + '_' + facet + '.freq'
        # To do:
            Expand this feature to individual frequency lists for individual
            catalysts.
"""
from os import listdir, environ, mkdir
import os.path
import numpy as np
import ase.db
from ase.atoms import string2symbols
from catmap.bee import BEEFEnsemble as BEE
import csv


class db2catmap(object):
    def __init__(self, mol_db=None, ads_db=None, mol_select=[], ads_select=[],
                 frequency_db=None, frequency_select=[],
                 ts_db=None, ts_selection=[]):
        self.state = BEE()
        self.epot = {}
        self.freq = {}
        self.de = {}
        self.dbid = {}
        self.reference_epot = {}
        self.reference_de = {}
        self.rxn_paths = {}
        self.formation_energies = {}
        self.std = {}
        if mol_db is not None:
            self.mol_db = mol_db
            self.get_molecules(mol_db, selection=mol_select)
            if ads_db is not None:
                self.ads_db = ads_db
                self.get_surfaces(ads_db, selection=ads_select)
                if ts_db is not None:
                    self.get_transition_states(ts_db, selection=ts_selection)
        self.ts_db = ts_db

    def get_molecules(self, fname, selection=[], frequency_db=None):
        """ Method for importing molecules.

        Parameters
        ----------
        fname : str
            path and filename of an ase database file containing molecules.
        """
        [mol_epot,
         mol_freq,
         mol_de,
         mol_dbid] = self._db2mol(fname, selection=selection,
                                  freq_path=frequency_db)
        self.epot.update(mol_epot)
        self.freq.update(mol_freq)
        self.de.update(mol_de)
        self.dbid.update(mol_dbid)

    def get_surfaces(self, fname, selection=[], frequency_db=None,
                     site_specific=False):
        """ Method for importing slabs and adsorbates.

        Parameters
        ----------
        fname : str
            path and filename of an ase database file containing slabs.
        """
        [surf_epot,
         surf_freq,
         surf_de,
         surf_dbid] = self._db2surf(fname, selection=selection,
                                    freq_path=frequency_db,
                                    site_specific=site_specific)
        self.epot.update(surf_epot)
        self.freq.update(surf_freq)
        self.de.update(surf_de)
        self.dbid.update(surf_dbid)

    def get_transition_states(self, fname, selection=[], frequency_db=None):
        """ Method for importing surface transition states.

        Parameters
        ----------
        fname : str
            path and filename of an ase database file
            containing reaction images.
        """

        rxn_paths = self._db2pes(fname, selection=selection)
        self.rxn_paths.update(rxn_paths)
        [surf_epot,
         surf_freq,
         surf_de,
         surf_dbid] = self.pes2ts(freq_path=frequency_db)
        self.epot.update(surf_epot)
        self.freq.update(surf_freq)
        self.de.update(surf_de)
        self.dbid.update(surf_dbid)

    def calc_formation_energies(self, references=(('H', 'H2_gas'),
                                                  ('O', 'H2O_gas'),
                                                  ('C', 'CH4_gas'),)):
        """ Method for generating formation energies.

        Parameters
        ----------
        references : list of tuples of strings.
            The first item in each tuple must be an atomic symbol, and the
            second item in each tuple must be a reference <species name>_gas
        """
        # Get atomic reference energies.
        self._mol2ref(references=references)
        # Make a dictionary with slab references and atomic references.
        self._get_refs()
        self._get_BEEstd()
        self._get_formation_energies()

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
        cmol = ase.db.connect(fname)
        smol = cmol.select(selection)
        if freq_path is not None:
            c_freq = ase.db.connect(freq_path)
        abinitio_energies = {}
        freq_dict = {}
        dbids = {}
        de = {}
        for d in smol:              # Get molecules from mol.db
            if 'energy' in d:
                abinitio_energy = float(d.energy)
            else:
                abinitio_energy = float(d.epot)
            species_name = str(d.formula)
            try:
                contribs = d.data.BEEFvdW_contribs
                ens = self.state.get_ensemble_perturbations(contribs)
            except AttributeError:
                ens = 0  # np.zeros(self.state.size)
            if species_name+'_gas' not in abinitio_energies:
                abinitio_energies[species_name+'_gas'] = abinitio_energy
                dbids[species_name+'_gas'] = int(d.id)
                de[species_name+'_gas'] = ens
                if freq_path is not None:
                    try:
                        d_freq = c_freq.get(['formula='+species_name])
                        frequencies = d_freq.data.frequencies
                        freq_dict.update({species_name + '_gas': frequencies})
                    except KeyError:
                        continue
            elif abinitio_energies[species_name+'_gas'] > abinitio_energy:
                abinitio_energies[species_name+'_gas'] = abinitio_energy
                dbids[species_name+'_gas'] = int(d.id)
                de[species_name+'_gas'] = contribs
        return abinitio_energies, freq_dict, de, dbids

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
            fname : str
                path/filname.
            selection : list of strings.
                Optional ase.db selection strings.
            site_specific : boolean
                If True: Dinstinguish sites using the site key value pair, and
                stores a the potential energy of adsorbates on each site.
                Else: Use the minimum ab initio energy, disregarding the site.
        """
        csurf = ase.db.connect(fname)
        if freq_path is not None:
            c_freq = ase.db.connect(freq_path)
        ssurf = csurf.select(selection)
        abinitio_energies = {}
        freq_dict = {}
        dbids = {}
        de = {}
        for d in ssurf:                     # Get slab and adsorbates from .db
            species = str(d.species)
            if '-' in species:
                continue
            ads = str(d.ads)
            name = str(d.name)
            if 'n' in d and int(d.n) > 1:  # Skip higher coverages for now.
                continue
            if 'energy' in d:
                abinitio_energy = float(d.energy)
            else:
                abinitio_energy = float(d.epot)
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
            if 'n' in d:
                n = int(d.n)
            else:
                n = 1
            # composition=str(d.formula)
            if species == '' or ('ads' in d and
                                 (ads == 'slab' or ads == 'clean')):
                species = ''
                site = 'slab'
                n = 0
            elif 'site' in d:
                site = str(d.site)
            else:
                site = 'site'
            cat = name + '_' + phase
            site_name = surf_lattice + '_' + facet + '_' + cell + '_' + site
            key = str(n) + '_' + species + '_' + cat + '_' + site_name
            try:
                contribs = d.data.BEEFvdW_contribs
                ens = self.state.get_ensemble_perturbations(contribs)
            except:
                ens = 0  # np.zeros(self.state.size)
            if key not in abinitio_energies:
                abinitio_energies[key] = abinitio_energy
                dbids[key] = int(d.id)
                de[key] = ens
                if species != '' and ads != 'slab' and freq_path is not None:
                    try:
                        d_freq = c_freq.get(['species='+species, 'name='+name])
                        frequencies = d_freq.data.frequencies
                        freq_dict.update({key: frequencies})
                    except KeyError:
                        continue
            elif abinitio_energies[key] > abinitio_energy:
                abinitio_energies[key] = abinitio_energy
                dbids[key] = int(d.id)
                de[key] = ens
        return abinitio_energies, freq_dict, de, dbids

    def _db2pes(self, fname, selection=[]):
        """ Returns a dictionary containing potential energy surfaces and metadata.

            File dependencies
            -----------------
            fname : ase.db file.

            Parameters
            ----------
            fname : str
                path/filname.
            selection : list of strings.
                Optional ase.db selection strings.
        """
        c = ase.db.connect(fname)
        s = c.select(selection)
        rxn_paths = {}
        # Get images from ase .db
        for d in s:
            species = str(d.species)
            if '-' not in species:
                continue
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
            surf = str(d.name) + '_' + phase + '_' + \
                surf_lattice + '_' + facet + '_' + cell
            if 'energy' in d:
                abinitio_energy = float(d.energy)
            else:
                abinitio_energy = float(d.epot)
            dbid = int(d.id)
            try:
                BEEFvdW_contribs = d.data.BEEFvdW_contribs
                ens = self.state.get_ensemble_perturbations(BEEFvdW_contribs)
            except AttributeError:
                ens = 0  # np.zeros(self.state.size)
            rxn_id = str(d.path_id)
            if rxn_id in rxn_paths:
                rxn_paths[rxn_id]['pes'].append(abinitio_energy)
                rxn_paths[rxn_id]['dbids'].append(dbid)
                rxn_paths[rxn_id]['de'].append(ens)
                try:
                    rxn_paths[rxn_id]['distance'].append(float(d.distance))
                except AttributeError:
                    # Do nothing.
                    continue
                if 'image' in d:
                    rxn_paths[rxn_id]['images'].append(int(d.image))
                else:
                    rxn_paths[rxn_id]['images'].append(int(d.step))
            else:
                if 'site' in d:
                    site = str(d.site)
                else:
                    site = 'site'
                rxn_paths[rxn_id] = {'surface_name': surf,
                                     'species': species,
                                     'pes': [abinitio_energy],
                                     'dbids': [dbid],
                                     'de': [ens],
                                     'site': site,
                                     'images': [int(d.step)],
                                     'distance': [float(d.distance)]}
        return rxn_paths

    def _mol2ref(self,
                 references=(('H', 'H2_gas'), ('O', 'H2O_gas'),
                             ('C', 'CH4_gas'),)):
        """ Returns two dictionaries containing:
            abinitio energy references for atoms
            32 non-selfconsitent perturbations for atoms.
        """
        atomic_e = {}
        atomic_de = {}
        print(references)
        for t in references:
            key = t[0]
            species = t[1]
            atomic_e[key] = self.epot[species]
            atomic_de[key] = np.array(self.de[species])
            composition = string2symbols(species.split('_')[0])
            n = 0.
            for symbol in composition:
                if symbol == key:
                    n += 1
                else:
                    atomic_e[key] -= atomic_e[symbol]
                    atomic_de[key] -= np.array(atomic_de[symbol])
            atomic_e[key] /= n
            atomic_de[key] /= n
        self.atomic_e = atomic_e
        self.atomic_de = atomic_de
        # return atomic_e, atomic_de

    def _get_refs(self):
        """ Returns dictionaries with referece energies of slabs and
        single atoms.
        """
        ref_dict = self.atomic_e
        ref_de = self.atomic_de
        for key in self.epot.keys():
            if 'slab' in key:
                ref_dict[key] = self.epot[key]
                ref_de[key] = self.de[key]
        self.reference_epot = ref_dict
        self.reference_de = ref_de

    def _get_formation_energies(self):
        """ Returns a dictionary with formation energies of adsorbates.
        Parameters
        ----------
        energy_dict : dictionary
            Each key is named in the format: adsorbate_name_phase_facet_site,
            and contains the potential energy of a slab an adsorbate or
            a molecule.
        ref_dict : dictionary
            Each key is either an atomic symbol and contains the reference
            potential energy of that atom,
            or the key is named in the format: _name_phase_facet_slab and it
            contains the reference potential energy of the slab.
        """
        formation_energies = {}
        for key in self.epot.keys():
            E0 = 0
            if 'gas' in key:
                species, site_name = key.split('_')  # Split key into name/site
            else:
                n, species, cat, pha, lattice, fac, cell, site = key.split('_')
                site_name = '0__' + cat + '_' + pha + '_' + lattice + '_' + \
                    fac + '_' + cell + '_slab'
                try:
                    E0 -= self.reference_epot[site_name]
                except KeyError:
                    continue
            if 'slab' not in key:
                composition = string2symbols(species.replace('-', ''))
                E0 += self.epot[key]
                for atom in composition:
                    E0 -= self.reference_epot[atom]
                formation_energies[key] = E0
        self.formation_energies.update(formation_energies)

    def _get_BEEstd(self):
        """ Returns dictionary with BEEF standard deviations on formation energies.
        """
        BEEstd = {}
        for key in self.de:
            de = np.zeros(self.state.size)
            if 'gas' in key:
                species, site_name = key.split('_')  # Split key into name/site
            else:
                n, species, cat, pha, lattice, fac, cell, site = key.split('_')
                site_name = '0__' + cat + '_' + pha + '_' + lattice + '_' + \
                    fac + '_' + cell + '_slab'
                try:
                    de -= self.reference_de[site_name]
                except KeyError:
                    continue
            if 'slab' not in key:
                composition = string2symbols(species.replace('-', ''))
                de += self.de[key]
                for atom in composition:
                    de -= self.reference_de[atom]
                BEEstd[key] = de.std()
        self.std.update(BEEstd)

    def get_BEE_PCA(self, ads_x, ads_y,
                    site_x='site', site_y='site'):
        """ Returns two dictionaries, BEE_PC with the principal component
        vectors and BEE_lambda with the corresponding eigenvalues of
        BEE pertubations.

        Parameters
        ----------
        de_dict  dict    contains beef perturbations of adsorbates on slabs,
                            where keys are named: adsorbate_name_phase_facet
        ref_de   dict    contains beef perturbations of references,
                            where keys are refernce elements, e.g: 'C','H',
                            and also slabs references: name_phase_facet
        ads_x    string  adsorbate first dimension
        ads_y    string  adsorbate second dimension
        """
        BEE_PC = {}
        BEE_lambda = {}
        # Loop over reference surfaces in ref_de.
        for slab in self.reference_epot.keys():
            # Ignore gas species.
            if 'gas' in slab or len(slab) <= 2:
                continue
            de_x = -self.reference_de[slab]
            de_y = -self.reference_de[slab]
            try:
                de_x += self.de[ads_x + '_' + slab[:-4] + site_x]
                de_y += self.de[ads_y + '_' + slab[:-4] + site_y]
            except KeyError:
                print('Missing BEEF ensemble for ', ads_x, ads_y, slab[:-4],
                      site_x, site_y)
                continue
            composition_x = string2symbols(ads_x)
            composition_y = string2symbols(ads_y)
            for atom in composition_x:
                de_x -= self.reference_de[atom]
            for atom in composition_y:
                de_y -= self.reference_de[atom]
            cov = np.cov(de_x, de_y)
            eigval, eigvec = np.linalg.eig(cov)
            BEE_PC[slab] = eigvec
            BEE_lambda[slab] = eigval / len(cov)
        return BEE_PC, BEE_lambda

    def pes2ts(self, freq_path=None):
        """ Returns dictionaries containing transition state ab initio energies.

            Input parameters
            ----------------
            freq_path : str
                path/folder where frequency files are located.
        """
        if freq_path is not None:
            c_freq = ase.db.connect(freq_path)
        abinitio_energies = {}
        freq_dict = {}
        dbids = {}
        de = {}
        calculate = []
        for rxn_id in self.rxn_paths:
            species = self.rxn_paths[rxn_id]['species']
            m = self.rxn_paths[rxn_id]['surface_name']
            site = self.rxn_paths[rxn_id]['site']
            key = '0_' + species + '_' + m + '_' + site
            images = self.rxn_paths[rxn_id]['images']
            if len(np.unique(images)) != len(images):
                print('non unique image number!', m, species, images)
                continue
            pes = np.array(self.rxn_paths[rxn_id]['pes'])
            s = np.argsort(images)
            # Look for local minima and maxima.
            localmins = np.where(np.r_[True, pes[s][1:] < pes[s][:-1]] &
                                 np.r_[pes[s][:-1] < pes[s][1:], True])[0]
            localmaxs = np.where(np.r_[True, pes[s][1:] > pes[s][:-1]] &
                                 np.r_[pes[s][:-1] > pes[s][1:], True])[0]
            differences = np.diff(pes[s])
            roughness = np.std(differences)
            if len(localmaxs) > 1 or len(localmins) > 2 or len(localmins) == 1:
                print('Warning!:')
            if np.min(self.rxn_paths[rxn_id]['distance']) > 1.6:
                print('Incomplete trajectory!', m, species, images)
                calculate.append(max(self.rxn_paths[rxn_id]['dbids']))
                continue
            tst = np.argmax(pes)
            print(species, m, roughness, len(localmaxs), len(localmins))
            if key not in abinitio_energies:
                abinitio_energies[key] = pes[tst]
                dbids[key] = self.rxn_paths[rxn_id]['dbids'][tst]
                de[key] = self.rxn_paths[rxn_id]['de'][tst]
                if freq_path is not None:
                    try:
                        d_freq = c_freq.get('path_id='+rxn_id)
                        frequencies = d_freq.data.frequencies
                        freq_dict.update({key: frequencies})
                    except KeyError:
                        continue
            elif abinitio_energies[key] > pes[tst]:
                abinitio_energies[key] = pes[tst]
                dbids[key] = self.rxn_paths[rxn_id]['dbids'][tst]
                de[key] = self.rxn_paths[rxn_id]['de'][tst]
        incomplete = ','.join([str(int(a)) for a in np.unique(calculate)])
        print(incomplete)
        return abinitio_energies, freq_dict, de, dbids

    def make_input_file(self, file_name, site_specific=False):
        """ Saves the catmap input file.

        Parameters
        ----------
        file_name : string
            path and name of the output file.
        energy_dict : dictionary
            Contains formation energies. Each species is represented by a field
            and the fields must be named in the format
            <species>_<surface>_<phase>_<surface lattice>_<facet>_<site>
            or
            <species>_gas
        freq_dict : dictionary
            Contains lists of frequencies. Each species is represented by a
            field and the fields must be named corresponding to energy_dict.
        bee_dict : dictionary
            Contains standard deviations on the energy.
            Each species is represented by a field and the fields must be named
            corresponding to energy_dict.
        """
        # Create a header.
        header = '\t'.join(['surface_name', 'phase', 'site_name',
                            'species_name', 'formation_energy',
                            'frequencies', 'reference', 'coverage', 'std'])
        # List of lines in the output.
        lines = []
        for key in self.formation_energies.keys():  # Iterate through keys
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
            if 'slab' not in key:  # Do not include empty site energy (0)
                try:
                    frequency = self.freq[key]
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
                    coverage = 0
                else:
                    surface = cat  # +'_'+pha #useful to include phase.
                    phase = pha
                    coverage = round(float(n) /
                                     (float(cell[0]) * float(cell[2])), 3)
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
                           E, list(frequency), environ['USER'], coverage, std]
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

    def make_nested_folders(self, project, reactions, surfaces=None,
                            site='site', ads_db=None, mol_db=None, ts_db=None,
                            publication='', url=''):
        """Saves a nested directory structure,
        compatible with the catapp project.

        Parameters
        ----------
        fname : string
            path and name of an ase database file.
        project : dictionary
        reactions : list
        """
        # Create a header
        spreadsheet = [['chemical_composition', 'facet', 'reactants',
                        'products', 'reaction_energy',
                        'beef_standard_deviation',
                        'activation_energy', 'DFT_code', 'DFT_functional',
                        'reference', 'url']]
        if surfaces is None:
            surfaces = [s for s in self.reference_epot.keys() if 'slab' in s]
        if mol_db is None:
            mol_db = self.mol_db
        if ads_db is None:
            ads_db = self.ads_db
        if ts_db is None:
            ts_db = self.ts_db
        c_ads = ase.db.connect(ads_db)
        c_mol = ase.db.connect(mol_db)
        if ts_db is not None:
            c_ts = ase.db.connect(ts_db)
        Nsurf = 0
        Nrxn = 0
        for i in range(len(reactions)):
            states = reactions[i].replace(' ', '').split('<->')
            if len(states) == 1:
                states = states[0].split('->')
                if len(states) == 1:
                    states = states[0].split('<-')
            elif len(states) < 3:
                states = [states[0]] + states[-1].split('->')
                if len(states) < 3:
                    states = states[0].split('<-') + states[1:]
            reactants = states[0].split('+')
            tstates = states[1].split('+')
            products = states[-1].split('+')
            rspecies = [r.replace('*',
                                  'star').split('_')[0] for r in reactants]
            pspecies = [p.replace('*',
                                  'star').split('_')[0] for p in products]
            rname = '_'.join(rspecies)
            pname = '_'.join(pspecies)
            reaction_name = '__'.join([rname, pname])
            path_reaction = project + '/' + reaction_name
            for slabkey in surfaces:
                totraj = {}
                [n, species, name, phase,
                 lattice, facet, cell, slab] = slabkey.split('_')
                path_surface = path_reaction + '/' + name.replace('/', '')
                path_facet = path_surface + '/' + facet.replace('x', '')
                DeltaE = 0.
                de = 0.  # np.zeros(self.state.size)
                ea = 0.
                intermediates_exist = True
                # Find reactant structures and energies
                for reactant in reactants:
                    species, sitesymbol = reactant.split('_')
                    if species[0].isdigit():
                        n = int(species[0])
                        species = species[1:]
                    else:
                        n = 1
                    if sitesymbol == 'g':
                        rkey = species + '_gas'
                        fname = path_reaction + '/' + rkey
                    elif species == '*':
                        rkey = slabkey
                        fname = path_facet + '/empty_surface'
                    else:
                        rkey = '1_' + species.replace('*', '') + \
                            '_' + name + '_' + phase + '_' + lattice + '_' + \
                            facet + '_' + cell + '_' + site
                        fname = path_facet + '/' + species
                    if rkey not in self.dbid:
                        intermediates_exist = False
                        break
                    totraj.update({rkey:
                                   {'dbid': self.dbid[rkey],
                                    'fname': fname}})
                    if species != '*':
                        DeltaE -= n * self.formation_energies[rkey]
                        de -= n * self.de[rkey]
                if not intermediates_exist:
                    continue
                # Find transition state structures and energies.
                if ts_db is not None:
                    for ts in tstates:
                        if '-' not in ts:
                            continue
                        species, sitesymbol = ts.split('_')
                        tskey = '0_' + species + '_' + name + '_' + phase + \
                            '_' + lattice + '_' + facet + '_' + cell + \
                            '_' + site
                        if tskey not in self.dbid:
                            continue
                        totraj.update({tskey:
                                       {'dbid': self.dbid[tskey],
                                        'fname': path_facet + '/' + ts}})
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
                        fname = path_reaction + '/' + pkey
                    elif species == '*':
                        pkey = slabkey
                        fname = path_facet + '/empty_surface'
                    else:
                        pkey = '1_' + species.replace('*', '') + \
                            '_' + name + '_' + phase + '_' + lattice + '_' + \
                            facet + '_' + cell + '_' + site
                        fname = path_facet + '/' + species
                    if pkey not in self.dbid:
                        intermediates_exist = False
                        break
                    totraj.update({pkey:
                                   {'dbid': self.dbid[pkey],
                                    'fname': fname}})
                    if species != '*':
                        DeltaE += n * self.formation_energies[pkey]
                        de += n * self.de[pkey]
                # If all states are found for this surface, write.
                if intermediates_exist:
                    if not os.path.isdir(project):
                        os.mkdir(project)
                    if reaction_name not in listdir(project):
                        os.mkdir(path_reaction)
                    if name.replace('/', '') not in listdir(path_reaction):
                        os.mkdir(path_surface)
                    if facet.replace('x', '') not in listdir(path_surface):
                        os.mkdir(path_facet)
                    for trajkey in totraj.keys():
                        if 'gas' in trajkey:
                            fname = totraj[trajkey]['fname'] + '.traj'
                            if fname in os.listdir(path_reaction):
                                continue
                            else:
                                atoms = c_mol.get_atoms(self.dbid[trajkey])
                        elif '-' in trajkey.split('_')[1]:
                            atoms = c_ts.get_atoms(self.dbid[trajkey])
                        else:
                            atoms = c_ads.get_atoms(self.dbid[trajkey])
                        atoms.write(totraj[trajkey]['fname'] + '.traj')
                    std = np.std(de)
                    print(DeltaE)
                    spreadsheet.append([name, facet, rname, pname,
                                        DeltaE, std, ea,
                                        'Quantum Espresso',
                                        'BEEF-vdW', publication, url])

                    Nsurf += 1
                else:
                    continue
            Nrxn += 1
        with open(project+'test.csv', 'wb') as f:
            writer = csv.writer(f)
            writer.writerows(spreadsheet)
        print(Nrxn, 'reactions imported.')
        print(Nsurf, 'surfaces saved.')
