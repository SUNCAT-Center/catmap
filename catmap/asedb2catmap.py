# -*- coding: utf-8 -*-
"""
Description:
    Module to convert ase.db files into a catmap input file

    File dependencies:
    ------------------
    <mol_name>.db: db file
        Contains molecular reference states.
        Mandatory key value pairs:
        --------------------------
                "epot" : float
                    DFT calculated potential energy.
        Optional key value pairs:
        --------------------------
                "data.BEEFens" : list
                    32 non-selfconsistent BEEF-vdW energies.

    <surface_name>.db: db file
        Contains all adsorbate states on
        surfaces and clean slabs for reference.
        Mandatory key value pairs:
        --------------------------
            "name" : str
                Value that identifies the catalyst composition.
            "phase" : str
                Value that identifies the catalyst phase.
            "facet" : str
                Name of the facet,
                preferably in hkl notation and separated by 'x', e.g. '0x0x0'.
            "surf_lattice" : str
                Name of the surface lattice geometry.
                E.g. hcp001 and fcc111 has "hexagonal" surface lattices.
            "species" : str
                value of the adsorbate chemical formula.
                '' should be the value for clean slabs.
                '-' should be inserted between seperate fragments.
            "epot" : float
                value with potential energy from DFT.
            "supercell" : str
                Supercell size separated by x, e.g. 2x2
            "layers" : int
                Number of atomic layers in slab.

        Recommended key value pairs:
        ----------------------------
            "site" : str
                name of adsorption site.

        Recommended keys in "data":
        ---------------------------
            "BEEFens" : list
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
from os import mkdir, listdir, environ
import numpy as np
import ase.db
from ase.atoms import string2symbols
import json
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
            self.attach_molecules(mol_db, select=mol_select)
            if ads_db is not None:
                self.attach_surfaces(ads_db, select=ads_select)
                if ts_db is not None:
                    self.rxn_paths = self._db2pes(ts_db,
                                                  select=ts_selection)

    def get_molecules(self, fname, selection=[]):
        [mol_epot,
         mol_freq,
         mol_de,
         mol_dbid] = self._db2mol(fname, selection=selection,
                                  freq_path='.')
        self.epot.update(mol_epot)
        self.freq.update(mol_freq)
        self.de.update(mol_de)
        self.dbid.update(mol_dbid)

    def get_surfaces(self, fname, selection=[]):
        [surf_epot,
         surf_freq,
         surf_de,
         surf_dbid] = self._db2surf(fname, selection=selection,
                                    freq_path='.')
        self.epot.update(surf_epot)
        self.freq.update(surf_freq)
        self.de.update(surf_de)
        self.dbid.update(surf_dbid)

    def get_transition_states(self, fname, selection=[]):
        rxn_paths = self._db2pes(fname, selection=selection)
        self.rxn_paths.update(rxn_paths)
        [surf_epot,
         surf_freq,
         surf_de,
         surf_dbid] = self.pes2ts()
        self.epot.update(surf_epot)
        self.freq.update(surf_freq)
        self.de.update(surf_de)
        self.dbid.update(surf_dbid)

    def calc_formation_energies(self,
                                references={'O': 'H2O_gas', 'C': 'CH4_gas'}):
        # Get atomic reference energies.
        self._mol2ref(references=references)
        # Make a dictionary with slab references and atomic references.
        self._get_refs()
        self._get_BEEstd()
        self._get_formation_energies()

    def _db2mol(self, fname, selection=[], freq_path='.'):
        """ Returns four dictionaries containing:
            ab initio energies,
            frequecies,
            non-selfconsistent BEEF perturbations,
            database ids.

            File dependencies
            -----------------
                fname : ase.db file.
        """
        cmol = ase.db.connect(fname)
        smol = cmol.select(selection)
        abinitio_energies = {}
        freq_dict = {}
        dbids = {}
        de = {}
        for d in smol:              # Get molecules from mol.db
            abinitio_energy = float(d.epot)
            species_name = str(d.formula)
            if species_name+'_gas' not in abinitio_energies:
                abinitio_energies[species_name+'_gas'] = abinitio_energy
                dbids[species_name+'_gas'] = int(d.id)
                de[species_name+'_gas'] = \
                    self.state.get_ensemble_perturbations(d.data.BEEFens)
                try:
                    freq = json.load(open(freq_path + '/' + species_name +
                                          '_gas.freq', 'r'))
                    freq_dict.update(freq)
                except IOError:
                    print('no frequencies for', species_name, '(g)')
            elif abinitio_energies[species_name+'_gas'] > abinitio_energy:
                abinitio_energies[species_name+'_gas'] = abinitio_energy
                dbids[species_name+'_gas'] = int(d.id)
                de[species_name+'_gas'] = \
                    self.state.get_ensemble_perturbations(d.data.BEEFens)
        return abinitio_energies, freq_dict, de, dbids

    def _db2surf(self, fname, selection=[], freq_path='.', sites=False):
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
            sites : boolean
                If True: Dinstinguish sites using the site key value pair,
                    in which case the site key value pair becomes mandatory.
                Else: Use the minimum ab initio energy regardless of site.
        """
        csurf = ase.db.connect(fname)
        ssurf = csurf.select(selection)
        abinitio_energies = {}
        freq_dict = {}
        dbids = {}
        de = {}
        for d in ssurf:                     # Get slab and adsorbates from .db
            ads = str(d.ads)
            species = str(d.species)
            if '-' in species:
                continue
            if 'n' in d and int(d.n) > 1:  # Skip higher coverages for now.
                continue
            cat = str(d.name)+'_'+str(d.phase)
            abinitio_energy = float(d.epot)
            surf_lattice = str(d.surf_lattice)
            cell = str(d.supercell) + 'x' + str(d.layers)
            n = int(d.n)
            # composition=str(d.formula)
            facet = str(d.facet)
            if species == '' or ('ads' in d and ads == 'slab'):
                species = ''
                site = 'slab'
            elif sites:
                site = str(d.site)
            else:
                site = 'site'
            site_name = surf_lattice + '_' + facet + '_' + cell + '_' + site
            key = str(n) + '_' + species + '_' + cat + '_' + site_name
            if key not in abinitio_energies:
                abinitio_energies[key] = abinitio_energy
                dbids[key] = int(d.id)
                de[key] = self.state.get_ensemble_perturbations(d.data.BEEFens)
                if species != '' and ads != 'slab':
                    try:
                        freq = json.load(open(freq_path + '/' + species + '_' +
                                              surf_lattice + '.freq', 'r'))
                        freq_dict.update({key: freq[species +
                                                    '_' + surf_lattice]})
                    except IOError:
                        print('no frequencies for',
                              species + '_' + site_name)
            elif abinitio_energies[key] > abinitio_energy:
                abinitio_energies[key] = abinitio_energy
                dbids[key] = int(d.id)
                de[key] = self.state.get_ensemble_perturbations(d.data.BEEFens)
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
            cell = str(d.supercell) + 'x' + str(d.layers)
            surf = str(d.name) + '_' + str(d.phase) + '_' + \
                str(d.surf_lattice) + '_' + str(d.facet) + '_' + cell
            abinitio_energy = float(d.epot)
            dbid = int(d.id)
            try:
                BEEFens = d.data.BEEFens
                de = self.state.get_ensemble_perturbations(BEEFens)
            except AttributeError:
                print(dbid, species, surf, 'missing BEEF perturbations.')
                de = np.zeros(2000)
            rxn_id = str(d.path_id)
            if rxn_id in rxn_paths:
                rxn_paths[rxn_id]['pes'].append(abinitio_energy)
                rxn_paths[rxn_id]['dbids'].append(dbid)
                rxn_paths[rxn_id]['de'].append(de)
                try:
                    rxn_paths[rxn_id]['distance'].append(float(d.distance))
                except AttributeError:
                    continue
                    # Do nothing.
                if 'image' in d:
                    rxn_paths[rxn_id]['images'].append(int(d.image))
                else:
                    rxn_paths[rxn_id]['images'].append(int(d.step))
            else:
                rxn_paths[rxn_id] = {'surface_name': surf,
                                     'species': species,
                                     'pes': [abinitio_energy],
                                     'dbids': [dbid],
                                     'de': [de],
                                     'images': [int(d.step)],
                                     'distance': [float(d.distance)]}
        return rxn_paths

    def _mol2ref(self,
                 references=(('H', 'H2_gas'), ('O', 'H2O_gas'),
                             ('C', 'CH4_gas'))):
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
            atomic_de[key] = self.de[species]
            composition = string2symbols(species.split('_')[0])
            n = 0.
            for symbol in composition:
                if symbol == key:
                    n += 1
                else:
                    atomic_e[key] -= atomic_e[symbol]
                    atomic_de[key] -= atomic_de[symbol]
            atomic_e[key] /= n
            atomic_de[key] /= n
        print(atomic_e)
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
            de = np.zeros(2000)
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

    def pes2ts(self, freq_path='.'):
        """ Returns dictionaries containing transition state ab initio energies.

            Input parameters
            ----------------
            rxn_paths : dictionary
                Created by the db2pes function.
        """
        abinitio_energies = {}
        freq_dict = {}
        dbids = {}
        de = {}
        calculate = []
        for rxn_id in self.rxn_paths:
            species = self.rxn_paths[rxn_id]['species']
            m = self.rxn_paths[rxn_id]['surface_name']
            key = '0_' + species + '_' + m + '_site'
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
                try:
                    freq = json.load(open(freq_path + '/' + species + '_' + m +
                                          '.freq', 'r'))
                    freq_dict.update({species + '_' + m + '_' + key:
                                      freq[species + '_' + m]})
                except IOError:
                    print('no frequencies for', key)
            elif abinitio_energies[key] > pes[tst]:
                abinitio_energies[key] = pes[tst]
                dbids[key] = self.rxn_paths[rxn_id]['dbids'][tst]
                de[key] = self.rxn_paths[rxn_id]['de'][tst]
        incomplete = ','.join([str(int(a)) for a in np.unique(calculate)])
        print(incomplete)
        return abinitio_energies, freq_dict, de, dbids

    def make_input_file(self, file_name):
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
                    std = self.std[key]
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
                    std = self.std[key]
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
                    site_name = lattice + '_' + site
                    coverage = round(float(n) /
                                     (float(cell[0]) * float(cell[2])), 3)
                E = float(E)
                outline = [surface, phase, site_name, name,
                           E, frequency, environ['USER'], coverage, std]
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

    def make_catapp_files(self, fname, project):
        """ Saves the catmap input file.

        Parameters
        ----------
        file_name : string
            path and name of an ase database file.
        energy_dict : dictionary
            Contains formation energies. Each species is represented by a field
            and the fields must be named in the format
            <species>_<surface>_<phase>_<surface lattice>_<facet>_<site>
            or
            <species>_gas
        freq_dict : dictionary
            Contains lists of frequencies.
            Each species is represented by a field and the fields must be named
            corresponding to energy_dict.
        bee_dict : dictionary
            Contains standard deviations on the energy.
            Each species is represented by a field and the fields must be named
            corresponding to energy_dict.
        """
        # Create a header
        c = ase.db.connect(fname)
        spreadsheet = [['chemical_composition', 'facet', 'AB', 'A', 'B',
                        'reaction_energy', 'beef_standard_deviation',
                        'activation_energy', 'DFT_code', 'DFT_functional',
                        'reference', 'url']]
        for key in self.dbid.keys():  # Iterate through keys
            if key in self.epot:
                Ef = self.epot[key]  # formation energy
                if 'gas' in key:
                    name, site = key.split('_')  # Split key into name/site
                    # try:
                    #     frequency = freq_dict[key]
                    # except KeyError:
                    #     frequency = []
                    # try:
                    #     bee = bee_dict[key]
                    # except KeyError:
                    #     bee = np.NaN
                else:
                    name, cat, pha, lattice, facet, site = key.split('_')
                if 'slab' not in key:  # Do not include empty site energy (0)
                    print(key)
                    # try:
                    #     frequency = freq_dict[key]
                    # except KeyError:
                    #     frequency = []
                    try:
                        bee = self.de[key]
                    except KeyError:
                        bee = np.NaN
                    if site == 'gas':
                        continue
                    else:
                        surface = cat  # +'_'+pha #useful to include phase.
                        # phase = pha
                        # site_name = lattice + '_' + site
                    composition = string2symbols(name)
                    nC = composition.count('C')
                    nH = composition.count('H')
                    AB = name
                    if nC is not 0:
                        A = (str(int(nC)) +
                             'CH4gas').replace('1.0', '').replace('1', '')
                    else:
                        A = ''
                    B = (str(nH/2. - 2 * nC) +
                         'H2gas').replace('1.0', '')
                    A_B = A + '_' + B
                    path_reaction = project + AB + 'star_' + A_B + '_star'
                    path_surface = path_reaction + '/' + surface + '_' + facet
                    atoms = c.get_atoms(self.epot[key])
                    slab = c.get_atoms(self.epot['_' + cat + '_' + pha +
                                       '_' + lattice + '_' + facet + '_slab'])
                    if not AB + 'star_' + A_B + '_star' in listdir(project):
                        print('Creating ' + path_reaction)
                        mkdir(path_reaction)
                    # CH4.write(path_reaction + '/CH4gas.traj')
                    # H2.write(path_reaction + '/H2gas.traj')
                    if not surface + '_' + facet in listdir(path_reaction):
                        print('Creating ' + path_surface)
                        mkdir(path_surface)
                    slab.write(path_surface + '/empty_surface.traj')
                    slab.write(path_surface + '/empty_surface.xyz')
                    atoms.write(path_surface + '/' + AB + '.traj')
                    atoms.write(path_surface + '/' + AB + '.xyz')
                    activation_energy = np.NaN
                    spreadsheet.append([surface, facet, AB, A, B,
                                        Ef, bee, activation_energy,
                                        'Quantum Espresso',
                                        'BEEF-vdW', 'M. H. Hansen, 2017',
                                        'suncat.stanford.edu'])
        with open(project+'test.csv', 'wb') as f:
            writer = csv.writer(f)
            writer.writerows(spreadsheet)
