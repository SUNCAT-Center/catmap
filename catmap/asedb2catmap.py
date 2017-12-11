# -*- coding: utf-8 -*-
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
    <fname>.db : db file
        Stores vibrational frequencies along with atomic structures
        and energies.
"""
from os import listdir, environ, mkdir
import os.path
import numpy as np
import ase.db
from ase.atoms import string2symbols
from ase.data import covalent_radii, atomic_numbers
from ase.calculators.singlepoint import SinglePointDFTCalculator
from catmap.bee import BEEFEnsemble as BEE
import csv


class db2catmap(object):
    def __init__(self, mol_db=None, ads_db=None, mol_select=[], ads_select=[],
                 frequency_db=None, frequency_select=[],
                 ts_db=None, ts_selection=[]):
        self.state = BEE()
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
         mol_ens,
         mol_dbid] = self._db2mol(fname, selection=selection,
                                  freq_path=frequency_db)
        self.epot.update(mol_epot)
        self.freq.update(mol_freq)
        self.ens.update(mol_ens)
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
         surf_ens,
         surf_dbid] = self._db2surf(fname, selection=selection,
                                    freq_path=frequency_db,
                                    site_specific=site_specific)
        self.epot.update(surf_epot)
        self.freq.update(surf_freq)
        self.ens.update(surf_ens)
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
         surf_ens,
         surf_dbid] = self.pes2ts(freq_path=frequency_db)
        self.epot.update(surf_epot)
        self.freq.update(surf_freq)
        self.ens.update(surf_ens)
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
        # Get atomic references.
        [self.atomic_e,
         self.atomic_ens] = self._mol2ref(references=references)
        # Get dictionaries with slab references and atomic references.
        [self.reference_epot,
         self.reference_ens] = self._get_refs()
        #
        self.de_dict, self.std = self._get_BEEstd()
        self.formation_energies = self._get_formation_energies()

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
        ens_dict = {}
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
                ens_dict[species_name+'_gas'] = ens
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
                ens_dict[species_name+'_gas'] = contribs
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
        ens_dict = {}
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
            elif 'site' in d and site_specific is True:
                site = str(d.site)
            else:
                site = 'site'
            cat = name + '_' + phase
            site_name = '_'.join([surf_lattice, facet, cell, site])
            key = '_'.join([str(n), species, cat, site_name])
            freq_key = key  # str(n) + '_' + species + '_' + surf_lattice
            try:
                contribs = d.data.BEEFvdW_contribs
                ens = self.state.get_ensemble_perturbations(contribs)
            except:
                ens = 0  # np.zeros(self.state.size)
            if key not in abinitio_energies:
                abinitio_energies[key] = abinitio_energy
                dbids[key] = int(d.id)
                ens_dict[key] = ens
                if species != '' and ads != 'slab' and freq_path is not None:
                    try:
                        freqsearch = ['species='+species, 'name='+name]
                        if site_specific is True:
                            freqsearch.append('site='+site)
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
            surf = '_'.join([str(d.name), phase, surf_lattice, facet, cell])
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
                rxn_paths[rxn_id]['ens'].append(ens)
                # try:
                rxn_paths[rxn_id]['distance'].append(float(d.distance))
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
                if 'site' in d:
                    site = str(d.site)
                else:
                    site = 'site'
                rxn_paths[rxn_id] = {'surface_name': surf,
                                     'species': species,
                                     'pes': [abinitio_energy],
                                     'dbids': [dbid],
                                     'ens': [ens],
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
        atomic_ens = {}
        print(references)
        for t in references:
            key = t[0]
            species = t[1]
            atomic_e[key] = self.epot[species]
            atomic_ens[key] = np.array(self.ens[species])
            composition = string2symbols(species.split('_')[0])
            n = 0.
            for symbol in composition:
                if symbol == key:
                    n += 1
                else:
                    atomic_e[key] -= atomic_e[symbol]
                    atomic_ens[key] -= np.array(atomic_ens[symbol])
            atomic_e[key] /= n
            atomic_ens[key] /= n
        return atomic_e, atomic_ens

    def _get_refs(self):
        """ Returns dictionaries with referece energies of slabs and
        single atoms.
        """
        ref_dict = self.atomic_e
        ref_ens = self.atomic_ens
        for key in self.epot.keys():
            if 'slab' in key:
                ref_dict[key] = self.epot[key]
                ref_ens[key] = self.ens[key]
        return ref_dict, ref_ens

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
                site_name = '_'.join(['0_', cat, pha, lattice, fac, cell,
                                      'slab'])
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
        return formation_energies

    def _get_BEEstd(self):
        """ Returns a dictionary with BEEF ensembles and with
        BEEF standard deviations on formation energies.
        """
        de_dict = {}
        std_dict = {}
        for key in self.ens:
            de = np.zeros(self.state.size)
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
                std_dict[key] = np.std(de)
        return de_dict, std_dict

    def get_ellipses(self, ads_x, ads_y,
                     site_x=None, site_y=None):
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
        widths = {}
        heights = {}
        angles = {}
        # Loop over reference surfaces in ref_de.
        for slab in self.reference_epot.keys():
            # Ignore gas species.
            if 'slab' not in slab:
                continue
            key_x = self._insert_species(slab, ads_x, site_x)
            if isinstance(key_x, list):
                de_x = np.zeros(self.state.size)
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
            key_y = self._insert_species(slab, ads_y, site_y)
            if key_y in self.de_dict:
                de_y = self.de_dict[key_y]
            else:
                continue
            if np.isclose(de_x, 0).all() or np.isclose(de_y, 0).all():
                continue
            cov = np.cov(de_x, de_y)
            eigval, eigvec = np.linalg.eig(cov)
            width, height = 2*np.sqrt(eigval)
            angle = np.rad2deg(np.arccos(eigvec[0, 0]))
            widths[slab] = width
            heights[slab] = height
            angles[slab] = angle
        return widths, heights, angles

    def pes2ts(self, freq_path=None, rtol=1.03):
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
            dbond = covalent_radii[atomic_numbers[g1[0]]] + \
                covalent_radii[atomic_numbers[g2[0]]]
            if len(np.unique(images)) != len(images):
                warn = True
                print('non unique image number!')
                print('Warning!', species, m, roughness,
                      len(localmaxs), len(localmins), images)
                continue
            if len(localmaxs) > 1 or len(localmins) > 2 or len(localmins) == 1:
                warn = True
            shortest = np.min(self.rxn_paths[rxn_id]['distance'])
            if shortest > dbond * rtol:
                warn = True
                assert np.argmax(self.rxn_paths[rxn_id]['dbids']) == \
                    np.argmax(self.rxn_paths[rxn_id]['images'])
                calculate.append(max(self.rxn_paths[rxn_id]['dbids']))
                continue
            tst = np.argmax(pes)
            if warn:
                print('Warning!', species, m, round(dbond * rtol, 3),
                      round(shortest, 3), roughness,
                      len(localmaxs), len(localmins), len(images))
            if key not in abinitio_energies:
                abinitio_energies[key] = pes[tst]
                dbids[key] = self.rxn_paths[rxn_id]['dbids'][tst]
                ens_dict[key] = self.rxn_paths[rxn_id]['ens'][tst]
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
                ens_dict[key] = self.rxn_paths[rxn_id]['ens'][tst]
        incomplete = ','.join([str(int(a)) for a in np.unique(calculate)])
        print('Incomplete:', incomplete)
        return abinitio_energies, freq_dict, ens_dict, dbids

    def scaling_analysis(self, x, y, lattice=None, site_x=None, site_y=None):
        """Returns the scaling relation information between the species
        x and y on the surface geometry 'lattice' and site.
        The outputs are in order: slope, intercept, which are floats,
        followed by X and Y, which are lists of formation energies.

        Parameters
        ----------
        x : str or list
            species x
        y : str
            species y
        lattice : str or None
            surface lattice of y
        site='site' : str
            Site of y. If _db2surf was run with site_specific=False, use the
            default value, 'site'.
        """
        X = []
        Y = []
        l = []
        # Loop over slabs.
        for slab in self.reference_epot.keys():
            if 'slab' not in slab:
                continue
            # Filter by surface lattice.
            if lattice is not None and lattice not in slab:
                continue
            # Get formation energy if it is stored.
            key_y = self._insert_species(slab, y, site_y)
            if key_y not in self.formation_energies:
                continue
            key_x = self._insert_species(slab, x, site_x)
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
            l.append(slab)
        # Get the scaling relation.
        slope, intercept = np.polyfit(X, Y, deg=1)
        return slope, intercept, X, Y, l

    def _insert_species(self, slab, species, site=None):
        """Return a key for a hypothetical adsorbate state. This is useful
        for reading or filling in formation energies of states related
        by scaling relations.

        Parameters
        ----------
        slab : str
            key of the slab
        species : str
            adsorbate formula
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

    def insert_interpolated_states(self, x, y, site_y=None,
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
        for key_x in x:
            X = self.formation_energies[key_x]
            Y = slope * X + intercept
            fields = key_x.split('_')
            fields[1] = y
            if site_y is not None:
                fields[-1] = site_y
            key_y = '_'.join(fields)
            self.formation_energies.update({key_y: Y})

    def make_input_file(self, file_name, site_specific=False,
                        covariance=None):
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
                    coverage = 0
                else:
                    surface = cat
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
                if covariance is not None and covariance[1] == name:
                    slab = '_'.join(['0_', cat, pha, lattice, facet, cell,
                                     'slab'])
                    if slab in self.width:
                        width = round(self.width[slab], 4)
                        height = round(self.height[slab], 4)
                        angle = round(self.angle[slab], 4)
                        reference = covariance[0]
                    else:
                        width = ''
                        height = ''
                        angle = ''
                        reference = ''
                    outline += [width, height, angle, reference]
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
                            site='site', mol_db=None,
                            slab_db=None, ads_db=None, ts_db=None,
                            publication='', url=''):
        """Saves a nested directory structure,
        compatible with the catapp project.

        Parameters
        ----------
        project : str
            parent folder name.
        reactions : list
            catmap's rxn_expressions. A list of strings.
        """
        # Create a header
        spreadsheet = [['chemical_composition', 'facet', 'reactants',
                        'products', 'reaction_energy',
                        'beef_standard_deviation',
                        'activation_energy', 'DFT_code', 'DFT_functional',
                        'reference', 'url']]
        # width, height, angle, covariance
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
                facet_name = facet.replace('(', '').replace(')', '')
                path_facet = path_surface + '/' + facet_name
                DeltaE = 0.
                de = np.zeros(self.state.size)
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
                        rkey = '_'.join(['1', species.replace('*', ''), name,
                                         phase, lattice, facet, cell, site])
                        fname = path_facet + '/' + species
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
                    for ts in tstates:
                        if '-' not in ts:
                            continue
                        species, sitesymbol = ts.split('_')
                        tskey = '_'.join(['0', species, name, phase, lattice,
                                          facet, cell, site])
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
                        pkey = '_'.join(['1', species.replace('*', ''), name,
                                         phase, lattice, facet, cell, site])
                        fname = path_facet + '/' + species
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
                    if not os.path.isdir(project):
                        mkdir(project)
                    if reaction_name not in listdir(project):
                        mkdir(path_reaction)
                    if name.replace('/', '') not in listdir(path_reaction):
                        mkdir(path_surface)
                    if facet_name not in listdir(path_surface):
                        mkdir(path_facet)
                    for trajkey in totraj.keys():
                        fname = totraj[trajkey]['fname'] + '.traj'
                        if 'gas' in trajkey:
                            if fname in os.listdir(path_reaction):
                                continue
                            else:
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
                            calc = SinglePointDFTCalculator(atoms)
                            calc.results['energy'] = float(d.epot)
                            atoms.set_calculator(calc)
                        fname = fname.replace('(', '').replace(')', '')
                        atoms.write(fname)
                    std = np.std(de)
                    spreadsheet.append([name, facet, rname, pname,
                                        DeltaE, std, ea,
                                        'Quantum Espresso',
                                        'BEEF-vdW', publication, url])
                    # width, height, angle, covariance
                    Nsurf += 1
                else:
                    continue
            Nrxn += 1
        with open(project + '/data.csv', 'wb') as f:
            writer = csv.writer(f)
            writer.writerows(spreadsheet)
        print(Nrxn, 'reactions imported.')
        print(Nsurf, 'surfaces saved.')
