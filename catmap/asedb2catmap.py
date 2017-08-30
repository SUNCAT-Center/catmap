# -*- coding: utf-8 -*-
"""
Description:
    Functions to convert ase.db files into a catmap input file

    File dependencies:
    ------------------
    <mol_name>.db: db file
        Contains molecular reference states.
        Mandatory key value pairs:
        --------------------------
                "enrgy" : float
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
            "adsorbate" : str
                value of the adsorbate chemical formula.
                '' should be the value for clean slabs.
                '-' should be inserted between seperate fragments.
            "enrgy" : float
                value with potential energy from DFT.

        Optional key value pairs:
        -------------------------
            "site" : str
                name of adsorption site.
            "data.BEEFens" : list
                32 non-selfconsistent BEEF-vdW energies.

    Optional file dependencies:
    ---------------------------
        .json files in the working directory containing a dictionary with
        frequencies. The file name and the key name should be
            adsorbate+'_'+facet
        # To do:
            Expand this feature to individual frequency lists for individual
            catalysts.
"""
import numpy as np
import ase.db
from ase.atoms import string2symbols
import json
from catmap.bee import BEEFEnsemble as BEE

state = BEE()


def get_refs(energy_dict, energy_mols, de_dict, de_mols):
    """ Returns dictionaries with referece energies of slabs and single atoms.
    """
    ref_dict = energy_mols
    ref_de = de_mols
    for key in energy_dict.keys():
        if 'slab' in key:
            ser, cat, pha, lattice, fac, site = key.split('_')
            name = cat + '_' + pha + '_' + lattice + '_' + fac + '_' + site
            ref_dict[name] = energy_dict[key]
            ref_de[name] = de_dict[key]
    return ref_dict, ref_de


def get_formation_energies(energy_dict, ref_dict):
    """ Returns a dictionary with formation energies of adsorbates.
    """
    formation_energies = {}
    for key in energy_dict.keys():
        E0 = 0
        if 'gas' in key:
            series, site_name = key.split('_')  # Split key into name/site
        else:
            try:
                series, cat, pha, lattice, fac, site = key.split('_')
            except ValueError as err:
                err.message += 'key=' + key
                raise
            site_name = cat + '_' + pha + '_' + lattice + '_' + fac + '_slab'
            try:
                E0 -= ref_dict[site_name]
            except KeyError:
                Warning('no slab reference ' + site_name)
                continue
        if 'slab' not in key:
            composition = string2symbols(series.replace('-', ''))
            E0 += energy_dict[key]
            for atom in composition:
                E0 -= ref_dict[atom]
            formation_energies[key] = round(E0, 4)
    return formation_energies


def get_BEEstd(de_dict, ref_de):
    """ Returns dictionary with BEEF standard deviations on formation energies.
    """
    BEEstd = {}
    for key in de_dict.keys():
        de = np.zeros(2000)
        if 'gas' in key:
            series, site = key.split('_')  # Split key into name/site
        else:
            series, cat, pha, lattice, fac, site = key.split('_')
            surface = cat + '_' + pha + '_' + lattice + '_' + fac + '_slab'
            try:
                de -= ref_de[surface]
            except KeyError:
                print('no slab BEE perturbation '+surface)
                continue
        if 'slab' not in key:
            composition = string2symbols(series.replace('-', ''))
            de += de_dict[key]
            for atom in composition:
                de -= ref_de[atom]
            BEEstd[key] = de.std()
    return BEEstd


def get_BEE_PCA(de_dict, ref_de, ads_x, ads_y, site_x, site_y):
    """ Returns two dictionaries, BEE_PC with the principal component vectors
    and BEE_lambda with the corresponding eigenvalues of BEE pertubations.

    Input:
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
    for slab in ref_de.keys():
        # Ignore gas species.
        if 'gas' in slab or len(slab) <= 2:
            continue
        de_x = -ref_de[slab]
        de_y = -ref_de[slab]
        try:
            de_x += de_dict[ads_x + '_' + slab[:-4] + site_x]
            de_y += de_dict[ads_y + '_' + slab[:-4] + site_y]
        except KeyError:
            print('Missing BEEF ensemble for ', ads_x, ads_y, slab[:-4],
                  site_x, site_y)
            continue
        composition_x = string2symbols(ads_x)
        composition_y = string2symbols(ads_y)
        for atom in composition_x:
            de_x -= ref_de[atom]
        for atom in composition_y:
            de_y -= ref_de[atom]
        cov = np.cov(de_x, de_y)
        eigval, eigvec = np.linalg.eig(cov)
        BEE_PC[slab] = eigvec
        BEE_lambda[slab] = eigval / len(cov)
    return BEE_PC, BEE_lambda


def make_input_file(file_name, energy_dict, frequency_dict={}, bee_dict={}):
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
    frequency_dict : dictionary
        Contains lists of frequencies. Each species is represented by a field
        and the fields must be named corresponding to energy_dict.
    bee_dict : dictionary
        Contains standard deviations on the energy.
        Each species is represented by a field and the fields must be named
        corresponding to energy_dict.
    """
    # Create a header
    header = '\t'.join(['surface_name', 'phase', 'site_name',
                        'species_name', 'formation_energy',
                        'frequencies', 'reference', 'bee'])
    lines = []  # List of lines in the output
    for key in energy_dict.keys():  # Iterate through keys
        E = energy_dict[key]  # Ab initio energy
        if 'gas' in key:
            name, site = key.split('_')  # Split key into name/site
            try:
                frequency = frequency_dict[key]
            except KeyError:
                frequency = []
            try:
                bee = bee_dict[key]
            except KeyError:
                bee = np.NaN
        else:
            name, cat, pha, lattice, facet, site = key.split('_')
        if 'slab' not in key:  # Do not include empty site energy (0)
            try:
                frequency = frequency_dict[key]
            except KeyError:
                frequency = []
            try:
                bee = bee_dict[key]
            except KeyError:
                bee = np.NaN
            if site == 'gas':
                surface = None
                phase = ''
                lattice = ''
                site_name = 'gas'
            else:
                surface = cat  # +'_'+pha #useful to include phase.
                phase = pha
                site_name = lattice + '_' + site
            outline = [surface, phase, site_name, name,
                       E, frequency, 'MHH_DFT', bee]
            line = '\t'.join([str(w) for w in outline])
            lines.append(line)
    lines.sort()  # The file is easier to read if sorted (optional)
    lines = [header] + lines  # Add header to top
    input_file = '\n'.join(lines)  # Join the lines with a line break
    input = open(file_name, 'w')  # Open the file name in write mode
    input.write(input_file)  # Write the text
    input.close()  # Close the file


def db2mol(fname, selection=[], freq_path='.'):  # fname for db with molecules.
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
    frequency_dict = {}
    dbids = {}
    de = {}
    for d in smol:              # Get molecules from mol.db
        abinitio_energy = float(d.enrgy)
        species_name = str(d.formula)
        if species_name+'_gas' not in abinitio_energies:
            abinitio_energies[species_name+'_gas'] = abinitio_energy
            dbids[species_name+'_gas'] = int(d.id)
            de[species_name+'_gas'] = \
                state.get_ensemble_perturbations(d.data.BEEFens)
            try:
                freq = json.load(open(freq_path + '/' + species_name +
                                      '_gas.freq', 'r'))
                frequency_dict.update(freq)
            except IOError:
                print('no frequencies for', species_name, '(g)')
        elif abinitio_energies[species_name+'_gas'] > abinitio_energy:
            abinitio_energies[species_name+'_gas'] = abinitio_energy
            dbids[species_name+'_gas'] = int(d.id)
            de[species_name+'_gas'] = \
                state.get_ensemble_perturbations(d.data.BEEFens)
    return abinitio_energies, frequency_dict, de, dbids


def mol2ref(abinitio_energies, de={}):  # , references={'H': 'H2_g',
                                        #          'O': 'H2O_g',
                                        #          'C': 'CH4_g', }):
    """ Returns two dictionaries containing:
        abinitio energy references for atoms
        32 non-selfconsitent perturbations for atoms.
    """
    mol_e = {}
    mol_de = {}
    mol_e['H'] = 0.5 * abinitio_energies['H2_gas']
    mol_e['O'] = abinitio_energies['H2O_gas'] - 2 * mol_e['H']
    mol_e['C'] = abinitio_energies['CH4_gas'] - 4 * mol_e['H']
    try:
        mol_de['H'] = 0.5 * de['H2_gas']
        mol_de['O'] = de['H2O_gas'] - 2*mol_de['H']
        mol_de['C'] = de['CH4_gas'] - 4*mol_de['H']
    except KeyError as err:
        err.message += ' Missing BEE perturbations for molecules!'
        raise
    return mol_e, mol_de


def db2surf(fname, selection=[], sites=False):
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
    frequency_dict = {}
    dbids = {}
    de = {}
    for d in ssurf:                     # Get slab and adsorbates from .db
        series = str(d.series)
        adsorbate = str(d.adsorbate)
        if 'FBL' in series or 'NEB' in series or '-' in adsorbate:
            continue
        if 'n' in d and int(d.n) > 1:  # Skip higher coverages for now.
            continue
        cat = str(d.name)+'_'+str(d.phase)
        abinitio_energy = float(d.enrgy)
        surf_lattice = str(d.surf_lattice)
        # composition=str(d.formula)
        facet = str(d.facet)
        if adsorbate == '' or series == 'slab':
            adsorbate = ''
            site = 'slab'
        elif sites:
            site = str(d.site)
        else:
            site = 'site'
        site_name = surf_lattice + '_' + facet + '_' + site
        key = adsorbate + '_' + cat + '_' + site_name
        if key not in abinitio_energies:
            abinitio_energies[key] = abinitio_energy
            dbids[key] = int(d.id)
            de[key] = state.get_ensemble_perturbations(d.data.BEEFens)
            if not series == 'slab':
                try:
                    freq = json.load(open(adsorbate + '_' + surf_lattice +
                                          '.freq', 'r'))
                    frequency_dict.update({adsorbate + '_' + cat + '_' +
                                           site_name: freq[adsorbate +
                                                           '_' + surf_lattice]}
                                          )
                except IOError:
                    print('no frequencies for', adsorbate + '_' + site_name)
        elif abinitio_energies[key] > abinitio_energy:
            abinitio_energies[key] = abinitio_energy
            dbids[key] = int(d.id)
            de[key] = state.get_ensemble_perturbations(d.data.BEEFens)
    return abinitio_energies, frequency_dict, de, dbids


def db2pes(fname, selection=[]):
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
    surfaces = {}
    # Get images from ase .db
    for d in s:
        adsorbate = str(d.adsorbate)
        if '-' not in adsorbate:
            continue
        surf = str(d.name) + '_' + str(d.phase) + '_' + \
            str(d.surf_lattice) + '_' + str(d.facet)
        abinitio_energy = float(d.enrgy)
        dbid = int(d.id)
        try:
            BEEFens = d.data.BEEFens
            de = state.get_ensemble_perturbations(BEEFens)
        except AttributeError:
            print(dbid, adsorbate, surf, 'missing BEEF perturbations.')
            de = np.zeros(2000)
        if surf in surfaces:
            if adsorbate in surfaces[surf]:
                surfaces[surf][adsorbate]['pes'].append(abinitio_energy)
                surfaces[surf][adsorbate]['dbids'].append(dbid)
                surfaces[surf][adsorbate]['de'].append(de)
            else:
                surfaces[surf][adsorbate] = {'pes': [abinitio_energy],
                                             'dbids': [dbid],
                                             'de': [de],
                                             'images': []}
        else:
            surfaces[surf] = {adsorbate: {'pes': [abinitio_energy],
                                          'dbids': [dbid],
                                          'de': [de],
                                          'images': []}}
        if 'image' in d:
            surfaces[surf][adsorbate]['images'].append(int(d.image))
        else:
            surfaces[surf][adsorbate]['images'].append(int(d.step))
    return surfaces


def pes2ts(surfaces):
    """ Returns dictionaries containing transition state ab initio energies.

        Input parameters
        ----------------
        surfaces : dictionary
            Created by the db2pes function.
    """
    abinitio_energies = {}
    frequency_dict = {}
    dbids = {}
    de = {}
    for m in surfaces:
        for adsorbate in surfaces[m]:
            key = adsorbate + '_' + m + '_site'
            images = surfaces[m][adsorbate]['images']
            if len(np.unique(images)) != len(images):
                print('non unique image number!', m, adsorbate, images)
                break
            pes = surfaces[m][adsorbate]['pes']
            if len(pes) < 4:
                continue
            #s = np.argsort(images)
            # Look for local minima and maxima.
            #localmins = np.where(np.r_[True, pes[1:] < pes[:-1]] &
            #                    np.r_[pes[:-1] < pes[1:], True])[0]
            #localmaxs = np.where(np.r_[True, pes[1:] > pes[:-1]] &
            #                    np.r_[pes[:-1] > pes[1:], True])[0]
            #if len(localmaxs) > 1:
            #   print(len(localmaxs), 'local maxima in', key)
            #if len(localmins) > 2:
            #   print(len(localmins), 'local minima in', key)
            tst = np.argmax(pes)
            if key not in abinitio_energies:
                abinitio_energies[key] = pes[tst]
                dbids[key] = surfaces[m][adsorbate]['dbids'][tst]
                de[key] = surfaces[m][adsorbate]['de'][tst]
                try:
                    freq = json.load(open(adsorbate + '_' + m +
                                          '.freq', 'r'))
                    frequency_dict.update({adsorbate + '_' + m + '_' +
                                           key: freq[adsorbate + '_' + m]})
                except IOError:
                    print('no frequencies for', key)
            elif abinitio_energies[key] > pes[tst]:
                abinitio_energies[key] = pes[tst]
                dbids[key] = surfaces[m][adsorbate]['dbids'][tst]
                de[key] = surfaces[m][adsorbate]['de'][tst]
    return abinitio_energies, frequency_dict, de, dbids
