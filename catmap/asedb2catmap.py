# -*- coding: utf-8 -*-
"""

Description:
    Functions to convert ase.db files into a catmap input file

Dependencies:
    mol.db: db file which must include molecular reference states and
        the key-value fields:
        "enrgy" with float value of DFT calculated potential energy.
        "vacuum" with float value specifying amount of vacuum.
        "PW" float with plane wave cutoff.

    <name>.db: db file, which must include adsorbate states and
        the key-value fields:
        "name" with str value that identifies the catalyst composition.
        "phase" with str value that identifies the catalyst phase.
        "facet" with str value that identifies the site.
            Indices should be separated by 'x', e.g. '1x1x1'
        # To do:
            Distinguish facet from site by adding a new required key.
        "adsorbate" with string value of the adsorbate chemical formula.
            '' should be the value for clean slabs.
        "enrgy" with float value with potential energy from DFT.
        "PW" float with plane wave cutoff.
        "kpts" str with k-points in the plane dimensions separated by 'x',
            e.g. '4x4'
        Optional keys:
            "data.BEEFens" containing a list of 32 non-selfconsistent BEEF-vdW
            energies.

    Optional dependencies:
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
    ref_dict = energy_mols
    ref_de = de_mols
    for key in energy_dict.keys():
        if 'slab' in key:
            ser, cat, pha, fac, lattice, site = key.split('_')
            name = cat + '_' + pha + '_' + fac + '_' + lattice + '_' + site
            ref_dict[name] = energy_dict[key]
            ref_de[name] = de_dict[key]
    return ref_dict, ref_de


def get_formation_energies(energy_dict, ref_dict):
    formation_energies = {}
    for key in energy_dict.keys():
        E0 = 0
        if 'gas' in key:
            series, site = key.split('_')  # Split key into name/site
        else:
            try:
                series, cat, pha, fac, lattice, site = key.split('_')
            except ValueError as err:
                err.message += 'key=' + key
                raise
            surface = cat + '_' + pha + '_' + fac + '_' + lattice + '_slab'
            try:
                E0 -= ref_dict[surface]
            except KeyError:
                Warning('no slab reference ' + surface)
                continue
        if 'slab' not in key:
            try:
                composition = string2symbols(series)
            except ValueError:
                series = series[:-2]
                composition = string2symbols(series)
            E0 += energy_dict[key]
            for atom in composition:
                E0 -= ref_dict[atom]
            formation_energies[key] = round(E0, 4)
    return formation_energies


def get_BEEstd(de_dict, ref_de):
    BEEstd = {}
    for key in de_dict.keys():
        de = np.zeros(2000)
        if 'gas' in key:
            series, site = key.split('_')  # Split key into name/site
        else:
            series, cat, pha, fac, lattice, site = key.split('_')
            surface = cat + '_' + pha + '_' + fac + '_' + lattice + '_slab'
            try:
                de -= ref_de[surface]
            except KeyError:
                print('no slab BEE perturbation '+surface)
                continue
        if 'slab' not in key:
            try:
                composition = string2symbols(series)
            except ValueError:
                series = series[:-2]
                composition = string2symbols(series)
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
            print('Found BEEF ensemble for ', ads_x, ads_y, slab[:-4],
                  site_x, site_y)
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
            name, cat, pha, fac, lattice, site = key.split('_')
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
                else:
                    surface = cat  # +'_'+pha #useful to include phase.
                    phase = pha
                outline = [surface, phase, lattice + '_' + site, name,
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


def mol2ref(abinitio_energies, de={}):
    mol_e = {}
    mol_de = {}
    mol_e['H'] = 0.5 * abinitio_energies['H2_gas']
    mol_e['O'] = abinitio_energies['H2O_gas'] - 2 * mol_e['H']
    mol_e['C'] = abinitio_energies['CH4_gas'] - 4 * mol_e['H']
    # mol_dict['C'] = abinitio_energies['CO_gas'] - mol_dict['O']
    try:
        mol_de['H'] = 0.5 * de['H2_gas']
        mol_de['O'] = de['H2O_gas'] - 2*mol_de['H']
        mol_de['C'] = de['CH4_gas'] - 4*mol_de['H']
    except KeyError as err:
        err.message += ' Missing BEE perturbations for molecules!'
        raise
    return mol_e, mol_de


def db2surf(fname, selection=[]):
    csurf = ase.db.connect(fname)
    ssurf = csurf.select(selection)
    abinitio_energies = {}
    frequency_dict = {}
    dbids = {}
    de = {}
    for d in ssurf:                     # Get slab and adsorbates from .db
        series = str(d.series)
        if 'FBL' in series or 'NEB' in series:
            continue
        if 'n' in d and int(d.n) > 1:   # Skip higher coverages
            continue
        if '-' in series:
            continue
        cat = str(d.name)+'_'+str(d.phase)
        abinitio_energy = float(d.enrgy)
        surf_lattice = str(d.surf_lattice)
        # composition=str(d.formula)
        adsorbate = str(d.adsorbate)
        facet = str(d.facet)
        if adsorbate == '' or series == 'slab':
            adsorbate = ''
            site = 'slab'
        else:
            site = str(d.site)
        site_name = facet + '_' + surf_lattice + '_' + site
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
