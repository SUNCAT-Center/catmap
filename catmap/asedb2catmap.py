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

def get_refs(energy_dict,energy_mols,de_dict,de_mols): #adapted from CATMAP wiki
    ref_dict = energy_mols
    ref_de = de_mols
    for key in energy_dict.keys():
        if 'slab' in key:
            ser,cat,pha,fac = key.split('_')
            name = cat+'_'+pha+'_'+fac
            ref_dict[name] = energy_dict[key]
            ref_de[name] = de_dict[key]
    return ref_dict, ref_de

def get_formation_energies(energy_dict,ref_dict): #adapted from CATMAP wiki
    formation_energies = {}
    for key in energy_dict.keys():
        E0 = 0
        if 'gas' in key:
            name,site = key.split('_') #split key into name/site
        else:
            try:
                name,cat,pha,fac = key.split('_')
            except ValueError as err:
                err.message += 'key='+key
                raise
            site = cat+'_'+pha+'_'+fac
            try:
                E0 -= ref_dict[site]
            except KeyError:
                Warning('no slab reference '+site)
                continue
        if name != 'slab':
            try:
                composition = string2symbols(name)
            except ValueError:
                name = name[:-2]
                composition = string2symbols(name)
            E0 += energy_dict[key]
            for atom in composition:
                E0 -= ref_dict[atom]
            formation_energies[key] = round(E0, 4)
    return formation_energies

def get_formation_energies_surftype(energy_dict,ref_dict): #adapted from CATMAP wiki
    formation_energies = {}
    for key in energy_dict.keys():
        E0 = 0
        if 'gas' in key:
            name,site = key.split('_') #split key into name/site
        else:
            name,cat,pha,fac = key.split('_')
            site = cat+'_'+pha+'_'+fac
            try:
                E0 -= ref_dict[site]
            except KeyError:
                print('no slab reference '+site)
                continue
        if name != 'slab':
            try:
                composition = string2symbols(name)
            except ValueError:
                name = name.replace('-','')
                composition = string2symbols(name)
            E0 += energy_dict[key]
            for atom in composition:
                E0 -= ref_dict[atom]
            formation_energies[key] = round(E0, 4)
    return formation_energies

def get_BEEstd(de_dict, ref_de): #adapted from CATMAP wiki and GPAW wiki.
    BEEstd = {}
    for key in de_dict.keys():
        de = np.zeros(2000)
        if 'gas' in key:
            name,site = key.split('_') #split key into name/site
        else:
            name,cat,pha,fac = key.split('_')
            site = cat+'_'+pha+'_'+fac
            try:
                de -= ref_de[site]
            except KeyError:
                print('no slab BEE perturbation '+site)
                continue
        if name != 'slab':
            try:
                composition = string2symbols(name)
            except ValueError:
                name = name[:-2]
                composition = string2symbols(name)
            de += de_dict[key]
            for atom in composition:
                de -= ref_de[atom]
            BEEstd[key] = de.std()
    return BEEstd

def get_BEE_PCA(de_dict, ref_de, ads_x, ads_y): #adapted from CATMAP wiki and GPAW wiki.
    """ Returns two dictionaries, BEE_PC with the principal component vectors 
    and BEE_lambda with the corresponding eigenvalues of BEE pertubations.
    
    Input:
        de_dict  dict    contains abinito energies of adsorbates on slabs,
                            where keys are named: adsorbate_name_phase_facet
        ref_de   dict    contains abinitio energies of references,
                            where keys are refernce elements, e.g: 'C','H',
                            and also slabs references: name_phase_facet
        ads_x    string  adsorbate first dimension
        ads_y    string  adsorbate second dimension
    
    """
    BEE_PC = {}
    BEE_lambda = {}
    for site in ref_de.keys():
        if 'gas' in site or len(site) <= 2:
            continue
        de_x =-ref_de[site]
        de_y =-ref_de[site]
        try:
            de_x += de_dict[ads_x + '_' + site]
            de_y += de_dict[ads_y + '_' + site]
        except KeyError:
            print('Missing BEEF ensemble for '+site)
            continue
        composition_x = string2symbols(ads_x)
        composition_y = string2symbols(ads_y)
        for atom in composition_x:
            de_x -= ref_de[atom]
        for atom in composition_y:
            de_y -= ref_de[atom]
        cov = np.cov(de_x,de_y)
        eigval, eigvec = np.linalg.eig(cov)
        BEE_PC[site] = eigvec
        BEE_lambda[site] = eigval
    return BEE_PC, BEE_lambda

def make_input_file(file_name,energy_dict,frequency_dict={},bee_dict={}): #adapted from CATMAP
    #create a header
    header = '\t'.join(['surface_name','site_name',
                        'species_name','formation_energy',
                        'frequencies','reference','bee'])
    lines = [] #list of lines in the output
    for key in energy_dict.keys(): #iterate through keys
        E = energy_dict[key] #raw energy
        if 'gas' in key:
            name,site = key.split('_') #split key into name/site
            try:
                frequency = frequency_dict[key]
            except KeyError:
                frequency = []
            try:
                bee = bee_dict[key]
            except KeyError:
                bee = np.NaN
        else:
            name,cat,pha,facet = key.split('_') #split key into name/site
            if (facet == '1x1x1' and pha == 'fcc') or (facet == '0x0x1' and pha == 'hcp'):
                site = 'hexagonal'
            else:
                site = facet
        if 'slab' not in name: #do not include empty site energy (0)          
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
            else:
                surface = cat #+'_'+pha #useful to include phase in some cases.
            outline = [surface,site,name,E,frequency,'MHH_DFT',bee]
            line = '\t'.join([str(w) for w in outline])
            lines.append(line)
    lines.sort() #The file is easier to read if sorted (optional)
    lines = [header] + lines #add header to top
    input_file = '\n'.join(lines) #Join the lines with a line break
    input = open(file_name,'w') #open the file name in write mode
    input.write(input_file) #write the text
    input.close() #close the file

def db2mol(fname,selection=[],freq_path='.'): #fname must be path/filename of db containing molecules
    cmol = ase.db.connect(fname)
    smol = cmol.select(selection)
    #mol_dict = {}
    abinitio_energies = {}
    frequency_dict = {}
    dbids = {}
    de = {}
    for d in smol:              #get molecules from mol.db
        #cat = 'None'
        #site_name = 'gas'
        abinitio_energy = float(d.enrgy)
        species_name=str(d.formula)
        if species_name+'_gas' not in abinitio_energies:
            abinitio_energies[species_name+'_gas'] = abinitio_energy
            dbids[species_name+'_gas'] = int(d.id)
            de[species_name+'_gas'] = state.get_ensemble_perturbations(d.data.BEEFens)
            try:
                freq = json.load(open(freq_path+'/'+species_name+'_gas.freq','r'))
                frequency_dict.update(freq)
            except IOError:
                print('no frequencies for',species_name,'(g)')
        elif abinitio_energies[species_name+'_gas'] > abinitio_energy:
            abinitio_energies[species_name+'_gas'] = abinitio_energy
            dbids[species_name+'_gas'] = int(d.id)
            de[species_name+'_gas'] = state.get_ensemble_perturbations(d.data.BEEFens)
    return abinitio_energies, frequency_dict, de, dbids

def mol2ref(abinitio_energies, de={}):
    mol_e = {}
    mol_de = {}
    mol_e['H'] = 0.5*abinitio_energies['H2_gas']
    mol_e['O'] = abinitio_energies['H2O_gas'] - 2*mol_e['H']
    mol_e['C'] = abinitio_energies['CH4_gas'] - 4*mol_e['H']
    #mol_dict['C'] = abinitio_energies['CO_gas'] - mol_dict['O']
    try:
        mol_de['H'] = 0.5*de['H2_gas']
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
    for d in ssurf:                     #get slab and adsorbates from .db
        series = str(d.series)
        if 'FBL' in series or 'NEB' in series:
            continue
        if 'n' in d and int(d.n) > 1:   #skip higher coverages
            continue
        if '-' in series:
            continue
        cat = str(d.name)+'_'+str(d.phase)
        site_name = str(d.facet)
        abinitio_energy = float(d.enrgy)
        surf_lattice = str(d.surf_lattice)
        #composition=str(d.formula)
        adsorbate = str(d.adsorbate)
        if adsorbate == '' and series == 'slab':
            adsorbate = 'slab'
        if adsorbate+'_'+cat+'_'+site_name not in abinitio_energies:
            abinitio_energies[adsorbate+'_'+cat+'_'+site_name] = abinitio_energy
            dbids[adsorbate+'_'+cat+'_'+site_name] = int(d.id)
            de[adsorbate+'_'+cat+'_'+site_name] = state.get_ensemble_perturbations(d.data.BEEFens)
            if not series == 'slab':
                try:
                    freq = json.load(open(adsorbate+'_'+surf_lattice+'.freq','r'))
                    frequency_dict.update({adsorbate+'_'+cat+'_'+site_name: freq[adsorbate+'_'+surf_lattice]})
                except IOError:
                    print('no frequencies for',adsorbate+'_'+site_name)
        elif abinitio_energies[adsorbate+'_'+cat+'_'+site_name] > abinitio_energy:
            abinitio_energies[adsorbate+'_'+cat+'_'+site_name] = abinitio_energy
            dbids[adsorbate+'_'+cat+'_'+site_name] = int(d.id)
            de[adsorbate+'_'+cat+'_'+site_name] = state.get_ensemble_perturbations(d.data.BEEFens)
    return abinitio_energies, frequency_dict, de, dbids






















































