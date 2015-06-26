from ase.atoms import string2symbols

abinitio_energies = {
        'CO_gas': -626.611970497, 
        'H2_gas': -32.9625308725, 
        'CH4_gas': -231.60983421,
        'H2O_gas': -496.411394229, 
        'CO_111': -115390.445596, 
        'C_111': -114926.212205,
        'O_111': -115225.106527,
        'H_111': -114779.038569,
        'CH_111': -114943.455431,
        'OH_111': -115241.861661,
        'CH2_111': -114959.776961,
        'CH3_111': -114976.7397,
        'C-O_111': -115386.76440668429,
        'H-OH_111': -115257.78796158083,
        'H-C_111': -114942.25042955727,
        'slab_111': -114762.254842,
        }

ref_dict = {}
ref_dict['H'] = 0.5*abinitio_energies['H2_gas']
ref_dict['O'] = abinitio_energies['H2O_gas'] - 2*ref_dict['H']
ref_dict['C'] = abinitio_energies['CH4_gas'] - 4*ref_dict['H']
ref_dict['111'] = abinitio_energies['slab_111']

def get_formation_energies(energy_dict,ref_dict):
    formation_energies = {}
    for key in energy_dict.keys(): #iterate through keys
        E0 = energy_dict[key] #raw energy
        name,site = key.split('_') #split key into name/site
        if 'slab' not in name: #do not include empty site energy (0)
            if site == '111':
                E0 -= ref_dict[site] #subtract slab energy if adsorbed
            #remove - from transition-states
            formula = name.replace('-','')
            #get the composition as a list of atomic species
            composition = string2symbols(formula)
            #for each atomic species, subtract off the reference energy
            for atom in composition:
                E0 -= ref_dict[atom]
            #round to 3 decimals since this is the accuracy of DFT
            E0 = round(E0,3)
            formation_energies[key] = E0
    return formation_energies

formation_energies = get_formation_energies(abinitio_energies,ref_dict)

for key in formation_energies:
    print key, formation_energies[key]

frequency_dict = {
                'CO_gas': [2170],
                'H2_gas': [4401], 
                'CH4_gas':[2917,1534,1534,3019,3019,3019,1306,
                            1306,1306],
                'H2O_gas': [3657, 1595, 3756],
                'CO_111': [60.8, 230.9, 256.0, 302.9, 469.9, 1747.3],
                'C_111': [464.9, 490.0, 535.9],
                'O_111': [359.5, 393.3, 507.0],
                'H_111': [462.8, 715.9, 982.5],
                'CH_111': [413.3, 437.5, 487.6, 709.6, 735.1, 3045.0],
                'OH_111': [55, 340.9, 396.1, 670.3, 718.0, 3681.7],
                'CH2_111': [55, 305.5, 381.3, 468.0, 663.4, 790.2, 1356.1, 
                            2737.7, 3003.9],
                'CH3_111': [55, 113.5, 167.4, 621.8, 686.0, 702.5, 1381.3, 
                            1417.5, 1575.8, 3026.6, 3093.2, 3098.9],
                'C-O_111': [],
                'H-OH_111': [],
                'H-C_111': []
                }


def make_input_file(file_name,energy_dict,frequency_dict):

    #create a header
    header = '\t'.join(['surface_name','site_name',
                        'species_name','formation_energy',
                        'frequencies','reference'])

    lines = [] #list of lines in the output
    for key in energy_dict.keys(): #iterate through keys
        E = energy_dict[key] #raw energy
        name,site = key.split('_') #split key into name/site
        if 'slab' not in name: #do not include empty site energy (0)
            frequency = frequency_dict[key]
            if site == 'gas':
                surface = None
            else:
                surface = 'Rh'
            outline = [surface,site,name,E,frequency,'Input File Tutorial.']
            line = '\t'.join([str(w) for w in outline])
            lines.append(line)

    lines.sort() #The file is easier to read if sorted (optional)
    lines = [header] + lines #add header to top
    input_file = '\n'.join(lines) #Join the lines with a line break

    input = open(file_name,'w') #open the file name in write mode
    input.write(input_file) #write the text
    input.close() #close the file

    print 'Successfully created input file'

file_name = 'energies.txt'
make_input_file(file_name,formation_energies,frequency_dict)

#Test that input is parsed correctly
from catmap.model import ReactionModel
from catmap.parsers import TableParser
rxm = ReactionModel()
#The following lines are normally assigned by the setup_file
#and are thus not usually necessary.
rxm.surface_names = ['Rh']
rxm.adsorbate_names = ('CO','C','O','H','CH','OH','CH2','CH3') 
rxm.transition_state_names = ('C-O','H-OH','H-C')
rxm.gas_names = ('CO_g','H2_g','CH4_g','H2O_g')
rxm.site_names = ('s',)
rxm.species_definitions = {'s':{'site_names':['111']}}
#Now we initialize a parser instance (also normally done by setup_file)
parser = TableParser(rxm)
parser.input_file = file_name
parser.parse()
#All structured data is stored in species_definitions; thus we can
#check that the parsing was successful by ensuring that all the
#data in the input file was collected in this dictionary.
for key in rxm.species_definitions:
    print key, rxm.species_definitions[key]
