# hydrogen bonding dictionary used in 'hori' based on Pt(111) from Peterson(2010)
def generate_hbond_dict():
    """Returns a dictionary of generic, surface-agnostic hydrogen bond stabilizations, in eV."""
    OH = -0.50  # OH directly on the surface
    ROH = -0.25  # a 'floppy' OH group
    CO = -0.1  # carbon monoxide
    d = {'COOH'   : ROH,
         'OCHO'   : 0.,
         'CO'     : CO,
         'CHO'    : CO,
         'CH2O'   : 0.,
         'OCH3'   : 0.,
         'O'      : 0.,
         'OH'     : OH,
         'H'      : 0.,
         'COH'    : ROH,
         'C'      : 0.,
         'CH'     : 0.,
         'CH2'    : 0.,
         'CH3'    : 0.,
         'CHOH'   : ROH,
         'COHOH'  : ROH,
         'OCH2O'  : 0.,
         'CH2OH'  : ROH,
         'OCHCH2' : 0.,
         'OCHCHO' : 0.,
        }
    return d

hbond_dict = generate_hbond_dict()
