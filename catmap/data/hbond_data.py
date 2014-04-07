# hydrogen bonding dictionary used in 'hori' based on Pt(111) from Peterson(2010)
def generate_hbond_dict():
    """Returns a dictionary of hydrogen bond stabilizations, in eV."""
    OH = -0.50  # OH directly on the surface
    ROH = -0.25  # a 'floppy' OH group
    CO = -0.1  # carbon monoxide
    d = {'COOH_s'   : ROH,
         'OCHO_s'   : 0.,
         'CO_s'     : CO,
         'CHO_s'    : CO,
         'CH2O_s'   : 0.,
         'OCH3_s'   : 0.,
         'O_s'      : 0.,
         'OH_s'     : OH,
         'H_s'      : 0.,
         'COH_s'    : ROH,
         'C_s'      : 0.,
         'CH_s'     : 0.,
         'CH2_s'    : 0.,
         'CH3_s'    : 0.,
         'CHOH_s'   : ROH,
         'COHOH_s'  : ROH,
         'OCH2O_s'  : 0.,
         'CH2OH_s'  : ROH,
         'OCHCH2_s' : 0.,
         'OCHCHO_s' : 0.,
        }
    return d

hbond_dict = generate_hbond_dict()
