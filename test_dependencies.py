import numpy
try:
    import scipy
except:
    print("Warning: scipy could not be imported! Some of the advanced features such as higher-order adsorbate-adsorbate interaction relying on arbitrary precision arithmethic may not for properly. Please install scipy for full feature support.")
import pylab
import mpmath
import ase
import catmap
