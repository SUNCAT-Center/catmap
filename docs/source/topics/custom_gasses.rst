Additional Gas Species Properties
================================

There are two main ways of using free energies in CatMAP:

1. **Frozen everything** - Useful if you've already calculated free energies for all species yourself, you can tell CatMAP
to use no free energy corrections by specifying ``frozen_gas`` and ``frozen_adsorbate`` as the ``gas_thermo_mode`` and ``adsorbate_thermo_mode``,
respectively.  If you are using scaling relationships with this method, keep in mind that linear scaling relationships are only
formulated for formation energies rather than free energies.

2. **Let CatMAP do it** - CatMAP has built-in functions to estimate free energies for gasses and adsorbates based on data you
provide and its own internal database of parameters.  The rest of this tutorial will go over how you can customize these
methods to your liking.  Currently, all corrections to electronic energies to get to free energies are located in the ``thermodynamics``
module.

**Ideal Gas Properties**

``ideal_gas`` is the default mode for correcting the energies of gas phase species.  It relies on `ASE's thermochemistry 
package <https://wiki.fysik.dtu.dk/ase/ase/thermochemistry/thermochemistry.html>`__, user-provided vibrational frequencies,
and a dictionary of ideal gas parameters, ``ideal_gas_params``, which stores ``[symmetry_number, geometry, spin]`` values for
each gas species key.  CatMAP provides and uses a default ``ideal_gas_params`` dictionary in ``catmap.data.ideal_gas_params``,
but you can provide your own parameters by adding something like the following to your .mkm file:

.. code:: python

  ideal_gas_params = {'Cl2_g':[2, 'linear', 0],
                      'F2_g':[2, 'linear', 0],}

and so on for each species in your system that CatMAP does not already have parameters for.

**Shomate Gas Properties**

In the tutorials, we use the Shomate equation to estimate free energy corrections for the gas phase molecules.  Shomate
parameters can be found in databases such as `NIST <http://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Type=JANAFL&Table=on#JANAFL>`__
or they can be fitted from experimental data with the ``fit_shomate`` function in ``catmap.thermodynamics``, which takes in lists of
temperatures, heat capacities, enthalpies, entropies, and initial guesses for the Shomate parameters and returns the least-squares fitted Shomate
parameters A, B, C, D, E, F, G, and H.  See ``fit_shomate.py`` in the ``custom_gasses`` tutorial for examples on how to do this.  This tutorial
also generates a plot that compares how the free energies generated from the Shomate equation and Ideal Gas equations can deviate as a function
of temperature.

If you have your own Shomate parameters you wish to use or you have calculated them using ``fit_shomate``, you can use them in your microkinetic model
by modifying your model's ``shomate_params`` attribute - which is a dictionary of ``gas_name:temperature_range`` keys to a list of 
Shomate parameters as values.  You can input this dictionary in your .mkm file like the example below:

.. code:: python

  shomate_params = {'CH3OH_g:298-1500':[-1.0846, 153.2464, -79.5305, 16.4713, 0.5220, -4.8974, 199.1894, 0.0],
                    'CH3CH2OH_g:298-1200':[-4.7368, 271.9618, -169.3495, 43.7386, 0.2464, -9.8283, 203.3326, 0.0],}


**Contributing Data to CatMAP**

If you have Shomate data or additional Ideal Gas parameter data that you'd like to see included in CatMAP, feel free to send us a `pull
request <https://help.github.com/articles/using-pull-requests/>`__ with your updated ``parameter_data.py`` file.  We will most likely only
include such new data in the main CatMAP repository only if the Shomate parameters are fit with well-established thermodynamic data (such as
NIST or the CRC) or if the Ideal Gas parameters are consistent with the geometries in `ASE's database <https://wiki.fysik.dtu.dk/ase/ase/structure.html>`__.

**Writing Your Own Corrections**

Modifying CatMAP source code to include your own custom free energy corrections may seem daunting, but it may be much easier than you think.
Check out the documentation for ``enthalpy_entropy.py``, specifically how ``ideal_gas`` and ``shomate_gas`` are implemented.  All a thermodynamic
correction needs to do is to return a dictionary of corrections it is reponsible for.  At the point your correction function is called, the ReactionModel
instance is already fully initialized, so you have access to all attributes of the reaction model from within your correction function.