Electrochemistry
=================

Electrochemistry in CatMAP is implemented as a set of thermodynamic
corrections in the enthalpy\_entropy module. To run models containing
electrochemical reactions in CatMAP, the user needs to provide the
following in addition to all of the typical components of a CatMAP model
described in the tutorials:

-  A fictitious gas molecule called 'pe' that represents the free energy
   of a proton-electron pair at 0V vs RHE. As part of the Computational
   Hydrogen Electrode, this species should have the same free energy of
   1/2 of a molecule of H2 gas. This species's energy (usually 0 if an
   H2 reference is used) should be defined in the input file of
   energies.
-  In the reaction definitions, any reaction involving the use of a
   proton-electron pair should have the ``pe_g`` gas molecule in the
   initial (if reductive) or final (if oxidative) state. In addition, in
   the reaction definitions and input file, transition states for these
   reactions should contain ``pe`` in the species name. For example, the
   reduction of adsorbed oxygen to adsorbed OH can be written as
   ``'O* + pe_g <-> O-pe* -> OH*'``.
-  A voltage (vs RHE) defined as a parameter in the .mkm file e.g.
   ``voltage = -0.2``
-  A transfer coefficient defined as a parameter in the .mkm file e.g.
   ``beta = 0.5``
-  In the .mkm file, ``electrochemical_thermo_mode`` can be specified to
   one of three possible values (defaults to
   ``simple_electrochemical``): ``simple_electrochemical``,
   ``hbond_electrochemical``, and
   ``hbond_with_estimates_electrochemical``. "simple" only adds free
   energy corrections to adsorbates and transition states to account for
   voltage and beta. "hbond" will take a default or user-provided
   ``hbond_dict`` to correct each species for hydrogen bonding
   stabilization (see ``catmap/data/hbond_data.py`` for the default
   hbond\_dict). "hbond\_with\_estimates" attempts to estimate a
   hydrogen bonding correction for a given species based on its chemical
   formula. This is a very crude process that is described in
   ``catmap/thermodynamics/enthalpy_entropy.py``.

There are two provided examples of using electrochemistry in CatMAP in
the tutorials/electrochemistry directory. The HER example shows the
hydrogen evolution reaction in the low-coverage limit, and the provided
README file explains the details of the simplest of electrochemical
reactions as done in CatMAP. The ORR example is meant to reproduce the
results discussed in Hansen et al. (2014) DOI: 10.1021/jp4100608, which
used a home-spun microkinetic model. The README in that directory goes
into depth on how you can replicate some of the major features of this
work with the provided scripts.
