Electrochemistry
=================

Electrochemistry in CatMAP is implemented as a set of thermodynamic
corrections in the enthalpy\_entropy module. To run models containing
electrochemical reactions in CatMAP, the user needs to provide the
following in addition to all of the typical components of a CatMAP model
described in the tutorials:

-  Some sort of voltage-dependent species.  Currently, ``pe`` and ``ele`` are supported as
   fictitious gas molecules that represents the free energy
   of a proton-electron pair at 0V vs RHE or an electron at 0V vs SHE, respectively.
   As part of the Computational
   Hydrogen Electrode, the ``pe`` gas species should have the same free energy of
   1/2 of a molecule of H2 gas. This species's energy (usually 0 if an
   H2 reference is used) can be defined manually in the input file, but defaults
   to half the free energy of the H2 gas molecule in your system (if using ``pe``)
   or 0 if using "ele".  These are special electrochemistry-specific gas-phase species
   that should ignore your choice of free energy correction scheme.
-  In the reaction definitions, any potential-dependent reaction
   should have the ``pe_g`` gas molecule or ``ele_g`` in the
   initial (if reductive) or final (if oxidative) state. In addition, in
   the reaction definitions and input file, transition states for these
   reactions should contain ``pe`` or ``ele`` in the species name. For example, the
   reduction of adsorbed oxygen to adsorbed OH can be written as
   ``'O* + pe_g <-> O-pe* -> OH*'``.
-  An alternative formulation for electrochemical transition states is available
   via a special notation shown in the ORR tutorials.  Instead of explicitly defining
   a transition state species in your input file and using it in a reaction expression,
   you may instead write something like ``^0.26eV_a`` where the 0.26eV represents
   the free energy barrier of elementary step at that step's equilibrium/limiting potential
   on site ``a``.
-  A voltage (vs RHE) defined as a parameter in the .mkm file e.g.
   ``voltage = -0.2``.  This is also a global thermodynamic variable, and can be manipulated
   as a descriptor as shown in the ORR\_thermo tutorial.
-  A transfer coefficient defined as a parameter in the .mkm file e.g.
   ``beta = 0.5``.  This value can be defined for each elementary step through
   a reaction-specific flag like ``;prefactor=None, beta=0.45`` at the end of
   the corresponding reaction expression string.  See the ORR tutorials for an example
   of this notation.
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

There are a few provided examples of using electrochemistry in CatMAP in
the tutorials/electrochemistry directory. The HER example shows the
hydrogen evolution reaction in the low-coverage limit, and the provided
README file explains the details of the simplest of electrochemical
reactions as done in CatMAP. The ORR example is meant to reproduce the
results discussed in Hansen et al. (2014) DOI: 10.1021/jp4100608, which
used a home-spun microkinetic model. The README in that directory goes
into depth on how you can replicate some of the major features of this
work with the provided scripts.  The electron example shows how to equivalently
use CatMAP for an SHE potential versus an RHE one.
