# ASE Database

This tutorial will introduce you to the CatMAP API for ASE databases.
The ASE-db API is made to interpret the atomic structures in your ase-db files as points in a potential energy landscape, on which the kinetic model is built. The `catmap.api.energy_landscape` class handles imports from the databases, calculation of formation energies for unique stable molecules and adsorbate states, and finally exports them to catmap energy txt files, as presented in tutorial 1.

A general tutorial in using ase-db can be found [here](https://wiki.fysik.dtu.dk/ase/tutorials/db/db.html),
and documentation can be found [here](https://wiki.fysik.dtu.dk/ase/ase/db/db.html).

Go to [this Jupyter Notebook](https://wiki.fysik.dtu.dk/ase/ase/db/db.html) to proceed with CatMAP's tutorial on importing ase-db data, calculation formation energies, and exporting the CatMAP input data file.

## How to make an ASE-db readable by the `EnergyLandscape` module in CatMAP

In the end of your calculation script, you can add the lines:
    
    import ase.db
    c = ase.db.connect('my_path_and_filename.db')
    c.write(atoms, species='OH', name='Pt',
            n=1, data={'BEEFvdW_contribs': contribs})  # contribs are the 32 non-selfconsistent energies.

which will write a row to your db file. If the db is not there already, it will be created. 
This will contain is all the basic information you need.

If the structure is a clean slab, put an empty string `''` in `species` and/or the string `'clean'` in `ads`.

For more advanced features the below keys are needed. The `EnergyLandscape` module recognizes and uses the following keys:

 - `energy` (immutable key retreived from an attached calculator)
 - `n`
 - `species`
 - `name`
 - `phase`, `crystal` (if both are present, `phase` will apply.)
 - `facet`
 - `surf_lattice`
 - `supercell`
 - `layers`
 - `site`
 - `data['BEEFvdW_contribs']`
 - `data['frequencies']`

`site` is not recognized by default, but can be switched on by a parameter `site_specific` as seen further below.

## How to make the CatMAP energy data file from the ASE-db.

An example jupyter notebook is provided in this folder. You can run it by

These will save a catmap input file in `my_input.txt`. Look at the code and you will see the procedure can be itemized in four tasks:
1) Define search filters. This is needed if some data has to be filtered out.
    E.g. if some data was calculated with different parameters than other data.
2) Import data from ase databases.
3) Store references and calculate formation energies.
4) Export catmap input file.

It is important to pay attention to the search filters. If you get garbage
results, it is likely because the search filters are not
sufficient for your dataset. Make sure you filter calculator parameters such as
XC-functional, basis sets cutoffs, k-point sampling, ect., when necessary.
Importing data from correctly formatted .db files is done like so:
    
    from catmap.api.ase_data import energy_landscape
    project = energy_landscape()
    project.get_molecules('molecules.db', selection=['fmax<0.05'])
    project.get_surfaces('surfaces.db', selection=['fmax<0.05'], site_specific=False)

The `site_specific` option accepts `True`, `False` or a string. In the latter case, the `site` key is recognized only if the value matches the string, while all other sites are treated as identical.

Your data is now stored in dictionaries that are attached to your `energy_landscape` object.

## Get formation energies and export to catmap format.

Formation energies are calculated like so:

    references = (('H', 'H2_gas'), ('O', 'H2O_gas'), ('C', 'CO_gas'),)
    project.calc_formation_energies(references)

`references` is a required parameter, that should contain gas phase references. If a gas phase reference is dependent on another, order the dependent one after the latter.

    project.make_input_file(file_name)

finally saves the input file, compatible with catmap.

## How to store and import frequencies.

A convenient way of storing the frequencies is along with the atoms object in the `data` field. Using ase, you calculate frequencies as described [here](https://wiki.fysik.dtu.dk/ase/ase/vibrations/vibrations.html), and store them in the database like so:

    vib = Vibrations(atoms, indices=vibrateatoms, delta=0.03)
    vib.run()
    vib.summary(method='standard')
    frequencies = vib.get_frequencies()
    c_out.write(atoms, data={'frequencies': frequencies})

Importing frequencies is handled by the methods `get_surfaces` and `get_molecules`, which we have already used. It is necessary to pass the parameter `frequency_db` to it to import frequencies along with atomic structures like so:

    project.get_molecules('molecules.db', frequency_db='frequencies.db', selection=['fmax<0.05'])

## How to store and import transition states and pathways.

The module expects you to leave all images in a database file and distunguish
paths using the key `path_id`, which can be generated using the uuid module.
Doing this after a NEB calculation can be done like so:

    path_id = uuid4().hex
    for im, atoms in enumerate(images):
        c_out.write(atoms, path_id=path_id, image=im)

Transition states and paths have the mandatory key value pairs:

 - `path_id`
 - `step` or `image`

`step` or `image` is used to order the images.
There is one additional recommended key value pair:

 - `distance`

which is useful for making plots of the energy versus a reaction coordinate.
To add formation energies of transition states to your catmap input, you can use the method:

    project.get_transition_states('neb.db')

## Issues with ASE-incompatible calculators

jvoss/ase-espresso is not fully compatible with ASE and it's database module.
A workaround is to store the energy and forces in a `SinglePointDFTCalculator` like so:

    from ase.calculators.singlepoint import SinglePointDFTCalculator
    epot = atoms.get_potential_energy()
    forces = atoms.get_forces()
    relaxed_atoms = atoms.copy()
    spcalc = SinglePointDFTCalculator(relaxed_atoms)
    spcalc.results['energy'] = epot
    spcalc.results['forces'] = forces
    relaxed_atoms.set_calculator(spcalc)
    c.write(relaxed_atoms, ...)

Currently, the `EnergyLandscape` module also recognizes the energy from the key `epot`, if a calculator is not attached.

## Retrieving data from catalysis-hub.org.

Please see tutorial 2.
