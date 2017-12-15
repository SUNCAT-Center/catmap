# ASE Database

A general tutorial in using ase-db can be found [here](https://wiki.fysik.dtu.dk/ase/tutorials/db/db.html),
and documentation can be found [here](https://wiki.fysik.dtu.dk/ase/ase/db/db.html).

## How to make an ase-db readable by the db2catmap module in Catmap

In the end of your calculation script, you can add the lines:
    
    c = ase.db.connect('my_path_and_filename.db')
    c.write(atoms, species='OH', ads='OH-1', name='Pt', phase='fcc',
            facet='(111)', supercell='2x2', layers=3, site='top',
            n=1, data={'BEEFvdW_contribs': contribs})

which will write a row to your db file. If the db is not there already,
it will be created. The above keys are recommended for relaxed slab and adsorbate structures.
If the structure is a clean slab, put an empty string `''` in `species` and/or the string `'clean'` in `ads`.

For surfaces, the db2catmap module recognizes the following keys:

 - `energy` (immutable key retreived from the calculator)
 - `n`
 - `species`
 - `name`
 - `phase`
 - `facet`
 - `surf_lattice`
 - `supercell`
 - `layers`
 - `site`
 - `data['BEEFvdW_contribs']`
 - `data['frequencies']`

`site` is not recognized by default, but can be switched on by the option `site_specific` as seen further below.

## How to import the data

An example script is provided in this folder, which you can run like so:

    python generate_input.py my_input.txt

This will save a catmap input file in `my_input.txt`.
Look at the code in `generate_input.py`, and you will see it carries out four tasks:
1) Define search filters. This is needed if some data has to be filtered out.
    E.g. if some data was calculated with different parameters than other data. 
2) Import data from ase databases.
3) Store references and calculate formation energies.
4) Export catmap input file.

It is important to pay attention to the search filters. If you get garbage
results, it is likely because the search filters are not
sufficient for your dataset. Make sure you filter calculator parameters such as
XC-funcional, basis sets cutoffs, k-point sampling, ect, when necessary.
Importing data from correctly formatted .db files is done like so:
    
    project = db2catmap()
    project.get_molecules('molecules.db', selection=['fmax<0.05'])
    project.get_surfaces('surfaces.db', selection=['fmax<0.05'], site_specific=False)

The data is now attached to your db2catmap object.

The `site_specific` option accepts `True`, `False` or a string, in which case the site is recognized only if it matches the string.

## Get formation energies and export to catmap format.

Formation energies are calculated like so:

    references = (('H', 'H2_gas'), ('O', 'H2O_gas'), ('C', 'CO_gas'),))
    project.calc_formation_energies(references=references)

If you dont specify any references, the defaults will use H2 (g), H2O (g) and CH4 (g) and assume you have the potential energies of those molecules. 

    project.make_input_file(file_name)

finally saves the input file, compatible with catmap.

## How to store and import frequencies.

A convenient way of storing the frequencies is along with the atoms object in the `data` field. Using ase, you calculate frequencies as described [here](https://wiki.fysik.dtu.dk/ase/ase/vibrations/vibrations.html), and store them in the database like so:

    vib = Vibrations(atoms, indices=vibrateatoms, delta=0.03)
    vib.run()
    vib.summary(method='standard')
    frequencies = vib.get_frequencies()
    c_out.write(atoms, data={'frequencies': frequencies})

The field `data` always contains a dictionary. Use it with the key `frequencies` to make them accessible to the asedb2catmap module.
Importing frequencies is handled by the methods `get_surfaces` and `get_molecules`, which we have already used. It is necessary to pass the parameter `frequency_db` to it to import frequencies along with atomic structures like so:

    project.get_molecules('molecules.db', frequency_db='frequencies.db', selection=['fmax<0.05'])

## How to store and import transition states and pathways.

The module expects you to leave all images in a database file and distunguish
paths using the key `path_id`, which can be generated using the uuid module.
Doing so from a NEB can be done like so:

    path_id = uuid4().hex
    im = 0
    for atoms in images:
        im += 1
        epot = atoms.get_potential_energy()
        c_out.write(atoms, path_id=path_id, image=im)

Transition states and paths have the mandatory key value pairs:

 - `path_id`
 - `step` or `image`

`step` or `image` is used to order the images.
There is one additional recommended key:

 - `distance`

which is useful for making plots of the energy versus a reaction coordinate.
To add formation energies of transition states to your catmap input, you can use the method:

    project.get_transition_states('neb.db')

## Issues with ASE-incompatible calculators

ase-espresso is not fully compatible with the ase database module.
A workaround is to store the energy and forces in a `SinglePointDFTCalculator` like so:

    epot = atoms.get_potential_energy()
    forces = atoms.get_forces()
    relaxed_atoms = atoms.copy()
    spcalc = SinglePointDFTCalculator(relaxed_atoms)
    spcalc.results['energy'] = epot
    spcalc.results['forces'] = forces
    relaxed_atoms.set_calculator(spcalc)
    c.write(relaxed_atoms, ...)

`SinglePointDFTCalculator` can be imported from `ase.calculators.singlepoint`

Currently, the db2catmap module also recognizes the energy from the key `epot`, if a calculator is not attached.

