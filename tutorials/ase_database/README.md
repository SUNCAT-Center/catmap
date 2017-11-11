# ASE Database

## How to make an ase-db readable by the db2catmap module in Catmap

In the end of your calculation script, you can add the lines:
    
    epot = atoms.get_potential_energy()
    contribs = calc.get_nonselfconsistent_energies()
    c = ase.db.connect('my_path_and_filename.db')
    c.write(atoms, epot=epot, species='OH', name='Pt', phase='fcc',
            facet='1x1x1', supercell='2x2', layers=3, site='top',
            n=1, data={'BEEFvdW_contribs': contribs})

which will write a row to your db file. If the file is not there already,
it will be created. The above keys are suitable for relaxed slab and adsorbate structures.
If the structure is a clean slab, put an empty string `''` in species.

## How to import the data

    python generate_input.py my_input.txt

will save the input file in `my_input.txt`.
In the script `generate_input.py`, you will see it carries out four steps:
1) Define search filters. Needed if some data has to be filtered out.
    E.g. if some data was calculated with different parameters than other data. 
2) Import data from ase databases.
3) Calculate references and formation energies.
4) Export catmap input file.

It is important to pay attention to the search filters. If you get garbage
results, it is likely because the search filters are not
sufficient for your dataset. Make sure you filter calculator parameters such as
XC-funcional, basis sets cutoffs, k-point sampling, ect.
For surfaces, the module have built-in distinctions that uses the keys:

 - n
 - species
 - name
 - phase
 - facet
 - surf_lattice
 - supercell
 - layers
 - site
 - data.BEEFvdW_contribs

`site` is optional and turned off by default. The potential energy is retrieved
from the key `energy` or `epot`, whichever is present.

Importing data from correctly formatted .db files is done like so:
    
    project = db2catmap()
    project.get_molecules('molecules.db')
    project.get_surfaces('surfaces.db')

The data is now attached to your asedb2catmap object.

## Get formation energies and export to catmap format.

Formation energies are calculated like so:

    references = (('H', 'H2_gas'), ('O', 'H2O_gas'), ('C', 'CO_gas'),))
    project.calc_formation_energies(references=references)

If you dont specify any references, the defaults will use H, H2O and CH4 and assume you have the potential energies of those molecules. 
To finally save the catmap input file do:

    project.make_input_file(file_name)

## How to store and import frequencies.

A convenient way of storing the frequencies is along with the atoms object in the `data` field. Using ase, you calculate and store the frequencies like so:

    vib = Vibrations(atoms, indices=vibrateatoms, delta=0.03)
    vib.run()
    vib.summary(method='standard')
    frequencies = vib.get_frequencies()
    c_out.write(atoms, epot=epot, data={'frequencies': frequencies})

Use the key `frequencies` to make it accessible to the asedb2catmap module.
Importing frequencies is handled by `get_surfaces`, which we have already used. 
It is necessary to pass the parameter `frequency_db` to it to import frequencies along with atomic structures.

## How to store and import transition states and pathways.

The module expects you to leave all images in a database file and distunguish
paths using the key `path_id`, which can be generated using the uuid module.
Doing so from a NEB can be done like so:

    path_id = uuid4().hex
    im = 0
    for atoms in images:
        im += 1
        epot = atoms.get_potential_energy()
        # ase-espresso calculator cannot be written to db.
        relaxed_atoms = atoms.copy()
        relaxed_atoms.set_calculator(None)
        c_out.write(relaxed_atoms, epot=epot, path_id=path_id, image=im)

Transition states and paths have the mandatory keys:

 - energy or epot
 - path_id
 - step or image

`step` or `image` is used to order the images.
There is one additional recommended key:

 - distance

which is useful for making plots of the energy versus reaction coordinate.
To add formation energies of transition states to your catmap input, you can simply do:

