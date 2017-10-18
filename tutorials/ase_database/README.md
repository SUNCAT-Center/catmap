# ASE Database

## How to make an ase-db readable by the db2catmap module in Catmap
In the end of your calculation script, you can add the lines:
    
    epot = get_potential_energy()
    c = ase.db.connect('my_path_and_filename.db')
    c.write(atoms, epot=epot, species='OH', name='Pt', phase='fcc',
            facet='1x1x1', supercell='2x2', layers=3, site='top')

which will write a row to your db file. If the file is not there already,
it will be created. The above keys are suitable for relaxed geometries.

## How to import the data

    python generate_input.py my_input.txt

will save the input file in `my_input.txt`.
In `generate_input.py`, you will see it carries out four steps:
1) Define search filters. Needed if some data has to be filtered out.
    E.g. if some data was calculated with different parameters than other data. 
2) Import data from ase databases.
3) Calculate references and formation energies.
4) Export catmap input file.

It is important to pay attention to the search filters. If you get garbage
results, it is likely because the search filters are not
sufficient for your dataset. Make sure you filter calculator parameters such as
XC-funcional, basis sets cutoffs, k-point sampling, ect.
For surfaces, the module have built-in distinctions for the keys:

 - name
 - phase
 - facet
 - surf_lattice
 - supercell
 - layers
 - site

Site is optional and turned off by default. The potential energy is retrieved
from the key `energy` or `epot`, whichever is present.

## How to store transition states and pathways.

The module expects you to leave all images in a database file and distunguish
paths using the key `path_id`, which can be generated using the uuid module.
Doing so from a NEB can be done like so:

    path_id = uuid4().hex
    for atoms in images:
        im += 1
        epot = atoms.get_potential_energy()
        forces = atoms.get_forces()
        fmaxout = np.sqrt((forces**2).sum(axis=1).max())
        relaxed_atoms = atoms.copy()
        relaxed_atoms.set_calculator(None)
        c_out.write(relaxed_atoms, epot=epot, path_id=path_id)

Transition states and paths have the mandatory keys:

 - path_id
 - step or image


'step' or `image` is used to order the images.
There is one recommended key:

 - distance

which is useful for making plots of the energy versus reaction coordinate.
