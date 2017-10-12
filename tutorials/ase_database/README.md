# ASE Database

## How to make an ase db readable by the db2catmap module in Catmap
In the end of your calculation script, you can add the lines:
    
    epot = get_potential_energy()
    c = ase.db.connect('my_database.db')
    c.write(atoms, epot=epot, species='OH', name='Pt', phase='fcc', supercell='2x2', layers=3, site='top')

which will write a row to your db file. If the file is not there already, it will be created.

## How to import the data

    python generate_input.py my_input.txt

will save the input file in `my_input.txt`.
