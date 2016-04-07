These scripts were written by Jordan Muscatello.

#----------------------------Description of files--------------------------------------------

graphene_sheet_gen.py generates a LAMMPS data file. The simulation box contains water molecules of a specified density with a given number of graphene sheets at the centre of the box. The sheets have a slit of the specified width with the potential for an offset in consecutive sheets.

vec.py, class_test.py and config_methods.py are all necessary to run graphene_sheet_gen.py

The in files are sample LAMMPS input files for this kind of setup.

#----------------------------Contribution by FJ----------------------------------------------

config_methods.py was amended by Frederike Jaeger to change the functionalisation on the graphene sheet edges (added COO, COOH, made sure COH was working) and to add a parameter determining the percentage of functionalisation.

graphene_sheet_gen.py was amended by FJ to include argparse parameters and facilitate the file structure of the generated input files.

