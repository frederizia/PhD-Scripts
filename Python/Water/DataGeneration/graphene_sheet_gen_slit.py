#!/usr/bin/python

## -- G R A P E N E _ S H E E T _ O H _ G E N . P Y -- ##
## Short script for generating functionalised graphene sheet configurations 
## with spce water for use with LAMMPS. 
##  
## Modified to functionalise surface with OH groups. 
## Slit size is fixed (try 5A before functionalisation)
##  
## Uses config_methods.py and "class_test.py" objects.

import math
import random
import argparse
from class_test import atom, bond, angle
from config_methods import *
    


## Parse arguments form command line
parser = argparse.ArgumentParser()
parser.add_argument("-nxy", type=int, nargs=2, \
                    help="Graphene unit cells in each direction",required=True)
parser.add_argument("-z", type=float, \
                    help="z length in Angtrom",required=True) 
parser.add_argument("-r", "--density", type=float, \
                    help="SPC/E Density in g/cm^3",required=True)
parser.add_argument("-o", "--offset", type=float, \
                    help="Offset l_off as fraction of y-length", default=0.5)
parser.add_argument("-n", "--nsheets", type=int, \
                    help="Number of sheets", default=2)
parser.add_argument("-s", "--spacing", type=float, \
                    help="Spacing between sheets", default=6.0)
parser.add_argument("-p", "--prob", type=float, \
                    help="Probability of functionalisation of sheet", default=0.5)
parser.add_argument("-delx", "--delx", type=float, \
                    help="Slit width", default=5.0)
parser.add_argument("-pe", "--probe", type=float, \
                    help="Probability of functionalisation of edge", default=0.5)
parser.add_argument("-f", "--f", type=str, \
                    help="Type of edge functionalisation", default='none')
parser.add_argument("-mol", "--molecule", type=str, \
                    help="Molecule or atom filling the box", default='spce')



args = parser.parse_args()
extents = args.nxy
DELX = args.delx
func = args.f
MOL = args.molecule

nx = extents[0]
ny = extents[1]

Lz = args.z
rho = args.density
off = float(args.offset)
if off < 0 or off > 1:
    print "Offset must be between 0 and 1" 
    exit()

nsheets = args.nsheets

p = args.prob
pe = args.probe

NA = 6.0221415e23

a = 1.42            # Spacing between carbon atoms (Angstrom)


## Dictionaries

# Just test for functionalisation 'CH' and 'CH_only'

if func == 'none':
    mass_index = {1:16.00, 2:1.00, 3:12.01}
    
    # Dictionary of types
    type_index = {'O':1, 'H':2, 'C':3}


if func == 'CH' or func == 'CH_only':

	if MOL == 'spce':
    		# Dictionary of masses
    		mass_index = {3:12.01, 1:16.00, 2:1.00, 4:1.00, 5:12.01, \
        	6:12.01, 7:16.00, 8:1.00, 9:12.01, 10:12.01, 11:16.00, 12:16.00}

    		# Dictionary of types
    		type_index = {'C':3, 'O':1, 'H':2, 'H_CH':4, 'C_CH':5,\
        	'C_COH':6, 'O_COH':7, 'H_COH':8, 'C1_CCOO':9, 'C2_CCOO':10, 			'O1_CCOO':11, 'O2_CCOO':12}


	elif MOL == 'He':
    		# Dictionary of masses
    		mass_index = {2:12.01, 1:4.00, 3:1.00, 4:12.01, \
        	5:12.01, 6:16.00, 7:1.00, 8:12.01, 9:12.01, 10:16.00, 11:16.00}

    		# Dictionary of types
    		type_index = {'He': 1, 'C': 2, 'H_CH':3, 'C_CH':4,\
        	'C_COH':5, 'O_COH':6, 'H_COH':7, 'C1_CCOO':8, 'C2_CCOO':9, 			'O1_CCOO':10, 'O2_CCOO':11}


if func == 'COH':
    # Dictionary of masses
    mass_index = {1:16.00, 2:1.00, 3:12.01, 4:12.01, 5:16.00, \
        6:1.00, 7:1.00, 8:12.01, 9:12.01, 10:12.01, 11:16.00, 12:16.00}
    
    # Dictionary of types
    type_index = {'O':1, 'H':2, 'C':3,'C_COH':4, 'O_COH':5, \
        'H_COH':6, 'H_CH':7, 'C_CH':8,'C1_CCOO':9, 'C2_CCOO':10, 'O1_CCOO':11,'O2_CCOO':12}

if func == 'CCOO':
    # Dictionary of masses
    mass_index = {1:16.00, 2:1.00, 3:12.01, 4:12.01, 5:12.01, \
        6:16.00, 7:16.00, 8:12.01, 9:16.00, \
        10:1.00, 11:1.00, 12:12.01, }
    
    # Dictionary of types
    type_index = {'O':1, 'H':2, 'C':3,'C1_CCOO':4, 'C2_CCOO':5, \
        'O1_CCOO':6, 'O2_CCOO':7, 'C_COH':8, 'O_COH':9, \
        'H_COH':10, 'H_CH':11, 'C_CH':12}

if func == 'CCOOH':
    # Dictionary of masses
    mass_index = {1:16.00, 2:1.00, 3:12.01, 4:12.01, 5:12.01, \
        6:16.00, 7:16.00, 8:1.0, 9:12.01, 10:16.00, \
        11:1.00, 12:1.00, 13:12.01, }
    
    # Dictionary of types
    type_index = {'O':1, 'H':2, 'C':3,'C1_CCOOH':4, 'C2_CCOOH':5, \
        'O1_CCOOH':6, 'O2_CCOOH':7, 'H_CCOOH':8,'C_COH':9, 'O_COH':10, \
        'H_COH':11, 'H_CH':12, 'C_CH':13}



## Lists for atoms, types and angles
atom_list = []
mol_count = 0
bond_list = []
bond_count = 0
angle_list = []
angle_count = 0
bond_types = []
angle_types = []
type_list = []
improper_list = []
improper_types = []
dihedral_list = []
dihedral_types= []

## System extents
delx = DELX  # slit width?
delz = args.spacing

## Define graphene cell
graphene_cell, a1, a2 = graphene_unit_cell(a)
Lx = nx*a1 + delx
Ly = ny*a2
off = Lx*off
x_off = 0.0
y_off = 0.0
z_off = 0.0

## Generate graphene sheets
for i in range(0,nsheets):
    ## offset alternate sheets
    print i 
    if i % 2 == 0: 
        x_off = 0.0
    else: 
        x_off = off

    atom_list, mol_count, type_list = duplicate_hex_cell(graphene_cell, \
                        atom_list, nx, ny, x_off, y_off, z_off, a, a1, a2, \
                                  mol_count, type_list, "n")
    z_off += delz

print "Number of atoms in sheets = %d" % len(atom_list)

print "Box extents:\nLx = %f\nLy = %f\nLz = %f\n" %  (Lx, Ly, Lz)

if MOL == 'spce':

	## Define SPC/E molecule
	spce, spce_bond, spce_angle, type_list, bond_types, angle_types \
                       = spce_molecule(type_list, bond_types, angle_types)

	write_lammps_molecule(spce, spce_bond, spce_angle, [], [], 
                      type_index, mass_index, "SPCE", "spce molecule")

elif MOL == 'He':		
	helium, type_list, helium_bond, helium_angle \
			= helium_atom(type_list)


## Apply bonding to graphene
bond_list, bond_types = add_bonds_nn(atom_list, bond_list, bond_types,
                                     0.0, Lx, 0.0, Ly, 0.0, Lz, a + 0.01)

## Add functional groups to edges
y_dist = 2.0*a2/3.0

#atom_list, bond_list, angle_list, improper_list, dihedral_list, mol_count, \
#bond_types, angle_types, improper_types, dihedral_types, type_list = \
#functionalise_edge3(atom_list, bond_list, angle_list, improper_list, \
                    # dihedral_list, mol_count, bond_types, angle_types, \
                    # improper_types, dihedral_types, type_list, \
                    # a, pe, Lx, Ly, 'y', '%s'%(func), y_dist)



print "lx = %f" % Lx

atom_list = apply_pbc(atom_list, Lx, Ly, Lz)

write_xyz(atom_list, "graphene_sheet", "A graphene sheet with a slit!")

zlo = z_off -delz # place lattice zlo from carbon

## Functionalise graphene surface
'''atom_list, bond_list,  angle_list, dihedral_list, mol_count, \
bond_types, angle_types, dihedral_types, type_list = \
functionalise_sheet2(atom_list, bond_list, angle_list, dihedral_list,\
       mol_count, bond_types, angle_types, dihedral_types, type_list,\
       'COH', a, p, 1, Lx, Ly, Lz)'''

## Random Insertion

r_tol = 1.2 # Tolerance for random insertion


if MOL == 'spce':

    atom_list, bond_list, angle_list, improper_list, dihedral_list,\
    mol_count, Nm, m = random_insertion_grid(atom_list, bond_list, angle_list,
                                         improper_list, dihedral_list,
                                         mol_count, type_index, mass_index,
                                         spce, spce_bond, spce_angle, [], [],
                                         0.0, Lx, 0.0, Ly, z_off*0.7, Lz*0.98, rho,
                                         r_tol, 11, 11, 11)
    
    print 'test'
    density = (Nm*m/NA)/(Lx*Ly*Lz*1e-24)
    print "Water density in cell = %f  g/cm^3" % density
    density_block = (Nm*m/NA)/(Lx*Ly*(Lz-z_off/2)*1e-24)
    print "Water density between the sheets = %f  g/cm^3" % density_block

elif MOL == 'He':
	atom_list, bond_list, angle_list, improper_list, dihedral_list,\
	mol_count, Nm, m = random_insertion_grid(atom_list, bond_list, angle_list,
                                         improper_list, dihedral_list,
                                         mol_count, type_index, mass_index,
                                         helium, helium_bond, helium_angle, [], [],
                                         0.0, Lx, 0.0, Ly, z_off*0.5, Lz, rho,
                                         r_tol, 9, 7, 10)


	density = (Nm*m/NA)/(Lx*Ly*Lz*1e-24)
	print "Helium density in cell = %f  g/cm^3" % density

#atom_list = shift_z(atom_list, Lz)


atom_list = create_slit(atom_list, delz, Lz)

write_xyz(atom_list, "graphene_slit_%i_%i_%i_%.2f" % (nx, ny, Lz, rho), "A graphene sheet")


write_lammps_data3(atom_list, bond_list, angle_list, improper_list,  
                   dihedral_list, mass_index, type_index, type_list,  
                   bond_types, angle_types, improper_types, dihedral_types, 
                   Lx, Ly, Lz, True,  
                   "graphene_slit_%i_%i_%i_%.2f" % (nx, ny, Lz, rho), "LAMMPS data for Graphene Sheet", delz)

