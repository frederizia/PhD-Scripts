# Contains all functions for generating graphene/water configurations
# for LAMMPS
import math
import random
import numpy
from class_test import *
from vec import *
from numpy import matrix
import numpy as np


def add_bonds_nn(atom_list, bond_list, bond_types, 
                 xlo, xhi, ylo, yhi, zlo, zhi, r_c):
    """
    Adds bonds to neighbouring monatomic species.
    To be used with graphene
    """
    Lx = xhi - xlo
    Ly = yhi - ylo
    Lz = zhi - zlo
 
    print "Lx = {}, Ly = {}, Lz = {}".format(Lx, Ly, Lz)

    l_x = matrix([[Lx, 0.0, 0.0]])
    l_y = matrix([[0.0, Ly, 0.0]])
    l_z = matrix([[0.0, 0.0, Lz]])

    new_bond_list = []
 
    r_c2 = r_c*r_c

    b_type = len(bond_types) + 1
    bond_types.append(b_type)
        
    bond_count = 0

    for i in range(0, len(atom_list)):        
        r_i = atom_list[i].return_r_matrix()
        for j in range(0, len(atom_list)):
            r_j = atom_list[j].return_r_matrix()
            if j != i:
                # check new_bondlist for duplicates
                dpl_flg = 0
                b = 0
                while not dpl_flg and b < len(new_bond_list):
                    
    #               if ((new_bond_list[b].atom1 == i + 1
     #                  and new_bond_list[b].atom2 == j + 1) or 
      #                 (new_bond_list[b].atom2 == i + 1 
       #                and new_bond_list[b].atom1 == j + 1)): 
                    if ((new_bond_list[b].atom1 == atom_list[i].atom 
                        and new_bond_list[b].atom2 == atom_list[j].atom) or 
                        (new_bond_list[b].atom2 == atom_list[i].atom 
                        and new_bond_list[b].atom1 == atom_list[j].atom)): 
                           
                        dpl_flg = 1  

                    b += 1  
                    
                r_ij = r_i - r_j
                #PBCs
                if r_ij[0,0] > Lx/2: r_ij -= l_x
                if r_ij[0,0] < -Lx/2: r_ij += l_x
                if r_ij[0,1] > Ly/2: r_ij -= l_y
                if r_ij[0,1] < -Ly/2: r_ij += l_y
                if r_ij[0,2] > Lz/2: r_ij -= l_z
                if r_ij[0,2] < -Lz/2: r_ij += l_z               
   #             print "i = {}, j = {}, tol = {}, r_ij = {}" .format(i, j, r_c, math.sqrt(float(r_ij*r_ij.T)))  
                if float(r_ij*r_ij.T) < r_c2 and not dpl_flg:
                    bond_count += 1
                    new_bond = bond(b_type, atom_list[i].atom, 
                                            atom_list[j].atom)
                    new_bond_list.append(new_bond)
                    
    print "{} bonds were created".format(bond_count)

    bond_list.extend(new_bond_list)
    return bond_list, bond_types 


def graphene_unit_cell(a):
    "Returns unit cell for graphene sheet"
    rt3 = math.sqrt(3.0)

    a1 = 3*a            # lattice constant x direction
    a2 = 3*rt3*a/2      # lattice constant y direction

    # Define graphene unit cell
    graphene_unit = []
    graphene_unit.append(atom(a/2, 0, 0, 'C', 0.0, 1))
    graphene_unit.append(atom(3*a/2, 0, 0, 'C', 0.0, 1))
    graphene_unit.append(atom(0, rt3*a/2, 0, 'C', 0.0, 1))
    graphene_unit.append(atom(2*a, rt3*a/2, 0, 'C', 0.0, 1))
    graphene_unit.append(atom(a/2, rt3*a, 0, 'C', 0.0 , 1))
    graphene_unit.append(atom(3*a/2, rt3*a, 0, 'C', 0.0, 1))

    return graphene_unit, a1, a2

def helium_atom(type_list):
    """
    Returns parameters for helium atom
    '6:He' or '9:He'
    Defines no angles or bonds
    """

    helium = []

    helium_types = ['He']

    type_list.extend(helium_types)

    helium.append(atom(0, 0, 0, 'He', 0.0, 1))


    # create empty lists as random_insertion_grid requires molecule bond list as input
    helium_bonds = []
    helium_angles = []

    return helium, type_list, helium_bonds, helium_angles

def argon_atom(type_list):
    """
    Returns parameters (?) associated with a single argon atom.
    """
   
    argon = []

    argon_types = ['Ar']   

    type_list.extend(argon_types)
  
    argon.append(atom(0.0, 0.0, 0.0, 'Ar', 0.0, 1)) 

    return argon, type_list
  

    type_list.extend(argon_types)

def Na_ion(type_list):
    """
    Returns parameters (Smith and Dang 1994) associated with sodium ion.
    """
    
    Na = []

    Na_types = ['Na']

    type_list.extend(Na_types)

    Na.append(atom(0.0, 0.0, 0.0, 'Na', 1.00, 1))

    return Na, type_list

def Cl_ion(type_list):
    """
    Returns parameters (Smith and Dang 1994) associated with chloride ion.
    """
    
    Cl = []

    Cl_types = ['Cl']

    type_list.extend(Cl_types)

    Cl.append(atom(0.0, 0.0, 0.0, 'Cl', 1.00, 1))

    return Cl, type_list

def ethanol_molecule(type_list, bond_types, angle_types):
    """
    Returns list of atoms, bond and angles corresponding 
    to ethanol molecule: CH3CH2OH
    Defines atoms: H2,C2,H1,C1,O_e,H_e

    """

    pi = 4.0*math.atan(1.0)

    # ANGLES
    ang_H1C1Oe = 109.44*pi/180 # bond angle in radians  #1
    ang_H1C1C2 = 109.46*pi/180 #2
    ang_H1C1H1 = 120.0*pi/180 #3
    ang_OeC1C2 = 109.0*pi/180 #4
    ang_HeOeC1 = 109.5*pi/180 #5
    ang_H2C2H2 = 109.5*pi/180 #6
    ang_H2C2C1 = 109.46*pi/180 #7

    # BOND LENGTHS
    l_C1C2 = 1.52             # bond length in Angstrom
    l_CH = 1.11
    l_C1Oe = 1.42 
    l_O2H2 = 0.94

    ethanol = []
    ethanol_bonds = []
    ethanol_angles = []

    ethanol_types = ['H2', 'C2', 'H1', 'C1', 'O_e', 'H_e']

    angles = [7]    # should this be a len(angle_types) + 7?
    bonds = [len(bond_types) + 4]
    

    type_list.extend(ethanol_types)

    q_C1 = 0.05
    q_H1 = 0.09 #check
    q_C2 = -0.27
    q_H2 = 0.09 #check
    q_O_e = -0.66
    q_H_e = 0.43

    # distances for hydrogen atoms wrt to oxygen 
    h1x = l*math.cos(ang/2.0)
    h1y = l*math.sin(ang/2.0)
    h2x = h1x
    h2y = -h1y


    # Define ethanol molecule
    ethanol.append(atom(0, 0, 0, 'C2', q_C2, 1)) #1
    ethanol.append(atom(0, 0, 0, 'H21', q_H2, 1)) #2
    ethanol.append(atom(0, 0, 0, 'H22', q_H2, 1)) #3
    ethanol.append(atom(0, 0, 0, 'H23', q_H2, 1)) #4
    ethanol.append(atom(0, 0, 0, 'C1', q_C1, 1)) #5
    ethanol.append(atom(0, 0, 0, 'H11', q_H1, 1)) #6
    ethanol.append(atom(0, 0, 0, 'H12', q_H1, 1)) #7
    ethanol.append(atom(0, 0, 0, 'O_e', q_O_e, 1)) #8
    ethanol.append(atom(0, 0, 0, 'H_e', q_H_e, 1)) #9

    # Bonds
    ethanol_bonds.append(bond(1, 1, 2)) # C2H21
    ethanol_bonds.append(bond(1, 1, 3)) # C2H22
    ethanol_bonds.append(bond(1, 1, 4)) # C2H23
    ethanol_bonds.append(bond(2, 1, 5)) # C2C1
    ethanol_bonds.append(bond(1, 5, 6)) # C1H11
    ethanol_bonds.append(bond(1, 5, 7)) # C1H12
    ethanol_bonds.append(bond(2, 5, 8)) # C1Oe
    ethanol_bonds.append(bond(3, 8, 9)) # OeHe

    # Angles
    ethanol_angles.append(angle(6, 2, 1,3)) #H21C2H22
    ethanol_angles.append(angle(6, 2, 1,4)) #H21C2H23
    ethanol_angles.append(angle(6, 4, 1,3)) #H23C2H22
    ethanol_angles.append(angle(3, 6, 5,7)) #H11C1H12
    ethanol_angles.append(angle(2, 1, 5,6)) #C2C1H11
    ethanol_angles.append(angle(2, 1, 5,7)) #C2C1H12
    ethanol_angles.append(angle(7, 2, 1,5)) #H21C2C1
    ethanol_angles.append(angle(7, 3, 1,5)) #H22C2C1
    ethanol_angles.append(angle(7, 4, 1,5)) #H23C2C1
    ethanol_angles.append(angle(1, 6, 5,8)) #H11C1Oe
    ethanol_angles.append(angle(1, 7, 5,8)) #H12C1Oe
    ethanol_angles.append(angle(4, 1, 5,8)) #C2C1Oe
    ethanol_angles.append(angle(5, 6, 8,9)) #C1OeHe

    #Dihedrals
    ethanol_dihedrals.append(dihedral(1,1,5,8,9)) # C2C1OeHe

# Impropers?



    bond_types.extend(bonds)
    angle_types.extend(angles)


    return ethanol, ethanol_bonds, ethanol_angles, type_list, bond_types, angle_types


def spce_molecule(type_list, bond_types, angle_types):
    """
    Returns list of atoms, bond and angles corresponding 
    to spc/e water molecule:
    1:O,2:H,3:H
    Defines angle and bond type 1
    """

    pi = 4.0*math.atan(1.0)

    ang = 109.47*pi/180 # bond angle in radians
    l = 1.0             # bond length in Angstrom

    spce = []
    spce_bonds = []
    spce_angles = []

    spce_types = ['O', 'H']

    angles = [1]
    bonds = [len(bond_types) + 1]
    

    type_list.extend(spce_types)

    q_O = -0.8476
    q_H = 0.4238

    # distances for hydrogen atoms wrt to oxygen 
    h1x = l*math.cos(ang/2.0)
    h1y = l*math.sin(ang/2.0)
    h2x = h1x
    h2y = -h1y
    spce.append(atom(0, 0, 0, 'O', q_O, 1))
    spce.append(atom(h1x, h1y, 0, 'H', q_H, 1))
    spce.append(atom(h2x, h2y, 0, 'H', q_H, 1))

    spce_bonds.append(bond(1, 1, 2))
    spce_bonds.append(bond(1, 1, 3))
    spce_angles.append(angle(1, 2, 1,3))

    bond_types.extend(bonds)
    angle_types.extend(angles)


    return spce, spce_bonds, spce_angles, type_list, bond_types, angle_types


def spce_unit_cell(a):
    "Returns fcc unit cell for spc/e water. Takes lattice constant argument"
    pi = 4.0*math.atan(1.0)

    ang = 109.47*pi/180 # bond angle in radians
    l = 1.0             # bond length in Angstrom

    # distances for hydrogen atoms wrt to oxygen 
    h1x = l*math.cos(ang/2.0)
    h1y = l*math.sin(ang/2.0)
    h2x = h1x
    h2y = -h2x

    spce_unit = []

    spce_unit.append(atom(0, 0, 0, 'O', -0.8476, 1))
    spce_unit.append(atom(0+h1x, 0+h1y, 0, 'H', 0.4238, 1))
    spce_unit.append(atom(0+h2x, 0+h2x, 0, 'H', 0.4238, 1))
      
    spce_unit.append(atom(a/2.0, 0, a/2.0, 'O', -0.8476, 2))
    spce_unit.append(atom(a/2.0+h1x, 0+h1y, a/2.0, 'H', 0.4238, 2))
    spce_unit.append(atom(a/2.0+h2x, 0+h2y, a/2.0, 'H', 0.4238, 2))
        
    spce_unit.append(atom(a/2.0, a/2.0, 0, 'O', -0.8476, 3))
    spce_unit.append(atom(a/2.0+h1x, a/2.0+h1y, 0, 'H', 0.4238, 3))
    spce_unit.append(atom(a/2.0+h2x, a/2.0+h2y, 0, 'H', 0.4238, 3))

    spce_unit.append(atom(0, a/2.0, a/2.0, 'O', -0.8276, 4))
    spce_unit.append(atom(0+h1x, a/2.0+h1y, a/2.0, 'H', 0.4238, 4))
    spce_unit.append(atom(0+h2x, a/2.0+h2y, a/2.0, 'H', 0.4238, 4))

    return spce_unit     

def generate_bilattice_nxyz_2(atom_list, bond_list, angle_list, improper_list,
                              dihedral_list, mol_count, 
                              type_index, mass_index, 
                              molecule1, bonds1, angles1, 
                              impropers1, dihedrals1,
                              molecule2, bonds2, angles2, 
                              impropers2, dihedrals2,
                              nx, ny, nz, rho):
    """
    Generates FCC lattice of two molecules with density rho
    Now with added impropers and dihedrals!
    """

    NA = 6.0221415e23

    rho = rho/1e24*NA # density in amu/A^3

    
    m = 0
    
    for a in molecule1:  #Get mass from dictionary
        m += mass_index[type_index[a.spec]]

    for a in molecule2:   
        m += mass_index[type_index[a.spec]]

    count = len(atom_list)
#   print "Mass of molecule = %f" % m
 
    a = (4*m/rho)**(1.0/3.0)

    # Set up 1st FCC cell (molecule on each site)
    fcc_unit1 = []

    for a1 in molecule1:
        fcc_unit1.append(atom(a1.x+0.0,a1.y+0.0,a1.z+0.0,\
                             a1.spec,a1.q,1))


    for a3 in molecule1:
        fcc_unit1.append(atom(a3.x+a/2.0,a3.y+a/2.0,a3.z+0.0,\
                             a3.spec,a3.q,2))

    # Overlap with 2nd FCC cell 
    fcc_unit2 = []

    for a2 in molecule2:
        fcc_unit2.append(atom(a2.x+a/2.0,a2.y+0.0,a2.z+a/2.0,\
                             a2.spec,a2.q,1))

    for a4 in molecule2:
        fcc_unit2.append(atom(a4.x+0.0,a4.y+a/2.0,a4.z+a/2.0,\
                             a4.spec,a4.q,2))

    # Duplicate cell              
    for i in range(0,nx): 
        for j in range(0,ny): 
            for k in range(0,nz):
                for site in fcc_unit1:
                    new_atom = atom(site.x + i*a, site.y + j*a, site.z + k*a,\
                                    site.spec, site.q, mol_count + site.mol)   
                    new_atom.atom = count + 1 
                    atom_list.append(new_atom)
                    count += 1
                mol_count += 2

    for i in range(0,nx): 
        for j in range(0,ny): 
            for k in range(0,nz):
                for site in fcc_unit2:
                    new_atom = atom(site.x + i*a, site.y + j*a, site.z + k*a,\
                                    site.spec, site.q, mol_count + site.mol)   
                    new_atom.atom = count + 1 
                    atom_list.append(new_atom)
                    count += 1
                mol_count += 2
 


    N1 = 2*nx*ny*nz*(len(molecule1)) #total number of created atoms of type1
    N2 = 2*nx*ny*nz*(len(molecule2)) #total number of created atoms of type1
    count = len(atom_list) - (N1 + N2)
  
#   print "N = %d" % N   

    # Add bonds and angles to respective lists
    for i in range(0,N1/(len(molecule1))):
        for c in angles1:
            newangle = angle(c.angle_type, c.atom1+count, c.atom2+count, 
                             c.atom3+count)
            angle_list.append(newangle)
                   
        for b in bonds1:
            newbond = bond(b.bond_type, b.atom1+count, b.atom2+count) 
            bond_list.append(newbond)          

        for d in dihedrals1:
            newdihedral = dihedral(d.dihedral_type, 
                                   d.atomi+count, d.atomj+count, 
                                   d.atomk+count, d.atoml+count)
            dihedral_list.append(newdihedral)                      

        for e in impropers1:
            newimproper = improper(e.improper_type, 
                                   e.atomi+count, e.atomj+count, 
                                   e.atomk+count, e.atoml+count)
            improper_list.append(newimproper)                      

        count += len(molecule1) 

    for i in range(0,N2/(len(molecule2))):
        for c in angles2:
            newangle = angle(c.angle_type, c.atom1+count, c.atom2+count, 
                             c.atom3+count)
            angle_list.append(newangle)
                   
        for b in bonds2:
            newbond = bond(b.bond_type, b.atom1+count, b.atom2+count) 
            bond_list.append(newbond)          

        for d in dihedrals2:
            newdihedral = dihedral(d.dihedral_type, 
                                   d.atomi+count, d.atomj+count, 
                                   d.atomk+count, d.atoml+count)
            dihedral_list.append(newdihedral)                      

        for e in impropers2:
            newimproper = improper(e.improper_type, 
                                   e.atomi+count, e.atomj+count, 
                                   e.atomk+count, e.atoml+count)
            improper_list.append(newimproper)                      

        count += len(molecule2) 

    return atom_list, bond_list, angle_list, \
           improper_list, dihedral_list, mol_count

def generate_bilattice_nxyz(atom_list, bond_list, angle_list, mol_count, \
                            type_index, mass_index, \
                            molecule1, bonds1, angles1,\
                            molecule2, bonds2, angles2,\
                            nx, ny, nz, rho):
    "Generates FCC lattice of two molecules with density rho"
    NA = 6.0221415e23

    rho = rho/1e24*NA # density in amu/A^3

    
    m = 0
    
    for a in molecule1:  #Get mass from dictionary
        m += mass_index[type_index[a.spec]]

    for a in molecule2:   
        m += mass_index[type_index[a.spec]]

#   print "Mass of molecule = %f" % m
 
    a = (4*m/rho)**(1.0/3.0)

    # Set up 1st FCC cell (molecule on each site)
    fcc_unit1 = []

    for a1 in molecule1:
        fcc_unit1.append(atom(a1.x+0.0,a1.y+0.0,a1.z+0.0,\
                             a1.spec,a1.q,1))


    for a3 in molecule1:
        fcc_unit1.append(atom(a3.x+a/2.0,a3.y+a/2.0,a3.z+0.0,\
                             a3.spec,a3.q,3))

    # Overlap with 2nd FCC cell 
    fcc_unit2 = []

    for a2 in molecule2:
        fcc_unit2.append(atom(a2.x+a/2.0,a2.y+0.0,a2.z+a/2.0,\
                             a2.spec,a2.q,2))

    for a4 in molecule2:
        fcc_unit2.append(atom(a4.x+0.0,a4.y+a/2.0,a4.z+a/2.0,\
                             a4.spec,a4.q,4))

    # Duplicate cell              
    for i in range(0,nx): 
        for j in range(0,ny): 
            for k in range(0,nz):
                for site in fcc_unit1:
                    new_atom = atom(site.x + i*a, site.y + j*a, site.z + k*a,\
                                    site.spec, site.q, mol_count + site.mol)   
                    atom_list.append(new_atom)
                mol_count += 2

    for i in range(0,nx): 
        for j in range(0,ny): 
            for k in range(0,nz):
                for site in fcc_unit2:
                    new_atom = atom(site.x + i*a, site.y + j*a, site.z + k*a,\
                                    site.spec, site.q, mol_count + site.mol)   
                    atom_list.append(new_atom)
                mol_count += 2
 


    N1 = 2*nx*ny*nz*(len(molecule1)) #total number of created atoms of type1
    N2 = 2*nx*ny*nz*(len(molecule2)) #total number of created atoms of type1
    count = len(atom_list) - (N1 + N2)
  
#   print "N = %d" % N   

    # Add bonds and angles to respective lists
    for i in range(0,N1/(len(molecule1))):
        for c in angles1:
            newangle = angle(c.angle_type, c.atom1+count, c.atom2+count, 
                             c.atom3+count)
            angle_list.append(newangle)
                   
        for b in bonds1:
            newbond = bond(b.bond_type, b.atom1+count, b.atom2+count) 
            bond_list.append(newbond)          

        count += len(molecule1) 

    for i in range(0,N2/(len(molecule2))):
        for c in angles2:
            newangle = angle(c.angle_type, c.atom1+count, c.atom2+count, 
                             c.atom3+count)
            angle_list.append(newangle)
                   
        for b in bonds2:
            newbond = bond(b.bond_type, b.atom1+count, b.atom2+count) 
            bond_list.append(newbond)          

        count += len(molecule2) 

    return atom_list, bond_list, angle_list, mol_count

def generate_lattice_nxyz(atom_list, bond_list, angle_list, mol_count, \
                     type_index, mass_index, molecule, bonds, angles,\
                     nx, ny, nz, rho):
    "Generates FCC lattice of desired molecule with density rho"
    NA = 6.0221415e23

    rho = rho/1e24*NA # density in amu/A^3

    
    m = 0
    
    for a in molecule:  #Get mass from dictionary
        m += mass_index[type_index[a.spec]]
   
    print "Mass of molecule = %f" % m
 
    a = (4*m/rho)**(1.0/3.0)

    # Set up FCC cell (molecule on each site)
    fcc_unit = []
    
    for a1 in molecule:
        fcc_unit.append(atom(a1.x+0.0,a1.y+0.0,a1.z+0.0,\
                             a1.spec,a1.q,1))

    for a2 in molecule:
        fcc_unit.append(atom(a2.x+a/2.0,a2.y+0.0,a2.z+a/2.0,\
                             a2.spec,a2.q,2))

    for a3 in molecule:
        fcc_unit.append(atom(a3.x+a/2.0,a3.y+a/2.0,a3.z+0.0,\
                             a3.spec,a3.q,3))

    for a3 in molecule:
        fcc_unit.append(atom(a3.x+0.0,a3.y+a/2.0,a3.z+a/2.0,\
                             a3.spec,a3.q,4))



    # Duplicate cell              
    for i in range(0,nx): 
        for j in range(0,ny): 
            for k in range(0,nz):
                for site in fcc_unit:
                    new_atom = atom(site.x + i*a, site.y + j*a, site.z + k*a,\
                                    site.spec, site.q, mol_count + site.mol)   
                    atom_list.append(new_atom)
                mol_count += 4
 
    N = 4*nx*ny*nz*len(molecule) # total number of created atoms
    count = len(atom_list) - N
  
    print "N = %d" % N   

    # Add bonds and angles to respective lists
    for i in range(0,N/len(molecule)):
        for c in angles:
            newangle = angle(c.angle_type, c.atom1+count, c.atom2+count, 
                             c.atom3+count)
            angle_list.append(newangle)
                   
        for b in bonds:
            newbond = bond(b.bond_type, b.atom1+count, b.atom2+count) 
            bond_list.append(newbond)          

        count += len(molecule) 

    return atom_list, bond_list, angle_list, mol_count

def generate_lattice(atom_list, bond_list, angle_list, mol_count, \
                     type_index, mass_index, molecule, bonds, angles,\
                     xlo, xhi, ylo, yhi, zlo, zhi, rho):
    "Generates FCC lattice in specified region with density closest to rho"
    NA = 6.0221415e23

    rho = rho/1e24*NA # density in amu/A^3

    lx = (xhi - xlo)
    ly = (yhi - ylo)
    lz = (zhi - zlo)

    m = 0
    
    for a in molecule:  #Get mass from dictionary
        m += mass_index[type_index[a.spec]]
   
    print "Mass of molecule = %f" % m
 
    a = (4*m/rho)**(1.0/3.0)
    nx = int(lx/a)            
    ny = int(ly/a)            
    nz = int(lz/a)

    print "area = %f\n" % (lx*ly)
    print "volume = %f\n" % (lx*ly*lz)

    Nm = 4*nx*ny*nz

    rho_real = Nm*m/(lx*ly*lz) # actual density after lattice
    rho_real = rho_real/NA*1e24

    print "Density in region is %f g/cm^3" % rho_real   

    # Set up FCC cell (molecule on each site)
    fcc_unit = []
    
    for a1 in molecule:
        fcc_unit.append(atom(a1.x+0.0+xlo,a1.y+0.0+ylo,a1.z+0.0+zlo,\
                             a1.spec,a1.q,1))

    for a2 in molecule:
        fcc_unit.append(atom(a2.x+a/2.0+xlo,a2.y+0.0+ylo,a2.z+a/2.0+zlo,\
                             a2.spec,a2.q,2))

    for a3 in molecule:
        fcc_unit.append(atom(a3.x+a/2.0+xlo,a3.y+a/2.0+ylo,a3.z+0.0+zlo,\
                             a3.spec,a3.q,3))

    for a3 in molecule:
        fcc_unit.append(atom(a3.x+0.0+xlo,a3.y+a/2.0+ylo,a3.z+a/2.0+zlo,\
                             a3.spec,a3.q,4))

    # Duplicate cell              
    for i in range(0,nx): 
        for j in range(0,ny): 
            for k in range(0,nz):
                for site in fcc_unit:
                    new_atom = atom(site.x + i*a, site.y + j*a, site.z + k*a,\
                                    site.spec, site.q, mol_count + site.mol)   
                    atom_list.append(new_atom)
                mol_count += 4
 
    N = 4*nx*ny*nz*len(molecule) # total number of created atoms
    count = len(atom_list) - N
  
    print "N = %d" % N   

    # Add bonds and angles to respective lists
    for i in range(0,N/len(molecule)):
        for c in angles:
            newangle = angle(c.angle_type, c.atom1+count, c.atom2+count, 
                             c.atom3+count)
            angle_list.append(newangle)
                   
        for b in bonds:
            newbond = bond(b.bond_type, b.atom1+count, b.atom2+count) 
            bond_list.append(newbond)          

        count += len(molecule) 

    return atom_list, bond_list, angle_list, mol_count, Nm, m

def generate_fcc_lattice_2(atom_list, bond_list, angle_list, improper_list,
                           dihedral_list, mol_count, 
                           type_index, mass_index, molecule, 
                           bonds, angles, impropers, dihedrals,
                           nx, ny, nz, x_off, y_off, z_off, a):
    """"
    Generates FCC lattice in 3 dimensions starting at point 
    (x_off, y_off, z_off) repeated for nx, ny, nz units. Cell size is 
    specified by lattice constant a. Now includes impropers and dihedrals.
    """ 

    Lx = nx*a  
    Ly = ny*a  
    Lz = nz*a  
 
    print "volume = %f\n" % (Lx*Ly*Lz)

    # Set up FCC cell (molecule on each site)
    fcc_unit = []
    
    for a1 in molecule:
        fcc_unit.append(atom(a1.x+0.0,a1.y+0.0,a1.z+0.0,\
                             a1.spec,a1.q,1))

    for a2 in molecule:
        fcc_unit.append(atom(a2.x+a/2.0,a2.y+0.0,a2.z+a/2.0,\
                             a2.spec,a2.q,2))

    for a3 in molecule:
        fcc_unit.append(atom(a3.x+a/2.0,a3.y+a/2.0,a3.z+0.0,\
                             a3.spec,a3.q,3))

    for a3 in molecule:
        fcc_unit.append(atom(a3.x+0.0,a3.y+a/2.0,a3.z+a/2.0,\
                             a3.spec,a3.q,4))

    # Duplicate cell              
    for i in range(0,nx): 
        for j in range(0,ny): 
            for k in range(0,nz):
                for site in fcc_unit:
                    x_set = site.x + i*a + x_off
                    y_set = site.y + j*a + y_off
                    z_set = site.z + k*a + z_off
                    new_atom = atom(x_set, y_set, z_set,
                                    site.spec, site.q, mol_count + site.mol)   
                    atom_list.append(new_atom)
                mol_count += 4
 
    N = 4*nx*ny*nz*len(molecule) # total number of created atoms
    count = len(atom_list) - N
  
    print "N = %d" % N   

    # Add bonds and angles to respective lists
    for i in range(0,N/len(molecule)):
        for c in angles:
            newangle = angle(c.angle_type, c.atom1+count, c.atom2+count, 
                             c.atom3+count)
            angle_list.append(newangle)
                   
        for b in bonds:
            newbond = bond(b.bond_type, b.atom1+count, b.atom2+count) 
            bond_list.append(newbond)          

        for d in dihedrals:
            newdihedral = dihedral(d.dihedral_type,
                                   d.atomi+count, d.atomj+count,
                                   d.atomk+count, d.atoml+count)
            dihedral_list.append(newdihedral)

        for e in impropers:
            newimproper = improper(e.improper_type,
                                   e.atomi+count, e.atomj+count,
                                   e.atomk+count, e.atoml+count)
            improper_list.append(newimproper)


        count += len(molecule) 

    return atom_list, bond_list, angle_list, \
           improper_list, dihedral_list, mol_count

def generate_fcc_lattice(atom_list, bond_list, angle_list, mol_count, \
                     type_index, mass_index, molecule, bonds, angles,\
                     nx, ny, nz, x_off, y_off, z_off, a):
    """"
    Generates FCC lattice in 3 dimensions starting at point 
    (x_off, y_off, z_off) repeated for nx, ny, nz units. Cell size is 
    specified by lattice constant a.
    """ 

    Lx = nx*a  
    Ly = ny*a  
    Lz = nz*a  
 
    print "volume = %f\n" % (Lx*Ly*Lz)

    # Set up FCC cell (molecule on each site)
    fcc_unit = []
    
    for a1 in molecule:
        fcc_unit.append(atom(a1.x+0.0,a1.y+0.0,a1.z+0.0,\
                             a1.spec,a1.q,1))

    for a2 in molecule:
        fcc_unit.append(atom(a2.x+a/2.0,a2.y+0.0,a2.z+a/2.0,\
                             a2.spec,a2.q,2))

    for a3 in molecule:
        fcc_unit.append(atom(a3.x+a/2.0,a3.y+a/2.0,a3.z+0.0,\
                             a3.spec,a3.q,3))

    for a3 in molecule:
        fcc_unit.append(atom(a3.x+0.0,a3.y+a/2.0,a3.z+a/2.0,\
                             a3.spec,a3.q,4))

    # Duplicate cell              
    for i in range(0,nx): 
        for j in range(0,ny): 
            for k in range(0,nz):
                for site in fcc_unit:
                    x_set = site.x + i*a + x_off
                    y_set = site.y + j*a + y_off
                    z_set = site.z + k*a + z_off
                    new_atom = atom(x_set, y_set, z_set,\
                                    site.spec, site.q, mol_count + site.mol)   
                    atom_list.append(new_atom)
                mol_count += 4
 
    N = 4*nx*ny*nz*len(molecule) # total number of created atoms
    count = len(atom_list) - N
  
    print "N = %d" % N   

    # Add bonds and angles to respective lists
    for i in range(0,N/len(molecule)):
        for c in angles:
            newangle = angle(c.angle_type, c.atom1+count, c.atom2+count, 
                             c.atom3+count)
            angle_list.append(newangle)
                   
        for b in bonds:
            newbond = bond(b.bond_type, b.atom1+count, b.atom2+count) 
            bond_list.append(newbond)          

        count += len(molecule) 

    return atom_list, bond_list, angle_list, mol_count

def spce_unit_cell(a):
    "Returns fcc unit cell for spc/e water. Takes lattice constant argument"
    pi = 4.0*math.atan(1.0)

    ang = 109.47*pi/180 # bond angle in radians
    l = 1.0             # bond length in Angstrom

    # distances for hydrogen atoms wrt to oxygen 
    h1x = l*math.cos(ang/2.0)
    h1y = l*math.sin(ang/2.0)
    h2x = h1x
    h2y = -h2x

    spce_unit = []

    spce_unit.append(atom(0, 0, 0, 'O', -0.8476, 1))
    spce_unit.append(atom(0+h1x, 0+h1y, 0, 'H', 0.4238, 1))
    spce_unit.append(atom(0+h2x, 0+h2x, 0, 'H', 0.4238, 1))
      
    spce_unit.append(atom(a/2.0, 0, a/2.0, 'O', -0.8476, 2))
    spce_unit.append(atom(a/2.0+h1x, 0+h1y, a/2.0, 'H', 0.4238, 2))
    spce_unit.append(atom(a/2.0+h2x, 0+h2y, a/2.0, 'H', 0.4238, 2))
        
    spce_unit.append(atom(a/2.0, a/2.0, 0, 'O', -0.8476, 3))
    spce_unit.append(atom(a/2.0+h1x, a/2.0+h1y, 0, 'H', 0.4238, 3))
    spce_unit.append(atom(a/2.0+h2x, a/2.0+h2y, 0, 'H', 0.4238, 3))

    spce_unit.append(atom(0, a/2.0, a/2.0, 'O', -0.8276, 4))
    spce_unit.append(atom(0+h1x, a/2.0+h1y, a/2.0, 'H', 0.4238, 4))
    spce_unit.append(atom(0+h2x, a/2.0+h2y, a/2.0, 'H', 0.4238, 4))

    return spce_unit     

def SAFT_MPD_molecule_4_SHAKE(type_list, bond_types, 
                              angle_types, dihedral_types):
    """
    Returns list of atoms, bonds, angles and dihedrals corresponding 
    to meta-pheyline diamine (MPD) molecule based on
    SAFT 3-segment benzene molecule, from 
    T. Lafitte et al. Molecular Physics, 
    Vol. 110, Nos. 11-12, 10-20 June 2012, 1189-1203 

    ** This version also includes dihedral interaction and correct
       bond distances from DFT structure 

    ** Version 4 contains two sticky sites labelled such that 
       correct adhesion to TMC is formed 
    
    ** SHAKE version includes two bonds and one angle per benzene triangle
       so that SHAKE can be used to retain molecular geometry instead of 
       the rigid integrator in lammps 

    1:MPD_B, 2:MPD_B, 3:MPD_B, 4:MPD_N, 5:MPD_N, 6:MPD_SA, 7:MPD_SB,
    8:MPD_SA, 9:MPD_SB
     
    Defines bond type 1, 2, 3, 4
            angle type 1, 2, 3, 4
            dihedral type 1
    """ 

    pi = 4.0*math.atan(1.0)

    ang = 60*pi/180     # bond angle in radians
    sigma = 2.42           # bond length in Angstrom

    d = 0.56 # Offset parameter for sigma_SS = 1 Angstrom
    theta = 30*pi/180 # angle of sticky sites wrt bond direction


    sigma_BB = 2.42   # bond between SAFT benzene beads in Angstrom
    sigma_BN = 1.42   # bond between SAFT benzene beads and carbon in Angstrom
    sigma_NS = 0.69   # bond between carbon and sticky site in Angstrom

    #sigma_CS -= 0.784

    saft_mpd = []
    saft_mpd_bonds = []
    saft_mpd_angles = []
    saft_mpd_dihedrals = []

    saft_mpd_types = ['MPD_B', 'MPD_N', 'MPD_SA', 'MPD_SB']

    angles = [1, 2, 3, 4]
    bonds = [1, 2, 3, 4]
    dihedrals = [1]
            
    type_list.extend(saft_mpd_types)

    rt3 = math.sqrt(3.0)

    B1x = 0.0
    B1y = 0.0
    B2x = sigma_BB
    B2y = 0.0
    B3x = 0.5*sigma_BB
    B3y = rt3/2.0*sigma_BB

    N1x = -rt3/2.0*sigma_BN
    N1y = -0.5*sigma_BN
    N2x = sigma_BB + rt3/2.0*sigma_BN
    N2y = -0.5*sigma_BN

#   S1x = -rt3/2.0*sigma_BN
#   S1y = -0.5*sigma_BN - sigma_NS
#   S2x = sigma_BB + rt3/2.0*sigma_BN + rt3/2.0*sigma_NS
#   S2y = -0.5*sigma_BN + 0.5*sigma_NS

    S1Ax = -rt3/2.0*sigma_BN - sigma_NS*math.sin(theta)
    S1Ay = -0.5*sigma_BN - sigma_NS*math.cos(theta)
    S1Bx = -rt3/2.0*sigma_BN + sigma_NS*math.sin(theta)
    S1By = -0.5*sigma_BN - sigma_NS*math.cos(theta)

    S2Ax = sigma_BB + rt3/2.0*sigma_BN + sigma_NS*math.cos(pi/6 - theta)
    S2Ay = -0.5*sigma_BN + sigma_NS*math.sin(pi/6 - theta)
    S2Bx = sigma_BB + rt3/2.0*sigma_BN + sigma_NS*math.cos(pi/6 + theta)
    S2By = -0.5*sigma_BN + sigma_NS*math.sin(pi/6 + theta)


    saft_mpd.append(atom(B1x, B1y, 0, 'MPD_B', 0.0, 1))  #1
    saft_mpd.append(atom(B2x, B2y, 0, 'MPD_B', 0.0, 1))  #2
    saft_mpd.append(atom(B3x, B3y, 0, 'MPD_B', 0.0, 1))  #3

    saft_mpd.append(atom(N1x, N1y, 0, 'MPD_N', 0.0, 1))  #4
    saft_mpd.append(atom(N2x, N2y, 0, 'MPD_N', 0.0, 1))  #5

    saft_mpd.append(atom(S1Ax, S1Ay, 0, 'MPD_SA', 0.0, 1))  #6
    saft_mpd.append(atom(S1Bx, S1By, 0, 'MPD_SB', 0.0, 1))  #7
    saft_mpd.append(atom(S2Ax, S2Ay, 0, 'MPD_SA', 0.0, 1))  #8
    saft_mpd.append(atom(S2Bx, S2By, 0, 'MPD_SB', 0.0, 1))  #9

    saft_mpd_bonds.append(bond(1, 1, 2))
    saft_mpd_bonds.append(bond(1, 2, 3))
#   saft_mpd_bonds.append(bond(1, 1, 3)) removed for SHAKE
        
    saft_mpd_bonds.append(bond(2, 1, 4))
    saft_mpd_bonds.append(bond(2, 2, 5))
 
    saft_mpd_bonds.append(bond(3, 4, 6))
    saft_mpd_bonds.append(bond(4, 4, 7))

    saft_mpd_bonds.append(bond(3, 5, 8))
    saft_mpd_bonds.append(bond(4, 5, 9))

#   saft_mpd_angles.append(angle(1, 3, 1, 2)) removed for SHAKE
#   saft_mpd_angles.append(angle(1, 2, 3, 1)) removed for SHAKE
    saft_mpd_angles.append(angle(1, 1, 2, 3))
    
    saft_mpd_angles.append(angle(2, 3, 1, 4))
    saft_mpd_angles.append(angle(2, 1, 2, 5))

    saft_mpd_angles.append(angle(3, 1, 4, 6)) # short angle
    saft_mpd_angles.append(angle(4, 1, 4, 7)) # wide angle
    saft_mpd_angles.append(angle(3, 2, 5, 8))
    saft_mpd_angles.append(angle(4, 2, 5, 9))

    saft_mpd_dihedrals.append(dihedral(1, 3, 1, 4, 6))
    saft_mpd_dihedrals.append(dihedral(1, 3, 1, 4, 7))
    saft_mpd_dihedrals.append(dihedral(1, 1, 2, 5, 8))
    saft_mpd_dihedrals.append(dihedral(1, 1, 2, 5, 9))

    bond_types.extend(bonds)
    angle_types.extend(angles)
    dihedral_types.extend(dihedrals)

    return saft_mpd, saft_mpd_bonds, saft_mpd_angles, saft_mpd_dihedrals, \
            type_list, bond_types, angle_types, dihedral_types

def SAFT_TMC_molecule_4_SHAKE(type_list, bond_types, 
                              angle_types, dihedral_types):
    """
    Returns list of atoms, bonds, angles and dihgedrals corresponding 
    to Benzene-1,3,5-tricarboxylicphenylamide (TMC) molecule based on
    SAFT 3-segment benzene molecule, from 
    T. Lafitte et al. Molecular Physics, 
    Vol. 110, Nos. 11-12, 10-20 June 2012, 1189-1203 

    ** This version also includes dihedral interaction and correct
       bond distances from DFT structure 

    ** Version 4 contains two sticky sites labelled such that 
       correct adhesion to MPD is formed 
    
    ** SHAKE version includes two bonds and one angle per benzene triangle
       so that SHAKE can be used to retain molecular geometry instead of 
       the rigid integrator in lammps 
    
    1:TMC_B, 2:TMC_B, 3:TMC_B, 4:TMC_C, 5:TMC_C, 6:TMC_C,
    7:TMS_SA, 8:TMC_SB, 9:TMC_SA, 10:TMC_SB, 11:TMC_SA, 12:TMC_SB  
    Defines bond type 1, 5, 6, 7
            angle type 1, 5, 6, 7 
            dihedral type 1
    """ 

    pi = 4.0*math.atan(1.0)

    ang = 60*pi/180     # bond angle in radians
    sigma = 2.42           # bond length in Angstrom

    d = 0.56 # Offset parameter for sigma_SS = 1 Angstrom

    sigma_BB = 2.42   # bond between SAFT benzene beads in Angstrom
    sigma_BC = 1.51   # bond between SAFT benzene beads and carbon in Angstrom
    sigma_CS = 0.69   # bond between carbon and sticky site in Angstrom

    theta = 30*pi/180 # angle of sticky sites wrt bond direction

    #sigma_CS -= 0.784

    saft_tmc = []
    saft_tmc_bonds = []
    saft_tmc_angles = []
    saft_tmc_dihedrals = []

    saft_tmc_types = ['TMC_B', 'TMC_C', 'TMC_SA', 'TMC_SB']

    angles = [1, 5, 6, 7]
    bonds = [1, 5, 6, 7]
    dihedrals = [1]
            
    type_list.extend(saft_tmc_types)

    rt3 = math.sqrt(3.0)

    B1x = 0.0
    B1y = 0.0
    B2x = sigma_BB
    B2y = 0.0
    B3x = 0.5*sigma_BB
    B3y = rt3/2.0*sigma_BB

    C1x = -rt3/2.0*sigma_BC
    C1y = -0.5*sigma_BC
    C2x = sigma_BB + rt3/2.0*sigma_BC
    C2y = -0.5*sigma_BC
    C3x = 0.5*sigma_BB
    C3y = rt3/2.0*sigma_BB + sigma_BC

#   S1x = -rt3/2.0*sigma_BC
#   S1y = -0.5*sigma_BC - sigma_CS
#   S2x = sigma_BB + rt3/2.0*sigma_BC + rt3/2.0*sigma_CS
#   S2y = -0.5*sigma_BC + 0.5*sigma_CS
#   S3x = 0.5*sigma_BB - rt3/2.0*sigma_CS
#   S3y = rt3/2.0*sigma_BB + sigma_BC + 0.5*sigma_CS

    S1Bx = -rt3/2.0*sigma_BC - sigma_CS*math.sin(theta)
    S1By = -0.5*sigma_BC - sigma_CS*math.cos(theta)
    S1Ax = -rt3/2.0*sigma_BC + sigma_CS*math.sin(theta)
    S1Ay = -0.5*sigma_BC - sigma_CS*math.cos(theta)

    S2Bx = sigma_BB + rt3/2.0*sigma_BC + sigma_CS*math.cos(pi/6 - theta)
    S2By = -0.5*sigma_BC + sigma_CS*math.sin(pi/6 - theta)
    S2Ax = sigma_BB + rt3/2.0*sigma_BC + sigma_CS*math.cos(pi/6 + theta)
    S2Ay = -0.5*sigma_BC + sigma_CS*math.sin(pi/6 + theta)

    S3Bx = 0.5*sigma_BB - sigma_CS*math.cos(pi/6 + theta)
    S3By = rt3/2.0*sigma_BB + sigma_BC + sigma_CS*math.sin(pi/6 + theta)
    S3Ax = 0.5*sigma_BB - sigma_CS*math.cos(pi/6 - theta)
    S3Ay = rt3/2.0*sigma_BB + sigma_BC + sigma_CS*math.sin(pi/6 - theta)


    saft_tmc.append(atom(B1x, B1y, 0, 'TMC_B', 0.0, 1))  #1
    saft_tmc.append(atom(B2x, B2y, 0, 'TMC_B', 0.0, 1))  #2
    saft_tmc.append(atom(B3x, B3y, 0, 'TMC_B', 0.0, 1))  #3

    saft_tmc.append(atom(C1x, C1y, 0, 'TMC_C', 0.0, 1))  #4
    saft_tmc.append(atom(C2x, C2y, 0, 'TMC_C', 0.0, 1))  #5
    saft_tmc.append(atom(C3x, C3y, 0, 'TMC_C', 0.0, 1))  #6

    saft_tmc.append(atom(S1Ax, S1Ay, 0, 'TMC_SA', 0.0, 1))  #7
    saft_tmc.append(atom(S1Bx, S1By, 0, 'TMC_SB', 0.0, 1))  #8
    saft_tmc.append(atom(S2Ax, S2Ay, 0, 'TMC_SA', 0.0, 1))  #9
    saft_tmc.append(atom(S2Bx, S2By, 0, 'TMC_SB', 0.0, 1))  #10
    saft_tmc.append(atom(S3Ax, S3Ay, 0, 'TMC_SA', 0.0, 1))  #11
    saft_tmc.append(atom(S3Bx, S3By, 0, 'TMC_SB', 0.0, 1))  #12

    saft_tmc_bonds.append(bond(1, 1, 2))
    saft_tmc_bonds.append(bond(1, 2, 3))
#   saft_tmc_bonds.append(bond(1, 1, 3)) removed for SHAKE
        
    saft_tmc_bonds.append(bond(5, 1, 4))
    saft_tmc_bonds.append(bond(5, 2, 5))
    saft_tmc_bonds.append(bond(5, 3, 6))
 
    saft_tmc_bonds.append(bond(6, 4, 7))
    saft_tmc_bonds.append(bond(7, 4, 8))
    saft_tmc_bonds.append(bond(6, 5, 9))
    saft_tmc_bonds.append(bond(7, 5, 10))
    saft_tmc_bonds.append(bond(6, 6, 11))
    saft_tmc_bonds.append(bond(7, 6, 12))

#   saft_tmc_angles.append(angle(1, 3, 1, 2)) removed for SHAKE
#   saft_tmc_angles.append(angle(1, 2, 3, 1)) removed for SHAKE
    saft_tmc_angles.append(angle(1, 1, 2, 3))
    
    saft_tmc_angles.append(angle(5, 4, 1, 3))
    saft_tmc_angles.append(angle(5, 1, 2, 5))
    saft_tmc_angles.append(angle(5, 6, 3, 1))

    saft_tmc_angles.append(angle(6, 1, 4, 7)) # short angle A
    saft_tmc_angles.append(angle(7, 1, 4, 8)) # wide angle B
    saft_tmc_angles.append(angle(6, 2, 5, 9))
    saft_tmc_angles.append(angle(7, 2, 5, 10))
    saft_tmc_angles.append(angle(6, 3, 6, 11))
    saft_tmc_angles.append(angle(7, 3, 6, 12))

    saft_tmc_dihedrals.append(dihedral(1, 3, 1, 4, 7))
    saft_tmc_dihedrals.append(dihedral(1, 3, 1, 4, 8))
    saft_tmc_dihedrals.append(dihedral(1, 1, 2, 5, 9))
    saft_tmc_dihedrals.append(dihedral(1, 1, 2, 5, 10))
    saft_tmc_dihedrals.append(dihedral(1, 2, 3, 6, 11))
    saft_tmc_dihedrals.append(dihedral(1, 2, 3, 6, 12))

    bond_types.extend(bonds)
    angle_types.extend(angles)
    dihedral_types.extend(dihedrals)

    return saft_tmc, saft_tmc_bonds, saft_tmc_angles, saft_tmc_dihedrals,\
           type_list, bond_types, angle_types, dihedral_types

def SAFT_MPD_molecule_4_imp(type_list, bond_types, 
                            angle_types, dihedral_types, improper_types):
    """
    Returns list of atoms, bonds, angles and dihedrals corresponding 
    to meta-pheyline diamine (MPD) molecule based on
    SAFT 3-segment benzene molecule, from 
    T. Lafitte et al. Molecular Physics, 
    Vol. 110, Nos. 11-12, 10-20 June 2012, 1189-1203 

    ** This version also includes dihedral interaction and correct
       bond distances from DFT structure 

    ** Version 4 contains two sticky sites labelled such that 
       correct adhesion to TMC is formed 

    ** This version now includes improper dihedrals so that the structure 
       can be maintained without the assumption of rigid bodies
    
    1:MPD_B, 2:MPD_B, 3:MPD_B, 4:MPD_N, 5:MPD_N, 6:MPD_SA, 7:MPD_SB,
    8:MPD_SA, 9:MPD_SB
     
    Defines bond type 1, 2, 3, 4
            angle type 1, 2, 3, 4
            dihedral type 1
            improper type 1 
    """ 

    pi = 4.0*math.atan(1.0)

    ang = 60*pi/180     # bond angle in radians
    sigma = 2.42           # bond length in Angstrom

    d = 0.56 # Offset parameter for sigma_SS = 1 Angstrom
    theta = 30*pi/180 # angle of sticky sites wrt bond direction


    sigma_BB = 2.42   # bond between SAFT benzene beads in Angstrom
    sigma_BN = 1.42   # bond between SAFT benzene beads and carbon in Angstrom
    sigma_NS = 0.69   # bond between carbon and sticky site in Angstrom

    #sigma_CS -= 0.784

    saft_mpd = []
    saft_mpd_bonds = []
    saft_mpd_angles = []
    saft_mpd_dihedrals = []
    saft_mpd_impropers = []

    saft_mpd_types = ['MPD_B', 'MPD_N', 'MPD_SA', 'MPD_SB']

    angles = [1, 2, 3, 4]
    bonds = [1, 2, 3, 4]
    dihedrals = [1]
    impropers = [1]
            
    type_list.extend(saft_mpd_types)

    rt3 = math.sqrt(3.0)

    B1x = 0.0
    B1y = 0.0
    B2x = sigma_BB
    B2y = 0.0
    B3x = 0.5*sigma_BB
    B3y = rt3/2.0*sigma_BB

    N1x = -rt3/2.0*sigma_BN
    N1y = -0.5*sigma_BN
    N2x = sigma_BB + rt3/2.0*sigma_BN
    N2y = -0.5*sigma_BN

#   S1x = -rt3/2.0*sigma_BN
#   S1y = -0.5*sigma_BN - sigma_NS
#   S2x = sigma_BB + rt3/2.0*sigma_BN + rt3/2.0*sigma_NS
#   S2y = -0.5*sigma_BN + 0.5*sigma_NS

    S1Ax = -rt3/2.0*sigma_BN - sigma_NS*math.sin(theta)
    S1Ay = -0.5*sigma_BN - sigma_NS*math.cos(theta)
    S1Bx = -rt3/2.0*sigma_BN + sigma_NS*math.sin(theta)
    S1By = -0.5*sigma_BN - sigma_NS*math.cos(theta)

    S2Ax = sigma_BB + rt3/2.0*sigma_BN + sigma_NS*math.cos(pi/6 - theta)
    S2Ay = -0.5*sigma_BN + sigma_NS*math.sin(pi/6 - theta)
    S2Bx = sigma_BB + rt3/2.0*sigma_BN + sigma_NS*math.cos(pi/6 + theta)
    S2By = -0.5*sigma_BN + sigma_NS*math.sin(pi/6 + theta)


    saft_mpd.append(atom(B1x, B1y, 0, 'MPD_B', 0.0, 1))  #1
    saft_mpd.append(atom(B2x, B2y, 0, 'MPD_B', 0.0, 1))  #2
    saft_mpd.append(atom(B3x, B3y, 0, 'MPD_B', 0.0, 1))  #3

    saft_mpd.append(atom(N1x, N1y, 0, 'MPD_N', 0.0, 1))  #4
    saft_mpd.append(atom(N2x, N2y, 0, 'MPD_N', 0.0, 1))  #5

    saft_mpd.append(atom(S1Ax, S1Ay, 0, 'MPD_SA', 0.0, 1))  #6
    saft_mpd.append(atom(S1Bx, S1By, 0, 'MPD_SB', 0.0, 1))  #7
    saft_mpd.append(atom(S2Ax, S2Ay, 0, 'MPD_SA', 0.0, 1))  #8
    saft_mpd.append(atom(S2Bx, S2By, 0, 'MPD_SB', 0.0, 1))  #9

    saft_mpd_bonds.append(bond(1, 1, 2))
    saft_mpd_bonds.append(bond(1, 2, 3))
    saft_mpd_bonds.append(bond(1, 1, 3))
        
    saft_mpd_bonds.append(bond(2, 1, 4))
    saft_mpd_bonds.append(bond(2, 2, 5))
 
    saft_mpd_bonds.append(bond(3, 4, 6))
    saft_mpd_bonds.append(bond(4, 4, 7))

    saft_mpd_bonds.append(bond(3, 5, 8))
    saft_mpd_bonds.append(bond(4, 5, 9))

    saft_mpd_angles.append(angle(1, 3, 1, 2))
    saft_mpd_angles.append(angle(1, 2, 3, 1))
    saft_mpd_angles.append(angle(1, 1, 2, 3))
    
    saft_mpd_angles.append(angle(2, 3, 1, 4))
    saft_mpd_angles.append(angle(2, 1, 2, 5))

    saft_mpd_angles.append(angle(3, 1, 4, 6)) # short angle
    saft_mpd_angles.append(angle(4, 1, 4, 7)) # wide angle
    saft_mpd_angles.append(angle(3, 2, 5, 8))
    saft_mpd_angles.append(angle(4, 2, 5, 9))

    saft_mpd_dihedrals.append(dihedral(1, 3, 1, 4, 6))
    saft_mpd_dihedrals.append(dihedral(1, 3, 1, 4, 7))
    saft_mpd_dihedrals.append(dihedral(1, 1, 2, 5, 8))
    saft_mpd_dihedrals.append(dihedral(1, 1, 2, 5, 9))

    saft_mpd_impropers.append(improper(1, 4, 2, 3, 1))
    saft_mpd_impropers.append(improper(1, 5, 1, 3, 2))

    bond_types.extend(bonds)
    angle_types.extend(angles)
    dihedral_types.extend(dihedrals)
    improper_types.extend(impropers)

    return saft_mpd, saft_mpd_bonds, saft_mpd_angles, saft_mpd_dihedrals, \
           saft_mpd_impropers, type_list, bond_types, angle_types, \
           dihedral_types, improper_types

def SAFT_TMC_molecule_4_imp(type_list, bond_types, angle_types, 
                            dihedral_types, improper_types):
    """
    Returns list of atoms, bonds, angles and dihgedrals corresponding 
    to Benzene-1,3,5-tricarboxylicphenylamide (TMC) molecule based on
    SAFT 3-segment benzene molecule, from 
    T. Lafitte et al. Molecular Physics, 
    Vol. 110, Nos. 11-12, 10-20 June 2012, 1189-1203 

    ** This version also includes dihedral interaction and correct
       bond distances from DFT structure 

    ** Version 4 contains two sticky sites labelled such that 
       correct adhesion to MPD is formed 
    
    ** This version now includes improper dihedrals so that the structure 
       can be maintained without the assumption of rigid bodies

    1:TMC_B, 2:TMC_B, 3:TMC_B, 4:TMC_C, 5:TMC_C, 6:TMC_C,
    7:TMS_SA, 8:TMC_SB, 9:TMC_SA, 10:TMC_SB, 11:TMC_SA, 12:TMC_SB  
 
   Defines bond type 1, 5, 6, 7
            angle type 1, 5, 6, 7 
            dihedral type 1
            improper type 1   
    """ 

    pi = 4.0*math.atan(1.0)

    ang = 60*pi/180     # bond angle in radians
    sigma = 2.42           # bond length in Angstrom

    d = 0.56 # Offset parameter for sigma_SS = 1 Angstrom

    sigma_BB = 2.42   # bond between SAFT benzene beads in Angstrom
    sigma_BC = 1.51   # bond between SAFT benzene beads and carbon in Angstrom
    sigma_CS = 0.69   # bond between carbon and sticky site in Angstrom

    theta = 30*pi/180 # angle of sticky sites wrt bond direction

    #sigma_CS -= 0.784

    saft_tmc = []
    saft_tmc_bonds = []
    saft_tmc_angles = []
    saft_tmc_dihedrals = []
    saft_tmc_impropers = []

    saft_tmc_types = ['TMC_B', 'TMC_C', 'TMC_SA', 'TMC_SB']

    angles = [1, 5, 6, 7]
    bonds = [1, 5, 6, 7]
    dihedrals = [1]
    impropers = [1]
            
    type_list.extend(saft_tmc_types)

    rt3 = math.sqrt(3.0)

    B1x = 0.0
    B1y = 0.0
    B2x = sigma_BB
    B2y = 0.0
    B3x = 0.5*sigma_BB
    B3y = rt3/2.0*sigma_BB

    C1x = -rt3/2.0*sigma_BC
    C1y = -0.5*sigma_BC
    C2x = sigma_BB + rt3/2.0*sigma_BC
    C2y = -0.5*sigma_BC
    C3x = 0.5*sigma_BB
    C3y = rt3/2.0*sigma_BB + sigma_BC

#   S1x = -rt3/2.0*sigma_BC
#   S1y = -0.5*sigma_BC - sigma_CS
#   S2x = sigma_BB + rt3/2.0*sigma_BC + rt3/2.0*sigma_CS
#   S2y = -0.5*sigma_BC + 0.5*sigma_CS
#   S3x = 0.5*sigma_BB - rt3/2.0*sigma_CS
#   S3y = rt3/2.0*sigma_BB + sigma_BC + 0.5*sigma_CS

    S1Bx = -rt3/2.0*sigma_BC - sigma_CS*math.sin(theta)
    S1By = -0.5*sigma_BC - sigma_CS*math.cos(theta)
    S1Ax = -rt3/2.0*sigma_BC + sigma_CS*math.sin(theta)
    S1Ay = -0.5*sigma_BC - sigma_CS*math.cos(theta)

    S2Bx = sigma_BB + rt3/2.0*sigma_BC + sigma_CS*math.cos(pi/6 - theta)
    S2By = -0.5*sigma_BC + sigma_CS*math.sin(pi/6 - theta)
    S2Ax = sigma_BB + rt3/2.0*sigma_BC + sigma_CS*math.cos(pi/6 + theta)
    S2Ay = -0.5*sigma_BC + sigma_CS*math.sin(pi/6 + theta)

    S3Bx = 0.5*sigma_BB - sigma_CS*math.cos(pi/6 + theta)
    S3By = rt3/2.0*sigma_BB + sigma_BC + sigma_CS*math.sin(pi/6 + theta)
    S3Ax = 0.5*sigma_BB - sigma_CS*math.cos(pi/6 - theta)
    S3Ay = rt3/2.0*sigma_BB + sigma_BC + sigma_CS*math.sin(pi/6 - theta)


    saft_tmc.append(atom(B1x, B1y, 0, 'TMC_B', 0.0, 1))  #1
    saft_tmc.append(atom(B2x, B2y, 0, 'TMC_B', 0.0, 1))  #2
    saft_tmc.append(atom(B3x, B3y, 0, 'TMC_B', 0.0, 1))  #3

    saft_tmc.append(atom(C1x, C1y, 0, 'TMC_C', 0.0, 1))  #4
    saft_tmc.append(atom(C2x, C2y, 0, 'TMC_C', 0.0, 1))  #5
    saft_tmc.append(atom(C3x, C3y, 0, 'TMC_C', 0.0, 1))  #6

    saft_tmc.append(atom(S1Ax, S1Ay, 0, 'TMC_SA', 0.0, 1))  #7
    saft_tmc.append(atom(S1Bx, S1By, 0, 'TMC_SB', 0.0, 1))  #8
    saft_tmc.append(atom(S2Ax, S2Ay, 0, 'TMC_SA', 0.0, 1))  #9
    saft_tmc.append(atom(S2Bx, S2By, 0, 'TMC_SB', 0.0, 1))  #10
    saft_tmc.append(atom(S3Ax, S3Ay, 0, 'TMC_SA', 0.0, 1))  #11
    saft_tmc.append(atom(S3Bx, S3By, 0, 'TMC_SB', 0.0, 1))  #12

    saft_tmc_bonds.append(bond(1, 1, 2))
    saft_tmc_bonds.append(bond(1, 2, 3))
    saft_tmc_bonds.append(bond(1, 1, 3))
        
    saft_tmc_bonds.append(bond(5, 1, 4))
    saft_tmc_bonds.append(bond(5, 2, 5))
    saft_tmc_bonds.append(bond(5, 3, 6))
 
    saft_tmc_bonds.append(bond(6, 4, 7))
    saft_tmc_bonds.append(bond(7, 4, 8))
    saft_tmc_bonds.append(bond(6, 5, 9))
    saft_tmc_bonds.append(bond(7, 5, 10))
    saft_tmc_bonds.append(bond(6, 6, 11))
    saft_tmc_bonds.append(bond(7, 6, 12))

    saft_tmc_angles.append(angle(1, 3, 1, 2))
    saft_tmc_angles.append(angle(1, 2, 3, 1))
    saft_tmc_angles.append(angle(1, 1, 2, 3))
    
    saft_tmc_angles.append(angle(5, 4, 1, 3))
    saft_tmc_angles.append(angle(5, 1, 2, 5))
    saft_tmc_angles.append(angle(5, 6, 3, 1))

    saft_tmc_angles.append(angle(6, 1, 4, 7)) # short angle A
    saft_tmc_angles.append(angle(7, 1, 4, 8)) # wide angle B
    saft_tmc_angles.append(angle(6, 2, 5, 9))
    saft_tmc_angles.append(angle(7, 2, 5, 10))
    saft_tmc_angles.append(angle(6, 3, 6, 11))
    saft_tmc_angles.append(angle(7, 3, 6, 12))

    saft_tmc_dihedrals.append(dihedral(1, 3, 1, 4, 7))
    saft_tmc_dihedrals.append(dihedral(1, 3, 1, 4, 8))
    saft_tmc_dihedrals.append(dihedral(1, 1, 2, 5, 9))
    saft_tmc_dihedrals.append(dihedral(1, 1, 2, 5, 10))
    saft_tmc_dihedrals.append(dihedral(1, 2, 3, 6, 11))
    saft_tmc_dihedrals.append(dihedral(1, 2, 3, 6, 12))

    saft_tmc_impropers.append(improper(1, 4, 2, 3, 1))
    saft_tmc_impropers.append(improper(1, 5, 1, 3, 2))
    saft_tmc_impropers.append(improper(1, 6, 1, 2, 3))

    bond_types.extend(bonds)
    angle_types.extend(angles)
    dihedral_types.extend(dihedrals)
    improper_types.extend(impropers)

    return saft_tmc, saft_tmc_bonds, saft_tmc_angles, saft_tmc_dihedrals, \
           saft_tmc_impropers, type_list, bond_types, angle_types, \
           dihedral_types, improper_types

def SAFT_MPD_molecule_4(type_list, bond_types, angle_types, dihedral_types):
    """
    Returns list of atoms, bonds, angles and dihedrals corresponding 
    to meta-pheyline diamine (MPD) molecule based on
    SAFT 3-segment benzene molecule, from 
    T. Lafitte et al. Molecular Physics, 
    Vol. 110, Nos. 11-12, 10-20 June 2012, 1189-1203 

    ** This version also includes dihedral interaction and correct
       bond distances from DFT structure 

    ** Version 4 contains two sticky sites labelled such that 
       correct adhesion to TMC is formed 
    
    1:MPD_B, 2:MPD_B, 3:MPD_B, 4:MPD_N, 5:MPD_N, 6:MPD_SA, 7:MPD_SB,
    8:MPD_SA, 9:MPD_SB
     
    Defines bond type 1, 2, 3, 4
            angle type 1, 2, 3, 4
            dihedral type 1
    """ 

    pi = 4.0*math.atan(1.0)

    ang = 60*pi/180     # bond angle in radians
    sigma = 2.42           # bond length in Angstrom

    d = 0.56 # Offset parameter for sigma_SS = 1 Angstrom
    theta = 30*pi/180 # angle of sticky sites wrt bond direction


    sigma_BB = 2.42   # bond between SAFT benzene beads in Angstrom
    sigma_BN = 1.42   # bond between SAFT benzene beads and carbon in Angstrom
    sigma_NS = 0.69   # bond between carbon and sticky site in Angstrom

    #sigma_CS -= 0.784

    saft_mpd = []
    saft_mpd_bonds = []
    saft_mpd_angles = []
    saft_mpd_dihedrals = []

    saft_mpd_types = ['MPD_B', 'MPD_N', 'MPD_SA', 'MPD_SB']

    angles = [1, 2, 3, 4]
    bonds = [1, 2, 3, 4]
    dihedrals = [1]
            
    type_list.extend(saft_mpd_types)

    rt3 = math.sqrt(3.0)

    B1x = 0.0
    B1y = 0.0
    B2x = sigma_BB
    B2y = 0.0
    B3x = 0.5*sigma_BB
    B3y = rt3/2.0*sigma_BB

    N1x = -rt3/2.0*sigma_BN
    N1y = -0.5*sigma_BN
    N2x = sigma_BB + rt3/2.0*sigma_BN
    N2y = -0.5*sigma_BN

#   S1x = -rt3/2.0*sigma_BN
#   S1y = -0.5*sigma_BN - sigma_NS
#   S2x = sigma_BB + rt3/2.0*sigma_BN + rt3/2.0*sigma_NS
#   S2y = -0.5*sigma_BN + 0.5*sigma_NS

    S1Ax = -rt3/2.0*sigma_BN - sigma_NS*math.sin(theta)
    S1Ay = -0.5*sigma_BN - sigma_NS*math.cos(theta)
    S1Bx = -rt3/2.0*sigma_BN + sigma_NS*math.sin(theta)
    S1By = -0.5*sigma_BN - sigma_NS*math.cos(theta)

    S2Ax = sigma_BB + rt3/2.0*sigma_BN + sigma_NS*math.cos(pi/6 - theta)
    S2Ay = -0.5*sigma_BN + sigma_NS*math.sin(pi/6 - theta)
    S2Bx = sigma_BB + rt3/2.0*sigma_BN + sigma_NS*math.cos(pi/6 + theta)
    S2By = -0.5*sigma_BN + sigma_NS*math.sin(pi/6 + theta)


    saft_mpd.append(atom(B1x, B1y, 0, 'MPD_B', 0.0, 1))  #1
    saft_mpd.append(atom(B2x, B2y, 0, 'MPD_B', 0.0, 1))  #2
    saft_mpd.append(atom(B3x, B3y, 0, 'MPD_B', 0.0, 1))  #3

    saft_mpd.append(atom(N1x, N1y, 0, 'MPD_N', 0.0, 1))  #4
    saft_mpd.append(atom(N2x, N2y, 0, 'MPD_N', 0.0, 1))  #5

    saft_mpd.append(atom(S1Ax, S1Ay, 0, 'MPD_SA', 0.0, 1))  #6
    saft_mpd.append(atom(S1Bx, S1By, 0, 'MPD_SB', 0.0, 1))  #7
    saft_mpd.append(atom(S2Ax, S2Ay, 0, 'MPD_SA', 0.0, 1))  #8
    saft_mpd.append(atom(S2Bx, S2By, 0, 'MPD_SB', 0.0, 1))  #9

    saft_mpd_bonds.append(bond(1, 1, 2))
    saft_mpd_bonds.append(bond(1, 2, 3))
    saft_mpd_bonds.append(bond(1, 1, 3))
        
    saft_mpd_bonds.append(bond(2, 1, 4))
    saft_mpd_bonds.append(bond(2, 2, 5))
 
    saft_mpd_bonds.append(bond(3, 4, 6))
    saft_mpd_bonds.append(bond(4, 4, 7))

    saft_mpd_bonds.append(bond(3, 5, 8))
    saft_mpd_bonds.append(bond(4, 5, 9))

    saft_mpd_angles.append(angle(1, 3, 1, 2))
    saft_mpd_angles.append(angle(1, 2, 3, 1))
    saft_mpd_angles.append(angle(1, 1, 2, 3))
    
    saft_mpd_angles.append(angle(2, 3, 1, 4))
    saft_mpd_angles.append(angle(2, 1, 2, 5))

    saft_mpd_angles.append(angle(3, 1, 4, 6)) # short angle
    saft_mpd_angles.append(angle(4, 1, 4, 7)) # wide angle
    saft_mpd_angles.append(angle(3, 2, 5, 8))
    saft_mpd_angles.append(angle(4, 2, 5, 9))

    saft_mpd_dihedrals.append(dihedral(1, 3, 1, 4, 6))
    saft_mpd_dihedrals.append(dihedral(1, 3, 1, 4, 7))
    saft_mpd_dihedrals.append(dihedral(1, 1, 2, 5, 8))
    saft_mpd_dihedrals.append(dihedral(1, 1, 2, 5, 9))

    bond_types.extend(bonds)
    angle_types.extend(angles)
    dihedral_types.extend(dihedrals)

    return saft_mpd, saft_mpd_bonds, saft_mpd_angles, saft_mpd_dihedrals, \
           type_list, bond_types, angle_types, dihedral_types

def SAFT_TMC_molecule_4(type_list, bond_types, angle_types, dihedral_types):
    """
    Returns list of atoms, bonds, angles and dihgedrals corresponding 
    to Benzene-1,3,5-tricarboxylicphenylamide (TMC) molecule based on
    SAFT 3-segment benzene molecule, from 
    T. Lafitte et al. Molecular Physics, 
    Vol. 110, Nos. 11-12, 10-20 June 2012, 1189-1203 

    ** This version also includes dihedral interaction and correct
       bond distances from DFT structure 

    ** Version 4 contains two sticky sites labelled such that 
       correct adhesion to MPD is formed 
    
    1:TMC_B, 2:TMC_B, 3:TMC_B, 4:TMC_C, 5:TMC_C, 6:TMC_C,
    7:TMS_SA, 8:TMC_SB, 9:TMC_SA, 10:TMC_SB, 11:TMC_SA, 12:TMC_SB  
    Defines bond type 1, 5, 6, 7
            angle type 1, 5, 6, 7 
            dihedral type 1
    """ 

    pi = 4.0*math.atan(1.0)

    ang = 60*pi/180     # bond angle in radians
    sigma = 2.42           # bond length in Angstrom

    d = 0.56 # Offset parameter for sigma_SS = 1 Angstrom

    sigma_BB = 2.42   # bond between SAFT benzene beads in Angstrom
    sigma_BC = 1.51   # bond between SAFT benzene beads and carbon in Angstrom
    sigma_CS = 0.69   # bond between carbon and sticky site in Angstrom

    theta = 30*pi/180 # angle of sticky sites wrt bond direction

    #sigma_CS -= 0.784

    saft_tmc = []
    saft_tmc_bonds = []
    saft_tmc_angles = []
    saft_tmc_dihedrals = []

    saft_tmc_types = ['TMC_B', 'TMC_C', 'TMC_SA', 'TMC_SB']

    angles = [1, 5, 6, 7]
    bonds = [1, 5, 6, 7]
    dihedrals = [1]
            
    type_list.extend(saft_tmc_types)

    rt3 = math.sqrt(3.0)

    B1x = 0.0
    B1y = 0.0
    B2x = sigma_BB
    B2y = 0.0
    B3x = 0.5*sigma_BB
    B3y = rt3/2.0*sigma_BB

    C1x = -rt3/2.0*sigma_BC
    C1y = -0.5*sigma_BC
    C2x = sigma_BB + rt3/2.0*sigma_BC
    C2y = -0.5*sigma_BC
    C3x = 0.5*sigma_BB
    C3y = rt3/2.0*sigma_BB + sigma_BC

#   S1x = -rt3/2.0*sigma_BC
#   S1y = -0.5*sigma_BC - sigma_CS
#   S2x = sigma_BB + rt3/2.0*sigma_BC + rt3/2.0*sigma_CS
#   S2y = -0.5*sigma_BC + 0.5*sigma_CS
#   S3x = 0.5*sigma_BB - rt3/2.0*sigma_CS
#   S3y = rt3/2.0*sigma_BB + sigma_BC + 0.5*sigma_CS

    S1Bx = -rt3/2.0*sigma_BC - sigma_CS*math.sin(theta)
    S1By = -0.5*sigma_BC - sigma_CS*math.cos(theta)
    S1Ax = -rt3/2.0*sigma_BC + sigma_CS*math.sin(theta)
    S1Ay = -0.5*sigma_BC - sigma_CS*math.cos(theta)

    S2Bx = sigma_BB + rt3/2.0*sigma_BC + sigma_CS*math.cos(pi/6 - theta)
    S2By = -0.5*sigma_BC + sigma_CS*math.sin(pi/6 - theta)
    S2Ax = sigma_BB + rt3/2.0*sigma_BC + sigma_CS*math.cos(pi/6 + theta)
    S2Ay = -0.5*sigma_BC + sigma_CS*math.sin(pi/6 + theta)

    S3Bx = 0.5*sigma_BB - sigma_CS*math.cos(pi/6 + theta)
    S3By = rt3/2.0*sigma_BB + sigma_BC + sigma_CS*math.sin(pi/6 + theta)
    S3Ax = 0.5*sigma_BB - sigma_CS*math.cos(pi/6 - theta)
    S3Ay = rt3/2.0*sigma_BB + sigma_BC + sigma_CS*math.sin(pi/6 - theta)


    saft_tmc.append(atom(B1x, B1y, 0, 'TMC_B', 0.0, 1))  #1
    saft_tmc.append(atom(B2x, B2y, 0, 'TMC_B', 0.0, 1))  #2
    saft_tmc.append(atom(B3x, B3y, 0, 'TMC_B', 0.0, 1))  #3

    saft_tmc.append(atom(C1x, C1y, 0, 'TMC_C', 0.0, 1))  #4
    saft_tmc.append(atom(C2x, C2y, 0, 'TMC_C', 0.0, 1))  #5
    saft_tmc.append(atom(C3x, C3y, 0, 'TMC_C', 0.0, 1))  #6

    saft_tmc.append(atom(S1Ax, S1Ay, 0, 'TMC_SA', 0.0, 1))  #7
    saft_tmc.append(atom(S1Bx, S1By, 0, 'TMC_SB', 0.0, 1))  #8
    saft_tmc.append(atom(S2Ax, S2Ay, 0, 'TMC_SA', 0.0, 1))  #9
    saft_tmc.append(atom(S2Bx, S2By, 0, 'TMC_SB', 0.0, 1))  #10
    saft_tmc.append(atom(S3Ax, S3Ay, 0, 'TMC_SA', 0.0, 1))  #11
    saft_tmc.append(atom(S3Bx, S3By, 0, 'TMC_SB', 0.0, 1))  #12

    saft_tmc_bonds.append(bond(1, 1, 2))
    saft_tmc_bonds.append(bond(1, 2, 3))
    saft_tmc_bonds.append(bond(1, 1, 3))
        
    saft_tmc_bonds.append(bond(5, 1, 4))
    saft_tmc_bonds.append(bond(5, 2, 5))
    saft_tmc_bonds.append(bond(5, 3, 6))
 
    saft_tmc_bonds.append(bond(6, 4, 7))
    saft_tmc_bonds.append(bond(7, 4, 8))
    saft_tmc_bonds.append(bond(6, 5, 9))
    saft_tmc_bonds.append(bond(7, 5, 10))
    saft_tmc_bonds.append(bond(6, 6, 11))
    saft_tmc_bonds.append(bond(7, 6, 12))

    saft_tmc_angles.append(angle(1, 3, 1, 2))
    saft_tmc_angles.append(angle(1, 2, 3, 1))
    saft_tmc_angles.append(angle(1, 1, 2, 3))
    
    saft_tmc_angles.append(angle(5, 4, 1, 3))
    saft_tmc_angles.append(angle(5, 1, 2, 5))
    saft_tmc_angles.append(angle(5, 6, 3, 1))

    saft_tmc_angles.append(angle(6, 1, 4, 7)) # short angle A
    saft_tmc_angles.append(angle(7, 1, 4, 8)) # wide angle B
    saft_tmc_angles.append(angle(6, 2, 5, 9))
    saft_tmc_angles.append(angle(7, 2, 5, 10))
    saft_tmc_angles.append(angle(6, 3, 6, 11))
    saft_tmc_angles.append(angle(7, 3, 6, 12))

    saft_tmc_dihedrals.append(dihedral(1, 3, 1, 4, 7))
    saft_tmc_dihedrals.append(dihedral(1, 3, 1, 4, 8))
    saft_tmc_dihedrals.append(dihedral(1, 1, 2, 5, 9))
    saft_tmc_dihedrals.append(dihedral(1, 1, 2, 5, 10))
    saft_tmc_dihedrals.append(dihedral(1, 2, 3, 6, 11))
    saft_tmc_dihedrals.append(dihedral(1, 2, 3, 6, 12))

    bond_types.extend(bonds)
    angle_types.extend(angles)
    dihedral_types.extend(dihedrals)

    return saft_tmc, saft_tmc_bonds, saft_tmc_angles, saft_tmc_dihedrals,\
           type_list, bond_types, angle_types, dihedral_types

def SAFT_MPD_molecule_3(type_list, bond_types, angle_types, dihedral_types):
    """
    Returns list of atoms, bonds, angles and dihedrals corresponding 
    to meta-pheyline diamine (MPD) molecule based on
    SAFT 3-segment benzene molecule, from 
    T. Lafitte et al. Molecular Physics, 
    Vol. 110, Nos. 11-12, 10-20 June 2012, 1189-1203 

    ** This version also includes dihedral interaction and correct
       bond distances from DFT structure 
    
    1:MPD_B, 2:MPD_B, 3:MPD_B, 4:MPD_N, 5:MPD_N, 6:MPD_S, 7:MPD_S
    Defines bond type 1, 2
            angle type 1, 2
            dihedral type 2
    """ 

    pi = 4.0*math.atan(1.0)

    ang = 60*pi/180     # bond angle in radians
    sigma = 2.42           # bond length in Angstrom

    d = 0.56 # Offset parameter for sigma_SS = 1 Angstrom

    sigma_BB = 2.42   # bond between SAFT benzene beads in Angstrom
    sigma_BN = 1.42   # bond between SAFT benzene beads and carbon in Angstrom
    sigma_NS = 0.69 - d  # bond between carbon and sticky site in Angstrom

    #sigma_CS -= 0.784

    saft_mpd = []
    saft_mpd_bonds = []
    saft_mpd_angles = []
    saft_mpd_dihedrals = []

    saft_mpd_types = ['MPD_B', 'MPD_N', 'MPD_S']

    angles = [1, 4, 5]
    bonds = [1, 4, 5]
    dihedrals = [1]
            
    type_list.extend(saft_mpd_types)

    rt3 = math.sqrt(3.0)

    B1x = 0.0
    B1y = 0.0
    B2x = sigma_BB
    B2y = 0.0
    B3x = 0.5*sigma_BB
    B3y = rt3/2.0*sigma_BB

    N1x = -rt3/2.0*sigma_BN
    N1y = -0.5*sigma_BN
    N2x = sigma_BB + rt3/2.0*sigma_BN
    N2y = -0.5*sigma_BN

    S1x = -rt3/2.0*sigma_BN
    S1y = -0.5*sigma_BN - sigma_NS
    S2x = sigma_BB + rt3/2.0*sigma_BN + rt3/2.0*sigma_NS
    S2y = -0.5*sigma_BN + 0.5*sigma_NS

    saft_mpd.append(atom(B1x, B1y, 0, 'MPD_B', 0.0, 1))  #1
    saft_mpd.append(atom(B2x, B2y, 0, 'MPD_B', 0.0, 1))  #2
    saft_mpd.append(atom(B3x, B3y, 0, 'MPD_B', 0.0, 1))  #3

    saft_mpd.append(atom(N1x, N1y, 0, 'MPD_N', 0.0, 1))  #4
    saft_mpd.append(atom(N2x, N2y, 0, 'MPD_N', 0.0, 1))  #5

    saft_mpd.append(atom(S1x, S1y, 0, 'MPD_S', 0.0, 1))  #7
    saft_mpd.append(atom(S2x, S2y, 0, 'MPD_S', 0.0, 1))  #8

    saft_mpd_bonds.append(bond(1, 1, 2))
    saft_mpd_bonds.append(bond(1, 2, 3))
    saft_mpd_bonds.append(bond(1, 1, 3))
        
    saft_mpd_bonds.append(bond(4, 1, 4))
    saft_mpd_bonds.append(bond(4, 2, 5))
 
    saft_mpd_bonds.append(bond(5, 4, 6))
    saft_mpd_bonds.append(bond(5, 5, 7))

    saft_mpd_angles.append(angle(1, 3, 1, 2))
    saft_mpd_angles.append(angle(1, 2, 3, 1))
    saft_mpd_angles.append(angle(1, 1, 2, 3))
    
    saft_mpd_angles.append(angle(4, 3, 1, 4))
    saft_mpd_angles.append(angle(4, 1, 2, 5))

    saft_mpd_angles.append(angle(5, 1, 4, 6))
    saft_mpd_angles.append(angle(5, 2, 5, 7))

    saft_mpd_dihedrals.append(dihedral(1, 3, 1, 4, 6))
    saft_mpd_dihedrals.append(dihedral(1, 1, 2, 5, 7))

    bond_types.extend(bonds)
    angle_types.extend(angles)
    dihedral_types.extend(dihedrals)

    return saft_mpd, saft_mpd_bonds, saft_mpd_angles, saft_mpd_dihedrals, \
           type_list, bond_types, angle_types, dihedral_types

def SAFT_TMC_molecule_3(type_list, bond_types, angle_types, dihedral_types):
    """
    Returns list of atoms, bonds, angles and dihgedrals corresponding 
    to Benzene-1,3,5-tricarboxylicphenylamide (TMC) molecule based on
    SAFT 3-segment benzene molecule, from 
    T. Lafitte et al. Molecular Physics, 
    Vol. 110, Nos. 11-12, 10-20 June 2012, 1189-1203 

    ** This version also includes dihedral interaction and correct
       bond distances from DFT structure 
    
    1:TMC_B, 2:TMC_B, 3:TMC_B, 4:TMC_C, 5:TMC_C, 6:TMC_C,
    7:TMS_S, 8:TMC_S, 9:TMC_S 
    Defines bond type 1,2
            angle type 1, 3, 4
            dihedral type 1
    """ 

    pi = 4.0*math.atan(1.0)

    ang = 60*pi/180     # bond angle in radians
    sigma = 2.42           # bond length in Angstrom

    d = 0.56 # Offset parameter for sigma_SS = 1 Angstrom

    sigma_BB = 2.42   # bond between SAFT benzene beads in Angstrom
    sigma_BC = 1.51   # bond between SAFT benzene beads and carbon in Angstrom
    sigma_CS = 0.69 - d  # bond between carbon and sticky site in Angstrom

    #sigma_CS -= 0.784

    saft_tmc = []
    saft_tmc_bonds = []
    saft_tmc_angles = []
    saft_tmc_dihedrals = []

    saft_tmc_types = ['TMC_B', 'TMC_C', 'TMC_S']

    angles = [1, 2, 3]
    bonds = [1, 2, 3]
    dihedrals = [1]
            
    type_list.extend(saft_tmc_types)

    rt3 = math.sqrt(3.0)

    B1x = 0.0
    B1y = 0.0
    B2x = sigma_BB
    B2y = 0.0
    B3x = 0.5*sigma_BB
    B3y = rt3/2.0*sigma_BB

    C1x = -rt3/2.0*sigma_BC
    C1y = -0.5*sigma_BC
    C2x = sigma_BB + rt3/2.0*sigma_BC
    C2y = -0.5*sigma_BC
    C3x = 0.5*sigma_BB
    C3y = rt3/2.0*sigma_BB + sigma_BC

    S1x = -rt3/2.0*sigma_BC
    S1y = -0.5*sigma_BC - sigma_CS
    S2x = sigma_BB + rt3/2.0*sigma_BC + rt3/2.0*sigma_CS
    S2y = -0.5*sigma_BC + 0.5*sigma_CS
    S3x = 0.5*sigma_BB - rt3/2.0*sigma_CS
    S3y = rt3/2.0*sigma_BB + sigma_BC + 0.5*sigma_CS

    saft_tmc.append(atom(B1x, B1y, 0, 'TMC_B', 0.0, 1))  #1
    saft_tmc.append(atom(B2x, B2y, 0, 'TMC_B', 0.0, 1))  #2
    saft_tmc.append(atom(B3x, B3y, 0, 'TMC_B', 0.0, 1))  #3

    saft_tmc.append(atom(C1x, C1y, 0, 'TMC_C', 0.0, 1))  #4
    saft_tmc.append(atom(C2x, C2y, 0, 'TMC_C', 0.0, 1))  #5
    saft_tmc.append(atom(C3x, C3y, 0, 'TMC_C', 0.0, 1))  #6

    saft_tmc.append(atom(S1x, S1y, 0, 'TMC_S', 0.0, 1))  #7
    saft_tmc.append(atom(S2x, S2y, 0, 'TMC_S', 0.0, 1))  #8
    saft_tmc.append(atom(S3x, S3y, 0, 'TMC_S', 0.0, 1))  #9

    saft_tmc_bonds.append(bond(1, 1, 2))
    saft_tmc_bonds.append(bond(1, 2, 3))
    saft_tmc_bonds.append(bond(1, 1, 3))
        
    saft_tmc_bonds.append(bond(2, 1, 4))
    saft_tmc_bonds.append(bond(2, 2, 5))
    saft_tmc_bonds.append(bond(2, 3, 6))
 
    saft_tmc_bonds.append(bond(3, 4, 7))
    saft_tmc_bonds.append(bond(3, 5, 8))
    saft_tmc_bonds.append(bond(3, 6, 9))

    saft_tmc_angles.append(angle(1, 3, 1, 2))
    saft_tmc_angles.append(angle(1, 2, 3, 1))
    saft_tmc_angles.append(angle(1, 1, 2, 3))
    
    saft_tmc_angles.append(angle(2, 4, 1, 3))
    saft_tmc_angles.append(angle(2, 1, 2, 5))
    saft_tmc_angles.append(angle(2, 6, 3, 1))

    saft_tmc_angles.append(angle(3, 1, 4, 7))
    saft_tmc_angles.append(angle(3, 2, 5, 8))
    saft_tmc_angles.append(angle(3, 3, 6, 9))

    saft_tmc_dihedrals.append(dihedral(1, 3, 1, 4, 7))
    saft_tmc_dihedrals.append(dihedral(1, 1, 2, 5, 8))
    saft_tmc_dihedrals.append(dihedral(1, 2, 3, 6, 9))

    bond_types.extend(bonds)
    angle_types.extend(angles)
    dihedral_types.extend(dihedrals)

    return saft_tmc, saft_tmc_bonds, saft_tmc_angles, saft_tmc_dihedrals,\
           type_list, bond_types, angle_types, dihedral_types

def SAFT_TMC_molecule_2(type_list, bond_types, angle_types):
    """
    Returns list of atoms, bonds and angles corresponding 
    to Benzene-1,3,5-tricarboxylicphenylamide molecule based on
    SAFT 3-segment benzene molecule, from 
    T. Lafitte et al. Molecular Physics, 
    Vol. 110, Nos. 11-12, 10-20 June 2012, 1189-1203 
    1:SAFT_B,2:SAFT_B,3:SAFT_B,4:SAFT_B,5:SAFT_B,6:SAFT_B,
    7:SAFT_O,8:SAFT_O,9:SAFT_O 
    Defines bond type 1,2
            angle type 1, 3, 4
    """ 

    pi = 4.0*math.atan(1.0)

    ang = 60*pi/180     # bond angle in radians
    sigma = 2.42           # bond length in Angstrom

    sigma_BB = 2.42   # bond between SAFT benzene beads in Angstrom
    sigma_BC = 1.51   # bond between SAFT benzene beads and carbon in Angstrom
    sigma_CS = 1.37   # bond between carbon and sticky site in Angstrom

    sigma_CS -= 0.784

    saft_tmc = []
    saft_tmc_bonds = []
    saft_tmc_angles = []

    saft_tmc_types = ['SAFT_B_TMC', 'SAFT_C_TMC', 'SAFT_S_TMC']

    angles = [1,3,4]
    bonds = [1,3]

              

    type_list.extend(saft_tmc_types)

    rt3 = math.sqrt(3.0)

    B1x = 0.0
    B1y = 0.0
    B2x = sigma_BB
    B2y = 0.0
    B3x = 0.5*sigma_BB
    B3y = rt3/2.0*sigma_BB

    C1x = -rt3/2.0*sigma_BC
    C1y = -0.5*sigma_BC
    C2x = sigma_BB + rt3/2.0*sigma_BC
    C2y = -0.5*sigma_BC
    C3x = 0.5*sigma_BB
    C3y = rt3/2.0*sigma_BB + sigma_BC

    S1x = -rt3/2.0*sigma_BC
    S1y = -0.5*sigma_BC - sigma_CS
    S2x = sigma_BB + rt3/2.0*sigma_BC + rt3/2.0*sigma_CS
    S2y = -0.5*sigma_BC + 0.5*sigma_CS
    S3x = 0.5*sigma_BB - rt3/2.0*sigma_CS
    S3y = rt3/2.0*sigma_BB + sigma_BC + 0.5*sigma_CS


    saft_tmc.append(atom(B1x, B1y, 0, 'SAFT_B_TMC', 0.0, 1))  #1
    saft_tmc.append(atom(B2x, B2y, 0, 'SAFT_B_TMC', 0.0, 1))  #2
    saft_tmc.append(atom(B3x, B3y, 0, 'SAFT_B_TMC', 0.0, 1))  #3

    saft_tmc.append(atom(C1x, C1y, 0, 'SAFT_C_TMC', 0.0, 1))  #4
    saft_tmc.append(atom(C2x, C2y, 0, 'SAFT_C_TMC', 0.0, 1))  #5
    saft_tmc.append(atom(C3x, C3y, 0, 'SAFT_C_TMC', 0.0, 1))  #6

    saft_tmc.append(atom(S1x, S1y, 0, 'SAFT_S_TMC', 0.0, 1))  #7
    saft_tmc.append(atom(S2x, S2y, 0, 'SAFT_S_TMC', 0.0, 1))  #8
    saft_tmc.append(atom(S3x, S3y, 0, 'SAFT_S_TMC', 0.0, 1))  #9

    saft_tmc_bonds.append(bond(1, 1, 2))
    saft_tmc_bonds.append(bond(1, 2, 3))
    saft_tmc_bonds.append(bond(1, 1, 3))
        
    saft_tmc_bonds.append(bond(1, 1, 4))
    saft_tmc_bonds.append(bond(1, 2, 5))
    saft_tmc_bonds.append(bond(1, 3, 6))
 
    saft_tmc_bonds.append(bond(3, 4, 7))
    saft_tmc_bonds.append(bond(3, 5, 8))
    saft_tmc_bonds.append(bond(3, 6, 9))

    saft_tmc_angles.append(angle(1, 3, 1, 2))
    saft_tmc_angles.append(angle(1, 2, 3, 1))
    saft_tmc_angles.append(angle(1, 1, 2, 3))
    
    saft_tmc_angles.append(angle(3, 4, 1, 3))
    saft_tmc_angles.append(angle(3, 1, 2, 5))
    saft_tmc_angles.append(angle(3, 6, 3, 1))

    saft_tmc_angles.append(angle(4, 1, 4, 7))
    saft_tmc_angles.append(angle(4, 2, 5, 8))
    saft_tmc_angles.append(angle(4, 3, 6, 9))

    bond_types.extend(bonds)
    angle_types.extend(angles)

    return saft_tmc, saft_tmc_bonds, saft_tmc_angles, \
           type_list, bond_types, angle_types

def SAFT_TMC_molecule(type_list, bond_types, angle_types):
    """
    Returns list of atoms, bonds and angles corresponding 
    to Benzene-1,3,5-tricarboxylicphenylamide molecule based on
    SAFT 3-segment benzene molecule, from 
    T. Lafitte et al. Molecular Physics, 
    Vol. 110, Nos. 11-12, 10-20 June 2012, 1189-1203 
    1:SAFT_B,2:SAFT_B,3:SAFT_B,4:SAFT_B,5:SAFT_B,6:SAFT_B,
    7:SAFT_O,8:SAFT_O,9:SAFT_O 
    Defines bond type 1,2
            angle type 1, 3, 4
    """ 

    pi = 4.0*math.atan(1.0)

    ang = 60*pi/180     # bond angle in radians
    sigma = 3.490           # bond length in Angstrom

    saft_tmc = []
    saft_tmc_bonds = []
    saft_tmc_angles = []

    saft_tmc_types = ['SAFT_B_TMC', 'SAFT_O']

    angles = [1,3,4]
    bonds = [1,3]

              

    type_list.extend(saft_tmc_types)

    rt3 = math.sqrt(3.0)

    B1x = 0.0
    B1y = 0.0
    B2x = sigma
    B2y = 0.0
    B3x = sigma/2.0
    B3y = sigma*rt3/2.0

    B4x = -sigma*rt3/2.0
    B4y = -sigma/2.0
    B5x = sigma*(1.0+rt3/2.0)
    B5y = -sigma/2.0
    B6x = sigma/2.0
    B6y = sigma*(1.0+rt3/2.0)

    O1x = -0.5*sigma*(1.0+rt3)
    O1y = -sigma/2.0
    O2x = sigma*(1.0+rt3/2.0)
    O2y = -sigma
    O3x = 0.5*sigma*(1.0+rt3/2.0)
    O3y = sigma*(1.25+rt3/2.0) 

    saft_tmc.append(atom(B1x, B1y, 0, 'SAFT_B_TMC', 0.0, 1))  #1
    saft_tmc.append(atom(B2x, B2y, 0, 'SAFT_B_TMC', 0.0, 1))  #2
    saft_tmc.append(atom(B3x, B3y, 0, 'SAFT_B_TMC', 0.0, 1))  #3

    saft_tmc.append(atom(B4x, B4y, 0, 'SAFT_B_TMC', 0.0, 1))  #4
    saft_tmc.append(atom(B5x, B5y, 0, 'SAFT_B_TMC', 0.0, 1))  #5
    saft_tmc.append(atom(B6x, B6y, 0, 'SAFT_B_TMC', 0.0, 1))  #6

    saft_tmc.append(atom(O1x, O1y, 0, 'SAFT_O', 0.0, 1))  #7
    saft_tmc.append(atom(O2x, O2y, 0, 'SAFT_O', 0.0, 1))  #8
    saft_tmc.append(atom(O3x, O3y, 0, 'SAFT_O', 0.0, 1))  #9

    saft_tmc_bonds.append(bond(1, 1, 2))
    saft_tmc_bonds.append(bond(1, 2, 3))
    saft_tmc_bonds.append(bond(1, 1, 3))
        
    saft_tmc_bonds.append(bond(1, 1, 4))
    saft_tmc_bonds.append(bond(1, 2, 5))
    saft_tmc_bonds.append(bond(1, 3, 6))
 
    saft_tmc_bonds.append(bond(3, 4, 7))
    saft_tmc_bonds.append(bond(3, 5, 8))
    saft_tmc_bonds.append(bond(3, 6, 9))

    saft_tmc_angles.append(angle(1, 3, 1, 2))
    saft_tmc_angles.append(angle(1, 2, 3, 1))
    saft_tmc_angles.append(angle(1, 1, 2, 3))
    
    saft_tmc_angles.append(angle(3, 4, 1, 3))
    saft_tmc_angles.append(angle(3, 1, 2, 5))
    saft_tmc_angles.append(angle(3, 6, 3, 1))

    saft_tmc_angles.append(angle(4, 1, 4, 7))
    saft_tmc_angles.append(angle(4, 2, 5, 8))
    saft_tmc_angles.append(angle(4, 3, 6, 9))

    bond_types.extend(bonds)
    angle_types.extend(angles)

    return saft_tmc, saft_tmc_bonds, saft_tmc_angles, \
           type_list, bond_types, angle_types


def SAFT_MPD_molecule_2(type_list, bond_types, angle_types):
    """
    Returns list of atoms, bond and angles corresponding 
    to Meta Phenyline Diamine molecule based on
    SAFT 3-segment benzene molecule, from 
    T. Lafitte et al. Molecular Physics, 
    Vol. 110, Nos. 11-12, 10-20 June 2012, 1189-1203 
    1:SAFT_B,2:SAFT_B,3:SAFT_B,4:SAFT_NH,5:SAFT_NH 
    Defines angle and bond type 1 and 2
    """ 

    pi = 4.0*math.atan(1.0)

    ang = 60*pi/180     # bond angle in radians
    sigma_BB = 2.42           # bond length in Angstrom
    sigma_BN = 1.42 - 0.336

    saft_mpd = []
    saft_mpd_bonds = []
    saft_mpd_angles = []

    saft_mpd_types = ['SAFT_B_MPD', 'SAFT_NH']

    angles = [1,2]
    bonds = [1,2]
    saft_tmc_dihedrals(dihedral(1, 3, 1, 4, 7))
    saft_tmc_dihedrals(dihedral(1, 1, 2, 5, 8))
    saft_tmc_dihedrals(dihedral(1, 2, 3, 6, 9))

    bond_types.extend(bonds)
    angle_types.extend(angles)
    dihedral_types.extend(dihedrals)

    return saft_tmc, saft_tmc_bonds, saft_tmc_angles, \
           type_list, bond_types, angle_types
def SAFT_TMC_molecule_2(type_list, bond_types, angle_types):
    """
    Returns list of atoms, bonds and angles corresponding 
    to Benzene-1,3,5-tricarboxylicphenylamide molecule based on
    SAFT 3-segment benzene molecule, from 
    T. Lafitte et al. Molecular Physics, 
    Vol. 110, Nos. 11-12, 10-20 June 2012, 1189-1203 
    1:SAFT_B,2:SAFT_B,3:SAFT_B,4:SAFT_B,5:SAFT_B,6:SAFT_B,
    7:SAFT_O,8:SAFT_O,9:SAFT_O 
    Defines bond type 1,2
            angle type 1, 3, 4
    """ 

    pi = 4.0*math.atan(1.0)

    ang = 60*pi/180     # bond angle in radians
    sigma = 2.42           # bond length in Angstrom

    sigma_BB = 2.42   # bond between SAFT benzene beads in Angstrom
    sigma_BC = 1.51   # bond between SAFT benzene beads and carbon in Angstrom
    sigma_CS = 1.37   # bond between carbon and sticky site in Angstrom

    sigma_CS -= 0.784

    saft_tmc = []
    saft_tmc_bonds = []
    saft_tmc_angles = []

    saft_tmc_types = ['SAFT_B_TMC', 'SAFT_C_TMC', 'SAFT_S_TMC']

    angles = [1,3,4]
    bonds = [1,3]

              

    type_list.extend(saft_tmc_types)

    rt3 = math.sqrt(3.0)

    B1x = 0.0
    B1y = 0.0
    B2x = sigma_BB
    B2y = 0.0
    B3x = 0.5*sigma_BB
    B3y = rt3/2.0*sigma_BB

    C1x = -rt3/2.0*sigma_BC
    C1y = -0.5*sigma_BC
    C2x = sigma_BB + rt3/2.0*sigma_BC
    C2y = -0.5*sigma_BC
    C3x = 0.5*sigma_BB
    C3y = rt3/2.0*sigma_BB + sigma_BC

    S1x = -rt3/2.0*sigma_BC
    S1y = -0.5*sigma_BC - sigma_CS
    S2x = sigma_BB + rt3/2.0*sigma_BC + rt3/2.0*sigma_CS
    S2y = -0.5*sigma_BC + 0.5*sigma_CS
    S3x = 0.5*sigma_BB - rt3/2.0*sigma_CS
    S3y = rt3/2.0*sigma_BB + sigma_BC + 0.5*sigma_CS


    saft_tmc.append(atom(B1x, B1y, 0, 'SAFT_B_TMC', 0.0, 1))  #1
    saft_tmc.append(atom(B2x, B2y, 0, 'SAFT_B_TMC', 0.0, 1))  #2
    saft_tmc.append(atom(B3x, B3y, 0, 'SAFT_B_TMC', 0.0, 1))  #3

    saft_tmc.append(atom(C1x, C1y, 0, 'SAFT_C_TMC', 0.0, 1))  #4
    saft_tmc.append(atom(C2x, C2y, 0, 'SAFT_C_TMC', 0.0, 1))  #5
    saft_tmc.append(atom(C3x, C3y, 0, 'SAFT_C_TMC', 0.0, 1))  #6

    saft_tmc.append(atom(S1x, S1y, 0, 'SAFT_S_TMC', 0.0, 1))  #7
    saft_tmc.append(atom(S2x, S2y, 0, 'SAFT_S_TMC', 0.0, 1))  #8
    saft_tmc.append(atom(S3x, S3y, 0, 'SAFT_S_TMC', 0.0, 1))  #9

    saft_tmc_bonds.append(bond(1, 1, 2))
    saft_tmc_bonds.append(bond(1, 2, 3))
    saft_tmc_bonds.append(bond(1, 1, 3))
        
    saft_tmc_bonds.append(bond(1, 1, 4))
    saft_tmc_bonds.append(bond(1, 2, 5))
    saft_tmc_bonds.append(bond(1, 3, 6))
 
    saft_tmc_bonds.append(bond(3, 4, 7))
    saft_tmc_bonds.append(bond(3, 5, 8))
    saft_tmc_bonds.append(bond(3, 6, 9))

    saft_tmc_angles.append(angle(1, 3, 1, 2))
    saft_tmc_angles.append(angle(1, 2, 3, 1))
    saft_tmc_angles.append(angle(1, 1, 2, 3))
    
    saft_tmc_angles.append(angle(3, 4, 1, 3))
    saft_tmc_angles.append(angle(3, 1, 2, 5))
    saft_tmc_angles.append(angle(3, 6, 3, 1))

    saft_tmc_angles.append(angle(4, 1, 4, 7))
    saft_tmc_angles.append(angle(4, 2, 5, 8))
    saft_tmc_angles.append(angle(4, 3, 6, 9))

    bond_types.extend(bonds)
    angle_types.extend(angles)

    return saft_tmc, saft_tmc_bonds, saft_tmc_angles, \
           type_list, bond_types, angle_types

def SAFT_TMC_molecule(type_list, bond_types, angle_types):
    """
    Returns list of atoms, bonds and angles corresponding 
    to Benzene-1,3,5-tricarboxylicphenylamide molecule based on
    SAFT 3-segment benzene molecule, from 
    T. Lafitte et al. Molecular Physics, 
    Vol. 110, Nos. 11-12, 10-20 June 2012, 1189-1203 
    1:SAFT_B,2:SAFT_B,3:SAFT_B,4:SAFT_B,5:SAFT_B,6:SAFT_B,
    7:SAFT_O,8:SAFT_O,9:SAFT_O 
    Defines bond type 1,2
            angle type 1, 3, 4
    """ 

    pi = 4.0*math.atan(1.0)

    ang = 60*pi/180     # bond angle in radians
    sigma = 3.490           # bond length in Angstrom

    saft_tmc = []
    saft_tmc_bonds = []
    saft_tmc_angles = []

    saft_tmc_types = ['SAFT_B_TMC', 'SAFT_O']

    angles = [1,3,4]
    bonds = [1,3]

              

    type_list.extend(saft_tmc_types)

    rt3 = math.sqrt(3.0)

    B1x = 0.0
    B1y = 0.0
    B2x = sigma
    B2y = 0.0
    B3x = sigma/2.0
    B3y = sigma*rt3/2.0

    B4x = -sigma*rt3/2.0
    B4y = -sigma/2.0
    B5x = sigma*(1.0+rt3/2.0)
    B5y = -sigma/2.0
    B6x = sigma/2.0
    B6y = sigma*(1.0+rt3/2.0)

    O1x = -0.5*sigma*(1.0+rt3)
    O1y = -sigma/2.0
    O2x = sigma*(1.0+rt3/2.0)
    O2y = -sigma
    O3x = 0.5*sigma*(1.0+rt3/2.0)
    O3y = sigma*(1.25+rt3/2.0) 

    saft_tmc.append(atom(B1x, B1y, 0, 'SAFT_B_TMC', 0.0, 1))  #1
    saft_tmc.append(atom(B2x, B2y, 0, 'SAFT_B_TMC', 0.0, 1))  #2
    saft_tmc.append(atom(B3x, B3y, 0, 'SAFT_B_TMC', 0.0, 1))  #3

    saft_tmc.append(atom(B4x, B4y, 0, 'SAFT_B_TMC', 0.0, 1))  #4
    saft_tmc.append(atom(B5x, B5y, 0, 'SAFT_B_TMC', 0.0, 1))  #5
    saft_tmc.append(atom(B6x, B6y, 0, 'SAFT_B_TMC', 0.0, 1))  #6

    saft_tmc.append(atom(O1x, O1y, 0, 'SAFT_O', 0.0, 1))  #7
    saft_tmc.append(atom(O2x, O2y, 0, 'SAFT_O', 0.0, 1))  #8
    saft_tmc.append(atom(O3x, O3y, 0, 'SAFT_O', 0.0, 1))  #9

    saft_tmc_bonds.append(bond(1, 1, 2))
    saft_tmc_bonds.append(bond(1, 2, 3))
    saft_tmc_bonds.append(bond(1, 1, 3))
        
    saft_tmc_bonds.append(bond(1, 1, 4))
    saft_tmc_bonds.append(bond(1, 2, 5))
    saft_tmc_bonds.append(bond(1, 3, 6))
 
    saft_tmc_bonds.append(bond(3, 4, 7))
    saft_tmc_bonds.append(bond(3, 5, 8))
    saft_tmc_bonds.append(bond(3, 6, 9))

    saft_tmc_angles.append(angle(1, 3, 1, 2))
    saft_tmc_angles.append(angle(1, 2, 3, 1))
    saft_tmc_angles.append(angle(1, 1, 2, 3))
    
    saft_tmc_angles.append(angle(3, 4, 1, 3))
    saft_tmc_angles.append(angle(3, 1, 2, 5))
    saft_tmc_angles.append(angle(3, 6, 3, 1))

    saft_tmc_angles.append(angle(4, 1, 4, 7))
    saft_tmc_angles.append(angle(4, 2, 5, 8))
    saft_tmc_angles.append(angle(4, 3, 6, 9))

    bond_types.extend(bonds)
    angle_types.extend(angles)

    return saft_tmc, saft_tmc_bonds, saft_tmc_angles, \
           type_list, bond_types, angle_types


def SAFT_MPD_molecule_2(type_list, bond_types, angle_types):
    """
    Returns list of atoms, bond and angles corresponding 
    to Meta Phenyline Diamine molecule based on
    SAFT 3-segment benzene molecule, from 
    T. Lafitte et al. Molecular Physics, 
    Vol. 110, Nos. 11-12, 10-20 June 2012, 1189-1203 
    1:SAFT_B,2:SAFT_B,3:SAFT_B,4:SAFT_NH,5:SAFT_NH 
    Defines angle and bond type 1 and 2
    """ 

    pi = 4.0*math.atan(1.0)

    ang = 60*pi/180     # bond angle in radians
    sigma_BB = 2.42           # bond length in Angstrom
    sigma_BN = 1.42 - 0.336

    saft_mpd = []
    saft_mpd_bonds = []
    saft_mpd_angles = []

    saft_mpd_types = ['SAFT_B_MPD', 'SAFT_NH']

    angles = [1,2]
    bonds = [1,2]

    type_list.extend(saft_mpd_types)

    rt3 = math.sqrt(3.0)

    B1x = 0.0
    B1y = 0.0
    B2x = sigma_BB
    B2y = 0.0
    B3x = 0.5*sigma_BB
    B3y = rt3/2.0*sigma_BB

    NH1x = -rt3/2.0*sigma_BN
    NH1y = -0.5*sigma_BN

    NH2x = sigma_BB + rt3/2.0*sigma_BN
    NH2y = -0.5*sigma_BN

    saft_mpd.append(atom(B1x, B1y, 0, 'SAFT_B_MPD', 0.0, 1))    #1
    saft_mpd.append(atom(B2x, B2y, 0, 'SAFT_B_MPD', 0.0, 1))    #2
    saft_mpd.append(atom(B3x, B3y, 0, 'SAFT_B_MPD', 0.0, 1))    #3
    saft_mpd.append(atom(NH1x, NH1y, 0, 'SAFT_NH', 0.0, 1)) #4 
    saft_mpd.append(atom(NH2x, NH2y, 0, 'SAFT_NH', 0.0, 1)) #5 

    saft_mpd_bonds.append(bond(1, 1, 2))
    saft_mpd_bonds.append(bond(1, 2, 3))
    saft_mpd_bonds.append(bond(1, 1, 3))

    saft_mpd_bonds.append(bond(2, 1, 4))
    saft_mpd_bonds.append(bond(2, 2, 5))
    
    saft_mpd_angles.append(angle(1, 3, 1, 2))
    saft_mpd_angles.append(angle(1, 2, 3, 1))
    saft_mpd_angles.append(angle(1, 1, 2, 3))
    
    saft_mpd_angles.append(angle(2, 3, 1, 4))
    saft_mpd_angles.append(angle(2, 3, 2, 5))
 

    bond_types.extend(bonds)
    angle_types.extend(angles)

    return saft_mpd, saft_mpd_bonds, saft_mpd_angles, \
           type_list, bond_types, angle_types

def SAFT_MPD_molecule(type_list, bond_types, angle_types):
    """
    Returns list of atoms, bond and angles corresponding 
    to Meta Phenyline Diamine molecule based on
    SAFT 3-segment benzene molecule, from 
    T. Lafitte et al. Molecular Physics, 
    Vol. 110, Nos. 11-12, 10-20 June 2012, 1189-1203 
    1:SAFT_B,2:SAFT_B,3:SAFT_B,4:SAFT_NH,5:SAFT_NH 
    Defines angle and bond type 1 and 2
    """ 

    pi = 4.0*math.atan(1.0)

    ang = 60*pi/180     # bond angle in radians
    sigma = 3.490           # bond length in Angstrom
    sigma2 = sigma - 0.2


    saft_mpd = []
    saft_mpd_bonds = []
    saft_mpd_angles = []

    saft_mpd_types = ['SAFT_B_MPD', 'SAFT_NH']

    angles = [1,2]
    bonds = [1,2]

    type_list.extend(saft_mpd_types)

    B1x = 0.0
    B1y = 0.0
    B2x = sigma
    B2y = 0.0
    B3x = sigma/2.0
    B3y = sigma*math.sqrt(3.0)/2.0

    NH1x = -sigma2*math.sqrt(3.0)/4.0
    NH1y = -sigma2/4.0

    NH2x = sigma2*(1+(math.sqrt(3.0)/4.0))
    NH2y = -sigma2/4.0

    saft_mpd.append(atom(B1x, B1y, 0, 'SAFT_B_MPD', 0.0, 1))    #1
    saft_mpd.append(atom(B2x, B2y, 0, 'SAFT_B_MPD', 0.0, 1))    #2
    saft_mpd.append(atom(B3x, B3y, 0, 'SAFT_B_MPD', 0.0, 1))    #3
    saft_mpd.append(atom(NH1x, NH1y, 0, 'SAFT_NH', 0.0, 1)) #4 
    saft_mpd.append(atom(NH2x, NH2y, 0, 'SAFT_NH', 0.0, 1)) #5 

    saft_mpd_bonds.append(bond(1, 1, 2))
    saft_mpd_bonds.append(bond(1, 2, 3))
    saft_mpd_bonds.append(bond(1, 1, 3))

    saft_mpd_bonds.append(bond(2, 1, 4))
    saft_mpd_bonds.append(bond(2, 2, 5))
    
    saft_mpd_angles.append(angle(1, 3, 1, 2))
    saft_mpd_angles.append(angle(1, 2, 3, 1))
    saft_mpd_angles.append(angle(1, 1, 2, 3))
    
    saft_mpd_angles.append(angle(2, 3, 1, 4))
    saft_mpd_angles.append(angle(2, 3, 2, 5))
 

    bond_types.extend(bonds)
    angle_types.extend(angles)

    return saft_mpd, saft_mpd_bonds, saft_mpd_angles, \
           type_list, bond_types, angle_types

def SAFTbenzene_molecule_rigid(type_list, bond_types, angle_types):
    """
    Returns list of atoms, bond and angles corresponding 
    to SAFT 3-segment benzene molecule, from 
    T. Lafitte et al. Molecular Physics, 
    Vol. 110, Nos. 11-12, 10-20 June 2012, 1189-1203 
    1:SAFT_B,2:SAFT_B,3:SAFT_B
    Defines angle and bond type 3
    """ 

    pi = 4.0*math.atan(1.0)

    ang = 60*pi/180     # bond angle in radians
    sigma = 3.490           # bond length in Angstrom

    saft_benzene = []
    saft_benzene_bonds = []
    saft_benzene_angles = []

    saft_benzene_types = ['SAFT_B']

    angles = [1]
    bonds = [1]

    type_list.extend(saft_benzene_types)

    B1x = 0.0
    B1y = 0.0
    B2x = sigma
    B2y = 0.0
    B3x = sigma/2.0
    B3y = sigma*math.sqrt(3.0)/2.0

    saft_benzene.append(atom(B1x, B1y, 0, 'SAFT_B', 0.0, 1))  
    saft_benzene.append(atom(B2x, B2y, 0, 'SAFT_B', 0.0, 1))  
    saft_benzene.append(atom(B3x, B3y, 0, 'SAFT_B', 0.0, 1))  

    saft_benzene_bonds.append(bond(1, 1, 2))
    saft_benzene_bonds.append(bond(1, 2, 3))
    saft_benzene_bonds.append(bond(1, 1, 3))
    saft_benzene_angles.append(angle(1, 3, 1, 2))
    saft_benzene_angles.append(angle(1, 2, 3, 1))
    saft_benzene_angles.append(angle(1, 1, 2, 3))
    
    bond_types.extend(bonds)
    angle_types.extend(angles)

    return saft_benzene, saft_benzene_bonds, saft_benzene_angles, \
           type_list, bond_types, angle_types

def SAFTbenzene_molecule_SHAKE(type_list, bond_types, angle_types):
    """
    Returns list of atoms, bond and angles corresponding 
    to SAFT 3-segment benzene molecule, from 
    T. Lafitte et al. Molecular Physics, 
    Vol. 110, Nos. 11-12, 10-20 June 2012, 1189-1203 
    1:SAFT_B,2:SAFT_B,3:SAFT_B
    Defines angle and bond type 3
    """ 

    pi = 4.0*math.atan(1.0)

    ang = 60*pi/180     # bond angle in radians
    sigma = 3.490           # bond length in Angstrom

    saft_benzene = []
    saft_benzene_bonds = []
    saft_benzene_angles = []

    saft_benzene_types = ['SAFT_B']

    angles = [1]
    bonds = [1]

    type_list.extend(saft_benzene_types)

    B1x = 0.0
    B1y = 0.0
    B2x = sigma
    B2y = 0.0
    B3x = sigma/2.0
    B3y = sigma*math.sqrt(3.0)/2.0

    saft_benzene.append(atom(B1x, B1y, 0, 'SAFT_B', 0.0, 1))  
    saft_benzene.append(atom(B2x, B2y, 0, 'SAFT_B', 0.0, 1))  
    saft_benzene.append(atom(B3x, B3y, 0, 'SAFT_B', 0.0, 1))  

    saft_benzene_bonds.append(bond(1, 1, 2))
#   saft_benzene_bonds.append(bond(1, 2, 3))
    saft_benzene_bonds.append(bond(1, 1, 3))
    saft_benzene_angles.append(angle(1, 3, 1, 2))
 #  saft_benzene_angles.append(angle(1, 2, 3, 1))
 #  saft_benzene_angles.append(angle(1, 1, 2, 3))
    
    bond_types.extend(bonds)
    angle_types.extend(angles)

    return saft_benzene, saft_benzene_bonds, saft_benzene_angles, \
           type_list, bond_types, angle_types


def duplicate_hex_cell(hex_cell, atom_list, nx, ny, x_off, y_off, z_off, \
                               a, a1, a2, mol_counter, type_list, hole):
    "Duplicates hexagonal unit cell"
    print "x_off = %f" % x_off 
    print "y_off = %f" % y_off
    print "z_off = %f" % z_off
    small = 1.0e-9
    new_list = []
    for i in range(0,nx):  # in x      
        for j in range(0,ny): # in y
            unit_count = 0
            for c  in hex_cell:

                if hole == 'y' and j == 1 and i ==1: # Put hole in sheet
                    pass
                else:
                    if j % 2 == 0:   # Define offset
                        b = 0 
                    else:
                        b = 3*a/2
                    x_set = c.x + (i*a1) + b + x_off
                    y_set = c.y + (j*a2) + y_off
                    z_set = c.z + z_off
                    ## Apply Periodic Boundaries

                    if x_set - small > nx*a1 + x_off - a:
                        x_set -= nx*a1
                    if x_set + small < x_off:
                        x_set += nx*a1
         #          if y_set + small > nx*a1 + y_off:
          #             y_set -= ny*a2 
           #        if y_set - small < y_off:
            #           y_set += ny*a2
                    new_atom = atom(x_set, y_set, z_set,\
                                    c.spec, c.q, c.mol)# + mol_counter)
                    ## UPDATE - label atom
                    new_atom.atom = c.mol + mol_counter
                    new_list.append(new_atom)
                    mol_counter += 1
               
    atom_list.extend(new_list)
 
    #mol_counter += 1

    # Add carbon to the list of present types if not already present
    
    flag = False
    for t in type_list:
         if t == 'C':
             flag = True
    if not flag:
        type_list.append('C')    
          
    return atom_list, mol_counter, type_list

def replace_molecule(atom_list, bond_list, angle_list, improper_list, 
                     dihedral_list, mol_count, type_index, mass_index, 
                     add_molecule, bonds, angles, impropers, dihedrals,
                     exclude_list, N):
    """
    Replace random molecules in atom_list with N copies of add_molecule. 
    Molecules with atoms found in exclude list are excluded are not replaced. 
    """
    i = 0
    print "{} molecules to replace!".format(N)  
    ## Assume last atom in atom_list has highest value
    count = atom_list[len(atom_list)-1].atom 
    while i < N:
       
        k = int(random.uniform(0, mol_count))

        print "Replacing molecule {}".format(k)

        exclude_flag = False
 
        ## use original length of atom list when assigning new atom ids
 
        j = 0
        remove_mol = []
        index_list = []
        mol_flag = False 
        for a in atom_list:
            

            if a.mol == k:
                mol_flag = True
                for t in exclude_list:
                    if a.spec == t:
                        print "Molecule excluded!"  
                        exclude_flag = True

                if not exclude_flag:       
                    index_list.append(a.atom) 
                    remove_mol.append(a)
        
           
        if not exclude_flag and mol_flag:
            # Find centre of mass coordinate of molecule
            R = ([[0.0, 0.0, 0.0]])
            M = 0.0
             
            for a in remove_mol:
                
                m = mass_index[type_index[a.spec]]
                M += m   
                R += m*a.return_r_matrix()

            R /= M
            xcoord = R[0,0]
            ycoord = R[0,1]
            zcoord = R[0,2]

 
            for a in remove_mol:
                ## iteratively remove atoms from list
                print "removing atom {}".format(a.atom)
                atom_list.remove(a)

            bond_remove = [] 
            for b in bond_list:
                for a in remove_mol:
                    if b.atom1 == a.atom or b.atom2 == a.atom:
                        bond_remove.append(b)

            angle_remove = []
            for c in angle_list:
                for a in remove_mol:
                    if (c.atom1 == a.atom or c.atom2 == a.atom or 
                        c.atom3 == a.atom):
                        angle_remove.append(c)

#           dihedral_remove = []
 #          for d in dihedral_list:
  #             if d.mol == k:
   #                dihedral_remove.append(d)

   #        improper_remove = []
    #       for e in improper_list:
     #          if e.mol == k: 
      #             improper_remove.append(e) 

            ## Remove tagged bonds, angles etc from respective lists:
            ## Remove duplicate bonds
            bond_remove = list(set(bond_remove))
            for b in bond_remove:
                print "removing bonds for {} and {}".format(b.atom1, 
                                                            b.atom2) 
                bond_list.remove(b)
            ## Remove duplicate angles
            angle_remove = list(set(angle_remove))
            for c in angle_remove:
                print "removing angles for {}, {} and {}".format(c.atom1, 
                                                             c.atom2,
                                                             c.atom3) 
                angle_list.remove(c)
      #     for d in dihedral_remove:
       #        dihedral_list.remove(d)
        #   for e in improper_remove:
         #      improper_list.remove(e)

            ## Add new atoms

                
            for c in angles:
                newangle = angle(c.angle_type, c.atom1+count, c.atom2+count, 
                             c.atom3+count)
                newangle.mol = k  
                angle_list.append(newangle)

            for b in bonds:
                newbond = bond(b.bond_type, b.atom1+count, b.atom2+count) 
                newbond.mol = k  
                bond_list.append(newbond)          

            for d in dihedrals:
                newdihedral = dihedral(d.dihedral_type,
                                       d.atomi+count, d.atomj+count,
                                       d.atomk+count, d.atoml+count)
                newdihedral.mol = k  
                dihedral_list.append(newdihedral)

            for e in impropers:
                newimproper = improper(e.improper_type,
                                       e.atomi+count, e.atomj+count,
                                       e.atomk+count, e.atoml+count)
                newimproper.mol = k  
                improper_list.append(newimproper)

            for a in add_molecule:
                newatom = atom(a.x + xcoord, a.y + ycoord, a.z + zcoord, 
                           a.spec, a.q, k)
                newatom.atom = count + 1   
                atom_list.append(newatom)
                count += 1

            i += 1                 

    return atom_list, bond_list, angle_list, improper_list, dihedral_list, \
           mol_count 


def random_insertion_grid(atom_list, bond_list, angle_list, 
                     improper_list, dihedral_list, mol_count, 
                     type_index, mass_index, molecule, bonds, angles,
                     impropers, dihedrals,
                     xlo, xhi, ylo, yhi, zlo, zhi, rho, r_tol, gx, gy, gz):
    """"
    Inserts \'molecules\' randomly into cell via grid to avoid extraneous 
    searching for overlaps. Appends atom list. Routine assumes empty cell at
    the moment. 
    -Density given in g/cm^3
    -Mass in a.m.u
    -Cell extents in Angstrom
    -Stop at 90% rejection rate
    """
    NA = 6.0221415e23

    m = 0.0

    for a in molecule:  #Get mass from dictionary
        m += mass_index[type_index[a.spec]]

    print "Mass of molecule = %f" % m
   
    V = (xhi-xlo)*(yhi-ylo)*(zhi-zlo)  
    N = int(V*rho*1e-24/(m/NA)) # Number of molecules
    rho = N*(m/NA)/(V*1e-24)
    print "Density = {} g/cm^3".format(rho)
    print "Molecules to be added = %d" % N

    # Check that there are more cells in grid than molecules
    while N > gx*gy*gz:
        gx += 1
        gy += 1
        gz += 1
        print "Resizing grid... gx = {}, gy = {}, gz = {} ".format(gx,gy,gz)
 
    grid = numpy.zeros((gx,gy,gz))


    if atom_list:
        count = atom_list[len(atom_list)-1].atom
    

    else:
        count = 0
    accept_flg = True
    critical_rate = 0.90
    rate_flg = False
    reject = 0
    i = 0   

    Lx = xhi - xlo
    Ly = yhi - ylo
    Lz = zhi - zlo

    l_x = matrix([[Lx, 0.0, 0.0]])
    l_y = matrix([[0.0, Ly, 0.0]])
    l_z = matrix([[0.0, 0.0, Lz]])


 
    while i <= N and not rate_flg:

        xcoord = random.uniform(0.0, Lx)
        ycoord = random.uniform(0.0, Ly)
        zcoord = random.uniform(0.0, Lz)

        xgrid = int(gx*xcoord/Lx)
        ygrid = int(gy*ycoord/Ly)
        zgrid = int(gz*zcoord/Lz)
        
        
        if grid[xgrid, ygrid, zgrid] == 0.0:
            # if nothing on grid point, insert
            mol_count += 1

            grid[xgrid, ygrid, zgrid] += 1.0

            xcoord = xgrid*Lx/gx + xlo
            ycoord = ygrid*Ly/gy + ylo
            zcoord = zgrid*Lz/gz + zlo

            for c in angles:
                newangle = angle(c.angle_type, c.atom1+count, c.atom2+count, 
                             c.atom3+count)
                newangle.mol = mol_count 
                angle_list.append(newangle)

            for b in bonds:
                newbond = bond(b.bond_type, b.atom1+count, b.atom2+count)
                newbond.mol = mol_count 
                bond_list.append(newbond)          

            for d in dihedrals:
                newdihedral = dihedral(d.dihedral_type,
                                       d.atomi+count, d.atomj+count,
                                       d.atomk+count, d.atoml+count)
                newdihedral.mol = mol_count
                dihedral_list.append(newdihedral)

            for e in impropers:
                newimproper = improper(e.improper_type,
                                       e.atomi+count, e.atomj+count,
                                       e.atomk+count, e.atoml+count)
                newimproper.mol = mol_count
                improper_list.append(newimproper)

            for a in molecule:
                newatom = atom(a.x + xcoord, a.y + ycoord, a.z + zcoord, 
                           a.spec, a.q, mol_count) 
                newatom.atom = count + 1 
                atom_list.append(newatom)
                count += 1
             

            i += 1  
        elif grid[xgrid, ygrid, zgrid]:
            reject += 1

    print "Number of tries = {}".format(i + reject -1)
    print "Number of rejections = {}".format(reject)
    print "Rejection rate = {}".format(float(reject)/float(i + reject))

    return atom_list, bond_list, angle_list, improper_list, dihedral_list, \
           mol_count, float(i), m

def random_insertion(atom_list, bond_list, angle_list, 
                     improper_list, dihedral_list, mol_count, 
                     type_index, mass_index, molecule, bonds, angles,
                     impropers, dihedrals,
                     xlo, xhi, ylo, yhi, zlo, zhi, rho, r_tol):
    """"
    Inserts \'molecules\' randomly into cell. Appends atom list.
    -Density given in g/cm^3
    -Mass in a.m.u
    -Cell extents in Angstrom
    -Stop at 90% rejection rate
    """
    NA = 6.0221415e23

    m = 0.0

    for a in molecule:  #Get mass from dictionary
        m += mass_index[type_index[a.spec]]

    print "Mass of molecule = %f" % m
   
    V = (xhi-xlo)*(yhi-ylo)*(zhi-zlo)  
    N = int(V*rho*1e-24/(m/NA)) # Number of molecules
    rho = N*(m/NA)/(V*1e-24)
    print "Density = {} g/cm^3".format(rho)
    print "Molecules to be added = %d" % N

    count = len(atom_list)

    accept_flg = True
    critical_rate = 0.90
    rate_flg = False
    reject = 0
    i = 0   

    Lx = xhi - xlo
    Ly = yhi - ylo
    Lz = zhi - zlo

    l_x = matrix([[Lx, 0.0, 0.0]])
    l_y = matrix([[0.0, Ly, 0.0]])
    l_z = matrix([[0.0, 0.0, Lz]])


 
    while i <= N and not rate_flg:
        accept_flg = True

        xcoord = random.uniform(xlo,xhi)
        ycoord = random.uniform(ylo,yhi)
        zcoord = random.uniform(zlo,zhi)

        r_i = ([[xcoord, ycoord, zcoord]])
        
        j = 0

        while j < len(atom_list) and accept_flg:

            r_ij = r_i - atom_list[j].return_r_matrix()

            if r_ij[0,0] > Lx/2: r_ij -= l_x
            if r_ij[0,0] < -Lx/2: r_ij += l_x
            if r_ij[0,1] > Ly/2: r_ij -= l_y
            if r_ij[0,1] < -Ly/2: r_ij += l_y
            if r_ij[0,2] > Lz/2: r_ij -= l_z
            if r_ij[0,2] < -Lz/2: r_ij += l_z

            r_ij2 = float(r_ij*r_ij.T)   

            if r_ij2 < r_tol*r_tol:
                accept_flg = False
                reject += 1
                rejection_rate = reject/(i + reject)
                if rejection_rate >= critical_rate: rate_flg = True

             
            j += 1
             
            
        if accept_flg: 

            mol_count += 1

            for c in angles:
                newangle = angle(c.angle_type, c.atom1+count, c.atom2+count, 
                             c.atom3+count)
                angle_list.append(newangle)

            for b in bonds:
                newbond = bond(b.bond_type, b.atom1+count, b.atom2+count) 
                bond_list.append(newbond)          

            for d in dihedrals:
                newdihedral = dihedral(d.dihedral_type,
                                       d.atomi+count, d.atomj+count,
                                       d.atomk+count, d.atoml+count)
                dihedral_list.append(newdihedral)

            for e in impropers:
                newimproper = improper(e.improper_type,
                                       e.atomi+count, e.atomj+count,
                                       e.atomk+count, e.atoml+count)
                improper_list.append(newimproper)

            for a in molecule:
                newatom = atom(a.x + xcoord, a.y + ycoord, a.z + zcoord, 
                           a.spec, a.q, mol_count)  
                atom_list.append(newatom)
                count += 1
             

            i += 1        
          #  print "Added molecules: {}".format(i)

    print "Number of tries = {}".format(i + reject)
    print "Number of rejections = {}".format(reject)
    print "Rejection rate = {}".format(float(reject)/float(i + reject))

    return atom_list, bond_list, angle_list, improper_list, dihedral_list, \
           mol_count, float(i), m

def offset_to_hole(atom_list, a, Lx, Ly):
    "Offsets lattice putting pore in centre"
    rt3 = math.sqrt(3.0)
    xc = 11*a/2
    yc = 2*rt3*a

    for i in atom_list:
        i.x -= xc - Lx/2
        i.y -= yc - Ly/2
        # apply pbc
        if i.x < 0:
            i.x += Lx
        if i.x > Lx:
            i.x -= Lx 
        if i.y < 0:
            i.y += Ly
        if i.y > Ly:
            i.y -= Ly 
    
    return atom_list

def apply_pbc(atom_list, Lx, Ly, Lz):
    """
    Applies periodic boundary conditions to all atoms in list
    """
    small = 1.0e-9

    for i in atom_list:
        # apply pbc
        if i.x - small < 0:
            i.x += Lx
        if i.x + small > Lx:
            i.x -= Lx 
        if i.y - small < 0:
            i.y += Ly
        if i.y + small > Ly:
            i.y -= Ly 
        if i.z - small < 0:
            i.z += Lz
        if i.z + small > Lz:
            i.z -= Lz 
    
    return atom_list

def create_hole(atom_list, l, shell, nx, ny, Lx, Ly):
    """
    Removes shells of atoms from graphene sheet with origin at centre 
    of cell (off-atom). Can remove up to 3 shells. 
    """

    N_n = 0    
    tol = 0.01*l

    if shell > 3:
        print "ERROR: Can only remove up to 3 carbon shells for hole.\n"
        exit()  

    if shell < 1: 
        print "ERROR: Shell number must be >=1\n"
        exit()  

    ## Check system size
    if (shell == 1 or shell == 2) and (nx < 2 or ny <2):
        print "ERROR: Too many atoms to remove. Hole too large!\n"
        exit()
    if (shell == 3) and (nx<3 or ny<3):
        print "ERROR: Too many atoms to remove. Hole too large!\n"

    if shell == 3:
        N_n =  12
    else:
        N_n = 6

    N_count = 0

    ## Assume that cell is already centred (off atom)
    if shell == 1: ## Remove single shell     
        r = l+tol 
    if shell == 2: ## Remove 2nd shell     
        r = 2*l+tol 
    if shell == 3: ## Remove 3rd shell     
        r = math.sqrt(7)*l+tol 

    r2 = r*r
    print "r^2 = %f" % r2

    i = 0
    dellist = []

    newlist = []
 #  while N_count < N_n:  
    for a in atom_list:
        x = a.x - Lx/2
        y = a.y - Ly/2
        s2 = x*x + y*y
#       print x, y, s2
                
        if s2 > r2:
            newlist.append(a)    
        else: # mark for deletion
            N_count += 1

    print "%d atoms removed" % N_count 

    return newlist

def functionalise_sheet2(atom_list, bond_list, angle_list, dihedral_list,\
                        mol_count, bond_types, angle_types, dihedral_types, \
                        type_list, group, a, p, seed, lx, ly, lz):
    """
    Functionalises graphene sheet. 
    Randomly selects carbon atoms for functionalisation with probability p.
    Then randomly selects side of (x-y) sheet with equal probability.  
    """
    tol = 0.01
    a += tol
    a2 = a*a
    pi = 4.0*math.atan(1.0)
    random.seed(seed)
    new_list = []
    newcount = len(atom_list)

    if group == 'COH': ## Test with COH groups
        ## Define bond type 3 and 4 angle type 2 and 3
        ## Values obtained from Mooney et al., Chemical Physics Letters,
        ## 294 (1998) 135-142                
        group_bonds = [3,4]
        group_angles = [3,4]
        group_dihedrals = [1]
        group_types = ['C_COH', 'O_COH', 'H_COH']
        l_CO = 1.364 
        l_OH = 0.960 
        ang_CCO = 90   # is this angle right?
#        ang_COH = 180 
        ang_COH = 110.5 
        
    if group == 'test': ## Test with test group
        ## Define bond type 3 
        group_bonds = [3]
        group_angles = [3]
        group_types = ['C_COH', 'O_COH']
        l_CO = 1.364 
        l_OH = 0.960 
        ang_CCO = 90
        ang_COH = 110.5 
 
    ## Define +/- normal vectors 
    n_p = vec3(0.0, 0.0, 1.0)
    n_m = vec3(0.0, 0.0, -1.0)
   
     
    icount = 0 
    for i in atom_list:
        icount += 1
        rand = random.uniform(0.0, 1.0)
        ## functionalise atom
        if i.spec == 'C' and rand <= p:
            side_rand = random.uniform(0.0, 1.0)
            if side_rand < 0.5:
                norm = n_p
            else:
                norm = n_m
            ## Test with OH group - defines angle tyes 2 and 3 
            if group == 'COH':

                CCO_angle_list = []
                
                r_CO = vec3(i.x, i.y, i.z + l_CO)

                ## find nearest neighbour for angle
                nn_count = 0
                while nn_count < 3:
                    jcount = 0
                    neighbour_list = []
                    for j in atom_list: 
                        x_i = i.x 
                        y_i = i.y
                        x_j = j.x
                        y_j = j.y 
                
                        x_ij = x_i - x_j 
                        y_ij = y_i - y_j

                        if x_ij < -lx/2: x_ij = x_i - x_j + lx              
                        if x_ij > lx/2: x_ij = x_i - x_j - lx              
                        if y_ij < -ly/2: y_ij = y_i - y_j + ly              
                        if y_ij > ly/2: y_ij = y_i - y_j - ly              
                        r_ij = vec3(x_ij, y_ij, i.z-j.z)                 
                        ## Apply PBC
                        jcount += 1
                        if abs(r_ij) < a and jcount != icount:
                            neighbour_list.append(jcount) 
                            nn_count += 1
          #     while nn_count == 0:
           #        jcount = 0
            #       for j in atom_list: 
             #          r_ij = vec3(i.x-j.x, i.y-j.y, i.z-j.z)
              #         jcount += 1
               #        if abs(r_ij) < a: 
                #           nn_count += 1

                ## C
                atom_list[icount-1] = atom(i.x, i.y, i.z, 'C_COH', 0.2, i.mol)
                atom_list[icount-1].atom = icount
                ## O
                newcount += 1
                O_COH = atom(i.x, i.y, i.z + l_CO*norm.z,'O_COH', -0.64, i.mol)
                O_COH.atom = newcount

                CCO_angle_list.append(angle(3, neighbour_list[0], \
                                      icount, newcount))
                CCO_angle_list.append(angle(3, neighbour_list[1], \
                                      icount, newcount))
         #      CCO_angle_list.append(angle(3, neighbour_list[2], \
          #                           icount, newcount))

                angle_list.extend(CCO_angle_list)
   
                CO_bond = bond(3, icount, newcount)
                bond_list.append(CO_bond)

                new_list.append(O_COH)
                 
                ## H                       
                newcount += 1
                rad_COH = ang_COH*pi/180 
                z_OH = l_OH*math.cos(pi-rad_COH)*norm.z
                y_OH = l_OH*math.sin(pi-rad_COH)
                H_COH = atom(i.x, i.y + y_OH, i.z + l_CO*norm.z + z_OH, \
                             'H_COH', 0.44, i.mol)
                H_COH.atom = newcount
                new_list.append(H_COH)

                OH_bond = bond(4, newcount-1, newcount)
                bond_list.append(OH_bond)
                
                COH_angle = angle(4, icount, newcount-1, newcount)   
                angle_list.append(COH_angle)

                ## Dihedrals
                CCOH_dihedral = dihedral(1, neighbour_list[0], icount, \
                                         newcount-1, newcount)
                dihedral_list.append(CCOH_dihedral)
            
            if group == 'test':
                
                CCO_angle_list = []
             
                r_CO = vec3(i.x, i.y, i.z + l_CO)

                ## find nearest neighbour for angle
                nn_count = 0
                while nn_count < 3:
                    jcount = 0
                    neighbour_list = []
                    for j in atom_list: 
                        r_ij = vec3(i.x-j.x, i.y-j.y, i.z-j.z)
                        jcount += 1
                        if abs(r_ij) < a and jcount != icount:
                            neighbour_list.append(jcount) 
                            nn_count += 1

                ## C
                atom_list[icount-1] = atom(i.x, i.y, i.z, 'C_COH', 0.0, i.mol)

                ## O
                newcount += 1
                O_COH = atom(i.x, i.y, i.z + l_CO*norm.z, 'O_COH', 0.0, i.mol)
                print "%d %d" % (icount, newcount)
                print "x : %f\ty : %f\tz: %f" % (i.x, i.y, i.z + l_CO*norm.z)
                new_list.append(O_COH)

                ## Bonds 
                CO_bond = bond(3, icount, newcount)
                bond_list.append(CO_bond)
                
                ## Angles
                CCO_angle_list.append(angle(3, neighbour_list[0], \
                                      icount, newcount))
                CCO_angle_list.append(angle(3, neighbour_list[1], \
                                      icount, newcount))
                angle_list.extend(CCO_angle_list)
               
                ## Dihedrals
                CCOH_dihedral = dihedral(1, neighbour_list[0], icount, \
                                         newcount-1, newcount)
                dihedral_list.append(CCOH_dihedral)
                 

    atom_list.extend(new_list)                                   
    bond_types.extend(group_bonds)
    angle_types.extend(group_angles)
    dihedral_types.extend(group_dihedrals)
    type_list.extend(group_types)
                
    return atom_list, bond_list, angle_list, dihedral_list, mol_count, \
           bond_types, angle_types, dihedral_types, type_list
 
def functionalise_sheet(atom_list, bond_list, angle_list, mol_count, \
                        bond_types, angle_types, type_list, group, a, p, seed,\
                        lx, ly, lz):
    """
    Functionalises graphene sheet. 
    Randomly selects carbon atoms for functionalisation with probability p.
    Then randomly selects side of (x-y) sheet with equal probability.  
    """
    tol = 0.01
    a += tol
    a2 = a*a
    pi = 4.0*math.atan(1.0)
    random.seed(seed)
    new_list = []
    newcount = len(atom_list)

    if group == 'COH': ## Test with COH groups
        ## Define bond type 3 and 4 angle type 2 and 3
        ## Values obtained from Mooney et al., Chemical Physics Letters,
        ## 294 (1998) 135-142                
        group_bonds = [3,4]
        group_angles = [3,4]
        group_types = ['C_COH', 'O_COH', 'H_COH']
        l_CO = 1.364 
        l_OH = 0.960 
        ang_CCO = 90
#        ang_COH = 180 
        ang_COH = 110.5 
        
    if group == 'test': ## Test with test group
        ## Define bond type 3 
        group_bonds = [3]
        group_angles = [3]
        group_types = ['C_COH', 'O_COH']
        l_CO = 1.364 
        l_OH = 0.960 
        ang_CCO = 90
        ang_COH = 110.5 
 
    ## Define +/- normal vectors 
    n_p = vec3(0.0, 0.0, 1.0)
    n_m = vec3(0.0, 0.0, -1.0)
   
     
    icount = 0 
    for i in atom_list:
        icount += 1
        rand = random.uniform(0.0, 1.0)
        ## functionalise atom
        if i.spec == 'C' and rand <= p:
            side_rand = random.uniform(0.0, 1.0)
            if side_rand < 0.5:
                norm = n_p
            else:
                norm = n_m
            ## Test with OH group - defines angle tyes 2 and 3 
            if group == 'COH':

                CCO_angle_list = []
                
                r_CO = vec3(i.x, i.y, i.z + l_CO)

                ## find nearest neighbour for angle
                nn_count = 0
                while nn_count < 3:
                    jcount = 0
                    neighbour_list = []
                    for j in atom_list: 
                        x_i = i.x 
                        y_i = i.y
                        x_j = j.x
                        y_j = j.y 
                
                        x_ij = x_i - x_j 
                        y_ij = y_i - y_j

                        if x_ij < -lx/2: x_ij = x_i - x_j + lx              
                        if x_ij > lx/2: x_ij = x_i - x_j - lx              
                        if y_ij < -ly/2: y_ij = y_i - y_j + ly              
                        if y_ij > ly/2: y_ij = y_i - y_j - ly              
                        r_ij = vec3(x_ij, y_ij, i.z-j.z)                 
                        ## Apply PBC
                        jcount += 1
                        if abs(r_ij) < a and jcount != icount:
                            neighbour_list.append(jcount) 
                            nn_count += 1
          #     while nn_count == 0:
           #        jcount = 0
            #       for j in atom_list: 
             #          r_ij = vec3(i.x-j.x, i.y-j.y, i.z-j.z)
              #         jcount += 1
               #        if abs(r_ij) < a: 
                #           nn_count += 1

                ## C
                atom_list[icount-1] = atom(i.x, i.y, i.z, 'C_COH', 0.2, i.mol)

                ## O
                newcount += 1
                O_COH = atom(i.x, i.y, i.z + l_CO*norm.z,'O_COH', -0.64, i.mol)

                CCO_angle_list.append(angle(3, neighbour_list[0], \
                                      icount, newcount))
                CCO_angle_list.append(angle(3, neighbour_list[1], \
                                      icount, newcount))
         #      CCO_angle_list.append(angle(3, neighbour_list[2], \
          #                           icount, newcount))

                angle_list.extend(CCO_angle_list)
   
                CO_bond = bond(3, icount, newcount)
                bond_list.append(CO_bond)

                new_list.append(O_COH)
                 
                ## H                       
                newcount += 1
                rad_COH = ang_COH*pi/180 
                z_OH = l_OH*math.cos(pi-rad_COH)*norm.z
                y_OH = l_OH*math.sin(pi-rad_COH)
                H_COH = atom(i.x, i.y + y_OH, i.z + l_CO*norm.z + z_OH, \
                             'H_COH', 0.44, i.mol)
                new_list.append(H_COH)

                OH_bond = bond(4, newcount-1, newcount)
                bond_list.append(OH_bond)
                
                COH_angle = angle(4, icount, newcount-1, newcount)   
                angle_list.append(COH_angle)
            
            if group == 'test':
                
                CCO_angle_list = []
             
                r_CO = vec3(i.x, i.y, i.z + l_CO)

                ## find nearest neighbour for angle
                nn_count = 0
                while nn_count < 3:
                    jcount = 0
                    neighbour_list = []
                    for j in atom_list: 
                        r_ij = vec3(i.x-j.x, i.y-j.y, i.z-j.z)
                        jcount += 1
                        if abs(r_ij) < a and jcount != icount:
                            neighbour_list.append(jcount) 
                            nn_count += 1

                ## C
                atom_list[icount-1] = atom(i.x, i.y, i.z, 'C_COH', 0.0, i.mol)

                ## O
                newcount += 1
                O_COH = atom(i.x, i.y, i.z + l_CO*norm.z, 'O_COH', 0.0, i.mol)
                print "%d %d" % (icount, newcount)
                print "x : %f\ty : %f\tz: %f" % (i.x, i.y, i.z + l_CO*norm.z)
                new_list.append(O_COH)

                ## Bonds 
                CO_bond = bond(3, icount, newcount)
                bond_list.append(CO_bond)
                
                ## Angles
                CCO_angle_list.append(angle(3, neighbour_list[0], \
                                      icount, newcount))
                CCO_angle_list.append(angle(3, neighbour_list[1], \
                                      icount, newcount))
                angle_list.extend(CCO_angle_list)

    atom_list.extend(new_list)                                   
    bond_types.extend(group_bonds)
    angle_types.extend(group_angles)
    type_list.extend(group_types)
                
    return atom_list, bond_list, angle_list, mol_count, \
           bond_types, angle_types, type_list 

def functionalise_edge3(atom_list, bond_list, angle_list, improper_list, \
                        dihedral_list, mol_count, bond_types, angle_types, \
                        improper_types, dihedral_types, type_list, \
                        a, p, lx, ly, pbc, group, ydist):
    """
    This version includes dihedrals in COH groups (Phenol)

    Identifies edge atoms using nearest neighbour search and uses vectors
    between neighbours to set up positions for functional groups
    """
  
    tol = 0.01
    pi = 4.0*math.atan(1.0)

    if group == 'COH': ## Test with COH groups
        ## Define bond type 3 and 4 angle type 2 and 3
        ## Values obtained from Mooney et al., Chemical Physics Letters,
        ## 294 (1998) 135-142
        ## and Konatham et al, Langmuir 2013, Supporting Information (angles)
        group_bonds = [3,4]
        group_angles = [2,3]
        group_improper = [1]
        group_dihedral = [1]
        group_types = ['C_COH', 'O_COH', 'H_COH']
        l_CO = 1.364
        l_OH = 0.960
        ang_CCO = 120
#        ang_COH = 180 
        ang_COH = 113*pi/180 # changed so that it's from same source as k vals
    
    if group == 'CCOO': ## Test with CCOO groups
        ## Define bond type 3,4,5 angle type 2,3,4,5 (these are just labels, shouldn't really matter)
        ## Values from Konatham et al, Langmuir 2013, Supporting Information
        group_bonds = [3,4,5]
        group_angles = [2,3,4,5]
        group_improper = [1]
        group_dihedral = [1,2]
        group_types = ['C1_CCOO', 'C2_CCOO', 'O1_CCOO', 'O2_CCOO']
        l_CC = 1.41
        l_CO1 = 1.22
        l_CO2 = 1.25
        ang_CCC = 120
        ang_OCO = 126  # what to do with this one?
        ang_OCO_r = 126*pi/180
        ang_CCO1 = 120*pi/180
        ang_CCO2 = 117*pi/180#246*pi/180#

    if group == 'CCOOH': ## Carboxyl functional group
        ## Define bond type 3,4,5,6 angle type 2,3,4,5,6
        ## Values from OPLS force fiels, BOSS distribution
        ## OPLS All-Atom Parameters for Organic Molecules, Ions, Peptides & Nucleic Acids: Jan  2007
        group_bonds = [3,4,5,6]
        group_angles = [2,3,4,5,6]
        group_improper = [1]
        group_dihedral = [1,2]
        group_types = ['C1_CCOOH', 'C2_CCOOH', 'O1_CCOOH', 'O2_CCOOH', 'H_CCOOH']
        l_CC = 1.529
        l_CO1 = 1.229
        l_CO2 = 1.364
        l_OH = 0.945
        ang_CCC = 111.1
        ang_OCO = 121
        ang_OCO_r = 121*pi/180
        ang_CCO1 = 120.4*pi/180
        ang_CCO2 = 108*pi/180
        ang_COH = 113*pi/180



    if group == 'O_test':
        # Replaces carbons with oxygens
        group_types = ['C_COH', 'O_COH', 'H_COH', 'C_CH', 'H_CH']
        group_bonds = []
        group_angles = []                  
        group_improper = []
        group_dihedral = []

    if group == 'CH': ## Test with hydrogenated groups
        ## Define bond type 3, angle type 2, improper type 1
        group_bonds = [3]
        group_angles = [2]
        group_improper = [1]
        group_dihedral = [1]
        group_types = ['C_CH', 'H_CH']
        l_CH = 1.08
        ang_CCH = 120

    if group == 'CH_only': ## For full functionalisation
        ## Define bond type 3, angle type 2, improper type 1
        group_bonds = [3]
        group_angles = [2]
        group_improper = [1]
        group_dihedral = [1]
        group_types = ['C_CH', 'H_CH']
        l_CH = 1.08
        ang_CCH = 120




    a2 = (a+tol)*(a+tol)
    icount = 0
    newcount = len(atom_list)
    new_list = []
    neigh_count = 1
    func_list = []

    for i in atom_list:
        rand = random.uniform(0.0, 1.0)
        icount += 1

        neighbour_vectors = []
        neighbour_list = []
        jcount = 0
        for j in atom_list:
            jcount += 1
            if icount != jcount:
                x_i = i.x 
                y_i = i.y
                x_j = j.x
                y_j = j.y 
                
                x_ij = x_i - x_j 
                y_ij = y_i - y_j
                

                if pbc == 'y': 
                    if x_ij < -lx/2: x_ij = x_i - x_j + lx              
                    if x_ij > lx/2: x_ij = x_i - x_j - lx              
                    if y_ij < -ly/2: y_ij = y_i - y_j + ly              
                    if y_ij > ly/2: y_ij = y_i - y_j - ly              
  
                r_ij = vec3(x_ij, y_ij, i.z-j.z)                 
                 
                if r_ij**2 < a2:
                    neighbour_vectors.append(r_ij)
                    neighbour_list.append(jcount)
    
    
        if len(neighbour_vectors) == 2: #if fewer than 3 nearest neighbours
            
            ## Define vector for functional group
            r_func = vec3(0.0,0.0,0.0)
            for r in neighbour_vectors:
                r_func += r
            r_func /= math.sqrt(r_func**2)
            
            if group == 'CH':
                
                
                ## Test with CH group
    

                
                # check if neighbouring edge carbon has been functionalised
                
                C_new_coords = np.around(np.array([i.x, i.y, i.z]),6)
                
                C_new_plus = np.around(np.array([i.x, i.y + ydist, i.z]),6)
                
                if i.y + ydist > ly:
                    C_new_plus = np.around(np.array([i.x, i.y + ydist - ly, i.z]),6)
                
                
                C_new_minus = np.around(np.array([i.x, i.y - ydist, i.z ]),6)
                if i.y - ydist < 0:
                    C_new_minus = np.around(np.array([i.x, i.y - ydist + ly, i.z]),6)


                
                # Check if neighbouring C has been functionalised
                if any((C_new_plus == x).all() for x in func_list) == True or any((C_new_minus == x).all() for x in func_list) == True:
                    print 'A neighbour has already been functionalised. Skip atom.'
            
                else:
                    func_list.append(C_new_coords)
                    newcount += 1
                
                    ## Replace Carbon
                    atom_list[icount-1] = atom(i.x, i.y, i.z, 'C_CH', -0.115, i.mol)
                    atom_list[icount-1].atom = icount
                    
                    
                    r_func *= 1.08
                    ## add hydrogen
                    H_CH = atom(r_func.x + i.x, r_func.y + i.y, r_func.z + i.z, \
                            'H_CH', 0.115, i.mol)
                    H_CH.atom = newcount
                    
                    print icount, newcount
                    
                    new_list.append(H_CH)
                    # CH Bond
                    bond_list.append(bond(3, icount, newcount))
                    # Angle
                    angle_list.append(angle(2, newcount, icount, neighbour_list[0]))
                    # Improper
                    newimproper = improper(1, newcount, neighbour_list[0], \
                                neighbour_list[1], icount)
                    improper_list.append(newimproper)


            if group == 'CH_only':
    
    
                ## CH on every carbon
    
    
                newcount += 1
        
                ## Replace Carbon
                atom_list[icount-1] = atom(i.x, i.y, i.z, 'C_CH', -0.115, i.mol)
                atom_list[icount-1].atom = icount
        
        
                r_func *= 1.08
                ## add hydrogen
                H_CH = atom(r_func.x + i.x, r_func.y + i.y, r_func.z + i.z, \
                    'H_CH', 0.115, i.mol)
                H_CH.atom = newcount
                    
 
                new_list.append(H_CH)
                # CH Bond
                bond_list.append(bond(3, icount, newcount))
                # Angle
                angle_list.append(angle(2, newcount, icount, neighbour_list[0]))
                # Improper
                newimproper = improper(1, newcount, neighbour_list[0], \
                                           neighbour_list[1], icount)
                improper_list.append(newimproper)


            if group == 'CCOOH':
    
                CCC_angle_list = []
                OCO_angle_list = []
    
                r_CC = r_func*l_CC
    
    
    
                # check if neighbouring edge carbon has been functionalised
    
                C_new_coords = np.around(np.array([i.x, i.y, i.z]),6)
    
                C_new_plus = np.around(np.array([i.x, i.y + ydist, i.z]),6)
    
                if i.y + ydist > ly:
                    C_new_plus = np.around(np.array([i.x, i.y + ydist - ly, i.z]),6)
    
    
                C_new_minus = np.around(np.array([i.x, i.y - ydist, i.z ]),6)
                if i.y - ydist < 0:
                    C_new_minus = np.around(np.array([i.x, i.y - ydist + ly, i.z]),6)
    
    
                # Check if neighbouring C has been functionalised
                if any((C_new_plus == x).all() for x in func_list) == True or any((C_new_minus == x).all() for x in func_list) == True:
                    print 'A neighbour has already been functionalised. Skip atom.'
    
                else:
                    func_list.append(C_new_coords)
        
                    ## C (GS)
                    atom_list[icount-1] = atom(i.x, i.y, i.z, 'C1_CCOOH', 0.08, i.mol)
                    atom_list[icount-1].atom = icount
        
                    ## C (FG)
                    newcount += 1
                    C2_CCOOH = atom(i.x + r_CC.x, i.y + r_CC.y, i.z + r_CC.z ,\
                       'C2_CCOOH', 0.55, i.mol)
                    C2_CCOOH.atom = newcount
                       
                    # Improper, what kind of impropers do I need?
                    newimproper = improper(1, newcount, neighbour_list[0], \
                                              neighbour_list[1], icount)
                    improper_list.append(newimproper)
                                              
                    # Angle CCC, what about neighbour list?
                    CCC_angle = angle(2, neighbour_list[1], icount, newcount)
                    angle_list.append(CCC_angle)
                                              
                    
                                              
                    # Bond CC
                    CC_bond = bond(3, icount, newcount)
                    bond_list.append(CC_bond)
                                              
                                           
                    new_list.append(C2_CCOOH)
                                              
                    #-----------------------------
                    ## O1 (double bond)
                    newcount += 1
                    sign = 1.0
                                              
                    x_CO1 = l_CO1*math.cos(pi-ang_CCO1)
                    y_CO1 = l_CO1*math.sin(pi-ang_CCO1)*sign
                    z_CO1 = 0.0
                    r_CO1 = vec3(x_CO1, y_CO1, z_CO1)
                                              
                    #determine angle of CC bond with x axis
                    n_x = vec3(1, 0, 0)
                                              
                    ang_CCZ = math.acos((r_CC*n_x)/abs(r_CC))
                                              
                    if r_CC.x > 0 + tol and r_CC.y > 0 + tol:
                        ang_CCZ = 2*pi - ang_CCZ
                                              
                    if r_CC.x < 0 - tol and r_CC.y > 0 + tol:
                        ang_CCZ = 2*pi - ang_CCZ
                                              
                                              
                    ## z-axis rotation matrix
                    R_z = numpy.matrix([\
                            [math.cos(ang_CCZ), -math.sin(ang_CCZ), 0.0],\
                            [math.sin(ang_CCZ), math.cos(ang_CCZ), 0.0],\
                            [0.0, 0.0, 1.0]])
                    r_CO1_mat = numpy.matrix([[r_CO1.x, r_CO1.y, r_CO1.z]])
                                                                  
                    r_CO1_mat = r_CO1_mat*R_z
                                                                  
                    r_CO1 = vec3(r_CO1_mat[0,0], r_CO1_mat[0,1], r_CO1_mat[0,2])
                                                                  
                    O1_CCOOH = atom(i.x + r_CC.x + r_CO1.x,
                                i.y + r_CC.y + r_CO1.y,
                                i.z + r_CC.z + r_CO1.z , \
                                'O1_CCOOH', -0.50, i.mol)
                    O1_CCOOH.atom = newcount
                                                                                 
                    new_list.append(O1_CCOOH)
                                                                                 
                                                                                 
                    # --------------
                                                                                 
                    ## O2 (H)
                    newcount += 1
                    sign = 1.0
                        
                    # try rotating CO1 vector by OCO angle
                    R_OCO = numpy.matrix([\
                                [math.cos(ang_OCO_r), -math.sin(ang_OCO_r), 0.0],\
                                [math.sin(ang_OCO_r), math.cos(ang_OCO_r), 0.0],\
                                [0.0, 0.0, 1.0]])
                                                                                                       
                    r_CO2_mat = r_CO1_mat*R_OCO
                    r_CO2 = vec3(r_CO2_mat[0,0], r_CO2_mat[0,1], r_CO2_mat[0,2])
                                                                                                       
                    O2_CCOOH = atom(i.x + r_CC.x + r_CO2.x,
                                   i.y + r_CC.y + r_CO2.y,
                                   i.z + r_CC.z + r_CO2.z , \
                                   'O2_CCOOH', -0.58, i.mol)
                                                                                                                      
                    O2_CCOOH.atom = newcount
                                                                                                                      
                    new_list.append(O2_CCOOH)
                                                                                                                      
                    #-----------
                    
                    # Angle OCO
                    OCO_angle = angle(3, newcount-1, newcount-2, newcount)
                    angle_list.append(OCO_angle)
                                                                                                                      
                    # Bond CO1 (double)
                    CO1_bond = bond(4, newcount-2, newcount-1)
                    bond_list.append(CO1_bond)
                                                                                                                      
                    # Bond CO2 (H)
                    CO2_bond = bond(5, newcount-2, newcount)
                    bond_list.append(CO2_bond)
                                                                                                                      
                    CCO1_angle = angle(4, icount, newcount-2, newcount-1)
                    angle_list.append(CCO1_angle)
                                                                                                                      
                    CCO2_angle = angle(5, icount, newcount-2, newcount)
                    angle_list.append(CCO2_angle)
                    
                    #------------
                    
                    ## H
                    
                    newcount += 1
                    sign = 1.0
                    
                    x_OH = l_OH*math.cos(pi-ang_COH)
                    y_OH = l_OH*math.sin(pi-ang_COH)*sign
                    z_OH = 0.0
                    r_OH = vec3(x_OH, y_OH, z_OH)
                    
                    #determine angle of CC bond with x axis
                    n_x = vec3(1, 0, 0)
                    
                    ang_COZ = math.acos((r_CO2*n_x)/abs(r_CO2))
                    
                    if r_CO2.x > 0 + tol and r_CO2.y > 0 + tol:
                        ang_COZ = 2*pi - ang_COZ
                    
                    if r_CO2.x < 0 - tol and r_CO2.y > 0 + tol:
                        ang_COZ = 2*pi - ang_COZ
                    
                    
                    ## z-axis rotation matrix
                    R_z = numpy.matrix([\
                            [math.cos(ang_COZ), -math.sin(ang_COZ), 0.0],\
                            [math.sin(ang_COZ), math.cos(ang_COZ), 0.0],\
                            [0.0, 0.0, 1.0]\
                                        ])
                    r_OH_mat = numpy.matrix([[r_OH.x, r_OH.y, r_OH.z]])
                                        
                    r_OH_mat = r_OH_mat*R_z
                                        
                    r_OH = vec3(r_OH_mat[0,0], r_OH_mat[0,1], r_OH_mat[0,2])
                                        
                    H_CCOOH = atom(i.x + r_CC.x + r_CO2.x + r_OH.x,
                                    i.y + r_CC.y + r_CO2.y + r_OH.y,
                                    i.z + r_CC.z + r_CO2.z + r_OH.z, \
                                    'H_CCOOH', 0.45, i.mol)
                    H_CCOOH.atom = newcount
                                                        
                    new_list.append(H_CCOOH)
                    
                    
                    #------------
                    
                    # OH bond
                    
                    OH_bond = bond(6, newcount-1, newcount)
                    bond_list.append(OH_bond)
                    
                    # COH angle
                    
                    COH_angle = angle(6, newcount-3, newcount-1, newcount)
                    angle_list.append(COH_angle)
                    
                    #------------
                                                                                                                      
                    # Dihedral (CHECK ATOMS)
                    
                    CCCO1_dihedral = dihedral(1, neighbour_list[0], icount,
                                              newcount-3, newcount-2)
                    dihedral_list.append(CCCO1_dihedral)
                                                                                                                                                
                    CCCO2_dihedral = dihedral(2, neighbour_list[1], icount,
                                            newcount-3, newcount-1)
                    dihedral_list.append(CCCO2_dihedral)



            if group == 'CCOO':
                
                CCC_angle_list = []
                OCO_angle_list = []
                
                r_CC = r_func*l_CC
                

                
                # check if neighbouring edge carbon has been functionalised
                
                C_new_coords = np.around(np.array([i.x, i.y, i.z]),6)
                
                C_new_plus = np.around(np.array([i.x, i.y + ydist, i.z]),6)
                
                if i.y + ydist > ly:
                    C_new_plus = np.around(np.array([i.x, i.y + ydist - ly, i.z]),6)
                
                
                C_new_minus = np.around(np.array([i.x, i.y - ydist, i.z ]),6)
                if i.y - ydist < 0:
                    C_new_minus = np.around(np.array([i.x, i.y - ydist + ly, i.z]),6)

                
                # Check if neighbouring C has been functionalised
                if any((C_new_plus == x).all() for x in func_list) == True or any((C_new_minus == x).all() for x in func_list) == True:
                    print 'A neighbour has already been functionalised. Skip atom.'
                
                else:
                    func_list.append(C_new_coords)
                    
                    ## C (GS)
                    atom_list[icount-1] = atom(i.x, i.y, i.z, 'C1_CCOO', 0.1, i.mol)
                    atom_list[icount-1].atom = icount
                
                    ## C (FG)
                    newcount += 1
                    C2_CCOO = atom(i.x + r_CC.x, i.y + r_CC.y, i.z + r_CC.z ,\
                                   'C2_CCOO', 0.70, i.mol)
                    C2_CCOO.atom = newcount
                             
                    # Improper, what kind of impropers do I need?
                    newimproper = improper(1, newcount, neighbour_list[0], \
                            neighbour_list[1], icount)
                    improper_list.append(newimproper)
                                                    
                    # Angle CCC, what about neighbour list?
                    CCC_angle = angle(2, neighbour_list[1], icount, newcount)
                    angle_list.append(CCC_angle)
                
                    # Angle OCO
                    OCO_angle = angle(3, newcount+2, newcount, newcount+1)
                    angle_list.append(OCO_angle)
                                                    
                    # Bond CC
                    CC_bond = bond(3, icount, newcount)
                    bond_list.append(CC_bond)
                
                                                    
                    new_list.append(C2_CCOO)
                
                    #-----------------------------
                    ## O1 (double bond)
                    newcount += 1
                    sign = 1.0
                    
                    x_CO1 = l_CO1*math.cos(pi-ang_CCO1)
                    y_CO1 = l_CO1*math.sin(pi-ang_CCO1)*sign
                    z_CO1 = 0.0
                    r_CO1 = vec3(x_CO1, y_CO1, z_CO1)
                                                    
                    #determine angle of CC bond with x axis
                    n_x = vec3(1, 0, 0)
                                                    
                    ang_CCZ = math.acos((r_CC*n_x)/abs(r_CC))
                                                    
                    if r_CC.x > 0 + tol and r_CC.y > 0 + tol:
                        ang_CCZ = 2*pi - ang_CCZ
                                                    
                    if r_CC.x < 0 - tol and r_CC.y > 0 + tol:
                        ang_CCZ = 2*pi - ang_CCZ
                                                    
                                                    
                    ## z-axis rotation matrix
                    R_z = numpy.matrix([\
                            [math.cos(ang_CCZ), -math.sin(ang_CCZ), 0.0],\
                            [math.sin(ang_CCZ), math.cos(ang_CCZ), 0.0],\
                            [0.0, 0.0, 1.0]])
                    r_CO1_mat = numpy.matrix([[r_CO1.x, r_CO1.y, r_CO1.z]])
                                                                        
                    r_CO1_mat = r_CO1_mat*R_z
                                                                        
                    r_CO1 = vec3(r_CO1_mat[0,0], r_CO1_mat[0,1], r_CO1_mat[0,2])
                                                                        
                    O1_CCOO = atom(i.x + r_CC.x + r_CO1.x,
                             i.y + r_CC.y + r_CO1.y,
                             i.z + r_CC.z + r_CO1.z , \
                            'O1_CCOO', -0.80, i.mol)
                    O1_CCOO.atom = newcount

                    new_list.append(O1_CCOO)
                
                
                    # --------------
                
                    ## O2 (negative)
                    newcount += 1
                    sign = 1.0
                
                    '''x_CO2 = l_CO2*math.cos(pi-ang_CCO2)
                    y_CO2 = l_CO2*math.sin(pi-ang_CCO2)*sign
                    z_CO2 = 0.0
                    r_CO2 = vec3(x_CO2, y_CO2, z_CO2)'''
                
                    '''#determine angle of CO bond with x axis
                    n_x = vec3(1, 0, 0)
                
                    ang_COZ = math.acos((r_CO*n_x)/abs(r_CO))
                
                    if r_CO.x > 0 + tol and r_CO.y > 0 + tol:
                        ang_COZ = 2*pi - ang_COZ
                
                    if r_CO.x < 0 - tol and r_CO.y > 0 + tol:
                        ang_COZ = 2*pi - ang_COZ
                
                
                    ## z-axis rotation matrix
                    R_z = numpy.matrix([\
                        [math.cos(ang_COZ), -math.sin(ang_COZ), 0.0],\
                        [math.sin(ang_COZ), math.cos(ang_COZ), 0.0],\
                        [0.0, 0.0, 1.0]])'''
                    '''r_CO2_mat = numpy.matrix([[r_CO2.x, r_CO2.y, r_CO2.z]])
                                    
                    r_CO2_mat = r_CO2_mat*R_z
                                    
                    r_CO2 = vec3(r_CO2_mat[0,0], r_CO2_mat[0,1], r_CO2_mat[0,2])
                
                    print "o1", i.x + r_CC.x + r_CO1.x
                    print "o2", i.x + r_CC.x + r_CO2.x
                                    
                    O2_CCOO = atom(i.x + r_CC.x + r_CO2.x,
                            i.y + r_CC.y + r_CO2.y,
                            i.z + r_CC.z + r_CO2.z , \
                            'O2_CCOO', -1.0, i.mol)'''
                
                    # try rotating CO1 vector by OCO angle
                    R_OCO = numpy.matrix([\
                            [math.cos(ang_OCO_r), -math.sin(ang_OCO_r), 0.0],\
                            [math.sin(ang_OCO_r), math.cos(ang_OCO_r), 0.0],\
                            [0.0, 0.0, 1.0]])
                                    
                    r_CO2_mat = r_CO1_mat*R_OCO
                    r_CO2 = vec3(r_CO2_mat[0,0], r_CO2_mat[0,1], r_CO2_mat[0,2])
                
                    O2_CCOO = atom(i.x + r_CC.x + r_CO2.x,
                               i.y + r_CC.y + r_CO2.y,
                               i.z + r_CC.z + r_CO2.z , \
                               'O2_CCOO', -1.0, i.mol)
                
                    O2_CCOO.atom = newcount
                                                 
                    new_list.append(O2_CCOO)

                    #-----------

                    # Bond CO1 (double)
                    CO1_bond = bond(4, newcount-2, newcount-1)
                    bond_list.append(CO1_bond)

                    # Bond CO2 (negative)
                    CO2_bond = bond(5, newcount-2, newcount)
                    bond_list.append(CO2_bond)

                    CCO1_angle = angle(4, icount, newcount-2, newcount-1)
                    angle_list.append(CCO1_angle)
                
                    CCO2_angle = angle(5, icount, newcount-2, newcount)
                    angle_list.append(CCO2_angle)
                                                                                     
                    # Dihedral
                    #print neighbour_list[0]
                    CCCO1_dihedral = dihedral(1, neighbour_list[0], icount,
                                        newcount-2, newcount-1)
                    dihedral_list.append(CCCO1_dihedral)
            
                    CCCO2_dihedral = dihedral(2, neighbour_list[0], icount,
                                     newcount-2, newcount)
                    dihedral_list.append(CCCO2_dihedral)

            if group == 'COH':

                CCO_angle_list = []
                
                r_CO = r_func*l_CO


                # check if neighbouring edge carbon has been functionalised

                C_new_coords = np.around(np.array([i.x, i.y, i.z]),6)

                C_new_plus = np.around(np.array([i.x, i.y + ydist, i.z]),6)

                if i.y + ydist > ly:
                    C_new_plus = np.around(np.array([i.x, i.y + ydist - ly, i.z]),6)


                C_new_minus = np.around(np.array([i.x, i.y - ydist, i.z ]),6)
                if i.y - ydist < 0:
                    C_new_minus = np.around(np.array([i.x, i.y - ydist + ly, i.z]),6)


                # Check if neighbouring C has been functionalised
                if any((C_new_plus == x).all() for x in func_list) == True or any((C_new_minus == x).all() for x in func_list) == True:
                    print 'A neighbour has already been functionalised. Skip atom.'

                else:
                    func_list.append(C_new_coords)

                    ## C

                    atom_list[icount-1] = atom(i.x, i.y, i.z, 'C_COH', 0.265, i.mol)
                    atom_list[icount-1].atom = icount

                    ## O
                    newcount += 1
                    O_COH = atom(i.x + r_CO.x, i.y + r_CO.y, i.z + r_CO.z ,\
                                 'O_COH', -0.70, i.mol)
                    O_COH.atom = newcount

                    # Improper
                    newimproper = improper(1, newcount, neighbour_list[0], \
                                       neighbour_list[1], icount) 
                    improper_list.append(newimproper)

                    # Angle
                    CCO_angle = angle(2, neighbour_list[1], icount, newcount)
                    angle_list.append(CCO_angle)

                    # Bond
                    CO_bond = bond(3, icount, newcount)
                    bond_list.append(CO_bond)

                    new_list.append(O_COH)
                 
                    ## H
                    newcount += 1
                    sign = 1.0
                    #rand = random.uniform(0.0,1.0)
                    #  if rand < 0.5: sign *= -1 # Random sign
                    x_OH = l_OH*math.cos(pi-ang_COH)
                    y_OH = l_OH*math.sin(pi-ang_COH)*sign
                    z_OH = 0.0
                    r_OH = vec3(x_OH, y_OH, z_OH)

                    #determine angle of CO bond with x axis
                    n_x = vec3(1, 0, 0)

                    ang_COZ = math.acos((r_CO*n_x)/abs(r_CO))
                
                    if r_CO.x > 0 + tol and r_CO.y > 0 + tol:
                        ang_COZ = 2*pi - ang_COZ

                    if r_CO.x < 0 - tol and r_CO.y > 0 + tol:
                        ang_COZ = 2*pi - ang_COZ

   
                    ## z-axis rotation matrix
                    R_z = numpy.matrix([\
                              [math.cos(ang_COZ), -math.sin(ang_COZ), 0.0],\
                              [math.sin(ang_COZ), math.cos(ang_COZ), 0.0],\
                              [0.0, 0.0, 1.0]\
                              ])
                    r_OH_mat = numpy.matrix([[r_OH.x, r_OH.y, r_OH.z]])

                    r_OH_mat = r_OH_mat*R_z

                    r_OH = vec3(r_OH_mat[0,0], r_OH_mat[0,1], r_OH_mat[0,2])
  
                    H_COH = atom(i.x + r_CO.x + r_OH.x,
                             i.y + r_CO.y + r_OH.y, 
                             i.z + r_CO.z + r_OH.z , \
                             'H_COH', 0.435, i.mol)
                    H_COH.atom = newcount

                    new_list.append(H_COH)
                

                    OH_bond = bond(4, newcount-1, newcount)
                    bond_list.append(OH_bond)
                
                    COH_angle = angle(3, icount, newcount-1, newcount)
                    angle_list.append(COH_angle)

                    # Dihedral
                    #print neighbour_list[0]
                    CCOH_dihedral = dihedral(1, neighbour_list[0], icount,
                                         newcount-1, newcount)
                    dihedral_list.append(CCOH_dihedral)




            


            if group == 'O_test':
               
                ## Replace Carbon with oxygen
                atom_list[icount-1] = atom(i.x, i.y, i.z, 'O_COH', -0.64, i.mol)
                
                #atom_list = sorted(atom_list, key = lambda atom:atom.atom)

    bond_types.extend(group_bonds)
    angle_types.extend(group_angles)
    improper_types.extend(group_improper)
    dihedral_types.extend(group_dihedral)
    type_list.extend(group_types)
    atom_list.extend(new_list)                                   
    return atom_list, bond_list, angle_list, improper_list, dihedral_list, \
           mol_count, bond_types, angle_types, improper_types, \
           dihedral_types, type_list 

def functionalise_edge2(atom_list, bond_list, angle_list, improper_list, \
                       mol_count, bond_types, angle_types, improper_types, 
                       type_list, a, p, lx, ly, pbc, group):
    """
    Identifies edge atoms using nearest neighbour search and uses vectors
    between neighbours to set up positions for functional groups

    **Test with C-H groups** 
    """
  
    tol = 0.01
    pi = 4.0*math.atan(1.0)

    if group == 'COH': ## Test with COH groups
        ## Define bond type 3 and 4 angle type 2 and 3
        ## Values obtained from Mooney et al., Chemical Physics Letters,
        ## 294 (1998) 135-142
        ## and Konatham et al, Langmuir 2013, Supporting Information (angles)
        group_bonds = [2,3,4]
        group_angles = [2,3,4]
        group_improper = [1]
        group_types = ['C_COH', 'O_COH', 'H_COH', 'C_CH', 'H_CH']
        l_CO = 1.364 
        l_OH = 0.960 
        ang_CCO = 120
#        ang_COH = 180 
        ang_COH = 113*pi/180

    if group == 'CH': ## Test with hydrogenated groups
        ## Define bond type 2, angle type 2, improper type 1
        nbond_types = len(bond_types)                
        group_bonds = [nbond_types+1,nbond_types+2,nbond_types+3]
        group_angles = [2,3,4]                  
        group_improper = [1]
        group_types = ['C_COH', 'O_COH', 'H_COH','C_CH', 'H_CH']
        l_CH = 1.08
        ang_CCH = 120
    
    if group == 'CH_only': ## Test with hydrogenated groups
        ## Define bond type 2, angle type 2, improper type 1
        nbond_types = len(bond_types)
        group_bonds = [nbond_types+1]
        group_angles = [2]
        group_improper = [1]
        group_types = ['C_CH', 'H_CH']
        l_CH = 1.08
        ang_CCH = 120




    a2 = (a+tol)*(a+tol)
    icount = 0
    newcount = len(atom_list)
    new_list = []
    for i in atom_list:
        rand = random.uniform(0.0, 1.0)
        icount += 1
        ni = i.atom
        neighbour_vectors = []
        neighbour_list = []
        jcount = 0
        for j in atom_list:
            jcount += 1
            nj = j.atom
            if icount != jcount:
                x_i = i.x 
                y_i = i.y
                x_j = j.x
                y_j = j.y 
                
                x_ij = x_i - x_j 
                y_ij = y_i - y_j

                if pbc == 'y': 
                    if x_ij < -lx/2: x_ij = x_i - x_j + lx              
                    if x_ij > lx/2: x_ij = x_i - x_j - lx              
                    if y_ij < -ly/2: y_ij = y_i - y_j + ly              
                    if y_ij > ly/2: y_ij = y_i - y_j - ly              
           #    print icount, jcount, len(atom_list)   
                r_ij = vec3(x_ij, y_ij, i.z-j.z)                 
                 
                if r_ij**2 < a2:
                    neighbour_vectors.append(r_ij)
                    neighbour_list.append(j) 

        if len(neighbour_vectors) == 2 and rand <= p: #if less than 3 nearest neighbours
            ## Define vector for functional group
            r_func = vec3(0.0,0.0,0.0)
            for r in neighbour_vectors:
                r_func += r    
            r_func /= math.sqrt(r_func**2)    

            if group == 'COH':

                CCO_angle_list = []
                
                r_CO = r_func*l_CO

                ## C
                atom_list[icount-1] = atom(i.x, i.y, i.z, 'C_COH', 0.2, i.mol)
                atom_list[icount-1].atom = ni
                ## O
                newcount += 1
                O_COH = atom(i.x + r_CO.x, i.y + r_CO.y, i.z + r_CO.z ,\
                             'O_COH', -0.64, i.mol)
                O_COH.atom = newcount

                # Improper 
                newimproper = improper(1, newcount, neighbour_list[0].atom, \
                                       neighbour_list[1].atom, icount) 
                improper_list.append(newimproper)

                # Angle
                CCO_angle = angle(3, neighbour_list[1].atom, icount, newcount)
                angle_list.append(CCO_angle)

                # Bond
                CO_bond = bond(3, ni, newcount)
                bond_list.append(CO_bond)

                new_list.append(O_COH)
                 
                ## H                       
                newcount += 1
                sign = 1.0
                rand = random.uniform(0.0,1.0)
                if rand < 0.5: sign *= -1
                x_OH = l_OH*math.cos(pi-ang_COH)
                y_OH = 0.0
                z_OH = l_OH*math.sin(pi-ang_COH)*sign
                r_OH = vec3(x_OH, y_OH, z_OH)

                #determine angle of CO bond with x axis
                n_x = vec3(1, 0, 0)

                ang_COZ = math.acos((r_CO*n_x)/abs(r_CO))
                
                if r_CO.x > 0 + tol and r_CO.y > 0 + tol:
                    ang_COZ = 2*pi - ang_COZ     

                if r_CO.x < 0 - tol and r_CO.y > 0 + tol:
                    ang_COZ = 2*pi - ang_COZ     

                

                ## z-axis rotation matrix
                R_z = numpy.matrix([\
                              [math.cos(ang_COZ), -math.sin(ang_COZ), 0.0],\
                              [math.sin(ang_COZ), math.cos(ang_COZ), 0.0],\
                              [0.0, 0.0, 1.0]\
                              ])
                r_OH_mat = numpy.matrix([[r_OH.x, r_OH.y, r_OH.z]])

                r_OH_mat = r_OH_mat*R_z

                r_OH = vec3(r_OH_mat[0,0], r_OH_mat[0,1], r_OH_mat[0,2])

                H_COH = atom(i.x + r_CO.x + r_OH.x, 
                             i.y + r_CO.y + r_OH.y, 
                             i.z + r_CO.z + r_OH.z , \
                             'H_COH', 0.44, i.mol)
                H_COH.atom = newcount
                new_list.append(H_COH)
                

                OH_bond = bond(4, newcount-1, newcount)
                bond_list.append(OH_bond)
                
                COH_angle = angle(4, icount, newcount-1, newcount)   
                angle_list.append(COH_angle)

            if group == 'CH':
                newcount += 1
                ## Test with CH group
       
                ## Replace Carbon
                atom_list[icount-1] = atom(i.x, i.y, i.z, 'C_CH', -0.115, i.mol)
                atom_list[icount-1].atom = ni
                r_func *= 1.08
                ## add hydrogen
                H_CH = atom(r_func.x + i.x, r_func.y + i.y, r_func.z + i.z,
                     'H_CH', 0.115, i.mol) 
                new_list.append(H_CH)
                H_CH.atom = newcount
                # Bond
                bond_list.append(bond(nbond_types+1, ni, newcount))
                # Angle
                angle_list.append(angle(2, newcount, 
                                  ni, neighbour_list[0].atom))
                # Improper 
                newimproper = improper(1, newcount, neighbour_list[0].atom, 
                                       neighbour_list[1].atom, ni)
                improper_list.append(newimproper)

            if group == 'CH_only':
                newcount += 1
                ## Test with CH group
    
                ## Replace Carbon
                atom_list[icount-1] = atom(i.x, i.y, i.z, 'C_CH', -0.115, i.mol)
                atom_list[icount-1].atom = ni
                r_func *= 1.08
                ## add hydrogen
                H_CH = atom(r_func.x + i.x, r_func.y + i.y, r_func.z + i.z,
                        'H_CH', 0.115, i.mol)
                new_list.append(H_CH)
                H_CH.atom = newcount
                # Bond
                bond_list.append(bond(nbond_types+1, ni, newcount))
                # Angle
                angle_list.append(angle(2, newcount,
                                        ni, neighbour_list[0].atom))
                # Improper
                newimproper = improper(1, newcount, neighbour_list[0].atom,
                                        neighbour_list[1].atom, ni)
                improper_list.append(newimproper)

    improper_types.extend(group_improper)
    bond_types.extend(group_bonds)
    angle_types.extend(group_angles)
    type_list.extend(group_types)
    atom_list.extend(new_list)
    
    return atom_list, bond_list, angle_list, improper_list, mol_count, \
           bond_types, angle_types, improper_types, type_list 

def functionalise_edge(atom_list, bond_list, angle_list, \
                       mol_count, bond_types, angle_types, type_list, \
                       a, lx, ly, pbc, group):
    """
    Identifies edge atoms using nearest neighbour search and uses vectors
    between neighbours to set up positions for functional groups

    **Test with C-H groups** 
    """
  
    tol = 0.01


    if group == 'CH': ## Test with hydrogenated groups
        ## Define bond type 2, angle type 2                
        group_bonds = [2]
        group_angles = [2]
        group_types = ['C_CH', 'H_CH']
        l_CH = 1.08
    


    a2 = (a+tol)*(a+tol)
    icount = 0
    newcount = len(atom_list)
    new_list = []
    for i in atom_list:
        icount += 1
        neighbour_vectors = []
        neighbour_list = []
        jcount = 0
        for j in atom_list:
            jcount += 1
            if icount != jcount:
                x_i = i.x 
                y_i = i.y
                x_j = j.x
                y_j = j.y 
                
                x_ij = x_i - x_j 
                y_ij = y_i - y_j

                if pbc == 'y': 
                    if x_ij < -lx/2: x_ij = x_i - x_j + lx              
                    if x_ij > lx/2: x_ij = x_i - x_j - lx              
                    if y_ij < -ly/2: y_ij = y_i - y_j + ly              
                    if y_ij > ly/2: y_ij = y_i - y_j - ly              
           #    print icount, jcount, len(atom_list)   
                r_ij = vec3(x_ij, y_ij, i.z-j.z)                 
                 
                if r_ij**2 < a2:
            #       print r_ij**2
                    neighbour_vectors.append(r_ij)
                    neighbour_list.append(jcount) 
  #     print icount, neighbour_vectors
        if len(neighbour_vectors) < 3: #if less than 3 nearest neighbours
            ## Define vector for functional group
            r_func = vec3(0.0,0.0,0.0)
            for r in neighbour_vectors:
                r_func += r    
            r_func /= math.sqrt(r_func**2)    

            if group == 'CH':
                newcount += 1
                ## Test with CH group
       
                ## Replace Carbon
                atom_list[icount-1] = atom(i.x, i.y, i.z, 'C_CH', -0.115, i.mol)
                r_func *= 1.08
                ## add hydrogen
                H_CH = atom(r_func.x + i.x, r_func.y + i.y, r_func.z + i.z, \
                     'H_CH', 0.115, i.mol) 
                new_list.append(H_CH)
                bond_list.append(bond(2, icount, newcount))
                newangle = angle(2, neighbour_list[0], icount, newcount)
                angle_list.append(newangle)

    angle_types.extend(group_angles)
    bond_types.extend(group_bonds)
    type_list.extend(group_types)
    atom_list.extend(new_list)                                   
    return atom_list, bond_list, angle_list, mol_count, \
           bond_types, angle_types, type_list 


def write_xyz(atom_list, title, header):
    "Writes configuration in .xyz format"
    xyz_out = open(title+".xyz", "w")
    n = len(atom_list)
    xyz_out.write("%d\n" % n)
    xyz_out.write(header+"\n")

    for c in atom_list:
        xyz_out.write("%s\t%f\t%f\t%f\n" % (c.spec, c.x, c.y, c.z))

    xyz_out.close()

def write_lammps_molecule(mol, mol_bond, mol_angle, mol_improper, mol_dihedral,
                          type_index, mass_index, fname, header):
    """
    23/5/14
    Creates lammps "molecule" file for a given molecular data
    """
    
    mol_out = open(fname+".mol", "w")
    mol_out.write(header+"\n\n")
    
    mol_out.write("%d atoms\n"%len(mol))
    if mol_bond: mol_out.write("%d bonds\n"%len(mol_bond))
    if mol_angle: mol_out.write("%d angles\n"%len(mol_angle))
    if mol_improper: mol_out.write("%d impropers\n"%len(mol_improper))
    if mol_dihedral: mol_out.write("%d dihedrals\n"%len(mol_dihedral))
    mol_out.write("\n\n")
         
    # Coords
    mol_out.write("Coords\n")
    i = 1
    for a in mol:
        mol_out.write("{0:<6d} {1:<.10e} {2:<.10e} {3:<.10e}\n"
                      .format(i, a.x, a.y, a.z))
        i += 1 
    mol_out.write("\n\n")

    # Types
    mol_out.write("Types\n")
    i = 1
    for a in mol:
        mol_out.write("{0:<6d} {1:<3d}\n"
                      .format(i, type_index[a.spec]))
        i += 1 
    mol_out.write("\n\n")
   
    # Charges 
    mol_out.write("Charges\n")
    i = 1
    for a in mol:
        mol_out.write("{0:<6d} {1:<12.6f}\n"
                      .format(i, a.q))
        i += 1 
    mol_out.write("\n\n")

    # Masses 
    mol_out.write("Masses\n")
    i = 1
    for a in mol:
        mol_out.write("{0:<6d} {1:<12.6f}\n"
                      .format(i, mass_index[type_index[a.spec]]))
        i += 1 
    mol_out.write("\n\n")

    # Bonds
    if mol_bond:
       mol_out.write("Bonds\n")
       i = 1 
       for b in mol_bond:
           mol_out.write("{0:<6d} {1:<6d} {2:<6d} {3:<6d}\n"
                         .format(i, b.bond_type, b.atom1, b.atom2))
           i += 1 
       mol_out.write("\n\n")

    # Angles
    if mol_angle:
       mol_out.write("Angles\n")
       i = 1 
       for c in mol_angle:
           mol_out.write("{0:<6d} {1:<6d} {2:<6d} {3:<6d} {4:<6d}\n"
                         .format(i, c.angle_type, c.atom1, c.atom2, c.atom3))
           i += 1 
       mol_out.write("\n\n")

    # Impropers
    if mol_improper:
       mol_out.write("Impropers\n")
       i = 1 
       for d in mol_improper:
           mol_out.write("{0:<6d} {1:<6d} {2:<6d} {3:<6d} {4:<6d} {5:<6d}\n"
                         .format(i, d.improper_type, d.atomi, d.atomj, 
                                                     d.atomk, d.atoml))
           i += 1 
       mol_out.write("\n\n")

    # Dihedrals
    if mol_dihedral:
       mol_out.write("Dihedrals\n")
       i = 1 
       for e in mol_dihedral:
           mol_out.write("{0:<6d} {1:<6d} {2:<6d} {3:<6d} {4:<6d} {5:<6d}\n"
                         .format(i, e.dihedral_type, e.atomi, e.atomj, 
                                                     e.atomk, e.atoml))
           i += 1 
       mol_out.write("\n\n")




def write_lammps_data3(atom_list, bond_list, angle_list, improper_list,
                      dihedral_list, mass_index, type_index, type_list, \
                      bond_types, angle_types, improper_types, dihedral_types,\
                      lx, ly, lz, charge_flg, title, header):
    "Writes configuration in LAMMPS file format"

    lammps_out = open(title+".dat", "w")

    ## Remove duplicate entries from type lists
    type_list = list(set(type_list))  
    bond_types = list(set(bond_types))  
    angle_types = list(set(angle_types))  
    improper_types = list(set(improper_types))  
    dihedral_types = list(set(dihedral_types))
    

    lammps_out.write(header+"\n\n")
    lammps_out.write("%d atoms\n"%len(atom_list))
    if bond_list: lammps_out.write("%d bonds\n"%len(bond_list))
    if angle_list: lammps_out.write("%d angles\n"%len(angle_list))
    if improper_list: lammps_out.write("%d impropers\n"%len(improper_list))
    if dihedral_list: lammps_out.write("%d dihedrals\n"%len(dihedral_list))
    lammps_out.write("0.0 %f xlo xhi\n" % lx)
    lammps_out.write("0.0 %f ylo yhi\n" % ly)
    lammps_out.write("0.0 %f zlo zhi\n" % lz)
    lammps_out.write("%d atom types\n" % len(type_list))

    if bond_list:
        lammps_out.write("%d bond types\n" % len(bond_types))
    if angle_list:
        lammps_out.write("%d angle types\n" % len(angle_types))
    if improper_list:
        lammps_out.write("%d improper types\n" % len(improper_types))
    if dihedral_list:
        lammps_out.write("%d dihedral types\n" % len(dihedral_types))

    # Masses
    lammps_out.write("\n Masses\n\n") 
    for t in type_list:       
        lammps_out.write("%d %f\n" % (type_index[t], mass_index[type_index[t]]))

    # Atoms
    lammps_out.write("\n Atoms\n\n")
    count = 1
    for a in atom_list:
        lammps_out.write("%d %d %d " % (a.atom, a.mol, type_index[a.spec]))
        if charge_flg:
            lammps_out.write("%f %f %f %f\n" % (a.q, a.x, a.y, a.z)) 
        else:
            lammps_out.write("%f %f %f\n" % (a.x, a.y, a.z))
        count += 1

    # Bonds
    if bond_list:
        count = 1
        lammps_out.write("\n Bonds\n\n")
        for b in bond_list:
            lammps_out.write("%d %d %d %d\n" % \
                             (count, b.bond_type, b.atom1, b.atom2))
            count += 1

    # Angles
    if angle_list:
        count = 1
        lammps_out.write("\n Angles\n\n")
        for c in angle_list:
            lammps_out.write("%d %d %d %d %d\n" % \
                             (count, c.angle_type, c.atom1, c.atom2, c.atom3))
            count += 1
    
    # Impropers
    if improper_list:
        count = 1
        lammps_out.write("\n Impropers\n\n")
        for d in improper_list:
             lammps_out.write("%d %d %d %d %d %d\n" % \
                             (count, d.improper_type, d.atomi, d.atomj, \
                                                      d.atomk, d.atoml))  
             count += 1   

    # Dihedrals
    if dihedral_list:
        count = 1
        lammps_out.write("\n Dihedrals\n\n")
        for e in dihedral_list:
             lammps_out.write("%d %d %d %d %d %d\n" % \
                             (count, e.dihedral_type, e.atomi, e.atomj, \
                                                      e.atomk, e.atoml))  
             count += 1   
    lammps_out.close()

def write_lammps_data2(atom_list, bond_list, angle_list, improper_list,\
                      mass_index, type_index, type_list, \
                      bond_types, angle_types, improper_types,\
                      lx, ly, lz, charge_flg, title, header):
    "Writes configuration in LAMMPS file format"

    lammps_out = open(title+".dat", "w")



    lammps_out.write(header+"\n\n")
    lammps_out.write("%d atoms\n"%len(atom_list))
    if bond_list: lammps_out.write("%d bonds\n"%len(bond_list))
    if angle_list: lammps_out.write("%d angles\n"%len(angle_list))
    if improper_list: lammps_out.write("%d impropers\n"%len(improper_list))
    lammps_out.write("0.0 %f xlo xhi\n" % lx)
    lammps_out.write("0.0 %f ylo yhi\n" % ly)
    lammps_out.write("0.0 %f zlo zhi\n" % lz)
    lammps_out.write("%d atom types\n" % len(type_list))

    if bond_list:
        lammps_out.write("%d bond types\n" % len(bond_types))
    if angle_list:
        lammps_out.write("%d angle types\n" % len(angle_types))
    if improper_list:
        lammps_out.write("%d improper types\n" % len(improper_types))

    # Masses
    lammps_out.write("\n Masses\n\n") 
    for t in type_list:       
        lammps_out.write("%d %f\n" % (type_index[t], mass_index[type_index[t]]))

    # Atoms
    lammps_out.write("\n Atoms\n\n")
    count = 1
    for a in atom_list:
        lammps_out.write("%d %d %d " % (a.atom, a.mol, type_index[a.spec]))
        if charge_flg:
            lammps_out.write("%f %f %f %f\n" % (a.q, a.x, a.y, a.z)) 
        else:
            lammps_out.write("%f %f %f\n" % (a.x, a.y, a.z))
        count += 1

    # Bonds
    if bond_list:
        count = 1
        lammps_out.write("\n Bonds\n\n")
        for b in bond_list:
            lammps_out.write("%d %d %d %d\n" % \
                             (count, b.bond_type, b.atom1, b.atom2))
            count += 1

    # Angles
    if angle_list:
        count = 1
        lammps_out.write("\n Angles\n\n")
        for c in angle_list:
             lammps_out.write("%d %d %d %d %d\n" % \
                             (count, c.angle_type, c.atom1, c.atom2, c.atom3))  
             count += 1   
    
    # Impropers
    if improper_list:
        count = 1
        lammps_out.write("\n Impropers\n\n")
        for d in improper_list:
             lammps_out.write("%d %d %d %d %d %d\n" % \
                             (count, d.improper_type, d.atomi, d.atomj, \
                                                      d.atomk, d.atoml))  
             count += 1   
    lammps_out.close()

def write_lammps_data(atom_list, bond_list, angle_list, \
                      mass_index, type_index, type_list, \
                      bond_types, angle_types, \
                      lx, ly, lz, charge_flg, title, header):
    "Writes configuration in LAMMPS file format"

    lammps_out = open(title+".dat", "w")

    lammps_out.write(header+"\n\n")
    lammps_out.write("%d atoms\n"%len(atom_list))
    lammps_out.write("%d bonds\n"%len(bond_list))
    lammps_out.write("%d angles\n"%len(angle_list))
    lammps_out.write("0.0 %f xlo xhi\n" % lx)
    lammps_out.write("0.0 %f ylo yhi\n" % ly)
    lammps_out.write("0.0 %f zlo zhi\n" % lz)
    lammps_out.write("%d atom types\n" % len(type_list))


    if bond_list:
        lammps_out.write("%d bond types\n" % len(bond_types))
    if angle_list:
        lammps_out.write("%d angle types\n" % len(angle_types))

    # Masses
    lammps_out.write("\n Masses\n\n") 
    for t in type_list:       
        lammps_out.write("%d %f\n" % (type_index[t], mass_index[type_index[t]]))

    # Atoms
    lammps_out.write("\n Atoms\n\n")
    count = 1
    for a in atom_list:
        lammps_out.write("%d %d %d " % (count, a.mol, type_index[a.spec]))
        if charge_flg:
            lammps_out.write("%f %f %f %f\n" % (a.q, a.x, a.y, a.z)) 
        else:
            lammps_out.write("%f %f %f\n" % (a.x, a.y, a.z))
        count += 1

    # Bonds
    if bond_list:
        count = 1
        lammps_out.write("\n Bonds\n\n")
        for b in bond_list: 
            lammps_out.write("%d %d %d %d\n" % \
                             (count, b.bond_type, b.atom1, b.atom2))
            count += 1

    # Angles
    if angle_list:
        count = 1
        lammps_out.write("\n Angles\n\n")
        for c in angle_list:
             lammps_out.write("%d %d %d %d %d\n" % \
                             (count, c.angle_type, c.atom1, c.atom2, c.atom3))  
             count += 1   
    
    lammps_out.close()

def insert_test(atom_list, x, y, z):
    atom_list.append(atom(x,y,z,'O',0.0,0))
    return atom_list

def shift_z(atom_list, lz):
    
    for a in atom_list:
        a.z += lz/2
        if a.z > lz:
            a.z -= lz
                  
    return atom_list
