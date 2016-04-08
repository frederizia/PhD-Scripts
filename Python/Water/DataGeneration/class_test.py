# Testing classes in python - a simple 'atom' class
from numpy import matrix

class atom:
    def __init__(self, xcoord, ycoord, zcoord, species, charge, molecule):
        self.x = xcoord
        self.y = ycoord
        self.z = zcoord
        self.spec = species
        self.q = charge
        self.mol = molecule
        self.atom = 0
    def return_r_matrix(self):
        r = matrix([[self.x, self.y, self.z]])
        return r

class bond:
    def __init__(self, Bond_type, first_atom, second_atom):
        self.bond_type = Bond_type
        self.atom1 = first_atom
        self.atom2 = second_atom
        self.mol = 0  

class angle:
    def __init__(self, Angle_type, first_atom, second_atom, third_atom):
        self.angle_type = Angle_type
        self.atom1 = first_atom
        self.atom2 = second_atom   # middle atom
        self.atom3 = third_atom
        self.mol = 0
class improper:
    def __init__(self, Improper_type, first_atom, second_atom, third_atom, fourth_atom):
        self.improper_type = Improper_type
        self.atomi = first_atom
        self.atomj = second_atom   # middle atom
        self.atomk = third_atom
        self.atoml = fourth_atom
        self.mol = 0 
class dihedral:
    def __init__(self, Dihedral_type, first_atom, second_atom, third_atom, fourth_atom):
        self.dihedral_type = Dihedral_type
        self.atomi = first_atom
        self.atomj = second_atom 
        self.atomk = third_atom
        self.atoml = fourth_atom
        self.mol = 0 
