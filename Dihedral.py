#! usr/bin/python

# Open relevant modules

import numpy

class Dihedral(object):
    """
    Class defining an angle between 3 atom objects
    Instance Variables
        Type = int
        Dihedral_Master1 = Atom object
        Dihedral_Master2 = Atom object
        Dihedral_Slave1 = Atom object
        Dihedral_Slave2= Atom object
        Dihedral_Eq = float
        Coeffs = [list]
    """

    def __init__(self, Dihedral_Master1, Dihedral_Master2, Dihedral_Slave1, Dihedral_Slave2, Dihedral_Eq):
        self.Dihedral_Master1 = Dihedral_Master1 # Atom object
        self.Dihedral_Master2 = Dihedral_Master2 # Atom object
        self.Dihedral_Slave1 = Dihedral_Slave1
        self.Dihedral_Slave2 = Dihedral_Slave2
        self.Dihedral_Eq = Dihedral_Eq # Assigned in OPLS.Assign_OPLS
        self.Coeffs = numpy.zeros(4,dtype=float) # Assigned in OPLS.Assign_OPLS
        self.Dihedral_ID = 0 # Assigned in OPLS.Assign_OPLS
        self.LAMMPS_Type = 0 # Assigned in Molecule.Assign_Lammps
        self.System_ID = 0 # Assigned in System.Write_LAMMPS_Data()
        return
        

