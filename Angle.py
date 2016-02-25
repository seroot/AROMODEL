#! usr/bin/python

# Open relevant modules

class Angle(object):
    """ 
    Class defining an angle between 3-atom objects
    Instance Variables:
        Type = int
        Angle_Master = Atom object
        Angle_Slave1 = Atom object
        Angle_Slave2 = Atom object
        Angle_Normal = np[3,float]
        Angle_Th = float (radians0
        Angle_Eq = float (degrees)
        ka = float
    """
    def __init__(self, Angle_Master, Angle_Slave1, Angle_Slave2, Angle_Eq):
        self.Angle_Master = Angle_Master #Atom Object (Atom in the middle)
        self.Angle_Slave1 = Angle_Slave1 #Atom Object
        self.Angle_Slave2 = Angle_Slave2 #Atom Object
        self.Angle_Eq = Angle_Eq
        self.ka = 0.0 # Defined in OPLS.Assign_OPLS
        self.Angle_ID = 0 # Defined in OPLS.Assign_OPLS
        self.LAMMPS_Type = 0 # Defined in Molecule.Assign_Lammps
        self.System_ID = 0 # Defined System.Write_LAMMPS_Data()
        return
