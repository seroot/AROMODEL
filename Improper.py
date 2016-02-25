#! usr/bin/python

# Improper Class

class Improper(object):
    """
    Class defining an improper interactions between 4 atoms
    This is mainly just for aromatic carbons to enforce planarity
    instance variables
        Type = int
        Improper_Master = Atom object
        Improper_Slave1 = Atom object
        Improper_Slave2 = Atom object
        Improper_Slave3 = Atom object ( This one is the "odd one out")
        Ki = float (Force/degree)
        Improper_Eq = float (degrees)  
    """

    def __init__(self, M, S1, S2, S3, Ki, Improper_Eq, Improper_ID):
        self.Improper_Master = M
        self.Improper_Slave1 = S1
        self.Improper_Slave2 = S2
        self.Improper_Slave3 = S3 # outside the ring
        self.Ki = Ki
        self.Improper_Eq = Improper_Eq
        self.Improper_ID = Improper_ID
        self.LAMMPS_Type = 1
        self.System_ID = 0
        return