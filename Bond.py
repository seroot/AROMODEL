#! usr/bin/python


# Import relevant modules


# Bond Class

class Bond(object):
    """
    Class defining a bond between atom objects
    instance variables:
        Type = int
        Bond_Master = Atom object
        Bond_Slave = Atom object
        Bond_Length = float
        Bond_Vector = np[3, float]
        req = float
        kb = float
    """
    def __init__(self, Bond_Master, Bond_Slave, req):
        self.Bond_Master = Bond_Master # Atom object
        self.Bond_Slave = Bond_Slave # Atom object
        self.req = req # distance
        self.kb = 0.0
        self.Bond_ID = 0
        self.Lammps_Type = 0
        self.System_ID = 0
        return




