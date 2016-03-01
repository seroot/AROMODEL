#! usr/bin/python


# Import relevant modules
import numpy as np
# Atom Class

Mass_Dict = {"C":12.011, "H":1.008, "Cl":35.453, "Si":28.0855, "O":15.9994, "S":32.06, "F":18.998}

class Atom(object):
    """
    Class defining an atom
    instance variables:
        Position = Position Vector numpy[3,float]
        Charge = Float
        Element = String
        Bond_List = list of atom objects
        Atom_ID = 0
    """
    def __init__(self, Position, Element, Atom_ID):
        self.Position = Position # Numpy[3,float]
        self.Element = Element # String
        self.Charge = 0.0 # Float
        self.Bond_List = [] # List of atom objects
        self.OPLS_Type = 0
        self.OPLS_Class = 0
        self.Mass = Mass_Dict[self.Element]
        self.Sigma = 0.0
        self.Epsilon = 0.0
        self.Atom_ID = Atom_ID
        self.LAMMPS_Type = 0
        self.System_ID = 0
        self.COM_Position = np.zeros(3,dtype=float)
        return


# Functions operating on sets of Atom objects

def Find_OPLS_ID(Atom):
    # Funcition for assigning OPLS Types and Classes based on connectivity of the atom
    if Atom.Element == "C":
        Temp_Bond_List = sorted([ Atomobj.Element for Atomobj in Atom.Bond_List ])
        if len(Temp_Bond_List) == 3:
            if Temp_Bond_List == ['C', 'C', 'H']:
                # Aromatic Carbon
                Atom.OPLS_Type = 90
                Atom.OPLS_Class = 48
            if Temp_Bond_List == ['C', 'C', 'Cl']:
                # Chlorobenzene >C-Cl
                Atom.OPLS_Type = 205
                Atom.OPLS_Class = 48
            if Temp_Bond_List == ['C', 'C', 'C']:
                # Aromatic Carbon
                Atom.OPLS_Type = 90
                Atom.OPLS_Class = 48
            if Temp_Bond_List == ['C', 'C', 'Si']:
                Atom.OPLS_Type = 907
                Atom.OPLS_Class = 48
        if len(Temp_Bond_List) == 4:
            if Temp_Bond_List == ['C', 'C', 'H', 'H']:
                # Alkane C-CH2-C
                Atom.OPLS_Type = 81
                Atom.OPLS_Class = 13
            if Temp_Bond_List == ['C', 'H', 'H', 'H']:
                # Alkane C-CH3
                Atom.OPLS_Type = 80
                Atom.OPLS_Class = 13
            if Temp_Bond_List == ['C', 'C', 'C', 'H']:
                # Alkane C-CH3
                Atom.OPLS_Type = 82
                Atom.OPLS_Class = 13
            if Temp_Bond_List == ['C', 'C', 'H', 'Si']:
                Atom.OPLS_Type = 873
                Atom.OPLS_Class = 13
            if Temp_Bond_List == ['Cl', 'Cl', 'Cl', 'H']:
                # CH-Cl3
                Atom.OPLS_Type = 48
                Atom.OPLS_Class = 13
    if Atom.Element == "H":
        Bonded_Atom = Atom.Bond_List[0]
        if Bonded_Atom.OPLS_Type == 90:
            # Aromatic Hydrogen C-H
            Atom.OPLS_Type = 91
            Atom.OPLS_Class = 49
        if Bonded_Atom.OPLS_Class == 13:
            Atom.OPLS_Type = 85
            Atom.OPLS_Class = 46
        else:
            Atom.OPLS_Type = 85
            Atom.OPLS_Class = 46
    if Atom.Element == "Cl":
        Bonded_Atom = Atom.Bond_List[0]
        if Bonded_Atom.OPLS_Class == 48:
            # ChloroBenzene
            Atom.OPLS_Type = 206
            Atom.OPLS_Class = 21
        if Bonded_Atom.OPLS_Class == 13:
            # CH-Cl3
            Atom.OPLS_Type = 49
            Atom.OPLS_Class = 21

    if Atom.Element == "Si":
        Temp_Bond_List = sorted([ Atomobj.Element for Atomobj in Atom.Bond_List])
        if Temp_Bond_List == ['C', 'C', 'C', 'C']:
            # Alkyl Silane
            Atom.OPLS_Type = 866
            Atom.OPLS_Class = 108
    return






