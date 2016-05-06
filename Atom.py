#! usr/bin/python


# Import relevant modules
import numpy as np
# Atom Class

Mass_Dict = {"C":12.011, "H":1.008, "Cl":35.453, "Si":28.086, "O":15.9994, "S":32.06, "F":18.998, "N":14.007}

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
        #self.COM_Position = np.zeros(3,dtype=float)
        self.Image_Flags = np.zeros(3,dtype=int)
        self.Unwrapped_Position = np.zeros(3,dtype=float)
        return


# Functions operating on sets of Atom objects


def Find_OPLS_ID(Atom, Fullerene):
    # Funcition for assigning OPLS Types and Classes based on connectivity of the atom
    if Atom.Element == "C":
        Temp_Bond_List = sorted([ Atomobj.Element for Atomobj in Atom.Bond_List ])
        if len(Temp_Bond_List) == 3:
            if Temp_Bond_List == ['C', 'C', 'H']:
                # Aromatic Carbon
                Atom.OPLS_Type = 90
                Atom.OPLS_Class = 48
            elif Temp_Bond_List == ['C', 'C', 'Cl']:
                # Chlorobenzene >C-Cl
                Atom.OPLS_Type = 205
                Atom.OPLS_Class = 48
            elif Temp_Bond_List == ['C', 'C', 'C']:
                if not Fullerene:
                    Fullerene = raw_input("Is this a fullerene? True or False: ")
                if Fullerene:
                    Atom.OPLS_Type = 907
                    Atom.OPLS_Class = 48
                # Aromatic Carbon
                else:
                    Atom.OPLS_Type = 90
                    Atom.OPLS_Class = 48
            elif Temp_Bond_List == ['C', 'C', 'Si']:
                Atom.OPLS_Type = 907
                Atom.OPLS_Class = 48
            elif Temp_Bond_List == ['C', 'O', 'O']:
                # Ester -COOR"
                Atom.OPLS_Type = 406
                Atom.OPLS_Class = 3
            elif Temp_Bond_List == ['C','C','S']:
                Atom.OPLS_Type =909
                Atom.OPLS_Class = 48
            elif Temp_Bond_List == ['C','H','S']:
                Atom.OPLS_Type = 911
                Atom.OPLS_Class =48
            elif Temp_Bond_List == ['C','N','O']:
                Atom.OPLS_Type = 853
                Atom.OPLS_Class = 3
            elif Temp_Bond_List == ['C','C', 'N']:
                # Quinoline C2
                Atom.OPLS_Type = 545
                Atom.OPLS_Class = 48
            elif Temp_Bond_List == ['C','C','O']:
                # Furan C2
                Atom.OPLS_Type = 508
                Atom.OPLS_Class = 84
            elif Temp_Bond_List == ['C', 'H', 'O']:
                # Furan C3
                Atom.OPLS_Type = 509
                Atom.OPLS_Class = 87
    
        elif len(Temp_Bond_List) == 4:
            if Temp_Bond_List == ['C', 'C', 'H', 'H']:
                # Alkane C-CH2-C
                Atom.OPLS_Type = 81
                Atom.OPLS_Class = 13
            elif Temp_Bond_List == ['C', 'H', 'H', 'H']:
                # Alkane C-CH3
                Atom.OPLS_Type = 80
                Atom.OPLS_Class = 13
            elif Temp_Bond_List == ['C', 'C', 'C', 'H']:
                # Alkane C-CH3
                Atom.OPLS_Type = 82
                Atom.OPLS_Class = 13
            elif Temp_Bond_List == ['C', 'C', 'H', 'Si']:
                Atom.OPLS_Type = 873
                Atom.OPLS_Class = 13
            elif Temp_Bond_List == ['Cl', 'Cl', 'Cl', 'H']:
                # CH-Cl3
                Atom.OPLS_Type = 48
                Atom.OPLS_Class = 13
            elif Temp_Bond_List == ['H', 'H', 'H', 'O']:
                # Methyl Ester CH3-O-R
                Atom.OPLS_Type = 409
                Atom.OPLS_Class = 13
            elif Temp_Bond_List == ['C', 'C', 'C', 'C']:
                # Alkcane >C<
                Atom.OPLS_Type = 84
                Atom.OPLS_Class = 13
            elif Temp_Bond_List == ['C', 'C', 'H', 'Si']:
                # Alkyl silane
                Atom.OPLS_Type = 873
                Atom.OPLS_Class = 13
            elif Temp_Bond_List == ['H','H','H', 'N']:
                # 1-methylimidazole
                Atom.OPLS_Type = 603
                Atom.OPLS_Class = 13
            elif Temp_Bond_List == ['H','H','H', 'H']:
                # 1-methylimidazole
                Atom.OPLS_Type = 83
                Atom.OPLS_Class = 13
                    
                    

    elif Atom.Element == "H":
        Bonded_Atom = Atom.Bond_List[0]
        if Bonded_Atom.OPLS_Type == 90:
            # Aromatic Hydrogen C-H
            Atom.OPLS_Type = 91
            Atom.OPLS_Class = 49
        elif Bonded_Atom.OPLS_Class == 13:
            Atom.OPLS_Type = 85
            Atom.OPLS_Class = 46
        elif Bonded_Atom.OPLS_Class == 48:
            Atom.OPLS_Type = 91
            Atom.OPLS_Class = 49
        else:
            Atom.OPLS_Type = 85
            Atom.OPLS_Class = 46
    elif Atom.Element == "Cl":
        Bonded_Atom = Atom.Bond_List[0]
        if Bonded_Atom.OPLS_Class == 48:
            # ChloroBenzene
            Atom.OPLS_Type = 206
            Atom.OPLS_Class = 21
        if Bonded_Atom.OPLS_Class == 13:
            # CH-Cl3
            Atom.OPLS_Type = 49
            Atom.OPLS_Class = 21

    elif Atom.Element == "Si":
        Temp_Bond_List = sorted([ Atomobj.Element for Atomobj in Atom.Bond_List])
        if Temp_Bond_List == ['C', 'C', 'C', 'C']:
            # Alkyl Silane
            Atom.OPLS_Type = 866
            Atom.OPLS_Class = 108

    elif Atom.Element == "O":
        Temp_Bond_List = sorted([ Atomobj.Element for Atomobj in Atom.Bond_List])
        if len(Temp_Bond_List) == 2:
            if Temp_Bond_List == ['C', 'C']:
                # Dialkyl Ether
                Atom.OPLS_Type = 408
                Atom.OPLS_Class = 20
        elif len(Temp_Bond_List) ==  1:
            if Temp_Bond_List == ['C']:
                #   Ester C=O
                Atom.OPLS_Type = 407
                Atom.OPLS_Class = 4
    elif Atom.Element == "S":
        Temp_Bond_List = sorted([ Atomobj.Element for Atomobj in Atom.Bond_List])
        if len(Temp_Bond_List) == 2:
            if Temp_Bond_List == ['C', 'C']:
                #Sulfide Thiophene
                Atom.OPLS_Type = 908
                Atom.OPLS_Class = 16
    elif Atom.Element == "N":
        Temp_Bond_List = sorted([ Atomobj.Element for Atomobj in Atom.Bond_List])
        if len(Temp_Bond_List) == 3:
            if Temp_Bond_List == ['C','C','C']:
                # DPP
                Atom.OPLS_Type = 847
                Atom.OPLS_Class = 107



    return Fullerene






