#! usr/bin/python

# OPLS

import Atom
import Bond
import Angle
import Dihedral
import Improper
import Molecule
import Configure


def Assign_OPLS(Molecule, ChelpG = True):
    # Assign OPLS TYPES and CLASSES
    print "Finding OPLS Types and Classes"
    print "------------------------------"
    print "Element Type Class"
    print "______________________________"
    Fullerene = False
    for Atom_Obj in Molecule.Atom_List:
        Fullerene = Atom.Find_OPLS_ID(Atom_Obj, Fullerene)
        print Atom_Obj.Atom_ID, Atom_Obj.Element, Atom_Obj.OPLS_Type, Atom_Obj.OPLS_Class
    
    # Open up OPLS FILE
    OPLS_Path = Configure.Template_Path + "oplsaa.prm.txt"
    OPLS_FILE = open(OPLS_Path, 'r')
    File_Lines = OPLS_FILE.readlines()

    # Initialize Parameters Lists
    OPLS_Bonds = []
    OPLS_Angles = []
    OPLS_Dihedrals = []
    OPLS_Impropers = []
    OPLS_VDW = []
    OPLS_CHARGE = []
    
    # Fill in parameter list
    for Line in File_Lines:
        Line = Line.split()
        try:
            if Line[0] == "bond":
                OPLS_Bonds.append(Line)
            if Line[0] == "angle":
                OPLS_Angles.append(Line)
            if Line[0] == "torsion":
                OPLS_Dihedrals.append(Line)
            if Line[0] == "imptors":
                OPLS_Impropers.append(Line)
            if Line[0] == "vdw":
                OPLS_VDW.append(Line)
            if Line[0] == "charge":
                OPLS_CHARGE.append(Line)
        except:
            continue


    print "------------------------------------"
    print "Finding VDW parameters"
    print "------------------------------------"
    print "Element, OPLS TYPE, Sigma, Epsilon"
    print "------------------------------------"

    i = 0
    for Atom_Obj in Molecule.Atom_List:
        for VDW in OPLS_VDW:
            if int(VDW[1]) == Atom_Obj.OPLS_Type:
                i += 1
                Atom_Obj.Sigma = float(VDW[2])
                Atom_Obj.Epsilon = float(VDW[3])
                print Atom_Obj.Atom_ID, Atom_Obj.Element, Atom_Obj.OPLS_Type, Atom_Obj.Sigma, Atom_Obj.Epsilon

    if i == len(Molecule.Atom_List):
        print "All atom types present and accounted for :)"
    else:
        print "Missing OPLS Type: Edit Find_OPLS_ID function in Atom.py"

    if ChelpG == False:
        print "--------------------------------"
        print "Finding partial charges"
        print "--------------------------------"
        for Atom_Obj in Molecule.Atom_List:
            for CHARGE in OPLS_CHARGE:
                if int(CHARGE[1]) == Atom_Obj.OPLS_Type:
                    print Atom_Obj.OPLS_Type
                    print "CHelpG:", Atom_Obj.Charge
                    print "OPLS:", float(CHARGE[2])
                    Atom_Obj.Charge = float(CHARGE[2])
        print "Neutralizing Molecule"
        Q = 0.0
        Num_Atoms = 0.0
        for Atom_Obj in Molecule.Atom_List:
            if Atom_Obj.Charge != 0.0:
                Q += Atom_Obj.Charge
                Num_Atoms += 1.0
        
        print "Total Charge on Molecule is", Q
        print "N is", Num_Atoms
        dQ = Q/Num_Atoms
        print "dQ =", dQ
        for Atom_Obj in Molecule.Atom_List:
            if Atom_Obj.Charge != 0.0:
                if Atom_Obj.Charge < 0.0:
                    Atom_Obj.Charge -= dQ
                    print "Adjusting Charge +"
                if Atom_Obj.Charge > 0.0:
                    Atom_Obj.Charge -= dQ
                    print "Adjusting Charge -"
        Q = 0.0
        Num_Atoms = 0.0
        for Atom_Obj in Molecule.Atom_List:
            if Atom_Obj.Charge != 0.0:
                Q += Atom_Obj.Charge
                Num_Atoms += 1.0
                                                
        print "Total Charge on Molecule is", Q


    print "------------------------------------"
    print "Finding Bonds"
    print "------------------------------------"
    print "Bond #, Orca eq, OPLS eq"
    i = 0
    for Bond_Obj in Molecule.Bond_List:
        Bonded_Atoms = sorted([Bond_Obj.Bond_Master.OPLS_Class, Bond_Obj.Bond_Slave.OPLS_Class])
        for B_OPLS in OPLS_Bonds:
            if int(B_OPLS[1]) == Bonded_Atoms[0] and int(B_OPLS[2]) == Bonded_Atoms[1]:
                i += 1
                Bond_Obj.kb = float(B_OPLS[3])
                if Molecule.UnConverged:
                    Bond_Obj.req = float(B_OPLS[4])
                Bond_Obj.Bond_ID = i
                print i, Bond_Obj.req, B_OPLS[4], Bond_Obj.kb

    print "-------------------------------------"
    if i == len(Molecule.Bond_List):
        print "All bonds present and accounted for :)"
    else:
        print "MISSING BOND PARAMETER!!!!!"
        for Bond_Obj in Molecule.Bond_List:
            if Bond_Obj.kb == 0:
                print Bond_Obj.Bond_Master.OPLS_Class, Bond_Obj.Bond_Slave.OPLS_Class
                # Set arbitrary value for bond constant
                Bond_Obj.Kb = 300.0


    print "-------------------------------------"
    print "Finding Angles"
    print "-------------------------------------"
    print "Angle #, Orca eq, OPLS eq"
    i = 0
    for Angle_Obj in Molecule.Angle_List:
        Angle_Slaves = sorted([Angle_Obj.Angle_Slave1.OPLS_Class, Angle_Obj.Angle_Slave2.OPLS_Class])
        for A_OPLS in OPLS_Angles:
            if int(A_OPLS[2]) == Angle_Obj.Angle_Master.OPLS_Class and Angle_Slaves[0] == int(A_OPLS[1]) and Angle_Slaves[1] == int(A_OPLS[3]):
                i += 1
                Angle_Obj.ka = float(A_OPLS[4])
                if Molecule.UnConverged:
                    Angle_Obj.Angle_Eq = float(A_OPLS[5])
                Angle_Obj.Angle_ID = i
                print i, Angle_Obj.Angle_Eq, float(A_OPLS[5])
    print "-------------------------------------"

    if i == len(Molecule.Angle_List):
        print " All angles present and accounted for :)"
    else:
        print "MISSING ANGLE PARAMETER!!!"
        for Angle_Obj in Molecule.Angle_List:
            if Angle_Obj.ka == 0:
                print Angle_Obj.Angle_Master.OPLS_Class, Angle_Obj.Angle_Slave1.OPLS_Class, Angle_Obj.Angle_Slave2.OPLS_Class
                Angle_Obj.ka = 35.0

    print "-------------------------------------"
    print "Finding Dihedrals"
    print "-------------------------------------"

    i = 0
    for Dihedral_Obj in Molecule.Dihedral_List:
        Dihedral_Masters = sorted([Dihedral_Obj.Dihedral_Master1.OPLS_Class, Dihedral_Obj.Dihedral_Master2.OPLS_Class])
        Dihedral_Slaves = sorted([Dihedral_Obj.Dihedral_Slave1.OPLS_Class, Dihedral_Obj.Dihedral_Slave2.OPLS_Class])
        #print Dihedral_Slaves[0], Dihedral_Masters, Dihedral_Slaves[1]
        for D_OPLS in OPLS_Dihedrals:
            M1 = int(D_OPLS[2])
            M2 = int(D_OPLS[3])
            S1 = int(D_OPLS[1])
            S2 = int(D_OPLS[4])
            if Dihedral_Masters[0] == M1 and Dihedral_Masters[1] == M2  and Dihedral_Slaves[0] == S1 and Dihedral_Slaves[1] == S2:
                i = i+1
                Dihedral_Obj.Coeffs[0] = float(D_OPLS[5])
                Dihedral_Obj.Coeffs[1] = float(D_OPLS[8])
                Dihedral_Obj.Coeffs[2] = float(D_OPLS[11])
                Dihedral_Obj.Dihedral_ID = i
                #print "Check"


    if i == len(Molecule.Dihedral_List):
        print "All dihedrals present and accounted for :)"
    else:
        print "MISSING DIHEDRAL PARAMETERS!!!"
        print "i =", i
        print "Num Dih=", len(Molecule.Dihedral_List)
        Molecule.Missing_Dihedrals = len(Molecule.Dihedral_List) - i 
        for Dihedral_Obj in Molecule.Dihedral_List:
            if Dihedral_Obj.Dihedral_ID == 0:
                print Dihedral_Obj.Dihedral_Slave1.OPLS_Class, Dihedral_Obj.Dihedral_Master1.OPLS_Class, Dihedral_Obj.Dihedral_Master2.OPLS_Class, Dihedral_Obj.Dihedral_Slave2.OPLS_Class
                print Dihedral_Obj.Dihedral_Slave1.Atom_ID, Dihedral_Obj.Dihedral_Master1.Atom_ID, Dihedral_Obj.Dihedral_Master2.Atom_ID, Dihedral_Obj.Dihedral_Slave2.Atom_ID
    print "Finding Improper interactions (Aromatic Carbons OPLS_CLASS=48)"
    i = 0
    for Atom_Obj in Molecule.Atom_List:
        if Atom_Obj.OPLS_Class == 48:
            i += 1
            Temp_Bond_List = sorted(Atom_Obj.Bond_List, key=lambda x: x.Element)
            Molecule.Improper_List.append(Improper.Improper( Atom_Obj, Temp_Bond_List[0], Temp_Bond_List[1], Temp_Bond_List[2], 5.0, 180.0, i ))
    print "Found ", i , "Improper interactions"
    print "---------------------------------------"
    print " Force Field Building is complete"
    print "---------------------------------------"


    return