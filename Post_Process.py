#!/usr/bin/python

import System
import Atom
import numpy as np

class Trajectory(object):
    """
        Class defining a complete trajectory outputted from an MD simulation
        instance variables: Num_Snaps, Snap_List
        """
    
    def __init__(self, File_Name):
        File = open(File_Name,'r')
        File_Lines = File.readlines()
        Temp_Snap_List = []
        for i in range(len(File_Lines)):
            line = File_Lines[i]
            if line == "ITEM: TIMESTEP\n":
                Time_Step = int(File_Lines[i+1])
                N = int(File_Lines[i+3])
                Atom_List = np.empty(N, dtype=object)
                XBounds = [ float(File_Lines[i+5].split()[0]), float(File_Lines[i+5].split()[1])]
                YBounds = [ float(File_Lines[i+6].split()[0]), float(File_Lines[i+6].split()[1])]
                ZBounds = [ float(File_Lines[i+7].split()[0]), float(File_Lines[i+7].split()[1])]
                Box_Dim = [ abs(XBounds[1]- XBounds[0]), abs(YBounds[1] - YBounds[0]), abs(ZBounds[1]- ZBounds[0])]
                
                # Define the unstrained Box Length
                if (Box_Dim[0] == Box_Dim[1]) and (Box_Dim[1] == Box_Dim[2]):
                    self.L0 = Box_Dim[1]
                    print "Setting L0 to ", self.L0
                
                print "Timestep = %d, N = %d\n" % (Time_Step, N)
                print Box_Dim
                # Extract Snapshot as list of atoms
                print "Extracting Atoms...\n"
                for j in range(N):
                    Atom_Line = File_Lines[i + 9 + j].split(' ')
                    ID = int(Atom_Line[0])
                    TYPE = int(Atom_Line[1])
                    MOL = int(Atom_Line[2])
                    POS= [ float(Atom_Line[3]), float(Atom_Line[4]), float(Atom_Line[5])]
                    IMAGE = [ int(Atom_Line[6]), int(Atom_Line[7]), int(Atom_Line[8])]


