#! usr/bin/python


# Import relevant modules
import numpy as np
import os
import subprocess
import shutil
import time
# import Class structure
import Atom
import Bond
import Angle
import Dihedral
import Configure


class Molecule(object):
    """
    Class defining a molecule
    instance variables (self) :
        N = Number of Atoms
        Name = Molecule Name
        Atom_list[N] = List of Atom objects
        Bond_List[Num_Bonds] = list of bond objects
        Angle_List[Num_Angles] = list of angle objects
        Dihedral_List[Num_Dihedrals] = list of Dihedral objects
        Improper_List[Num_Impropers] = list of Improper objects
        COM = Center of Mass
        MOL_ID = ID Number for molecule
        unConverged = Flag for bypassing Orca convergence (Default = False)
    """

    def __init__(self, File_Name):
        File = open(File_Name,'r') # File_Name is the name of an .xyz file outputted by Avogadro
        File_Lines = File.readlines()
        
        # Define Instance Variables
        self.Name = File_Name.split('.xyz')[0]
        self.N = int(File_Lines[0].strip('\n')) # Integer
        self.Atom_List = np.empty(self.N, dtype=object) # Numpy object array
        self.Bond_List = []
        self.Angle_List = []
        self.Dihedral_List = []
        self.Improper_List = []
        self.MW = 0.0
        self.COM = np.zeros(3,float)
        self.Mol_ID = 0
        self.Missing_Dihedrals = 0
        self.unConverged = False # Unconverged Orca Optimization
        
        print "Setting up molecule"
        print "Molecule Name:", self.Name
        print self.N, "Atoms in ", self.Name
        print "----------------------------------"
        for i in range(self.N):
            Line = File_Lines[2+i]
            Line = Line.strip('\n').split()
            Element = Line[0]
            Position = np.array( [ float(Line[1]), float(Line[2]), float(Line[3]) ], dtype=float )
            self.Atom_List[i] = Atom.Atom(Position, Element, i+1) # Instantiate Atom_List with Atom objects
            self.MW += Atom_List[i].Mass
        
        print "Initial XYZ Coordinates:\n"
        for Atom_Obj in self.Atom_List:
            print Atom_Obj.Element, Atom_Obj.Position
        print "----------------------------------"
        return

    def Adjust_COM(self):
        for Atom_Obj in self.Atom_List:
            Atom_Obj.Position += self.COM
        return

    def Set_Up_FF(self, run_orca=True):
        if run_orca:
            print "Setting up Orca input script"
        
            # Write Orca Input File
            File_Name = self.Name + ".inp"
            File = open(File_Name, 'w')
            File.write('! RKS B3LYP 6-31+G** TightSCF Opt CHELPG\n\n')
            File.write('*xyz 0 1\n')
            for Atom_Obj in self.Atom_List:
                File.write('%s %.5f %.5f %.5f\n' % ( Atom_Obj.Element, Atom_Obj.Position[0], Atom_Obj.Position[1],Atom_Obj.Position[2]))
            File.write('*')
            File.close()
            
            """
            #Run subprocess on local machine
            os.system('mkdir Orca')
            os.system('mv %s ./Orca' % File_Name)
            os.chdir('./Orca')
            File_Out = self.Name + ".out"
            os.system('orca %s > %s' %(File_Name, File_Out)) # Run Orca Job
            """
            
            File_Out = self.Name + ".out"
            print "Running Orca Geometry Optimization on Comet"
            cmd = "mkdir" + Configure.Comet_Path % self.Name
            subprocess.call(["ssh", Configure.Comet_Login, cmd])
            subtemp = Configure.Template_Path + "sub_orca_temp"
            submit = "submit_orca"
            Path = Configure.Comet_Path % self.Name
            # Write submit script
            with open(subtemp) as f:
                template = f.read()
            s = template.format(Comet_Path=Path, Orca_Path = Configure.Orca_Path, name = self.Name )
            with open(submit ,"w") as f:
                f.write(s)
            
            # Copy Files over  to Comet
            os.system( Configure.c2c % (submit, self.Name))
            os.system( Configure.c2c % (File_Name, self.Name))
            # Run job
            os.system( Configure.c2l % (self.Name, File_Out))
            try:
                File = open(File_Out,'r')
            except:
                subprocess.call(["ssh",Configure.Comet_Login, Configure.SBATCH % (self.Name, submit)])

            Finished = False
            i = 0
            while not Finished:
                os.system( Configure.c2l % (self.Name, File_Out))
                try:
                    File = open(File_Out,'r')
                    File_Lines = File.readlines()
                    print File_Lines[-1]
                    if File_Lines[-1].split(' ')[0] == "TOTAL":
                        Finished = True
                    else:
                        print "Not Finished"
                        i += 10
                        print "Sleeping process", i, "Minutes"
                        time.sleep(600)
                except:
                    print "Sleeping process", i, "miniutes"
                    time.sleep(600)
                    i += 10
    
                
            
    
    
        else:
            os.chdir('./Orca')
        
        File_Out = self.Name + ".out"
        
        # Extract info from Orca output file
        Orca_File = open(File_Out, 'r')
        File_Lines = Orca_File.readlines()
        print "Extracting Redundant Coordinates..."
        for i in range(len(File_Lines)):
            Line = File_Lines[i].strip('\n').split()
            try:
                if Line[0] == "Redundant" and (File_Lines[i+2].strip('\n').split()[1] == "Optimized" or self.UnConverged == True):
                    Redundant = True
                    j = i + 7
                    while Redundant:
                        R_Line = File_Lines[j]
                        R_Line = R_Line.split()
                        try:
                            if int(R_Line[0].strip('.')) >= 1:
                                j = j + 1
                                Type = R_Line[1].split('(')[0]
                                # Extract Bonds
                                if Type == "B":
                                    Master_ID = int(R_Line[2].split(',')[0])
                                    Slave_ID = int(R_Line[3].split(')')[0])
                                    req = float(R_Line[-1])
                                    self.Atom_List[Master_ID].Bond_List.append(self.Atom_List[Slave_ID])
                                    self.Atom_List[Slave_ID].Bond_List.append(self.Atom_List[Master_ID])
                                    self.Bond_List.append( Bond.Bond(self.Atom_List[Master_ID], self.Atom_List[Slave_ID], req))
                                # Extract Angles
                                if Type == "A":
                                    Master_ID = int(R_Line[3].split(',')[0])
                                    Slave1_ID = int(R_Line[2].split(',')[0])
                                    Slave2_ID = int(R_Line[4].split(')')[0])
                                    Angle_Eq = float(R_Line[-1])
                                    self.Angle_List.append( Angle.Angle(self.Atom_List[Master_ID], self.Atom_List[Slave1_ID], self.Atom_List[Slave2_ID], Angle_Eq))
                                if Type == "D":
                                    Master1_ID = int(R_Line[3].split(',')[0])
                                    Master2_ID = int(R_Line[4].split(',')[0])
                                    Slave1_ID = int(R_Line[2].split(',')[0])
                                    Slave2_ID = int(R_Line[5].split(')')[0])
                                    Dihedral_Eq = float(R_Line[-1])
                                    self.Dihedral_List.append(Dihedral.Dihedral(self.Atom_List[Master1_ID],self.Atom_List[Master2_ID], self.Atom_List[Slave1_ID], self.Atom_List[Slave2_ID], Dihedral_Eq))
                        except:
                            Redundant = False
    
                if Line[0] == "CHELPG" and len(Line) == 2:
                    for j in range(self.N):
                        Chelp_Line = File_Lines[i+2+j].split()
                        index = int(Chelp_Line[0])
                        charge = float(Chelp_Line[3])
                        self.Atom_List[index].Charge = charge
                
                if Line[0] == "CARTESIAN" and Line[2] == "(ANGSTROEM)":
                    for j in range(self.N):
                        XYZ_Line = File_Lines[i+2+j].split()
                        Position_Temp = np.array( [float(XYZ_Line[1]) ,float(XYZ_Line[2]) ,float(XYZ_Line[3])], dtype = float)
                        self.Atom_List[j].Position = Position_Temp
            except:
                continue
        print "Redundant Internal Coordinates and ChelpG partial charges Extracted"
        
        if not run_orca:
            os.chdir('..')
        
        print "Optimized XYZ Coordinates:\n"
        for Atom_Obj in self.Atom_List:
            # Finds OPLS Types and Classes
            print Atom_Obj.Atom_ID, Atom_Obj.Element, Atom_Obj.Position, sorted([ Atomobj.Element for Atomobj in Atom_Obj.Bond_List ])
        print "----------------------------------"

        return


# Functions Operating on sets of Molecule objects

def Assign_Lammps(Moltemp_List):
    """
        Function that inputs a list of molecule templates. It searches through all the atoms, bonds, angles etc. to find the unique types of interactions
        present in an arbitrary system object made of molecule templates.
        
        returns a list of unique params for writing to a LAMMPS data file. these are lists defined such that the i-1 element corresponds to the ith LAMMPS_Type
    """
    
    print "Finding unique Atoms"
    Unique_Atoms = []
    Atom_Params = []
    for Moltemp_Obj in Moltemp_List:
        for Atom_Obj in Moltemp_Obj.Atom_List:
            i = 1
            for Type in Unique_Atoms:
                if Atom_Obj.OPLS_Type != Type:
                    i += 1
                if Atom_Obj.OPLS_Type == Type:
                    Atom_Obj.LAMMPS_Type = i
            if i > len(Unique_Atoms):
                Unique_Atoms.append(Atom_Obj.OPLS_Type)
                Atom_Params.append([Atom_Obj.Mass, Atom_Obj.Sigma, Atom_Obj.Epsilon])
                Atom_Obj.LAMMPS_Type = i
                print "Found new atom type:", i, "OPLS_ID =", Atom_Obj.OPLS_Type, Atom_Params[i-1]

    print "Finding unique Bonds"
    Bond_Params = []
    for Moltemp_Obj in Moltemp_List:
        for Bond_Obj in Moltemp_Obj.Bond_List:
            i = 1
            for Params in Bond_Params:
                if Params[0] != Bond_Obj.kb or Params[1] != Bond_Obj.req:
                    i += 1
                if Params[0] == Bond_Obj.kb and Params[1] == Bond_Obj.req:
                    Bond_Obj.LAMMPS_Type = i
            if i > len(Bond_Params):
                Bond_Params.append([Bond_Obj.kb, Bond_Obj.req])
                Bond_Obj.LAMMPS_Type = i
                print "Found new bond type:", i, Bond_Obj.kb, Bond_Obj.req

    print "Finding unique Angles"
    Angle_Params = []
    for Moltemp_Obj in Moltemp_List:
        for Angle_Obj in Moltemp_Obj.Angle_List:
            i = 1
            for Params in Angle_Params:
                if Params[0] != Angle_Obj.ka or Params[1] != Angle_Obj.Angle_Eq:
                    i += 1
                if Params[0] == Angle_Obj.ka and Params[1] == Angle_Obj.Angle_Eq:
                    Angle_Obj.LAMMPS_Type = i
            if i > len(Angle_Params):
                Angle_Params.append([Angle_Obj.ka, Angle_Obj.Angle_Eq])
                Angle_Obj.LAMMPS_Type = i
                print "Found new angle type:", i, Angle_Params[i-1]


    print "Finding unique Dihedrals"
    Dihedral_Params = []
    for Moltemp_Obj in Moltemp_List:
        for Dihedral_Obj in Moltemp_Obj.Dihedral_List:
            i = 1
            for Params in Dihedral_Params:
                if Dihedral_Obj.Coeffs[0] != Params[0] or Dihedral_Obj.Coeffs[1] != Params[1] or Dihedral_Obj.Coeffs[2] != Params[2] or Dihedral_Obj.Coeffs[3] != Params[3]:
                    i += 1
                if Dihedral_Obj.Coeffs[0] == Params[0] and Dihedral_Obj.Coeffs[1] == Params[1] and Dihedral_Obj.Coeffs[2] == Params[2] and Dihedral_Obj.Coeffs[3] == Params[3]:
                    Dihedral_Obj.LAMMPS_Type = i
            if i > len(Dihedral_Params):
                Dihedral_Params.append(Dihedral_Obj.Coeffs)
                Dihedral_Obj.LAMMPS_Type = i
                print "Found new dihedral type:", i, Dihedral_Params[i-1]


    print "Finding unique Impropers"
    Improper_Params = []
    for Moltemp_Obj in Moltemp_List:
        for Improper_Obj in Moltemp_Obj.Improper_List:
            i = 1
            for Params in Improper_Params:
                if Improper_Obj.Ki != Params[0] or Improper_Obj.Improper_Eq != Params[1]:
                    i += 1
                if Improper_Obj.Ki == Params[0] and Improper_Obj.Improper_Eq == Params[1]:
                    Improper_Obj.LAMMPS_Type = i
        try:
            if i > len(Improper_Params):
                Improper_Params.append([Improper_Obj.Ki, Improper_Obj.Improper_Eq])
                Improper_Obj.LAMMPS_Type = i
                print "Found improper type:", i, Improper_Obj.Ki, Improper_Obj.Improper_Eq
        except:
            continue



    return Atom_Params, Bond_Params, Angle_Params, Dihedral_Params, Improper_Params










