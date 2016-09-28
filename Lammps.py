#!usr/bin/python
import numpy as np
import Configure
import subprocess
import os
from lammps import lammps
import System

def Run_Sim_Anneal( File_List, Name):
    
    # Make new directory
    cmd = "mkdir " + Configure.Comet_Path % Name
    subprocess.call(["ssh", Configure.Comet_Login, cmd])

    for File in File_List:
        Data_File = File
        Sim_Name = Data_File.split('.')[1] + "_Fold"
        Data_Out = Data_File + "_Folded"
        In_File = "in." + Sim_Name
        
        Rand = np.random.randint(10000, high = 100000)
        In_Temp = Configure.Template_Path + "in.quench.txt"
        Path = Configure.Comet_Path % Name
        
        # Prepare input script
        with open(In_Temp) as f:
            template = f.read()
            s = template.format(Data_In= Data_File, rand = Rand, Data_Out = Data_Out)
        with open(In_File,'w') as f:
            f.write(s)
        
        # Prepare submission script
        sub_temp = Configure.Template_Path + "sub_Fold"
        NProcs = 8
        Nodes=1
        submit = "sub_%s" % Sim_Name
        with open(sub_temp) as f:
            template = f.read()
            s = template.format(Sim_Name = Sim_Name, path = Path, NProcs = Nodes*NProcs, Nodes=Nodes, tpn = NProcs)
        with open(submit,'w') as f:
            f.write(s)

        # Copy over to Comet
        os.system( Configure.c2c % (submit, Name))
        os.system( Configure.c2c % (In_File, Name))
        os.system( Configure.c2c % (Data_File, Name))
        subprocess.call(["ssh", Configure.Comet_Login, Configure.SBATCH % (Name, submit)])
    return


def Run_Dihedral_Scan( Molecule, File_List):
    Energy = []

    for file in File_List:
        file = open(file, 'r')
        file = file.readlines()
        i = 0
        for line in file:
            try:
                line = line.split()
                coords = np.array([float(line[1]), float(line[2]), float(line[3])])
                Molecule.Atom_List[i].Position = coords
                i += 1
            except:
                continue
        D_Sys = System.System([Molecule], [1], 100.0, "Dihedral")
        D_Sys.Gen_Rand()
        D_Sys.Write_LAMMPS_Data( Dihedral=True)
        os.system("cp %sin.Dihedral_Energy ./" % Configure.Template_Path)
        lmp = lammps()
        lmp.file("in.Dihedral_Energy")
        Output_File = open("log.lammps", 'r')
        Output_File = Output_File.readlines()
        for Line in Output_File:
            try:
                if len(Line.split()) == 2 and float(Line.split()[0]) == float(Line.split()[1]):
                    Energy.append(float(Line.split()[0]))
            except:
                continue

    Energy = np.asarray(Energy)

    return Energy- Energy[0]




