import Atom
import Molecule
import OPLS
import System
import sys
import DA_Polymer
import glob
import numpy as np
import os

# This script is specific to the DA project, it takes in lammps data files
# and generates a polymer system

def main():
    Script, Mult = sys.argv
    Density = 0.1
    Mult = int(Mult)
    File_List = glob.glob('data.*')
    print File_List
    Mol_Temp_List = []
    i = 0
    Total_Mass = 0.0
    MW = 0.0
    for File in File_List:
        Mol_Temp_List.append(DA_Polymer.DA_Polymer(File))
        Total_Mass += Mol_Temp_List[i].MW
        MW = Mol_Temp_List[i].MW
        i += 1
    Avogadro = 6.0221413e23
    Comp_List = np.ones(i+1, dtype=int)
    Comp_List = Comp_List*Mult
    Total_Mols = float(i*Mult)/Avogadro
    Molar_Volume = MW/Density
    Volume = Molar_Volume*Total_Mols
    Box_Length_CM = Volume**(1./3.)
    Box_Length = Box_Length_CM*100000000
    #Name = File_List[0].split('.')[1].split('_')[0] + "_%s" % i*Mult
    
    Name = os.getcwd().split('/')[-1] + "_%s" % int(i*Mult)
    print Name
    DA_System = System.System(Mol_Temp_List, Comp_List, Box_Length, Name)
    DA_System.Gen_Rand()
    DA_System.Write_LAMMPS_Data()
    DA_System.Current_Restart = 'restart.%s_%d_%d' % (DA_System.Name, 420, 25)
    DA_System.Temperature = 420
    DA_System.Run_Lammps_NPT(Temp_Out = 420.0)
    DA_System.Run_Lammps_NPT()
    DA_System.Run_Lammps_NPT()
    DA_System.Run_Lammps_NPT()
    DA_System.Run_Lammps_NPT()
    DA_System.Run_Lammps_NPT()
    DA_System.Run_Lammps_NPT()
    DA_System.Run_Lammps_NPT()
    #System.Uniaxial_Strain(DA_System, Restart_In, steps = 50)


if __name__=='__main__': main()


