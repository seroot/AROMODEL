#!/usr/bin/python



import /Users/seroot/Desktop/Research/AROMODEL/Atom
import /Users/seroot/Desktop/Research/AROMODEL/Molecule
import /Users/seroot/Desktop/Research/AROMODEL/OPLS
import /Users/seroot/Desktop/Research/AROMODEL/System
import /Users/seroot/Desktop/Research/AROMODEL/sys
import /Users/seroot/Desktop/Research/AROMODEL/DA_Polymer

def main():
    Script, File_Name, N = sys.argv
    N = int(N)
    print File_Name
    Name = File_Name.split('.')[0] + "_%s" % N
    Solvent = Molecule.Molecule(File_Name)
    Solvent.UnConverged = True
    Solvent.Set_Up_FF(run_orca=True, local = False)
    OPLS.Assign_OPLS(Solvent, ChelpG = False)
    Solvent_System = System.System([Solvent], [N], 200.0, Name)
    #Solvent_System.Gen_Rand()
    Solvent_System.Write_LAMMPS_Data()
    #Solvent_System.Run_Lammps_Init()
    

if __name__=='__main__': main()


