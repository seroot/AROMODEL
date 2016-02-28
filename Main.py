import Atom
import Molecule
import OPLS
import System
import sys
import DA_Polymer

def main():
    Script, File_Name, N = sys.argv
    N = int(N)
    print File_Name
    Solvent = Molecule.Molecule(File_Name)
    Solvent.Set_Up_FF(run_orca=True)
    OPLS.Assign_OPLS(Solvent, ChelpG = False)
    Solvent_System = System.System([Solvent], [N], 50.0)
    Solvent_System.Gen_Rand()
    Solvent_System.Write_LAMMPS_Data()
    Solvent_System.Run_Lammps_Init()
    

if __name__=='__main__': main()


