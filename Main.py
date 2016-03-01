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
    Name = File_Name.split('.')[0] + "_%s" % N
    Solvent = Molecule.Molecule(File_Name)
    Solvent.UnConverged = True
    Solvent.Set_Up_FF(run_orca=True)
    OPLS.Assign_OPLS(Solvent, ChelpG = False)
    Solvent_System = System.System([Solvent], [N], 100.0, Name)
    Solvent_System.Gen_Rand()
    Solvent_System.Write_LAMMPS_Data()
    Solvent_System.Run_Lammps_Init()
    

if __name__=='__main__': main()


