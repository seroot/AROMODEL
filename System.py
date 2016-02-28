#! usr/bin/python

import Molecule
import random
from copy import deepcopy
import os
import subprocess
import time
import Configure
# Class defining an MD system for simulation with LAMMPS


class System(object):
    """
    Class defining an MD system for simulation with LAMMPS
    instance variables:
        Molecule_List = List of Molecule objects ( d
        Composition_List = List of integers defining the number of each molecule type in the system
        Box_Size = float (Angstrom)
        Atom_Params = [[Mass, Sigma, Epsilon]] List of floats
        Bond_Params = [[KB, RO]] List of floats
        Angle_Params = [[KA, A0]] List of floats
        Dihedral_Params = [[V1,V2,V3,V4]] List of floats
        Improper_Params = [[KI, AI]] List of floats
        Num_Atoms = float
        
    """

    def __init__(self, Moltemp_List, Composition_List, Box_Size, Name):
        self.Name = Name
        
    
        self.Moltemp_List = Moltemp_List
        self.Atom_Params, self.Bond_Params, self.Angle_Params, self.Dihedral_Params, self.Improper_Params = Molecule.Assign_Lammps(Moltemp_List)
        self.Composition_List = Composition_List
        self.Box_Size = Box_Size
        self.Molecule_List = []
        

        return

    def Gen_Rand(self):
        i = 0
        for Molecule_Obj in self.Moltemp_List:
            print "Genenerated", self.Composition_List[i], Molecule_Obj.Name, "Molecules"
            for j in range(self.Composition_List[i]):
                Temp_Mol = deepcopy(Molecule_Obj)
                Temp_Mol.COM += [random.random()*self.Box_Size, random.random()*self.Box_Size, random.random()*self.Box_Size]
                Temp_Mol.Mol_ID = j + i + 1
                Temp_Mol.Adjust_COM()
                self.Molecule_List.append(Temp_Mol)
            i += 1
    
        Q = 0
        Num_Atoms = 0
        for Mol_Obj in self.Molecule_List:
            for Atom_Obj in Mol_Obj.Atom_List:
                Q += Atom_Obj.Charge
                Num_Atoms += 1.0
        print "The total charge of the system is ", Q
        
        
        dQ = Q/Num_Atoms
        print "dQ =", dQ

        for Mol_Obj in self.Molecule_List:
            for Atom_Obj in Mol_Obj.Atom_List:
                if Atom_Obj.Charge < 0.0:
                    Atom_Obj.Charge += dQ
                if Atom_Obj.Charge > 0.0:
                    Atom_Obj.Charge -= dQ
    
        for Molecule_Obj in self.Molecule_List:
            print Molecule_Obj.Mol_ID, Molecule_Obj.Name, Molecule_Obj.COM


    def Write_LAMMPS_Data(self):
        """
        Function for writing LAMMPS data file 
        """
        self.Data_File = self.Name + ".data"
        File = open(self.Data_File, 'w')
        File.write('LAMMPS data file via System.Write_LAMMPS_Data()\n\n')
        
        # Find number of atoms, bonds, dihedrals, impropers in the system
        self.Num_Atoms = 0
        self.Num_Bonds = 0
        self.Num_Angles = 0
        self.Num_Dihedrals = 0
        self.Num_Impropers = 0
        i = 0
        for Moltemp_Obj in self.Moltemp_List:
            self.Num_Atoms += len(Moltemp_Obj.Atom_List)*self.Composition_List[i]
            self.Num_Bonds += len(Moltemp_Obj.Bond_List)*self.Composition_List[i]
            self.Num_Angles += len(Moltemp_Obj.Angle_List)*self.Composition_List[i]
            self.Num_Dihedrals += len(Moltemp_Obj.Dihedral_List)*self.Composition_List[i]
            self.Num_Impropers += len(Moltemp_Obj.Improper_List)*self.Composition_List[i]
            i += 1
        
        
        File.write('%d atoms\n' % self.Num_Atoms)
        File.write('%d atom types\n' % len(self.Atom_Params))
        File.write('%d bonds\n' % self.Num_Bonds)
        File.write('%d bond types\n' % len(self.Bond_Params))
        File.write('%d angles\n' % self.Num_Angles)
        File.write('%d angle types\n' % len(self.Angle_Params))
        if self.Num_Dihedrals > 0:
            File.write('%d dihedrals\n' % self.Num_Dihedrals)
            File.write('%d dihedral types\n' % len(self.Dihedral_Params))
    
        if self.Num_Impropers > 0:
            File.write('%d impropers\n' % self.Num_Impropers)
            File.write('%d improper types\n' % len(self.Improper_Params))

        File.write('\n\n0.0000 %.4f xlo xhi\n' % self.Box_Size)
        File.write('0.0000 %.4f ylo yhi\n' % self.Box_Size)
        File.write('0.0000 %.4f zlo zhi\n' % self.Box_Size)

        File.write('\n\nMasses\n\n')
        i = 1
        for Params in self.Atom_Params:
            File.write('%d %.3f\n' % ( i, Params[0]))
            i += 1

        File.write('\n\nPair Coeffs # lj/cut/coul/long\n\n')
        i = 1
        for Params in self.Atom_Params:
            File.write('%d %.3f %.3f\n' % (i, Params[2], Params[1]))
            i += 1

        File.write('\n\nBond Coeffs # harmonic\n\n')
        i = 1
        for Params in self.Bond_Params:
            File.write('%d %.4f %.4f\n' % (i, Params[0], Params[1])) # Its possible that im missing a factor of 2 here
            i += 1

        File.write('\n\nAngle Coeffs # harmonic\n\n')
        i = 1
        for Params in self.Angle_Params:
            File.write('%d %.4f %.4f\n' % (i, Params[0], Params[1])) # Its possible that im missing a factor of 2 here
            i += 1
        if self.Num_Dihedrals > 0:
            File.write('\n\nDihedral Coeffs # opls\n\n')
            i = 1
            for Params in self.Dihedral_Params:
                File.write('%d %.4f %.4f %.4f %.4f\n' % (i, Params[0], Params[1], Params[2], Params[3]))
                i += 1

        if self.Num_Impropers > 0:
            File.write('\n\nImproper Coeffs # harmonic\n\n')
            i = 1
            for Params in self.Improper_Params:
                File.write('%d %.4f -1 2\n' % (i, Params[0]))
                i += 1

        File.write('\n\nAtoms # full\n\n')
        i = 1
        j = 1
        for Molecule_Obj in self.Molecule_List:
            for Atom_Obj in Molecule_Obj.Atom_List:
                Atom_Obj.System_ID = i
                File.write('%d %d %d %.8f %.4f %.4f %.4f\n' % (Atom_Obj.System_ID, Molecule_Obj.Mol_ID, Atom_Obj.LAMMPS_Type, Atom_Obj.Charge, Atom_Obj.Position[0], Atom_Obj.Position[1], Atom_Obj.Position[2]))
                i += 1
            j += 1
        

        File.write('\n\nBonds\n\n')
        i = 1
        for Molecule_Obj in self.Molecule_List:
            for Bond_Obj in Molecule_Obj.Bond_List:
                Bond_Obj.System_ID = i
                File.write('%d %d %d %d\n' % ( Bond_Obj.System_ID, Bond_Obj.LAMMPS_Type, Bond_Obj.Bond_Master.System_ID, Bond_Obj.Bond_Slave.System_ID))
                i += 1


        File.write('\n\nAngles\n\n')
        i = 1
        for Molecule_Obj in self.Molecule_List:
            for Angle_Obj in Molecule_Obj.Angle_List:
                Angle_Obj.System_ID = i
                File.write('%d %d %d %d %d\n' % (Angle_Obj.System_ID, Angle_Obj.LAMMPS_Type, Angle_Obj.Angle_Slave1.System_ID, Angle_Obj.Angle_Master.System_ID,  Angle_Obj.Angle_Slave2.System_ID))
                i += 1
        if self.Num_Dihedrals > 0:
            File.write('\n\nDihedrals\n\n')
            i = 1
            for Molecule_Obj in self.Molecule_List:
                for Dihedral_Obj in Molecule_Obj.Dihedral_List:
                    Dihedral_Obj.System_ID = i
                    File.write('%d %d %d %d %d %d\n' % (Dihedral_Obj.System_ID, Dihedral_Obj.LAMMPS_Type, Dihedral_Obj.Dihedral_Slave1.System_ID, Dihedral_Obj.Dihedral_Master1.System_ID, Dihedral_Obj.Dihedral_Master2.System_ID, Dihedral_Obj.Dihedral_Slave2.System_ID))
                    i += 1

        if self.Num_Impropers > 0:
            File.write('\n\nImpropers\n\n')
            i = 1
            for Molecule_Obj in self.Molecule_List:
                for Improper_Obj in Molecule_Obj.Improper_List:
                    Improper_Obj.System_ID = i
                    File.write('%d %d %d %d %d %d\n' % (Improper_Obj.System_ID, Improper_Obj.LAMMPS_Type, Improper_Obj.Improper_Master.System_ID, Improper_Obj.Improper_Slave1.System_ID, Improper_Obj.Improper_Slave2.System_ID, Improper_Obj.Improper_Slave3.System_ID))
                    i += 1

    def Run_Lammps_Init(self, GPU = False):
        cmd = "mkdir" + Configure.Comet_Path % self.Name
        subprocess.call(["ssh", Configure.Comet_Login, cmd])
        
        # Set up input file
        In_Temp = Configure.Template_Path + "in.init_temp"
        In_File = "in.init_%s" % self.Name
        with open(In_Temp) as f:
            template = f.read()
        s = template.format(System_Name = self.Name)
        with open(In_File,'w') as f:
            f.write(s)
        """
        Current Rule of thumb for Parallel Job Submission:
        MPI: 1 processor per 1000 atoms
        MPI + GPU: 4 GPU and 24 Processors --> Only use for >100K particles
        Need to do more benchmarks with current system --> good task for CJ
        """
        
        NProcs = int(self.Num_Atoms/1000.)
        if NProcs > 48:
            GPU = True
        # Set up submit script
        if not GPU:
            if NProcs < 24:
                sub_temp = Configure.Template_Path + "sub_Lammps"
                submit = "sub_%s" % self.Name
                with open(sub_temp) as f:
                    template = f.read()
                    s = template.format(System_Name = self.Name, path = Configure.Comet_Path % self.Name, NProcs = NProcs, Nodes=1, tpn = NProcs)
                with open(submit,'w') as f:
                    f.write(s)
            if NProcs > 24 and NProcs < 48:
                sub_temp = Configure.Template_Path + "sub_Lammps"
                submit = "sub_%s" % self.Name
                tpn = NProcs/2
                NProcs = tpn*2
                with open(sub_temp) as f:
                    template = f.read()
                    s = template.format(System_Name = self.Name, path = Configure.Comet_Path % self.Name, NProcs = NProcs, Nodes=2, tpn = tpn)
                with open(submit,'w') as f:
                    f.write(s)
    
        elif GPU:
            sub_temp = Configure.Template_Path + "GPU_Sub"
            submit = "GPU_sub_%s" % self.Name
            with open(sub_temp) as f:
                template = f.read()
            s = template.format(System_Name = self.Name, path = Configure.Comet_Path % self.Name)
            with open(submit,'w') as f:
                f.write(s)
            

        File_Out1 = 'log.Init_%s' % self.Name
        File_Out = 'restart.Condensed_%s' % self.Name

        # Copy over to Comet
        os.system( Configure.c2c % (submit, self.Name))
        os.system( Configure.c2c % (In_File, self.Name))
        os.system( Configure.c2c % (self.Data_File, self.Name))
        os.system( Configure.c2l % (self.Name, File_Out1))
        try:
            File = open(File_Out1,'r')
        except:
            subprocess.call(["ssh", Configure.Comet_Login, Configure.SBATCH % (self.Name, submit)])
        Finished = False
        i = 0
        while not Finished:
            os.system( Configure.c2l % (self.Name, File_Out))
            try:
                File = open(File_Out,'r')
                Finished = True
            except:
                print "Sleeping process", i, "miniutes"
                time.sleep(600)
                i += 10
        os.system( 'rm %s' % File_Out)
        return





