#! usr/bin/python

import Molecule
import random
import numpy as np
from copy import deepcopy
import os
import subprocess
import time
import Configure
import Parallel
import pickle
# Class defining an MD system for simulation with LAMMPS


class System(object):
    """
    Class defining an MD system for simulation with LAMMPS
    instance variables:
        Molecule_List = List of Molecule objects
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
        Num_Mol = 0
        for Comp in self.Composition_List:
            Num_Mol += Comp
        print "Number of Molecules = ", Num_Mol
        self.Molecule_List = []
        self.Current_Restart = ""
        self.Temperature = 800
        

        return

    def Gen_Rand(self):
        i = 0
        k = 1
        for Molecule_Obj in self.Moltemp_List:
            print "Genenerated", self.Composition_List[i], Molecule_Obj.Name, "Molecules"
            for j in range(self.Composition_List[i]):
                Temp_Mol = deepcopy(Molecule_Obj)
                #Temp_Mol = Molecule_Obj
                """
                Temp_Mol.COM = np.asarray([random.random()*self.Box_Size, random.random()*self.Box_Size, random.random()*self.Box_Size], dtype = float)
                Temp_Mol.Mol_ID = k
                k += 1
                Temp_Mol.Adjust_COM()
                self.Molecule_List.append(Temp_Mol)
                """
                
                Deposited = False
                while not Deposited:
                    
                    Temp_Mol.COM = np.asarray([random.random()*self.Box_Size, random.random()*self.Box_Size, random.random()*self.Box_Size], dtype = float)
                    Temp_Mol.Mol_ID = k
                    
                    
                    
                    if len(self.Molecule_List) == 0:
                        Temp_Mol.Adjust_COM()
                        print "No Overlap"
                        self.Molecule_List.append(Temp_Mol)
                        Deposited = True
                        k+= 1
                    
                
                    e = 0
                    for Mol_Obj in self.Molecule_List:
                        Distance = np.linalg.norm(Temp_Mol.COM - Mol_Obj.COM)
                        if Distance > 30:
                            e += 1

                        else:
                            print "Overlap"
                            print Distance
                            break

        
                    if e == len(self.Molecule_List):
                        Temp_Mol.Adjust_COM()
                        self.Molecule_List.append(Temp_Mol)
                        print "Depositing", j
                        Deposited = True
                        k += 1
                    else:
                        print "Didn't Deposit"
                
        
            i += 1
        
    
        Q = 0.0
        Num_Atoms = 0
        for Mol_Obj in self.Molecule_List:
            for Atom_Obj in Mol_Obj.Atom_List:
                Q += Atom_Obj.Charge
                if Atom_Obj.Charge != 0.0:
                    Num_Atoms += 1.0
        print "The total charge of the system is ", Q
        
        dQ = abs(Q/Num_Atoms)
        print "dQ =", dQ

        for Mol_Obj in self.Molecule_List:
            for Atom_Obj in Mol_Obj.Atom_List:
                if Atom_Obj.Charge < 0.0:
                    Atom_Obj.Charge -= dQ
                if Atom_Obj.Charge > 0.0:
                    Atom_Obj.Charge -= dQ
    
        for Molecule_Obj in self.Molecule_List:
            print Molecule_Obj.Mol_ID, Molecule_Obj.Name, Molecule_Obj.COM
        return

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

    def Run_Lammps_Init(self, Nodes = 1):
        cmd = "mkdir " + Configure.Comet_Path % self.Name
        subprocess.call(["ssh", Configure.Comet_Login, cmd])
        
        # Set up input file
        In_Temp = Configure.Template_Path + "in.init_temp"
        In_File = "in.init_%s" % self.Name
        Sim_Name = "init_%s" % self.Name
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
        sub_temp = Configure.Template_Path + "sub_Lammps"
        NProcs = 24
        submit = "sub_%s" % self.Name
        with open(sub_temp) as f:
            template = f.read()
            s = template.format(Sim_Name = Sim_Name, path = Configure.Comet_Path % self.Name, NProcs = Nodes*NProcs, Nodes=Nodes, tpn = NProcs)
        with open(submit,'w') as f:
            f.write(s)
        
        

        File_Out1 = 'log.%s' % Sim_Name
        File_Out = 'restart.%s_800_1' % self.Name

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
                print "Sleeping process", i, "minutes"
                time.sleep(600)
                i += 10
        os.system( 'rm %s' % File_Out)
        self.Current_Restart = File_Out
        return

    def Run_Lammps_NPT(self, GPU = False, Temp_Out = 0.0, time_steps = 1000000, Nodes = 1):
        """
            Function for running NPT dynamics for the system in lammps
            
        """
        Temp_In = self.Temperature
        count = int(self.Current_Restart.split('_')[-1])
        count += 1
        NPT_Temp = Configure.Template_Path + "in.NPT_Temp"
        NPT_In = "in.NPT_%s_%d_%d" % (self.Name, count, Temp_In)
        Sim_Name = NPT_In.split('.')[-1]
        if Temp_Out != 0.0:
            NPT_In += "_%d" % Temp_Out
            Sim_Name = NPT_In.split('.')[-1]
            self.Temperature = Temp_Out
        if Temp_Out == 0.0:
            Temp_Out = Temp_In
        New_Restart = 'restart.%s_%d_%d' % (self.Name, Temp_Out, count)

        with open(NPT_Temp) as f:
            template = f.read()
        s = template.format(Name = self.Name, count = count, Temp_In = Temp_In, Temp_Out = Temp_Out, Restart_In = self.Current_Restart, Restart_Out = New_Restart, Steps = time_steps)
        with open(NPT_In,'w') as f:
            f.write(s)
            #f.insert(33, 'fix def1 all print 100 "${temperature} ${volume} ${dens} ${Enthalpy}" append Thermo_{Temp_In}_{Temp_Out}.txt screen no')
        
        
        sub_temp = Configure.Template_Path + "sub_Lammps"
        NProcs = 24
        submit = "sub_%s" % self.Name
        with open(sub_temp) as f:
            template = f.read()
            s = template.format(Sim_Name = Sim_Name, path = Configure.Comet_Path % self.Name, NProcs = Nodes*NProcs, Nodes=Nodes, tpn = NProcs)
        with open(submit,'w') as f:
            f.write(s)

        File_Out1 = 'log.%s' % Sim_Name
        File_Out = New_Restart
        Traj_File = self.Name + "_%d.lammpstrj" % count

        # Copy over to Comet
        os.system( Configure.c2c % (submit, self.Name))
        os.system( Configure.c2c % (NPT_In, self.Name))
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
                print "Sleeping process", i, "minutes"
                time.sleep(600)
                i += 10

        os.system( 'rm %s' % File_Out)
        self.Current_Restart = File_Out
        os.system( Configure.c2l % (self.Name,  Traj_File))
        return


    def Run_Lammps_Strain(self, Restart_In, Restart_Out, Index, Nodes = 1):
        " Function for running uniaxial straining simulations on the system in LAMMPS"

        Strain_Temp = Configure.Template_Path + "in.strain_Temp"
        Strain_In = "in.%s_strain_%s" % (self.Name, Index)
        Sim_Name = Strain_In.split('.')[-1]
        
        # Prepare input script
        with open(Strain_Temp) as f:
            Template = f.read()
        s = Template.format( Name= self.Name, index = Index, Restart_In= Restart_In, Restart_Out = Restart_Out)
        with open(Strain_In, 'w') as f:
            f.write(s)

        # Prepare Submission Script
        sub_temp = Configure.Template_Path + "sub_Lammps"
        NProcs = 24
        submit = "sub_%s_strain" % self.Name
        with open(sub_temp) as f:
            template = f.read()
        s = template.format(Sim_Name = Sim_Name, path = Configure.Comet_Path % self.Name, NProcs = Nodes*NProcs, Nodes=Nodes, tpn = NProcs)
        with open(submit,'w') as f:
            f.write(s)

        File_Out1 = 'log.%s' % Sim_Name
        File_Out = Restart_Out
        Traj_File1 = "Strain_Equil_%s_%d.lammpstrj" % (self.Name, Index)
        Traj_File2 = "Strain_%s_%d.lammpstrj" % (self.Name, Index)
        Strain_File1 = "%s_Strain_%d.txt" % (self.Name, Index)

        # Copy over to Comet
        os.system( Configure.c2c % (submit, self.Name))
        os.system( Configure.c2c % (Strain_In, self.Name))
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
                print "Sleeping process", i, "minutes"
                time.sleep(600)
                i += 10

        os.system( 'rm %s' % File_Out)
        self.Current_Restart = File_Out
        os.system( Configure.c2l % (self.Name,  Traj_File1))
        os.system( Configure.c2l % (self.Name,  Traj_File2))
        os.system( Configure.c2l % (self.Name,  Strain_File1))
        
        
        
        return

def Run_Glass_Transition(system, Interval, Ramp_Steps = 50000, Equil_Steps = 50000, T_End = 100, Nodes = 1):
    T_init = system.Temperature
    Range = T_init - T_End
    Steps = Range/Interval
    T_out = T_init
    for i in range(Steps):
        T_in = T_out
        T_out = T_out - Interval
        system.Run_Lammps_NPT(Temp_Out= T_out, time_steps = Ramp_Steps,Nodes=Nodes)
        system.Run_Lammps_NPT( time_steps = Equil_Steps, Nodes=Nodes)
        File_Out = "Thermo_%d_%d" % (T_out, T_out) + ".txt"
        os.system( Configure.c2l % (system.Name, File_Out))
    return


def Uniaxial_Strain(system, Restart_In, steps = 20):
    Restart_In_Temp = Restart_In
    for i in range(steps):
        index = i+1
        Restart_Out = "restart." + system.Name + "_Strain_%d" % (i+1)
        system.Run_Lammps_Strain(Restart_In_Temp, Restart_Out, index )
        Restart_In_Temp = Restart_Out

    return







