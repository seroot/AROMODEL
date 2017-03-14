#!usr/bin/python

# Import relevant modules
import numpy as np
import os
import subprocess
import shutil
import time
import random
import glob

# Import classes
import Molecule
import Configure

def Run_Dihedral_Scan(Molecule, dihedral_id, dih_init, dih_end, local = False):
    # Function for running a dihedral scan using Orca on Comet
    File_Name = self.Name + "_Dih_%d.inp" % dihedral_id
    File = open(File_Name, 'w')
    File.write('! RKS B3LYP 6-31+G** NormalSCF Opt NOSOSCF CHELPG PAL8\n\n')
    File.write('%scf MaxIter 500 end\n')
    File.write('%geom Scan\n')
    return
