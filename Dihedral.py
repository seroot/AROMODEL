#! usr/bin/python

# Open relevant modules

import numpy
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def Fourier_Series(Phi, A1, A2, A3, A4):
    return 0.5*(A1*(1+numpy.cos(Phi)) + A2*(1-numpy.cos(2.0*Phi)) + A3*(1+numpy.cos(3.0*Phi)) + A4*(1-numpy.cos(4.0*Phi)))

def Fourier_Series_Val(Phi, Pop):
    return 0.5*(Pop[0]*(1+numpy.cos(Phi)) + Pop[1]*(1-numpy.cos(2.0*Phi)) + Pop[2]*(1+numpy.cos(3.0*Phi)) + Pop[3]*(1-numpy.cos(4.0*Phi)))

def Multi_Harmonic(Phi, A1, A2, A3, A4, A5):
    return A1 + A2*numpy.cos(Phi) + A3*numpy.cos(Phi)**2 + A4*numpy.cos(Phi)**3 + A5*numpy.cos(Phi)**4

def Multi_Harmonic_Val(Phi, Pop):
    return Pop[0] + Pop[1]*numpy.cos(Phi) + Pop[2]*numpy.cos(Phi)**2 + Pop[3]*numpy.cos(Phi)**3 + Pop[4]*numpy.cos(Phi)**4

class Dihedral(object):
    """
    Class defining a dihedral angle between 4 atom objects
    Instance Variables
        Type = int
        Dihedral_Master1 = Atom object
        Dihedral_Master2 = Atom object
        Dihedral_Slave1 = Atom object
        Dihedral_Slave2= Atom object
        Dihedral_Eq = float
        Coeffs = [list]
    """

    def __init__(self, Dihedral_Master1, Dihedral_Master2, Dihedral_Slave1, Dihedral_Slave2, Dihedral_Eq):
        self.Dihedral_Master1 = Dihedral_Master1 # Atom object
        self.Dihedral_Master2 = Dihedral_Master2 # Atom object
        self.Dihedral_Slave1 = Dihedral_Slave1
        self.Dihedral_Slave2 = Dihedral_Slave2
        self.Dihedral_Eq = Dihedral_Eq # Assigned in OPLS.Assign_OPLS
        self.Coeffs = numpy.zeros(4,dtype=float) # Assigned in OPLS.Assign_OPLS
        self.Dihedral_ID = 0 # Assigned in OPLS.Assign_OPLS
        self.LAMMPS_Type = 0 # Assigned in Molecule.Assign_Lammps
        self.System_ID = 0 # Assigned in System.Write_LAMMPS_Data()
        return
        

    def Fit_Parameters(self, Dihedral_Energy, Dihedral_Angles):
        Pop, Pco = curve_fit(Fourier_Series, Dihedral_Angles, Dihedral_Energy)
        plt.plot(Dihedral_Angles, Dihedral_Energy, 'o', label= "MP2 Derived")
        plt.plot(Dihedral_Angles, Fourier_Series_Val(Dihedral_Angles, Pop), '-', label="Fourier Fit")
        plt.ylabel('Relative Energy (kcal/mol)', fontsize=25)
        plt.xlabel('Dihedral Angle (Radians)', fontsize=25)
        plt.xlim((Dihedral_Angles[0],Dihedral_Angles[-1]))
        plt.axhline(y=0.0, linestyle='--')
        plt.legend(frameon=False)
        plt.savefig('MP2_Dihedral_%d.png' % self.Dihedral_ID)
        plt.close()
        self.Coeffs = Pop
        print self.Coeffs
        return





