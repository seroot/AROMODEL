#! usr/bin/python

# Open relevant modules

import Bond
import numpy
import time

class Ring(object):
    """
    Class defining a ring object containing 5 or 6 atom objects.
    Instance Variables
    Atom_List = List of atom objects
    ring_type = integer
    Fused = Boolean
    Ring_Name = string
    """
    
   
    def __init__(self, Atom_List):
        # Initialize from list of atom objects
        self.Atom_List = sorted(Atom_List)
        for Atom in self.Atom_List:
            Atom.Ring = True
            print Atom.Element 
        ring_type = len(Atom_List)
        self.Fused = False
        self.Ring_Name = ""
        self.Classify_Ring()
        return
    
    def Classify_Ring(self):
        """
        Function for classifying ring objects according to their alphabetically sorted list of atoms
        """
        Element_List = []
        for Atom in self.Atom_List:
            Element_List.append(Atom.Element)
    
        print "Atom List = ", Element_List
        
        # iterate through commonly occuring ring structures
        
        if Element_List == ['C', 'C', 'C', 'C', 'S']:
            self.Ring_Name = "Thiophene"
        elif Element_List == ['C', 'C', 'C', 'C', 'O']:
            self.Ring_Name = "Furan"
        else: 
            print "This ring type is not currently defined, add it to Ring.py"
        
        print "This ring is a", self.Ring_Name
        return
    
    def Find_OPLS_ID(self):
        """Function for finding the OPLS ID's for atoms embedded within ring objects, this is done seperately and before
        Atom.Find_OPLS_ID() because OPLS ID's depend not only on the connectivity but also, which type ring the atom is in.
        """
    
        return
    
    
    
    
    
    
    
    
    
    

def create_rings(atoms):
    """ Creates Ring object through the use of following each atom's bonds

        Keyword Arguments:
        atoms - The list of atom objects to find the rings
    """
    ring = []
    ringlist = []
    for i in range(len(atoms)):
        layer1 = atoms[i]
        if layer1.ring:
            continue
        path = []
        path.append(layer1)
        for j in range(len(layer1.Bond_List)):
            layer2 = layer1.Bond_List[j]
            if layer2 in path:
                continue
            path.append(layer2)
            for k in range(len(layer2.Bond_List)):
                layer3 = layer2.Bond_List[k]
                if layer3 in path:
                    continue
                path.append(layer3)
                for l in range(len(layer3.Bond_List)):
                    layer4 = layer3.Bond_List[l]
                    if layer4 in path:
                        continue
                    path.append(layer4)
                    for m in range(len(layer4.Bond_List)):
                        layer5 = layer4.Bond_List[m]
                        if layer5 in path:
                            continue
                        for n in range(len(layer5.Bond_List)):
                            layer6 = layer5.Bond_List[n]
                            if layer6 == path[0]:
                                ringlist.append(Ring([layer1,layer2,layer3,layer4,layer5]))
                                layer1.ring = True
                                layer2.ring = True
                                layer3.ring = True
                                layer4.ring = True
                                layer5.ring = True
                            if layer6 in path:
                                continue
                            path.append(layer6)
                            for o in range(len(layer6.Bond_List)):
                                layer7 = layer6.Bond_List[o]
                                if layer7 == path[0]:
                                    ringlist.append(Ring([layer1,layer2,layer3,layer4,layer5,layer6]))
                                    layer1.ring = True
                                    layer2.ring = True
                                    layer3.ring = True
                                    layer4.ring = True
                                    layer5.ring = True
                                    layer6.ring = True
    return ringlist




        
        
        
        
        
        
        
        
        
        
        
        
        
        
          
          