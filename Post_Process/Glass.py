#!usr/bin/python


from sys import argv
import matplotlib.pyplot as plt
import glob
import numpy
from scipy.optimize import curve_fit
import pylab

def Extract_Therm( Filename):
    File = open(Filename)
    TEMP = []
    temp = []
    DENS = []
    dens = []
    VOL = []
    vol = []
    ENTH = []
    enth = []
    for line in File:
        line = line.split()
        try:
            if (float(line[0]>0.0)):
                temp.append(float(line[0]))
                dens.append(float(line[2]))
                vol.append(float(line[1]))
                enth.append(float(line[3]))
        except:
            TEMP.append(numpy.mean(numpy.asarray(temp)))
            DENS.append(numpy.mean(numpy.asarray(dens)))
            VOL.append(numpy.mean(numpy.asarray(vol)))
            ENTH.append(numpy.mean(numpy.asarray(enth)))
            temp = []
            dens = []
            vol = []
            enth = []
            continue

    return TEMP, DENS, VOL, ENTH



