# -*- coding: utf-8 -*-
"""
Created on Thu May 21 14:02:54 2015

@author: geography
"""

#  Normalises active layer GSD to get rid of mud fraction

import os
import numpy as np

path = str(os.path.join(os.pardir, "Output"))

def normaliseGSD(path):

    mudlist = []
    
    file = os.path.join(path, "Out_ActiveLayer.GSD.F[0]")
    file = open(file, 'r')
    file = file.readlines()
    del file[0:3]
    
    for line in file:
        node = line.split()
        node = [1-float(i) for i in node]
        mudlist.append(node[1:])
        
    n = 1
    while n < 9:
        file = os.path.join(path, "Out_ActiveLayer.GSD.F[" + str(n) + "]")
        file = open(file, 'r')
        file = file.readlines()
        str(file)
        lines = file[3:]
        
        header0 = file[0]
        header1 = file[1]
        header2 = file[2]
        
        outputfile = os.path.join(path, "Out_NormalisedActiveLayer.GSD.F[" + str(n) + "]")
        outputfile = open(outputfile, 'w')
        outputfile.write(str(header0).strip('[').strip(']').strip("'"))
        outputfile.write(str(header1).strip('[').strip(']').strip("'"))
        outputfile.write(str(header2).strip('[').strip(']').strip("'"))
        
        j = 0
        for line in lines:
            splitset = line.split()
            node = splitset[0]
            mud = np.array(mudlist[j])
            bedsize = [float(i) for i in splitset[1:]]
            bedsize = np.array(bedsize)
            fraction = bedsize/mud
            fraction = fraction.tolist()
            outputfile.write(str(node) + '\t' + str(fraction).strip('[').strip(']').replace(', ','\t') + '\n')        
            j = j + 1
            
        outputfile.close()                
        n = n + 1
        
normaliseGSD(path)        