# -*- coding: utf-8 -*-
"""
Created on Tue May 19 11:35:48 2015

@author: geography

Plots Degradation/Active Layer D50 against migration rate for various floodplain scenarios
"""
import os
import matplotlib.pyplot as plt
from size_spec_load_plot import extractGSDdata
from ActiveLayerGSDnormaliser import normaliseGSD

miglist = [0.0, 1.0, 1.5, 2.0]
#pblist = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
pblist = [0.0, 0.6, 1.0]
#labellist = ['0.0 (3.85 mm)', '0.2 (5.06 mm)','0.4 (7.01 mm)','0.6 (12.04 mm)','0.8 (24.45 mm)','1.0 (45.25 mm)']
labellist = ['0.0 (3.85 mm)','0.6 (12.04 mm)','1.0 (45.25 mm)']
colorlist = ['k','b','r','g','c','m']

bigdeglist = []
bigDlist = []

# Load data
mainpath = os.path.join(os.pardir, "migration_vs_alphabar")
outputfolders = os.listdir(mainpath)

for p in pblist:
    deglist = []
    Dlist = []
    for m in miglist:
        name = str(m)+"_"+str(p)
        for folder in outputfolders:
            if name in folder:
                path = str(os.path.join(mainpath, folder))
                degfile = open(path + '//' + 'Out_CumulativeBedChange', 'r')
                degfile = degfile.readlines()
                str(degfile)
                Dfile = open(path + '//' + 'Out_ActiveLayer.GSD.D50', 'r')
                Dfile = Dfile.readlines()
                str(Dfile)
                #normaliseGSD(path)
                #data = extractGSDdata(path)
                #cobmob = data[3]
                degnode = degfile[3].split()
                Dnode = Dfile[3].split()
                degvalue = float(degnode[-1])
                #Dvalue = float(cobmob[-1])
                Dvalue = float(Dnode[-1])
                deglist.append(degvalue)
                Dlist.append(Dvalue)
    bigdeglist.append(deglist)
    bigDlist.append(Dlist)
    
# Plot data
    
fig = plt.figure()
degsub = fig.add_subplot(2,1,1)
Dsub = fig.add_subplot(2,1,2)

i = 0
handlelist = []
for p in pblist:
    deg, = degsub.plot(miglist,bigdeglist[i], label=labellist[i], color = colorlist[i])
    degsub.set_ylabel('Bed elevation change (m)')
    Dsub.plot(miglist,bigDlist[i],color = colorlist[i])
    Dsub.set_ylabel('Active Layer D50 (mm)')
    plt.xlabel('Migration rate (m/yr)')
    plt.hold(True)
    handlelist.append(deg)
    i = i + 1
fig.tight_layout()
plt.legend(handles=handlelist, loc="center", bbox_to_anchor=(1.4,1), title="AlphaBar (Initial Floodplain D50)")
    
plt.show()
