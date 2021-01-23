# -*- coding: utf-8 -*-
"""
Created on Tue May 19 13:45:32 2015

@author: geography
"""

import os
import matplotlib.pyplot as plt

folder = str(os.path.join(os.pardir, "Output"))
folder0 = str(os.path.join(os.pardir, "Output_0.0mig"))

def extractGSDdata(folder):
    inmiglist = []
    outmiglist = []
    bedlist = []
    loadlist = []
    sandplot = []
    gravelplot = []
    cobbleplot = []
    boulderplot = []
    sandplotfp = []
    gravelplotfp = []
    cobbleplotfp = []
    boulderplotfp = []
    
    mobilitylist = []
    
    qtotfile = open(folder + '//' + "Out_Load.QsavBedTot")
    qtotfile = qtotfile.readlines()
    str(qtotfile)
    qtotnode = qtotfile[3].split()
    del qtotnode[0]
    
    time = qtotfile[1]
    times = time.split()
    
    i = 1
    while i < 9:
        inmigfile = open(folder + '//' + "Out_ActiveLayer.ExSed.InMigration[" + str(i) + "]", 'r')
        inmigfile = inmigfile.readlines()
        str(inmigfile)
        inmignode = inmigfile[3].split()
        del inmignode[0]
        
        outmigfile = open(folder + '//' + "Out_ActiveLayer.ExSed.OutMigration[" + str(i) + "]", 'r')
        outmigfile = outmigfile.readlines()
        str(outmigfile)
        outmignode = outmigfile[3].split()
        del outmignode[0]
        
        bedfile = open(folder + '//' + "Out_NormalisedActiveLayer.GSD.F[" + str(i) + "]", 'r')
        bedfile = bedfile.readlines()
        str(bedfile)
        bednode = bedfile[3].split()
        del bednode[0]
        
        loadfile = open(folder + '//' + "Out_Load.GSDBedloadAv.F[" + str(i) + "]", 'r')
        loadfile = loadfile.readlines()
        str(loadfile)
        loadnode = loadfile[3].split()
        del loadnode[0]
        
        mobilitylist.append(float(loadnode[0])/float(bednode[0])) 
        
        j = 0
        
        innerinmiglist = []
        inneroutmiglist = []
        innerbedlist = []
        innerloadlist = []
        
        for t in times:
            
            inmigvalue = float(inmignode[j])
            innerinmiglist.append(inmigvalue)
            
            outmigvalue = float(outmignode[j])
            inneroutmiglist.append(outmigvalue)
            
            bedvalue = float(bednode[j])
            innerbedlist.append(bedvalue)
            
            loadvalue = float(loadnode[j])*float(qtotnode[j])
            innerloadlist.append(loadvalue)
            
            j = j + 1
            
        inmiglist.append(innerinmiglist)
        outmiglist.append(inneroutmiglist)
        bedlist.append(innerbedlist)
        loadlist.append(innerloadlist)
        
        i = i + 1
    
    j = 0    
    for t in times:
        sandinmig = inmiglist[0][j]
        gravelinmig = inmiglist[1][j] + inmiglist[2][j] + inmiglist[3][j]
        cobbleinmig = inmiglist[4][j] + inmiglist[5][j]
        boulderinmig = inmiglist[6][j] + inmiglist[7][j]
        
        sandoutmig = outmiglist[0][j]
        graveloutmig = outmiglist[1][j] + outmiglist[2][j] + outmiglist[3][j]
        cobbleoutmig = outmiglist[4][j] + outmiglist[5][j]
        boulderoutmig = outmiglist[6][j] + outmiglist[7][j]
    
        sandbed = bedlist[0][j]
        sandbed0 = bedlist[0][0]
        gravelbed = bedlist[1][j] + bedlist[2][j] + bedlist[3][j]
        gravelbed0 = bedlist[1][0] + bedlist[2][0] + bedlist[3][0]
        cobblebed = bedlist[4][j] + bedlist[5][j]
        cobblebed0 = bedlist[4][0] + bedlist[5][0]
        boulderbed = bedlist[6][j] + bedlist[7][j]
        boulderbed0 = bedlist[6][0] + bedlist[7][0]
        
        sandload = loadlist[0][j]
        gravelload = loadlist[1][j] + loadlist[2][j] + loadlist[3][j]
        cobbleload = loadlist[4][j] + loadlist[5][j] 
        boulderload = loadlist[6][j] + loadlist[7][j]

  #Outputs/inputs divided 
#        sandplot.append((sandload)/sandinmig)
#        sandplotfp.append((sandoutmig)/sandinmig)
#        gravelplot.append((gravelload)/gravelinmig)
#        gravelplotfp.append((graveloutmig)/gravelinmig)
#        cobbleplot.append((cobbleload)/cobbleinmig)
#        cobbleplotfp.append((cobbleoutmig)/cobbleinmig)
#        boulderplot.append((boulderload)/boulderinmig)
#        boulderplotfp.append((boulderoutmig)/boulderinmig)
        
  #Outputs/inputs total 
        sandplot.append((sandload+sandoutmig)/sandinmig)
        gravelplot.append((gravelload+graveloutmig)/gravelinmig)
        cobbleplot.append((cobbleload+cobbleoutmig)/cobbleinmig)
        boulderplot.append((boulderload+boulderoutmig)/boulderinmig)
    
  #Mobility
#        sandplot.append((sandload/float(qtotnode[j]))/sandbed)
#        gravelplot.append((gravelload/float(qtotnode[j]))/gravelbed)
#        cobbleplot.append((cobbleload/float(qtotnode[j]))/cobblebed)
#        boulderplot.append((boulderload/float(qtotnode[j]))/boulderbed)
    
#  Active Layer GSD
#        sandplot.append(sandbed)
#        gravelplot.append(gravelbed)
#        cobbleplot.append(cobblebed)
#        boulderplot.append(boulderbed)    
    
        j = j + 1

    return [times, sandplot, gravelplot, cobbleplot, boulderplot]
    #return [times, sandplot, gravelplot, cobbleplot, boulderplot,sandplotfp, gravelplotfp, cobbleplotfp, boulderplotfp]

data = extractGSDdata(folder)
times = data[0]
sandplot = data[1]
gravelplot = data[2]
cobbleplot = data[3]
boulderplot = data[4]  

#data0 = extractGSDdata(folder0)
#times0 = data0[0]
#sandplot0 = data0[1]
#gravelplot0 = data0[2]
#cobbleplot0 = data0[3]
#boulderplot0 = data0[4]      

# Plot data
    
fig = plt.figure()

"""
Over time
"""
#sub = fig.add_subplot(1,1,1)
#sub.plot(times[:],sandplot[:], label="Sand (fully mobile)")
#sub.set_yscale('log')
#sub.set_ylabel('pi/fi')
#sub.set_xlabel('Time (yrs)')
#plt.hold(True)
#sub.plot(times[:],gravelplot[:], label="Gravel (fully mobile)")
#plt.hold(True)
#sub.plot(times[:],cobbleplot[:], label="Cobbles (partially mobile)")
#plt.hold(True)
#sub.plot(times[:],boulderplot[:], label="Boulders (immobile)")
#fig.tight_layout()
#plt.legend(loc="center", bbox_to_anchor=(1.4,.5))
"""
Initial mobility
"""
#sub = fig.add_subplot(1,1,1)
#sizelist = [0.355,4.0,16.0,45.25,128.0,362.0,724.1,1.448e+03]
#line = [1,1,1,1,1,1,1,1]
#sub.plot(sizelist, mobilitylist)
#plt.hold(True)
#sub.plot(sizelist, line)
#sub.set_xlim(sizelist[0],sizelist[-1])
#sub.set_xscale('log')
#sub.set_ylabel('pi/fi')
#sub.set_xlabel('D (mm)')
#sub.set_title('Initial grain mobility')

"""
vs. No migration
"""
#sandplotfp = data[5]
#gravelplotfp = data[6]
#cobbleplotfp = data[7]
#boulderplotfp = data[8] 

sub = fig.add_subplot(2,2,1)
sub.plot(times[2:],sandplot[2:], label="Load")
plt.hold(True)
#sub.plot(times0[:],sandplot0[:], label="m = 0.0 m/yr")
#sub.plot(times[:],sandplotfp[:], label="Floodplain")
sub.set_ylabel('Output/Input')
#sub.set_xlabel('t (yrs)')
sub.set_title('Sand')

sub2 = fig.add_subplot(2,2,2)
sub2.plot(times[2:],gravelplot[2:], label="Load")
plt.hold(True)
#sub2.plot(times0[:],gravelplot0[:], label="m = 0.0 m/yr")
#sub2.plot(times[:],gravelplotfp[:], label="Floodplain")
#sub2.set_ylabel('Output/Input')
#sub2.set_xlabel('t (yrs)')
sub2.set_title('Gravel')

sub3 = fig.add_subplot(2,2,3)
sub3.plot(times[2:],cobbleplot[2:], label="Load")
plt.hold(True)
#sub3.plot(times0[:],cobbleplot0[:], label="m = 0.0 m/yr")
#sub3.plot(times[:],cobbleplotfp[:], label="Floodplain")
sub3.set_ylabel('Output/Input')
sub3.set_xlabel('t (yrs)')
sub3.set_title('Cobbles')

sub4 = fig.add_subplot(2,2,4)
sub4.plot(times[2:],boulderplot[2:], label="Load")
plt.hold(True)
#sub4.plot(times0[:],boulderplot0[:], label="m = 0.0 m/yr")
#sub4.plot(times[:],boulderplotfp[:], label="Floodplain")
#sub4.set_ylabel('Output/Input')
sub4.set_xlabel('t (yrs)')
sub4.set_title('Boulders')
#plt.legend(loc="center", bbox_to_anchor=(1.5,1.2))

fig.tight_layout()
  
plt.show()    
    
    
    
    