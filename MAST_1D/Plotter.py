# -*- coding: utf-8 -*-
"""
This script contains a function that converts MAST-1D output into a form usable by ggplot.
It then makes simple plots of either an output variable with respect to time at a user-defined
node or an output variable with respect to distance downstream at a user-defined time.

Database Format:

Node    Time    Variable    Value
x       x       x           x
x       x       x           x
x       x       x           x
...

The code can be modified to add more columns later (e.g. feed, migration)

"""

import os
import numpy as npy
import pandas
from ggplot import*

"""
Inputs
"""

folder = str(os.path.join(os.pardir, "Output")) #  Folder where model output text files are stored
variable1 = 'ActiveLayer.GSD.D50' #  Variable to plot (must be a file with the prefix 'Out_'; do not include 'Out_' in name here)
xaxis = 'Time' #  Variable you want as the x-axis; 'Node' or 'Time'
locid = '0.0' #  Point in space/time--either distance downstream in meters or time in years
reachlength = 7500 #  Reach length in meters

"""
Converts input files to database
"""

def datamaker(path):

#  Import and read model output data (files with prefix 'Out_')

    filelist = []
    datalist = []
    
    variables = os.listdir(path)

    for x in variables:
        if 'Out_' in x:
            filelist.append(x)
            
    for x in filelist:
        file = open(path + '//' + x, 'r')
        file = file.readlines()
        str(file)

        time = file[1]
        times = time.split()

        del file[0:3]
 
#  Forms data into ggplot-style array form
        
        i = 0
        for t in times:
            for row in file:
                splitset = row.split()
                
#                node = (reachlength-float(splitset[0]))/1.2
#                node = (reachlength-float(splitset[0]))
                node = float(splitset[0])
                time = float(t)
                variable = x
                variable = variable[4:]
                value = float(splitset[i+1])
                
                out = [node, time, variable, value]
                
                datalist.append(out)

            i = i + 1


    da = npy.array(datalist)
    df = pandas.DataFrame(data=da, columns = ['Node', 'Time', 'Variable', 'Value'])
    return df

def gsdata(textfile):

    datalist = []    
    
    file = open(textfile, 'r')
    file = file.readlines()
    
    del file[0]
    del file[-12:-1] # Get rid of Middle reach
    del file[-1]
    
    for row in file:
        splitset = row.split()
        node = float(splitset[1])
        D50 = float(splitset[3])
        D90 = float(splitset[5])
        out = [node, 'D50', D50]
        datalist.append(out)
        out = [node, 'D90', D90]
        datalist.append(out)
        
    da = npy.array(datalist)
    df = pandas.DataFrame(data=da, columns = ['Node', 'Variable', 'Value'])
    return df
    
def longdata(textfile):

    datalist = []    
    
    file = open(textfile, 'r')
    file = file.readlines()
    
    del file[0]
    
    for row in file:
        splitset = row.split()
        node = float(splitset[6])
        minp = float(splitset[1])
        maxp = float(splitset[2])
        meanp = float(splitset[3])
        out = [node, 'min', minp]
        datalist.append(out)
        out = [node, 'max', maxp]
        datalist.append(out)
        out = [node, 'mean', meanp]
        datalist.append(out)
        
    da = npy.array(datalist)
    df = pandas.DataFrame(data=da, columns = ['Node', 'Variable', 'Value'])
    return df
        
"""
Plotting
"""

df = datamaker(folder)

sindex = ''
xlabel = ''

if xaxis == 'Node':
    sindex = 'Time'
    slabel = 'Distance Upstream (m)'
    stitle = 't = '
else:
    sindex = 'Node'
    slabel = 't (yr)'
    stitle = 'x = '
#variable2 = 'Floodplain.GSD.D50'
#variable2 = 'Floodplain.GSD.D50'

df.set_index([sindex,'Variable'], drop='False', inplace = 'True')

#  Set node/timestep to plot

df1 = df.loc[locid,variable1,:]

#  The plot

print ggplot(aes(x=xaxis, y='Value'), data=df1) +\
    geom_line(color = 'blue') + xlab(slabel) + ylab(variable1) \
    + ggtitle(stitle + locid)

"""
For field data
"""

#fdf = gsdata('D:\MAST-1D_version_K4\Field_data\Wolman_grainsize.txt')
#fdf.set_index(['Variable'], drop='False', inplace = 'True')
#fdf1 = fdf.loc['D90',:]

#fdf = longdata('D:\MAST-1D_version_K5\Field_data\Long_prof.txt')
#fdf.set_index(['Variable'], drop='False', inplace = 'True')
#fdf1 = fdf.loc['mean',:]

#print ggplot(aes(x=xaxis, y='Value'), data=df1) +\
#    geom_line(color = 'blue') + xlab(slabel) + ylab(variable1) \
#    + ggtitle(stitle + locid) + geom_line(aes(x='Node', y='Value'), data=fdf1)




