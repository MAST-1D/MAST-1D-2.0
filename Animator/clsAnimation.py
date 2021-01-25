# -*- coding: utf-8 -*-
"""
Created on Mon May 08 15:45:43 2017

@author: geography
"""
import os
import numpy as np
import pandas
from clsOutputdata import clsOutputdata
import matplotlib.pyplot as plt
import tkinter as T

class clsAnimation():
    
    def __init__(self):
        self.figure = plt.Figure()
        self.figure.ax = self.figure.add_subplot(1,1,1)
        for line in self.figure.ax.yaxis.get_ticklines():
            # line is a matplotlib.lines.Line2D instance
            line.set_color('white') 
        for line in self.figure.ax.xaxis.get_ticklines():
            # line is a matplotlib.lines.Line2D instance
            line.set_color('white') 
        
        self.line, = self.figure.ax.plot([], [])
        self.annotation = "k"
        self.data = [] # pandas dataframe
        emptyatt = clsOutputdata()
        self.attributes = {'empty':emptyatt} # Outputdata objects
        self.variables = ['',''] # variables
        self.var = T.StringVar() # current variable
        self.var.set('empty')
        self.t = "" # time
        self.j = 0 # time tracker
        self.nnodes = 1
        self.ntimes = 1
        self.timestamp = [['0'],['0']]
        self.values = []
        self.pause = False
        self.CheckPlot = False
        self.AdjustPlot()
    
    def datamaker(self, path):
    
    #  Import and read model output data (files with prefix 'Out_')
    
        filelist = []
        datalist = []
        varlist = []
        objectdict = {}
        
        try:
            variables = os.listdir(path)
        except:
            BadPathMessage = 'Please choose an output folder.'
            return BadPathMessage

        for x in variables:
            if 'Out_' in x:
                filelist.append(x)
        
        if len(filelist) == 0:
            BadPathMessage = 'Improper folder. Please choose a MAST-1D Output folder.'
            return BadPathMessage
        
        for x in filelist:
    
            file = open(path + '//' + x, 'r')
            file = file.readlines()
            str(file)
    
            time = file[1]
            times = time.split()
            
            variable = x
            variable = variable[4:]
    
            del file[0:3]
     
    #  Forms data into ggplot-style array form
            
            i = 0
            for t in times:
                for row in file:            
                    splitset = row.split()
                    
                    if '-' not in splitset:                
                        node = float(splitset[0])
                        time = float(t)
                        value = float(splitset[i+1])                
                        out = [node, time, variable, value]                
                        datalist.append(out)
    
                i = i + 1    
            varlist.append(variable)
    
    
        da = np.array(datalist)        
        df = pandas.DataFrame(data=da, columns = ['Node', 'Time', 'Variable', 'Value'])
        df['Node'] = df['Node'].astype(float)
        df['Value'] = df['Value'].astype(float)
        
        df.set_index(['Variable'], drop=False, inplace = True)
        vals = df.loc[varlist[0],:]
        minx = min(vals['Node'].tolist())
        maxx = max(vals['Node'].tolist())
        
        for var in varlist:
            vals = df.loc[var,:]
            miny = min(vals['Value'].tolist())
            maxy = max(vals['Value'].tolist())
            
            varobject = clsOutputdata()
            varobject.minx = minx
            varobject.maxx = maxx
            varobject.miny = miny
            varobject.maxy = maxy
            varobject.xlabel = 'Distance downstream (m)'
            varobject.ylabel = var
            
            objectdict[var] = varobject
    
        #df.reset_index(inplace=True)
        self.data = df
        self.attributes = objectdict
        self.variables = objectdict.keys()
        
        return 'fine'
            
    ####################################################
    def SetPlot(self, attribute):        

        self.figure.ax = self.figure.add_subplot(1,1,1)
        self.figure.ax.set_xlim([attribute.minx, attribute.maxx])
        self.figure.ax.set_ylim([attribute.miny,attribute.maxy])
#        for line in self.figure.ax.yaxis.get_ticklines():
#            # line is a matplotlib.lines.Line2D instance
#            line.set_color('white') 
#        for line in self.figure.ax.xaxis.get_ticklines():
#            # line is a matplotlib.lines.Line2D instance
#            line.set_color('white') 

        self.figure.ax.set_xlabel(attribute.xlabel)
        self.figure.ax.set_ylabel(attribute.ylabel)
        self.line.set_data(0,0)
        #self.AdjustPlot
        
        
        annotation = self.figure.ax.annotate('t = ' + "" + ' yrs', xy=(attribute.minx+0.75\
            *((attribute.maxx-attribute.minx)),attribute.miny+0.5*((attribute.maxy-\
            attribute.miny)))) 
        self.annotation = annotation
        self.annotation.set_animated(True)
    
    # Only adjusts plot if it hasn't been already
    def AdjustPlot(self):
        if self.CheckPlot == False:
            #self.CheckPlot = True
            self.figure.subplots_adjust(left=0.2)
    
    #Init only required for blitting to give a clean slate.
    def Aninit(self):
        self.line.set_data([],[])
        return self.line, self.annotation

        
    def calc_t_and_nodes(self):
        # Find out how many records there are
        self.data.set_index(['Node', 'Variable'], drop=False, inplace = True)
        values = self.data.loc[0, self.var.get(),:]
        self.ntimes = len(values.index)
        #self.data.reset_index(inplace=True)
        
        # Find out how many nodes there are
        self.data.set_index(['Time', 'Variable'], drop=False, inplace = True)
        values = self.data.loc['0.0', self.var.get(),:]
        self.nnodes = len(values.index)
        self.timestamp = list(self.data.index.values)