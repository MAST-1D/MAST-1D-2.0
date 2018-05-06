# -*- coding: utf-8 -*-
"""
Created on Fri Feb 13 12:21:25 2015

@author: geography
"""

from Tkinter import *
import os

from clsInputs import clsInputs
from ReachTester import clsModel
import numpy as np
import pickle as pic

root = Tk()
root.wm_title("MAST-1D vK1")

font = 'Times'
size = '12'

class clsFrame(object):
    """
    Attributes:
        pf--str (root frame)
        f--frame (the object's frame)
        row--int (row index of frame)
        column--int (column index of frame)
        nrows--int (number of rows in frame)
        ncolumn--int (number of columns in frame)
        ent--clsEntry
        cr--frame (child frame with f as parent)
        var--str (variable)
        i--int (left-hand column index)
        box--entry box
        par--frame (parent frame for entry)
        ind--int (order index for frame)
        buttonlist--[[str],[int]] (list of variable options and indexes)
        reflist--[widget references] (list tracking references of all entries, buttons,etc.)
    """

    def __init__(self, pf, receiver, row, column, nrows, ncolumns):
        self.pf=pf
        self.receiver = receiver # class receiver instance indexing where links should go        
        self.f=Frame(pf, bd=2, relief = 'groove')
        self.f.grid(row=row, column=column, sticky=N+S+E+W)
        self.row=row
        self.column=column
        self.nrows=nrows
        self.ncolumns=ncolumns
        self.loc=["a"]*(nrows*ncolumns)
        self.i = 0
        self.c = 0
        self.reflist = []
        
    def addFrame(self, row, column, nrows, ncolumns, ind): 
        self.loc[ind]=clsFrame(self.f,self.receiver, row, column, nrows, ncolumns)
        self.loc[ind].f.grid(row=row, column=column)
        self.i = self.i+1
        
    def addEntry(self, var, alias):
        ent = clsEntry(self.f, var, self.i,self.c)
        ent.labelmaker()
        ent.entrymaker()
        self.i = self.i+1
        self.receiver.entrylist.append([alias, ent.box])
 
        
    def addRadiobutton(self, var, alias, buttonlist):
        ent = clsEntry(self.f, var, self.i, self.c)
        ent.labelmaker()
        j = 0
        for b in buttonlist:
            ent.radiomaker(var, buttonlist,j)
            #self.i=self.i+1
            self.c = 2 #  Add if you want to indent
            exec b[2]
            self.i = self.i+1
            self.c = 0  #  Add if you want to indent
            ent.i=ent.i+2            
            j = j + 1
            self.receiver.radiolist.append([alias, ent.vartag])

            
    def addButton(self, alias, var):
        ent = clsEntry(self.f, var, self.i, self.c)
        if self.c == 0:
            ent.labelmaker()
        ent.buttonmaker(alias, var, self.receiver)
        self.receiver.windowlist.append([alias,ent])
        self.i = self.i+1

        
    def addHeading(self, var):
        ent = clsEntry(self.f, var, self.i, self.c)
        ent.headingmaker()
        self.i = self.i+1
        
    def configurecolumns(self):
        k = 0
        for col in range(0, self.ncolumns):
            self.f.columnconfigure(k,weight=1)
            k = k + 1

class clsEntry(object):
    
    """
    Attributes:
        f--frame
        var--variable
        i--left-hand column index
        box--entry box
        vartag--object naming the variable
        buttonlist--[[str],[int]] (list of variable options and indexes)
    """
    def __init__(self, f, var, i, c):
        """
        Arguments:
            f--frame
            var--variable
            i--left-hand column index
            box--the widget
        """
        self.f=f
        self.var=var
        self.i=i
        self.c=c
        self.vartag = BooleanVar()
        self.box = StringVar()
        self.active = False
        self.text = ""
        
    
    def labelmaker(self):
        Label(self.f, text=self.var).grid(row=self.i, column=self.c,sticky=E)
 
    def entrymaker(self):        
        self.box = Entry(self.f)
        self.box.grid(row = self.i, column=self.c+1, sticky=W)
        
    def radiomaker(self, alias, buttonlist,j):
        self.box = Radiobutton(self.f, text=buttonlist[j][0], variable=self.vartag, value=buttonlist[j][1])
        self.box.grid(row=self.i, column=1, sticky=W)
        
    def buttonmaker(self, alias, var, receiver):
        self.box = Button(self.f, text='Configure', command = lambda alias=alias: self.messageWindow(receiver, alias, var))
        self.box.grid(row = self.i, column = self.c+1, sticky=W)
        receiver.buttonlist.append([self.f,var, self.i, self.c+1])
            
    def headingmaker(self):
        Label(self.f, text=self.var, justify=CENTER, font=(font + " " + size + " " + 'bold')).grid(row=self.i,columnspan=2)        

    def messageWindow(self, receiver, alias, text): # Make save function that both fills out input class and reloads text in new window.
        win = Toplevel()
        message = text
        Label(win, text=message).pack()
        self.y = Text(win, height=20, width=30) 
        self.y.pack()
        if self.active == True:
            self.y.insert(END,self.text)
        Button(win, text='Ok', command=lambda alias=alias:receiver.windowsave(alias, self)).pack()
        
        
class clsReceiver(object):
    
    def __init__(self, inputs):
        self.entrylist = [] # [alias, widget] List that all entries will go into to run them to the model
        self.radiolist = []        
        self.inputs = inputs # Input class where the entries will be stored
        self.buttonlist = [] #  Where button instances will be stored
        self.windowlist = [] #  [[alias], widget] Where tab-separated data from pop-up windows will be saved
        
    def setvals(self): 
        for entry in self.entrylist:
            setattr(self.inputs, entry[0], eval(entry[1].get()) )
        for entry in self.radiolist:
            setattr(self.inputs, entry[0], entry[1].get())
        for entry in self.windowlist:
            splitset = entry[1].text.split('\n')
            numentries = len(splitset)
            alias0 = [0.]*numentries
            alias1 = [0.]*numentries
            n0 = 0
            n1 = 0
            for splits in splitset:
                row = splits.split('\t')
                num = len(row)
                if row == ['']:
                    del alias0[-1]
                    del alias1[-1]
                else:
                    alias0[n0] = float(row[0])
                    n0 = n0 + 1
                    if num > 1:
                        alias1[n1] = float(row[1])
                        n1 = n1 + 1
                    else:
                        del alias1[-1] # For Dbdy/F groupings
            
            setattr(self.inputs, entry[0][0], alias0)
            setattr(self.inputs, entry[0][1], alias1)
            
    def windowsave(self, alias, ent): #need to get these to entrylist  Save to inputs from entrylist
        tlist = ent.y.get("1.0",'end-1c')
        ent.active = True
        ent.text = tlist
        
    def fillinsaved(self, savedinputs):
        entries = savedinputs.__dict__
        for entry in self.entrylist:
            key = entry[0]
            ent = entry[1]
            ent.delete(0,END)
            ent.insert(0, entries.get(key))
        for entry in self.radiolist:
            key = entry[0]
            entry[1].set(entries.get(key))
        for entry in self.windowlist:
            entry[1].active=True
            text = ""
            key1 = entry[0][0]
            key2 = entry[0][1]
            var1 = entries.get(key1)
            var2 = entries.get(key2)
            n0 = 0
            n1 = 0
            while n0 < len(var1):
                if len(var1)>len(var2) and n0 == 0:
                    text = str(var1[n0]) + '\n'
                    n0 = n0 + 1
                else:
                    text = text + str(var1[n0]) + '\t' + str(var2[n1]) + '\n'
                    n0 = n0 + 1
                    n1 = n1 + 1
            text = text.rstrip()
            entry[1].text=text
                    

###############################################################################

inputs = clsInputs()
receiver = clsReceiver(inputs)

"""
Main frame
"""
mainframe = clsFrame(root,receiver,0,0,3,1)

"""
Top frame:  run parameters
"""
runparams = clsFrame(mainframe.f,receiver,0,0,3,2)
runparams.addHeading('Run parameters')
runparams.addEntry('Timestep (yr)', 'dt')
runparams.addEntry('Runtime (in timesteps)', 'MaxSteps')
runparams.addEntry('Number of nodes', 'Nnodes')
runparams.addEntry('Disturbances', 'LoadFactor')
runparams.addEntry('Disturbance time', 'LoadFactorCount')

runparams.configurecolumns()

"""
Middle frame:  Output specifications
"""
outputspecs = clsFrame(mainframe.f,receiver, 1,0,3,2)
outputspecs.addHeading('Output specifications')

outputspecs.addEntry('Output', 'Outputvars')
outputspecs.addEntry('Number of Printouts', 'NumberOfPrintouts')
outputspecs.addEntry('Number of Animations', 'NumberOfAnimations')

outputspecs.configurecolumns()

"""
Lower frame:  fluvial fluvial
"""
fluvial = clsFrame(mainframe.f,receiver,2,0,2,2)
fluvial.addHeading('Fluvial inputs')
fluvial.f.columnconfigure(0,pad=50)
fluvial.f.columnconfigure(1,pad=50)
fluvial.configurecolumns()

"""
Channel geometry
"""
fluvial.addFrame(0,0,7,4,0)
fluvial.loc[0].f.config(relief='flat')
fluvial.loc[0].addHeading('Channel geometry')
#fluvial.loc[0].addRadiobutton('Cross-section type', 'dummy alias', [['Complex',False, "fluvial.loc[0].addButton('dummy alias','Cross-section')"],['Rectangular',True,""]])
fluvial.loc[0].addEntry('Channel width (m)', 'Bc')
fluvial.loc[0].addEntry('Reach length', 'reachlength')
fluvial.loc[0].addEntry('Floodplain width (m)', 'Bf')
fluvial.loc[0].addEntry('Bed slope', 'Slope')
fluvial.loc[0].addEntry('Migration rate (m/yr)', 'migration')
fluvial.loc[0].addEntry('Sinuosity', 'ChSin')

"""
Hydraulics
"""
#fluvial.addFrame(1,0,5,4,2)
fluvial.loc[0].f.config(relief='flat')
fluvial.loc[0].addHeading('Hydraulics')
fluvial.loc[0].addButton(['Qw','p'],'Flow duration curve')
fluvial.loc[0].addRadiobutton('Velocity function', 'vfunc', [['Chezy',False, "fluvial.loc[0].i=fluvial.loc[0].i+1"],['Manning',True,"fluvial.loc[0].i=fluvial.loc[0].i+1"]])
fluvial.loc[0].addEntry('Cfc', 'Cfc')
fluvial.loc[0].addEntry('Cff', 'Cff')

fluvial.loc[0].configurecolumns

"""
Sediment
"""
fluvial.addFrame(0,1,9,3,1)
fluvial.loc[1].f.config(relief='flat')
fluvial.loc[1].addHeading('Sediment')
fluvial.loc[1].addRadiobutton('Sediment transport function', 'TrinityFit', [['Wilcock and Crowe',False, "fluvial.loc[1].i=fluvial.loc[1].i+1"],['Trinity River',True,"fluvial.loc[1].i=fluvial.loc[1].i+1"]])
fluvial.loc[1].addButton(['Dbdy','Fa'],'Active layer GSD')
fluvial.loc[1].addButton(['Dbdy','Fp'],'Floodplain GSD')
fluvial.loc[1].addEntry('Fraction of sand in suspension', 'FSandSusp')
fluvial.loc[1].addEntry('Floodplain depth (m)', 'FloodplainL')
fluvial.loc[1].addEntry('Active layer depth (m)', 'ActiveLayerL')
fluvial.loc[1].addEntry('Number of substrate layers', 'NLayers')
fluvial.loc[1].addEntry('Thickness of substrate layers (m)', 'LayerL')
fluvial.loc[1].addEntry('Point bar height (m)', 'Hpb')
fluvial.loc[1].addEntry('Sigma', 'sigma')
fluvial.loc[1].addEntry('Porosity', 'lambdap')

"""
Parameters
"""
#fluvial.addFrame(1,1,6,2,3)
fluvial.loc[1].f.config(relief='flat')
fluvial.loc[1].addHeading('Parameters')
fluvial.loc[1].addEntry('Kbar', 'Kbar')
fluvial.loc[1].addEntry('Alphabar', 'AlphaBar')
fluvial.loc[1].addEntry('Alphapartlyalluvial', 'AlphaPartlyAlluvial')
fluvial.loc[1].addEntry('ncAddons', 'ncAddons')
fluvial.loc[1].addEntry('ncMultiplier', 'ncMultiplier')
fluvial.loc[1].addEntry('nf', 'nf')

fluvial.loc[1].configurecolumns

"""
Lowest frame:  Tracers
"""

tracers = clsFrame(mainframe.f,receiver,3,0,1,2)
tracers.f.columnconfigure(0,weight=1)
tracers.f.columnconfigure(1,weight=1)
tracers.addHeading('Tracers')

tracers.addEntry('Process production rates', 'coj')
tracers.addEntry('Attenuation rates', 'Lcj')
tracers.addEntry('Nuclide name', 'Name')
tracers.addEntry('Decay constant', 'DecayConst')
tracers.addEntry('Total production rate', 'ProductionRate')
tracers.addEntry('Fallout rate', 'FalloutRate')

tracers.configurecolumns

"""
Action buttons
"""

def testcomb():
    savedinputs = pic.load(open('inputparams.inputs', "rb"))
    receiver.fillinsaved(savedinputs)
    
def setmodel():
    receiver.setvals()
    
def runmodel():
    run = clsModel(inputs)
    run.plotresponse.trace("w", instantplotter(run.plotresponse))
    run.main()
    print ('Model run complete!')
    
def instantplotter(thing):
    print 'Ping!'
    print thing

action = clsFrame(mainframe.f,receiver,4,0,1,3)

x=Button(action.f, text = 'Load inputs', command=testcomb)
x.grid(row=0, column=0)
y=Button(action.f, text = 'RUN', command=runmodel, bg='Green')
y.grid(row=0, column=2)
z = Button(action.f, text = 'set vals', command = setmodel)
z.grid(row = 0, column=1)

action.f.columnconfigure(0,weight=1)
action.f.columnconfigure(1,weight=1)
action.f.columnconfigure(2,weight=1)


mainloop()