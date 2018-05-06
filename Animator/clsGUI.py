# -*- coding: utf-8 -*-
"""
Created on Mon May 08 20:14:23 2017

@author: geography
"""
import os
import Tkinter as T
import tkFileDialog
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.animation as animation
from clsAnimation import clsAnimation
import matplotlib.pyplot as plt

class clsGUI(object):
    
    def __init__(self):
        self.root = T.Tk()
        self.fig = clsAnimation()
        self.var = T.StringVar()
        self.CurrentPath = ""
        self.FileName = T.StringVar() 
        self.f = T.Frame(self.root, bd=2)
        self.ft = T.Frame(self.root, bd = 2)
        self.fname = T.Frame(self.root, bd = 2)
        self.VarsToPlot = T.OptionMenu(self.ft, self.var, *sorted(set(self.fig.variables)), command=self.updateplot).grid(column=1, row=0, sticky=T.W)
        self.Filelabel = T.Label(self.fname, text='Model Run:', font = "Times 12 bold").grid(column=0, row=0, pady=10, sticky=T.E)        
        self.FileNameLabel = T.Label(self.fname, textvariable = self.FileName, font = 'Times 12').grid(column=1, row=0, pady=10, sticky=T.E)        
        self.Plotlabel = T.Label(self.ft, text='Variable to plot:', font = "Times 12 bold").grid(column=0, row=0, pady=10, sticky=T.E)
        self.dir_opt = {}        
        self.pathbutton = T.Button(self.fname, text = "Choose run", font="Times 12", bg='LIGHTSKYBLUE', command = self.chooserun)     
        self.pausebutton = T.Button(self.f, text = u"\u2016", font=('bold'), bg='RED', command = self.pauseorplay)
        self.rewindbutton = T.Button(self.f, text = u"\u23EA", font=('bold'), command = self.rewind)
        self.ffbutton = T.Button(self.f, text = u"\u23E9", font=('bold'), command = self.fforward)
        self.aniplot = FigureCanvasTkAgg(self.fig.figure, master=self.f)
        self.error = ""
    
    def SetUpGUI(self, path):
        
        self.root.wm_title("MAST-1D Animator")
        
        self.fname.pack()
        self.ft.pack()
        self.f.pack()
        
        self.pathbutton.config(height=1, width=10)
        self.pathbutton.grid(column=3, row=0, padx=20)
 
        self.pausebutton.config(height=2, width=5)
        self.pausebutton.grid(column=1, row=3, pady=10)        
        
        self.rewindbutton.config(height=1, width = 3)
        self.rewindbutton.grid(column=0, row=3)
                
        self.ffbutton.config(height=1, width = 3)
        self.ffbutton.grid(column=2, row=3)
                
        self.aniplot.get_tk_widget().grid(column=0, row=1, padx=10, columnspan=3)        
        self.BuildPlot(path)

    def PlotlabelOptions(self):
        
        self.dir_opt['initialdir'] = self.CurrentPath
        self.dir_opt['mustexist'] = False
        self.dir_opt['title'] = 'This is a title'
        
    def FindRunName(self, path):
        splitpath = path.split("/")
        name = splitpath[-1]
        
        return name
        
    def BuildPlot(self, path):
        self.CurrentPath = path
        name = self.FindRunName(path)
        result = self.fig.datamaker(path)
        self.error = result
        
        # Check for bad paths
        if result != 'fine':
            #self.fig.AdjustPlot()
            self.fig.figure.ax.set_xlim([-5, 5])
            self.fig.figure.ax.set_ylim([-5, 5])
            self.fig.variables = ['','']
            self.var.set("")
            self.FileName.set("")
            self.VarsToPlot = T.OptionMenu(self.ft, self.var, *sorted(set(self.fig.variables)), command=self.updateplot).grid(column=1, row=0, sticky=T.W)           
            self.fig.annotation = self.fig.figure.ax.annotate(result, xy=(0,0), ha='center',va='center',size=10)
            self.fig.annotation.set_animated(True)
            
        else:
            self.fig.SetPlot(self.fig.attributes[self.fig.variables[0]])
            self.var.set(self.fig.variables[0])
            self.updateplot(self.var.get())
            self.fig.var = self.var
            self.fig.calc_t_and_nodes()
            self.FileName.set(name)      
            self.VarsToPlot = T.OptionMenu(self.ft, self.var, *sorted(set(self.fig.variables)), command=self.updateplot).grid(column=1, row=0, sticky=T.W)       
            self.PlotlabelOptions()
    
    def chooserun(self):
        filename = tkFileDialog.askdirectory(**self.dir_opt)
        
        # Write new current path
        pathfile = open(os.path.join('Animator', 'PathLastOpened.txt'), 'w')
        pathfile.write(filename)
        pathfile.close()
        
        # Update plot and GUI
        self.BuildPlot(filename)
    
    # Sets the variable selected from the dropdown list
    def updateplot(self, var):
        axes = self.fig.attributes[var]

        self.fig.figure.ax.set_xlim([axes.minx, axes.maxx])
        self.fig.figure.ax.set_ylim([axes.miny,axes.maxy])
#        for line in self.fig.figure.ax.yaxis.get_ticklines():
#            # line is a matplotlib.lines.Line2D instance
#            line.set_color('white') 
#        for line in self.fig.figure.ax.xaxis.get_ticklines():
#            # line is a matplotlib.lines.Line2D instance
#            line.set_color('white') 

        self.fig.figure.ax.set_ylabel(axes.ylabel)
        self.fig.figure.ax.set_xlabel(axes.xlabel)

        self.fig.figure.canvas.draw()

    def pauseorplay(self):
        if self.fig.pause == False:
            self.pausebutton.configure(text = u"\u25B6", bg='GREEN')
            self.fig.pause = True
        else:
            self.pausebutton.configure(text = u"\u2016", font = ('bold'), bg = 'RED')
            self.fig.pause = False
            
    def rewind(self):
        self.fig.pause = True
        self.pausebutton.configure(text = u"\u25B6", bg='GREEN')
        self.fig.j = self.fig.j - 1
        
    def fforward(self):
        self.fig.pause = True
        self.pausebutton.configure(text = u"\u25B6", bg='GREEN')
        self.fig.j = self.fig.j + 1
