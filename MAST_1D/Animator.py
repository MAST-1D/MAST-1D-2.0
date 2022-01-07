# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 13:39:45 2015

@author: geography
"""

"""
A simple example of an animated plot
"""

import numpy as npy
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import os
import sys
sys.path.append("..")
#from Tkinter import *
import tkinter as tk
import pandas

if __name__ == "__main__":
    # User inputs output folder here.
    root = tk.Tk()
    root.withdraw()
    path = tk.filedialog.askdirectory(title="Select Folder Where Output Data are Located")
    root.destroy()

    ##################################################

    class clsOutputdata(object):
        
        def __init__(self):
            self.variablename = ''
            self.minx = 0.
            self.maxx = 0.
            self.miny = 0.
            self.maxy = 0.
            self.xlabel = ""
            self.ylabel = ""
            
    ##################################################
    def datamaker(path):

    #  Import and read model output data (files with prefix 'Out_')

        filelist = []
        datalist = []
        varlist = []
        objectdict = {}
        
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


        da = npy.array(datalist)        
        df = pandas.DataFrame(data=da, columns = ['Node', 'Time', 'Variable', 'Value'])
        df['Node'] = df['Node'].astype(float)
        df['Value'] = df['Value'].astype(float)
        
        df.set_index(['Variable'], drop='False', inplace = True)
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

        df.reset_index(inplace=True)
            
        return([df,objectdict]) # for matplotlib
    ####################################################

    #Init only required for blitting to give a clean slate.
    def init():
        line.set_data([],[])
        return line, annotation

    def animate(i):
        global pause
        global annotation
        global j
        axes = attributes[var.get()]
        t = timestamp[j*nnodes][0]
        values = data.loc[t, var.get(), :]
        x = values['Node'].tolist()
        y = values['Value'].tolist()
        
        if j > ntimes:
            j = 0
        line.set_data(x,y)  # update the data
        annotation = ax.annotate('t = ' + t + ' yrs', xy=(axes.minx+0.75*((axes.maxx-axes.minx)),axes.miny+0.5*((axes.maxy-axes.miny))))
        if not pause:        
            j = j + 1
        
        return line, annotation

    def updateplot(var):
        axes = attributes[var]
        line.axes.set_ylim([axes.miny,axes.maxy])
        ax.set_ylabel(axes.ylabel)
        ax.figure.canvas.draw()
        19
    def pauseorplay():
        global pause
        if pause == False:
            pausebutton.configure(text = u"\u25B6", bg='GREEN')
            pause = True
        else:
            pausebutton.configure(text = u"\u2016", font = ('bold'), bg = 'RED')
            pause = False
            
    def rewind():
        global pause
        global j
        pause = True
        pausebutton.configure(text = u"\u25B6", bg='GREEN')
        j = j - 1
        
    def fforward():
        global pause
        global j
        pause = True
        pausebutton.configure(text = u"\u25B6", bg='GREEN')
        j = j + 1    
            
    # Data extraction
    #extraction = datamaker(str(os.path.join(os.pardir, "Output"))) # Make it draw into inputs file
    #extraction = datamaker('C:\Users\Katie\Dropbox\MAST-1D_version_K11\Output\Removal_hydrograph\RemovalTestControlledWidening_newC')
    extraction = datamaker(path)
    data = extraction[0]
    attributes = extraction[1]
    #variablelist = attributes.keys()
    variablelist = list(attributes.keys())

    #  Animation plot
    figure = plt.Figure()
    ax = figure.add_subplot(1,1,1)

    axes = attributes[variablelist[0]]
    ax.set_xlim([axes.minx, axes.maxx])
    ax.set_ylim([axes.miny,axes.maxy])
    ax.set_xlabel(axes.xlabel)
    ax.set_ylabel(axes.ylabel)
    figure.subplots_adjust(left=0.2)

    line, = ax.plot([], [])
    annotation = ax.annotate('t = ' + "" + ' yrs', xy=(axes.minx+0.75*((axes.maxx-axes.minx)),axes.miny+0.5*((axes.maxy-axes.miny))))
    annotation.set_animated(True)

    # GUI
    root = tk.Tk()
    root.wm_title("MAST-1D vK3 animator")
    var = tk.StringVar()
    var.set(variablelist[0])

    f = tk.Frame(root, bd=2)
    ft = tk.Frame(root, bd = 2)
    ft.pack()
    f.pack()
    VarsToPlot = tk.OptionMenu(ft, var, *sorted(set(variablelist)), command=updateplot).grid(column=1, row=0, sticky=tk.W) 
    Plotlabel = tk.Label(ft, text='Variable to plot:', font = ('12' + 'bold')).grid(column=0, row=0, pady=10, sticky=tk.E)  

    aniplot = FigureCanvasTkAgg(figure, master=f)
    aniplot.get_tk_widget().grid(column=0, row=1, padx=10, columnspan=3)

    pause = False
    pausebutton = tk.Button(f, text = u"\u2016", font=('bold'), bg='RED', command = pauseorplay)
    pausebutton.config(height=2, width=5)
    pausebutton.grid(column=1, row=3, pady=10)

    rewindbutton = tk.Button(f, text = u"\u23EA", font=('bold'), command = rewind)
    rewindbutton.config(height=1, width = 3)
    rewindbutton.grid(column=0, row=3)

    ffbutton = rewindbutton = tk.Button(f, text = u"\u23E9", font=('bold'), command = fforward)
    ffbutton.config(height=1, width = 3)
    ffbutton.grid(column=2, row=3)

    # Find out how many records there are
    data.set_index(['Node', 'Variable'], drop='False', inplace = True)
    values = data.loc[0, var.get(),:]
    ntimes = len(values.index)
    data.reset_index(inplace=True)

    # Find out how many nodes there are
    data.set_index(['Time', 'Variable'], drop='False', inplace = True)
    values = data.loc['0.0', var.get(),:]
    nnodes = len(values.index)
    timestamp = list(data.index.values)

    # The animation
    j = 0
    ani = animation.FuncAnimation(figure, animate, init_func=init,frames=ntimes, interval=50, blit=True)

    tk.mainloop()
