# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 13:39:45 2015

@author: geography
"""

"""
A simple example of an animated plot
"""

import matplotlib.animation as animation
import os
import sys
sys.path.append("..")
import Tkinter as T
from Animator.clsGUI import clsGUI

# Find initial output folder.
pathfile = open(os.path.join('Animator', 'PathLastOpened.txt')).readlines()
path = pathfile
if len(path) > 0:
    path = pathfile[0]
else:
    path = ""
    
def animate(i):

    if 'folder' not in GUI.error:
        axes = GUI.fig.attributes[GUI.fig.var.get()]

    try:
        GUI.fig.t = GUI.fig.timestamp[GUI.fig.j*GUI.fig.nnodes][0]
    except:
        GUI.fig.t = '0'
    values = 0
    x = []
    y = []
    try:
        values = GUI.fig.data.loc[GUI.fig.t, GUI.fig.var.get(), :]
        x = values['Node'].tolist()
        y = values['Value'].tolist()
    except:
        x = [0]
        y = [0]

    # This is sort of a repeat of above but if there is an improper folder
    if 'folder' in GUI.error:
        x = [0]
        y = [0]
    
    if GUI.fig.j > GUI.fig.ntimes:
        GUI.fig.j = 0
    
    GUI.fig.line.set_data(x,y)  # update the data
    
    if 'folder' not in GUI.error: # Or else error message disappears
        GUI.fig.annotation = GUI.fig.figure.ax.annotate('t = ' + GUI.fig.t + ' yrs', xy=(axes.minx+0.75*((axes.maxx-axes.minx)),axes.miny+0.5*((axes.maxy-axes.miny))))
    
    if not GUI.fig.pause:        
        GUI.fig.j = GUI.fig.j + 1

    return GUI.fig.line, GUI.fig.annotation

# The GUI
GUI = clsGUI()
GUI.SetUpGUI(path)

# The animation
ani = animation.FuncAnimation(GUI.fig.figure, animate, init_func=GUI.fig.Aninit,frames=GUI.fig.ntimes, interval=50, blit=True)

T.mainloop()
