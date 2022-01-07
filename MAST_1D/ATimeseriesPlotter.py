# -*- coding: utf-8 -*-
"""
Created on Tue Aug 24 19:27:09 2021

@author: lauerj
"""
import tkinter as tk
import pandas as pd
import os
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from functools import partial

class App:
    """
    A GUI for easily cycling through output timeseries variables.  
    
    Allows the user to interactively select any number of variables from 
    the output dataframe that are then plotted against time in a multi-
    part figure.  
    
    Parameters
    ----------
    root : object
        Tkinter root object.
    df : object
        Pandas dataframe containing time-series output. Index must be set to
        date-time column.
    plots : integer
        Number of sub-plots to graph.
    """
    def __init__(self, root, df, plots):
        # Create a container
        frame = tk.Frame(root)
        self.df = df
        variables = tuple(df.columns.values)
        self.fig = Figure()
        
        #initialize lists
        options_var, label, menu, self.ax, self.line = [None]*plots,[None]*plots,[None]*plots,[None]*plots,[None]*plots

        for i in range(plots):
            options_var[i] = tk.StringVar(root)
            options_var[i].set(variables[0])

            label[i] = tk.Label(root,  text='Select dataset ' + str(i+1)+ ' to plot')
            label[i].pack()       

            menu[i] = tk.OptionMenu(root, options_var[i], *variables, command=partial(self.callback,i))
            menu[i].pack()

            self.ax[i] = self.fig.add_subplot(plots,1,i+1)
            self.line[i], = self.ax[i].plot(df.index,df[variables[0]])

        self.canvas = FigureCanvasTkAgg(self.fig,master=root)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side='top', fill='both', expand=1)
        frame.pack()
    
    def callback(self,i,selection):
        self.line[i].set_ydata(self.df[selection])
        x, y = self.line[i].get_data()
        self.ax[i].set_ylim(min(y),max(y))
        self.ax[i].set_title(selection)
        self.fig.tight_layout()
        self.canvas.draw()

if __name__ == "__main__":
    root = tk.Tk()
    root.withdraw()
    path = tk.filedialog.askdirectory(title="Select Folder Where Output Data are Located")
    root.destroy()
    
    df = pd.read_csv(os.path.join(path, "DailyOutput.csv"), 
                              parse_dates=True, 
                              index_col='time_year')
    root = tk.Tk()
    plots = 3
    app = App(root,df,plots)
    root.mainloop()