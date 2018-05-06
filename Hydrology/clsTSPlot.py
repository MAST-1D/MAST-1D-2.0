# -*- coding: utf-8 -*-
"""
Created on Tue Mar 01 18:36:06 2016

@author: geography
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 10:40:16 2016

Handles hydrologic and climate data for Elwha

@author: geography
"""

import xlrd
import datetime
import matplotlib.pyplot as plt
import numpy as np
from copy import deepcopy
from clsTimeSeries import clsTimeSeries

class clsTSPlot(object):
    
    def __init__(self):
        
        """
        Attributes:

        Series--[clsTimeSeries] (List that stores time series objects)
        """
        
        self.Series = []

    def Add_Series(self,path, Dtype,label, smoothed=0,timechunks=[],text = ""):
        
        """
        Creates clsTimeSeries object, loads the data to it, and adds it to the Series list
        
        Attributes:
        
        path--str (path where file is stored)
        Dtype--str (type of data in order to reference the correct function to load the data)
            -'Q' = Daily Discharge
            -'P' = Daily Precipitation
            - smoothed--int (default is 0; if it is not 0, data is smoothed with smoothing function and plotted; value must be odd)
            - timechunks--[str] (list of periods to divide data)
            - text -- str (optional Label for subplot)
        """
        
        series = clsTimeSeries(smoothed,timechunks)
        
        if Dtype == 'Q':
            series.Load_Q_data(path,label,text)
            if smoothed != 0:
                series.Smooth_data(smoothed)
                self.smoothed = True
            self.Series.append(series)
            
        elif Dtype == 'P':
            series.Load_P_data(path,label=label,text=text)
            if smoothed != 0:
                series.Smooth_data(smoothed)
                self.smoothed = True
            self.Series.append(series)
            
        elif Dtype == 'Pk':
            series.Load_Pk_data(path,label=label,text=text)
            if smoothed != 0:
                series.Smooth_data(smoothed)
                self.smoothed = True
            self.Series.append(series)
            
        else:
            print 'Improper spreadsheet format'
            
            
    def Plot_series(self,shading=[],smoothed=0):
        
        """
        Creates plot of time series
        
        Attributes:
        
        shading--[[int],[int]] (For shading date ranges.  Contains lists of start-end date pairs with the format 'yyyymmdd')
        smoothed--int (default is 0; if it is not 0, data is smoothed with smoothing function and plotted; value must be odd)
        """
        
        fig, axes = plt.subplots(nrows=len(self.Series), ncols=1, sharex=True,sharey=False)
        
        i = 0        
        if len(self.Series) > 1:
            for ax in axes:
                data = self.Series[i]
                self.Plot_func(ax,data,ylabel=data.label,shading=shading)
                if smoothed != 0 and data.smoothed == False:
                    smoothdata = deepcopy(data)
                    smoothdata.Smooth_data(smoothed)
                    ax.plot(smoothdata.Dates,smoothdata.Values,color='r',linewidth=2)
                i = i + 1
        else:
            data = self.Series[i]
            self.Plot_func(axes,self.Series[i],ylabel=data.label,shading=shading)
            if smoothed != 0 and data.smoothed == False:
                    smoothdata = deepcopy(self.Series[i])
                    smoothdata.Smooth_data(smoothed)
                    ax.plot(smoothdata.Dates,smoothdata.Values)
        
        plt.xlabel('Date')
        plt.show
        
        return fig
        
    def Plot_func(self,ax,data,ylabel,shading):
        ax.plot(data.Dates,data.Values)
        ax.set_ylabel(ylabel)
        
        if len(shading) > 0:
            for x in shading:
                x1 = x[0]
                x2 = x[1]
                if type(x1)==int: #  Filters for whether date is in years or day/month/year
                    ax.axvspan(x1,x2, color='lightgrey', alpha=0.5)
                else:    
                    ax.axvspan(datetime.datetime(int(x1[0:4]),int(x1[4:6]),int(x1[6:])),datetime.datetime(int(x2[0:4]),int(x2[4:6]),int(x2[6:])), color='lightgrey', alpha=0.5)
            
        ax.text(0.01,.9,data.text,transform = ax.transAxes,fontweight='bold')

