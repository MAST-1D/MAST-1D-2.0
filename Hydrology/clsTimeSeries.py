# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 16:16:22 2016

@author: geography
"""
#import xlrd
import datetime
import numpy as np
from copy import deepcopy
import math

class clsTimeSeries(object):

    """
    Holds data for time chunks of time series
    
    Attributes:
    
    Dates--[date] (Date in time series, in matplotlib's datetime format)
    Values--[float] (Time series values)
    label--str (Label for y-axis)  
    averages--[[int,int,float]] (List of averages for defined periods with their date indexes)
    timechunks--[str] (list of periods to divide data)
    text -- str (optional Label for subplot)
    """
    
    def __init__(self,Dates,Values, label="", text=""):
        self.Dates = Dates
        self.Values = Values
        #self.average = np.mean(self.Values)
        self.average = np.mean(list(self.Values))
        self.DurationCurves = []
        self.label = label
        self.text = text
        
##    def Load_Pk_data(self,Path,Worksheet='ForAnalysis',label="",text=""):
##        
##        """
##        Loads in discharge data from the 'Peak_McDonaldGauge' Excel file and fills Dates and Values lists and label
##        
##        Attributes:
##        Worksheet--str (Optional name of Excel worksheet; default is 'Q')
##        label--str (name of series (for labeling y-axis of figure))
##        text -- str (optional Label for subplot)
##        """
##        
##        self.label = label
##        self.text = text
##        
##        book = xlrd.open_workbook(Path)
##        sheet = book.sheet_by_name(Worksheet)
##        
##        self.Dates = map(lambda x: int(x.value), sheet.col_slice(colx=1,start_rowx = 1,end_rowx = None))
##        self.Values = map(lambda x: x.value, sheet.col_slice(colx=3,start_rowx = 1,end_rowx = None))
    
    def CreateDurationCurve(self, bins, feet=False, minQ=40., logbin=False):
        
        """
        Creates a flow duration curve for daily discharge data given a number of bins to use.  If the attribute
        feet is set to True, the function will run a ft^3/s to m^3/s conversion after the 
        data has been binned.
        
        Attributes:
        
        bins--int (number of bins to use)
        feet--bool (denotes whether discharge data is in cubic meters per second (False)
            or cubic feet per second (True))
        """
        
        binQ = []
        
        
        #  Define bin edges
        maxbin = max(self.Values)
        binlist = list(range(0,int(maxbin),int(bins)))
        
        if logbin == True:
            binlist = [0,bins]
            while binlist[-1] < maxbin:
                newbin = binlist[-1]+bins*math.log(binlist[-1],10)
                binlist.append(newbin)
            binlist[-1]=maxbin
        
        if logbin == False:
            maxbinlist = max(binlist)
            binlist.append(maxbinlist+bins) # to get the last bin
            
        cfd, binedgesnp = np.histogram(self.Values,bins=binlist)
        DCraw = cfd.tolist()
        binedges = binedgesnp.tolist()
        
        #  Convert from histogram to frequency plot
        total = float(sum(DCraw))
        DC = list(map(lambda x: x/total, DCraw))

        #  Find average value of each bin
        i = 0
        while i < len(binedges)-1:
            Q = np.mean([binedges[i],binedges[i+1]])
            binQ.append(Q)
            i = i + 1
        
        # Get rid of bins with no discharges or else MAST will crash
        newbinQ = []
        newDC = []
        for i in range(len(binQ)):
            if DC[i] != 0.:
                newbinQ.append(binQ[i])
                newDC.append(DC[i])
         
#        #  Combine bins with discharges below a threshold--for model stability
#        Trigger = False
#        for Q in binQ:
#            if Q > minQ and Trigger == False:
#                Trigger = True
#                Qind = binQ.index(Q)
#                combinedBin = float(sum(binQ[0:Qind]))/len(binQ[0:Qind])
#                combinedp = sum(DC[0:Qind])
#                
#                binQ = [combinedBin] + binQ[Qind:]
#                DC = [combinedp] + DC[Qind:]
                
                
        
        #  If discharge is in ft^3/s, convert to m^3/s
        if feet == True:
            newbinQ = map(lambda x: x*.3048**3, newbinQ)            
            
        return newbinQ, newDC
        
    

   
         
