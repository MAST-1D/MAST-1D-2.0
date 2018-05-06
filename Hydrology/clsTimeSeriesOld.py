# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 16:16:22 2016

@author: geography
"""
import xlrd
import datetime
import matplotlib.pyplot as plt
import numpy as np
from copy import deepcopy

class clsTimeSeries(object):

    """
    Holds data for time series, e.g. discharge or precipitation record.
    
    Attributes:
    
    Dates--[date] (Date in time series, in matplotlib's datetime format)
    Values--[float] (Time series values)
    label--str (Label for y-axis)
    smoothed--bool (False by default; if smoothing is applied to data, becomes True)   
    averages--[[int,int,float]] (List of averages for defined periods with their date indexes)
    timechunks--[str] (list of periods to divide data)
    text -- str (optional Label for subplot)
    """
    
    def __init__(self,smoothed=0,timechunks=[]):
        self.Dates = []
        self.Values = []
        self.label = ""
        self.timechunks = timechunks
        self.averages = []
        self.text = ""
        self.Bins = []
        self.DurationCurves = []
        
        if smoothed != 0:
            self.smoothed = True
        else:
            self.smoothed = False
        
    def Load_Q_data(self,Qpath,label="",text=""):
        
        """
        Loads in discharge data from the 'Q_McDonaldBridge' Excel file and fills Dates and Values lists and label
        
        Attributes:
        Qpath--str (path to Excel file)
        label--str (name of series (for labeling y-axis of figure))
        text -- str (optional Label for subplot)
        """
        
        self.label = label 
        self.text = text
        
        book = xlrd.open_workbook(Qpath)
        sheet = book.sheet_by_name('Q')
        
        #  Converts string date to matplotlib-compatable date
        rawDates = map(lambda x: x.value,sheet.col_slice(colx=0,start_rowx = 1,end_rowx = None))
        tupleDates = map(lambda x: x.split('-'),rawDates)
        self.Dates = map(lambda x: datetime.datetime(int(x[0]),int(x[1]),int(x[2])), tupleDates)
                    
        #  Loads Q data and calculates averages for specified timechunks
        self.Values = map(lambda x: x.value, sheet.col_slice(colx=1,start_rowx=1,end_rowx = None))
        if len(self.timechunks) > 0:
            self.Calculate_time_averages()
            
    def Load_Pk_data(self,Pkpath,label,text):
        
        """
        Loads in discharge data from the 'Peak_McDonaldGauge' Excel file and fills Dates and Values lists and label
        
        Attributes:
        Qpath--str (path to Excel file)
        label--str (name of series (for labeling y-axis of figure))
        text -- str (optional Label for subplot)
        """
        
        self.label = label
        self.text = text
        
        book = xlrd.open_workbook(Pkpath)
        sheet = book.sheet_by_name('ForAnalysis')
        
        self.Dates = map(lambda x: int(x.value), sheet.col_slice(colx=1,start_rowx = 1,end_rowx = None))
        self.Values = map(lambda x: x.value, sheet.col_slice(colx=3,start_rowx = 1,end_rowx = None))

    def Load_P_data(self,Ppath,label,text):
        
        """
        Loads in precip data from 'Elwha_daily_climate_data' spreadsheet and fills Dates and Values lists and label
        
        Attributes:
        Qpath--str (path to Excel file)
        label--str (name of series (for labeling y-axis of figure))
        """
        
        self.label = label
        self.text = text
        
        book = xlrd.open_workbook(Ppath)
        sheet = book.sheet_by_name('646810')
        
        #  Converts string date to matplotlib-compatable date
        rawDates = map(lambda x: str(x.value),sheet.col_slice(colx=2,start_rowx = 2,end_rowx = None))
        self.Dates = map(lambda x: datetime.datetime(int(x[0:4]),int(x[4:6]),int(x[6:8])), rawDates)
        
        #  Loads values--deletes values of '-9999' (i.e. no data) and removes the date from the Dates list
        rawValues = map(lambda x: x.value/10., sheet.col_slice(colx=7,start_rowx=2,end_rowx = None)) # Division by 10 to conver to mm
        i = 0
        for val in rawValues:
            if val != -999.9:
                self.Values.append(val)
            else:
                self.Dates[i] = 'nan'
            i = i + 1
        self.Dates = filter(lambda x: x != 'nan', self.Dates)
        
        if len(self.timechunks) > 0:
            self.Calculate_time_averages()
        
    def Calculate_time_averages(self):
        
        """
        Calculates averages for chunks of time specified in self.timechunks
        """

        for x in self.timechunks:  
            startind, endind = self.Extract_timechunk(x[0],x[1])              
            
            meanVal = np.mean(self.Values[startind:endind])
            self.averages.append([startdate,enddate,meanVal])

    def Extract_timechunk(self,date1,date2):                    
        x1 = date1
        x2 = date2
        startdate = datetime.datetime(int(x1[0:4]),int(x1[4:6]),int(x1[6:]))
        enddate = datetime.datetime(int(x2[0:4]),int(x2[4:6]),int(x2[6:]))
        startind = self.Dates.index(startdate)
        endind = self.Dates.index(enddate)
        
        return startind, endind
    
    def Smooth_data(self,x):

        """
        Calculates a moving average of a smoothing value for a given Values dataset and
            adjusts Values and Date set        
        
        Attributes:
        x--int (Smoothing value--must be an odd number)
        """     
        
        smoothedvals = []        
        
        startind = (x-1)/2
        endind = len(self.Dates)-((x-1)/2)
        
        self.Dates = self.Dates[startind:endind]
        
        i = startind
        while i < endind:
            val = np.mean(self.Values[i-((x-1)/2):i+((x-1)/2)])
            smoothedvals.append(val)
            i = i + 1
            
        self.Values = smoothedvals
   
    def Cumu_departure(self,startyear,endyear):
        
        """
        Calculates the cumulative departure of total annual flow given a daily discharge record.
        -Sums flow for each year.
        -Calculates mean for all years.
        -Sets Values list to cumulative departure from this mean.
        Method is adapted from the USGS (Kresch, 1994)
        
        Attributes:
        startyear--int
        endyear--int
        """

        yearlist = []
        datelist = []
        
        #  Extract discharge values for each year and sum them.

        meanflow = sum(map(lambda x: x*86400*1,self.Values))/len(self.Dates) # Equation 1 in Kresch:  ((m^3/s)*(s/day)*(day)/((day)*(yr/d)) = m^3
        i = startyear
        Dxt_old = 0.
        
        while i <= endyear:
            year = i+1
            yearsum = 0.
            index1 = self.Dates.index(datetime.datetime(i,10,01)) # Water new-years
            if i == endyear:
                yearsum = sum(map(lambda x: x*86400*1,self.Values[index1:]))/len(self.Dates[index1:])
            else:
                index2 = self.Dates.index(datetime.datetime(i+1,10,01)) # One past Water new-years eve
                yearsum = sum(map(lambda x: x*86400*1,self.Values[index1:index2]))/len(self.Dates[index1:index2])
                        
            Dx = (yearsum-meanflow) # Equation 2 in Kresch (m^3)
            Dxt_new = Dxt_old + (365.25*Dx) #  Kresch Equation 3 (converting from m^3/day to m^3/yr)
            Dxt_old = Dxt_new            
            yearlist.append(Dxt_new)
            datelist.append(year)
            
            i = i + 1

        self.Values = yearlist
        self.Dates = datelist            