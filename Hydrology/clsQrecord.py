# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 13:54:29 2016

@author: geography
"""
from clsTimeSeries import clsTimeSeries
import xlrd
import datetime
import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib.lines as lines 

class clsQrecord(object):
    
    """
    Loads and holds data for entire time series, e.g. discharge or precipitation record.
    
    Attributes:
    
    Dates--[date] (Date in time series, in matplotlib's datetime format)
    Values--[float] (Time series values)
    smoothed--bool (False by default; if smoothing is applied to data, becomes True)   
    averages--[[int,int,float]] (List of averages for defined periods with their date indexes)
    recordbreaks--[str] (List of dates to divide data)
    timechunks--[clsTimeSeries] (list of divided data)
    CumDep--[float] (list of stored values for cumulative departure analysis)
    Years--[int] (list of years to correspond with CumDep)
    indexdict--{str:int} (Matches string of record chunk in form ddmmyyyyddmmyyyy to index in timechunks list)
    """
    
    def __init__(self, Path, recordbreaks=[],Qsheet='Q',smoothed=0):
        
        self.Path = Path
        self.recordbreaks = recordbreaks
        self.timechunks = []
        self.Dates = []
        self.Values = []
        self.CumDep = []
        self.Years = []
        self.indexdict = {}
        
        if smoothed != 0:
            self.smoothed = True
        else:
            self.smoothed = False
        
    def Load_Q_data(self,Worksheet = 'Q',label="",text=""):
        
        """
        Loads in discharge data from the 'Q_McDonaldBridge' Excel file and fills Dates and Values lists and label
        
        Attributes:
        Worksheet--str (Optional name of Excel worksheet; default is 'Q')
        label--str (name of series (for labeling y-axis of figure))
        text -- str (optional Label for subplot)
        """
        
        self.label = label 
        self.text = text
        
        book = xlrd.open_workbook(self.Path)
        sheet = book.sheet_by_name(Worksheet)
        
        #  Converts string date to matplotlib-compatable date
        rawDates = map(lambda x: x.value,sheet.col_slice(colx=0,start_rowx = 1,end_rowx = None))
        tupleDates = map(lambda x: x.split('-'),rawDates)
        self.Dates = map(lambda x: datetime.datetime(int(x[0]),int(x[1]),int(x[2])), tupleDates)
                    
        #  Loads Q data and calculates averages for specified timechunks
        self.Values = map(lambda x: x.value, sheet.col_slice(colx=1,start_rowx=1,end_rowx = None))

            
    def Load_Pk_data(self,Worksheet='ForAnalysis',label="",text=""):
        
        """
        Loads in discharge data from the 'Peak_McDonaldGauge' Excel file and fills Dates and Values lists and label
        
        Attributes:
        Worksheet--str (Optional name of Excel worksheet; default is 'Q')
        label--str (name of series (for labeling y-axis of figure))
        text -- str (optional Label for subplot)
        """
        
        self.label = label
        self.text = text
        
        book = xlrd.open_workbook(self.Path)
        sheet = book.sheet_by_name(Worksheet)
        
        self.Dates = map(lambda x: int(x.value), sheet.col_slice(colx=1,start_rowx = 1,end_rowx = None))
        self.Values = map(lambda x: x.value, sheet.col_slice(colx=3,start_rowx = 1,end_rowx = None))

    def Load_SD_data(self,Worksheet='Elwha_monthly_climate_data',label="",text=""):
        
        """
        Loads in discharge data from the 'Elwha_monthly_climate_data_with_snow_graphs' Excel file and fills Dates and Values lists and label
        Currently brings in max January snow depth
        
        Attributes:
        Worksheet--str (Optional name of Excel worksheet; default is 'Q')
        label--str (name of series (for labeling y-axis of figure))
        text -- str (optional Label for subplot)
        """
        
        self.label = label
        self.text = text
        
        book = xlrd.open_workbook(self.Path)
        sheet = book.sheet_by_name(Worksheet)
        
        self.Dates = map(lambda x: int(x.value), sheet.col_slice(colx=40,start_rowx = 2,end_rowx = 75))       
        self.Values = map(lambda x: x.value, sheet.col_slice(colx=41,start_rowx = 2,end_rowx = 75))

    def Load_P_data(self,Worksheet='646810',label="",text=""):
        
        """
        Loads in precip data from 'Elwha_daily_climate_data' spreadsheet and fills Dates and Values lists and label
        
        Attributes:
        Worksheet--str (Optional name of Excel worksheet; default is 'Q')
        label--str (name of series (for labeling y-axis of figure))
        text -- str (optional Label for subplot)
        """
        
        self.label = label
        self.text = text
        
        book = xlrd.open_workbook(self.Path)
        sheet = book.sheet_by_name(Worksheet)
        
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
            
    def Extract_timechunks(self): 
        
        if len(self.recordbreaks) > 0:
            i = 0
            for x in self.recordbreaks:                                  
                x1 = x[0]
                x2 = x[1]
                startdate = datetime.datetime(int(x1[0:4]),int(x1[4:6]),int(x1[6:]))
                enddate = datetime.datetime(int(x2[0:4]),int(x2[4:6]),int(x2[6:]))
                startind = self.Dates.index(startdate)
                endind = self.Dates.index(enddate)+1
                
                if enddate == self.Dates[-1]:
                    TS = clsTimeSeries(self.Dates[startind:],self.Values[startind:])
                    self.timechunks.append(TS) 
                else:
                    TS = clsTimeSeries(self.Dates[startind:endind],self.Values[startind:endind])
                    self.timechunks.append(TS) 
                
                self.indexdict[x1+x2]=i
                i = i + 1
            
        else:
            TS = clsTimeSeries(self.Dates,self.Values)
            self.timechunks.append(TS)
            
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
        index1 = self.Dates.index(datetime.date(startyear,10,01))
        index2 = self.Dates.index(datetime.date(endyear,10,01))
        meanflow = sum(map(lambda x: x*86400*1,self.Values[index1:index2+1]))/len(self.Dates[index1:index2+1]) # Equation 1 in Kresch:  ((m^3/s)*(s/day)*(day)/((day)*(yr/d)) = m^3
        i = startyear
        Dxt_old = 0.
        
        while i <= endyear:
            year = i+1
            yearsum = 0.
            
            if type(self.Dates[0])==int:
                index1 = self.Dates.index(i)
            else:
                index1 = self.Dates.index(datetime.date(i,10,01)) # Water new-years
            
            if i == endyear:
                yearsum = sum(map(lambda x: x*86400*1,self.Values[index1:]))/len(self.Dates[index1:])
            else:
                if type(self.Dates[0])==int:
                    index2 = self.Dates.index(i+1)
                else:
                    index2 = self.Dates.index(datetime.date(i+1,10,01)) # One past Water new-years eve
                yearsum = sum(map(lambda x: x*86400*1,self.Values[index1:index2]))/len(self.Dates[index1:index2])
                        
            Dx = (yearsum-meanflow) # Equation 2 in Kresch (m^3)
            Dxt_new = Dxt_old + (365.25*Dx) #  Kresch Equation 3 (converting from m^3/day to m^3/yr)
            Dxt_old = Dxt_new            
            yearlist.append(Dxt_new)
            datelist.append(year)
            
            i = i + 1

        self.CumDep = yearlist
        self.Years = datelist  
        
    def PlotDurationCurve(self,bins,feet,export=False):
        
        """
        Runs the CreateDurationCurve method in clsTimeSeries for each timechunk
        feet--bool (denotes whether discharge data is in cubic meters per second (False)
            or cubic feet per second (True))
        """
        fig = plt.figure(figsize = (10,8))
        ax = plt.subplot(1,1,1)
        legendinfo = []
        labellist = []
        colorlist = ['r','g','b']
        i=0
        
        exportDC = []  
        exportQ = []
        
        #  Add entrainment Qs
        ax.fill_between([8*10**-5,1],213.823,1000.,facecolor='.95',linewidth=0)
        #ax.fill_between([8*10**-5,1],27.,61.,facecolor='.85',linewidth=0)
        #ax.fill_between([8*10**-5,1],61.,338.,facecolor='.75',linewidth=0)
        #ax.fill_between([8*10**-5,1],338.,1000.,facecolor='.65',linewidth=0)

        #plt.text(10**-4,700,'Boulder transport',size=16,color='.35')
        #plt.text(10**-4,30.,'Gravel entrainment',size=16,color='.35')
        #plt.text(10**-4,65.,'Cobble entrainment',size=16,color='.35')
        #plt.text(10**-4,700.,'Boulder entrainment',size=16,color='.35')

        ax.axhline(107,color='.75',linewidth=3,ls='dashed') # Approx. 1 year flood
        ax.text(10**-4,80.,'Approx. 1 year flood',size=16, color='.55')

        def tick_function(D):
            Q = -4*10**-6*D**3+.003*D**2+.2841*D+11.594
            return Q
            #return ["%.3f" % z for z in Q]
        
        for x in self.timechunks:
            cumubin = []
            
            binQraw, DC = x.CreateDurationCurve(bins,feet)
            total = 0
            for val in DC:
                cumubin.append(1-total)
                total = total + val
                
            exportDC.append(cumubin) # for adding things to the plot in other modules
            
            # Convert binQ to m^3/s
            #binQ = map(lambda x: x*.0283, binQraw)
            binQ = binQraw
            #print max(binQ)
            exportQ.append(binQ)
            ax.plot(cumubin, binQ, color=colorlist[i],linewidth=1.5)
            plt.hold=True
            minyear = min(x.Dates)
            maxyear = max(x.Dates)
            label = str(minyear.year) + '-' + str(maxyear.year)
            labellist.append(label)
            legform = lines.Line2D([],[],color=colorlist[i],label=label,linewidth=1.5)
            legendinfo.append(legform)
            i=i+1
        
        ax.set_ylabel('Q (m$^3$/s)')
        ax.set_xlabel('Frequency Exceeded')
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_xlim(8*10**-5,1)
        ax.set_ylim(1,1000)
        ax2 = ax.twinx()
        
        Dtickvals = np.array([2, 16, 64, 128, 256, 512])
        Dticklocs = tick_function(Dtickvals)
        Dticklocs = map(lambda x: math.log(x,10), Dticklocs.tolist())

        ax2.set_ylim(0,3)
        ax2.set_yticks(Dticklocs)
        ax2.set_yticklabels(map(lambda x: str(x),Dtickvals))
        #ax2.set_yscale('log')

        ax2.set_ylabel('Maximum particle size entrained (mm)')        
  
        plt.rcParams.update({'font.size': 16})
        ax.legend(handles=legendinfo)
#        plt.show()
        
        if export == True:
            return fig, ax, exportDC, exportQ
        
            
            
            