# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 17:38:44 2016

@author: geography
"""
from clsQrecord import clsQrecord
from clsTimeSeries import clsTimeSeries

class clsDCfromQ(object):
    
    def __init__(self,Qpath, HydroBins, recordbreaks=[], feet=False, logbin=False):
        
        """
        Attributes:
        
        Qpath--str (Path to discharge record Excel file)
        Hydrobins--int (interval (in m^3/s) to divide discharge into)
        recordbreaks--[str,str] (optional list of strings contains start and end dates of time chunks.  If not specified, entire record will be used)
        Qrecord--clsQrecord (discharge record object; gets filled automatically)
        feet--bool (Optional tag that tells the code which unit the discharge data is in.  Default is False (m^3/s), True is ft^3/s.)        
        """
        self.Qpath = Qpath
        self.HydroBins = HydroBins
        self.recordbreaks = recordbreaks
        self.feet = feet
        self.logbin = logbin
        self.Qrecord = clsQrecord(self.Qpath,self.recordbreaks)
        self.LoadQ()
        
    def LoadQ(self):
        self.Qrecord.Load_Q_data()
        self.Qrecord.Extract_timechunks()
        """
        Get list of timechunks somewhere around here.  Need a way to individually load them to MAST.
        Chronology extension? (Spreadsheet of sediment supply, Q)
        """        
        
#        self.Qrecord.PlotDurationCurve(self.HydroBins,self.feet)
        
    def ExtractDC(self,recordperiod=[]):
        
        timechunk = []
        if len(recordperiod) > 0:
            TSindex = self.Qrecord.indexdict[recordperiod[0]+recordperiod[1]]
            timechunk = self.Qrecord.timechunks[TSindex]
        else:
            timechunk = clsTimeSeries([], self.Qrecord.Values)
        binQ, DC = timechunk.CreateDurationCurve(self.HydroBins, self.feet, self.logbin)
        
        return binQ, DC