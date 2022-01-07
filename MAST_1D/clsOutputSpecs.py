# -*- coding: utf-8 -*-
"""
Created on Tue May 10 16:37:05 2016

@author: geography

Stores daily data for output.  For use with hydrographs.
"""

import datetime
from clsNode import clsNode
import os


class clsOutputSpecs(object):
    
    def __init__(self, startdate = ""):
        
        # Daily variables        
        self.startdate = startdate
        self.Date = []
        self.Q = []
        self.QsavBedTot = []
        self.QsavTotAllFeed = []
        self.Qsk_1 = []
        self.F = []
        self.FpF = []
        self.SubF = []
        self.Bc = []
        self.CumuWiden = []
        self.CumuNarrow = []
        self.Eshear = [] 
        self.InVChange = []
        self.OutVChange = []   
        self.SinkLoadSed = []            
        
    def AddDay(self):
        if len(self.Date) == 0:
            self.Date.append(self.startdate)
        else:
            PriorDate = self.Date[-1]
            NewDate = PriorDate + datetime.timedelta(days=1)
            self.Date.append(NewDate)

    def PopulateLists(self,Node):
        
        # Adds daily variables to list--to be called at the end of each day.
        
        self.AddDay()
        self.Q.append(Node.DC.Qw[0])
        self.QsavBedTot.append(Node.Load.QsavBedTot)
        #self.QsavTotAllFeed.append(Node.Load.QsavTotAllFeed)
        x = Node.Load.QsAvkLoad[1]
        self.Qsk_1.append(x.tolist())
        #x = Node.ActiveLayer.GSD.F
        #self.F.append(x.tolist())
        #x = Node.Substrate[-1].C.GSD.F
        #self.SubF.append(x.tolist())
        #x = Node.Floodplain.GSD.F
        #self.FpF.append(x.tolist())
        x = Node.Load.Eshear
        self.Eshear.append(x.tolist())
        #x = Node.ActiveLayer.ExSed.InVerticalChange
        #self.InVChange.append(x.tolist())
        #x = Node.ActiveLayer.ExSed.OutVerticalChange
        #self.OutVChange.append(x.tolist())
        #x = Node.ActiveLayer.SinkLoadSed
        #self.SinkLoadSed.append(x.tolist())
        self.Bc.append(Node.Bc)
        #self.CumuWiden.append(Node.CumulativeWidening)
        #self.CumuNarrow.append(Node.CumulativeNarrowing)
      
      
      
      
     
      
      
      
      
      
      
      
      
      
      
