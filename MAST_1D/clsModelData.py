# -*- coding: utf-8 -*-
"""
Created on Wed May 25 13:42:31 2016

Makes numpy arrays of 'regular' model outputs.
Loads, sorts, calculates, and stores daily data from the MAST-1D hydrograph output.

PDO-based timechunks and log/full discharge data are hard-coded.

@author: geography
"""

import sys
sys.path.append("..")
from MAST_1D.clsOutputSpecs import clsOutputSpecs
from Sediment.clsEffectiveT import clsEffectiveT
from Sediment.clsFieldGSD import clsFieldGSD
import datetime
import json
import os
import numpy as np
from copy import deepcopy
import math

class clsModelData(object):
    
    def __init__(self, folder, gscategories = []):
        
        """
        Attributes:
        
        folder--str (folder containing run output)
        gscategories--[float] (list with grainsize bins)
        output--clsOutputSpecs (object with all of the daily output data)
        Capacity--clsEffectiveT (object with sediment capacity (non-normalized for flow duration))
        CapacityRC--[float] (total discharge-specific sediment capacity)
        D50 [float] (Median grainsize normalized to get rid of mud fraction)
        Q--{str:[float]} (discharge lists for all three periods)
        Qs--{str:[float]} (total sediment load lists for all three periods)
        A--{str:[float]} (armor ratio lists for all three periods)
        Qsk--{str:[float]} (size-specific sediment load nested lists for all three periods)
        EffQ--{str:clsEffectiveT} (object with log-binned duration-averaged sediment loads for all three periods)
        Eff50Q--{str:clsEffectiveT} (object with non-binned duration-averaged sediment loads for all three periods)
        Variables--{str:[[float],[float],nparray]} (dictionary with 'regular' outputs)        
        """
        
        self.folder = folder
        self.gscategories = gscategories
        self.output = clsOutputSpecs()
        self.Capacity = clsEffectiveT(folder)
        self.CapacityRC = []
        self.D50 = []
        self.GSD = {}
        self.Mobility = {}
        self.Date = {}
        self.WaterYear = {}
        self.Q = {}
        self.Qs = {}
        self.A = {}
        self.Qsk = {}
        self.EffQ = {}
        self.Eff50Q = {}
        self.Variables = {}
        
        
    def ReloadSpecs(self, Node):
        
        """
        Takes json files saved from the model run and reloads them into the object
        """
        
        Datepath = os.path.join(self.folder,'save.DailyDate')
        Qpath = os.path.join(self.folder,'save.DailyQ' + str(Node))
        QsTotpath = os.path.join(self.folder,'save.DailyQsTot'+ str(Node))
        Qskpath = os.path.join(self.folder,'save.DailyQsk'+ str(Node))
        SubGSDpath = os.path.join(self.folder,'save.DailySubGSD'+ str(Node))
        GSDpath = os.path.join(self.folder,'save.DailyGSD'+ str(Node))
        
        self.output.Date = map(lambda x: datetime.datetime.strptime(x, '%Y,%m,%d').date(), json.load(open(Datepath)))
        self.output.Q = json.load(open(Qpath))
        self.output.QsTot = json.load(open(QsTotpath))
        self.output.Qsk = json.load(open(Qskpath))
        self.output.SubGSD = json.load(open(SubGSDpath))
        self.output.GSD = json.load(open(GSDpath))
        
    def LoadCapacity(self):
        
        """
        Fills the Capacity object from the FeedRC file.
        """
        
        self.Capacity.SplitByIndividualGSDs(cumulativeSize=False)
        
        if len(self.gscategories) > 0:
            self.Capacity.categorizeGS(self.gscategories, cumulativeSize=False)
        
        DurCapacityRC = []
        
        j = 0
        while j < len(self.Capacity.QsMatrix[:,0]):
            CapacityQsj = sum(self.Capacity.QsMatrix[j,:])
            DurCapacityRC.append(CapacityQsj)
            j = j + 1
        
        DurCapacityRC = np.array(DurCapacityRC)
        p = np.array(self.Capacity.plist)*.01
        
        # Get rid of this eventually and just use sum of QsMatrix.        
        self.CapacityRC = (DurCapacityRC/p)/365.25 # Convert from duration-averaged m^3/yr to non-duration m^3/day
        
        i = 0        
        while i < len(self.Capacity.QsMatrix[0,:]):
            self.Capacity.QsMatrix[:,i] = self.Capacity.QsMatrix[:,i]/p
            i = i + 1
            
        self.Capacity.QsMatrix = self.Capacity.QsMatrix/365.25
        
    def PrepareLists(self):
        
        """
        Creates empty lists and objects to be filled by the data. Discharges are hard-coded.
        """
        
        # All (unlogarithmic) discharge values; too clunky to find them a better way.
        fullQ = ['1.416', '4.248', '7.079', '9.911', '12.74', '15.57', '18.41', '21.24', '24.07', '26.9', '29.73', '32.56', '35.4', '38.23', '41.06', '43.89', '46.72', '49.55', '52.39', '55.22', '58.05', '60.88', '63.71', '66.54', '69.38', '72.21', '75.04', '77.87', '80.7', '83.53', '86.37', '89.2', '92.03', '94.86', '97.69', '100.5', '103.4', '106.2', '109.0', '111.9', '114.7', '117.5', '120.3', '123.2', '126.0', '128.8', '131.7', '134.5', '137.3', '140.2', '143.0', '145.8', '148.7', '151.5', '154.3', '157.2', '160.0', '162.8', '165.7', '168.5', '171.3', '174.1', '177.0', '179.8', '182.6', '185.5', '188.3', '191.1', '194.0', '196.8', '199.6', '202.5', '205.3', '208.1', '211.0', '213.8', '216.6', '219.5', '222.3', '225.1', '228.0', '230.8', '233.6', '236.4', '239.3', '242.1', '244.9', '247.8', '250.6', '253.4', '256.3', '259.1', '261.9', '264.8', '267.6', '270.4', '273.3', '276.1', '278.9', '281.8', '284.6', '287.4', '290.2', '293.1', '295.9', '298.7', '301.6', '304.4', '307.2', '310.1', '312.9', '315.7', '318.6', '321.4', '324.2', '327.1', '329.9', '332.7', '335.6', '338.4', '341.2', '344.0', '346.9', '349.7', '352.5', '355.4', '358.2', '361.0', '363.9', '366.7', '369.5', '372.4', '375.2', '378.0', '380.9', '383.7', '386.5', '389.4', '392.2', '395.0', '397.9', '400.7', '403.5', '406.3', '409.2', '412.0', '414.8', '417.7', '420.5', '423.3', '426.2', '429.0', '431.8', '434.7', '437.5', '440.3', '443.2', '446.0', '448.8', '451.7', '454.5', '457.3', '460.1', '463.0', '465.8', '468.6', '471.5', '474.3', '477.1', '480.0', '482.8', '485.6', '488.5', '491.3', '494.1', '497.0', '499.8', '502.6', '505.5', '508.3', '511.1', '514.0', '516.8', '519.6', '522.4', '525.3', '528.1', '530.9', '533.8', '536.6', '539.4', '542.3', '545.1', '547.9', '550.8', '553.6', '556.4', '559.3', '562.1', '564.9', '567.8', '570.6', '573.4', '576.2', '579.1', '581.9', '584.7', '587.6', '590.4', '593.2', '596.1', '598.9', '601.7', '604.6']
        fullQ = map(lambda x: float(x), fullQ)
        #Dbdy = [1.0, 2.0,  4.0,  8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0]
        Dbdy = [1.0, 2.0, 2.8284271247461903, 4.0, 5.656854249492381, 8.0, 11.313708498984761, 16.0, 22.627416997969522, 32.0, 45.254833995939045, 64.0, 90.50966799187809, 128.0, 181.01933598375618, 256.0, 362.03867196751236, 512.0,1024.]
        
        self.Q['Period 1'] = []
        self.Q['Period 2'] = []
        self.Q['Period 3'] = []

        self.Qs['Period 1'] = []
        self.Qs['Period 2'] = []
        self.Qs['Period 3'] = []
        
        self.Mobility['Period 1'] = []
        self.Mobility['Period 2'] = []
        self.Mobility['Period 3'] = []
        
        self.A['Period 1'] = []
        self.A['Period 2'] = []
        self.A['Period 3'] = []
        
        self.Qsk['Period 1'] = []
        self.Qsk['Period 2'] = []
        self.Qsk['Period 3'] = []
        
        self.GSD['Period 1'] = []
        self.GSD['Period 2'] = []
        self.GSD['Period 3'] = []
        
        self.Date['Period 1'] = []
        self.Date['Period 2'] = []
        self.Date['Period 3'] = []
        
        self.WaterYear['Period 1'] = []
        self.WaterYear['Period 2'] = []
        self.WaterYear['Period 3'] = []
        
        # Fill EffectiveT objects with empty lists
        periods = ['Period 1', 'Period 2', 'Period 3']
        
        for period in periods:
            
            EffObject = clsEffectiveT('filler')
            #EffObject.QsMatrix = np.zeros((len(self.Capacity.Qlist), len(self.gscategories)))
            #EffObject.p = np.zeros(len(self.Capacity.Qlist))
            #EffObject.Qlist = deepcopy(self.Capacity.Qlist)
            #EffObject.QsMatrix = np.zeros((len(fullQ), len(self.gscategories)))
            EffObject.QsMatrix = np.zeros((len(fullQ), len(Dbdy)))
            EffObject.p = np.zeros(len(fullQ))
            EffObject.Qlist = deepcopy(fullQ)
            #EffObject.Dlist = self.gscategories
            EffObject.Dlist = Dbdy
            
            self.EffQ[period] = deepcopy(EffObject)
        
            Eff50Object = clsEffectiveT('filler') # will need to make a flow duration curve.
            #Eff50Object.QsMatrix = np.zeros((len(fullQ), len(self.gscategories)))
            Eff50Object.QsMatrix = np.zeros((len(fullQ), len(Dbdy)))
            Eff50Object.p = np.zeros(len(fullQ))
            Eff50Object.Qlist = fullQ
            #Eff50Object.Dlist = self.gscategories
            Eff50Object.Dlist = Dbdy
            
            self.Eff50Q[period] = deepcopy(Eff50Object)
         
    def find_nearest(self, array,value):
        
        """
        Find discharge bin for Q         
        """
        
        idx = np.searchsorted(array, value, side="left")
        if idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx]):
            return idx-1
        else:
            return idx
            
    def DivideData(self):
        
        """
        Loops through the model output and fills lists according to time period.  
        Grainsize, discharge, and time periods are hardcoded.
        """
        
        Qbins = deepcopy(map(lambda x: float(x), self.Capacity.Qlist))
        #Dbdy = [0.251, 1.414, 2.828, 5.657, 11.31, 22.63, 45.25, 90.51, 181.0, 362.0]
        Dbdy = [1.0, 2.0, 2.8284271247461903, 4.0, 5.656854249492381, 8.0, 11.313708498984761, 16.0, 22.627416997969522, 32.0, 45.254833995939045, 64.0, 90.50966799187809, 128.0, 181.01933598375618, 256.0, 362.03867196751236, 512.0,1024.]
        fullQ = ['1.416', '4.248', '7.079', '9.911', '12.74', '15.57', '18.41', '21.24', '24.07', '26.9', '29.73', '32.56', '35.4', '38.23', '41.06', '43.89', '46.72', '49.55', '52.39', '55.22', '58.05', '60.88', '63.71', '66.54', '69.38', '72.21', '75.04', '77.87', '80.7', '83.53', '86.37', '89.2', '92.03', '94.86', '97.69', '100.5', '103.4', '106.2', '109.0', '111.9', '114.7', '117.5', '120.3', '123.2', '126.0', '128.8', '131.7', '134.5', '137.3', '140.2', '143.0', '145.8', '148.7', '151.5', '154.3', '157.2', '160.0', '162.8', '165.7', '168.5', '171.3', '174.1', '177.0', '179.8', '182.6', '185.5', '188.3', '191.1', '194.0', '196.8', '199.6', '202.5', '205.3', '208.1', '211.0', '213.8', '216.6', '219.5', '222.3', '225.1', '228.0', '230.8', '233.6', '236.4', '239.3', '242.1', '244.9', '247.8', '250.6', '253.4', '256.3', '259.1', '261.9', '264.8', '267.6', '270.4', '273.3', '276.1', '278.9', '281.8', '284.6', '287.4', '290.2', '293.1', '295.9', '298.7', '301.6', '304.4', '307.2', '310.1', '312.9', '315.7', '318.6', '321.4', '324.2', '327.1', '329.9', '332.7', '335.6', '338.4', '341.2', '344.0', '346.9', '349.7', '352.5', '355.4', '358.2', '361.0', '363.9', '366.7', '369.5', '372.4', '375.2', '378.0', '380.9', '383.7', '386.5', '389.4', '392.2', '395.0', '397.9', '400.7', '403.5', '406.3', '409.2', '412.0', '414.8', '417.7', '420.5', '423.3', '426.2', '429.0', '431.8', '434.7', '437.5', '440.3', '443.2', '446.0', '448.8', '451.7', '454.5', '457.3', '460.1', '463.0', '465.8', '468.6', '471.5', '474.3', '477.1', '480.0', '482.8', '485.6', '488.5', '491.3', '494.1', '497.0', '499.8', '502.6', '505.5', '508.3', '511.1', '514.0', '516.8', '519.6', '522.4', '525.3', '528.1', '530.9', '533.8', '536.6', '539.4', '542.3', '545.1', '547.9', '550.8', '553.6', '556.4', '559.3', '562.1', '564.9', '567.8', '570.6', '573.4', '576.2', '579.1', '581.9', '584.7', '587.6', '590.4', '593.2', '596.1', '598.9', '601.7', '604.6']
        fullQ = map(lambda x: float(x), fullQ)
        
        i = 0
        for val in self.output.Date:            
            ActiveGSD = clsFieldGSD()
            SubstrateGSD = clsFieldGSD()
            SizeLoad = clsEffectiveT('filler')
            
            j = 0

            for x in Dbdy:
                ActiveGSD.totalmass[math.log(Dbdy[j], 2)] = self.output.GSD[i][j+1]
                SubstrateGSD.totalmass[math.log(Dbdy[j], 2)] = self.output.SubGSD[i][j+1]
                j = j + 1

            ActiveGSD.TotalBulk = sum(ActiveGSD.totalmass.values())
            SubstrateGSD.TotalBulk = sum(SubstrateGSD.totalmass.values())
            ActiveGSD.FractionMaker()
            ActiveGSD.calculate_stats()
            SubstrateGSD.FractionMaker()
            SubstrateGSD.calculate_stats()
            
            self.GSD[i] = ActiveGSD.GSD.values()          
            self.D50.append(ActiveGSD.D50)
            
            # Create mobility list
            TotalTrans = sum(self.output.Qsk[i])
            LoadF = map(lambda x: x/TotalTrans*100, self.output.Qsk[i])
            Mobility = (np.array(LoadF[1:])/np.array(self.GSD[i])).tolist()
            
            # For size-specific rating curves
            QsVals = map(lambda x: x*60*60*24, self.output.Qsk[i]) # m^3/day
            QsVals = QsVals[1:] # Gets rid of mud fraction
            SizeLoad.QsMatrix = np.matrix(QsVals)
            SizeLoad.Dlist = Dbdy
            SizeLoad.Qlist = [self.output.Q[i]]
            
            SizeCatLoad = deepcopy(SizeLoad)
            if len(self.gscategories) > 0:
                SizeCatLoad.categorizeGS(self.gscategories, cumulativeSize=False)
            
            # Figure out which bin to put Q into
            Qindex = self.find_nearest(fullQ, self.output.Q[i])
            #Qindex = self.find_nearest(Qbins, self.output.Q[i])
            fullQindex = self.find_nearest(fullQ, self.output.Q[i])
            WY = self.GetWaterYear(val)
            
            if WY >= 1927 and WY <1947:
                self.WaterYear['Period 1'].append(WY)
                self.Date['Period 1'].append(val)
                self.Q['Period 1'].append(self.output.Q[i])
                self.Qs['Period 1'].append(self.output.QsTot[i])
                self.A['Period 1'].append(ActiveGSD.D50/SubstrateGSD.D50)
                self.Mobility['Period 1'].append(Mobility)
                QskOut = SizeLoad.QsMatrix[0].tolist()
                QskCatOut = SizeCatLoad.QsMatrix[0].tolist()
                self.Qsk['Period 1'].append(QskCatOut)
                self.EffQ['Period 1'].QsMatrix[Qindex,:] = self.EffQ['Period 1'].QsMatrix[Qindex,:] + np.array(QskOut)
                self.Eff50Q['Period 1'].QsMatrix[fullQindex,:] = self.Eff50Q['Period 1'].QsMatrix[fullQindex,:] + np.array(QskOut)
                
            elif WY >= 1948 and WY < 1977:
                self.WaterYear['Period 2'].append(WY)
                self.Date['Period 2'].append(val)
                self.Q['Period 2'].append(self.output.Q[i])
                self.Qs['Period 2'].append(self.output.QsTot[i])
                self.A['Period 2'].append(ActiveGSD.D50/SubstrateGSD.D50)
                self.Mobility['Period 2'].append(Mobility)
                QskOut = SizeLoad.QsMatrix[0].tolist()
                QskCatOut = SizeCatLoad.QsMatrix[0].tolist()
                self.Qsk['Period 2'].append(QskCatOut)
                self.EffQ['Period 2'].QsMatrix[Qindex,:] = self.EffQ['Period 2'].QsMatrix[Qindex,:] + np.array(QskOut)
                self.Eff50Q['Period 2'].QsMatrix[fullQindex,:] = self.Eff50Q['Period 2'].QsMatrix[fullQindex,:] + np.array(QskOut)
                
            elif WY >=1977:
                self.WaterYear['Period 3'].append(WY)
                self.Date['Period 3'].append(val)
                self.Q['Period 3'].append(self.output.Q[i])
                self.Qs['Period 3'].append(self.output.QsTot[i])
                self.A['Period 3'].append(ActiveGSD.D50/SubstrateGSD.D50)
                self.Mobility['Period 3'].append(Mobility)
                QskOut = SizeLoad.QsMatrix[0].tolist()
                QskCatOut = SizeCatLoad.QsMatrix[0].tolist()
                self.Qsk['Period 3'].append(QskCatOut)
                self.EffQ['Period 3'].QsMatrix[Qindex,:] = self.EffQ['Period 3'].QsMatrix[Qindex,:] + np.array(QskOut)
                self.Eff50Q['Period 3'].QsMatrix[fullQindex,:] = self.Eff50Q['Period 3'].QsMatrix[fullQindex,:] + np.array(QskOut)
                
            i = i + 1
        
        # Convert Size-specific transport from model into annual average
        self.EffQ['Period 1'].QsMatrix = self.EffQ['Period 1'].QsMatrix/(1948-1925)
        self.EffQ['Period 2'].QsMatrix = self.EffQ['Period 2'].QsMatrix/(1977-1948)
        self.EffQ['Period 3'].QsMatrix = self.EffQ['Period 3'].QsMatrix/(2011-1977)
        
        self.Eff50Q['Period 1'].QsMatrix = self.Eff50Q['Period 1'].QsMatrix/(1948-1925)
        self.Eff50Q['Period 2'].QsMatrix = self.Eff50Q['Period 2'].QsMatrix/(1977-1948)
        self.Eff50Q['Period 3'].QsMatrix = self.Eff50Q['Period 3'].QsMatrix/(2011-1977)
        
        # Put grainsizes into categories if necessary
        if len(self.gscategories) > 0: 
            
            self.EffQ['Period 1'].categorizeGS(self.gscategories, cumulativeSize=False)
            self.EffQ['Period 2'].categorizeGS(self.gscategories, cumulativeSize=False)
            self.EffQ['Period 3'].categorizeGS(self.gscategories, cumulativeSize=False)
            
            self.Eff50Q['Period 1'].categorizeGS(self.gscategories, cumulativeSize=False)
            self.Eff50Q['Period 2'].categorizeGS(self.gscategories, cumulativeSize=False)
            self.Eff50Q['Period 3'].categorizeGS(self.gscategories, cumulativeSize=False)
        
        # Make cumulative
        self.EffQ['Period 1'].MakeCumulativeSizeMatrix()
        self.EffQ['Period 2'].MakeCumulativeSizeMatrix()
        self.EffQ['Period 3'].MakeCumulativeSizeMatrix()
        
        self.Eff50Q['Period 1'].MakeCumulativeFlowMatrix()
        self.Eff50Q['Period 2'].MakeCumulativeFlowMatrix()
        self.Eff50Q['Period 3'].MakeCumulativeFlowMatrix()
        
    def LoadObject(self, Node):
        
        """
        Runs through all the functions to fill in the object
        """
        
        self.ReloadSpecs(Node)
        self.LoadCapacity()
        self.PrepareLists()
        self.DivideData()
        
    def LoadNodeOutput(self, variablelist):

        """
        Takes specified variable and loads node/time output data into a numpy matrix.
        
        Attributes:
        
        path (str)--Output path
        variablelist ([str])--List of output variables to extract (without 'Out_' prefix)
        """

        #  Import and read model output data (files with prefix 'Out_')    
        nodelist = []
        valuelist = []
        
        for variable in variablelist:
            dfile = open(os.path.join(self.folder, 'Out_' + variable) ,'r')
            dfile = dfile.readlines()
            str(dfile)
    
            time = dfile[1]
            times = time.split()
    
            del dfile[0:3]
         
            #  Forms data into numpy array
            for row in dfile:
                splitset = row.split()
                
                node = float(splitset[0])
                nodelist.append(node)
    
                values = map(lambda x: float(x), splitset[1:])                
                valuelist.append(values)                
        
            valuearray = np.array(valuelist)     
                
            self.Variables[variable] = [times, nodelist, valuearray] # for matplotlib
            
    def GetWaterYear(self, date):
        
        """
        Uses date to come up with water year
        
        Attributes:
        date (datetime)--The date
        """
        
        wateryear = 0
        month = date.month
        if month <10:
            wateryear = date.year
        else:
            wateryear = date.year + 1
            
        return wateryear
            
            
            
            
            
            
            
            
            
            
            
            
            
                    
                    
                    
                    
                    