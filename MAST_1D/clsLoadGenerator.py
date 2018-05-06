# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 13:08:13 2016

@author: geography
"""

import math

class clsLoadGenerator(object):
    
    def __init__(self,C,tau):
        self.C = C
        self.tau = tau
        
    def create_feed(self,predamyears,damyears,removalyears,dtvals, dtcount):
        
        cumucount = 0 # in years
        counter = 0
        dtcounter = 0
        removalcumu = 0
        LoadFactor = [1,1]
        LoadFactorCount = [0,0]
        

        while cumucount <= (predamyears + damyears + removalyears):
            
            if counter == dtcount[dtcounter]:
                if dtcounter < len(dtcount) - 1: dtcounter += 1
                
            dt = dtvals[dtcounter]
            cumucount = cumucount+dt
                
            # Calculate timesteps for each feed class (pre-dam, dam, start of removal)
            if cumucount <= predamyears:
                LoadFactor[0] = 1.
                LoadFactorCount[0]=0
                LoadFactorCount[1]=counter
                
            if predamyears < cumucount <= predamyears + damyears:
                LoadFactor[1]=0.
                
                
            if predamyears + damyears < cumucount <= predamyears + damyears + removalyears:
                LoadFactorCount.append(counter)
                dtdays = dt*365.25
                removalcumu = removalcumu+dtdays
                
                Load = 1 + self.C*math.exp(-removalcumu/self.tau)
                LoadFactor.append(Load)
            print counter
            print LoadFactor[-1]
            counter = counter + 1
            
        return counter, LoadFactor, LoadFactorCount