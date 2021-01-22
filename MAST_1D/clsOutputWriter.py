# -*- coding: utf-8 -*-
"""
Created on Tue May 10 16:37:05 2016

@author: geography

Stores functions for writing output to files.
"""

import datetime
from clsNode import clsNode
import os
from clsOutputSpecs import clsOutputSpecs
import json


class clsOutputWriter(object):
    
    def __init__(self, Outputfolder, DailyNodes, startdate = ""):
        
        """
        Attributes:
        
        -Outputfolder--str (Name of folder (within parent folder) to write outputfiles)
        -DailyNodes--[int] (List of integers denoting the nodes for which daily data will be saved)
        """
        
        self.Outputfolder = Outputfolder 
        self.DailyNodes = {}
        
        for Node in DailyNodes:
            DailyNode = clsOutputSpecs(startdate)
            self.DailyNodes[Node] = DailyNode
        
    def Output(self, Filename, Reach, VariableName, Nprint, T):
        """
        Arguments:
            Filename -- str
            Reach -- clsReach
            VariableName -- str
            Nprint -- int
            T -- float
            toplot -- list # Katie add
        """
     
        outputpath = os.path.join(os.pardir, self.Outputfolder, Filename)
        # VariableName is a string that describes the reach property to be printed.
        # For example, if the property reach.node(i).etabav is to be printed, 
        # VariableName should be "etabav."
        if Nprint == 0 or os.path.isfile(outputpath) == False:
            with open(outputpath, 'w') as f:
                f.write('{:10}{:10}'.format('', 'Time(yr)') + '\n')
                f.write('{:10}{:<10.4}'.format('', T) + '\n')
                f.write('{:10}{:10}'.format('xc (m)', VariableName[:10]) + '\n')
                for Node in Reach.Node:
                    NextAttr = Node
                    for Name in VariableName.split('.'):
                        if len(Name.split('[')) == 2:
                            index = int(Name.split('[')[1][:-1])
                            NextAttr = getattr(NextAttr, Name.split('[')[0])[index]
                        else:
                            NextAttr = getattr(NextAttr, Name)
                    f.write('{:<10}{:<10.4}'.format(int(Node.xc), NextAttr) + ' ' + '\n')
        else:
            with open(outputpath, 'r+') as f:
                lines = f.readlines()
                lines[0] = lines[0][:-1] + '{:10}'.format('Time(yr)') + '\n'
                lines[1] = lines[1][:-1] + '{:<10.4}'.format(T) + '\n'
                lines[2] = lines[2][:-1] + '{:10}'.format(VariableName[:10]) + '\n'
                for i in range(Reach.nnodes()):
                    NextAttr = Reach.Node[i]
                    for Name in VariableName.split('.'):
                        if len(Name.split('[')) == 2:
                            index = int(Name.split('[')[1][:-1])
                            NextAttr = getattr(NextAttr, Name.split('[')[0])[index]
                        else:
                            NextAttr = getattr(NextAttr, Name)
                    lines[3 + i] = lines[3 + i][:-1] + '{:<10.4}'.format(NextAttr)\
                        + ' ' + '\n'                  
                    
                f.seek(0)
                f.writelines(lines)
                
    def OutputFlux(self, Reach, Nprint, T, tag):
        """
        Katie add!
        
        Outputs size-specific flux leaving bottom-most node or the feed into the top node.
        """
        name = ""
        outlist = ""        
        
        if tag == 'QsOut':
            name = 'QsOut'
            outlist = Reach.CumulativeOutput
        elif tag == 'QsIn':            
            name = 'QsIn'
            outlist = Reach.CumulativeFeed
        elif tag == 'BankIn':            
            name = 'BankIn'
            outlist = getattr(Reach, 'CumulativeBankSupply')
        elif tag == 'BankOut':            
            name = 'BankOut'
            outlist = getattr(Reach, 'CumulativeBankSink')
            
        outputpath = os.path.join(os.pardir, self.Outputfolder, name)
        if Nprint == 0 or os.path.isfile(outputpath) == False:
            with open(outputpath, 'w') as f:
                f.write('{:10}{:10}'.format('', 'Time(yr)') + '\n')
                f.write('{:10}{:<10}'.format('', T) + '\n')
                f.write('{:10}{:10}'.format('D (mm)', name) + '\n')
                
                i = 0
                for D in Reach.Node[-3].Load.GSDBedloadAv.D:
                    f.write('{:<10.5}{:<10.4}'.format(str(D), outlist[i]) + ' ' + '\n')
                    i = i + 1
                    
        else:
            with open(outputpath, 'r+') as f:
                lines = f.readlines()
                lines[0] = lines[0][:-1] + '{:10}'.format('Time(yr)') + '\n'
                lines[1] = lines[1][:-1] + '{:<10}'.format(T) + '\n'
                lines[2] = lines[2][:-1] + '{:10}'.format(name) + '\n'
                
                i = 0
                for D in Reach.Node[-1].Load.GSDBedloadAv.D:
                    lines[3+i] = lines[3+i][:-1] + '{:<10.4}'.format(outlist[i])\
                        + ' ' + '\n' 
                    i = i + 1
                        
                f.seek(0)
                f.writelines(lines) 
                
    def PopulateDailyLists(self, Reach):
        
        """
        Attributes:
        
        -Reach--clsReach (Reach with data to output)
        """
        
        for DailyNode in self.DailyNodes.keys():
            self.DailyNodes[DailyNode].PopulateLists(Reach.Node[DailyNode])
            
    def WriteDailyFiles(self):
        """
        Writes the daily data stored in the clsOutputSpecs object to a json file
        at the end of the run.  For hydrograph runs.
        """
        
        outputpath = os.path.join(os.pardir, self.Outputfolder)
        datelist = []
        for DailyNode in self.DailyNodes.keys():
            datelist = list(map(lambda x: x.strftime('%Y,%m,%d'), self.DailyNodes[DailyNode].Date))
            json.dump(self.DailyNodes[DailyNode].Q, open(os.path.join(outputpath, "save.DailyQ" + str(DailyNode)), "w", encoding="utf8" ))
            json.dump(self.DailyNodes[DailyNode].QsavBedTot, open(os.path.join(outputpath, "save.DailyQsavBedTot" + str(DailyNode)), "w", encoding="utf8" ))
            json.dump(self.DailyNodes[DailyNode].QsavTotAllFeed, open(os.path.join(outputpath, "save.DailyQsavTotAllFeed" + str(DailyNode)), "w" )) 
            json.dump(self.DailyNodes[DailyNode].Qsk, open(os.path.join(outputpath, "save.DailyQsk" + str(DailyNode)), "w" )) 
            json.dump(self.DailyNodes[DailyNode].F, open(os.path.join(outputpath, "save.DailyF" + str(DailyNode)), "w" )) 
            json.dump(self.DailyNodes[DailyNode].FpF, open(os.path.join(outputpath, "save.DailyFpF" + str(DailyNode)), "w" )) 
            json.dump(self.DailyNodes[DailyNode].SubF, open(os.path.join(outputpath, "save.DailySubF" + str(DailyNode)), "w" )) 
            json.dump(self.DailyNodes[DailyNode].Bc, open(os.path.join(outputpath, "save.DailyBc" + str(DailyNode)), "w" )) 
            json.dump(self.DailyNodes[DailyNode].CumuWiden, open(os.path.join(outputpath, "save.DailyCumuWiden" + str(DailyNode)), "w" )) 
            json.dump(self.DailyNodes[DailyNode].CumuNarrow, open(os.path.join(outputpath, "save.DailyCumuNarrow" + str(DailyNode)), "w" ))
            json.dump(self.DailyNodes[DailyNode].Eshear, open(os.path.join(outputpath, "save.DailyEshear" + str(DailyNode)), "w" )) 
            json.dump(self.DailyNodes[DailyNode].InVChange, open(os.path.join(outputpath, "save.DailyInVChange" + str(DailyNode)), "w" ))
            json.dump(self.DailyNodes[DailyNode].OutVChange, open(os.path.join(outputpath, "save.DailyOutVChange" + str(DailyNode)), "w" ))
            json.dump(self.DailyNodes[DailyNode].SinkLoadSed, open(os.path.join(outputpath, "save.DailySinkLoadSed" + str(DailyNode)), "w" ))
            
            
      
        #json.dump(datelist, open(os.path.join(outputpath, "save.DailyDate"), "wb" ))
        json.dump(datelist, open(os.path.join(outputpath, "save.DailyDate"), "w", encoding="utf8" ))
        #json.dump(datelist, open(os.path.join(outputpath, "save.DailyDate"), "w"))


      
     
      
      
      
      
      
      
      
