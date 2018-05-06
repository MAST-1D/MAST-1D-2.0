import numpy as np
from copy import deepcopy
import os
#from Tkinter import *

from clsReach import clsReach
from clsTracerProperties import clsTracerProperties
from clsInputs import clsInputs
from clsOutputSpecs import clsOutputSpecs
import datetime
import pickle as pickle
import pdb
import math
import numpy as np
import json



class clsModel(object):
    
    def __init__(self, inputs):
        
        self.inputs = inputs
        #self.plotresponse = StringVar()

    def SetupTracers(self, TracerProperties):
        """
        Setup as Cosmogenic 14C
        
        Arguments:
            TracerProperties -- clsTracerProperties
            
        Return: clsTracerProperties
        """
        TracerProperties.coj[0] = self.inputs.coj[0] # % production from this process 
            # (at surface, presumably--not integrated over depth)
        TracerProperties.coj[1] = self.inputs.coj[1] # % production from this process 
            # (at surface, presumably--not integrated over depth)
        TracerProperties.coj[2] = self.inputs.coj[2] # % production from this process 
            # (at surface, presumably--not integrated over depth)
        TracerProperties.Lcj[0] = self.inputs.Lcj[0] # g/cm2
        TracerProperties.Lcj[1] = self.inputs.Lcj[1] # g/cm2
        TracerProperties.Lcj[2] = self.inputs.Lcj[2] # g/cm2
        TracerProperties.Name = self.inputs.Name
        TracerProperties.DecayConst = self.inputs.DecayConst # 1/yr
        TracerProperties.ProductionRate = self.inputs.ProductionRate # Atoms/g Si02/yr
        TracerProperties.FalloutRate = self.inputs.FalloutRate # Atoms/cm2/yr
        return TracerProperties

    def SetupReach(self, TransFunc, TrinityFit, CalibrationFactor):
        """
        Arguments:
            TrinityFit -- bool
            CalibrationFactor -- float
            
        Return: clsReach
        """
        Dbdy = self.inputs.Dbdy

        Reach = clsReach(self.inputs.Nnodes, self.inputs.NLayers, 1, Dbdy, len(self.inputs.Qw), self.inputs.ManningStabilizer)        
        Reach.BoundaryHc = self.inputs.BoundaryHc
        Reach.SetBoundary = self.inputs.SetBoundary
        
        i = 0
        for Node in Reach.Node:
        
            Node.DC.Qw = self.inputs.Qw
            Node.DC.p = self.inputs.p

            Node.ActiveLayer.GSD.F = self.inputs.Fa 
            Node.ActiveLayer.GSD.UpdateStatistics()
            
            Node.Floodplain.GSD.F = self.inputs.Fp
            Node.Floodplain.GSD.UpdateStatistics()
            Node.Bc = self.inputs.Bc
            Node.Bf = self.inputs.Bf              
            Node.Cfc = self.inputs.Cfc # = 1/Cz^2 # Need to check value in Excel version
            Node.Cff = self.inputs.Cff # Need to check value in Excel version            
            Node.Slope = self.inputs.Slope                
            migration = self.inputs.migration # 0.8 m/yr
            Node.cbank = migration / 24. / 60. / 60. / 365.25 # m/s
            Node.ChSin = self.inputs.ChSin
            
            #  Only assigns a cross-section length if it is present in the inputs,
            #  or else calculates it from the reach length and number of nodes.
            if hasattr(self.inputs, 'dxf'):
                Node.dxc = self.inputs.dxf[i]*Node.ChSin
            else:
                Node.dxc = self.inputs.reachlength / Reach.nnodes()
            
            #  Katie add--assigns channel coordinate at beginning (upstream) end of node
            Node.xc = 0
            if i == 0:
                Node.xc = 0
            else:
                Node.xc = Reach.Node[i-1].xc + Reach.Node[i-1].dxc
                
            #Node.xc = Node.dxc * (i - 1)
                
            Node.FSandSusp = self.inputs.FSandSusp
            Node.MudFraction = self.inputs.MudFraction
            Node.Floodplain.L = self.inputs.FloodplainL                
            Node.ActiveLayer.L = self.inputs.ActiveLayerL
            Node.Floodplain.Volume = Node.Floodplain.L * Node.Bf * Node.dxc / \
                Node.ChSin
            Node.ActiveLayer.Volume = Node.ActiveLayer.L * Node.Bc * Node.dxc
            Node.Hpb = self.inputs.Hpb
            Node.sigma = self.inputs.sigma
            Node.lambdap = self.inputs.lambdap
            Node.Kbar = self.inputs.Kbar
            Node.AlphaBar = self.inputs.AlphaBar
            Node.AlphaPartlyAlluvial = self.inputs.AlphaPartlyAlluvial
            Node.ncAddons = self.inputs.ncAddons
            Node.ncMultiplier = self.inputs.ncMultiplier
            Node.nf = self.inputs.nf
            
            for m in range(Reach.NLayers):
                Node.Substrate[m].C.L = self.inputs.LayerL
                Node.Substrate[m].F.L = self.inputs.LayerL
                Node.Substrate[m].C.Volume = Node.Substrate[m].C.L * Node.Bc * \
                    Node.dxc
                Node.Substrate[m].F.Volume = Node.Substrate[m].F.L * Node.Bf * \
                    Node.dxc / Node.ChSin
                
                #  Substrate created automatically
                #Node.Substrate[m].C.GSD = deepcopy(Node.ActiveLayer.GSD)
                #Node.Substrate[m].F.GSD = deepcopy(Node.Floodplain.GSD)
                
                #  Katie changed--define substrate (able to handle lag deposits)
                DLag = self.inputs.DLag
                FLag = self.inputs.FLag
                totF = 1 + sum(FLag)
                
                #Node.Substrate[m].F.GSD = deepcopy(Node.Floodplain.GSD)
                Node.Substrate[m].F.GSD = deepcopy(Node.ActiveLayer.GSD) # Katie change floodplain substrate to same GSD as active layer
                Node.Substrate[m].C.GSD = deepcopy(Node.ActiveLayer.GSD)
                
                for k in range(1, Node.NSizes + 1):
                    Node.Substrate[m].F.GSD.F[k] = Node.Substrate[m].F.GSD.F[k]/totF # Katie change
                    Node.Substrate[m].C.GSD.F[k] = Node.Substrate[m].C.GSD.F[k]/totF
                    
                for D in DLag:
                    Node.Substrate[m].F.GSD.F[D]=FLag[D]
                    Node.Substrate[m].C.GSD.F[D]=FLag[D]

                Node.Substrate[m].F.GSD.UpdateStatistics()
                Node.Substrate[m].C.GSD.UpdateStatistics()
                
             #  Katie move this line from the top of the loop to here.
            Node.etabav = 0.
            if i == 0:
                Node.etabav = 33.9
            else:
                Node.etabav = Reach.Node[i-1].etabav-Reach.Node[i-1].Slope*Reach.Node[i-1].dxc
            #Node.etabav = 33.9 - Node.Slope * Node.dxc * i
            Node.InitialBedElev = Node.etabav
            i = i + 1

        # Katie add for reservoir
#        i = 1
#        for Node in Reach.Node[-6:]:
#            Node.Bc = Node.Bf #Node.Bc + 5*i
#            i = i +1
            
        # ***********************Set Up Hydraulics******************************
        # Find Normal Depth--probably can be removed
        for Node in Reach.Node:
            Node.UpdateDepthAndDischargeAtAllFlows(self.inputs.vfunc)
            for j in range(Reach.NFlows):
                Node.DC.WSE[j] = Node.DC.Hc[j] + Node.etabav

        # Find downstream boundary condition--Katie commented out for normal flow            
        Reach.Node[-1].DC.WSE = self.find_downstream_boundary(Reach, Reach.Node[-1], Reach.Node[-1].DC.Qw, setB = Reach.BoundaryHc)                
        # Compute backwater hydraulics
        Reach.UpdateManningDepthAtAllFlowAndNodes()

        # *****************************END HYDRAULICS***************************
        
        # *****************************SET FEED*********************************
        for Node in Reach.Node:
            Node.Load.UpdateSedimentLoadByDurationAndSize(Node.DC, \
                Node.ActiveLayer.GSD, 1000., 2.7, Node.Bc, \
                Node.FractionAlluvial, TransFunc, TrinityFit, CalibrationFactor)

        # Assume mud feed is 10 times sand feed--Katie:  now there is a variable for this.
        Node = Reach.Node[0]
        Node.Load.QsAvkFeed[0] = 0.
        for j in range(Node.DC.NFlows()):
            Node.Load.Qsjkfeed[j, 0] = Node.MudFraction * Node.Load.Qsjk[j, 1] # this may not be all sand sizes
            Node.Load.QsAvkFeed[0] = Node.Load.QsAvkFeed[0] + \
                Node.Load.Qsjkfeed[j, 0] * Node.DC.p[j]
        # Assume feed is 100% of capacity for all bed material
        for k in range(1, Node.NSizes + 1):
            Node.Load.QsAvkFeed[k] = 1. * Node.Load.QsAvkLoad[k]
            for j in range(Node.DC.NFlows()):
                Node.Load.Qsjkfeed[j, k] = 1. * Node.Load.Qsjk[j, k]

        #*********************Katie add--rating curve for shear stress********
        with open(os.path.join(os.pardir, self.inputs.Outputfolder, "TauRC"), 'w') as f:
            f.write('{:10}{:10}'.format('Q', 'Tau (rho*g*Hc*Sf)') + '\n')
            for k in range(Node.DC.NFlows()):
                shear = 1000*9.81*Node.DC.Hc[k]*Node.DC.Sf[k]
                f.write('{:<10.4}{:<10.4}'.format(Node.DC.Qw[k],shear) + '\n')
            f.close()
   
        # ********************Output gravel rating curve for feed*************
#        with open('Output' + '//' + 'FeedRC', 'w') as f:  
        with open(os.path.join(os.pardir, self.inputs.Outputfolder, "FeedRC"), 'w') as f:
            f.write('{:10}{:10}{}'.format('', '', 'Particle Diameter (mm)') + '\n'\
                + '{:10}{:10}'.format('', ''))
            for k in range(1, Node.NSizes + 1):
                f.write('{:<10.4}'.format(Node.Load.GSDBedloadAv.D[k]))
            f.write('\n\n' + '{:10}'.format(''))
    
            f.write('{:10}{}'.format('Q (cms)', 'Load (m3/yr)') + '\n')
            for j in range(Reach.NFlows):
                f.write('{:10}{:<10.4}'.format('', Node.DC.Qw[j]))
                for k in range(1, Node.NSizes + 1):
                    f.write('{:<10.4}'.format(Node.Load.Qsjkfeed[j, k] * 60. * 60.\
                        * 24. * 365.25))
                f.write('\n')
            
            f.write('\n')
            f.write('{:10}{:10}{}'.format('probabilit', 'Q(cms)', 'Load * p[j]' + \
                'm3/yr') + '\n')
            for j in range(Reach.NFlows):
                f.write('{:<10.4}{:<10.4}'.format(Node.DC.p[j], Node.DC.Qw[j]))
                for k in range(1, Node.NSizes + 1):
                    f.write('{:<10.4}'.format(Node.Load.Qsjkfeed[j, k] * \
                        Node.DC.p[j] * 60. * 60. * 24. * 365.25))
                f.write('\n')
            f.write('{:10}{:10}'.format('', 'Sum'))
            for k in range(1, Node.NSizes + 1):
                f.write('{:<10.4}'.format(Node.Load.QsAvkFeed[k] * 60. * 60. * 24.\
                    * 365.25))
    
        # ******************End Output Rating Curve for Feed***********************
        
        # *****************************END FEED************************************
        
        # *************************MAKE SURE SYSTEM IS AT GEOMORPHIC EQUILIBRIUM***
        # Compute and set the floodplain number required to achieve the correct 
        # bank height at equilibrium conditions
        Node = Reach.Node[0]
        Fl = Node.EquilibriumMudFloodplainNumber(Node.Floodplain.L - Node.Hpb - \
            Node.ActiveLayer.L, Node.cbank)
        for Node in Reach.Node:
            Node.Flmud = Fl
        
        # Determine migration factor based on mud deposition rate at upstream node.
        # This ensures that true equilibrium is possible.
        # Otherwise, migration at low flow exports most of the mud from the 
        # floodplain.
        
        # Katie note:  this is not accounted for in the hydrograph runs--does not matter if width change is on.
        Migrationfactor = [0.] * Reach.NFlows
        Node = Reach.Node[0]
        if Node.cbank > 0.:
            Node.UpdateDMudj()
            for j in range(Reach.NFlows):
                Migrationfactor[j] = Node.Dfjk[j, 0] * Node.DC.p[j] / Node.Dfav[0]
            for i in range(Reach.nnodes()):
                for j in range(Reach.NFlows):
                    Reach.Node[i].DC.Migrationfactor[j] = Migrationfactor[j]
    
        # Compute the size distribution of the floodplain at equilibrium
        FFpEquilib = [0.] * (Reach.NBedSizes + 1)
        FFpEquilib[0] = (Node.ActiveLayer.GSD.F[0] * Node.ActiveLayer.L + \
            Node.FkPointBarDepositAnnualAverage[0] * Node.Hpb + \
            (Node.Floodplain.L - Node.ActiveLayer.L - Node.Hpb)) / \
            (Node.Floodplain.L)

        for k in range(1, Reach.NBedSizes + 1):
            FFpEquilib[k] = (Node.ActiveLayer.GSD.F[k] * Node.ActiveLayer.L + \
                Node.FkPointBarDepositAnnualAverage[k] * Node.Hpb) / \
                Node.Floodplain.L
        for i in range(Reach.nnodes()):
            for k in range(Reach.NBedSizes + 1):
                Reach.Node[i].Floodplain.GSD.F[k] = FFpEquilib[k]
        
        for i in range(Reach.nnodes()):
            for m in range(Reach.NLayers):
                Reach.Node[i].Substrate[m].F.GSD = \
                    deepcopy(Reach.Node[i].Floodplain.GSD) # Katie change from Floodplain to ActiveLayer--then changed back.
                    
        # *********************************END GEOMORPHIC EQUILIBRIUM**************
                    
        
        # ***********************************SETUP TRACERS*************************
        Node = Reach.Node[0]
        for k in range(Node.NSizes + 1):
            if k == 0:
                # Set tracer concentration in mud feed to 1
                for j in range(Reach.NFlows):
                    Node.Load.TMudFeedj[j, 0] = 1.
            else:
                # Set tracer concentration in bed material feed to 1
                Node.Load.TBedFeedk[k, 0] = 1.
            
            # Set tracer concentration in all reservoirs to arbitrary value
            ArbitraryValue = 1.
            for i in range(Reach.nnodes()):
                Reach.Node[i].Floodplain.T[k, 0] = ArbitraryValue
                Reach.Node[i].ActiveLayer.T[k, 0] = ArbitraryValue
                for m in range(Reach.NLayers):
                    Reach.Node[i].Substrate[m].C.T[k, 0] = ArbitraryValue
                    Reach.Node[i].Substrate[m].F.T[k, 0] = ArbitraryValue
        # ********************************END SETUP TRACERS ***********************

        return Reach

    def find_downstream_boundary(self, Reach, Node, Qw, setB = 0):
        """
        Finds the downstream boundary condition given a Node object and list and
        number of discharges.
        """
        WSE = np.zeros(Reach.NFlows)
        # Find downstream boundary condition--Katie commented out for normal flow
        for j in range(Reach.NFlows):
            # this is a poor way to set the downstream boundary condition, but it
            # cannot be solved directly using the normal depth solution implemented
            # in the node class since the definition of friction slope is slightly
            # different.
            # This simply iterates through depth until water surface slope is 
            # parallel to the bed.
            
            depth = 0. # Katie add this and following if statement to allow user to set water surface elevation
            if Reach.SetBoundary == True:
                depth = setB
                
            else:
                depth = Reach.ManningNormalDepthForBoundary(Node.Slope, Node.Bc, \
                    Node.Bf, Node.nc(), Node.nf, Node.Floodplain.L - \
                    Node.ActiveLayer.L, Qw[j], Node.ChSin)
            
            WSE[j] = depth + Node.etabav
        return WSE

    def set_up_hydrograph(self, Reach):
        """
        After the floodplain initial conditions are set up, the flow and transport
        bins are reset to accomodate just one flow.
        """
        TyearOld = 0
        Reach.NFlows = 1
        for Node in Reach.Node:
            Node.Flmud = .2 # Here is where Katie hard-codes the floodplain number.
            Node.NFlows = 1                        
            Node.DC.Qw = np.zeros(Node.NFlows)
            #Node.DC.Qw[0]=self.inputs.Qlist[0]
            Node.DC.p = np.zeros(Node.NFlows)
            Node.DC.p[0] = 1.
            Node.DC.Uc = np.zeros(Node.NFlows)
            Node.DC.Uf = np.zeros(Node.NFlows)
            Node.DC.Hc = np.zeros(Node.NFlows)
            Node.DC.Hf = np.zeros(Node.NFlows)
            Node.DC.Qs = np.zeros(Node.NFlows)
            Node.DC.Qwc = np.zeros(Node.NFlows)
            Node.DC.Qwf = np.zeros(Node.NFlows)
            Node.DC.WSE = np.zeros(Node.NFlows)
            Node.Load.Qsjk = np.zeros((Node.NFlows, Node.NSizes + 1))
            Node.Load.Qsjkfeed = np.zeros((Node.NFlows, Node.NSizes + 1))
            Node.Load.QsjTot = np.zeros(Node.NFlows)
            #Node.DC.WSE[0] = Node.etabav
            Node.DC.Sf = np.zeros(Node.NFlows) 
            
    def CustomNodes(self, Reach):
        """
        Inputs variables for individual Nodes (to be run after initial floodplain
        conditions and Control Nodes are set). More custom variables can be added here.
        """
        
        for i in range(Reach.nnodes()):

            Node = Reach.Node[i]
            
            #  Clear some attributes from original node
            Node.Bcrate = 0.
            Node.Floodplain.Volume = 0.
            
            if hasattr(self.inputs, 'BcNodes'):
                Node.Bc = self.inputs.BcNodes[i]
            if hasattr(self.inputs, 'BfNodes'):
                Node.Bf = self.inputs.BfNodes[i]
                
            # *****************Katie add:  SET PARTLY ALLUVIAL REACHES*****************
        #  Partly-alluvial reaches are set after the floodplain number calculation
        #  is done so that the migration rate can be set to zero.  In theory, a 
        #  partly-alluvial node with a migrating floodplain can exist.  Here, however,
        #  for the purpose of the Elwha River, I defined 'Canyon' nodes in the input
        #  class, and these nodes will be set to have migration rates of zero, no
        #  width change, and very high floodplain heights to simulate a canyon.
               
            if hasattr(self.inputs, 'Canyon'):            
                if self.inputs.Canyon[i] == True:
                    Node.Canyon = True
                    Node.cbank = 0.
                    Node.PartlyAlluvial = True
                    Node.Floodplain.L = 4. # Arbitrarily set to be very high.
                    
            Node.Floodplain.Volume = Node.Floodplain.L * Node.Bf * Node.dxc / \
                Node.ChSin
        # *****************END SET PARTLY ALLUVIAL REACHES*****************
        

  
    def main(self):
            
        TracerProperties = [clsTracerProperties()]
        CalibrationFactor = self.inputs.CalibrationFactor # Will need to do something with this later
        
                
        if self.inputs.initialcond == True:
            Reach = pickle.load( open( self.inputs.priorReach, "rb" ) )
            for Node in Reach.Node:
                Node.DC.p = self.inputs.p
                Node.LayerL = 1.
                #Node.cbank = 0./ 24. / 60. / 60. / 365.25
            #Reach.Node = Reach.Node[0:30] # For dam removal to get rid of unstable backwater section              
                
        else:
            Reach = self.SetupReach(self.inputs.TransFunc, self.inputs.TrinityFit, CalibrationFactor) 
        
        TracerProperties[0] = self.SetupTracers(TracerProperties[0])
        counter = 0
        date = datetime.date(*self.inputs.StartDate)
        Tyear = 0.
        dt = self.inputs.dt[0] * 365.25*24*60*60 # seconds 631152 s = 0.02 yr or about 7 days
        dtcount = self.inputs.dtcount

        
        OutputSpecs1 = clsOutputSpecs(date)
        OutputSpecs2 = clsOutputSpecs(date)
        OutputSpecs3 = clsOutputSpecs(date)
        OutputSpecs4 = clsOutputSpecs(date)
        OutputSpecs5 = clsOutputSpecs(date)
        
        ValidateDates = map(lambda x: datetime.date(*x[:6]), self.inputs.ValidateDates)
                    
        # *********SETUP INITIAL CONDITIONS FOR PLOTTING AND FOR LATER USE ********
        # Get initial load conditions
        Reach.UpdateSlope()
        Reach.UpdateManningDepthAtAllFlowAndNodes()
        for Node in Reach.Node:
            Node.Load.UpdateSedimentLoadByDurationAndSize(Node.DC, \
                Node.ActiveLayer.GSD, 1000., 2.7, Node.Bc, \
                Node.FractionAlluvial, self.inputs.TransFunc, self.inputs.TrinityFit, CalibrationFactor) 

        InitialTotalBedloadFeed = Reach.Node[0].Load.QsavBedTotFeed
        InitialMigrationRate = Reach.Node[0].cbank
        # *************************************************************************
        
        MaxSteps = self.inputs.MaxSteps

        #  Katie modify:  Load factor count--if you provide integers, it will use
        #  the counter variable to determine when to switch sediment load.  If a
        #  tuple is provided, it will convert the tuple to a datetime object and use
        #  that.
        LoadFactor = self.inputs.LoadFactor
        LoadFactorCount = self.inputs.LoadFactorCount
        LoadType = ""
        if type(self.inputs.LoadFactorCount[0]) == int:
            LoadFactorCount = self.inputs.LoadFactorCount
            LoadType = 'counter'
        if type(self.inputs.LoadFactorCount[0]) == tuple:
            LoadFactorCount = map(lambda x: datetime.date(*x[:6]), self.inputs.LoadFactorCount)
            LoadType = 'date'
        
        #  Boundary change can be determined from either the counter or date, like
        #  the sediment load (see above)
        BoundaryChange = self.inputs.BoundaryChange
        BoundaryChangeCount = self.inputs.BoundaryChangeCount
        BoundaryType = ""
        if type(self.inputs.BoundaryChangeCount[0]) == int:
            BoundaryType = 'counter'
        if type(self.inputs.BoundaryChangeCount[0]) == tuple:
            BoundaryChangeCount = map(lambda x: datetime.date(*x[:6]), self.inputs.BoundaryChangeCount)
            BoundaryType = 'date'
        
        
        n = 0
        m = 0
        oldTyear = 0 # Katie add  
        Tmult = self.inputs.Tmult # Katie add
        subdaycount = 0 # Katie add
        
        NumberOfPrintouts = self.inputs.NumberOfPrintouts
        NumberOfAnimations = self.inputs.NumberOfAnimations
        Interval = int(MaxSteps / NumberOfPrintouts)
        AnimationInterval = int(MaxSteps / NumberOfAnimations)
        Printstep = 0
        VarPrintstep = 0
        #NextCount = 0
        NextCount = 0
        NextBoundary = 0
        NextAnimationCount = 0
        
        Reach.StepDownstream(0., self.inputs.AlphaBed, 0.2, 0.1, self.inputs.LayerL,\
            TracerProperties, self.inputs.TransFunc, self.inputs.TrinityFit,\
            CalibrationFactor, False, self.inputs.thresholdQ, self.inputs.W, self.inputs.ErodeT, self.inputs.Bp, self.inputs.vfunc) # Katie: may want to change other parameters to variables.  Width change variables are optional.
        
        #  Controlled variables for later
        ControlNode = deepcopy(Reach.Node[0]) # Katie add--for calculating feed from a rating curve (see update feed below)   
        zControlBoundary = deepcopy(Reach.Node[-1]) # Katie add--to set boundary condition with changing water depth        
        ReferenceAvkFeed = deepcopy(Reach.Node[0].Load.QsAvkFeed) # Katie add
        ReferenceQsjkFeed = deepcopy(Reach.Node[0].Load.Qsjkfeed)
                        
        #  Need to normalize  ReferenceQsjkFeed by duration--contains raw flux for duration curve--used in dam removal feed algorithm        
        NormalizedQsjkFeed = deepcopy(ReferenceQsjkFeed)        
        for k in range(ControlNode.NSizes + 1): # start at index k = 1 
            for j in range(ControlNode.DC.NFlows()):
                NormalizedQsjkFeed[j,k] = ReferenceQsjkFeed[j, k] * ControlNode.DC.p[j]        
                
        # Set up custom nodes here. 
        self.CustomNodes(Reach)
        
        if self.inputs.Hydrograph == True:
            self.set_up_hydrograph(Reach)
            Steps = len(self.inputs.Qlist)
        else:
            Steps = MaxSteps
            
        print 'Model setup complete!  Starting timesteps...'
        while counter < Steps:
            # ***************************OUTPUT RESULTS TO SPREADSHEET************* # Katie deleted a whole section here; see historic copies.
            if counter == NextCount: # Katie comment out
#            if counter % 1 == 0: # To print every timestep Katie add

                year = 0                
                if Tyear == 0:
                    year = 0.
                else:
                    year = Tyear-(dt/60./60./24./365.25)
                    
                for out in self.inputs.Outputvars:
                    self.Output('Out_'+ out, Reach, out, \
                        Printstep, year)
       
                self.OutputFlux(Reach, Printstep, year, 'QsOut') # Katie add
                self.OutputFlux(Reach, Printstep, year, 'QsIn') # Katie add

                with open(os.path.join(os.pardir, self.inputs.Outputfolder, 'NewControl'), 'w') as f:
                    f.write(str(counter))
                Printstep = Printstep + 1
                NextCount = NextCount + Interval
#                oldTyear = Tyear
            #print(counter)
#            if Tyear == 0 or Tyear > TyearOld + .01:
#                print Tyear
#                TyearOld = Tyear
                
            # Katie add--write output on certain dates in order to validate model.
            if date in ValidateDates and subdaycount == 1:
                for out in self.inputs.Validatevars:
                    self.Output('OutValidate_'+ out, Reach, out, \
                        VarPrintstep, float(date.year))
                VarPrintstep += 1
            
            if date.day == 1 and subdaycount == 1:
                print str(date)

            # **********************END OUTPUT*************************************
            
            # **********************CALCULATE HYDRAULICS***************    
            # For Hydrograph # Katie add
            if self.inputs.Hydrograph == True:                
                if subdaycount == Tmult or Tyear == 0:
                    subdaycount = 0
                    
                    #***********Output variables to object********************
                    OutputSpecs1.PopulateLists(Reach.Node[0])
                    OutputSpecs2.PopulateLists(Reach.Node[int(2*Reach.nnodes()/5.)])
                    OutputSpecs3.PopulateLists(Reach.Node[int(3*Reach.nnodes()/5.)])
                    OutputSpecs4.PopulateLists(Reach.Node[int(4*Reach.nnodes()/5.)])
                    OutputSpecs5.PopulateLists(Reach.Node[-1])
                    #*********************************************************
                    for Node in Reach.Node:
                        for j in range(Reach.NFlows):
                            Node.DC.Qw[j] = self.inputs.Qlist[counter]
                        Node.UpdateDepthAndDischargeAtAllFlows(self.inputs.vfunc)
                        for j in range(Reach.NFlows):
                            Node.DC.WSE[j] = Node.DC.Hc[j] + Node.etabav

                ####################################################################
#               # Find downstream boundary condition--will want to move this out of the hydrograph loop when simulating the removal with a flow duration curve.
                    Reach.Node[-1].DC.WSE = self.find_downstream_boundary(Reach, zControlBoundary, Reach.Node[-1].DC.Qw, setB = Reach.BoundaryHc)
                        
                    ####################################################################
            # **********************DO MASS CONSERVATION COMPUTATION***************    
                Reach.StepDownstream(dt, self.inputs.AlphaBed, 0.2, 0.1, self.inputs.LayerL, TracerProperties, self.inputs.TransFunc, self.inputs.TrinityFit, \
                    CalibrationFactor, self.inputs.WidthChange, self.inputs.thresholdQ, self.inputs.W, self.inputs.ErodeT, self.inputs.Bp, self.inputs.vfunc)                                
                Reach.UpdateSlope()                
                Reach.UpdateOutput(dt) # Katie add

            Tyear = Tyear + dt/ 365.25 / 24. / 60. /60.
            # *********************************************************************
    
            # ********************CHECK FOR BEDROCK *******************************
             #This section of code is experimental and is intended to crudely 
             #represent the exposure of bedrock, which would result in a partly 
             #alluvial bed and would halt incision where it is exposed
            for Node in Reach.Node:
                if Node.PartlyAlluvial == True: # Katie add PartlyAlluvial attribute to Node object
                    if Node.CumulativeBedChange < -0.2 and not Node.FixedElev: # Katie comment out
                        Node.FixedElev = True
                    if Node.FractionAlluvial > 1.:
                        Node.FixedElev = False
                        Node.ActiveLayer.Volume = Node.ActiveLayer.L * Node.Bc * \
                            Node.dxc
            # *********************************************************************

            
            # **********************ACTIVATE DAM REMOVAL--KATIE ADD***************
            # There are two processes--emptying of Aldwell and the consequent boundary
            # lowering, and the pulse of sediment from Lake Mills.            
            
            # Simulate the emptying of the Aldwell reservoir.  The water depth
            # begins at 30m and decreases incrementally over ~6 months (see Warrick
            # et al., 2015).  Then, when the depth gets down to 2m (arbitrarily
            # selected to be close enough to normal conditions but not too low to 
            # screw up the backwater calculation if there is a high flow), the 
            # boundary condition is set to depend on discharge once more.
            if date >= datetime.date (2011, 9, 17):
                Reach.SetBoundary = True # Should technically start as True in the inputs page
                timedelta = date - datetime.date(2011, 9, 17) # datetime.timedelta object with number of days elapsed since removal began
                Reach.BoundaryHc = 30.-.17*timedelta.days # 181 days between start of removal and base level point; 30 m/181 days ~ .17 m/day               
                if Reach.BoundaryHc < 2.:
                    Reach.SetBoundary = False
            
            # Katie added this section for the dam removal--when bedload reaches the 
            # former Glines site, the feed changes (below).  Based on the date
            # published in East (2015).
            if date == datetime.date(2012, 10, 14):
                self.inputs.FeedType = 'DamRemoval'
                DamTyearStart = Tyear
            # **********************END ACTIVATE DAM REMOVAL***************
            
            # *****************************UPDATE FEED*****************************

            # Update sediment load tracker
            LoadGauge = ""
            if LoadType == 'counter':
                LoadGauge = counter
            elif LoadType == 'date':
                LoadGauge = date

            if LoadGauge == LoadFactorCount[n]: # Katie change counter to LoadGauge--allows user to specify date in which load changes.
                if n < len(LoadFactor):
                    n += 1
                if n == len(LoadFactor) - 1:
                    LoadFactorCount = LoadFactorCount + ['filler'] # Katie add to keep script from crashing--LoadFactorCount's length is originally one smaller than LoadFactor's


            # Other time-sensitive parameters such as flow duration curve can also 
            # be updated here
            # Assume feed is 0% of capacity for all bed material
            
            # Standard method--designed for the flow duration curve.  Feed is based on a duration-averaged value                   
                for k in range(0, Reach.Node[0].NSizes + 1):   # Katie change range starting from 1 to 0...or else I think the mud feed isn't altered                
                    Reach.Node[0].Load.QsAvkFeed[k] = (ReferenceAvkFeed[k])*LoadFactor[n]
                    #Reach.Node[0].Load.QsAvkFeed[k] = (.002135)*LoadFactor[n]# Katie change for when starting from 0 feed
                    for j in range(Reach.Node[0].DC.NFlows()):
                        Reach.Node[0].Load.Qsjkfeed[j, k] = LoadFactor[n] * \
                            ReferenceQsjkFeed[j, k]
                

    ################################################################################            
#            #  Dam removal method # Katie add!
#            print ReferenceAvkFeed
            #print ReferenceQsjkFeed
#            ReferenceAvkFeed = np.array([1.09612004e-03,5.48060020e-04,2.21086265e-04,2.03938754e-04,\
#            1.74955003e-04,1.73541104e-04,1.85895866e-04,2.43681377e-04,\
#            9.80422282e-05,8.13663274e-06,1.92993808e-07])
#            ReferenceQsjkFeed = np.array([[0.00026878132620208845, 0.00013439066310104423, 3.542103133212312e-05, 2.6138185793052566e-05, 1.7252182357076505e-05, 1.2354204996358027e-05, 8.272881158626023e-06, 4.621639146231713e-06, 2.973589772790212e-07, 3.2870173822028424e-09, 4.015553583855542e-11], [0.003973628866615778, 0.001986814433307889, 0.0008136480629372178, 0.0007477947277530212, 0.0006346646894085722, 0.0006171266960434657, 0.0006344308621932922, 0.0007254161649506475, 0.00012710646459301026, 1.410878705740759e-06, 1.7235865784883655e-08], [0.007878032693055904, 0.003939016346527952, 0.0017287838750075929, 0.001640014355930321, 0.0014403216075062446, 0.0014559173550096093, 0.001579909180273495, 0.0020471300388793076, 0.000642102994600644, 1.0317941303134297e-05, 1.2604815052741097e-07], [0.012703189642204982, 0.006351594821102491, 0.0029050579748032546, 0.002807475787146348, 0.002514891741298256, 0.002599081512904223, 0.0029074715653942436, 0.004032993159145689, 0.001677583664708771, 4.645895882275157e-05, 5.675614604686764e-07], [0.018268011852770487, 0.009134005926385243, 0.004296146769715382, 0.004203772090815557, 0.0038155802410322136, 0.004001443701647746, 0.004566107150260334, 0.006617164906425084, 0.0032636429666528287, 0.0001555466895668803, 1.9002213682501524e-06], [0.024444964682390314, 0.012222482341195157, 0.005867838390929279, 0.005793807257107721, 0.005309127556185841, 0.005626700176727855, 0.006512634672374437, 0.009734529394649522, 0.0053895628921964926, 0.0004146328985735156, 5.211290909272014e-06], [0.031140502905703187, 0.015570251452851593, 0.007594296300765042, 0.007550702652552253, 0.006969658444054861, 0.007446076920671486, 0.008711872740841743, 0.013328018852767156, 0.008029526006354173, 0.0008560314928100303, 1.238282017909386e-05], [0.03645570471262003, 0.018227852356310015, 0.008977604718320993, 0.008964146901866659, 0.00831130984199701, 0.008923024001830338, 0.010508544184490841, 0.01630416711477033, 0.010327063087311837, 0.0013183508909896168, 2.203854256010635e-05], [0.04012922082062649, 0.020064610410313245, 0.009939059187843944, 0.009948974506804013, 0.00924853840162256, 0.009957706577458487, 0.011772023377995657, 0.01841431675622052, 0.012004466586323601, 0.0016929309070636235, 3.1469806412042644e-05], [0.04389755647030302, 0.02194877823515151, 0.010929314679590645, 0.010965092355889735, 0.010217331967237495, 0.011029407427435219, 0.013084263308669708, 0.020618679851158953, 0.013793231296464384, 0.002121857612427033, 4.4057235504175164e-05], [0.047755187358133916, 0.023877593679066958, 0.011946771550069743, 0.012010799339668754, 0.01121601318268437, 0.012136209038950697, 0.014442825367476144, 0.022912903702882313, 0.015689710987901735, 0.0026059145113202353, 6.059661839091822e-05], [0.05169713066782097, 0.025848565333910486, 0.01298997999278105, 0.01308455031586817, 0.01224305489340196, 0.013276360951279226, 0.015845475127585475, 0.025292965689841015, 0.01769039298850074, 0.003145587969773449, 8.202585106581885e-05], [0.05571887208262665, 0.027859436041313326, 0.014057620776729542, 0.014184936527172691, 0.013297061790844594, 0.014448260001940279, 0.01729015939152818, 0.027755141201337853, 0.01979190683822719, 0.0037411227419446155, 0.00010944091510404673], [0.0597827504641249, 0.02989137523206245, 0.015139543648977714, 0.015301432329862781, 0.014367890413110049, 0.015640556823819585, 0.018762776669578603, 0.030275042211610932, 0.021972791767421677, 0.0043870545729459715, 0.0001437982621462234], [0.0633013674123513, 0.03165068370617565, 0.01607861228248238, 0.016271545517659282, 0.015299361983639884, 0.016678953960170835, 0.020047393460968566, 0.0324808173633173, 0.02390451678192469, 0.004980742768205176, 0.00017977025131998516], [0.06652682056313333, 0.033263410281566666, 0.016941188936573454, 0.01716342266636962, 0.01615650000245919, 0.01763544308114694, 0.021232259179032525, 0.03452106469126449, 0.025708629711170317, 0.005552010226466863, 0.00021849567250445086], [0.06959937487978346, 0.03479968743989173, 0.017764336402927922, 0.018015186229864505, 0.016975743291249838, 0.01855044558847291, 0.022367051865494404, 0.0364799033126146, 0.027455354697411935, 0.00611947718495718, 0.0002610501237606703], [0.07507514546831993, 0.03753757273415997, 0.019234600986346957, 0.019538038457839194, 0.01844193219306445, 0.020189818320864115, 0.024403190314209233, 0.040005493669313996, 0.03063239610674278, 0.007184986176986215, 0.00034804580734510397]])

            if self.inputs.FeedType == 'DamRemoval': 
                
                C = 150
                tau = 1095.*24*60*60 # Convert from days to seconds
                t = (Tyear-DamTyearStart)*365.25*24*60*60 # Convert from years to seconds

                for k in range(0, Reach.Node[0].NSizes + 1):                   
                    Reach.Node[0].Load.QsAvkFeed[k] = (ReferenceAvkFeed[k])*(1 + C*math.exp(-t/tau))
                    for j in range(Reach.Node[0].DC.NFlows()):
                        Reach.Node[0].Load.Qsjkfeed[j, k] = (1 + C*math.exp(-t/tau)) * \
                        ReferenceQsjkFeed[j, k]
    
################################################################################     
            # Rating curve method # Katie add!
            if counter < Steps-1:                  
                if self.inputs.FeedType == 'RatingCurve':
                    
                    NormFlow = ControlNode.NormalChannelDepthAndDischarge(self.inputs.Qlist[counter], .01, self.inputs.vfunc)                
                    U = NormFlow[1]/(NormFlow[0]*ControlNode.Bc)
                    
                    FeedArray = ControlNode.Load.WilcockCroweLoadBySize(U, ControlNode.DC.Sf[0], ControlNode.ActiveLayer.GSD, 1000., 2.7, ControlNode.Bc, self.inputs.TrinityFit, self.inputs.CalibrationFactor)
                    
                    for k in range(1, Reach.Node[0].NSizes + 1):                   
                        Reach.Node[0].Load.QsAvkFeed[k] = FeedArray[k]*LoadFactor[n]
                        #Reach.Node[0].Load.QsAvkFeed[k] = (.002135)*LoadFactor[n]# Katie change for when starting from 0 feed
                        for j in range(Reach.Node[0].DC.NFlows()):
                            Reach.Node[0].Load.Qsjkfeed[j, k] = LoadFactor[n] * \
                                FeedArray[k] # Only works when there is one flow bin.
                                
                # Calculate mud load using Konrad-modified Curran relation for Elwha
                    Reach.Node[0].Load.QsAvkFeed[0] = Reach.Node[0].Load.CurranLoad(self.inputs.Qlist[counter])*LoadFactor[n]
                    Reach.Node[0].Load.Qsjkfeed[0, 0] = LoadFactor[n] * \
                                Reach.Node[0].Load.CurranLoad(self.inputs.Qlist[counter])*LoadFactor[n] # Only works when there is one flow bin.
                
            # *********************************************************************
            
            # *******************UPDATE BOUNDARY CONDITION (KATIE ADD)*************
            BoundaryGauge = ""
            if BoundaryType == 'counter':
                BoundaryGauge = counter
            elif BoundaryType == 'date':
                BoundaryGauge = date

            if BoundaryGauge == BoundaryChangeCount[NextBoundary]:

                # Hydrograph method                
                ############################################################
                Reach.SetBoundary = True
                Reach.BoundaryHc = BoundaryChange[NextBoundary]                

                if NextBoundary < len(BoundaryChangeCount) - 1: NextBoundary += 1
            
            # ************************ADJUST LATERAL EXCHANGE PARAMETERS***********
            # In this seciton, the bank migration rate, width, etc can be adjusted.
            # For now, bank migration rate is a function of load.  It would also be
            # easy to compute width as a function of surface size and bankfull 
            # capacity, but this would require a function for finding bankfull 
            # capacity--easy since can use steady uniform flow for channel.
            
            # Update the bank migration rate as a function of the bed material load
            #for Node in Reach.Node: # Katie comment out
            #    Node.cbank = InitialMigrationRate * Node.Load.QsavBedTotFeed / \
            #        InitialTotalBedloadFeed
            #    if Node.cbank > InitialMigrationRate:
            #        Node.cbank = InitialMigrationRate
            # *********************************************************************
            
            #  Update time trackers.  Counter only updates for a new day if hydrograph is on.
            if len(dtcount) > 0 and counter == dtcount[m]:
                dt = self.inputs.dt[m+1]* 365.25*24*60*60
                if m < len(dtcount) - 1:
                    m = m + 1
            
            subdaycount = subdaycount + 1 # Katie add
            if subdaycount == Tmult or self.inputs.Hydrograph == False:
                counter = counter + 1 # Katie add 
                date = date + datetime.timedelta(days=1)
        
        #  Save final Reach object
        pickle.dump(Reach, open(os.path.join(os.pardir, self.inputs.Outputfolder, "save.Reach"), "wb" )) 
        pickle.dump(self.inputs, open(os.path.join(os.pardir, self.inputs.Outputfolder, "inputparams.inputs"), "wb")) 
        
        datelist = map(lambda x: x.strftime('%Y,%m,%d'), OutputSpecs1.Date) 
        
        json.dump(datelist, open(os.path.join(os.pardir, self.inputs.Outputfolder, "save.DailyDate"), "wb" ))

        json.dump(OutputSpecs1.Q, open(os.path.join(os.pardir, self.inputs.Outputfolder, "save.DailyQ1"), "wb" )) 
        json.dump(OutputSpecs1.QsavBedTot, open(os.path.join(os.pardir, self.inputs.Outputfolder, "save.DailyQsTot1"), "wb" ))
        json.dump(OutputSpecs1.Qsk, open(os.path.join(os.pardir, self.inputs.Outputfolder, "save.DailyQsk1"), "wb" )) 
        json.dump(OutputSpecs1.F, open(os.path.join(os.pardir, self.inputs.Outputfolder, "save.DailyGSD1"), "wb" )) 
        json.dump(OutputSpecs1.SubF, open(os.path.join(os.pardir, self.inputs.Outputfolder, "save.DailySubGSD1"), "wb" ))  
        json.dump(OutputSpecs1.Bc, open(os.path.join(os.pardir, self.inputs.Outputfolder, "save.DailyBc1"), "wb" ))
        
        json.dump(OutputSpecs2.Q, open(os.path.join(os.pardir, self.inputs.Outputfolder, "save.DailyQ2"), "wb" )) 
        json.dump(OutputSpecs2.QsavBedTot, open(os.path.join(os.pardir, self.inputs.Outputfolder, "save.DailyQsTot2"), "wb" ))
        json.dump(OutputSpecs2.Qsk, open(os.path.join(os.pardir, self.inputs.Outputfolder, "save.DailyQsk2"), "wb" )) 
        json.dump(OutputSpecs2.F, open(os.path.join(os.pardir, self.inputs.Outputfolder, "save.DailyGSD2"), "wb" )) 
        json.dump(OutputSpecs2.SubF, open(os.path.join(os.pardir, self.inputs.Outputfolder, "save.DailySubGSD2"), "wb" ))
        json.dump(OutputSpecs2.Bc, open(os.path.join(os.pardir, self.inputs.Outputfolder, "save.DailyBc2"), "wb" ))

        json.dump(OutputSpecs3.Q, open(os.path.join(os.pardir, self.inputs.Outputfolder, "save.DailyQ3"), "wb" )) 
        json.dump(OutputSpecs3.QsavBedTot, open(os.path.join(os.pardir, self.inputs.Outputfolder, "save.DailyQsTot3"), "wb" ))
        json.dump(OutputSpecs3.Qsk, open(os.path.join(os.pardir, self.inputs.Outputfolder, "save.DailyQsk3"), "wb" )) 
        json.dump(OutputSpecs3.F, open(os.path.join(os.pardir, self.inputs.Outputfolder, "save.DailyGSD3"), "wb" )) 
        json.dump(OutputSpecs3.SubF, open(os.path.join(os.pardir, self.inputs.Outputfolder, "save.DailySubGSD3"), "wb" )) 
        json.dump(OutputSpecs3.Bc, open(os.path.join(os.pardir, self.inputs.Outputfolder, "save.DailyBc3"), "wb" )) 

        json.dump(OutputSpecs4.Q, open(os.path.join(os.pardir, self.inputs.Outputfolder, "save.DailyQ4"), "wb" )) 
        json.dump(OutputSpecs4.QsavBedTot, open(os.path.join(os.pardir, self.inputs.Outputfolder, "save.DailyQsTot4"), "wb" ))
        json.dump(OutputSpecs4.Qsk, open(os.path.join(os.pardir, self.inputs.Outputfolder, "save.DailyQsk4"), "wb" )) 
        json.dump(OutputSpecs4.F, open(os.path.join(os.pardir, self.inputs.Outputfolder, "save.DailyGSD4"), "wb" )) 
        json.dump(OutputSpecs4.SubF, open(os.path.join(os.pardir, self.inputs.Outputfolder, "save.DailySubGSD4"), "wb" )) 
        json.dump(OutputSpecs4.Bc, open(os.path.join(os.pardir, self.inputs.Outputfolder, "save.DailyBc4"), "wb" )) 

        json.dump(OutputSpecs5.Q, open(os.path.join(os.pardir, self.inputs.Outputfolder, "save.DailyQ5"), "wb" )) 
        json.dump(OutputSpecs5.QsavBedTot, open(os.path.join(os.pardir, self.inputs.Outputfolder, "save.DailyQsTot5"), "wb" ))
        json.dump(OutputSpecs5.Qsk, open(os.path.join(os.pardir, self.inputs.Outputfolder, "save.DailyQsk5"), "wb" )) 
        json.dump(OutputSpecs5.F, open(os.path.join(os.pardir, self.inputs.Outputfolder, "save.DailyGSD5"), "wb" )) 
        json.dump(OutputSpecs5.SubF, open(os.path.join(os.pardir, self.inputs.Outputfolder, "save.DailySubGSD5"), "wb" )) 
        json.dump(OutputSpecs5.Bc, open(os.path.join(os.pardir, self.inputs.Outputfolder, "save.DailyBc5"), "wb" ))                 		

        # Produce done file
        donefile = os.path.join(os.pardir, self.inputs.Outputfolder, 'Model_finished')
        with open (donefile, 'w') as f:
            f.write('Model completed')
            f.close()		
           
           # ****************MAKE CHARTS WHEN ALL COMPUTATIONS ARE DONE***********        
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
     
        outputpath = os.path.join(os.pardir, self.inputs.Outputfolder, Filename)
        # VariableName is a string that describes the reach property to be printed.
        # For example, if the property reach.node(i).etabav is to be printed, 
        # VariableName should be "etabav."
        if Nprint == 0:
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
        else:
            name = 'QsIn'
            outlist = Reach.CumulativeFeed
            
        outputpath = os.path.join(os.pardir, self.inputs.Outputfolder, name)
        if Nprint == 0:
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

#if __name__ == '__main__':
#    main()
