from copy import deepcopy
import os
from clsReach import clsReach
from clsTracerProperties import clsTracerProperties
from clsOutputWriter import clsOutputWriter
from clsReservoir import clsReservoir
from clsSubstratePairClass import clsSubstratePairClass
import datetime
import pickle as pickle
import numpy as np
from math import exp
import math

class clsModel(object):
    
    def __init__(self, inputs):
        
        self.inputs = inputs
        
    def SetInitialConditions(self, inputs):
        """
        Section I:  Load the Reach object from specifications in 'InputsMainpage...'
        and perform initial hydraulic and sediment transport computations to 
        fully define the initial condition.
        
        Parameters
        ----------
        inputs : :obj:`MAST_1D.clsInputs`
            The class defining all in put parameters.
        Returns
        -------
        Reach : :obj:`MAST_1D.clsReach`
            The initialized reach.
        ControlNode : :obj:`MAST_1D.clsNode`
            A node with basic default input parameters.
        ControlGSD : :obj:`MAST_1D.clsGSD`
            A basic grain size distribution.
        zControlBoundary : :obj:`MAST_1D.clsNode`
            A node representing the downstreawm boundary
        TracerProperties : :obj:`MAST_1D.clsTracerProperties`
        CalibrationFactor : float
            Sediment transport calibration factor. Used to adjust
            reference Shields stress in Wilcock Crowe type calculation.
            (UNSURE IF IT IS BEING SET HERE).
        
        """
            
        #  Sets reach parameters by either loading existing conditions from a 
        #  prior run or creating a reach object using the user-defined inputs.
        Reach = clsReach(inputs)
        Reach.SetupReach(inputs)
        
        #  Set the calibration factor for the sediment transport equation
        CalibrationFactor = self.inputs.CalibrationFactor # Sets sediment transport calibration

        #  Sets initial tracer conditions.  These functions are untested in 
        #  this version of MAST-1D.
        TracerProperties = [clsTracerProperties()]
        TracerProperties[0] = Reach.SetupTracers(self.inputs, TracerProperties[0])

        #  Set the initial hydraulic and sediment transport conditions
        Reach.UpdateSlope()
        Reach.UpdateManningDepthAtAllFlowAndNodes()
        for Node in Reach.Node:
            Node.Load.UpdateSedimentLoadByDurationAndSize(Node.DC, \
                Node.ActiveLayer.GSD, 1000., 2.7, Node.Bc, \
                Node.FractionAlluvial, self.inputs.TransFunc,\
                self.inputs.TrinityFit, CalibrationFactor) 
               
        Reach.StepDownstream(0., self.inputs.AlphaBed, 0.2, 0.1, self.inputs.LayerL,\
            TracerProperties, self.inputs.TransFunc, self.inputs.TrinityFit,\
            CalibrationFactor, False, self.inputs.W,\
            self.inputs.ErodeT, self.inputs.alphatau, self.inputs.BcMin, self.inputs.vfunc, self.inputs.AvulsionThreshold, Reach.Node[0].ActiveLayer.GSD, self.inputs.AvulsionExchange, self.inputs.MobilityThreshold) # Katie: may want to change other parameters to variables.  Width change variables are optional.
        
        #  Preserve control nodes to be used to set boundary conditions later
        ControlNode = deepcopy(Reach.Node[0]) # For calculating feed from a rating curve (see update feed below)
        ControlGSD = deepcopy(Reach.Node[0].ActiveLayer.GSD)
        zControlBoundary = deepcopy(Reach.Node[-1]) # To set boundary condition with changing water depth                       
        #  Set reach to existing conditions if desired after control nodes have been set
        if inputs.initialcond == True:
            Reach = pickle.load(open(inputs.priorReach, "rb" ))
            for Node in Reach.Node: # Resets widening and narrowing in case width change is turned off.
                Node.NarrowRate = 0.
                Node.WidenRate = 0.
            for Node in Reach.Node:
                setattr(Node, 'TotalSubstrateVolume', 0.)
            zControlBoundary.etabav = Reach.Node[-1].etabav
            for j in range(Reach.NFlows):
                Reach.Node[-1].DC.WSE[j] = Reach.Node[-1].etabav + inputs.BoundaryFactor[0]
            #Temp
            #Reach.SetBoundary = False
        
        #  Sets up node-specific values for select variables, if provided by the
        #  user in the CustomNodes function in Section 2.
        if inputs.initialcond == False:
            self.CustomNodes(Reach)
            for j in range(Reach.NFlows):
                zControlBoundary.DC.WSE[j] = Reach.Node[-1].DC.WSE[j]
            # Reset controlNode bed elevation
            zControlBoundary.etabav = Reach.Node[-1].etabav
            #print Reach.Node[-1].etabav
                        
        return Reach, ControlNode, ControlGSD, zControlBoundary, TracerProperties, CalibrationFactor
            
    def CustomNodes(self, Reach):
        """
        Section II:  Control the behavior of individual nodes if desired.
        Specify which input variables have node-specific initial or
        boundary conditions.  The user can add/remove custom variables.
        
        Parameters
        ----------
        Reach : :obj:`MAST_1D.clsReach`
            The initialized reach.
        """
        
        #  This function is run after initial floodplain conditions and Control 
        #  Nodes are set.
        
        for i in range(Reach.nnodes()):

            Node = Reach.Node[i]
            
            #  Clear some attributes from original node--not sure if this is 
            #  necessary--may be relic from previous version.
            Node.Bcrate = 0.
            Node.NarrowRate = 0.
            Node.WidenRate = 0.
            Node.Floodplain.Volume = 0.
            
            #  Custom node variables are set here.
            
            #  Checks to see if the node-specific variable exists in the provided inputs
            #  Here, BcNodes is the channel width
            if hasattr(self.inputs, 'BcNodes'):
            
            #  If the variable exists, it is changed to the node specific value.
            #  (Note:  if it doesn't exist, the value will keep the non-node-specific
            #  value also provided in the inputs.)
                Node.Bc = self.inputs.BcNodes[i]
                
            #  Repeat for all 'custom' variables
            #  Here, BfNodes is the valley width and ChSinNodes is sinuosity.
            if hasattr(self.inputs, 'BfNodes'):
                Node.Bf = self.inputs.BfNodes[i]-Node.Bc  
                Node.ValleyMargin = self.inputs.BfNodes[i]
            if hasattr(self.inputs, 'ChSinNodes'):
                Node.ChSin = self.inputs.ChSinNodes[i]
            if hasattr(self.inputs, 'SlopeNodes'):
                Node.Slope = self.inputs.SlopeNodes[i]
                
                # Redoes bed elevation to incorporate new slope.
                Node.etabav = 0.
                if i == 0:
                    Node.etabav = 127.
                else:
                    Node.etabav = Reach.Node[i-1].etabav-Reach.Node[i-1].Slope*Reach.Node[i-1].dxc
                #Node.etabav = 33.9 - Node.Slope * Node.dxc * i
                Node.InitialBedElev = Node.etabav

            #  Redoes channel distance/elevation for custom sinuosity
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

                
        #  *****************Example: ELWHA PARTLY ALLUVIAL REACHES*************
        #  This section is an example of how the CustomNodes function can be used
        #  to customize the behavior of particular reaches.  Here, particular
        #  nodes are set as bedrock reaches.  They are not allowed to widen or 
        #  narrow, and the maximum degradation is set using the partially-alluvial
        #  functions in MAST-1D.
                
        #  Partly-alluvial reaches are set after the floodplain number calculation
        #  is done so that the migration rate can be set to zero.  In theory, a 
        #  partly-alluvial node with a migrating floodplain can exist.  Here, however,
        #  for the purpose of the Elwha River, I defined 'Canyon' nodes in the input
        #  class, and these nodes will be set to have migration rates of zero, no
        #  width change, and very high floodplain heights to simulate a canyon.
               
            if hasattr(self.inputs, 'Canyon'):            
                if self.inputs.Canyon[i] == True:
                    Node.Canyon = True
                    Node.cbank = 0. # Set migration rate to zero.
                    Node.PartlyAlluvial = True
                    Node.Floodplain.L = 15. # Arbitrarily set to be very high.
                    #Node.ActiveLayer.GSD.F = np.zeros(Node.NSizes + 1)
                    #Node.ActiveLayer.GSD.F[-1] = 1.
            #setattr(Node, 'ValleyMargin', 500.)        
            Node.Floodplain.Volume = Node.Floodplain.L * Node.Bf * Node.dxc / \
                Node.ChSin
            Node.ActiveLayer.Volume = Node.ActiveLayer.L*Node.Bc*Node.dxc
            TotalVolume = 0
            for m in range(Node.NLayers()):
                Node.Substrate[m].F.Volume = Node.Substrate[m].F.L*Node.Bf*Node.dxc/Node.ChSin
                Node.Substrate[m].C.Volume = Node.Substrate[m].C.L*Node.Bc*Node.dxc
                TotalVolume = TotalVolume + Node.Substrate[m].C.Volume + Node.Substrate[m].F.Volume
            Node.TotalSubstrateVolume = TotalVolume
        # *****************END ELWHA PARTLY ALLUVIAL REACHES*******************
        if self.inputs.SetBoundary==True:
            for j in range(Reach.NFlows):
                Reach.Node[-1].DC.WSE[j] = Reach.Node[-1].etabav + self.inputs.BoundaryFactor[0]
        
    def RunModel(self):
        """
        Perform the model run.
        """        
        #**********************************************************************
        #Section III:  Set initial conditions.
                
        #  Loads reach inputs and sets initial bed and floodplain conditions.
        Reach, ControlNode, ControlGSD, zControlBoundary, TracerProperties, CalibrationFactor = \
            self.SetInitialConditions(self.inputs)
        
        ## temp workaround for old reach object
        #setattr(Reach, 'CumulativeBankSupply', np.zeros(Reach.NBedSizes + 1))
        #setattr(Reach, 'CumulativeBankSink', np.zeros(Reach.NBedSizes + 1))
        #for Node in Reach.Node:
            #setattr(Node, 'ValleyMargin', 500.)
            #setattr(Node, 'TotalSubstrateVolume', 0.)
        
        ##print self.CumulativeBankSupply
        ##print self.CumulativeBankSupply
        
        #  Reset hydraulics to accommodate hydrograph if necessary
        if self.inputs.Hydrograph == True:
            Reach.set_up_hydrograph()
            
        # Hard-code the floodplain number
        #for Node in Reach.Node:
         #   Node.Flmud = .05
            
        #**********************************************************************
        #Section IV:  Set boundary conditions.  Only simple changes to the feed rate
        #and downstream water surface elevation are set here.  More complicated 
        #algorithms can be hard-coded in Section VI, although some may require
        #global variables to be defined here.
         
        #  Load factor sets the timing of changes in upstream sediment feed.
        #  If you provide integers, it will use the timestep number to determine
        #  when to switch sediment load.  If a tuple is provided, it will convert 
        #  the tuple to a datetime object and determine the feed change using the
        #  date.
        LoadFactor = self.inputs.LoadFactor
        LoadFactorCount = self.inputs.LoadFactorCount
        LoadType = ""
        if type(self.inputs.LoadFactorCount[0]) == int:
            LoadFactorCount = self.inputs.LoadFactorCount
            LoadType = 'counter'
        if type(self.inputs.LoadFactorCount[0]) == tuple:
            LoadFactorCount = map(lambda x: datetime.date(*x[:6]),\
                self.inputs.LoadFactorCount)
            LoadType = 'date'
        n = 0 # Keeps track of which item in LoadFactorCount is the current trigger
        NewAvkFeed = Reach.Node[0].Load.QsAvkFeed #Katie add: Updated feed goes here, and feed in Node objects are filled after results are written to file
        NewJKFeed = Reach.Node[0].Load.Qsjkfeed
        
        #  Boundary change can be determined from either the counter or date, like
        #  the sediment load (see above)        
        BoundaryFactor = self.inputs.BoundaryFactor
        BoundaryFactorCount = self.inputs.BoundaryFactorCount
        BoundaryType = ""
        if self.inputs.Hydrograph == False:
            BoundaryType = 'counter'
        else:
            BoundaryFactorCount = list(map(lambda x: datetime.date(*x[:6]),\
                self.inputs.BoundaryFactorCount))
            BoundaryType = 'date'
        NextBoundary = 0 #  Keeps track of which item in BoundaryChangeCount is the current trigger      
        
        #  This customized global variable is needed to trigger a more sophisticated 
        #  boundary condition change in Section VI. 
        BoundaryTrigger = False # Trigger that determines whether nodes are added in a Dam Removal situation      
        NodeTrigger = False        
        
        #**********************************************************************
        #Section V:  Define time-keeping counters and set up output
        
        
        #  Timesteps elapsed for number flow duration curves or the number of 
        #  days elapsed for hydrographs
        counter = 0 
        
        #  This is a custom counter for monitoring dam removal feed
        startcounter = 0
        
        #  The date.
        date = datetime.date(*self.inputs.StartDate)
        
        #  The time elapsed in years.
        Tyear = 0
        
        #  Set the length of the timestep and determine when it changes
        #  when a flow duration curve is used
        dt = self.inputs.dt[0] * 365.25*24*60*60 # initial timestep in seconds
        dtlist = self.inputs.dt # list of timesteps
        dtcount = self.inputs.dtcount # list of when timestep changes occur
        m = 0 # keeps track of which integer in dtcount is the current trigger
        MaxSteps = self.inputs.MaxSteps # Time counter--Total number of timesteps
        
        #  Set the length of the timestep and determine when it changes
        #  when a hydrograph is used
        Tmult = self.inputs.Tmult # Number of timesteps in a day
        subdaycount = 0 # Timestep counter within a day
        if self.inputs.Hydrograph == True:
            dt = 1./365.25/Tmult*365.25*24*60*60 # In seconds
            MaxSteps = len(self.inputs.Qlist) # Time counter--Total number of days
        LowFlow = False # Sets low flow tag
        
        try:
            self.inputs.CyclingHydrograph
        except:
            self.inputs.CyclingHydrograph = False
        if self.inputs.CyclingHydrograph == True:
            dt = 1./365.25/27*365.25*24*60*60
            Reach.Node[-1].Canyon = True

        #  Set up output
        OutputObj = clsOutputWriter(self.inputs.Outputfolder, self.inputs.DailyNodes, date)
        
        #  'Standard' output files that are printed at consistent time intervals
        NumberOfPrintouts = self.inputs.NumberOfPrintouts # Number of times during run to record output        
        Interval = MaxSteps / NumberOfPrintouts # Interval in time counters between output recordings
        Printstep = 0 # Locator for standard output files
        NextCount = 0 # Keeps track of what value of counter will trigger the next standard print interval
        
        #  Output files that are written at user-defined dates for model performance testing      
        ValidateDates = list(map(lambda x: datetime.date(*x[:6]), self.inputs.ValidateDates)) # List of dates to output data
        VarPrintstep = 0 # Locator for performance-testing output files  

        #**********************************************************************
        #Section VI:  Iterate through timesteps.        
        
        if self.inputs.Hydrograph == True:
            Reach.Node[0].Load.QsAvkFeed, Reach.Node[0].Load.Qsjkfeed = Reach.Node[0].Load.UpdateElwhaFeedRatingCurve(ControlNode, LoadFactor[n],\
                self.inputs.MudFraction, self.inputs.Qlist[counter],\
                self.inputs.vfunc, self.inputs.TrinityFit, CalibrationFactor, ControlGSD, Section = 'Middle')       
        
        print('Model setup complete!  Starting timesteps...')

        while counter < MaxSteps and Reach.Node[3].ActiveLayer.GSD.D50 == Reach.Node[3].ActiveLayer.GSD.D50:
            #************************************************
            #VI.A:  Record output if at the proper timestep

            #  If time tracker is at the defined interval, write 'standard' output
            #  to file and save a copy of Reach in case of model crashing.
            #if Tyear == 0 or date.month == 10 and date.day == 1 and subdaycount == 1:
            #if counter % 1 == 0: # To print every timestep for debugging purposes
            #if subdaycount == Tmult: # To output every day for debugging purposes
            if counter == int(NextCount):
                NextCount = NextCount + Interval
                year = 0                
                if Tyear == 0:
                    year = 0.
                else:
                    year = Tyear-(dt/60./60./24./365.25)
                
                print('Saving regular output at count = %s' % counter) 
                for out in self.inputs.Outputvars:
                    OutputObj.Output('Out_'+ out, Reach, out, \
                        Printstep, year)

                OutputObj.OutputFlux(Reach, Printstep, year, 'QsOut') # Katie add
                OutputObj.OutputFlux(Reach, Printstep, year, 'QsIn') # Katie add
                OutputObj.OutputFlux(Reach, Printstep, year, 'BankIn') # Katie add
                OutputObj.OutputFlux(Reach, Printstep, year, 'BankOut') # Katie add
                
                Printstep += 1
                #print 'nc'
                #print Reach.Node[0].nc()
                #print Reach.Node[-3].nc()
                #pickle.dump(Reach, open(os.path.join(os.pardir, self.inputs.Outputfolder, "save.Reach"), "wb" )) 
                #pickle.dump(self.inputs, open(os.path.join(os.pardir, self.inputs.Outputfolder, "inputparams.inputs"), "wb"))
                
            #  If date is one specified for performance testing, write output.
            if date in ValidateDates and subdaycount == 1:
                print('Saving validation output at specified date of: %s ' % date)
                for out in self.inputs.Validatevars:
                    OutputObj.Output('OutValidate_'+ out, Reach, out, \
                        VarPrintstep, float(date.year))
                VarPrintstep += 1
            
            #************************************************
            #VI.B:  Perform mass conservation (for more details, see the 'Step Downstream'
            #function in the Reach class)
            
            #  If the model is running a hydrograph, the discharge is updated
            #  here at the start of each day and daily output variables are 
            #  recorded in a list.
            if self.inputs.Hydrograph == True:                
                if subdaycount == Tmult or Tyear == 0: #  Tmult denotes the number of timesteps in each day.
                    subdaycount = 0
                    FeedTrigger = False
                    
                    #  Record output variables
                    OutputObj.PopulateDailyLists(Reach)
                    LowFlow = False
                    #  Update the discharge
                    #print self.inputs.Qlist[counter]
                    for Node in Reach.Node:
                        for j in range(Reach.NFlows):
                            Node.DC.Qw[j] = self.inputs.Qlist[counter]
                        Node.UpdateDepthAndDischargeAtAllFlows(self.inputs.vfunc)
                        for j in range(Reach.NFlows):
                            Node.DC.WSE[j] = Node.DC.Hc[j] + Node.etabav

                    # Find downstream boundary condition for the new discharge.                   
                    Reach.find_downstream_boundary(zControlBoundary)   
                    
                    # Make it so that sediment transport only occurs if flow is above
                    # threshold.  Increase timestep in that case. Only works with hydrograph
                    if Reach.Node[0].DC.Qw[0] < 40.:
                        LowFlow = True
                        dt = dt*Tmult
                        Tmult = 1
                    elif Reach.Node[0].DC.Qw[0] > 300. and self.inputs.CyclingHydrograph == False:
                        Tmult = 150
                        dt = 1./365.25/Tmult*365.25*24*60*60
                    else:
                        Tmult = self.inputs.Tmult
                        if self.inputs.CyclingHydrograph == True:
                            dt = 1./365.25/27*365.25*24*60*60
                        else:
                            dt = 1./365.25/Tmult*365.25*24*60*60

            else:
                subdaycount = 1
                

            # Apply feed update after output has been written to file
            Reach.Node[0].Load.QsAvkFeed, Reach.Node[0].Load.Qsjkfeed = NewAvkFeed, NewJKFeed            
            
            #  Run through the hydraulics, sediment transport, and mass conservation
            #  computations. See 'Step Downstream' in the Reach class for more
            #  details.
            if LowFlow == True:
                Reach.AbridgedStepDownstream(self.inputs.BcMin, self.inputs.W, dt, TracerProperties, self.inputs.WidthChange, self.inputs.TransFunc, self.inputs.TrinityFit, self.inputs.CalibrationFactor, self.inputs.alphatau, self.inputs.AlphaBed)
                
            else:
                Reach.StepDownstream(dt, self.inputs.AlphaBed, 0.2, 0.1, self.inputs.LayerL,\
                    TracerProperties, self.inputs.TransFunc, self.inputs.TrinityFit, \
                    CalibrationFactor, self.inputs.WidthChange,\
                    self.inputs.W, self.inputs.ErodeT, self.inputs.alphatau, self.inputs.BcMin, self.inputs.vfunc,self.inputs.AvulsionThreshold, ControlGSD, self.inputs.AvulsionExchange, self.inputs.MobilityThreshold)                                
                Reach.UpdateSlope()                
            
            #  Update some cumulative variables and the time elapsed            
            Reach.UpdateOutput(dt)
            Tyear = Tyear + dt/ 365.25 / 24. / 60. /60.

            #************************************************
            #VI.C:  Check for bedrock and perform partly-alluvial computations if
            #activated
    
            #This section of code is experimental and is intended to crudely 
            #represent the exposure of bedrock, which would result in a partly 
            #alluvial bed and would halt incision where it is exposed            
            for Node in Reach.Node:
                if Node.PartlyAlluvial == True:
                    if Node.CumulativeBedChange < self.inputs.PartlyAlluvialMin and not Node.FixedElev:
                        Node.FixedElev = True
                    if Node.FractionAlluvial > 1.:
                        Node.FixedElev = False
                        Node.ActiveLayer.Volume = Node.ActiveLayer.L * Node.Bc * \
                            Node.dxc
            #if Reach.Node[-1].CumulativeBedChange > 1.: # Katie add
            #    Reach.Node[-1].etabav = Reach.Node[-1].InitialBedElev + 1.
            #if Reach.Node[-1].CumulativeBedChange < .5: # Katie add
            #    Reach.Node[-1].etabav = Reach.Node[-1].InitialBedElev -.5

            #************************************************
            #VI.D:  EXAMPLE CUSTOMIZATION:  Activate dam removal.  This piece of
            #the model was added to trigger the sediment feed and downstream boundary
            #changes associated with dam removal in the Middle Elwha River. A boolean
            #variable called 'Removal' was added to the inputs class to trigger 
            #this section of code if set to True.  If this section is unwanted,
            #the user can either set inputs.Removal to False or comment out this
            #section.
              
            # There are two processes--emptying of Lake Aldwell and the consequent 
            # boundary lowering, and the pulse of sediment from Lake Mills.            
            
            if self.inputs.Removal == True:
                #  Prevent the alluvial node between the canyons from widening
#                Reach.Node[1].Canyon = True
#                Reach.Node[1].NarrowRate = 0.
#                Reach.Node[1].WidenRate = 0.

                #  Trigger removal at the proper time and set up time variables
                if date >= datetime.date (2011, 9, 17) or self.inputs.Hydrograph == False and counter >= 1860:
                    deltat = 0  # Time elaspsed since removal                  
                    if self.inputs.Hydrograph == True:                    
                        timedelta = date - datetime.date(2011, 9, 17) # datetime.timedelta object with number of days elapsed since removal began
                        deltat = timedelta.days
                    else:
                        timeseconds = (counter-1860)*dt#9300
                        deltat = timeseconds/(60.*60.*24.)
                    Reach.BoundaryHc = 30.-.17*deltat # 181 days between start of removal and base level point; 30 m/181 days ~ .17 m/day                            
                   # Simulate the emptying of the Aldwell reservoir.  The water depth
                    # begins at 30m and decreases incrementally over ~6 months (see Warrick
                    # et al., 2015).  Then, when the depth gets down to 2m (arbitrarily
                    # selected to be close enough to normal conditions but not too low to 
                    # screw up the backwater calculation if there is a high flow), the 
                    # boundary condition is set to depend on discharge once more.
                    # In runs with a hydrograph, nodes are added to the bottom of
                    # the reach for stability.     
    
                    if Reach.BoundaryHc < 3. and BoundaryTrigger == False:
                        if self.inputs.Hydrograph == True:                        
                            BoundaryTrigger = True                        
                            Reach.SetBoundary = False
#                            NewNodes = clsReach(self.inputs)
#                            NewNodes.SetupReach(self.inputs)
#                            Reach.Add_Nodes(NewNodes.Node[0], 5, zControlBoundary.etabav, zControlBoundary.xc, self.inputs.Outputvars,self.inputs.Validatevars,self.inputs.ValidateDvars, self.inputs.Outputfolder)
#                            zControlBoundary = deepcopy(Reach.Node[-1]) # Set downstream boundary to new downstreammost node                       
#                            Reach.Node[-1].PartlyAlluvial = True
#                            Reach.Node[-1].Canyon = True
#                            Reach.Node[-1].Bc = 40.
#                            Reach.Node[-2].PartlyAlluvial = True
#                            Reach.Node[-2].Canyon = True
#                            Reach.Node[-2].Bc = 40.
#                            Reach.Node[-3].PartlyAlluvial = True
#                            Reach.Node[-3].Canyon = True
#                            Reach.Node[-3].Bc = 40.
#                            Reach.Node[-4].PartlyAlluvial = True
#                            Reach.Node[-4].Canyon = True
#                            Reach.Node[-4].Bc = 40.
#                            Reach.Node[-5].PartlyAlluvial = True
#                            Reach.Node[-5].Canyon = True
#                            Reach.Node[-5].Bc = 40.
                            
                            # Reset hydraulics
                            Reach.set_up_hydrograph()
                            for Node in Reach.Node:
                                for j in range(Reach.NFlows):
                                    Node.DC.Qw[j] = self.inputs.Qlist[counter]
                                Node.UpdateDepthAndDischargeAtAllFlows(self.inputs.vfunc)
                                for j in range(Reach.NFlows):
                                    Node.DC.WSE[j] = Node.DC.Hc[j] + Node.etabav
                            Reach.find_downstream_boundary(zControlBoundary)
                        
                        else:
                            BoundaryTrigger = True
                            Reach.SetBoundary = False
                            #Reach.Add_Nodes(zControlBoundary, 5, zControlBoundary.etabav, zControlBoundary.xc, self.inputs.Outputvars, self.inputs.Outputfolder)
                            #zControlBoundary = deepcopy(Reach.Node[-1]) # Set downstream boundary to new downstreammost node                       
                            for Node in Reach.Node:
                                Node.UpdateDepthAndDischargeAtAllFlows(self.inputs.vfunc)
                                for j in range(Reach.NFlows):
                                    Node.DC.WSE[j] = Node.DC.Hc[j] + Node.etabav
                            Reach.find_downstream_boundary(zControlBoundary)

                    # Find downstream boundary condition--will want to move this out of the hydrograph loop when simulating the removal with a flow duration curve.
                    Reach.find_downstream_boundary(zControlBoundary)

                # This section simulates the upstream sediment feed change
                # for the dam removal--when bedload reaches the 
			 # former Glines site, the feed changes (below).  Based on the date
			 # published in East (2015).                        
                if date == datetime.date(2012, 10, 14) or self.inputs.Hydrograph == False and counter == 1880:                   
                    
                    if self.inputs.Hydrograph == True:
                        if NodeTrigger == False:
                            NodeTrigger = True
                            self.inputs.FeedType = 'DamRemovalFirstPulseRatingCurve'
                            # Reduce the thickness of the active layer--only for 
                            # hydrograph and add remaining active layer as new 
                            # substrate layer
#                            for Node in Reach.Node:
#                                NewALL = .1
#                                Node.Substrate.append(clsSubstratePairClass())
#                                Node.Substrate[-1].C = clsReservoir(range(Node.Substrate[-2].C.NSizes + 2), Node.NTracers)
#                                Node.Substrate[-1].F = clsReservoir(range(Node.Substrate[-2].F.NSizes + 2), Node.NTracers)
#                                Node.Substrate[-1].C.GSD = deepcopy(Node.ActiveLayer.GSD)
#                                Node.Substrate[-1].C.L = Node.ActiveLayer.L-NewALL
#                                Node.Substrate[-1].C.Volume = Node.Substrate[-1].C.L*Node.Bc*Node.dxc
#                                Node.Substrate[-1].F.GSD = deepcopy(Node.Floodplain.GSD)
#                                Node.Substrate[-1].F.L = Node.ActiveLayer.L-NewALL
#                                Node.Substrate[-1].F.Volume = Node.Substrate[-1].F.L*Node.Bf*Node.dxc/Node.ChSin
#                                Node.ActiveLayer.L = NewALL                            
#                                Node.ActiveLayer.Volume = Node.ActiveLayer.L*Node.Bc*Node.dxc
                                
                    else:
                        self.inputs.FeedType = 'DamRemovalDuration'
                    DamTyearStart = Tyear
                    
                    startcounter = counter
                    
                if date == datetime.date(2013, 10, 1):
                    if self.inputs.Hydrograph == True:                    
                        self.inputs.FeedType = 'DamRemovalFirstPulseRatingCurve'
                    #DamTyearStart = Tyear
                    # Counter at feed change for feed code below
                    #startcounter = counter
                    
            #************************************************
            #VI.E:  Update sediment feed.  Custom sediment feed regimes can be called
            #here, but they should be written in the Load class.
            
            # Trigger a change in feed at the proper time--used for the 'default'
            # feed function.
            LoadGauge = ""
            if LoadType == 'counter':
                LoadGauge = counter
            elif LoadType == 'date':
                LoadGauge = date
            if n < len(LoadFactorCount):
                if self.inputs.Hydrograph == False and str(LoadGauge) == str(LoadFactorCount[n]):
                    n += 1
                elif self.inputs.Hydrograph == True and str(LoadGauge) == str(LoadFactorCount[n]) and FeedTrigger == False:#subdaycount == 1: # Katie change counter to LoadGauge--allows user to specify date in which load changes.
                    FeedTrigger = True                    
                    n += 1
            
            # A method for calculating feed is called here based on the FeedType
            # variable specified by the user on the inputs page.            
            
            # Standard method--designed for the flow duration curve.  
            # Feed is a user-defined fraction of sediment transport capacity                   
            if self.inputs.FeedType == 'DurationCurve':            
                #Node =  Reach.Node[0]   
                NewAvkFeed, NewJKFeed = Node.Load.UpdateFeedDurationCurve(\
                    ControlNode, LoadFactor[n])  
      
            # Rating curve method--the standard for hydrograph runs.
            # Feed is a user-defined fraction of sediment transport capacity
            # for a given flow.  Mud feed is a set proportion of feed for the
            # next finest size class.
            if counter < MaxSteps-1:                                  
                if self.inputs.FeedType == 'RatingCurve' and subdaycount == Tmult - 1 and self.inputs.Qlist[counter + 1]>40.:                
                    Node = Reach.Node[0]
                    NewAvkFeed, NewJKFeed = \
                        Node.Load.UpdateFeedRatingCurve(ControlNode, LoadFactor[n],\
                        self.inputs.MudFraction, self.inputs.Qlist[counter + 1],\
                        self.inputs.vfunc, self.inputs.TrinityFit, CalibrationFactor, ControlGSD) 
            
            # ******************Example: ELWHA RATING CURVE *******************
            # This is an example of a customized method of updating the sediment
            # feed.            
            
            # It is the same as 'RatingCurve' method above when calculating feed 
            # for the bed material load, but uses an emperical sediment rating 
            # curve provided by Curran/Konrad for the suspended load.  
            if counter < MaxSteps-1:
                if self.inputs.FeedType == 'RatingCurveMiddleElwha' and subdaycount == Tmult - 1 and self.inputs.Qlist[counter + 1]>40.:
                #print self.inputs.Qlist                                  
                    Node = Reach.Node[0]
                    NewAvkFeed, NewJKFeed = \
                        Node.Load.UpdateElwhaFeedRatingCurve(ControlNode, LoadFactor[n],\
                        self.inputs.MudFraction, self.inputs.Qlist[counter + 1],\
                        self.inputs.vfunc, self.inputs.TrinityFit, CalibrationFactor, ControlGSD, Section = 'Middle')
                        
            if counter < MaxSteps-1:
                if self.inputs.FeedType == 'RatingCurveMiddleElwhaStochastic' and subdaycount == Tmult - 1 and self.inputs.Qlist[counter + 1]>40.:
                #print self.inputs.Qlist
                    multiplier = np.random.lognormal(0., .9)
                    #print multiplier                                  
                    Node = Reach.Node[0]
                    NewAvkFeed, NewJKFeed = \
                        Node.Load.UpdateElwhaFeedRatingCurve(ControlNode, multiplier,\
                        self.inputs.MudFraction, self.inputs.Qlist[counter + 1],\
                        self.inputs.vfunc, self.inputs.TrinityFit, CalibrationFactor, ControlGSD, Section = 'Middle')

            if counter < MaxSteps-1:      
                if self.inputs.FeedType == 'RatingCurveUpperElwha' and subdaycount == Tmult - 1:
                    Node = Reach.Node[0]
                    NewAvkFeed, NewJKFeed = \
                        Node.Load.UpdateElwhaFeedRatingCurve(ControlNode, LoadFactor[n],\
                        self.inputs.MudFraction, self.inputs.Qlist[counter + 1],\
                        self.inputs.vfunc, self.inputs.TrinityFit, CalibrationFactor, ControlGSD, Section = 'Upper')

            # ******************End ELWHA RATING CURVE ************************
            
            # *********Example: ELWHA DAM REMOVAL:  DURATION CURVE ************
            # This is an example of a customized method of updating the sediment
            # feed. 
             
            # Simulates exponential decay of the LoadFactor variable and was 
            # calibrated to data provided by Tim Randle.  It is written to be
            # used with a flow duration curve and activated partway through
            # a run (see Activate Dam Removal, Section VI.D).
            if self.inputs.FeedType == 'DamRemovalDuration': 
                Node = Reach.Node[0]                
                C = self.inputs.C1
                tau = self.inputs.tau
                t = (Tyear-DamTyearStart)*365.25 # In days
                PercSus = self.inputs.PercSus
                multiplier = ''
                    
                if self.inputs.DecayType == 'Exp':
                    multiplier = 1 + C*exp(-t/tau)

                NewAvkFeed, NewJKFeed = \
                    Node.Load.UpdateElwhaDamRemovalDurationCurve(ControlNode,\
                    multiplier)                
            # ************End ELWHA DAM REMOVAL:  DURATION CURVE **************
                    
            # *********Example: ELWHA DAM REMOVAL:  RATING CURVE **************
            # This is an example of a customized method of updating the sediment
            # feed. 
             
            # Simulates exponential decay of the LoadFactor variable and was 
            # calibrated to data provided by Tim Randle.  It is written to be
            # used with a hydrograph and activated partway through
            # a run (see Activate Dam Removal, Section VI.D).  The feed is determined
            # by the daily flow, and the multiplier by the decay function.
            if counter < MaxSteps-1: 
                if self.inputs.FeedType == 'DamRemovalRatingCurve' and subdaycount == Tmult - 1 and self.inputs.Qlist[counter + 1]>40.:                
#                    DamTyearStart = Tyear
                    Node = Reach.Node[0]
                    C = self.inputs.C2
                    tau = self.inputs.tau
                    t = (Tyear-DamTyearStart)*365.25
                    multiplier = ''
                    FixedCapacity = False
                    
                    if self.inputs.DecayType == 'Exp':
                        multiplier = 1 + C*exp(-t/tau)
                    elif self.inputs.DecayType == 'Pow':
                        v = 500 # Hard coded here for now--can make it an input variable later.
                        multiplier = 1 + C*(t+1)**(-v/tau)

                    if self.inputs.FixedCapacity == True: # Determines whether feed is based on duration-averaged capacity or capacity for daily discharge
                        FixedCapacity = True

                    NewAvkFeed, NewJKFeed = \
                        Node.Load.UpdateElwhaFeedRatingCurve(ControlNode,\
                        multiplier, self.inputs.MudFraction,\
                        self.inputs.Qlist[counter + 1], self.inputs.vfunc,\
                        self.inputs.TrinityFit, CalibrationFactor, ControlGSD,\
                        Removal = True, PercSus = self.inputs.PercSus, FixedCapacity = FixedCapacity)

                    feedcounter = counter - startcounter
                    
            if counter < MaxSteps-1:         
                if self.inputs.FeedType == 'DamRemovalFirstPulseRatingCurve' and subdaycount == Tmult - 1 and self.inputs.Qlist[counter + 1]>40.:                
#                    DamTyearStart = Tyear
                    Node = Reach.Node[0]
                    C = self.inputs.C1
                    tau = self.inputs.tau
                    t = (Tyear-DamTyearStart)*365.25 # In days
                    PercSus = self.inputs.PercSus
                    multiplier = ''
                    FixedCapacity = False
                    
                    if self.inputs.DecayType == 'Exp':
                        multiplier = 1 + C*exp(-t/tau)
                    elif self.inputs.DecayType == 'Pow':
                        v = 500 # Hard coded here for now--can make it an input variable later.
                        multiplier = 1 + C*(t+1)**(-v/tau)

                    if self.inputs.FixedCapacity == True: # Determines whether feed is based on duration-averaged capacity or capacity for daily discharge
                        FixedCapacity = True
                    #print date
                    #print multiplier    
                    NewAvkFeed, NewJKFeed = \
                        Node.Load.UpdateElwhaFeedFineRatingCurve(ControlNode,\
                        multiplier, self.inputs.MudFraction,\
                        self.inputs.Qlist[counter + 1], self.inputs.vfunc,\
                        self.inputs.TrinityFit, CalibrationFactor, ControlGSD, ControlGSD,\
                        Removal = True, PercSus = PercSus, FixedCapacity = FixedCapacity, FirstPulse = True)
                    #print sum(NewAvkFeed)*60*60*24
                    feedcounter = counter - startcounter
                
            # ************End ELWHA DAM REMOVAL:  Rating CURVE ****************
            
            #************************************************
            #VI.F:  Update the downstream boundary condition if it is user-defined.

            BoundaryGauge = ""
            if BoundaryType == 'counter':
                BoundaryGauge = counter
            elif BoundaryType == 'date':
                BoundaryGauge = date
            
            if n < len(BoundaryFactorCount):
                if BoundaryGauge == BoundaryFactorCount[NextBoundary]:
                    Reach.SetBoundary = True
                    Reach.BoundaryHc = BoundaryFactor[NextBoundary]                
                    if NextBoundary < len(BoundaryFactorCount) - 1: NextBoundary += 1
                        
            #************************************************
            #VI.G:  Adjust lateral exchange parameters.  This section of code is a 
            #relic from a previous version and is a way to adjust the migration 
            #rate as a function of the sediment supply rate.  It is currently
            #commented out but could be used as an alternative to the width change
            #function.         

            # In this section, the bank migration rate, width, etc can be adjusted.
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
            
            #************************************************
            #VI.H:  Update time trackers

            #  Print record of time--every 10 timesteps for duration curve and 
            #  every month for hydrographs
            if date.day == 1 and subdaycount == 0 and self.inputs.Hydrograph == True and self.inputs.CyclingHydrograph == False:
                print(str(date))
            if self.inputs.Hydrograph == False and counter % 10 == 0:
                print(counter)
            
            #  Change the length of the timestep if at the user-specified point in the model            
            if len(dtcount) > 0 and counter == dtcount[m]:
                dt = dtlist[m+1]* 365.25*24*60*60
                if m < len(dtcount) - 1:
                    m = m + 1
            
            subdaycount = subdaycount + 1

            #  Counter only updates for a new day if hydrograph is on.  If using
            #  a duration curve, it updates at each time step.
            if subdaycount == Tmult or self.inputs.Hydrograph == False:
                counter = counter + 1
                date = date + datetime.timedelta(days=1)
                
            # Filler if you want option to cycle endlessly thorough hydrograph.
            # Not currently connected to anything
            #filler = False 
                
            if self.inputs.CyclingHydrograph == True and counter == MaxSteps: # Cycles endlessly through the hydrograph.
                counter = 0
                NextCount = Interval
        
        #**********************************************************************
        #VII:  Save final output

        #  Save final objects
#        pickle.dump(Reach, open(os.path.join(os.pardir, self.inputs.Outputfolder, "save.Reach"), "wb" )) 
#        pickle.dump(self.inputs, open(os.path.join(os.pardir, self.inputs.Outputfolder, "inputparams.inputs"), "wb"))         
        
        #  Save objects with daily record (for hydrographs)
        OutputObj.WriteDailyFiles()

        # Produce done file
        donefile = os.path.join(os.pardir, self.inputs.Outputfolder, 'Model_finished')
        with open (donefile, 'w') as f:
            f.write('Model completed')
            f.close()		
           


