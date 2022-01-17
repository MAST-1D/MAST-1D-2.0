# -*- coding: utf-8 -*-
"""
Created on Mon Feb 23 15:44:47 2015

@author: Kathryn De Rego, UBC Geography
"""


"""
REQUIRED PACKAGES
User does not modify.
"""

import os
import sys
#import csv
import datetime
import numpy as np
import matplotlib.pyplot as plt
import pickle


sys.path.append("..")
from Hydrology.clsTimeSeries import clsTimeSeries
#from MAST_1D.clsModel import clsModel
#from MAST_1D.clsInputs import clsInputs
#from MAST_1D.clsOutputWriter import clsOutputWriter
from clsInputs import clsInputs
from clsOutputWriter import clsOutputWriter
from clsReach import clsReach
from clsTracerProperties import clsTracerProperties
#from clsModel import clsModel
from math import exp
from copy import deepcopy


def load_hydrographs(inputs):
    #Initialize lists
    inputs.Qlist = []
    inputs.Qw = []
    inputs.p = []
    
    for f in inputs.DischargeFiles:
        # Load discharge files as a lists
        DischargeFile = os.path.join(os.pardir,"Discharge_Files", f)
        print(DischargeFile)
        Qlist = open(DischargeFile).readlines()
        Qlist = list(map(lambda x: float(x), Qlist))
        binsize = (max(Qlist)-min(Qlist))/inputs.NumberofHydroBins
        print('binsize = %s' %(binsize))
        # Create a duration curve from the list for setting up equilibrium floodplain
        # conditions and feed.  Other parameters can be customized (see ExtractDC function).
        Q = clsTimeSeries([],Qlist) 
        Qw, p = Q.CreateDurationCurve(binsize)
        print('number of bins made = %s' % len(Qw))
        inputs.Qlist.append(Qlist)
        inputs.Qw.append(Qw)
        inputs.p.append(p)
      
    return inputs
               
inputs = clsInputs() 

#**********************************BEGIN INPUTS*****************************************

"""
I.A  IMPORT CSV WITH NODE-SPECIFIC INPUTS FOR SOME ATTRIBUTES

Here you can import a .csv file with spacially-explicit initial conditions
"""

#inputfile = os.path.join(os.pardir, 'Elwha_spatial_data', 'UR_Valley_Nodes_sinuosity.csv')
#
## open the file in universal line ending mode 
#with open(inputfile, 'rU') as infile: # This loop copied from Stack Overflow.
#  # read the file as a dictionary for each row ({header : value})
#  reader = csv.DictReader(infile)
#  data = {}
#  for row in reader:
#    for header, value in row.items():
#      try:
#        data[header].append(value)
#      except KeyError:
#        data[header] = [value]
#        
#inputs.BcNodes = map(lambda x: float(x), data['InitialBc'])
#inputs.BfNodes = map(lambda x: float(x), data['Avg_width'])
#inputs.dxf = map(lambda x: float(x), data['Valley_length']) # Optional length of valley to determine node length with sinuosity
#inputs.Canyon = map(lambda x: bool(int(x)), data['Canyon'])
#inputs.ReachwideBedrock = False # If true, every reach will only be allowed to degrade a user-specified amount and bed will become 'partly-alluvial'
#inputs.PartlyAlluvialMin = -.2
#inputs.ChSinNodes = map(lambda x: float(x), data['Avg_Sinuosity'])

"""
I.B  IMPORT REACH CONDITIONS FROM A PREVIOUS RUN
"""
inputs.initialcond = False # If you are using a prior run as initial conditions, True; if starting from scratch, False
inputs.priorReach = 'D:\MAST-1D_version_K6\Output\Pre_vs_post_dam\WholePeriod'+ '//' +'save.Reach' # File with initial conditions Reach object

"""
II:  SET CHANNEL GEOMETRY.

Note that you must fill these variables even if you are importing them in Section
I.A or I.B.  MAST-1D will set the feed and floodplain number based on these values.
"""
inputs.Nnodes = 66 #  Number of nodes
inputs.Bc = 94. # Channel width (m)
inputs.reachlength = 13673. # Length of reach (channel length) (m)
inputs.Bf = 500 # Total valley width (m)
inputs.ChSin = 1.02 # Channel sinuosity
inputs.Slope = 0.0074 # Bed gradient

"""
III:  CHOOSE HYDROGRAPH OR FLOW DURATION CURVE
"""
inputs.Hydrograph = True # True if supplying hydrograph instead of flow duration curve--not used for this demo
inputs.CyclingHydrograph = False #Flag to determine if hydrograph is repeated
"""
IV.A  TIMESTEP PARAMETERS FOR FLOW DURATION CURVE
"""    
inputs.MaxSteps = 200 #  Total number of timesteps to run.
# for hydrograph runs, ends at smaller of maxsteps or the number of days in the hydrograph.
inputs.dt = [0] # Timestep (years)--can add multiple values for timestep adjustment during run
inputs.dtcount = [] # List of timesteps to instigate timestep interval change (length of list should be one fewer than dt)

"""
IV.B  TIMESTEP PARAMETERS FOR HYDROGRAPH

MAST-1D is set up to run hydrographs with a daily resolution.  It will keep track
of the date and instigate user-specified boundary condition changes based on the date.
"""
inputs.Tmult = 1 # Number of timesteps per day
inputs.TmultHighFlow = 8 #150 # Number of timesteps for days with flow above High Flow Timestep Threshold
inputs.TmultCyclingHydrograph = 27 #Number of timesteps per day for Cycling Hydrograph runs.
inputs.StartDate = (2011,9,15) # (year, month, day).
inputs.LowFlowThreshold = 40 # threshold discharge below which hydrograph runs ignore sediment transport
inputs.HighFlowTimestepThreshold = 120 #300 # threshold discharge above which hydrograph runs use a shorter timestep.


"""
V:  BOUNDARY CONDITIONS
"""
inputs.Removal = False # This is a custom variable that is used to set some boundary
    # condition behavior for the Elwha River Dam removal (see Section VI.D in 
    # clsModel).

#  Upstream sediment feed
inputs.LoadFactor = [0.8,0.8]  #  List of upstream sediment feed (in fraction of equilibrium capacity)
inputs.LoadFactorCount = [900] #  List of times to instigate feed change--should
    # be an int (# of timesteps) for duration curves or a tuple ((yyyy, m, dd))
    # for hydrographs.  Number of entries should be one less than LoadFactor.
inputs.FeedType = "RatingCurveUpperElwha" # Choose 'DurationCurve' if using a duration 
    # curve or 'RatingCurve' if using a hydrograph.  You can also create a custom
    # method for applying feed--see Section VI.E in clsModel. Note that for many runs,
    # supply is specified at the capacity of the upstream node, UNSURE IF USED.

#  Downstream water surface elevation
inputs.SetBoundary = False # True if you want to set a downstream boundary condition; false, model calculates it
inputs.BoundaryFactor = [30.]# The downstream WSE--is the same for all flows
inputs.BoundaryFactorCount = [] #  List of times to instigate WSE change--should
    # be an int (# of timesteps) for duration curves or a tuple ((yyyy, m, dd))
    # for hydrographs.  Number of entries should be one less than BoundaryFactor.

"""
VI:  HYDRAULICS
"""
inputs.ManningStabilizer = 3 # Controls how much the water surface elevation changes
    #with each iteration in the backwater calculation
inputs.vfunc = True # Velocity function, True for Manning, False for Chezy
inputs.Cfc = 0.0075 # 1/Cz^2 for channel
inputs.Cff = 0.69 # 1/Cz^2 for floodplain
inputs.ncAddons = 0.0066 # Form roughness component for Manning's n; 
    # will be added to calculated n value
inputs.ncMultiplier = 1.15 # Sinuosity multiplier for Manning's n
inputs.nf = 0.1 # Manning's n for floodplain

"""
VII.A:  DISCHARGE RECORD FOR FLOW DURATION CURVE

Note that it is not clear if this is used for the hydrograph runs.
"""
# List of discharges in flow duration curve (m^3/s)
# This distribution was created from the Elwha River near McDonald Bridge, 1927-1989
inputs.Qw = [\
30.,\
50.,\
70.,\
90.,\
150.,\
275.,\
455.,\
540.] 

# List of flow frequencies for flow duration curve
inputs.p = [\
0.3,\
0.5,\
0.1,\
0.05,\
0.045,\
0.004,\
0.0007,\
0.0003]

"""
VII.B:  DISCHARGE RECORD FOR HYDROGRAPH
"""
inputs.DischargeFileID = [0] # ID of each of the discharge timeseries.  Each
#inputs.DischargeFileID = [0,1] # ID of each of the discharge timeseries.  Each
    # node is given this id, which correspond with the given timeseries files.
inputs.DischargeFileCoords = [0] #downstream coordinate at upper end of a
    # given discharge timeseries
inputs.DischargeFiles = ['Qfile15Sep2011pres'] # Names 
    # of discharge files.  They should be
    # stored in the "Discharge_Files" folder in the parent directory. 
    #  See examples there for formatting.   
#inputs.DischargeFile = inputs.DischargeFiles[0] # Names of discharge file.  It should be
    # stored in the "Discharge_Files" folder in the parent directory.  See examples
    # there for formatting.
#inputs.HydroBins = 2.832 # The increment for bins for the flow duration curve.
    # MAST-1D calculates initial floodplain and feed parameters based on a flow
    # duration curve representing the modeled period, even when the hydrograph 
    # function is turned on.  A flow duration curve will be calculated automatically
    # using DischargeFile and the bin increment.  A good initial rule of thumb 
    # for the bin increment is the resolution of the discharge record.
inputs.NumberofHydroBins = 15
    # Number of bins for all flow duration curves.  Note that this number
    # cannot be too large if multiple 
"""
VIII:  GRAINSIZE AND SEDIMENT TRANSPORT
"""
# Bounds of grainsize classes (mm).  The finest size class should be the silt/clay
# class.  There must be a lower bound (here it is estimated to be 0.002).
inputs.Dbdy = [\
0.002,\
.063,\
1.0,\
2.0,\
4.0,\
8.0,\
16.0,\
32.0,\
64.0,\
128.0,\
256.0,\
512.0,\
1024.]

# List of size frequencies for Active Layer (can be in any unit and does not
# need to add up to 1 or 100--MAST-1D will normalize)
inputs.Fa = [
0.,\
7.0,\
4.3,\
4.7,\
4.7,\
5.7,\
7.4,\
15.,\
21.,\
20.,\
9.7,\
1.]

# Substrate GSD alteration--to add lag deposits to the substrate, choose a fraction
# for the grainsize of choice.  If the lists are left blank, the substrate is 
# composed of the same material as the Active Layer and Floodplain.
inputs.DLag = [] # List of indexes of grainsizes to alter
inputs.FLag = [] # Fraction of GSD

# Parameters for suspended load
inputs.FSandSusp = 0.12 # Fraction of sand in the suspension.  This parameter is 
    # not currently connected in the model.
inputs.MudFraction = 18. # Mud feed multiplier (multiple of next finest size class)
inputs.FlBed = .75 # Floodplain number for bed material

# Bedload sediment transport equation
inputs.TransFunc = 'WilcockCrowe' # Transport function; choices are 'WrightParker' and 'WilcockCrowe'
inputs.TrinityFit = True # True is Trinity River form (Gaeuman 2009), False is normal Wilcock and Crowe        
inputs.CalibrationFactor = 1.  # Calibration coefficient for critical shear stress

# Lateral Sediment Supply
#inputs.LateralSedimentSourceNodes = [20,60]  # indices of nodes with lateral supply
#inputs.LateralMultiplier = [0.1,0.1]  # fraction of upstream supply for lateral source nodes

inputs.LateralSedimentSourceNodes = []  # indices of nodes with lateral supply
inputs.LateralMultiplier = []  # fraction of upstream supply for lateral source nodes

"""
VIII. CHANNEL MIGRATION AND WIDTH CHANGE
"""
inputs.WidthChange = True # If true, turns off constant migration rate and 
    # channel-floodplain coupling is determined by width change functions
inputs.migration = 1.2 # Channel migration rate (m/yr).  Estimate a value even if 
    # the width change function is on; it is used to calculate the equilibrium 
    # floodplain number.
    
# Vegetation encroachment terms
inputs.BcMin = 40.  # Minimum channel width
inputs.W = .07 # Narrowing constant--used to calibrate narrowing function.  
    # As a starting guide, estimate the percentage of bar that is vegetated
    # annually and double it.
inputs.alphatau = 32. # Shear stress threshold for channel narrowing

# Bank erosion terms
inputs.ErodeT = .55 # Fraction of near-bank sediment sourced from the active layer
inputs.MobilityThreshold = 10**-7 # Mobility threshold for initiation of bank erosion

# Avulsion terms
inputs.AvulsionExchange = .1

inputs.AvulsionThreshold = -1.25 # Minimum bank height below which avulsion will occur 
# can be set negative to prevent all avulsions.

"""
X. RESERVOIR THICKNESSES AND EXCHANGE PARAMETERS
"""
# Reservoir characteristics
inputs.FloodplainL = 1.75 # Initial thickness of the active floodplain (m)
inputs.ActiveLayerL = .4 # Initial thickness of the Active Layer
inputs.LayerL = 1.5 # Thickness of substrate layers (m)
inputs.NLayers = 2 # Number of substrate layers
inputs.Hpb = 1.7 # Thickness of the point bar (constant through time)
inputs.lambdap = 0.5 # Porosity

# Reservoir exchange parameters
inputs.Kbar = 1/1000000. # Parameter controlling fraction washload in point bar deposits
inputs.AlphaBed = .9 # Proportion of new substrate composed of active layer material
    # (verses load material)
inputs.AlphaBar = 1. # Parameter controlling similarity between bed material load 
    # and bar deposition
inputs.AlphaPartlyAlluvial = 0.9 # Parameter controlling similarity between bed
    # material load and deposition in the active layer of a partly alluvial node

"""
XI. TRACER PROPERTIES 

The tracer component has not been tested in this version of MAST-1D, but should 
theoretically work.  They may need to be initialized in order to compute a run.
"""
# Set up as Cosmogenic 14C.
inputs.NTracers = 1 # Number of tracers
inputs.coj = [82.96, 1.98, 15.06] # production from different processes (at surface, presumably--not integrated over depth)
inputs.Lcj = [160., 738., 2688.] # Attenuation rates (g/cm^2)
inputs.Name = "'14C'"  # Name of radioisotopic tracer
inputs.DecayConst = 0.000121 # Decay constant (can be zero for a cnservative tracer)
inputs.ProductionRate = 15.1 # Nuclide production rate
inputs.FalloutRate = 0. # Nuclide fallout rate--most relevant for fine sediment.

inputs.TracerICFloodplain = 1.
inputs.TracerICActiveLayer = 0.
inputs.TracerICSubstrate = 0.
inputs.TracerBCFeed = 0.
        
"""
XII. OUTPUT SPECIFICATIONS 

There are three types of output:
A. Text files of node attributes for each node over time periods of equal intervals
    (for hydrographs and duration curves)
B. Text files of node attributes for each node on user-specified dates (currently 
    for hydrographs only)
C. JSON files with daily output at a given location. User specifies which 
    nodes and variables to write out. For hydrograph runs only.  
"""
# Output folder
RunName = "ElwhaHydrographDemo_Simple_short" # Name of the run (a folder of outputs files will be created
    # under this name).
inputs.Outputfolder = os.path.join("Output","Demo",RunName) # Parent folder for 
    # the output folder 

# A.  General output

# Number of dataslices to save.  Will be written at equal time intervals.
inputs.NumberOfPrintouts = 50 
 
# List variables to save at regular intervals (strings of clsNode attributes)
inputs.Outputvars = ['Slope',
                     'Bf',
                     'Bc',
                     'DC.Qwf[-1]',
                     'DC.Sf[0]',
                     'DC.Hc[0]',
                     'DC.Hc[-1]',
                     'DC.Uc[-1]',
                     'DC.Uc[0]',
                     'DC.WSE[0]',
                     'DC.WSE[-1]',
                     'DC.Qw[0]',
                     'DC.Uc[0]',
                     'ActiveLayer.GSD.D50',
                     'ActiveLayer.GSD.D84',
                     'Load.QsavBedTot',
                     'Load.QsavBedTotFeed',
                     'Floodplain.GSD.D50',
                     'etabav',
                     'Substrate[-1].C.GSD.D50',
                     'CumulativeBedChange',
                     'Floodplain.L',
                     'Floodplain.GSD.F[0]',
                     'Floodplain.GSD.D84',
                     'CumulativeTotBedMaterialFeed',
                     'SLatSourceAv[-1]',
                     'CobbleMobility',
                     'WidenRate',
                     'NarrowRate',
                     'PointBarSubsurfaceGSD.D50',
                     ]

# List variables for daily output (strings of clsReach attributes)    
inputs.DailyOutputVars=['Node[0].DC.Qw[0]',
                        'Node[-1].DC.Qw[0]',
                        'Node[10].Load.QsjTot[0]',
                        'Node[0].CumulativeTotFeed',
                        'Node[-1].CumulativeTotFeed',
                        'Node[10].Load.QsAvkLoad[0]',
                        'Node[10].Load.QsAvkLoad[7]',
                        'Node[10].ActiveLayer.GSD.D50',
                        'Node[10].ActiveLayer.ExSed.InWidthChange[1]',
                        'Node[10].ActiveLayer.ExSed.InWidthChange[7]',
                        'Node[10].Bc',
                        'Node[10].NarrowRate',
                        'Node[10].WidenRate',
                        'Node[10].CumulativeNarrowing',
                        'Node[10].CumulativeWidening', 
                        'Node[10].ActiveLayer.GSD.F[0]',
                        'Node[10].ActiveLayer.GSD.F[1]',
                        'Node[10].ActiveLayer.GSD.F[7]',
                        'Node[10].ActiveLayer.GSD.F[-1]',
                        'Node[9].ActiveLayer.T[1,0]',
                        'Node[10].ActiveLayer.T[1,0]',
                        'Node[-1].ActiveLayer.T[1,0]',
                        'Node[-1].ActiveLayer.T[7,0]',
                        'Node[10].ActiveLayer.T[7,0]',
                        'Node[0].Floodplain.T[1,0]',
                        'Node[0].Load.TBedFeedk[0,0]',
                        'Node[1].Load.TBedFeedk[0,0]',
                        'Node[1].ActiveLayer.T[0,0]',
                        'Node[0].ActiveLayer.T[0,0]',
                        'Node[0].ActiveLayer.T[1,0]',
                        'Node[10].FkPointBarDepositAnnualAverage[1]',
                        'Node[10].Dfav[0]',
                        'Node[10].Dfav[1]'
                        ]
    
# B. Output on specific dates (for hydrograph runs)

# Dates (yyyy, m, dd) in which to output variables for model validation 
inputs.ValidateDates = [(2011, 9, 30),
                        (2012, 9, 30),
                        (2013, 9, 30),
                        (2014, 9, 30),
                        (2015, 9, 30),
                        (2016, 9, 30)] 
    
# Variables (attributes of clsNode) in which to output on specific dates for model validation
inputs.Validatevars = ['CumulativeTotFeed',
                       'Bc',
                       'CumulativeNarrowing',
                       'CumulativeWidening',
                       'CumulativeTotalAvulsionWidth']

# C. Output daily (for hydrograph runs)

# Nodes that will be output daily.  The variables to save are hard-coded in
# clsOutputSpecs. Note that the file sizes of daily
# output can be up to several megabytes, so use judiciosly.  It is better to use
# inputs.DailyOutputVars to specify daily output.
inputs.DailyNodes = [0,35] 


#**********************************END INPUTS*******************************************


#  Checks to see if the specified output folder exists and creates it if it doesn't  
directory = str(os.path.join(os.pardir,inputs.Outputfolder))
if not os.path.exists(directory):
    os.makedirs(directory)

# Creates a duration curve for the discharge record for a hydrograph run and
# creates a model object with the proper inputs.
run = ""
if inputs.Hydrograph == True:
    newinputs = load_hydrographs(inputs)
    #run = clsModel(newinputs)
    #inputs = inputs
#else:
    #run = clsModel(inputs)
    #inputs = inputs
    
#**********************************************************************
#Section III:  Set initial conditions.
#Initialize the reach, control node, grain size distribution, boundaries, tracers, and calibration factor

def CustomNodes(Reach):
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
        if hasattr(inputs, 'BcNodes'):
        
        #  If the variable exists, it is changed to the node specific value.
        #  (Note:  if it doesn't exist, the value will keep the non-node-specific
        #  value also provided in the inputs.)
            Node.Bc = inputs.BcNodes[i]
            
        #  Repeat for all 'custom' variables
        #  Here, BfNodes is the valley width and ChSinNodes is sinuosity.
        if hasattr(inputs, 'BfNodes'):
            Node.Bf = inputs.BfNodes[i]-Node.Bc  
            Node.ValleyMargin = inputs.BfNodes[i]
        if hasattr(inputs, 'ChSinNodes'):
            Node.ChSin = inputs.ChSinNodes[i]
        if hasattr(inputs, 'SlopeNodes'):
            Node.Slope = inputs.SlopeNodes[i]
            
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
        if hasattr(inputs, 'dxf'):
            Node.dxc = inputs.dxf[i]*Node.ChSin
        else:
            Node.dxc = inputs.reachlength / Reach.nnodes()
        
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
           
        if hasattr(inputs, 'Canyon'):            
            if inputs.Canyon[i] == True:
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
    if inputs.SetBoundary==True:
        for j in range(Reach.NFlows):
            Reach.Node[-1].DC.WSE[j] = Reach.Node[-1].etabav + inputs.BoundaryFactor[0]
    


def SetInitialConditions(inputs):
    """
    Section I:  Load the Reach ipnuts 
    and perform initial hydraulic and sediment transport computations to 
    fully define the initial condition.
    
    Parameters
    ----------
    inputs : :obj:`MAST_1D.clsInputs`
        The class defining all input parameters.
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
    CalibrationFactor = inputs.CalibrationFactor # Sets sediment transport calibration

    #  Sets initial tracer conditions.  These functions are untested in 
    #  this version of MAST-1D.
    TracerProperties = [clsTracerProperties()]
    TracerProperties[0] = Reach.SetupTracers(inputs, TracerProperties[0])

    #  Set the initial hydraulic and sediment transport conditions
    Reach.UpdateSlope()
    Reach.UpdateManningDepthAtAllFlowAndNodes()
    for Node in Reach.Node:
        Node.Load.UpdateSedimentLoadByDurationAndSize(Node.DC,
                                                      Node.ActiveLayer.GSD,
                                                      1000.,
                                                      2.7,
                                                      Node.Bc,
                                                      Node.FractionAlluvial,
                                                      inputs.TransFunc,
                                                      inputs.TrinityFit,
                                                      CalibrationFactor) 
           
    Reach.StepDownstream(0.,
                         inputs.AlphaBed,
                         0.2,
                         0.1,
                         inputs.LayerL,
                         TracerProperties,
                         inputs.TransFunc,
                         inputs.TrinityFit,
                         CalibrationFactor,
                         False,
                         inputs.W,
                         inputs.ErodeT,
                         inputs.alphatau,
                         inputs.BcMin,
                         inputs.vfunc,
                         inputs.AvulsionThreshold,
                         Reach.Node[0].ActiveLayer.GSD,
                         inputs.AvulsionExchange,
                         inputs.MobilityThreshold) # Katie: may want to change other parameters to variables.  Width change variables are optional.
    
    #  Preserve control nodes to be used to set boundary conditions later
    
    ControlNode = deepcopy(Reach.Node[0]) # For calculating feed from a rating curve (see update feed below)
    ControlGSD = deepcopy(Reach.Node[0].ActiveLayer.GSD)
    # Note that for now, there is only one control node representing the upstream
    # end of the system.  However, multiple control nodes could be used
    # to compute a sediment supply from tributaries.  If this approach is 
    # taken, at a minimum, the parameters in the control notes used for
    # computing sediment transport capacity would need to be defined.  A
    # non-complete list of such parameters includes width, depth, slope,
    # grain size distribution, channel roughness parameters, and all the
    # flags related to the choice of roughness formula and sediment
    # transport equation. An alternative approach for determining lateral
    # sediment supply would simply be to specify a size-specific supply
    # based on a rating curve.  This would need to be done at each time 
    # step and is not implemented.
    
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
        CustomNodes(Reach)
        for j in range(Reach.NFlows):
            zControlBoundary.DC.WSE[j] = Reach.Node[-1].DC.WSE[j]
        # Reset controlNode bed elevation
        zControlBoundary.etabav = Reach.Node[-1].etabav
        #print Reach.Node[-1].etabav
                    
    return Reach, ControlNode, ControlGSD, zControlBoundary, TracerProperties, CalibrationFactor

Reach, ControlNode, ControlGSD, zControlBoundary, TracerProperties, CalibrationFactor = \
    SetInitialConditions(inputs)
    
    
    
    
#Set up the hydrograph
if inputs.Hydrograph == True:
    Reach.set_up_hydrograph()


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
LoadFactor = inputs.LoadFactor
LoadFactorCount = inputs.LoadFactorCount
LoadType = ""
if type(inputs.LoadFactorCount[0]) == int:
    LoadFactorCount = inputs.LoadFactorCount
    LoadType = 'counter'
if type(inputs.LoadFactorCount[0]) == tuple:
    LoadFactorCount = map(lambda x: datetime.date(*x[:6]),\
        inputs.LoadFactorCount)
    LoadType = 'date'
n = 0 # Keeps track of which item in LoadFactorCount is the current trigger
NewAvkFeed = Reach.Node[0].Load.QsAvkFeed #Katie add: Updated feed goes here, and feed in Node objects are filled after results are written to file
NewJKFeed = Reach.Node[0].Load.Qsjkfeed

#  Boundary change can be determined from either the counter or date, like
#  the sediment load (see above)        
BoundaryFactor = inputs.BoundaryFactor
BoundaryFactorCount = inputs.BoundaryFactorCount
BoundaryType = ""
if inputs.Hydrograph == False:
    BoundaryType = 'counter'
else:
    BoundaryFactorCount = list(map(lambda x: datetime.date(*x[:6]),\
        inputs.BoundaryFactorCount))
    BoundaryType = 'date'
NextBoundary = 0 #  Keeps track of which item in BoundaryChangeCount is the current trigger      

#  This customized global variable is needed to trigger a more sophisticated 
#  boundary condition change in Section VI. 
BoundaryTrigger = False # Trigger that determines whether nodes are added in a Dam Removal situation      
NodeTrigger = False        

#Define a function to plot a size distribution
def plot_GSD(GSD,graph_title,exclude_bins=None):  
    FF = np.zeros(len(GSD.F))
    fig, ax = plt.subplots(1,1)
    labels = []
    for i in range(len(FF)):
        labels.append(f'{np.format_float_positional(GSD.D_lower[i])} - {np.format_float_positional(GSD.D_upper[i])}')
        if exclude_bins is not None and i in exclude_bins:
            FF[i]=0
        else:
            FF[i]=GSD.F[i]
    FF = FF/sum(FF)
    x_pos = [i for i, _ in enumerate(FF)]
    ax.bar(x_pos, FF, color='green')
    plt.xticks(rotation = 90)
    plt.xlabel("Grain Size Bin (mm)")
    plt.ylabel("Fraction")
    plt.title(graph_title)
    plt.xticks(x_pos, labels)
    return fig

#plot the size distribution of the bedload in the upstream node
fig4 = plot_GSD(Reach.Node[0].Load.GSDBedloadAv,
                "Bedload Size Distribution At Upstream Node")
fig5 = plot_GSD(Reach.Node[0].Load.GSDBedloadAv,
                "Gravel Portion of Bedload Distribution At Upstream Node",
                exclude_bins=[1])

#plot the size distribution of the active layer in the control node
fig6 = plot_GSD(ControlNode.ActiveLayer.GSD, 
                "Bed Material Distribution in Control Node")


#HERE is where one could update a bunch of geometric or node parameters

def plot_profile(x,y,graph_title):  
    fig, ax = plt.subplots(1,1)
    ax.plot(x, y)
    plt.xticks(rotation = 90)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(graph_title)
    return fig

x = [Node.xc for Node in Reach.Node]
y = [Node.etabav for Node in Reach.Node]
fig7 = plot_profile(x,y,"Longitudinal Profile")

y = [Node.ActiveLayer.GSD.D50 for Node in Reach.Node]
fig8 = plot_profile(x,y,"Grain Size Profile")

#**********************************************************************
#Section V:  Define time-keeping counters and set up output


#  Timesteps elapsed for number flow duration curves or the number of 
#  days elapsed for hydrographs
counter = 0 

#  This is a custom counter for monitoring dam removal feed
startcounter = 0

#  The date.
date = datetime.date(*inputs.StartDate)

#  The time elapsed in years.
Tyear = 0

#  Set the length of the timestep and determine when it changes
#  when a flow duration curve is used
dt = inputs.dt[0] * 365.25*24*60*60 # initial timestep in seconds
dtlist = inputs.dt # list of timesteps
dtcount = inputs.dtcount # list of when timestep changes occur
m = 0 # keeps track of which integer in dtcount is the current trigger
MaxSteps = inputs.MaxSteps # Time counter--Total number of timesteps

#  Set the length of the timestep and determine when it changes
#  when a hydrograph is used
Tmult = inputs.Tmult # Number of timesteps in a day
subdaycount = 0 # Timestep counter within a day
if inputs.Hydrograph == True:
    dt = 1./365.25/Tmult*365.25*24*60*60 # In seconds
    MaxSteps = min(len(inputs.Qlist[0]),inputs.MaxSteps) # Time counter--Total number of days
LowFlow = False # Sets low flow tag

try:
    inputs.CyclingHydrograph
except:
    inputs.CyclingHydrograph = False
if inputs.CyclingHydrograph == True:
    dt = 1./365.25/27*365.25*24*60*60
    Reach.Node[-1].Canyon = True

#  Set up output--note that the variables saved are specified in clsOutputSpecs
OutputObj = clsOutputWriter(inputs.Outputfolder, inputs.DailyNodes, date)

#  'Standard' output files that are printed at consistent time intervals
NumberOfPrintouts = inputs.NumberOfPrintouts # Number of times during run to record output        
Interval = MaxSteps / NumberOfPrintouts # Interval in time counters between output recordings
Printstep = 0 # Locator for standard output files
NextCount = 0 # Keeps track of what value of counter will trigger the next standard print interval

#  Output files that are written at user-defined dates for model performance testing      
ValidateDates = list(map(lambda x: datetime.date(*x[:6]), inputs.ValidateDates)) # List of dates to output data
VarPrintstep = 0 # Locator for performance-testing output files  

#**********************************************************************
#Section VI:  Iterate through timesteps.        

if inputs.Hydrograph == True:
    Reach.Node[0].Load.QsAvkFeed, Reach.Node[0].Load.Qsjkfeed = Reach.Node[0].Load.UpdateElwhaFeedRatingCurve(ControlNode, 
                                                                                                              LoadFactor[n],
                                                                                                              inputs.MudFraction,
                                                                                                              inputs.Qlist[0][counter],
                                                                                                              inputs.vfunc,
                                                                                                              inputs.TrinityFit,
                                                                                                              CalibrationFactor,
                                                                                                              ControlGSD,
                                                                                                              Section = 'Middle')       

print('Model setup complete!  Starting timesteps...')

def update_count_and_save():
    global counter,NextCount,Interval,Printstep,year,OutputObj,Reach
    if counter == int(NextCount):
        NextCount = NextCount + Interval
        year = 0                
        if Tyear == 0:
            year = 0.
        else:
            year = Tyear-(dt/60./60./24./365.25)
        
        print('Saving regular output at count = %s' % counter) 
        for out in inputs.Outputvars:
            OutputObj.Output('Out_'+ out, Reach, out, \
                Printstep, year)

        OutputObj.OutputFlux(Reach, Printstep, year, 'QsOut')
        OutputObj.OutputFlux(Reach, Printstep, year, 'QsIn') 
        OutputObj.OutputFlux(Reach, Printstep, year, 'BankIn') 
        OutputObj.OutputFlux(Reach, Printstep, year, 'BankOut') 
        
        Printstep += 1

def print_validate_data():
    global ValidateDates, subdaycount, inputs, OutputObj, Reach, VarPrintstep
    if date in ValidateDates and subdaycount == 1:
        print('Saving validation output at specified date of: %s ' % date)
        for out in inputs.Validatevars:
            OutputObj.Output('OutValidate_'+ out, Reach, out, \
                VarPrintstep, float(date.year))
        VarPrintstep += 1

def check_hydrograph_and_adjust_dt():
    global inputs, subdaycount, FeedTrigger, Tmult, Tyear, OutputObj, Reach
    global date, LowFlow, dt
    
    if inputs.Hydrograph == True:                
        if subdaycount == Tmult or Tyear == 0: #  Tmult denotes the number of timesteps in each day.
            subdaycount = 0
            FeedTrigger = False
            
            #  Record output variables
            OutputObj.PopulateDailyLists(Reach)
            OutputObj.PopulateDailyDictList(Reach, inputs.DailyOutputVars, date)
            LowFlow = False
            #  Update the discharge
            #print self.inputs.Qlist[counter]
            for Node in Reach.Node:
                for j in range(Reach.NFlows):
                    Node.DC.Qw[j] = inputs.Qlist[Node.Q_ID][counter]
                Node.UpdateDepthAndDischargeAtAllFlows(inputs.vfunc)
                for j in range(Reach.NFlows):
                    Node.DC.WSE[j] = Node.DC.Hc[j] + Node.etabav

            # Find downstream boundary condition for the new discharge.                   
            Reach.find_downstream_boundary(zControlBoundary)   
            
            # Sediment transport only occurs if flow is above a
            # threshold.  Increase timestep if above even higher 
            # threshold. Only works with hydrograph
            #if Reach.Node[0].DC.Qw[0] < 40.:
            if Reach.Node[0].DC.Qw[0] < inputs.LowFlowThreshold:
                LowFlow = True
                dt = dt*Tmult
                Tmult = 1
            #elif Reach.Node[0].DC.Qw[0] > 300. and self.inputs.CyclingHydrograph == False:
            elif Reach.Node[0].DC.Qw[0] > inputs.HighFlowTimestepThreshold and inputs.CyclingHydrograph == False:
                #Tmult = 150
                Tmult = inputs.TmultHighFlow
                dt = 1./365.25/Tmult*365.25*24*60*60
                print ('dt in model set to high flow dt of %s' % dt)
            else:
                Tmult = inputs.Tmult
                if inputs.CyclingHydrograph == True:
                    #dt = 1./365.25/27*365.25*24*60*60
                    dt = 1./365.25/inputs.TmultCyclingHydrograph*365.25*24*60*60
                else:
                    dt = 1./365.25/Tmult*365.25*24*60*60
                    #print ('dt in model = %s' % dt)

    else:
        subdaycount = 1

def update_sediment_supply():
    global Reach, NewAvkFeed, NewJKFeed, inputs, LowFlow, multiplier, Qsjk, Qsav
    
    # Apply feed update after output has been written to file
    Reach.Node[0].Load.QsAvkFeed, Reach.Node[0].Load.Qsjkfeed = NewAvkFeed, NewJKFeed

    # Define any lateral sediment supply here.  In this implementation, the
    # lateral sources are simply a multiple of the sediment supplied
    # to the upper end of the reach, but it should be possible to compute
    # the lateral sources for other geometries using different control
    # nodes. These nodes would need to have transport rates computed
    # in advance of this call, which can be done by 
    # 1) defining the flow duration curve for each control node,
    # which for hydrograph runs is a single discharge and probability,
    # 2) running ControlNode.UpdateDepthAndDischargeAtAllFlows, 
    # 3) calling ControlNode.UpdateSedimentLoadByDurationAndSize, and
    # 4) using ControlNode.Qsjk and ControlNode.QsAvkLoad as the load, 
    # potentially after multiplying by the appropriate multiplier.
    
    # Note that lateral sources are not computed for low flow conditions.
    for i in range(len(inputs.LateralSedimentSourceNodes)):
        index = inputs.LateralSedimentSourceNodes[i]
        if LowFlow == True:
            multiplier = 0
        else:
            multiplier = inputs.LateralMultiplier[i]
        Qsjk = Reach.Node[0].Load.Qsjkfeed
        Qsav = Reach.Node[0].Load.QsAvkFeed
        Reach.Node[index].SLatSourcejk = np.array(Qsjk)*multiplier
        Reach.Node[index].SLatSourceAv = np.array(Qsav)*multiplier        

def sediment_transport_and_exner():
    global LowFlow, Reach, inputs, dt, TracerProperties, CalibrationFactor 
    global ControlGSD
    
    if LowFlow == True:
        Reach.AbridgedStepDownstream(inputs.BcMin, 
                                     inputs.W, 
                                     dt, 
                                     TracerProperties, 
                                     inputs.WidthChange, 
                                     inputs.TransFunc, 
                                     inputs.TrinityFit, 
                                     inputs.CalibrationFactor, 
                                     inputs.alphatau, 
                                     inputs.AlphaBed)
        
    else:
        Reach.StepDownstream(dt, 
                             inputs.AlphaBed, 
                             0.2, 
                             0.1, 
                             inputs.LayerL,
                             TracerProperties,
                             inputs.TransFunc,
                             inputs.TrinityFit,
                             CalibrationFactor,
                             inputs.WidthChange,
                             inputs.W,
                             inputs.ErodeT,
                             inputs.alphatau,
                             inputs.BcMin,
                             inputs.vfunc,
                             inputs.AvulsionThreshold,
                             ControlGSD,
                             inputs.AvulsionExchange,
                             inputs.MobilityThreshold)                                
        Reach.UpdateSlope()                

def partly_alluvial_calcs():
    #This section of code is experimental and is intended to crudely 
    #represent the exposure of bedrock, which would result in a partly 
    #alluvial bed and would halt incision where it is exposed            
    global Node, Reach, inputs
    for Node in Reach.Node:
        if Node.PartlyAlluvial == True:
            if Node.CumulativeBedChange < inputs.PartlyAlluvialMin and not Node.FixedElev:
                Node.FixedElev = True
            if Node.FractionAlluvial > 1.:
                Node.FixedElev = False
                Node.ActiveLayer.Volume = Node.ActiveLayer.L * Node.Bc * \
                    Node.dxc
    #if Reach.Node[-1].CumulativeBedChange > 1.: # Katie add
    #    Reach.Node[-1].etabav = Reach.Node[-1].InitialBedElev + 1.
    #if Reach.Node[-1].CumulativeBedChange < .5: # Katie add
    #    Reach.Node[-1].etabav = Reach.Node[-1].InitialBedElev -.5

def apply_customizations():
    global inputs, date, counter, deltat, timedelta, timeseconds, deltat
    global Reach, BoundaryTrigger, zControlBoundary, DamTyearStart
    global NodeTrigger, startcounter
        
    if inputs.Removal == True:
        #  Prevent the alluvial node between the canyons from widening
#                Reach.Node[1].Canyon = True
#                Reach.Node[1].NarrowRate = 0.
#                Reach.Node[1].WidenRate = 0.

        #  Trigger removal at the proper time and set up time variables
        if date >= datetime.date (2011, 9, 17) or inputs.Hydrograph == False and counter >= 1860:
            deltat = 0  # Time elaspsed since removal                  
            if inputs.Hydrograph == True:                    
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
                if inputs.Hydrograph == True:                        
                    BoundaryTrigger = True                        
                    Reach.SetBoundary = False
                    
                    # Reset hydraulics
                    Reach.set_up_hydrograph()
                    Reach.find_downstream_boundary(zControlBoundary)
                
                    for Node in Reach.Node:
                        for j in range(Reach.NFlows):
                            Node.DC.Qw[j] = inputs.Qlist[Node.Q_ID][counter]
                        Node.UpdateDepthAndDischargeAtAllFlows(inputs.vfunc)
                        for j in range(Reach.NFlows):
                            Node.DC.WSE[j] = Node.DC.Hc[j] + Node.etabav
                    #Reach.find_downstream_boundary(zControlBoundary)
                
                else:
                    BoundaryTrigger = True
                    Reach.SetBoundary = False
                    #Reach.Add_Nodes(zControlBoundary, 5, zControlBoundary.etabav, zControlBoundary.xc, self.inputs.Outputvars, self.inputs.Outputfolder)
                    #zControlBoundary = deepcopy(Reach.Node[-1]) # Set downstream boundary to new downstreammost node                       
                    Reach.find_downstream_boundary(zControlBoundary)
                    
                    for Node in Reach.Node:
                        Node.UpdateDepthAndDischargeAtAllFlows(inputs.vfunc)
                        for j in range(Reach.NFlows):
                            Node.DC.WSE[j] = Node.DC.Hc[j] + Node.etabav
                    

            # Find downstream boundary condition--will want to move this out of the hydrograph loop when simulating the removal with a flow duration curve.
            Reach.find_downstream_boundary(zControlBoundary)

        # This section simulates the upstream sediment feed change
        # for the dam removal--when bedload reaches the 
			 # former Glines site, the feed changes (below).  Based on the date
			 # published in East (2015).                        
        if date == datetime.date(2012, 10, 14) or inputs.Hydrograph == False and counter == 1880:                   
            
            if inputs.Hydrograph == True:
                if NodeTrigger == False:
                    NodeTrigger = True
                    inputs.FeedType = 'DamRemovalFirstPulseRatingCurve'
                    # Reduce the thickness of the active layer--only for 
                    # hydrograph and add remaining active layer as new 
                    # substrate layer

            else:
                inputs.FeedType = 'DamRemovalDuration'
            DamTyearStart = Tyear
            
            startcounter = counter
            
        if date == datetime.date(2013, 10, 1):
            if inputs.Hydrograph == True:                    
                inputs.FeedType = 'DamRemovalFirstPulseRatingCurve'

def change_sediment_feed():
    global LoadGauge, LoadType, date, LoadFactorCount, inputs, FeedTrigger
    global NewAvkFeed, NewJKFeed, Node, ControlNode, LoadFactor, counter
    global MaxSteps, subdaycount, Tmult, Reach, CalibrationFactor, ControlGSD
    global Section, C, tau, t, PercSus, FixedCapacity, feedcounter, startcounter
    global v, n
    LoadGauge = ""
    if LoadType == 'counter':
        LoadGauge = counter
    elif LoadType == 'date':
        LoadGauge = date
    if n < len(LoadFactorCount):
        if inputs.Hydrograph == False and str(LoadGauge) == str(LoadFactorCount[n]):
            n += 1
        elif inputs.Hydrograph == True and str(LoadGauge) == str(LoadFactorCount[n]) and FeedTrigger == False:
            FeedTrigger = True                    
            n += 1
    
    # A method for calculating feed is called here based on the FeedType
    # variable specified by the user on the inputs page.            
    
    # Standard method--designed for the flow duration curve.  
    # Feed is a user-defined fraction of sediment transport capacity                   
    if inputs.FeedType == 'DurationCurve':            
        #Node =  Reach.Node[0]   
        NewAvkFeed, NewJKFeed = Node.Load.UpdateFeedDurationCurve(\
            ControlNode, LoadFactor[n])  
  
    # Rating curve method--the standard for hydrograph runs.
    # Feed is a user-defined fraction of sediment transport capacity
    # for a given flow.  Mud feed is a set proportion of feed for the
    # next finest size class.
    if counter < MaxSteps-1:                                  
        if inputs.FeedType == 'RatingCurve' and subdaycount == Tmult - 1 and inputs.Qlist[0][counter + 1]>40.:                
            Node = Reach.Node[0]
            NewAvkFeed, NewJKFeed = \
                Node.Load.UpdateFeedRatingCurve(ControlNode, LoadFactor[n],\
                inputs.MudFraction, inputs.Qlist[0][counter + 1],\
                inputs.vfunc, inputs.TrinityFit, CalibrationFactor, ControlGSD) 
    
    # ******************Example: ELWHA RATING CURVE *******************
    # This is an example of a customized method of updating the sediment
    # feed.            
    
    # It is the same as 'RatingCurve' method above when calculating feed 
    # for the bed material load, but uses an emperical sediment rating 
    # curve provided by Curran/Konrad for the suspended load.  
    if counter < MaxSteps-1:
        if inputs.FeedType == 'RatingCurveMiddleElwha' and subdaycount == Tmult - 1 and inputs.Qlist[0][counter + 1]>40.:
        #print self.inputs.Qlist                                  
            Node = Reach.Node[0]
            NewAvkFeed, NewJKFeed = \
                Node.Load.UpdateElwhaFeedRatingCurve(ControlNode, LoadFactor[n],\
                inputs.MudFraction, inputs.Qlist[0][counter + 1],\
                inputs.vfunc, inputs.TrinityFit, CalibrationFactor, ControlGSD, Section = 'Middle')
                
    if counter < MaxSteps-1:
        if inputs.FeedType == 'RatingCurveMiddleElwhaStochastic' and subdaycount == Tmult - 1 and inputs.Qlist[0][counter + 1]>40.:
        #print self.inputs.Qlist
            multiplier = np.random.lognormal(0., .9)
            #print multiplier                                  
            Node = Reach.Node[0]
            NewAvkFeed, NewJKFeed = \
                Node.Load.UpdateElwhaFeedRatingCurve(ControlNode, multiplier,\
                inputs.MudFraction, inputs.Qlist[0][counter + 1],\
                inputs.vfunc, inputs.TrinityFit, CalibrationFactor, ControlGSD, Section = 'Middle')

    if counter < MaxSteps-1:      
        if inputs.FeedType == 'RatingCurveUpperElwha' and subdaycount == Tmult - 1:
            Node = Reach.Node[0]
            NewAvkFeed, NewJKFeed = \
                Node.Load.UpdateElwhaFeedRatingCurve(ControlNode, LoadFactor[n],\
                inputs.MudFraction, inputs.Qlist[0][counter + 1],\
                inputs.vfunc, inputs.TrinityFit, CalibrationFactor, ControlGSD, Section = 'Upper')

    # ******************End ELWHA RATING CURVE ************************
    
    # *********Example: ELWHA DAM REMOVAL:  DURATION CURVE ************
    # This is an example of a customized method of updating the sediment
    # feed. 
     
    # Simulates exponential decay of the LoadFactor variable and was 
    # calibrated to data provided by Tim Randle.  It is written to be
    # used with a flow duration curve and activated partway through
    # a run (see Activate Dam Removal, Section VI.D).
    if inputs.FeedType == 'DamRemovalDuration': 
        Node = Reach.Node[0]                
        C = inputs.C1
        tau = inputs.tau
        t = (Tyear-DamTyearStart)*365.25 # In days
        PercSus = inputs.PercSus
        multiplier = ''
            
        if inputs.DecayType == 'Exp':
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
        if inputs.FeedType == 'DamRemovalRatingCurve' and subdaycount == Tmult - 1 and inputs.Qlist[0][counter + 1]>40.:                
#                    DamTyearStart = Tyear
            Node = Reach.Node[0]
            C = inputs.C2
            tau = inputs.tau
            t = (Tyear-DamTyearStart)*365.25
            multiplier = ''
            FixedCapacity = False
            
            if inputs.DecayType == 'Exp':
                multiplier = 1 + C*exp(-t/tau)
            elif inputs.DecayType == 'Pow':
                v = 500 # Hard coded here for now--can make it an input variable later.
                multiplier = 1 + C*(t+1)**(-v/tau)

            if inputs.FixedCapacity == True: # Determines whether feed is based on duration-averaged capacity or capacity for daily discharge
                FixedCapacity = True

            NewAvkFeed, NewJKFeed = Node.Load.UpdateElwhaFeedRatingCurve(ControlNode, 
                                                                         multiplier,
                                                                         inputs.MudFraction,
                                                                         inputs.Qlist[0][counter + 1],
                                                                         inputs.vfunc,
                                                                         inputs.TrinityFit,
                                                                         CalibrationFactor,
                                                                         ControlGSD,
                                                                         Removal = True,
                                                                         PercSus = inputs.PercSus,
                                                                         FixedCapacity = FixedCapacity)

            feedcounter = counter - startcounter
            
    if counter < MaxSteps-1:         
        if inputs.FeedType == 'DamRemovalFirstPulseRatingCurve' and subdaycount == Tmult - 1 and inputs.Qlist[0][counter + 1]>40.:                
#                    DamTyearStart = Tyear
            Node = Reach.Node[0]
            C = inputs.C1
            tau = inputs.tau
            t = (Tyear-DamTyearStart)*365.25 # In days
            PercSus = inputs.PercSus
            multiplier = ''
            FixedCapacity = False
            
            if inputs.DecayType == 'Exp':
                multiplier = 1 + C*exp(-t/tau)
            elif inputs.DecayType == 'Pow':
                v = 500 # Hard coded here for now--can make it an input variable later.
                multiplier = 1 + C*(t+1)**(-v/tau)

            if inputs.FixedCapacity == True: # Determines whether feed is based on duration-averaged capacity or capacity for daily discharge
                FixedCapacity = True
            #print date
            #print multiplier    
            NewAvkFeed, NewJKFeed = Node.Load.UpdateElwhaFeedFineRatingCurve(ControlNode,
                                                                             multiplier,
                                                                             inputs.MudFraction,
                                                                             inputs.Qlist[0][counter + 1],
                                                                             inputs.vfunc,
                                                                             inputs.TrinityFit,
                                                                             CalibrationFactor,
                                                                             ControlGSD,
                                                                             ControlGSD,
                                                                             Removal = True,
                                                                             PercSus = PercSus,
                                                                             FixedCapacity = FixedCapacity,
                                                                             FirstPulse = True)
            #print sum(NewAvkFeed)*60*60*24
            feedcounter = counter - startcounter
        
    # ************End ELWHA DAM REMOVAL:  Rating CURVE ****************

def update_downstream_boundary():
    global BoundaryGauge, BoundaryType, counter, date
    global BoundaryFactorCount, NextBoundary, Reach    
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

def update_time_trackers():
    global date, subdaycount, inputs, counter, dtcount, m, subdaycount 
    global Tmult, MaxSteps, NextCount, Interval, dt, dtlist
    
    if date.day == 1 and subdaycount == 0 and inputs.Hydrograph == True and inputs.CyclingHydrograph == False:
        print(str(date))
    if inputs.Hydrograph == False and counter % 10 == 0:
        print(counter)
    
    #  Change the length of the timestep if at the user-specified point in the model            
    if len(dtcount) > 0 and counter == dtcount[m]:
        dt = dtlist[m+1]* 365.25*24*60*60
        if m < len(dtcount) - 1:
            m = m + 1
    
    subdaycount = subdaycount + 1

    #  Counter only updates for a new day if hydrograph is on.  If using
    #  a duration curve, it updates at each time step.
    if subdaycount == Tmult or inputs.Hydrograph == False:
        counter = counter + 1
        date = date + datetime.timedelta(days=1)
        
    # Filler if you want option to cycle endlessly thorough hydrograph.
    # Not currently connected to anything
    #filler = False 
        
    if inputs.CyclingHydrograph == True and counter == MaxSteps: # Cycles endlessly through the hydrograph.
        counter = 0
        NextCount = Interval

            
    #************************************************
    

while counter < MaxSteps and Reach.Node[3].ActiveLayer.GSD.D50 == Reach.Node[3].ActiveLayer.GSD.D50:
    #************************************************
    #VI.A:  Record output if at the proper timestep
    update_count_and_save()
        
    #  If date is one specified for performance testing, write output.
    print_validate_data()
    
    #************************************************
    #VI.B:  Perform mass conservation (for more details, see the 'Step Downstream'
    #function in the Reach class)
    
    #  If the model is running a hydrograph, the discharge is updated
    #  here at the start of each day and daily output variables are 
    #  recorded in a list.
    check_hydrograph_and_adjust_dt()

    update_sediment_supply()

    
    #  Run through the hydraulics, sediment transport, and mass conservation
    #  computations. See 'Step Downstream' in the Reach class for more
    #  details.
    sediment_transport_and_exner()
    
    #  Update some cumulative variables and the time elapsed            
    Reach.UpdateOutput(dt)
    Tyear = Tyear + dt/ 365.25 / 24. / 60. /60.
    
    #************************************************
    #VI.C:  Check for bedrock and perform partly-alluvial computations if
    #activated
    partly_alluvial_calcs()
    
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
    apply_customizations()
    
    #VI.E:  Update sediment feed.  Custom sediment feed regimes can be called
    #here, but they should be written in the Load class.
    
    # Trigger a change in feed at the proper time--used for the 'default'
    # feed function.
    change_sediment_feed()
    
    #************************************************
    #VI.F:  Update the downstream boundary condition if it is user-defined.
    update_downstream_boundary()
                
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
    update_time_trackers()
    
#**********************************************************************
#VII:  Save final output

#  Save final objects
#        pickle.dump(Reach, open(os.path.join(os.pardir, self.inputs.Outputfolder, "save.Reach"), "wb" )) 
#        pickle.dump(self.inputs, open(os.path.join(os.pardir, self.inputs.Outputfolder, "inputparams.inputs"), "wb"))         

#  Save objects with daily record (for hydrographs)
OutputObj.WriteDailyFiles()
OutputObj.SimpleWriteDailyFile()

# Produce done file
donefile = os.path.join(os.pardir, inputs.Outputfolder, 'Model_finished')
with open (donefile, 'w') as f:
    f.write('Model completed')
    f.close()	
# Starts the model
#run.RunModel()