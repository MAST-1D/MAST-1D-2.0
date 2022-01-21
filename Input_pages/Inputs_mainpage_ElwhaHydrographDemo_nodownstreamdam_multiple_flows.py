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
import csv
sys.path.append("..")
from Hydrology.clsTimeSeries import clsTimeSeries
from MAST_1D.clsModel import clsModel
from MAST_1D.clsInputs import clsInputs

"""
OBJECT AND FUNCTION DEFINITION
User does not modify.
"""                  
def main(inputs):
    
    #  Checks to see if the specified output folder exists and creates it if it doesn't  
    directory = str(os.path.join(os.pardir,inputs.Outputfolder))
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    # Creates a duration curve for the discharge record for a hydrograph run and
    # creates a model object with the proper inputs.
    run = ""
    if inputs.Hydrograph == True:
        newinputs = load_hydrographs(inputs)
        run = clsModel(newinputs)
    else:
        run = clsModel(inputs)
        
    # Starts the model
    run.RunModel()

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

"""
IV.A  TIMESTEP PARAMETERS FOR FLOW DURATION CURVE
"""    
inputs.MaxSteps = 20000 #  Total number of timesteps to run
inputs.dt = [0] # Timestep (years)--can add multiple values for timestep adjustment during run
inputs.dtcount = [] # List of timesteps to instigate timestep interval change (length of list should be one fewer than dt)

"""
IV.B  TIMESTEP PARAMETERS FOR HYDROGRAPH

MAST-1D is set up to run hydrographs with a daily resolution.  It will keep track
of the date and instigate user-specified boundary condition changes based on the date.
"""
inputs.Tmult = 1 # Number of timesteps per day
inputs.StartDate = (2011,9,15) # (year, month, day).

"""
V:  BOUNDARY CONDITIONS
"""
inputs.Removal = False # This is a custom variable that is used to set some boundary
    # condition behavior for the Elwha River Dam removal (see Section VI.D in 
    # clsModel).

#  Upstream sediment feed
inputs.LoadFactor = [1.,1.]  #  List of upstream sediment feed (in fraction of equilibrium capacity)
inputs.LoadFactorCount = [900] #  List of times to instigate feed change--should
    # be an int (# of timesteps) for duration curves or a tuple ((yyyy, m, dd))
    # for hydrographs.  Number of entries should be one less than LoadFactor.
inputs.FeedType = "RatingCurveUpperElwha" # Choose 'DurationCurve' if using a duration 
    # curve or 'RatingCurve' if using a hydrograph.  You can also create a custom
    # method for applying feed--see Section VI.E in clsModel.

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
inputs.DischargeFileID = [0,1] # ID of each of the discharge timeseries.  Each
    # node is given this id, which correspond with the given timeseries files.
inputs.DischargeFileCoords = [0,8000] #downstream coordinate at upper end of a
    # given discharge timeseries
inputs.DischargeFiles = ['Qfile15Sep2011pres','Qfile15Sep2011pres_b'] # Names of discharge file.  It should be
    # stored in the "Discharge_Files" folder in the parent directory.  See examples
    # there for formatting.   
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
    # Number of bins for all flow duration curves
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
#inputs.AvulsionThreshold = 1.25 # Minimum bank height below which avulsion will occur 

inputs.AvulsionThreshold = -1.25 # Minimum bank height below which avulsion will occur 

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
theoretically work.
"""
# Set up as Cosmogenic 14C.
inputs.NTracers = 1 # Number of tracers
inputs.coj = [82.96, 1.98, 15.06] # production from different processes (at surface, presumably--not integrated over depth)
inputs.Lcj = [160., 738., 2688.] # Attenuation rates? (g/cm^2)
inputs.Name = "'14C'"
inputs.DecayConst = 0.000121 # Decay constant
inputs.ProductionRate = 15.1 # Nuclide production rate
inputs.FalloutRate = 0. # Nuclide fallout rate
        
"""
XII. OUTPUT SPECIFICATIONS

There are three types of output:
A. Text files of node attributes for each node over time periods of equal intervals
    (for hydrographs and duration curves)
B. Text files of node attributes for each node on user-specified dates (currently 
    for hydrographs only)
C. JSON files with daily output for 1 node.  The user specifies which nodes to 
    write output for (for hydrographs only)  
"""
# Output folder
RunName = "ElwhaHydrographDemo" # Name of the run (a folder of outputs files will be created
    # under this name).
inputs.Outputfolder = os.path.join("Output","Demo",RunName) # Parent folder for 
    # the output folder 

# A.  General output

# Number of dataslices to save.  Will be written at equal time intervals.
inputs.NumberOfPrintouts = 100 
 
# List of variables to save (as strings of clsNode attributes)
inputs.Outputvars = ['Slope','Bf','Bc',\
'DC.Qwf[-1]','DC.Sf[0]','DC.Hc[0]','DC.Hc[-1]','DC.Uc[-1]','DC.Uc[0]','DC.WSE[0]','DC.WSE[-1]','DC.Qw[0]',\
'DC.Uc[0]','ActiveLayer.GSD.D50',\
'ActiveLayer.GSD.D84','Load.QsavBedTot','Load.QsavBedTotFeed','Floodplain.GSD.D50','etabav',\
'Substrate[-1].C.GSD.D50','CumulativeBedChange','Floodplain.L',\
'Floodplain.GSD.F[0]', 'Floodplain.GSD.D84']

# B. Output on specific dates (for hydrograph runs)

# Dates (yyyy, m, dd) in which to output variables for model validation 
inputs.ValidateDates = [(2011, 9, 30),(2012, 9, 30),(2013, 9, 30),(2014, 9, 30),(2015, 9, 30),(2016, 9, 30)]     
# Variables (attributes of clsNode) in which to output on specific dates for model validation
inputs.Validatevars = ['CumulativeTotFeed','Bc','CumulativeNarrowing','CumulativeWidening','CumulativeTotalAvulsionWidth']

# C. Output daily (for hydrograph runs)

# Nodes that will be output daily.  The outputted variables are hard-coded in
# clsOutputSpecs and more can be added there.  Note that the file sizes of daily
# output can be up to several megabytes.
inputs.DailyNodes = [0,3] 

#**********************************END INPUTS*******************************************

"""
RUN THE MODEL
User does not modify.
"""

if __name__ == '__main__':
    main(inputs)
