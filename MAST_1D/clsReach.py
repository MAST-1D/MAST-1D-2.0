
from clsSubstratePairClass import clsSubstratePairClass
from clsReservoir import clsReservoir
from clsNode import clsNode
import numpy as np
from copy import deepcopy
import os

class clsReach(object):
    """
    A class defining a river reach.
    
    Defines the properties of a river reach within which sediment will be
    routed downstream.  A reach consists of a set of nodes connected in
    series.  Each node represents multiple bends of the channel.  The reach
    class includes methods for computing water level in each node given
    flow in the reach and downstream boundary conditions.  It also includes
    methods for specifying a size-specific sediment feed rate.
    
    Note
    ----
    Because the reach object requires many input parameters upon instantiation,
    all parameters are passed as class attributes of a user-defined "inputs" object.
    Because the nature of the inputs expected by clsReach are quite specific, they
    are documented here rather than in the clsInputs class definition (which is
    left relativly undefined so that the user can add model features
    as needed.  A helpful addition to the code would be to break the 
    clsInput class up so that the variables required by all reaches 
    are explicitly called out, perhaps as a separate sub-class or
    individual well-defined attributes.
    
    Parameters
    ----------
    inputs.Dbdy : array_like(float)
        Sediment grain size at each bin boundary. Referred to as
        BinBdySizes elsewhere.
    inputs.Qw : array_like(float)
        Flows in flow duration distribution applied to reach.
    inputs.NLayers : int
        Number of substrate layers to consider.
    inputs.NTracers : int
        Number of tracers being considered by the model
    inputs.Nnodes : int
        Number of nodes in the reach.
    inputs.ManningStabilizer : ???
        NEEDS DOCUMENTATION
    inputs.p : array_like(float)
        Fraction of time represented by each discharge bin.
    inputs.Fa : array_like(float)
        Initial fractions for each sediment size class in active layer.

        Note that Presently the initial floodplain size fractions are solved for to 
        ensure that exchange between floodplain and channel is initially in
        equilibrium.

    inputs.Bc : float
        Initial channel width (m).
    inputs.Bf : float
        Initial valley width (m).

        Note that it appears that the user enters intial valley width inputs.Bf as a parameter, but 
        the node attribute Node.Bf represents the width of the floodplain, not the
        entire valley, whose width is the sum of Node.Bf + Node.Bc.

    inputs.Cfc : float
        Friction factor for channel. Cfc = 1/Czc^2.
    inputs.Cff : float
        Friction factor for floodplain. Cff = 1/Czf^2.
    inputs.Slope : float
        Initial slope of reach.
    inputs.migration : float
        Reach-average lateral migration rate (m/s) (yes, it's a small number).
    inputs.ChSin : float
        Reach-averaged sinusity (channel length / valley length).
    inputs.dxf : float (optional)
        Floodplain spacing between nodes.
    inputs.FSandSusp : float
        Fraction sand moving in suspension.  This load is in
        in addition to any sand transport computed from the surface-based
        bed material transport equation.
    inputs.FlBed : float
        Floodplain number for suspended bed material.  Determines the
        fraction of bed material flux moving across the floodplain zone
        that deposits on the floodplain as overbank deposition.
    inputs.MudFraction : float
        ???

        It appears that inputs.MudFraction still needs documentation.

    inputs.FloodplainL : float
        Initial floodplain thickness (m).
    inputs.ActiveLayerL : float
        Thickness of active layer (m).
    inputs.Hpb : float
        Point bar thickness (m).
    inputs.lambdap : float
        Porosity of all sediment deposits.
    inputs.Kbar : float
        Parameter influencing fraction mud stored in bars.
        Motivated by Hoey and Ferguson, who created sediment deposits
        in substrate as mixture of load and sediment surface.
    inputs.AlphaBar : float
        Another parameter influencing fraction mud stored in bars.
    inputs.AlphaPartlyAlluvial : float
        Parameter influencing sediment passing through partly alluvial reaches.
    inputs.ncAddons : float
        Total additional form roughness added to grain roughness.
    inputs.ncMultiplier : float
        Multiplier for grain roughness.
    inputs.nf : float
        Floodplain roughness parameter.
    inputs.LayerL : float
        Initial thickness of substrate layers (m).
    inputs.DLag : float
        Relates to lag deposit in nodes. NEED TO DEFINE
    inputs.FLag : float
        Relates to lag deposit in nodes. NEED TO DEFINE

        Note that presently, there is an elevation hard coded into node.

    inputs.vfunc : UNKNOWN
        NEED TO DOCUMENT: PROBABLY A VELOCITY FUNCTION.
    inputs.ReachwideBedrock : bool
        Flag that allows bedrock to be modeled for all nodes.  If False,
        the reach is assumed to always be fully alluvial.
    inputs.SetBoundary : bool
        NEED TO DOCUMENT: PROBABLY HELPS SET DOWNSTREAM BOUNDARY CONDITION.
    inputs.BoundaryFactor : array_like(float)???
        NEED TO DOCUMENT: PROBABLY ADJUSTS BOUNDARY CONDITION
    inputs.TransFunc : UNKOWN
        NEED TO DOCUMENT
    inputs.TrinityFit : Bool
        Flag that determines if Gaeuman fit to Wilcock & Crowe is used.  Regular
        Wilcock and Crowe used if false.
    inputs.CalibrationFactor : Float
        Multiplier applied to all computed bed material transport rates.

    Attributes
    ----------
    Node : array_like(:obj:`MAST_1D.clsNode`, length = NNodes)
        An array of interconnected node objects.
    NFlows : int
        Number of discharge bins considered in flow duration distribution 
    NBedSizes : int 
        Number of bed material sediment size bins
     
    Notes
    -----
    Lots more to document here:
    Qsjkfeed -- [float]
    QsAvkFeed -- [float]
    BoundaryConditionDownstreamEtaBed -- float
    NTracers -- int (number of sizes)
    NFlows -- int (number of flows in FDC)
    NLayers -- int
    CumulativeOutput -- [float] (Cumulative amount of sediment for each size
        class transported from the last node) # Katie add to track reservoir filling!
    Initialized -- bool
    
    """
    Initialized = False
    
    def nnodes(self):
        return len(self.Node)
        
    def QsavBedTotFeed(self):
        return np.sum(self.QsAvkFeed[1:self.NBedSizes + 1])
    
    def __init__(self, inputs): # Katie add Manning Stabilizer
    
        if not self.Initialized:
            self.NBedSizes = len(inputs.Dbdy) - 2
            self.NFlows = len(inputs.Qw)
            self.NTracers = inputs.NTracers
            self.NLayers = inputs.NLayers
            
            self.Node = [clsNode(self.NLayers, self.NTracers, inputs.Dbdy, self.NFlows) \
                for i in range(inputs.Nnodes)]
            
            self.Qsjkfeed = np.zeros((self.NFlows, self.NBedSizes + 1))
            self.QsAvkFeed = np.zeros(self.NBedSizes + 1)
            
            self.Initialized = True
            self.ManningStabilizer = inputs.ManningStabilizer # Katie add
            self.CumulativeOutput = np.zeros(self.NBedSizes+1) # Katie add
            self.CumulativeFeed = np.zeros(self.NBedSizes + 1) # Katie add
            self.CumulativeBankSupply = np.zeros(self.NBedSizes + 1) # Katie add
            self.CumulativeBankSink = np.zeros(self.NBedSizes + 1) # Katie add
            self.BoundaryHc = "" # Katie add
            self.SetBoundary = False
        else:
            raise RuntimeError('Tried to initiate clsReach twice.')
    
    def SetupReach(self, inputs):
        """
        Subroutine that uses the clsInputs object to set up initial 
        properties of the reach, especially the properties of all 
        nodes within the reach.  
        """
        
        i = 0
        for Node in self.Node:
            Node.DC.Qw = inputs.Qw
            Node.DC.p = inputs.p

            Node.ActiveLayer.GSD.F = inputs.Fa 
            Node.ActiveLayer.GSD.UpdateStatistics()
            
            #Node.Floodplain.GSD.F = inputs.Fp
            #Node.Floodplain.GSD.UpdateStatistics()
            Node.Bc = inputs.Bc
            Node.Bf = inputs.Bf-inputs.Bc # Katie change       
            Node.ValleyMargin = inputs.Bf #Katie add
            Node.Cfc = inputs.Cfc # = 1/Cz^2 # Need to check value in Excel version
            Node.Cff = inputs.Cff # Need to check value in Excel version            
            Node.Slope = inputs.Slope                
            migration = inputs.migration # 0.8 m/yr
            Node.cbank = migration / 24. / 60. / 60. / 365.25 # m/s
            Node.ChSin = inputs.ChSin

            #  Only assigns a cross-section length if it is present in the inputs,
            #  or else calculates it from the reach length and number of nodes.
            if hasattr(inputs, 'dxf'):
                Node.dxc = inputs.dxf[i]*Node.ChSin
            else:
                Node.dxc = inputs.reachlength / self.nnodes()
            
            #  Katie add--assigns channel coordinate at beginning (upstream) end of node
            Node.xc = 0
            if i == 0:
                Node.xc = 0
            else:
                Node.xc = self.Node[i-1].xc + self.Node[i-1].dxc
                
            #Node.xc = Node.dxc * (i - 1)
                
            Node.FSandSusp = inputs.FSandSusp
            Node.Flbed = inputs.FlBed
            Node.MudFraction = inputs.MudFraction
            Node.Floodplain.L = inputs.FloodplainL                
            Node.ActiveLayer.L = inputs.ActiveLayerL
            Node.Floodplain.Volume = Node.Floodplain.L * Node.Bf * Node.dxc / \
                Node.ChSin
            Node.ActiveLayer.Volume = Node.ActiveLayer.L * Node.Bc * Node.dxc
            Node.Hpb = inputs.Hpb
            Node.lambdap = inputs.lambdap
            Node.Kbar = inputs.Kbar
            Node.AlphaBar = inputs.AlphaBar
            Node.AlphaPartlyAlluvial = inputs.AlphaPartlyAlluvial
            Node.ncAddons = inputs.ncAddons
            Node.ncMultiplier = inputs.ncMultiplier
            Node.nf = inputs.nf
            
            for m in range(self.NLayers):
                Node.Substrate[m].C.L = inputs.LayerL
                Node.Substrate[m].F.L = inputs.LayerL
                Node.Substrate[m].C.Volume = Node.Substrate[m].C.L * Node.Bc * \
                    Node.dxc
                Node.Substrate[m].F.Volume = Node.Substrate[m].F.L * Node.Bf * \
                    Node.dxc / Node.ChSin
                
                #  Substrate created automatically
                #Node.Substrate[m].C.GSD = deepcopy(Node.ActiveLayer.GSD)
                #Node.Substrate[m].F.GSD = deepcopy(Node.Floodplain.GSD)
                
                #  Katie changed--define substrate (able to handle lag deposits)
                DLag = inputs.DLag
                FLag = inputs.FLag
                totF = 1 + sum(FLag)
                
                Node.Substrate[m].F.GSD = deepcopy(Node.Floodplain.GSD)
                #Node.Substrate[m].F.GSD = deepcopy(Node.ActiveLayer.GSD) # Katie change floodplain substrate to same GSD as active layer
                Node.Substrate[m].C.GSD = deepcopy(Node.ActiveLayer.GSD)
                
                for k in range(1, Node.NSizes + 1):
                    Node.Substrate[m].F.GSD.F[k] = Node.Substrate[m].F.GSD.F[k]/totF # Katie change
                    Node.Substrate[m].C.GSD.F[k] = Node.Substrate[m].C.GSD.F[k]/totF
                    
                for D in DLag:
                    Node.Substrate[m].F.GSD.F[D]=FLag[D]
                    Node.Substrate[m].C.GSD.F[D]=FLag[D]

                Node.Substrate[m].F.GSD.UpdateStatistics()
                Node.Substrate[m].C.GSD.UpdateStatistics()
                
            #  Katie add this--keep track of original substrate attributes.
            #  Can add substrate layers if channel degrades.
            Node.ControlSubstrate = clsSubstratePairClass()
            Node.ControlSubstrate.C = clsReservoir(inputs.Dbdy, inputs.NTracers)
            Node.ControlSubstrate.F = clsReservoir(inputs.Dbdy, inputs.NTracers)
            Node.ControlSubstrate.C.L = inputs.LayerL
            Node.ControlSubstrate.F.L = inputs.LayerL
            Node.ControlSubstrate.C.Volume = Node.ControlSubstrate.C.L * Node.Bc * Node.dxc
            Node.ControlSubstrate.F.Volume = Node.ControlSubstrate.F.L * Node.Bf * \
                    Node.dxc / Node.ChSin
            Node.ControlSubstrate.C.GSD = deepcopy(Node.ActiveLayer.GSD)
            Node.ControlSubstrate.F.GSD = deepcopy(Node.Floodplain.GSD)
            for k in range(1, Node.NSizes + 1):
                Node.ControlSubstrate.F.GSD.F[k] = Node.ControlSubstrate.F.GSD.F[k]/totF # Katie change
                Node.ControlSubstrate.C.GSD.F[k] = Node.ControlSubstrate.C.GSD.F[k]/totF
            Node.ControlSubstrate.F.GSD.UpdateStatistics()
            Node.ControlSubstrate.C.GSD.UpdateStatistics()
             #  Katie move this line from the top of the loop to here.
            Node.etabav = 0.
            if i == 0:
                Node.etabav = 127.
            else:
                Node.etabav = self.Node[i-1].etabav-self.Node[i-1].Slope*self.Node[i-1].dxc
            #Node.etabav = 33.9 - Node.Slope * Node.dxc * i
            Node.InitialBedElev = Node.etabav
            i = i + 1

            #  Makes partly-alluvial function activatable in all nodes if desired
            if hasattr(inputs, 'ReachwideBedrock'):
                if inputs.ReachwideBedrock == True:
                    Node.PartlyAlluvial = True
        # temp workaround for old reach object
        #setattr(self, 'CumulativeBankSink', np.zeros(self.NBedSizes + 1))
        #print self.CumulativeBankSupply
            
        # ***********************Set Up Hydraulics******************************
        # Find Normal Depth--probably can be removed
        for Node in self.Node:
            Node.UpdateDepthAndDischargeAtAllFlows(inputs.vfunc)
            for j in range(self.NFlows):
                Node.DC.WSE[j] = Node.DC.Hc[j] + Node.etabav

        # Find downstream boundary condition
        self.SetBoundary = inputs.SetBoundary
        if self.SetBoundary == True:
            self.BoundaryHc = inputs.BoundaryFactor[0]         
        self.find_downstream_boundary(self.Node[-1])                
        # Compute backwater hydraulics
        self.UpdateManningDepthAtAllFlowAndNodes()

        # *****************************END HYDRAULICS***************************
        
        # *****************************SET FEED*********************************
        for Node in self.Node:
            #filler = Node.ActiveLayer.GSD.D65 # To initialize GSD--or else it changes
            Node.Load.UpdateSedimentLoadByDurationAndSize(Node.DC, \
                Node.ActiveLayer.GSD, 1000., 2.7, Node.Bc, \
                Node.FractionAlluvial, inputs.TransFunc, inputs.TrinityFit, inputs.CalibrationFactor)

        # Assume mud feed is 10 times sand feed--Katie:  now there is a variable for this.
        Node = self.Node[0]
        Node.Load.QsAvkFeed[0] = 0.
        if inputs.FeedType == 'RatingCurveMiddleElwha':
            for j in range(Node.DC.NFlows()):
                Node.Load.Qsjkfeed[j, 0] = Node.Load.CurranLoad(Node.DC.Qw[j])
                Node.Load.QsAvkFeed[0] = Node.Load.QsAvkFeed[0] + \
                    Node.Load.Qsjkfeed[j, 0] * Node.DC.p[j]
        else:        
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
        with open(os.path.join(os.pardir, inputs.Outputfolder, "TauRC"), 'w') as f:
            f.write('{:10}{:10}'.format('Q', 'Tau (rho*g*Hc*Sf)') + '\n')
            for k in range(Node.DC.NFlows()):
                shear = 1000*9.81*Node.DC.Hc[k]*Node.DC.Sf[k]
                f.write('{:<10.4}{:<10.4}'.format(Node.DC.Qw[k],shear) + '\n')
            f.close()
   
        # ********************Output gravel rating curve for feed*************
#        with open('Output' + '//' + 'FeedRC', 'w') as f:  
        with open(os.path.join(os.pardir, inputs.Outputfolder, "FeedRC"), 'w') as f:
            f.write('{:10}{:10}{}'.format('', '', 'Particle Diameter (mm)') + '\n'\
                + '{:10}{:10}'.format('', ''))
            for k in range(0, Node.NSizes + 1):
                f.write('{:<10.4}'.format(Node.Load.GSDBedloadAv.D[k]))
            f.write('\n\n' + '{:10}'.format(''))
    
            f.write('{:10}{}'.format('Q (cms)', 'Load (m3/yr)') + '\n')
            for j in range(self.NFlows):
                f.write('{:10}{:<10.4}'.format('', Node.DC.Qw[j]))
                for k in range(0, Node.NSizes + 1):
                    f.write('{:<10.4}'.format(Node.Load.Qsjkfeed[j, k] * 60. * 60.\
                        * 24. * 365.25))
                f.write('\n')
            
            f.write('\n')
            f.write('{:10}{:10}{}'.format('probabilit', 'Q(cms)', 'Load * p[j]' + \
                'm3/yr') + '\n')
            for j in range(self.NFlows):
                f.write('{:<10.4}{:<10.4}'.format(Node.DC.p[j], Node.DC.Qw[j]))
                for k in range(0, Node.NSizes + 1):
                    f.write('{:<10.4}'.format(Node.Load.Qsjkfeed[j, k] * \
                        Node.DC.p[j] * 60. * 60. * 24. * 365.25))
                f.write('\n')
            f.write('{:10}{:10}'.format('', 'Sum'))
            for k in range(0, Node.NSizes + 1):
                f.write('{:<10.4}'.format(Node.Load.QsAvkFeed[k] * 60. * 60. * 24.\
                    * 365.25))
    
        # ******************End Output Rating Curve for Feed***********************
        
        # *****************************END FEED************************************
        
        # *************************MAKE SURE SYSTEM IS AT GEOMORPHIC EQUILIBRIUM***
        # Compute and set the floodplain number required to achieve the correct 
        # bank height at equilibrium conditions
        Node = self.Node[0]
        Fl = Node.EquilibriumMudFloodplainNumber(Node.Floodplain.L - Node.Hpb - \
            Node.ActiveLayer.L, Node.cbank)

        for Node in self.Node:
            Node.Flmud = Fl


        # Determine migration factor based on mud deposition rate at upstream node.
        # This ensures that true equilibrium is possible.
        # Otherwise, migration at low flow exports most of the mud from the 
        # floodplain.
        
        # Katie note:  this is not accounted for in the hydrograph runs--does not matter if width change is on.
        Migrationfactor = [0.] * self.NFlows
        Node = self.Node[0]
        if Node.cbank > 0.:
            Node.UpdateDMudj()
            for j in range(self.NFlows):
                Migrationfactor[j] = Node.Dfjk[j, 0] * Node.DC.p[j] / Node.Dfav[0]
            for i in range(self.nnodes()):
                for j in range(self.NFlows):
                    self.Node[i].DC.Migrationfactor[j] = Migrationfactor[j]
       # print Node.DC.Migrationfactor
    
        # Compute the size distribution of the floodplain at equilibrium
        FFpEquilib = [0.] * (self.NBedSizes + 1)
        FFpEquilib[0] = (Node.ActiveLayer.GSD.F[0] * Node.ActiveLayer.L + \
            Node.FkPointBarDepositAnnualAverage[0] * Node.Hpb + \
            (Node.Floodplain.L - Node.ActiveLayer.L - Node.Hpb)) / \
            (Node.Floodplain.L)

        for k in range(1, self.NBedSizes + 1):
            FFpEquilib[k] = (Node.ActiveLayer.GSD.F[k] * Node.ActiveLayer.L + \
                Node.FkPointBarDepositAnnualAverage[k] * Node.Hpb) / \
                Node.Floodplain.L
        for i in range(self.nnodes()):
            for k in range(self.NBedSizes + 1):
                self.Node[i].Floodplain.GSD.F[k] = FFpEquilib[k]
        
        for i in range(self.nnodes()):
            for m in range(self.NLayers):
                self.Node[i].Substrate[m].F.GSD = \
                    deepcopy(self.Node[i].Floodplain.GSD) # Katie change from Floodplain to ActiveLayer--then changed back.
        Node.ControlSubstrate.F.GSD = deepcopy(Node.Floodplain.GSD)
        # Katie add--to try to conserve mud
        #for Node in self.Node:
         #   Node.ActiveLayer.GSD.F[0] = 0.
          #  Node.ActiveLayer.GSD.UpdateStatistics()
        # *********************************END GEOMORPHIC EQUILIBRIUM**************
                    
        
        # ***********************************SETUP TRACERS*************************
        Node = self.Node[0]
        for k in range(Node.NSizes + 1):
            if k == 0:
                # Set tracer concentration in mud feed to 1
                for j in range(self.NFlows):
                    Node.Load.TMudFeedj[j, 0] = 1.
            else:
                # Set tracer concentration in bed material feed to 1
                Node.Load.TBedFeedk[k, 0] = 1.
            
            # Set tracer concentration in all reservoirs to arbitrary value
            ArbitraryValue = 1.
            for i in range(self.nnodes()):
                self.Node[i].Floodplain.T[k, 0] = ArbitraryValue
                self.Node[i].ActiveLayer.T[k, 0] = ArbitraryValue
                for m in range(self.NLayers):
                    self.Node[i].Substrate[m].C.T[k, 0] = ArbitraryValue
                    self.Node[i].Substrate[m].F.T[k, 0] = ArbitraryValue
        # ********************************END SETUP TRACERS ***********************

    def SetupTracers(self, inputs, TracerProperties):
        """
        Sets up tracers using parameters appropriate for cosmogenic 14C
        
        Parameters
        ----------
        TracerProperties : :obj:`MAST_1D.clsTracerProperties`
            Tracer properties object to modify.
        inputs.coj : array_like(float, length = 3)
            Surface production rates from spallation and fast or slow muons
        inputs.Lcj :  array_like(float, length = 3)
            Attenuation lengths from spallation, fast, or slow muons.
        inputs.name : string
            Name of tracer under consideration.
        inputs.DecayConst : float
            Decay constant for tracer (1/yr)
        inputs.ProductionRate : float   
            Surface production rate of tracer, Atoms/g SiO2/yr.
        inputs.FalloutRate : float
            For fallout radionuclides, fallout rate in atoms/cm2/yr.
                
        Returns
        -------
        :obj:`MAST_1D.clsTracerProperties`
            Object formatted with specified tracer properties.
        """
        TracerProperties.coj[0] = inputs.coj[0] # % production from this process 
            # (at surface, presumably--not integrated over depth)
        TracerProperties.coj[1] = inputs.coj[1] # % production from this process 
            # (at surface, presumably--not integrated over depth)
        TracerProperties.coj[2] = inputs.coj[2] # % production from this process 
            # (at surface, presumably--not integrated over depth)
        TracerProperties.Lcj[0] = inputs.Lcj[0] # g/cm2
        TracerProperties.Lcj[1] = inputs.Lcj[1] # g/cm2
        TracerProperties.Lcj[2] = inputs.Lcj[2] # g/cm2
        TracerProperties.Name = inputs.Name
        TracerProperties.DecayConst = inputs.DecayConst # 1/yr
        TracerProperties.ProductionRate = inputs.ProductionRate # Atoms/g Si02/yr
        TracerProperties.FalloutRate = inputs.FalloutRate # Atoms/cm2/yr
        
        return TracerProperties

    def find_downstream_boundary(self, ControlNode):
        """
        Finds the downstream boundary condition given a Node object and list and
        number of discharges.
        
        Katie moved here.
        """

        # Find downstream boundary condition--Katie commented out for normal flow
        for j in range(self.NFlows):
            # this is a poor way to set the downstream boundary condition, but it
            # cannot be solved directly using the normal depth solution implemented
            # in the node class since the definition of friction slope is slightly
            # different.
            # This simply iterates through depth until water surface slope is 
            # parallel to the bed.
            
            depth = 0. # Katie add this and following if statement to allow user to set water surface elevation
            if self.SetBoundary == True:
                depth = self.BoundaryHc
                
            else:
                depth = self.ManningNormalDepthForBoundary(ControlNode.Slope, ControlNode.Bc, \
                    ControlNode.Bf, ControlNode.nc(), ControlNode.nf, ControlNode.Floodplain.L - \
                    ControlNode.ActiveLayer.L, self.Node[-1].DC.Qw[j], ControlNode.ChSin)
            self.Node[-1].DC.WSE[j] = depth + ControlNode.etabav
            #print ControlNode.etabav
   
    def set_up_hydrograph(self):
        """
        After the floodplain initial conditions are set up, the flow and transport
        bins are reset to accomodate just one flow. Katie add.
        """
        
        self.NFlows = 1
        for Node in self.Node:
            #Node.Flmud = .05 # Here is where Katie hard-codes the floodplain number.
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
            Node.DC.Migrationfactor = np.zeros(Node.NFlows)
            Node.DC.Migrationfactor[0] = 1.
            Node.Dfjk = np.zeros((Node.NFlows, Node.NSizes + 1))
            Node.Dfav = np.zeros((Node.NSizes + 1))
            #InitialMudFeed = sum(Node.Load.Qsjkfeed[:,0])
            #InitialMudLoad = sum(Node.Load.Qsjk[:,0])
            Node.Load.SinkSed = np.zeros((Node.NFlows, Node.NSizes + 1))
            Node.Load.LatSed = np.zeros((Node.NFlows, Node.NSizes + 1))
            Node.Load.Qsjkfeed = np.array([deepcopy(Node.Load.QsAvkFeed)])
            Node.Load.QsAvkLoad = deepcopy(Node.Load.QsAvkFeed)            
            #Node.Load.Qsjkfeed[0,0]= InitialMudFeed
            Node.Load.Qsjk = np.array([deepcopy(Node.Load.QsAvkFeed)])#Feed because we assume none gets left on fp initially. Shouldn't matter after first timestep
            #Node.Load.Qsjk[0,0] = InitialMudLoad
            Node.Load.QsjTot = np.zeros(Node.NFlows)
            #Node.DC.WSE[0] = Node.etabav
            Node.DC.Sf = np.zeros(Node.NFlows) 
            
    
    def BackwaterPredictor(self, WSE1, KC1, KF1, KC2, KF2, AC1, AF1, AC2, \
        AF2, dxc, dxf, Q1, Q2):
        """
        Applies energy equation between simplified cross sections.
        
        This function computes an upstream predicted water surface elevation 
        using the energy equation, assuming subcritical flow and ignoring
        expansion/contraction losses (which should be accounted for in the estimate 
        of the friction coefficients).  It is written using the conveyance form
        of the energy equation, which means that it does not require any 
        specific parameterization for friction as long as conveyance can be 
        estimated. The algorithm is based on the USACE HEC-RAS hydraulic 
        reference manual.
    
        Parameters
        ----------
        WSE1 : float 
            Water surface elevation in downstream section.
        KC1 : float 
            Channel conveyance in section 1. For manning,
            k = 1/n*A*R^(2/3))
        KF1 : float
            Floodplain conveyance in secion 1.
        KC2 : float
            Channel conveyance in section 2.
        KF2 : float 
            Floodplain conveyance in section 2.
        AC1 : float
            Flow area for channel in section 1.
        AF1 : float 
            Flow area for floodplain in section 1.
        AC2 : float 
            Flow area for channel in section 2.
        AF2 : float 
            Flow area for floodplain in section 2.
        dxc : float 
            Down channel length between sections.
        dxf : float 
            Down foodplain length between sections.
        Q1 : float 
            Total discharge in m3/s at section 2.
        Q2 : float 
            Total discharge in m3/s at section 2.
        WSE2Guess : float 
            Estimated water surface elevation in upstream section.       
        
        Returns
        -------
        float
            Water surface elevation in upstream section.
        """
        
        if AF1 > 0.:
            a1 = (AC1 + AF1) ** 2 / (KC1 + KF1) ** 3 * (KF1 ** 3 / AF1 ** 2 + \
                KC1 ** 3 / AC1 ** 2)
        else:
            a1 = 1.
        
        if AF2 > 0.:
            a2 = (AC2 + AF2) ** 2 / (KC2 + KF2) ** 3 * (KF2 ** 3 / AF2 ** 2 + \
                KC2 ** 3 / AC2 ** 2)
        else:
            a2 = 1.
            
        Vbar1 = Q1 / (AC1 + AF1)
        Vbar2 = Q2 / (AC2 + AF2)
        
        Sf = ((Q1 + Q2) / (KC1 + KF1 + KC2 + KF2)) ** 2
        
        QC1 = KC1 * Sf ** 0.5
        QF1 = KF1 * Sf ** 0.5
        QC2 = KC2 * Sf ** 0.5
        QF2 = KF2 * Sf ** 0.5
        
        dxbar = ((QC1 + QC2) * dxc + (QF1 + QF2) * dxf) / (QC1 + QF1 + QC2 + \
            QF2)
        
        WSE2 = WSE1 + a1 * Vbar1 ** 2 / 2. / 9.807 - a2 * Vbar2 ** 2 / 2. / \
            9.807 + Sf * dxbar
        return WSE2

    def ManningBackwaterDepthDischargeandSf(self, Q1, Q2, DSNode1, USNode2, \
        WSE1, WSE2guess,j,ManningStabilizer): 
        """
        Applies backwater equation between two nodes.

        Given the geometry of an upstream and downstram node, discharge, and
        water level at th downstream node, computes water level, channel discharge,
        and friction slope for the upstream node.

        Parameters
        ----------
        Q1 : float
            Discharge in the downstream section.
        Q2 : float
            Discharge in the upstream section.
        DSNode1 : :obj:`MAST_1D.clsNode`
            The downstream node.
        USNode2 : :obj:`MAST_1D.clsNode`
            The upstream node.
        WSE1 : float
            The water surface elevation in the downstream node.
        WSE2guess : float
            An initial estimate for the water surface elevation at the upstream node.
        j : int
            Number of discharge bins.
        ManningStabilizer : ???
            NEEDS TO BE DOCUMENTED.

        Returns
        -------
        HC2 : float 
            Channel depth at the upstream section.
        Q2 : float  
            Channel discharge in the upstream section.
        Sf : float
            Friction slope.
        """
        # Katie add ManningStabilizer to allow user to adjust amount of chanel in iteration for stabilization purposes
        
        # Compute flow area and conveyance for downstream section
        HC1 = WSE1 - DSNode1.etabav
        HF1 = WSE1 - (DSNode1.etabav + DSNode1.Floodplain.L - \
            DSNode1.ActiveLayer.L)
            
        ## Katie add for stabilization
        #if DSNode1.Floodplain.L <= DSNode1.ActiveLayer.L:
        #    HF1 = WSE1 - 1.

        HF1 = float(HF1)
        if HF1 < 0.:
            HF1 = 0.
        AC1 = HC1 * DSNode1.Bc
        AF1 = HF1 * DSNode1.Bf
        # Conveyance computed using Manning's eq. since we have Manning's 
        # roughness estimates on the Ain. Can be changed to be consistent with
        # Chezy eq.
        KC1 = 1 / DSNode1.nc() * AC1 * HC1 ** (2. / 3.)

        KF1 = 1 / DSNode1.nf * AF1 * HF1 ** (2. / 3.)
    
        dxc = USNode2.dxc
        dxf = USNode2.dxc * USNode2.ChSin

        RefinedWSEGuess = WSE2guess

        error = 1.
        Count = 0

        while abs(error) >= .001 and Count < 19: # Katie change error from 0.001
            # Update conveyance at section 2 based on WSE2Guess
            HC2 = RefinedWSEGuess - USNode2.etabav
            HF2 = RefinedWSEGuess - (USNode2.etabav + USNode2.Floodplain.L - \
                USNode2.ActiveLayer.L)
            
            ## Katie add for stabilization
            #if USNode2.Floodplain.L <= USNode2.ActiveLayer.L:
            #    HF2 = RefinedWSEGuess - 1.
                
            if HF2 < 0.:
                HF2 = 0.
            AC2 = HC2 * USNode2.Bc
            AF2 = HF2 * USNode2.Bf
            # Conveyance computed using Manning's eq. since we have Manning's
            # roughness estimates on the Ain. Can be changed to be consistent 
            # with Chezy eq.
            KC2 = 1. / USNode2.nc() * AC2 * HC2 ** (2. / 3.)
            KF2 = 1. / USNode2.nf * AF2 * HF2 ** (2. / 3.)
            
            WSEpredicted = self.BackwaterPredictor(WSE1, KC1, KF1, KC2, KF2, \
                AC1, AF1, AC2, AF2, dxc, dxf, Q1, Q2)
            error = RefinedWSEGuess - WSEpredicted
            if WSEpredicted-USNode2.etabav>3: #Katie add if statement
                RefinedWSEGuess = (RefinedWSEGuess) - (error / ManningStabilizer) # Katie change from error/3.
            else:
                RefinedWSEGuess = (RefinedWSEGuess) - (error / 1000000.)
            Count = Count + 1

        if Count == 20:
            print('Backwater iteration did not converge in 20 iterations')
        Sf = (Q2 / (KC2 + KF2)) ** 2
        Qc = KC2 / (KC2 + KF2) * Q2
        
        return np.array([HC2, Qc, Sf])
    
    def ManningNormalDepthForBoundary(self, Slope, Bc, Bf, nc, nf, \
        ChannelDepth, Q, Sinuosity):
        """
        Solves for steady uniform flow depth.
        
        This function iterates for normal depth (steady uniform flow) using 
        the backwater solver. Normal depth is assumed to have been found 
        when the calculated water surface is parallel to the specified bed 
        slope for an arbitrarily specified channel distance (100 m).
        
        Parameters
        ----------
        Slope : float
            Average channel slope.
        Bc : float
            Channel width (m).
        Bf : float
            Floodplain width (m).
        nc : float
            Manning's n for channel.
        nf : float
            Manning's n for floodplain.
        ChannelDepth : float
            Depth of flow in channel.
        Q : float
            Total discharge for cross-section.
        Sinuosity : float
            Channel sinuosity.
            
        Returns
        -------
        float
            Steady uniform flow depth in channel. 
        """
        FlowDepth = (nc * Q / Bc / Slope ** 0.5) ** (3. / 5.)
        error = 1.
        Count = 0
        if FlowDepth > ChannelDepth:
            MaxDepth = FlowDepth
            MinDepth = ChannelDepth
            while abs(error) >= 1e-05 and Count < 50:
                Ac = FlowDepth * Bc
                Kc = 1. / nc * Ac * FlowDepth ** (2. / 3.)
                Af = (FlowDepth - ChannelDepth) * Bf
                Kf = 1. / nf * Af * (FlowDepth - ChannelDepth) ** (2. / 3.)
                
                WsePred = self.BackwaterPredictor(FlowDepth, Kc, Kf, Kc, Kf, \
                    Ac, Af, Ac, Af, 100., 100. * Sinuosity, Q, Q)
                DepthPred = WsePred - Slope * 100.
                error = FlowDepth - DepthPred
                
                if DepthPred > FlowDepth:
                    MinDepth = FlowDepth
                else:
                    MaxDepth = FlowDepth
                
                FlowDepth = (MaxDepth + MinDepth) / 2.
                Count = Count + 1
            if Count > 50:
                print('Normal depth computation did not converge at \
                    downstream boundary')
        return FlowDepth
    
    def UpdateManningDepthAtAllFlowAndNodes(self):
        """
        Iterates through all nodes to perform hydraulic computations.
        """
        for i in range(self.nnodes() - 2, -1, -1): # Loop from downstream to
        # upstream (high i to low) For gradually-varied flow, should be self.nnodes()-2
            for j in range(self.NFlows):
                result = self.ManningBackwaterDepthDischargeandSf(self.\
                    Node[i + 1].DC.Qw[j], self.Node[i].DC.Qw[j], \
                    self.Node[i + 1], self.Node[i], \
                    self.Node[i + 1].DC.WSE[j], self.Node[i].DC.WSE[j],self.Node[-1].DC.WSE[j]-self.Node[-1].etabav,self.ManningStabilizer) # Katie add Manning Stabilizer # Katie change to node index from plus to minus
            
#                result = self.Node[i].NormalChannelDepthAndDischarge(self.Node[i].DC.Qw[j], .1, True)
                #print self.Node[i].DC.Qw[j]                
                #print result[0]
                self.Node[i].DC.WSE[j] = result[0] + self.Node[i].etabav
                self.Node[i].DC.Qwc[j] = result[1]
                
                self.Node[i].DC.Sf[j] = result[2] # Gradually-varied flow                
#                self.Node[i].DC.Sf[j] = self.Node[i].Slope # Uniform flow
                
                self.Node[i].DC.Qwf[j] = self.Node[i].DC.Qw[j] - \
                    self.Node[i].DC.Qwc[j]

                self.Node[i].DC.Hc[j] = self.Node[i].DC.WSE[j] - \
                    self.Node[i].etabav
                        
                self.Node[i].DC.Uc[j] = self.Node[i].DC.Qwc[j] / \
                    (self.Node[i].Bc * self.Node[i].DC.Hc[j])

#                # Katie add--to prevent crashing
#                if self.Node[i].DC.Hc[j] != self.Node[i].DC.Hc[j] or self.Node[i].DC.Hc[j]<0.:
#                    #print self.Node[i].nc()
#                    #print i
#                    Node = self.Node[i]
#                    output = Node.NormalChannelDepthAndDischarge(Node.DC.Qw[j], .01, True)
#                    Node.DC.Hc[j] = output[0]
#                    #print Node.DC.Hc
#                    Node.DC.Qwc[j] = output[0]
#                    Node.DC.Qwf[j] = 0.#Node.DC.Qw[j] - Node.DC.Qwc[j]
#                    Node.DC.Sf[j] = Node.Slope
#                    Node.DC.WSE[j] = Node.etabav + Node.DC.Hc[j]
#                    
#                                    
#                    self.Node[i].DC.Uc[j] = 0.
                    #print Node.DC.Hc

                if self.Node[i].DC.Hc[j] > self.Node[i].Floodplain.L - \
                    self.Node[i].ActiveLayer.L:
                    self.Node[i].DC.Hf[j] = self.Node[i].DC.Hc[j] - \
                        (self.Node[i].Floodplain.L - \
                        self.Node[i].ActiveLayer.L)
                    self.Node[i].DC.Uf[j] = self.Node[i].DC.Qwf[j] / \
                        (self.Node[i].DC.Hf[j] * self.Node[i].Bf)
                else:
                    self.Node[i].DC.Hf[j] = 0.
                    self.Node[i].DC.Uf[j] = 0.

        # Set parameters for downstream-most node
        for j in range(self.NFlows):
            self.Node[-1].DC.Hc[j] = self.Node[-1].DC.WSE[j] - \
                self.Node[-1].etabav
            self.Node[-1].DC.Hf[j] = self.Node[-1].DC.WSE[j] - \
                (self.Node[-1].etabav + self.Node[-1].Floodplain.L - \
                self.Node[-1].ActiveLayer.L)
            if self.Node[-1].DC.Hf[j] < 0.:
                self.Node[-1].DC.Hf[j] = 0.        
            Ac = self.Node[-1].DC.Hc[j] * self.Node[-1].Bc
            Af = self.Node[-1].DC.Hf[j] * self.Node[-1].Bf
            
            Kc = 1. / self.Node[-1].nc() * Ac * (self.Node[-1].DC.Hc[j]) ** \
                (2. / 3.)
            Kf = 1. / self.Node[-1].nf * Af * (self.Node[-1].DC.Hf[j]) ** \
                (2. / 3.)

            self.Node[-1].DC.Qwc[j] = self.Node[-1].DC.Qw[j] * Kc / (Kc + Kf)
            self.Node[-1].DC.Qwf[j] = self.Node[-1].DC.Qw[j] * Kf / (Kc + Kf)
            
            self.Node[-1].DC.Uc[j] = self.Node[-1].DC.Qwc[j] / Ac
            if Af > 0.:
                self.Node[-1].DC.Uf[j] = self.Node[-1].DC.Qwf[j] / Af
            else:
                self.Node[-1].DC.Uf[j] = 0.
            self.Node[-1].DC.Sf[j] = (self.Node[-1].DC.Qw[j] / (Kc + Kf)) ** 2
    
    def AbridgedStepDownstream(self, BcMin, W, dt, TracerProperties, WidthChange, TransFunc, TrinityFit, CalibrationFactor, alphatau, alphabed):
        """
        Katie add.  This function is used during low flow, when only narrowing
        is allowed to occur.  It makes all exchange functions zero except narrowing.
        
        Right now, lateral inputs are not 'turned off'--should fix this.
        """
        self.UpdateManningDepthAtAllFlowAndNodes()
        for i in range(self.nnodes()):
            Node = self.Node[i] 
            Node.Load.UpdateSedimentLoadByDurationAndSize(Node.DC, \
                Node.ActiveLayer.GSD, 1000., 2.7, Node.Bc, \
                Node.FractionAlluvial, TransFunc, TrinityFit, CalibrationFactor) # Only for determining point bar GSDs.
            if Node.Canyon == False and WidthChange == True:
                Node.Narrowing(BcMin, W, alphatau)
                #print Node.ActiveLayer.Volume
            Node.UpdateLateralSedFluxes()
            Node.UpdateLateralTracerConcentrations() # Note that lateral tracer 
                # concentrations in mud are set using the mud feed
            Node.WidenRate = 0.
            Node.cbank = 0.
            Node.DeltaEtaB = 0.
            
            for k in range(Node.NSizes + 1):
                Node.Load.QsAvkLoad[k] = 0.
                Node.Load.QsAvkFeed[k] = 0.
                Node.Dfav[k] = 0.
                for j in range(Node.NFlows):
                    Node.Load.Qsjkfeed[j,k] = 0.
                    Node.Load.Qsjk[j,k] = 0.
                    Node.Dfjk[j,k] = 0.
                    

            Node.UpdateNetSedimentSourcesAndSinks() # Overbank deposition is 
                # included here as floodplain source and load or activelayer 
                # sink (depending on size)
            Node.UpdateVerticalExchangeSedFluxes(Node.DeltaEtaB, alphabed)
            Node.UpdateVerticalExchangeTracerConcentrations()
            
            Node.UpdateExtraWidthFluxes(Node.WidenRate, Node.NarrowRate, Node.DeltaEtaB) # Katie add
            Node.UpdateVolumeSizeFractionsandTracersInAllReservoirs(dt, \
                TracerProperties)
            Node.UpdateGeometricParameters(dt)
            
    def StepDownstream(self, dt, alphabed, LMinAfterSplit, LMinBeforeRemove, \
        SubstrateSpacing, TracerProperties, TransFunc, TrinityFit, CalibrationFactor,\
        WidthChange, W, ErodeT, alphatau, BcMin, Manning, AvulsionThreshold, ControlGSD,\
        AvulsionExchange, MobilityThreshold):
        """
        Subroutine that performs the computations necessary to step forward in time.
        
        Parameters
        ----------
        dt : float
            Timestep.
        alphabed : float
            Hoey-Ferguson parameter determining the ratios of load to active later that
            get transferred to channel substrate when aggrading.
        LMinAfterSplit : float
            Minimum thickness (m) for a substrate node if a new substrate layer is to be spawned (when aggrading).
            This It is necessary to split substrate so as not to mix material too deeply.
        LMinBeforeRemove : float
            Minimum thickness (m) allowed for the top substrate layer if system is degrading.  
            If too thin, the layer is removed and its sediment is mixed with the next lower layer.             
        SubstrateSpacing : float
            Thickness of all but top-most substrate layers (m).
        TracerProperties : :obj:`MAST_1D.clsTracerProperties`
            Properties of tracer under consideration.
        TrinityFit : bool
            Flag that determines if Gaeuman fit to Wilcock & Crowe is used.  Regular
            Wilcock and Crowe used if false.
        CalibrationFactor : float
            Multiplier applied to all computed bed material transport rates.
        W : float (optional)
            Parameter for width change function.
        ErodeT : float (optional)
            Parameter for width change function.
        Bp : float (optional)
            Parameter for width change function.
        Manning : bool (optional)
            Parameter for width change function.
        """
        Dmudj = np.zeros(self.NFlows)
        
        # do hydraulics, sediment transport, morphodynamics and sediment mass
        # conservation
        
        #********
        #Calculate hydraulics
        
        self.UpdateManningDepthAtAllFlowAndNodes()
        for i in range(self.nnodes()):
            if i == self.nnodes()-1:
                x = 'filler'
            # Call sediment transport calculation.  This requires the 
            # specification of the equation to use (WC 2003 or Guiman 2009) and
            # any calibration necessary for reference shear stress.
            Node = self.Node[i]                          
            #print i
            
            #********
            #Calculate sediment transport            
            
            Node.Load.UpdateSedimentLoadByDurationAndSize(Node.DC, \
                Node.ActiveLayer.GSD, 1000., 2.7, Node.Bc, \
                Node.FractionAlluvial, TransFunc, TrinityFit, CalibrationFactor)
            
            #********
            #Determine width change
            
            if  Node.Canyon == False and WidthChange == True:
                Node.WidthChange(Manning, TrinityFit, CalibrationFactor, 2.7, 1000.,\
                    9.81, W, ErodeT, alphatau, BcMin, dt, ControlGSD, MobilityThreshold) # Katie add
            
            #********
            #Calculate lateral reservoir exchanges
            
            Node.UpdateLateralSedFluxes()
            Node.UpdateLateralTracerConcentrations() # Note that lateral tracer 
                # concentrations in mud are set using the mud feed
            
            Node.UpdateDMudj() # Note that overbank deposition is 
                # considered a source/sink term, not a lateral flux term, and 
                # is computed using mud feed

            Node.UpdateDSandj()
            for j in range(self.NFlows):
                Dmudj[j] = Node.Dfjk[j, 0]
            Node.Load.UpdateMudTracerAverages(Dmudj, Node.DC) # the deposition
                # averaged mud tracer concentration is needed for computing mud 
                # deposition concentration for floodplain
            
            Node.UpdateNetSedimentSourcesAndSinks() # Overbank deposition is 
                # included here as floodplain source and load or activelayer 
                # sink (depending on size)
            Node.UpdateNetTracerSourceAndSinkConcentrations() # Tracer 
                # concentration in overbank deposition is computed from 
                # deposition-averaged mud feed
            
            #********
            #Determine vertical bed change and reservoir exchanges
            
            Node.ExnerBed()
            #print 'reach' + str(Node.DeltaEtaB)          
            Node.UpdateVerticalExchangeSedFluxes(Node.DeltaEtaB, alphabed)
            Node.UpdateVerticalExchangeTracerConcentrations()
            
            Node.UpdateExtraWidthFluxes(Node.WidenRate, Node.NarrowRate, Node.DeltaEtaB) # Katie add
            
            Node.Load.ApplyMassConservationToMudLoad(Node.DC) # This cannot be 
                # applied until here since vertical exchange fluxes influence 
                # mud conservation
            Node.Load.ApplyTracerConservationToMudLoad(Node.DC)
            Node.Load.UpdateMudTracerAverages(Dmudj, Node.DC) # This is 
                # probably not required, but is included here so that the mud 
                # tracer properties can be printed out if desired
            Node.UpdateVolumeSizeFractionsandTracersInAllReservoirs(dt, \
                TracerProperties)
            
            #********
            #Determine if the channel avulses
            
            if Node.Canyon == False:
                Node.Avulsion(AvulsionThreshold, SubstrateSpacing, AvulsionExchange) # Katie add
            Node.UpdateGeometricParameters(dt)
            Node.SplitOrCombineSubstrate(LMinAfterSplit, LMinBeforeRemove, \
                SubstrateSpacing)

            # Still need to implement update of tracer production and decay
            # in all reservoirs
            
            #********
            #Set feed for downstream node
            
            # Set feed for next downstream node unless on the last node
            if i < self.nnodes() - 1:
                for k in range(self.NBedSizes + 1):
                    for j in range(self.NFlows):
                        self.Node[i + 1].Load.Qsjkfeed[j, k] = \
                            Node.Load.Qsjk[j, k]
                        for L in range(self.NTracers):
                            if self.NTracers > 0 and k == 0:
                                self.Node[i + 1].Load.TMudFeedj[j, L] = \
                                    Node.Load.TMudj[j, L]
                    self.Node[i + 1].Load.QsAvkFeed[k] = \
                        self.Node[i].Load.QsAvkLoad[k]
                    for L in range(self.NTracers):
                        if self.NTracers > 0 and k > 0:
                            self.Node[i + 1].Load.TBedFeedk[k, L] = \
                                Node.ActiveLayer.T[k, L]
                
    def UpdateSlope(self):
        """
        Subroutine to computes slope for each node from geometry of nodes in reach.
        """
        Node = self.Node
        for i in range(self.nnodes()):
            if i == self.nnodes() - 1:
                Node[i].Slope = (Node[i - 1].etabav - Node[i].etabav) / \
                    Node[i - 1].dxc
            else:
                Node[i].Slope = (Node[i].etabav - Node[i + 1].etabav) / \
                    Node[i].dxc
            #if Node[i].Slope < 0.: # Katie add to prevent crashing
             #   Node[i].Slope = .002
        
#        for i in range(self.nnodes()): # Katie add this calculation of slope to make bottommost node have same slope as second to last
#            if i < self.nnodes() - 1:
#                
#                Node[i].Slope = (Node[i].etabav - Node[i + 1].etabav) / \
#                    Node[i].dxc                
#            else:
#                Node[i].Slope = Node[i - 1].Slope
                
    def UpdateOutput(self, dt):
        """
        Calculates cumulative bed material output from bottom-most node and cumulative feed.
        Also calculates cumulative supply from banks.
        
        Katie add!
        
        Parameters
        ----------
        dt : float
            Timestep in seconds.
        """        
        
        Node = self.Node[3] # Elwha specific!
        #fluxrate = map(lambda x: x*dt, Node.Load.QsAvkLoad[0:self.NBedSizes + 1])
        fluxrate = list(map(lambda x: x*dt, Node.Load.QsAvkLoad[0:self.NBedSizes + 1]))
        self.CumulativeOutput[0:] = self.CumulativeOutput[0:] + fluxrate 
        
        Node = self.Node[0]
        #feedrate = map(lambda x: x*dt, Node.Load.QsAvkFeed[0:self.NBedSizes + 1])
        feedrate = list(map(lambda x: x*dt, Node.Load.QsAvkFeed[0:self.NBedSizes + 1]))
        self.CumulativeFeed[0:] = self.CumulativeFeed[0:] + feedrate

        # Sum bank supply from all nodes in Middle Elwha:  Elwha Specific!
        BankSup = np.zeros(self.NBedSizes + 1)
        BankSink = np.zeros(self.NBedSizes + 1)
        if self.nnodes() > 0:
            for i in range(4):
                Node = self.Node[i]
                BankSup = BankSup + Node.ActiveLayer.ExSed.InWidthChange*dt
                BankSink = BankSink + Node.ActiveLayer.ExSed.OutWidthChange*dt
        self.CumulativeBankSupply = self.CumulativeBankSupply + BankSup
        self.CumulativeBankSink = self.CumulativeBankSink + BankSink
        
    def Add_Nodes(self, NewNode, numnodes, etabav, xc, Outputvars, Validatevars, ValidateDvars, Outputfolder):    
        """
        Katie add this function to add nodes to the end of the reach.  Was
        developed for the dam removal, where aggradation on the final node
        rose above boundary condition for WSE, causing the model to crash.
        
        This function does two things:  it adds nodes to the bottom of the
        model, with the attributes of a control node.  Bed elevations are
        adjusted given the slope of the control node.  It also updates the
        output files with no data for the time elapsed before the nodes
        were added.
        
        Parameters
        ----------
        NewNode : :obj:`MAST_1D.clsNode`
            Control node object that whose attributes will be copied 
            into all of the new nodes.
        numnodes : int
            Number of new nodes to be added.
        etabav : float
            Initial etabav of the downstreammost node.
        xc : int
            Longitudinal coordinate of the downstreammost node.
        Outputvars : ???
            NEEDS TO BE DOCUMENTED
        Validatevars : ???
            NEEDS TO BE DOCUMENTED
        ValidateDvars : ???
            NEEDS TO BE DOCUMENTED
        Outputfolder : ???
            NEEDS TO BE DOCUMENTED
        """        
        
        #  Add specified number of nodes to downstream end of reach and 
        #  set the new bed elevations
        newnodes = []        
        for i in range(numnodes):
            Node = deepcopy(NewNode) # Need to fix this--can't deepcopy
            Node.etabav = etabav-Node.dxc*Node.Slope*(i+1)
            Node.xc = xc + Node.dxc*(i+1)
            Node.InitialBedElev = Node.etabav
            newnodes.append(Node)
        self.Node = self.Node + newnodes

        #  Add new nodes to output files
        for name in Outputvars:
            outputpath = os.path.join(os.pardir, Outputfolder, 'Out_' + name)
            if os.path.isfile(outputpath) == True:
                with open(outputpath, 'r+') as f:
                    lines = f.readlines()                
                    lenlines = len(lines[-1].split()) # Get number of entries already in file
                    for i in range(self.nnodes()):
                        if i >= self.nnodes()- numnodes:
                            outxc = int(self.Node[i].xc)
                            newrow = str(outxc) + '\t' + '\t' + ('-'+ '\t')*(lenlines-1) + '\n'
                            lines.append(newrow)
                    f.seek(0)
                    f.writelines(lines)

        #  Add new nodes to Validate files
        for name in Validatevars:
            outputpath = os.path.join(os.pardir, Outputfolder, 'OutValidate_' + name)
            if os.path.isfile(outputpath) == True:
                with open(outputpath, 'r+') as f:
                    lines = f.readlines()                
                    lenlines = len(lines[-1].split()) # Get number of entries already in file
                    for i in range(self.nnodes()):
                        if i >= self.nnodes()- numnodes:
                            outxc = int(self.Node[i].xc)
                            newrow = str(outxc) + '\t' + '\t' + ('-'+ '\t')*(lenlines-1) + '\n'
                            lines.append(newrow)
                    f.seek(0)
                    f.writelines(lines)

        #  Add new nodes to ValidateD files
        for name in ValidateDvars:
            outputpath = os.path.join(os.pardir, Outputfolder, 'OutValidateD_' + name)
            if os.path.isfile(outputpath) == True:
                with open(outputpath, 'r+') as f:
                    lines = f.readlines()                
                    lenlines = len(lines[-1].split()) # Get number of entries already in file
                    for i in range(self.nnodes()):
                        if i >= self.nnodes()- numnodes:
                            outxc = int(self.Node[i].xc)
                            newrow = str(outxc) + '\t' + '\t' + ('-'+ '\t')*(lenlines-1) + '\n'
                            lines.append(newrow)
                    f.seek(0)
                    f.writelines(lines)





