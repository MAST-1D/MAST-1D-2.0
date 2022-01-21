
from clsSubstratePairClass import clsSubstratePairClass
from clsReservoir import clsReservoir
from clsDurationCurve import clsDurationCurve
from clsLoad import clsLoad
from copy import deepcopy
import pdb
import numpy as np
from math import exp

class clsNode(object):
    """
    A segment of a river valley in which sediment transport is computed and
    and bed/bar sediment is conserved.
    
    The node object is the heart of MAST-1D. It include sediment
    storage reservoirs representing active layer of the channel bed, floodplain,
    and a series of substrate sediment storage reservoirs.  The 
    geometry of these reservoirs as well as the size distributions of their
    sediment are available within each node.  The node object 
    includes attributes used for representing bedrock or coarse, non-erodible
    lag layers.  It also stores parameters such as hydraulic roughness, a distribution
    of discharge, and hydraulic output.
    
    Parameters
    ----------
    NLayers : int
        Number of initial layers in substrate.
    NTracers : int
        Number of tracers.
    BinBdySizes : array_like(float)
        Sediment grain size at each bin boundary.
    NFlows : int
        Number of discharge bins in flow duration distribution. 
    
    Attributes
    ----------
    ActiveLayer : :obj:`MAST_1D.clsReservoir`
        A sediment storage reservoir representing the active layer of the channel.
    Substrate : array_like(:obj:`MAST_1D.clsSubstratePairClass`, length = NLayers)
        The sediment substrate.  There are two zones, one representing the average
        properties of sediment accessible to the channel if the channel incises and 
        another representing sediment that would be transferred to the floodplain and 
        could thus become available to the channel through lateral bank erosion. There
        are an arbitrary number of substrate layers. 
    Floodplain : :obj:`MAST_1D.clsReservoir`
        A sediment storage reservoir representing material stored adjacent to the
        channel and thus availble to the channel through lateral bank migration.  Sediment
        can be added to the floodplain through bar deposition or by overbank deposition.
    Load : :obj:`MAST_1D.clsLoad`
        Sediment moving downstream out of the node.
    Slope : float 
        Average down-channel slope. 
    H : float
        Channel depth.
    Bf : float 
        Floodplain width).
    Bc : float
        Channel width.
    cbank : float 
        Average lateral bank migration rate normal to down-channel direction for node.
    DeltaEtaB : float 
        Average rate of bed elevation change within channel.
    InitialBedElev : float 
        Initial bed elevation.
    Hpb : float
        Point bar thickness.  Used to compute lateral sediment flux to floodplain.
    Bcrate : float
        Channel widening rate.
    etabav : float
        Average channel bed elevation within node.
    ChSin : float 
        Channel sinusity.
    xc : float 
        Down channel coordinate.
    x : float
        Down valley coordinate.
    Dfjk : array_like(float) 
        Deposition rate on the floodplain in each bin of the
        low duration curve j and for each size k.  NOTE:  MAY 
        BE WORTH MOVING SO THIS IS A PROPERTY OF THE LOAD OBJECT.
    SLatSourcejk : array_like(float, dim = 2, length = (NFlows, NSizes))
        Lateral sources of sediment in m3/s in each 
        bin j of the flow duration in size class k.
    SLatSourceAv : array_like(float, length = NSizes)
        Mean annual lateral source of sediment in m3/s in size k.
    Dfav : array_like(float, length = NSizes) 
        Duration-averaged deposition rate on the floodplain in size k.  Usually 
        very small excepth in mud and sand sizes.
    FixedElev : bool 
        Flag to determine if the node's elevation is fixed 
        (e.g. bedrock or boundary condition).
    pTLatSource : ???
        Not sure if this attribute is still implemented--it should reprsent
        the probabilities of lateral erosion events.  There is a note saying that 
        pTLatSource could be adjusted as a function of time to 
        represent individual erosion events.  This would require changing the 
        duration curve at various times in the computation as well.
    TLatSourceAv : array_like(float, dim = 2, length = (NSizes, NTracers))
        Duration-averaged concentration of cosmogenic tracer nuclides in the 
        sediment coming from lateral sources in size k and tracer l.
    DC : :obj:`MAST_1D.clsDurationCurve` 
        Duration curve object used to store flow and sediment flux in each bin 
        of flow duration curve.
    Initialized : bool 
        Flag to determine if arrays have been redimensioned appropriately.
    lambdap : float
        Porosity of all sediment deposits.
    sigma : float
        Subsidence rate.
    dxc : float
        Down-channel distance to next node.
    Cfc : float 
        Friction coefficient for channel.
    Cff : float
        Friction coefficient for floodplain.
    FSandSusp : float
        Fraction sand load that is suspended.
    Kbar : float
        Parameter controlling fraction washload in point bar deposits.
    AlphaBar : float
        Parameter controlling similarity between bed material load and bar deposition.
    AlphaPartlyAlluvial : float
        Parameter controlling similarity between bed material load and deposition in 
        the active layer of a partly alluvial node.
    nc : float (read only)
        Manning's n for channel, including grain roughness, which is computed from
        the grain size distribution of the active layer, with form drag and sinuosity 
        multiplier included.
    ncAddons : float
        Form drag addition for Manning's n.
    ncMultiplier : float
        Sinuosity multiplier for Manning's n.
    nf : float
        Manning's n for floodplain)
    PointBarSubsurfaceGSD : :obj:`MAST_1D.clsGSD`
        Size distribution of subsurface of point bar.  Material with this distribution 
        is transferred into the floodplain by lateral channel shifting. It is a mixture
        of the active layer and the load.
    FkPointBarDepositAnnualAverage : ???
        Fraction in size class k in long-term average transfer of sediment
        to point bar.  COMPUTED IN __getattribute__. ADDITIONAL DOCUMENTATION 
        MAY BE APPROPRIATE.
    Flmud : float
        Floodplain number for mud.  Influences fraction of suspended mud flowing above
        floodplain level and across floodplain that gets stored in the floodplain through
        overbank deposition.
    Flbed : float
        Floodplain number for bed material sizes
    NLayers : int (read only)
        Number of layers of substrate.
    NSizes : int
        Number of sediment size bins.
    NTracers : int
        Number of tracers.
    NFlows : int
        Number of discharge bins in flow duration distribution.
    CumulativeBedChange : float (read only)
        Total channel bed elevation change since beginning of simulation.
    BarPavingRatio : float (read only)
        Ratio of D50 in active layer to D50 in subsurface of point bars.
    FractionAlluvial : float (read only)
        Ratio of sediment storage volume in active layer to maximum
        sediment storage volume of active layer.  
    CumulativeNarrowing : float
        Katie add--NEEDS DOCUMENTATION
    CumulativeWidening : float
        Katie add--NEEDS DOCUMENTATION
    CumulativeTotalAvulsionWidth : float
        Katie add--NEEDS DOCUMENTATION
    CumulativeTotFeed : float
        Katie add--NEEDS DOCUMENTATION
    CumulativeTotBedMaterialFeed : float
        Katie add--NEEDS DOCUMENTATION
    CumulativeWideningVolume : float
        Katie add--NEEDS DOCUMENTATION
    CumulativeWideningBedMaterialVolume : float
        Katie add--NEEDS DOCUMENTATION
    CumulativeChannelVolumeChange : float
        Katie add--NEEDS DOCUMENTATION
    PartlyAlluvial : boolean
        Katie add--NEEDS DOCUMENTATION
    Canyon : Boolean 
        Katie add--NEEDS DOCUMENTATION
    ControlSubstrate : :obj:`MAST_1D.clsSubstratePairClass`
        Stores original substrate in case incision is so large that a new layer is needed below 
        all the initial layers.    
    NarrowRate : float
        Channel narrowing rate.
    WideningRate : float
        Channel widening rate.
    ValleyMargin : float
        Katie add--it is the original fp + bc width.
    CumuTotVolume : float
        Katie add--NEEDS DOCUMENTATION 
    TotalSubstrateVolume : float. 
        To check mass conservation (Katie add).
    CumuTotalSubstrateCDeltaS : float
        Katie add--NEEDS DOCUMENTATION
    TotalCVol : float
        Katie add--NEEDS DOCUMENTATION
    CobbleMobility : float
        Katie add--NEEDS DOCUMENTATION  
    
    """
    
    Initialized = False
    g = 9.81
    
    def NLayers(self):
        return len(self.Substrate)
    
    def getCumulativeBedChange(self):
        return self.etabav - self.InitialBedElev
    
    CumulativeBedChange = property(getCumulativeBedChange)

    def __getattribute__(self, name):
        if name == 'PointBarSubsurfaceGSD':
            self.PointBarSubsurfaceGSD = deepcopy(self.ActiveLayer.GSD)
            
            for k in range(self.NSizes + 1):
                object.__getattribute__(self, 'PointBarSubsurfaceGSD').F[k] = \
                    self.FkPointBarDepositAnnualAverage[k]
            return object.__getattribute__(self, name)
        elif name == 'FkPointBarDepositAnnualAverage':
            """
            Note
            ----
            
            There are some interesting code comments in the __getattribute__ function.
            Basically, it looks like a call to return the 'PointBarSubsurfaceGSD'
            attribute will return the annual average of this attribute.
            
            
            In this function, Kbar defines the weighting of overall washload vs. bed 
            material load in new bar material and is defined as 
            kbar = pb/pw*(Fw,bar/Fb,bar), where pb and pw are the fractions bed
            material and wash load in the total load, respectively, and Fw,bar 
            and Fb,bar are the fractions wash load and bed material in the bar.
            For Kbar = 0, the bars are entirely bed material, and for 
            Kbar = infinity, the bars are entirely wash load. Note that in this 
            formulation, Kbar is equivalent to the ratio of BetaW/BetaB in the 
            Computers and Geosciences paper.  Also note that For bed material 
            load = 10% of total load and for bars that are 90% bed material, 
            k = 1/81. The ratio should be solved using observations of 
            mud content in the lower part of the floodplain (including in 
            oxbow-lake fill) along with estimates of annual average bed 
            material load (computed by the same equation used in the model) and
            annual mud load (computed from a rating curve).

            AlphaBar determines the size fractionation of the bed material 
            fraction of the new bar deposit, which is assumed to be a mixture 
            of material from the active layer and material from the load.  For 
            alpha = 1, the bed material (non-mud) portion of the point bar 
            deposit has the same size distribution as the active layer.  For 
            alpha = 0, the bed material (non-mud) portion of the point bar 
            deposit has the size distribution of the duration-averaged load.
        
            This is computed based on the size distributions within the load 
            for bed material but based on feed rate for mud (i.e. before mass 
            conservation for mud).  Might be possible to computed based only on
            feed or only on load.  In any case, may be sensitive to large 
            changes in feed size distribution, particularly if bed material 
            supply is reduced to zero.  Probably needs a way of reducing 
            migration rate in this case so that the channel does not run out of
            mud due to lateral exchange.
            """
            if self.Load.QsavBedTot > 0:
                FWashloadInPointBar = 1. - (1. / (1. + self.Kbar * \
                    self.Load.QsAvkFeed[0] / self.Load.QsavBedTot))
            else:
                FWashloadInPointBar = 1.
                #print('Bed material load is zero.  Washload fraction in ' + \
                #    'point bar set equal to 1.')
            
            Fkpointbar = [0.] * (self.NSizes + 1)
            Fkpointbar[0] = FWashloadInPointBar
            for k in range(1, self.NSizes + 1):
                Fkpointbar[k] = (1. - FWashloadInPointBar) * \
                    (self.ActiveLayer.GSD.F[k] * self.AlphaBar + \
                    self.Load.GSDBedloadAv.F[k] * (1. - self.AlphaBar))
                    
            return Fkpointbar
        else:
            return object.__getattribute__(self, name)
    
    def nGrain(self): # Grain roughness
        """
        Grain Roughness computed from size distribution in active layer.
        """
        # Katie make it so that mud is not in roughness 
        #calculation--for partly alluvial problem
        mud = self.ActiveLayer.GSD.F[0] 
        self.ActiveLayer.GSD.F[0] = 0.
        
        nGrain =  (self.ActiveLayer.GSD.D65) ** (1. / 6.) * 0.0146
        self.ActiveLayer.GSD.F[0] = mud
        #self.ActiveLayer.GSD.UpdateStatistics()

        return nGrain
    
    def nc(self):
        """
        Total channel roughness including form drag and multiplier for sinuosity.
        """
        #print 'nc'
        #print self.nGrain() + self.ncAddons
        return (self.nGrain() + self.ncAddons) * self.ncMultiplier
        
    def getBarPavingRatio(self):
        return self.ActiveLayer.GSD.D50 / self.PointBarSubsurfaceGSD.D50
    
    BarPavingRatio = property(getBarPavingRatio)
    
    def getBedPavingRatio(self):
        return self.ActiveLayer.GSD.D50 / self.Substrate[-1].C.GSD.D50
    
    BedPavingRatio = property(getBedPavingRatio)
    
    def getFractionAlluvial(self):
        return self.ActiveLayer.Volume / self.ActiveLayer.L / self.Bc / self.dxc
    
    FractionAlluvial = property(getFractionAlluvial)

    def SplitOrCombineSubstrate(self, LMinAfterSplit, LMinBeforeRemove, Spacing):
        """
        Method for updating the substrate if the top layer is too thick or too thin.
        
        If the top layer is thicker than allowed, the top layer is split into two 
        layers, both having the same size distribution.  If the top layer is thinner than
        allowed, it is combined with the layer below it.  Splitting or combining 
        substrate layers does not result in horizontal mixing between channel and floodplain
        zones of the substrate.  Also, if the number of subsrate layers drops to 1, adds
        a control substrate below the substrate.
        
        Parameters
        ----------
        LMinAfterSplit : float
            Minimum thickness (m) for a substrate node if a new substrate layer 
            is to be spawned (when aggrading). This It is necessary to split 
            substrate so as not to mix material too deeply.
        LMinBeforeRemove : float
            Minimum thickness (m) allowed for the top substrate layer if system is degrading.  
            If too thin, the layer is removed and its sediment is mixed with the next lower layer.             
        Spacing : float
            Thickness of all but top-most substrate layers (m).
        """
        
        # Katie--add this loop.  If the number of substrate layers gets down to
        # 1, a new layer is added to the bottom.  This could be a problem because
        # it will not have been exchanging with the floodplain constantly like the
        # other substrate layers, but I don't think it will make a big difference.
        if len(self.Substrate) == 1:
            self.Substrate = [deepcopy(self.ControlSubstrate)] + self.Substrate
        # Add a new layer if necessary
        if self.Substrate[-1].C.L > Spacing + LMinAfterSplit:
            OldL = self.Substrate[-1].C.L
            NewL = OldL - Spacing
            self.Substrate.append(clsSubstratePairClass())
            #self.Substrate[-1] = deepcopy(self.Substrate[-2])
            # Katie change--deepcopy keeps references for nested objects--need to define new object
            # and set it manually.
            self.Substrate[-1].C = clsReservoir(range(self.Substrate[-2].C.NSizes + 2), self.NTracers)
            self.Substrate[-1].F = clsReservoir(range(self.Substrate[-2].F.NSizes + 2), self.NTracers)
            self.Substrate[-1].C.GSD = deepcopy(self.Substrate[-2].C.GSD)
            self.Substrate[-1].F.GSD = deepcopy(self.Substrate[-2].F.GSD)
            self.Substrate[-2].C.L = Spacing
            self.Substrate[-2].F.L = Spacing
            self.Substrate[-2].C.Volume = Spacing * self.Bc * self.dxc
            self.Substrate[-2].F.Volume = Spacing * self.Bf * self.dxc / \
                self.ChSin
            self.Substrate[-1].C.L = NewL
            self.Substrate[-1].F.L = NewL
            self.Substrate[-1].C.Volume = NewL * self.Bc * self.dxc
            self.Substrate[-1].F.Volume = NewL * self.Bf * self.dxc / self.ChSin # Katie fix from 'Me.ChSin'

            # Katie add
            self.Substrate[-2].C.DeltaS = 0.
            self.Substrate[-2].F.DeltaS = 0.

            print('SubstrateSplit')
            #print self.NLayers()
        
        # Remove a layer if necessary
        if self.Substrate[-1].C.L < LMinBeforeRemove:
            NewL = self.Substrate[-1].C.L + self.Substrate[-2].C.L
            Substrate = self.Substrate[-2]
            for k in range(self.NSizes+1): # Katie add + 1--currently missing largest size class
                Substrate.C.GSD.F[k] = (self.Substrate[-1].C.GSD.F[k] * \
                    self.Substrate[-1].C.L + Substrate.C.GSD.F[k] * \
                    Substrate.C.L) / NewL
                Substrate.F.GSD.F[k] = (self.Substrate[-1].F.GSD.F[k] * \
                    self.Substrate[-1].F.L + Substrate.F.GSD.F[k] * \
                    Substrate.F.L) / NewL
                for m in range(self.NTracers):
                    Substrate.C.T[k, m] = (self.Substrate[-1].C.T[k, m] * \
                        self.Substrate[-1].C.L + Substrate.C.T[k, m] * \
                        Substrate.C.L) / NewL
                    Substrate.F.T[k, m] = (self.Substrate[-1].F.T[k, m] * \
                        self.Substrate[-1].F.L + Substrate.F.T[k, m] * \
                        Substrate.F.T[k, m] * Substrate.F.L) / NewL
            Substrate.C.L = NewL
            Substrate.F.L = NewL
            Substrate.C.Volume = NewL * self.Bc * self.dxc
            Substrate.F.Volume = NewL * self.Bf * self.dxc / self.ChSin
            del self.Substrate[-1]
            print('Substrate combined')
            #print self.NLayers()
    
    def __init__(self, NLayers, NTracers, BinBdySizes, NFlows):
        if not self.Initialized:
            self.NSizes = len(BinBdySizes) - 2
            self.NTracers = NTracers
            self.NFlows = NFlows
            self.Substrate = [clsSubstratePairClass() for i in range(NLayers)]
            self.ActiveLayer = clsReservoir(BinBdySizes, NTracers)
            self.Floodplain = clsReservoir(BinBdySizes, NTracers)
            self.Load = clsLoad(NFlows, BinBdySizes, NTracers)
            for Substrate in self.Substrate:
                Substrate.C = clsReservoir(BinBdySizes, NTracers)
                Substrate.F = clsReservoir(BinBdySizes, NTracers)
            self.DC = clsDurationCurve(NFlows)

            self.Dfjk = np.zeros((NFlows, self.NSizes + 1))
            self.Ssum = np.zeros((NFlows + 1, self.NSizes + 1))
            self.Sav = np.zeros(self.NSizes + 1)
            self.Dfav = np.zeros(self.NSizes + 1)
            self.STav = np.zeros((self.NSizes + 1, NTracers))
            self.Tlat = np.zeros((self.NSizes + 1, NTracers))
            self.Flatb = np.zeros(self.NSizes + 1)
            self.Flatu = np.zeros(self.NSizes + 1)
            self.SLatSourceAv = np.zeros(self.NSizes + 1)
            self.SLatSourcejk = np.zeros((NFlows, self.NSizes + 1))
            self.TLatSourceAv = np.zeros((self.NSizes + 1, NTracers))
            self.Bcrate = 0.
            self.DeltaEtaB = 0. # Katie add
            self.CumulativeNarrowing = 0. # Katie add
            self.CumulativeWidening = 0. # Katie add
            self.CumulativeTotalAvulsionWidth = 0.
            self.CumulativeTotFeed = 0. # Katie add
            self.CumulativeTotBedMaterialFeed = 0. # Katie add
            self.CumulativeWideningVolume = 0. # Katie add
            self.CumulativeWideningBedMaterialVolume = 0. # Katie add
            self.CumulativeChannelVolumeChange = 0. # Katie add
            self.PartlyAlluvial = False # Katie add
            self.Canyon = False # Katie add
            self.ControlSubstrate = clsSubstratePairClass() # Katie add to hold original substrate (see Add/Remove Substrate function)
            self.NarrowRate = 0. # Katie add these to separate BcRate into narrowing and widening
            self.WidenRate = 0.
            self.ValleyMargin = 0. # Katie add--it is the original fp + bc width.
            self.CumuTotVolume = 0. # Katie add            
            self.TotalSubstrateVolume = 0. # Katie add--to check mass conservation
            self.CumuTotalSubstrateCDeltaS = 0. # Katie add
            self.CumuTotalSubstrateFDeltaS = 0. # Katie add
            self.TotalCVol = 0. # Katie add
            self.CobbleMobility = 0.
            
            self.Initialized = True
        else:
            raise RuntimeError('Tried to initiate clsNode twice.')
    
    def EquilibriumMudFloodplainNumber(self, DeltaEta, cbank):
        """
        Computes floodplain number required to achieve long-term deposition
        rate that produces a given overbank sediment thickness.
        
        This function can be used to compute the floodplain number for mud at 
        perfect equilibrium, when the bed is not changing and there is no net 
        divergence in load.  It must be called after the hydraulics have been 
        updated and mud feed has been specified for the node.
        The variable DeltaEta should represent the mud fraction of the 
        thickness of the overbank sediment deposit.
        
        Parameters
        ----------
        DeltaEta : float
            Average thickness of overbank fines at eroding banks, 
            assuming long-term bed elevation is constant.
        cbank : float
            Long-term average bank migration rate. 
        
        Return: float
            Equilibrium floodplain number for mud.
        """
        
        # Store old values of Fl, Dfjk(j,0) and Dfav(0)
        try:
            OldvalFl = self.Flmud
        except AttributeError:
            OldvalFl = 0.
        
        OldvalDfav0 = self.Dfav[0]
        OldvalDfj0 = np.array(self.Dfjk[:, 0])
        
        # Compute the total overbank deposition rate if Flmud is equal to one,
        # then calculate the fraction this must be reduced in order to equal
        # the observed erosion rate of the overbank deposit.  
        # This is the equilibrium floodplain number.
        self.Flmud = 1.
        self.UpdateDMudj()

        EquilibriumMudFloodplainNumber = DeltaEta * cbank / self.Dfav[0]
        
        self.Flmud = OldvalFl

        self.Dfav[0] = OldvalDfav0
        self.Dfjk[:, 0] = OldvalDfj0

        return EquilibriumMudFloodplainNumber
    
    def NormalChannelDepthAndDischarge(self, Qw, threshold, Manning): #  Input variables are now the same between VBA and this Python version (Katie)
        """
        Calculates steady uniform flow hydraulics.
        
        This method computes the normal water depth as a function of 
        channel width, floodplain thickness, floodplain width, friction 
        coefficients in channel and on floodplain, and slope, and partitions
        flow into channel and floodplain zones. It can perform calculation
        either using Manning's equation (for Manning = True) or Chezy
        equation (For Manning = False)
        
        Parameters
        ----------
            Qw : float
                Discharge (m3/s)
            threshold : float
                allowable error threshold.
            Manning : bool
                Flag indicating whether to use Manning's equation (if True) or
                Chezy equation (if False).
            
        Return: float
            Normal depth in channel (m).
        """

        g = 9.81
        Czc = 1 / self.Cfc ** 0.5
        Czf = 1/ self.Cff ** 0.5


        if Manning:
            #print self.nc()
            #print Qw
            #print self.Bc
            #print self.Slope
            HcGuess = (self.nc() * Qw / (self.Bc * self.Slope ** 0.5)) ** \
                (3. / 5.)

        else:
            HcGuess = (Qw / (self.Bc * Czc * (g * self.Slope) ** 0.5)) ** \
                (1 / 1.5)
        if HcGuess < (self.Floodplain.L - self.ActiveLayer.L):
            Hc = HcGuess
            Qc = Qw

        else:
            #print 'oogala'
            # iterate to solution using Newton's method
            dHc = 1e-05
            error = threshold + 1.
            counter = 0
            while error >= threshold:
                if Manning:
                    QwGuess = self.ManningDischarge(HcGuess, self.Bc, \
                        self.nc(), self.Slope, self.Bf, self.nf, \
                        self.Floodplain.L - self.ActiveLayer.L, self.Slope * \
                        self.ChSin)
                else:
                    QwGuess = self.ChezyDischarge(HcGuess, self.Bc, Czc, \
                        self.Slope, self.Bf, Czf, self.Floodplain.L - \
                        self.ActiveLayer.L, self.Slope * self.Chsin)
                error = QwGuess - Qw
                
                if Manning:
                    QwGuessPlusDQw = self.ManningDischarge(HcGuess + dHc, \
                        self.Bc, self.nc(), self.Slope, self.Bf, self.nf, \
                        self.Floodplain.L - self.ActiveLayer.L, self.Slope * \
                        self.ChSin)
                else:
                    QwGuessPlusDQw = self.ChezyDischarge(HcGuess + dHc, \
                        self.Bc, Czc, self.Slope, self.Bf, Czf, \
                        self.Floodplain.L - self.ActiveLayer.L, self.Slope * \
                        self.ChSin)

                dQw = QwGuessPlusDQw - QwGuess
                HcGuess = HcGuess - error * dHc / dQw
                counter += 1
                if counter > 100:
                    raise RuntimeError('iteration for channel flow did not \
                        converge in clsNode')
            Hc = HcGuess

            if Manning:
                Qc = self.ManningDischarge(Hc, self.Bc, self.nc(), self.Slope,\
                    0., self.nf, 0., 0.)

            else:
                Qc = self.ChezyDischarge(Hc, self.Bc, Czc, self.Slope, 0., \
                    0., 0., 0.)

        return np.array([Hc, Qc])
    
    def ChezyDischarge(self, H, Bc, Czc, Sc, Bf, Czf, Tf, Sf):
        """
        Estimates overall discharge in channel/floodplain complex
        using Chezy equation.
                
        Parameters
        ----------
        H : float
            Flow Depth in channel (m)
        Bc : float
            Channel width (m)
        Czc : float
            Chezy coefficient in channel zone.
        Sc : float
            Average Channel Slope
        Bf : float
            Floodplain width (m)
        Czf : float
            Chezy coefficient in floodplain zone.
        Tf : float
            Floodplain thickness (same as channel depth).
        Sf : float
            Averate down-channel slope in floodplain zone.
        
        Returns
        -------
        float
            Discharge (m3/s)
        """
        g = 9.81
        if H > Tf:
            return Bc * Czc * H * (g * H * Sc) ** 0.5 + Bf * Czf * (H - Tf) * \
                (g * (H - Tf) * Sf) ** 0.5
        else:
            return Bc * Czc * H * (g * H * Sc) ** 0.5
    
    def ManningDischarge(self, H, Bc, nc, Sc, Bf, nf, Tf, Sf):
        """
        Estimates overall discharge in channel/floodplain complex
        using Manning equation.       
        
        Parameters
        ----------
        H : float
            Flow Depth in channel (m).
        Bc : float
            Channel width (m).
        nc : float
            Manning's n for channel zone.
        Sc : float
            Average Channel Slope.
        Bf : float
            Floodplain width (m).
        nf : float
            Mannin's n for floodplain zone.
        Tf : float
            Floodplain thickness (same as channel depth).
        Sf : float
            Averate down-channel slope in floodplain zone.
        
        Returns
        -------
        float
            Discharge (m3/s)
        """

        if H > Tf:
            return 1. / nc * Bc * H * H ** (2. / 3.) * Sc ** (1. / 2.) + 1. \
                / nf * Bf * (H - Tf) * (H - Tf) ** (2. / 3.) * Sf ** (1. / 2.)
        else:
            return 1. / nc * Bc * H * H ** (2. / 3.) * Sc ** (1. / 2.)
    
    def UpdateDepthAndDischargeAtAllFlows(self, Manning):
        """
        Method to compute and update hydraulic parameters using steady uniform
        flow approximation.  Performs computation for all discharges in discharge distribution.
        
        Parameters
        ----------
        Manning : bool
            Flag indicating whether to use Manning's equation (if True) or
            Chezy equation (if False).
        """
        
        for j in range(self.DC.NFlows()):
            
            ChannelHQ = self.NormalChannelDepthAndDischarge(self.DC.Qw[j], \
                0.01, Manning)
            DC = self.DC
            DC.Hc[j] = ChannelHQ[0]
            DC.Qwc[j] = ChannelHQ[1]
            DC.Qwf[j] = DC.Qw[j] - DC.Qwc[j]
            DC.Uc[j] = DC.Qwc[j] / (DC.Hc[j] * self.Bc)
            if DC.Hc[j] > self.Floodplain.L - self.ActiveLayer.L:
                DC.Hf[j] = DC.Hc[j] - (self.Floodplain.L - self.ActiveLayer.L)
                DC.Uf[j] = DC.Qwf[j] / (DC.Hf[j] * self.Bf)
            else:
                DC.Hf[j] = 0.
                DC.Uf[j] = 0.
            DC.Sf[j] = self.Slope
            
    def UpdateDepthAndDischargeAtOneFlow(self, Manning, Qw):
        """
        Function to compute velocity and friction slope using steady uniform
        flow approximation.  Performs computation for one discharge. Katie add for 
        width change function with flow duration curve.
        
        Parameters
        ----------
        Manning : bool
            Flag indicating whether to use Manning's equation (if True) or
            Chezy equation (if False).
        Qw : float  
            Discharge (m3/s)
            
        Returns
        -------
        Uc : float
            Channel Velocity
        Sf : float  
            Average slope of node (note this is simply the node's slope property)
        
        """
                    
        ChannelHQ = self.NormalChannelDepthAndDischarge(Qw, \
            0.01, Manning)

        Hc = ChannelHQ[0]
        Qwc = ChannelHQ[1]
        Qwf = Qw - Qwc
        Uc = Qwc / (Hc * self.Bc)
        if Hc > self.Floodplain.L - self.ActiveLayer.L:
            Hf = Hc - (self.Floodplain.L - self.ActiveLayer.L)
            Uf = Qwf / (Hf * self.Bf)
        else:
            Hf = 0.
            Uf = 0.
        Sf = self.Slope
        
        return Uc, Sf
    
    def UpdateLateralSedFluxes(self):
        """
        Method that computes and updates the size-specific lateral fluxes for the node.
        
        These are volume fluxes and do not include overbank deposition, which
        is handled as a net source term for the floodplian and as a net sink 
        term for the water column (for mud) or the active layer (for sand)
        """
        
        # Set up even outgoing fluxes for substrate layers due to lateral 
        # boundary movement of channel zone
        for m in range(self.NLayers()):
            for k in range(self.NSizes + 1):
                 #For channel zone # Katie:  comment this section out if you don't want active layer and 
                # floodplain substrates mixing.
                self.Substrate[m].C.ExSed.OutMigration[k] = \
                    (1. - self.lambdap) * self.Substrate[m].C.GSD.F[k] * \
                    self.Substrate[m].C.L * self.cbank * self.dxc
                
                # The following ensures that if channel is narrowing, the 
                # channel zone substrate is exporting sediment to floodplain 
                # zone
                #if self.Bcrate < 0.:

                if self.NarrowRate != 0.: # Katie add
                    self.Substrate[m].C.ExSed.OutWidthChange[k] = \
                        -(1. - self.lambdap) * self.Substrate[m].C.GSD.F[k] * \
                        self.Substrate[m].C.L * self.NarrowRate * self.dxc # Katie:  good:  this produces a positive number + change BcRate to NarrowRate
                else:
                    self.Substrate[m].C.ExSed.OutWidthChange[k] = 0.
                # For floodplain zone # Katie comment out
                self.Substrate[m].F.ExSed.OutMigration[k] = \
                    (1. - self.lambdap) * self.Substrate[m].F.GSD.F[k] * \
                    self.Substrate[m].F.L * self.cbank * self.dxc # Katie unindent beginning of this line bloc
                
                # The following ensures that if channel is widening, the
                # floodplain zone substrate is exporting sediment to channel
                # zone
                #if self.Bcrate > 0.:
                if self.WidenRate != 0.: # Katie add
                    self.Substrate[m].F.ExSed.OutWidthChange[k] = \
                        (1. - self.lambdap) * self.Substrate[m].F.GSD.F[k] \
                        * self.Substrate[m].F.L * self.WidenRate * self.dxc # Katie get rid of negative sign at the beginning of this code bloc--number needs to be positive for UpdateFandT function. + change BcRate to WidenRate
                else:
                    self.Substrate[m].F.ExSed.OutWidthChange[k] = 0.
        
        # Set up outgoing fluxes for floodplain, active layer, and water column
        # (i.e. load) due to lateral boundary movement of channel zone
        
        for k in range(self.NSizes + 1):
            # active layer flux includes flux due to point bar deposition 
            # and due to lateral movement of boundary through thickness of 
            # active layer fraction bed material in new point bar material
            
            if k == 0:
            # there is no mud moved to the floodplain from the active 
            # layer.  Mud in point bar deposits comes from water column and
            # is accounted for in the conservation of mud equation.
                self.ActiveLayer.ExSed.OutMigration[k] = 0.
                self.Load.ExSed.OutMigration[k] = (1. - self.lambdap) * \
                    self.FkPointBarDepositAnnualAverage[k] * self.Hpb * \
                    self.cbank * self.dxc
                if self.DC.Uc[-1] < .5: # Katie add so no migration occurs in reservoirs where there is negligible bedload transport
                    self.Load.ExSed.OutMigration[k] = 0.
                    
            else:
                self.ActiveLayer.ExSed.OutMigration[k] = \
                    (1. - self.lambdap) * (self.ActiveLayer.GSD.F[k] * \
                    self.ActiveLayer.L + \
                    self.FkPointBarDepositAnnualAverage[k] * self.Hpb) * \
                    self.cbank * self.dxc
                if self.DC.Uc[-1] < .5: # Katie add so no migration occurs in reservoirs where there is negligible bedload transport
                    self.ActiveLayer.ExSed.OutMigration[k] = 0.
        
            Floodplain = self.Floodplain # Katie:  I think the same mud is both going to the active layer (here) and the load (Load.MassConservationForMudload)
            Floodplain.ExSed.OutMigration[k] = (1. - self.lambdap) * \
                Floodplain.GSD.F[k] * Floodplain.L * self.cbank * self.dxc
            if self.DC.Uc[-1] < .5: # Katie add so no migration occurs in reservoirs where there is negligible bedload transport
                self.Floodplain.ExSed.OutMigration[k] = 0.
        
        # Set up outgoing fluxes for floodplain, active layer, and water column
        # (i.e. load) due to width change of channel zone
        
        for k in range(self.NSizes + 1):
            #if self.Bcrate < 0.:
            if self.NarrowRate != 0.:
                # Assumes narrowing only occurs by point bar deposition.  Flux
                # also accounts for movement of vertical boundary between 
                # active layer and floodplain.
                # This routine still needs to include mud in the floodplain 
                # since the (alpha) term only allows for bed material sieving.
    
                if k == 0:
                    # there is no mud moved to the floodplain from the active 
                    # layer.  Mud in point bar deposits comes from water column
                    # and is accounted for in the conservation of mud equation.
                    # This routine assumes all narrowing is by point bar 
                    # deposition.
    
                    # NOTE: Katie changed it so that mud is moved to the floodplain
                    # from the active layer.  This is because some mud now
                    # enters the active layer during widening when the active layer
                    # expands laterally.
                    self.ActiveLayer.ExSed.OutWidthChange[k] = 0.
                    self.Load.ExSed.OutWidthChange[k] = (1. - self.lambdap) * \
                        self.FkPointBarDepositAnnualAverage[k] * self.Hpb * \
                        self.NarrowRate * self.dxc # Katie change BcRate to NarrowRate
                else:
                    self.ActiveLayer.ExSed.OutWidthChange[k] = \
                        -(1 - self.lambdap) * (self.ActiveLayer.GSD.F[k] * \
                        self.ActiveLayer.L + \
                        self.FkPointBarDepositAnnualAverage[k] * self.Hpb) * \
                        self.NarrowRate * self.dxc # Katie add negative sign to make the whole value positive for the UpdateFandT function + change BcRate to NarrowRate

            
            else: 
                self.ActiveLayer.ExSed.OutWidthChange[k] = 0. #if widening,
                    # there is not flux from active layer to floodplain due to
                    # narrowing
    
            Floodplain = self.Floodplain
            #if self.Bcrate > 0.:
            if self.WidenRate != 0.: # Katie add
                # Assumes narrowing only occurs by point bar deposition. 
                # Flux also accounts for movement of vertical boundary 
                # between active layer and floodplain.
                Floodplain.ExSed.OutWidthChange[k] = (1. - self.lambdap) \
                    * Floodplain.GSD.F[k] * (Floodplain.L-self.ActiveLayer.L) * self.WidenRate * \
                    self.dxc # Katie:  this is already positive and doesn't need the negtive sign (see line 546) + change Bcrate to WidenRate
                    # Katie--have to add in lower part of floodplain supply after exner equation so that the active layer volume increases and
                    # isn't compensated for by bed aggradation.
            else:
                Floodplain.ExSed.OutWidthChange[k] = 0. # if narrowing, there 
                    # is no flux from floodplain to active layer due to
                    # widening

        # Set up incoming fluxes due to lateral boundary movement.  In all 
        # cases these are equal to one of the already computed outgoing fluxes.
        for k in range(self.NSizes + 1):
            if k == 0:
                self.Load.ExSed.InMigration[k] = \
                    self.Floodplain.ExSed.OutMigration[k]
                self.Load.ExSed.InWidthChange[k] = \
                    self.Floodplain.ExSed.OutWidthChange[k]
            else:
                self.ActiveLayer.ExSed.InMigration[k] = \
                    self.Floodplain.ExSed.OutMigration[k]
                self.ActiveLayer.ExSed.InWidthChange[k] = \
                    self.Floodplain.ExSed.OutWidthChange[k]
            self.Floodplain.ExSed.InMigration[k] = \
                self.ActiveLayer.ExSed.OutMigration[k]
            self.Floodplain.ExSed.InWidthChange[k] = \
                self.ActiveLayer.ExSed.OutWidthChange[k]
            for m in range(self.NLayers()): # Katie comment out
                self.Substrate[m].C.ExSed.InMigration[k] = \
                    self.Substrate[m].F.ExSed.OutMigration[k]
                self.Substrate[m].F.ExSed.InMigration[k] = \
                    self.Substrate[m].C.ExSed.OutMigration[k]
                self.Substrate[m].C.ExSed.InWidthChange[k] = \
                    self.Substrate[m].F.ExSed.OutWidthChange[k]
                self.Substrate[m].F.ExSed.InWidthChange[k] = \
                    self.Substrate[m].C.ExSed.OutWidthChange[k]

    def UpdateLateralTracerConcentrations(self):
        """
        Method that computes and updates concentrations of tracers in lateral fluxes.
        
        This subroutine is necessary because the subroutine that updates tracer
        concentrations based on mass conservation at each node operates as a
        method of clsReservoir and thus does not have access to the tracer 
        concentrations in the adjacent reservoirs. Consequently, these must be 
        set up at node level.
    
        Note
        ----
        This does not need to be called before calling mass conservaiton on the
        load because all washload tracer concentrations are computed using the 
        tracer concentration of the washload feed.
        """
    
        for L in range(self.NTracers):
            # Set up even outgoing tracer concentrations for substrate layers
            # due to lateral boundary movement of channel zone
            for m in range(self.NLayers()):
                for k in range(self.NSizes + 1):
                    self.Substrate[m].C.ExTracer[L].OutMigration[k] = \
                        self.Substrate[m].C.T[k, L]
                    self.Substrate[m].C.ExTracer[L].OutWidthChange[k] = \
                        self.Substrate[m].C.T[k, L]
                    self.Substrate[m].F.ExTracer[L].OutMigration[k] = \
                        self.Substrate[m].F.T[k, L]
                    self.Substrate[m].F.ExTracer[L].OutWidthChange[k] = \
                        self.Substrate[m].F.T[k, L]
        
            # Set up outgoing fluxes for floodplain and active layer.
            for k in range(self.NSizes + 1):
                if k == 0:
                    self.Load.ExTracer[L].OutMigration[0] = \
                        self.Load.TTemporalAvMudFeed[L]
                else:
                    self.ActiveLayer.ExTracer[L].OutMigration[k] = \
                        self.ActiveLayer.T[k, L]
                self.Floodplain.ExTracer[L].OutMigration[k] = \
                    self.Floodplain.T[k, L]
        
            # Set up outgoing fluxes for floodplain, active layer, and water
            # column (i.e. load) due to width change of channel zone
            for k in range(self.NSizes + 1):
                #if self.Bcrate < 0.:
                if self.NarrowRate != 0.: # Katie add
                    if k == 0:
                        self.Load.ExTracer[L].OutWidthChange[k] = \
                            self.Load.TTemporalAvMudFeed[L]
                    else:
                        self.ActiveLayer.ExTracer[L].OutWidthChange[k] = \
                            self.ActiveLayer.T[k, L]
                Floodplain = self.Floodplain
                #if self.Bcrate > 0.:
                if self.WidenRate != 0.: # Katie add
                    Floodplain.ExTracer[L].OutWidthChange[k] = \
                        Floodplain.T[k, L]
                else:
                    Floodplain.ExTracer[L].OutWidthChange[k] = 0. # if
                        # narrowing, there is no flux from floodplain to active
                        # layer due to widening
        
            # Set up incoming concentrations due to lateral boundary movement.
            # In all cases these are equal to one of the already computed 
            # outgoing fluxes.
            for k in range(self.NSizes + 1):
                if k == 0:
                    self.Load.ExTracer[L].InMigration[0] = \
                        self.Floodplain.ExTracer[L].OutMigration[0]
                    self.Load.ExTracer[L].InWidthChange[0] = \
                        self.Floodplain.ExTracer[L].OutWidthChange[0]
                    self.Floodplain.ExTracer[L].InMigration[0] = \
                        self.Load.ExTracer[L].OutMigration[0]
                    self.Floodplain.ExTracer[L].InWidthChange[0] = \
                        self.Load.ExTracer[L].OutWidthChange[0]
                
                    self.Substrate[-1].F.ExTracer[L].InMigration[0] = 0. 
                    self.Substrate[-1].F.ExTracer[L].InWidthChange[0] = 0.
                        # These could be set to the duration-averaged tracer
                        # concentration, but for now, assume no mud infiltrates
                        # through active layer.
                    
                    # Mud could be allowed to enter the substrate through the 
                    # channel due either to hyphoreic exchange in which case a
                    # hyphoreic flux would need to be computed and a trap 
                    # efficiency of some sort would need to be specified for 
                    # the gravel.  This would probably need to be applied for
                    # each bin of the duration curve.  Alternatively, it would
                    # be possible to allow the active layer to change size over
                    # time, as a function of the flow.  The difference between
                    # active layer thickness for the high and low flow could be
                    # used to approximate this exchange thickness.  However, it
                    # is not clear what mechanism would be responsible for
                    # movement of fines deeper into the substrate. Probably 
                    # best to use a power function distribution for hyphoreic
                    # residence time and see what can be done with that, 
                    # settling velocity, and pore size in the upper substrate.
                    # The residence time would need to depend on fines 
                    # content...
                else:
                    self.ActiveLayer.ExTracer[L].InMigration[k] = \
                        self.Floodplain.ExTracer[L].OutMigration[k]
                    self.ActiveLayer.ExTracer[L].InWidthChange[k] = \
                        self.Floodplain.ExTracer[L].OutWidthChange[k]
                    self.Floodplain.ExTracer[L].InMigration[k] = \
                        self.ActiveLayer.ExTracer[L].OutMigration[k]
                    self.Floodplain.ExTracer[L].InWidthChange[k] = \
                        self.ActiveLayer.ExTracer[L].OutWidthChange[k]
                for m in range(self.NLayers()):
                    self.Substrate[m].C.ExTracer[L].InMigration[k] = \
                        self.Substrate[m].F.ExTracer[L].OutMigration[k]
                    self.Substrate[m].F.ExTracer[L].InMigration[k] = \
                        self.Substrate[m].C.ExTracer[L].OutMigration[k]
                    self.Substrate[m].C.ExTracer[L].InWidthChange[k] = \
                        self.Substrate[m].F.ExTracer[L].OutWidthChange[k]
                    self.Substrate[m].F.ExTracer[L].InWidthChange[k] = \
                        self.Substrate[m].C.ExTracer[L].OutWidthChange[k]

    def UpdateVerticalExchangeSedFluxes(self, BedAggradationRate, alpha):
        """
        Method that computes and updates fluxes associated with vertical bed movement.
        
        This subroutine considers only fluxes associated with boundary movement
        of bottom of active layer/floodplain. Net fluxes associated with 
        overbank deposition are handled in UpdateSedimentSourcesAndSinks. Results 
        in a source of mud to the water column from the 
        bed if the channel is degrading into a substrate that contains mud. 
        Only the duration-averaged supply rate of mud to the water column is 
        specified.
        
        Parameters
        ----------
        BedAggredationRate : float
            Rate at which bed is aggrading.
        alpha : float
            Coefficient determining the
            mixture between bed material and load that gets transferred 
            to substrate when system is aggrading.
        """
        #print 'Agg' + str(BedAggradationRate)
        for k in range(self.NSizes + 1):
            
            # Outgoing flux from active layer
            # Katie comment out            
            if k == 0:
                # no mud infiltrates for now, but this should be updated if 
                # possible
                self.ActiveLayer.ExSed.OutVerticalChange[0] = 0.
            else:
            # Katie add normalization
                FNorm = self.ActiveLayer.GSD.F[k]/sum(self.ActiveLayer.GSD.F[1:])
                if BedAggradationRate > 0:
                    self.ActiveLayer.ExSed.OutVerticalChange[k] = (alpha * \
                        FNorm + (1. - alpha) * \
                        self.Load.GSDBedloadAv.F[k]) * BedAggradationRate * \
                        self.Bc * self.dxc * (1. - self.lambdap)
                else:
                    self.ActiveLayer.ExSed.OutVerticalChange[k] = 0.
        
            # Outgoing flux from floodplain
            if BedAggradationRate > 0.:
                self.Floodplain.ExSed.OutVerticalChange[k] = \
                    self.Floodplain.GSD.F[k] * BedAggradationRate * self.Bf * \
                    self.dxc / self.ChSin * (1. - self.lambdap)
            else:
                self.Floodplain.ExSed.OutVerticalChange[k] = 0. # Katie change from ActiveLayer to Floodplain...bug.
            
            # Outgoing flux from upper substrate layer floodplain zone # Katie comment out to get rid of floodplain/substrate coupling
            if BedAggradationRate < 0.:
                self.Substrate[-1].F.ExSed.OutVerticalChange[k] = \
                    -self.Substrate[-1].F.GSD.F[k] * \
                    BedAggradationRate * self.Bf * self.dxc / self.ChSin * \
                    (1. - self.lambdap)
            else:
                self.Substrate[-1].F.ExSed.OutVerticalChange[k] = 0.
            
            # Outgoing flux from upper substrate channel zone
            if BedAggradationRate < 0.:
                 # Katie add normalization
                FNorm = self.Substrate[-1].C.GSD.F[k]/sum(self.Substrate[-1].C.GSD.F[1:])
                self.Substrate[-1].C.ExSed.OutVerticalChange[k] = \
                    -FNorm * \
                    BedAggradationRate * self.Bc * self.dxc * \
                    (1. - self.lambdap)
            else:
                self.Substrate[-1].C.ExSed.OutVerticalChange[k] = 0.
            
            
            # Setup incoming fluxes
            
            self.Substrate[-1].F.ExSed.InVerticalChange[k] = \
                self.Floodplain.ExSed.OutVerticalChange[k]
            self.Substrate[-1].C.ExSed.InVerticalChange[k] = \
                deepcopy(self.ActiveLayer.ExSed.OutVerticalChange[k])
            self.Floodplain.ExSed.InVerticalChange[k] = \
                self.Substrate[-1].F.ExSed.OutVerticalChange[k]
            self.ActiveLayer.ExSed.InVerticalChange[k] = \
                deepcopy(self.Substrate[-1].C.ExSed.OutVerticalChange[k])

        # set up mud fluxes to water column due to incision
        self.Load.ExSed.InVerticalChange[0] = \
            self.Substrate[-1].C.ExSed.OutVerticalChange[0]
       
        if self.FixedElev: # Katie add
            if sum(self.ActiveLayer.SourceFeedSed) > sum(self.ActiveLayer.SinkLoadSed):
               self.Load.ExSed.InVerticalChange[0] = self.ActiveLayer.SinkLoadSed[0]-self.ActiveLayer.SourceFeedSed[0] # Katie add
            else: 
                self.Load.ExSed.InVerticalChange[0] = -self.ActiveLayer.GSD.F[0]*(sum(self.ActiveLayer.SourceFeedSed)-sum(self.ActiveLayer.SinkLoadSed))                
        # Katie add:
        #self.Load.MudEnt = self.Substrate[-1].C.ExSed.OutVerticalChange[0]
        self.ActiveLayer.ExSed.InVerticalChange[0] = 0.
        
        # this assumes mud does NOT infiltrate into bed if bed is aggrading
        

    def UpdateVerticalExchangeTracerConcentrations(self):
        """
        Method that computes and updates concentrations of tracers in vertical fluxes.
        
        This method is necessary because the method that updates tracer
        concentrations based on mass conservation at each node operates is a 
        method of clsReservoir and thus does not have access to the tracer 
        concentrations in the adjacent reservoirs.
        Consequently, these must be set up at node level.
        """
        
        for L in range(self.NTracers):
            for k in range(self.NSizes + 1):
                
                # Outgoing Concentrations
                
                self.ActiveLayer.ExTracer[L].OutVerticalChange[k] = \
                    self.ActiveLayer.T[k, L]
                self.Floodplain.ExTracer[L].OutVerticalChange[k] = \
                    self.Floodplain.T[k, L]
                self.Substrate[-1].F.ExTracer[L].OutVerticalChange[k] = \
                    self.Substrate[-1].F.T[k, L]
                self.Substrate[-1].C.ExTracer[L].OutVerticalChange[k] = \
                    self.Substrate[-1].C.T[k, L]
                
                # Setup incoming Concentrations
                
                self.Substrate[-1].F.ExTracer[L].InVerticalChange[k] = \
                    self.Floodplain.ExTracer[L].OutVerticalChange[k]
                self.Substrate[-1].C.ExTracer[L].InVerticalChange[k] = \
                    self.ActiveLayer.ExTracer[L].OutVerticalChange[k]
                self.Floodplain.ExTracer[L].InVerticalChange[k] = \
                    self.Substrate[-1].F.ExTracer[L].OutVerticalChange[k]
                self.ActiveLayer.ExTracer[0].InVerticalChange[k] = \
                    self.Substrate[-1].C.ExTracer[L].OutVerticalChange[k]
            
            # set up mud fluxes to water column due to incision
            self.Load.ExTracer[L].InVerticalChange[0] = \
                self.Substrate[-1].C.ExTracer[L].OutVerticalChange[0]
    
    def UpdateExtraWidthFluxes(self, WidenRate, NarrowRate, AggRate):
        """
        Katie add. Method that updates sediment fluxes in the floodplain and
        active layer reservoir objects--to reflect width change.
        """
        # Katie add        
        
        # Now add extra bit of floodplain adjacent to active layer to sediment
        # exchanges; after DeltaEtaB is calculated.  Ensures that the active layer
        # volume increases when the width increases.  Here, mud from fp can enter active layer
        
        # Corresponding tracer function not yet written.
        for k in range(self.NSizes + 1):
            if k == 0:
                NewALMud = self.Floodplain.GSD.F[k]*self.ActiveLayer.L * self.dxc * WidenRate*(1 - self.lambdap)
                self.Floodplain.ExSed.OutWidthChange[k] += NewALMud
                self.ActiveLayer.ExSed.InWidthChange[k] = NewALMud
            else:
                self.Floodplain.ExSed.OutWidthChange[k] += self.Floodplain.GSD.F[k]*self.ActiveLayer.L * self.dxc * WidenRate*(1 - self.lambdap)
                self.ActiveLayer.ExSed.InWidthChange[k] = self.Floodplain.ExSed.OutWidthChange[k]
                
            
            # These not needed.
            #self.Substrate[-1].F.ExSed.OutWidthChange[k] += self.Substrate[-1].F.GSD.F[k]*AggRate*self.dxc*WidenRate
            #self.Substrate[-1].C.ExSed.InWidthChange[k] = self.Substrate[-1].F.ExSed.OutWidthChange[k]
    
    def UpdateNetSedimentSourcesAndSinks(self):
        """
        Method that updates fluxes in all sediment storage reservoirs.
        
        Floodplain deposition is handled as a sink here.
        """
        
        # *********************************************************************
        # This code is required in order to make the transition from partly 
        # alluvial to fully alluvial occur without unnecessary fining of the
        # bed, which can increase rates in the partly alluvial sections 
        # unreasonably high.  To handle this, the total divergence in flux is 
        # first computed, summed across all grainsizes.  At this stage, it is 
        # assumed that no fine material bypasses the partly alluvial reach, so 
        # the only outflux to downstream is the transport computed for the 
        # alluvial section (this may underrepresent actual transport).  Then, 
        # after total flux is computed the size distribution of the material 
        # stored in the active layer of the partly alluvial channel is assumed
        # be a weighted fraction of the size distribution of material already 
        # in the active layer and the load.  It must be mostly the same as the 
        # activelayer itself in order to prevent unnreasonably fining.
        # A better way to track this might be to develop functions that 
        # identify which grain sizes are likely to bypass the partly alluvial 
        # reach.  Presumably coarse material is much more likely to be captured
        # by the partly alluvial bed, but it is not clear how to determine how
        # much of each size bypasses the node while traveling across the 
        # non-alluvial sections.
        try:
            self.FixedElev
        except AttributeError:
            self.FixedElev = False
        if self.FixedElev:
            NetInflux = np.zeros(self.NSizes + 1)
            NetInfluxTotal = 0.
            Outflux = 0.
            for k in range(self.NSizes + 1):
                NetInflux[k] = self.ActiveLayer.ExSed.InMigration[k] \
                    + self.ActiveLayer.ExSed.InWidthChange[k] \
                    - self.ActiveLayer.ExSed.OutMigration[k] \
                    - self.ActiveLayer.ExSed.OutWidthChange[k] \
                    + self.ActiveLayer.SourceLatSed[k] \
                    - self.ActiveLayer.SinkLatSed[k] \
                    + self.Load.QsAvkFeed[k]
                NetInfluxTotal += NetInflux[k]
                Outflux += deepcopy(self.Load.QsAvkLoad[k])
            Outflux += self.Load.QsAvkFeed[0] # Katie add mud to output
        # **********************End computation of local net flux *************

        # Net source/sink fluxes to all reservoirs except water column are 
        # computed based on averages across the duration curve
        for k in range(1, self.NSizes + 1):
            self.ActiveLayer.SourceLatSed[k] = self.SLatSourceAv[k]
            # net shaving included in the lateral exchange fluxes
            self.ActiveLayer.SinkLatSed[k] = self.Dfav[k] * self.dxc 
                # Additional sinks such as dredging could be added here.
        
        for k in range(self.NSizes + 1):
            self.ActiveLayer.SourceFeedSed[k] = self.Load.QsAvkFeed[k]
            # **********MASS CONSERVATION IF BED ELEVATION IS FIXED************
        
            #if self.FixedElev:
                # If the total influx in size k is greater than the transport 
                # capacity assuming fully alluvial conditions, this is where we 
                # have to figure out how to handle the extra sediment entering 
                # a node in the partly alluvial case.
                # If there is more entering the node than the transport 
                # capacity out of the node (which depends on the fraction of 
                # the bed that is covered) then some should be passed out and 
                # some should be stored.  Presumably there should not be more 
                # passed out than the total capacity for transport if the node 
                # were fully alluvial. For now, the assumption is that the 
                # material stored has a size distribution similar to the 
                # existing active layer.
            #    if NetInfluxTotal > Outflux:
                    # if the system is gaining sediment, make sure the sediment
                    # stored is similar to the size distribution of the active 
                    # layer
                    # The sediment passed downstream in the load is then equal 
                    # to the total influx in size k minus the sediment in size
                    # k that is stored.
                
                    # Katie modify: make it so that only a portion of the outflux
                    # leaves the system if outfluxes are greater than influxes
                    # and prevents the load from being completely zero (which
                    # was what the old method did).                    
            #        NormStay = (NetInfluxTotal - Outflux) * \
            #            (self.ActiveLayer.GSD.F[k] * \
            #            self.AlphaPartlyAlluvial + \
            #            self.Load.GSDBedloadAv.F[k] * \
            #            (1. - self.AlphaPartlyAlluvial))
                        
            #        self.ActiveLayer.SinkLoadSed[k] = NetInflux[k] - NormStay
                    
            #        x = 'fller'
            #        if self.ActiveLayer.SinkLoadSed[k] < 0.:
            #            self.ActiveLayer.SinkLoadSed[k] = 0. # Katie comment out
                        #self.ActiveLayer.SinkLoadSed[k] = .5*NetInflux[k] # Katie add
                        #print('negative size specific load in partly alluvial \
                            #computation in UpdateSedimentSourcesAndSinks of Node \
                            #was converted to zero')
            #        self.Load.QsAvkLoad[k] = self.ActiveLayer.SinkLoadSed[k] # Katie add to maintain mass conservation
                    #self.Load.ExSed.InVerticalChange[0] = self.ActiveLayer.SinkLoadSed[0]-self.ActiveLayer.SourceFeedSed[0] # Katie add                    
                    #self.ActiveLayer.SinkLoadSed[0] = deepcopy(self.Load.QsAvkFeed[0]) # Katie add
                                        
            #    else:
                    # If the system is sediment starved, assume the only outflux
                    # is that computed from the fraction of the bed that is 
                    # alluvial.
                    
#                    if self.FractionAlluvial < .00001: # Katie add--so that the active layer doesn't run out of sediment
#                        self.ActiveLayer.SinkLoadSed[k] = self.Load.QsAvkFeed[k]                    
#                    else:
            #        self.ActiveLayer.SinkLoadSed[k] = self.Load.QsAvkLoad[k]                    
            #        self.ActiveLayer.SinkLoadSed[0] = deepcopy(self.Load.QsAvkFeed[0]) # Katie add
                    # Katie:  this is where you would add the code to allow mud to 
                    # be entrained--so that the partly alluvial GSD doesn't get
                    # overtaken by mud.
                    
                # Katie try to make it so that partly-alluvial node isn't dominated by mud
            #        self.ActiveLayer.SinkLoadSed[0] = self.ActiveLayer.SourceFeedSed[0] - self.ActiveLayer.GSD.F[0]*(NetInfluxTotal-Outflux)
                    #self.Load.QsAvkLoad[0] = self.ActiveLayer.SinkLoadSed[0]
                    #self.Load.ExSed.InVerticalChange[0] = -self.ActiveLayer.GSD.F[0]*(NetInfluxTotal-sum(self.ActiveLayer.SinkLoadSed))
                    
            #else:
            self.ActiveLayer.SinkLoadSed[k] = deepcopy(self.Load.QsAvkLoad[k])
            self.ActiveLayer.SinkLoadSed[0] = deepcopy(self.Load.QsAvkFeed[0])#Katie add
            #This may not matter if we don't track mud in the active layer.
            #Mud (so size 0) should be accounted for in the water column mass conservation.

        for k in range(self.NSizes + 1):
            self.Floodplain.SourceLatSed[k] = self.Dfav[k] * self.dxc
    
        # Substrate does not have net source or sink terms.  However, 
        # size-specific source/sink terms could be added to any layer to 
        # represent weathering and/or abrasion

        # Net source/sink fluxes to/from water column for mud.  Assumes sand is
        # conserved entirely in active layer.
        # This means sand conservation is NOT applied to each bin of FDC and is
        # just based on annual averages.  Sand is always at capacity based on 
        # sand content in active layer since more sand can always be entrained 
        # as long as sand is present.

        for j in range(self.NFlows):
            self.Load.LatSed[j, 0] = self.SLatSourcejk[j, 0]
            self.Load.SinkSed[j, 0] = self.Dfjk[j, 0] * self.dxc
    
    def UpdateNetTracerSourceAndSinkConcentrations(self):
        """
        Method that updates tracer concentrations in sediment moving into/out of 
        sediment storage reserviors for the node. Also computes the 
        duration-averaged tracer concentration in the load.
        
        This method is necessary because the mass conservation at each node
        operates as a method of clsReservoir and thus does not have access to 
        the node's lateral sources and sinks.
        """
        for L in range(self.NTracers):
            # Set up tracer concentrations in source and sinks for washload 
            # sediment for each bin of duration curve.
            # Duration averaging across load should work for overbank, but
            # water column tracers are a problem.  Needs to be addressed in
            # clsLoad, where the mass conservation for tracers is performed.
            self.Floodplain.SourceLatTracer[0, L] = \
                self.Load.TDepositionAvMudFeed[L]
            self.ActiveLayer.SourceLatTracer[0, L] = 0.
            self.ActiveLayer.SinkLatTracer[0, L] = 0.
        
            # Set up tracer concentrations in sources and sinks for bed 
            # material sediment.
            for k in range(1, self.NSizes + 1):
                self.Floodplain.SourceLatTracer[k, L] = \
                    self.ActiveLayer.T[k, L]
                self.ActiveLayer.SourceFeedTracer[k, L] = \
                    self.Load.TBedFeedk[k, L]
                self.ActiveLayer.SinkLoadTracer[k, L] = \
                    self.ActiveLayer.T[k, L]
                self.ActiveLayer.SourceLatTracer[k, L] = \
                    self.TLatSourceAv[k, L]
    
    def FloodplainDeposition(self, C, Qwf, Fl, Bf):
        """
        Computes overbank deposition rate.
        
        This function computes overbank deposition rates (m/s), outgoing 
        suspended mud flux (m3/s), and fraction mud in point bars (no units) 
        for a single set of steady discharge and sediment supply terms at a 
        single node.  All terms are assumed constant over the interval of the 
        computation. If the computed outgoing suspended sediment flux is less 
        than zero, the deposition rate is reduced so that the outgoing sediment 
        flux is zero and mass is conserved. See Lauer and Parker, 2008, Water
        Resources Research.
        
        Parameters
        ----------
        C : float
            Average suspended sediment concentration above 
            floodplain level.
        Qwf : float
            Water discharge across floodplain (m3/s)
        Fl : float
            Floodplain number (can be size specific)
        Bf : float 
            Floodplain width (m)
    
        Returns
        -------
        float
            Volumetric overbank deposition per unit channel length (m2/s).
        
        Notes
        -----
        Note that the overbank deposition rate is a volume rate of sediment 
        particles deposited on floodplain per unit channel length, so has 
        units of L2/T.  It is not the vertical rate of change floodplain
        elevation. Also note that it does not include a porosity term. 
                
        """
        return Fl * C * Qwf / Bf

    
    def RouseFractionSuspendedLoadAboveFloodplainLevel(self, SettlingVel, \
        ustar, Hc, Lf, Intervals):
        """
        Computes fraction of load for particles of a given settling velocity that turbulence
        keeps suspendeed above the floodplain elevation. 
        
        Computed using Rouse suspended sediment profile.  
        
        Note
        ----
        Units of settling velocity and ustar must be consistent and
        units of Hc and Hf must be consistent.
        
        Parameters
        ----------
        SusSedLoad : float 
            Volumetric suspended sediment load in size class of interest.
        Qc : float 
            Water discharge (m3/s) in channel.
        SettlingVel : float
            Settling velocity for suspended sediment (same units as Ustar)
        Ustar : float 
            Shear velocity in channel (same units as SettlingVel)
        Hc : float
            Flow depth in channel
        Lf : float
            Height of floodplain with respect to channel bed (m).
        Intervals : float 
            Number of intervals over which integral is numerically integrated.
        
        Returns
        -------
        float
            Fraction of load moving above floodplain level.
        """
    
        # Test whether flow is above floodplain and return if flow is below
        # floodplain level

        if Lf / Hc > 1.:
            return 0.
        else:
            # Test whether floodplain elevation is above near bed elevation 
            # of zeta = 0.05.  If not, assume average concentration above
            # floodplain level is just total suspended load divided by channel 
            # discharge
            if Lf / Hc < 0.05:
                print('floodplain elevation is below near bed elevation so Rouse Profile can not be computed.  Concentration set to load/discharge.')

            # Integrate for lower part of water column up to floodplain level.
            # Assume near bed concentration is evaluated at zeta = 0.05 and
            # that no suspended transport occurs below this level.
            dzita = (Lf/Hc - 0.05) / Intervals
            sum1 = 0. # Katie change variable name from sum to sum1 to avoid using the name of a built-in function.
            for n in range(Intervals):
                ZitaDummy = 0.05 + n * dzita - dzita / 2.
                #print dzita                 

                sum1 = sum1 + dzita * self.RouseIntegrand(ZitaDummy, 0.05, \
                    SettlingVel, ustar)
            IntegralBelowFloodplainLevel = sum1
            
            # Integrate for upper part of water column above floodplain level.
            dzita = (1- Lf/Hc) / Intervals
            sum1 = 0.
            ZitaDummy = 0.
            for n in range(Intervals):
                ZitaDummy = Lf/Hc + n * dzita - dzita / 2.
                sum1 = sum1 + dzita * self.RouseIntegrand(ZitaDummy, 0.05, SettlingVel, ustar)

            IntegralAboveFloodplainLevel = sum1
            
            return 0.95 / (1. - (Lf/Hc)) * IntegralAboveFloodplainLevel / \
                (IntegralAboveFloodplainLevel + IntegralBelowFloodplainLevel)

    def RouseIntegrand(self, Zita, Zitab, Vs, ustar):
        """
        Computes integrand for use in Rouse integration.
        
        CONSIDER MAKING A PRIVATE FUNCTION IN RouseFractionSuspendedLoadAboveFloodplainLevel.
        
        Parameters
        ----------
            Zita : float
            ZItab : float
            Vs : float
            ustar : float
            
        Returns
            float
        """
        #print 'Zita ' + str(Zita)
        #print 'Zitab ' + str(Zitab)
        #print 'Vs ' + str(Vs)
        #print 'ustar ' + str(ustar)
        #print ((1. - Zita) * Zitab / ((1. - Zitab) * Zita)) ** (Vs / (0.4 * \
        #    ustar))
        return ((1. - Zita) * Zitab / ((1. - Zitab) * Zita)) ** (Vs / (0.4 * \
            ustar))
    
    def UpdateDMudj(self):
        """
        Method for updating mud deposition rate for each bin and
        computing flow-duration averaged value of overbank mud deposition rate (k = 0).
        """
        self.Dfav[0] = 0.
        for j in range(self.DC.NFlows()):
            C = self.Load.Qsjkfeed[j, 0] / self.DC.Qwc[j] # Could decided to 
                # use Qw rather than Qwc here to be consistent with rating 
                # curve

            self.Dfjk[j, 0] = self.FloodplainDeposition(C, self.DC.Qwf[j], \
                self.Flmud, self.Bf)
            #print self.Dfjk[j,0]
            #print 'Fp Q' + str(self.DC.Qwf[j])
            #print 'total Q' + str(self.DC.Qw[j])

            self.Dfav[0] += self.Dfjk[j, 0] * self.DC.p[j]

    def UpdateDSandj(self):
        """
        Method for updating sand deposition rate for each bin and
        computing flow-duration averaged value of overbank sand
        deposition rate (k = represntative of all sand sizes).
        """
        g = 9.81
        #print self.Floodplain.L-self.ActiveLayer.L
        #print self.DC.Hc
        for k in range(1, self.ActiveLayer.GSD.NBedSizes + 1):
            # Only used for sand size classes, but could be used for gravel too
            # if shear velocity was high enough. 
            # Actually, this routing could work for mud too as long as a 
            # reasonable settling velocity was available.
            # It is written as a separate subroutine since the floodplain 
            # numbers for sand and mud can be different.
            self.Dfav[k] = 0. # Katie move here--so coarse fractions have value of 0.
            if self.ActiveLayer.GSD.D[k] < 32. and self.ActiveLayer.GSD.D[k] > .063:
                #self.Dfav[k] = 0. # Katie comment out and move above loop
                for j in range(self.NFlows):
                    try:
                        self.Flbed
                    except AttributeError:
                        self.Flbed = 0.
                    if self.Flbed == 0.: # Do not perform rouse integration if
                    # we already know there is no sand deposition on floodplain
                        self.Dfjk[j, k] = 0.
                    else: # Perform rouse integration and compute overbank 
                    # deposition
                        ustar = (g * self.DC.Hc[j] * self.Slope) ** 0.5
                        # The ratio of suspended sediment above floodplain to 
                        # below is computed using only 20 intervals to save 
                        # time. This can be changed.
                        # This ratio could also be made a property of the node
                        # and only be recomputed periodically.
                        CRatio = self.RouseFractionSuspendedLoadAboveFloodplainLevel(\
                            self.ActiveLayer.GSD.Vs[k], ustar, self.DC.Hc[j], \
                            self.Floodplain.L - self.ActiveLayer.L, 20)
                        # Suspended sediment concentration above floodplain is computed assuming a constant fraction me.FSandSusp of the load is suspended.
                        # Could decided to use Qw rather than Qwc here to be consistent with rating curve
                        C = self.Load.Qsjkfeed[j, k] * self.FSandSusp * CRatio / self.DC.Qwc[j]                     
                        self.Dfjk[j, k] = self.FloodplainDeposition(C, self.DC.Qwf[j], self.Flbed, self.Bf)
                        self.Dfav[k] = self.Dfav[k] + self.Dfjk[j, k] * self.DC.p[j]
                        #print 'size' + str(self.ActiveLayer.GSD.D[k])                       
                        #print self.Dfav[k]
        #print 'sand' +  str(self.Dfav)
    def WidthChange(self, Manning, TrinityFit, CalibrationFactor, SG, rho_w, g, W, alphaF, alphatau, BcMin, dt, ControlGSD, MobilityThreshold):
        """
        Katie added this method.  It is an optional add-on to replace the constant
        migration rate (Node.cbank).  Bank erosion and floodplain sequestration are
        treated as separate, independent processes.
        
        Parameters
        ----------
        Manning : bool
            APPEARS TO BE UNUSED
        TrinityFit : bool
            Flag that determines if Gaeuman fit to Wilcock & Crowe is used.  Regular
            Wilcock and Crowe used if false.
        CalibrationFactor : float
            Sediment transport calibration factor. Used to adjust reference Shields
            stress in Wilcock Crowe type sediment transport computation.
        SG : float
            Specific gravity of sediment.
        rho_w : float
            Density of water
        g : float
            APPEARS TO BE UNUSED
        W : float
            Narrowing constant--used to calibrate narrowing function.  
            As a starting guide, estimate the percentage of bar that is 
            vegetated annually and double it.    
        alphaF : float
            NEEDS TO BE DOCUMENTED
        alphatau : float
            NEEDS TO BE DOCUMENTED
        BcMin : float
            Minimum channel width (m)       
        dt : float
            Timestep
        ControlGSD : ???
            APPEARS TO BE UNUSED
        MobilityThreshold : float   
            Parameter representing mobility of bank toe material.
        """
        index = -3
        self.cbank = 0. # Turns off constant migration rate so that ExchangeType values can be filled with this function        
        #self.NarrowRate = 0.
        self.WidenRate = 0.        
        
        # Bank erosion function.  This is based on the work from Sarah Davidson's
        # PhD thesis but modified to fulfill our mass conservation requirements
        # and smaller timestep.  The channel is widened incrementally until the
        # D84 is no longer fully mobile (full mobility is defined as twice the 
        # critical shear stress for the D84--see Sarah's thesis for details.)
        
        # In order to calculate the shear stress on the widened bed, we need to
        # recalculate hydraulics.  For simplicity, normal flow is assumed within
        # each node.
    
        oldBc = deepcopy(self.Bc) # Must save to use the old channel width for the conservation of mass later on.                
        
        TotalWidening = np.zeros(self.NFlows)
        #TotalNarrowing = []        
        
        for j in range(self.NFlows):
            self.Bc = oldBc
            #pi = self.Load.Qsjk[j,-1]/sum(self.Load.Qsjk[j,1:])
            pi = sum(self.Load.Qsjk[j,index:])/self.Bc
            fi = sum(self.ActiveLayer.GSD.F[index:])/sum(self.ActiveLayer.GSD.F[1:])
            
#            ControlQsk = self.Load.WilcockCroweLoadBySize(self.DC.Uc[j], self.DC.Sf[j], ControlGSD, rho_w, SG, self.Bc, TrinityFit, \
#                    CalibrationFactor)
                    
            self.CobbleMobility = (pi/fi)#*(sum(self.Load.QsAvkLoad[1:])/sum(ControlQsk))   # Normalize by ratio above 'reference load'--bar push mechanism     
            
            if self.CobbleMobility > MobilityThreshold:
        
                # Create composite GSD from active layer and floodplain
                HybridGSD = deepcopy(self.ActiveLayer.GSD)
                for k in range(self.NSizes + 1):
                    HybridGSD.F[k] = self.ActiveLayer.GSD.F[k]*alphaF + self.Floodplain.GSD.F[k]*(1-alphaF)
                HybridGSD.UpdateStatistics()
                
                # Calculate total transport rate of hybrid mixture.
#                FPcap = (self.Load.Qsjk[j,index]/self.ActiveLayer.GSD.F[index]/self.Bc)*self.Floodplain.GSD.F[index] # capacy of floodplain coarse grains
                qsk = self.Load.WilcockCroweLoadBySize(self.DC.Uc[j], self.DC.Sf[j], HybridGSD, rho_w, SG, self.Bc, TrinityFit, \
                    CalibrationFactor, IndD = index)
                qs = (sum(qsk[index:])/self.Bc)*((sum(self.Floodplain.GSD.F[index:])*(1-alphaF))/sum(HybridGSD.F[index:]))
                #qs = (qsk/self.Bc)*(self.Floodplain.GSD.D90*(1-alpha))/HybridGSD.D90
                
                TotalWidening[j] = (qs/sum((self.Floodplain.L*self.Floodplain.GSD.F[index:])))*self.DC.p[j]
                #TotalWidening[j] = qs/self.Floodplain.L*.1
                #TotalWidening[j] = (qs/self.Floodplain.L)*(1-alpha)
#                TotalWidening[j] = (FPcap/(self.Hpb*self.Floodplain.GSD.F[index]))*(1-alpha)
        
        if oldBc + sum(TotalWidening)*dt >= self.ValleyMargin: # Makes sure that width does not exceed valley width
            TotalWidening = [(self.ValleyMargin-oldBc)/dt] # Adjusts Bcrate to max widening

        self.WidenRate = sum(TotalWidening)

        self.Bc = oldBc
        
        self.Narrowing(BcMin, W, alphatau)

    
    def Narrowing(self, BcMin, W, alphatau):
        """
        Katie added this method.  It appears to complete the width change 
        computation performed by WidthChange method. CONSIDER INCLUDING 
        THIS AS PART OF THE WIDTHCHANGE METHOD.
        
        Parameters
        ----------
        BcMin : float
            Minimum channel width (m)  
        W : float
            NEEDS TO BE DOCUMENTED   
        alphatau : float
            NEEDS TO BE DOCUMENTED
        """
        
        self.NarrowRate = 0.
        TotalNarrowing = []
        oldHc = self.DC.Hc
        oldBc = self.Bc

        for j in range(self.NFlows):
            widthchangeN = 0.
            TauPrime = self.Load.findWC_TauPrime(9.81, self.DC.Sf[j], self.DC.Uc[j], self.ActiveLayer.GSD.D65)

            if TauPrime < alphatau and self.DC.Qwf[j] == 0.: # Don't want narrowing in flooded zones (i.e. reservoirs)            
#            if oldHc[j] < .25*self.Floodplain.L and oldBc + sum(TotalNarrowing) > BcMin: #.465 # Make narrowing a function of channel depth, not discharge.  Narrows when flow depth is quarter bankfull.
    
    #            Fs = 0. # Adjust narrowing rate by sand content--not currently in use
    #            i = 0            
    #            while self.ActiveLayer.GSD.D_upper[i] <= 2.: # Katie change from .002 (which is the lower boundary)
    #                Fs = self.ActiveLayer.GSD.FfinerThanD_upper[i] # Katie change from parentheses to brackets
    #                i = i + 1		
                widthchangeN = -((oldBc+sum(TotalNarrowing))-BcMin)*W/(365.25*24*60*60) #Katie multiply by percent sand # Change per second--simplest solution; only one calibration term which is percent widening per year.
            #else:
                #print('widthchangeN not set, Tauprime = %s and alphatau = %s' % (TauPrime,alphatau)) 
                #widthchangeN = -((oldBc+sum(TotalNarrowing))-BcMin)*W/(365.25*24*60*60)
            TotalNarrowing.append(widthchangeN*self.DC.p[j])                
        self.NarrowRate = sum(TotalNarrowing)
        #print self.NarrowRate*576.*150
        #self.Bcrate = sum(TotalWidening) + sum(TotalNarrowing)
        #self.Bc = oldBc
    
    def Avulsion(self, threshold, spacing, exchange):       
        """
        Katie added this method. Simulates avulsion on a node-by-node basis.
        
        When the channel elevation approaches that of the bank (how close is a
        user-defined threshold), the uppermost substrate layer will be incorporated
        back into the floodplain. The active layer will also join the Floodplain
        reservoir.  The bed elevation (etabav) will be adjusted so that it is
        lowered by an increment of 1 substrate spacing. If this occurs repeatedly
        in a given node, it should simulate valley filling, like that of a delta.
        
        Parameters
        ----------
        threshold : float
            difference between the bank (Floodplain.L-ActiveLayer.L)
            and the channel elevation after which avulsion will occur.
        spacing : ???
            NEEDS TO BE DOCUMENTED
        exchange : ???
            NEEDS TO BE DOCUMENTED
        """
        
        if self.Floodplain.L-self.ActiveLayer.L < threshold:
            print('Avulsion!')
            print(self.Floodplain.L)
            # This part of the code extracts a slice of substrate that is the 
            # thickness of the user-specified original substrate spacing.  This
            # is so that the depth of the new channel is consistent and does not
            # depend on the thickness of the upper substrate layer.  It is assumed
            # that no extra floodplain material enters the load when the avulsion
            # occurs.  This is not perfect because the floodplain height increases
            # more than the bed elevation decreases.

            # If topmost substrate layer is thicker than the spacing, split it.
            if self.Substrate[-1].C.L > spacing:
                self.SplitOrCombineSubstrate(0, 0, spacing)
                newLayer = deepcopy(self.Substrate[-2]) # this is the spacing-thick substrate to keep
                del self.Substrate[-2]
                self.Substrate.append(newLayer) # Now new layer is at the top.
            # If topmost substrate layer is thinner than the spacing, combine
            # substrate layers and extract the spacing-thick layer
            elif self.Substrate[-1].C.L < spacing:
                self.SplitOrCombineSubstrate(spacing, spacing, spacing)
                self.SplitOrCombineSubstrate(0, 0, spacing)
                newLayer = deepcopy(self.Substrate[-2]) # this is the spacing-thick substrate to keep
                del self.Substrate[-2]
                self.Substrate.append(newLayer) # Now new layer is at the top.
                            
            # The new floodplain volume includes the uppermost substrate layer
            NewFpVolume = self.Floodplain.Volume + self.Substrate[-1].F.Volume\
                + self.Substrate[-1].C.Volume

            # New floodplain grainsize distribution incorporates former floodplain,
            # part of the old active layer, and both substrate channel and floodplain reservoirs.
            # An part chunk is taken out to form part of the new active layer.
            NewFpF = np.zeros(self.NSizes+1)
            NewFpT = np.zeros((self.NSizes+1, self.NTracers))
            
            for k in range(self.NSizes+1):
                NewFpF[k] = (self.Floodplain.GSD.F[k]*self.Floodplain.Volume + \
                    self.Substrate[-1].C.GSD.F[k]*self.Substrate[-1].C.Volume +\
                    self.Substrate[-1].F.GSD.F[k]*self.Substrate[-1].F.Volume +\
                    exchange*self.ActiveLayer.GSD.F[k]*self.Bc*self.dxc*self.ActiveLayer.L -\
                    exchange*self.Floodplain.GSD.F[k]*self.Bc*self.dxc*self.ActiveLayer.L)/\
                    NewFpVolume
                    
                # Do the same thing for tracers                
                for m in range(self.NTracers):
                    NewFpT[k, m] = (self.Floodplain.T[k, m]*self.Floodplain.Volume + \
                    self.Substrate[-1].C.T[k, m]*self.Substrate[-1].C.Volume +\
                    self.Substrate[-1].F.T[k, m]*self.Substrate[-1].F.Volume +\
                    exchange*self.ActiveLayer.T[k, m]*self.Bc*self.dxc*self.ActiveLayer.L -\
                    exchange*self.Floodplain.T[k, m]*self.Bc*self.dxc*self.ActiveLayer.L)/\
                    NewFpVolume
                    
            # The active layer is given the GSD and tracer content of part of the (old) floodplain
            # Might think about changing this later--would be technically more
            # correct to have GSD a combination of Floodplain and Substrate,
            # but the difference is likely to be miniscule.
            for k in range(self.NSizes+1):
                self.ActiveLayer.GSD.F[k] = self.Floodplain.GSD.F[k]*exchange+self.ActiveLayer.GSD.F[k]*(1-exchange)
                for m in range(self.NTracers):
                    self.ActiveLayer.T[k,m] = self.Floodplain.T[k,m]*exchange + self.ActiveLayer.T[k,m]*(1-exchange)
            self.ActiveLayer.GSD.UpdateStatistics()
            
            # Update floodplain reservoir 
            self.Floodplain.Volume = NewFpVolume
            self.Floodplain.GSD.F = NewFpF
            self.Floodplain.T = NewFpT
            self.Floodplain.GSD.UpdateStatistics()
            
            # Update substrate reservoirs--mix
            for m in range(self.NLayers()):
                for k in range(self.NSizes+1):
                    self.Substrate[m].C.GSD.F[k] = self.dxc*self.Substrate[m].C.L*\
                        (self.Substrate[m].C.GSD.F[k]*self.Bc\
                        + self.Substrate[m].F.GSD.F[k]*self.Bc*exchange - \
                        self.Substrate[m].C.GSD.F[k]*self.Bc*exchange)/self.Substrate[m].C.Volume
                        
                    self.Substrate[m].F.GSD.F[k] = self.dxc*self.Substrate[m].F.L*\
                        (self.Substrate[m].F.GSD.F[k]*self.Bc\
                        + self.Substrate[m].C.GSD.F[k]*self.Bc*exchange - \
                        self.Substrate[m].F.GSD.F[k]*self.Bc*exchange)/self.Substrate[m].F.Volume
                        
                self.Substrate[m].C.GSD.UpdateStatistics()
                self.Substrate[m].F.GSD.UpdateStatistics()
            
            # Lower channel elevation to base of former substrate level 
            self.etabav = self.etabav-self.Substrate[-1].C.L
            self.Floodplain.L = self.Floodplain.Volume / (self.Bf * self.dxc / \
                self.ChSin)
            
            # Get rid of old uppermost substrate
            del self.Substrate[-1]
            self.CumulativeTotalAvulsionWidth = self.CumulativeTotalAvulsionWidth + self.Bc*exchange

    def ExnerBed(self):
        """
        Method for computing total net influx of sediment to active layer
        from all sources and then updating bed elevation. Calculation is summed across all sizes.
        """
    
        # Compute total net lateral influx to active layer due to migration,
        # width change, and lateral sources, summed across all bed material
        # sizes

        NetQsIn = np.sum(self.ActiveLayer.ExSed.InMigration[1:self.NSizes + 1]\
            - self.ActiveLayer.ExSed.OutMigration[1:self.NSizes + 1] + \
            self.ActiveLayer.ExSed.InWidthChange[1:self.NSizes + 1] - \
            self.ActiveLayer.ExSed.OutWidthChange[1:self.NSizes + 1] +\
            self.ActiveLayer.SourceLatSed[1:self.NSizes + 1] - \
            self.ActiveLayer.SinkLatSed[1:self.NSizes + 1] + \
            self.ActiveLayer.SourceFeedSed[1:self.NSizes + 1] - \
            self.ActiveLayer.SinkLoadSed[1:self.NSizes + 1])

        #NetQsIn = np.sum(self.ActiveLayer.ExSed.InMigration[1:self.NSizes + 1]\
        #    - self.ActiveLayer.ExSed.OutMigration[1:self.NSizes + 1] + \
        #    self.ActiveLayer.ExSed.InWidthChange[1:self.NSizes + 1] + \
        #    self.ActiveLayer.SourceLatSed[1:self.NSizes + 1] - \
        #    self.ActiveLayer.SinkLatSed[1:self.NSizes + 1] + \
        #    self.ActiveLayer.SourceFeedSed[1:self.NSizes + 1] - \
        #    self.ActiveLayer.SinkLoadSed[1:self.NSizes + 1])
             
        # This is where a formulation accounts for the possibility that the bed
        # is partially alluvial:
        if self.FixedElev:
            self.DeltaEtaB = 0.
        else:
            self.DeltaEtaB = NetQsIn / self.dxc / self.Bc / (1. - self.lambdap)
            
        # Katie add variables for cumulative channel volume change--as an alternative
        # to CumulativeBedChange when the Avulsion fuction is working.
        self.CumulativeChannelVolumeChange = self.CumulativeChannelVolumeChange + NetQsIn


    def UpdateVolumeSizeFractionsandTracersInAllReservoirs(self, dt, \
        TracerProperties):
        """
        Method for mixing sediment and tracer concentrations into reservoirs.
        
        Parameters
        ----------
            dt : float
                Timestep.
            TracerProperties : array_like(:obj:clsTracerProperties, length = NTracers)
                Production and decay properties for each tracer type.
        """
        self.ActiveLayer.UpdateFandT(dt, self.lambdap, TracerProperties)
        self.Floodplain.UpdateFandT(dt, self.lambdap, TracerProperties)
        #for m in range(self.NLayers()): # Katie take out of loop
        self.Substrate[-1].C.UpdateFandT(dt, self.lambdap, TracerProperties)
        self.Substrate[-1].F.UpdateFandT(dt, self.lambdap, TracerProperties)

        # Katie add--need to artificially add point bar volume to active layer
        # if channel is narrowing to maintain correct fraction alluvial.
        self.ActiveLayer.Volume += -self.NarrowRate*self.Hpb*self.dxc*dt

        # Should include tracer production and decay either here or in a
        # separate sub
    
    def UpdateGeometricParameters(self, dt):
        """
        Method for updating the geometric parameters of the node.
        
        Parameters
        ----------
            dt : float
                Timestep.
                
        Note
        ----
        For simulating a node with partly alluvial cover on bedrock, the 
        subroutine may need some adjustment.  See comments in the code
        for ideas.        
        """

        self.etabav = self.etabav + self.DeltaEtaB * dt
            
        #self.etabav = round(self.etabav,4) # Katie add
        #self.Bc = self.Bc + self.Bcrate * dt # Katie:  take this and move it to WidthChange function to put it in correct place before Exner equation is done.        
        #self.Bf = self.Bf - self.Bcrate * self.ChSin * dt # Katie:  take this and move it to WidthChange function to put it in correct place before Exner equation is done.
        
        #self.Bc = self.Bc + self.Bcrate * dt
        self.Bc = self.Bc + (self.NarrowRate + self.WidenRate) * dt # Katie add
        #self.Bf = self.Bf - self.Bcrate * self.ChSin * dt
        self.Bf = self.Bf - (self.NarrowRate + self.WidenRate) * self.ChSin * dt # Katie add
        
                    
        # Katie add to keep track of how much narrowing and bank erosion occurs
        #if self.NarrowRate != 0.:
        self.CumulativeNarrowing = self.CumulativeNarrowing + self.NarrowRate*dt

        self.CumulativeWidening = self.CumulativeWidening + self.WidenRate*dt
        self.CumulativeWideningVolume = self.CumulativeWideningVolume + (self.WidenRate*dt*self.Floodplain.L*self.dxc)
        self.CumulativeWideningBedMaterialVolume = self.CumulativeWideningBedMaterialVolume + (self.WidenRate*dt*self.Floodplain.L*self.dxc*(1-self.Floodplain.GSD.F[0]))
            
       # Katie add--keep track of total feed
        self.CumulativeTotFeed = self.CumulativeTotFeed + (np.sum(self.Load.QsAvkFeed)*dt)
        oldFeed = self.CumulativeTotBedMaterialFeed
        self.CumulativeTotBedMaterialFeed = oldFeed + (np.sum(self.Load.QsAvkFeed[1:])*dt)
        
        self.Floodplain.L = self.Floodplain.Volume / (self.Bf * self.dxc / \
            self.ChSin)
        self.H = self.Floodplain.L

        self.Substrate[-1].C.L = self.Substrate[-1].C.L + self.DeltaEtaB * dt
        self.Substrate[-1].F.L = self.Substrate[-1].F.L + self.DeltaEtaB * dt
        #self.Substrate[-1].C.Volume = self.Substrate[-1].C.L * self.Bc * \
        #    self.dxc # Katie comment out--unnecessary code as volume is updated in UpdateFandT
        #self.Substrate[-1].F.Volume = self.Substrate[-1].F.L * self.Bf * \
        #    self.dxc/self.ChSin # Katie add ChSin
        if abs(self.Substrate[-1].F.L - self.Substrate[-1].C.L) > 0.001:
            print('Problem with substrate thickness')

        TotalSubstrateVolume = 0.
        TotalSubstrateFDeltaS = 0.
        TotalSubstrateCDeltaS = 0.
        TotalCVol = 0

        self.CumuTotVolume = self.CumuTotVolume + self.Floodplain.DeltaS + self.ActiveLayer.DeltaS
        for SubClass in self.Substrate:
            TotalCVol = TotalCVol + SubClass.C.Volume
            TotalSubstrateVolume = TotalSubstrateVolume + SubClass.C.Volume + SubClass.F.Volume
            self.CumuTotVolume = self.CumuTotVolume + SubClass.C.DeltaS + SubClass.F.DeltaS          
            self.CumuTotalSubstrateFDeltaS = self.CumuTotalSubstrateFDeltaS + SubClass.F.DeltaS
            self.CumuTotalSubstrateCDeltaS = self.CumuTotalSubstrateCDeltaS = SubClass.C.DeltaS
        self.TotalSubstrateVolume = TotalSubstrateVolume
        self.TotalCVol = TotalCVol

        #self.ActiveLayer.Volume = self.ActiveLayer.L*self.Bc*self.dxc*self.FractionAlluvial


        # This method should probably be modified for a partly alluvial case.
        # If dt is so long that channel incises through all substrate, then 
        # need to reduce alluvial cover fraction.  Conversely, for partly 
        # alluvial case, if bed cover fraction would increase to greater than
        # 1, then need to allow bed elevation to change and also to compute 
        # appropriate vertical boundary exchange fluxes.  Should be easy to 
        # write if statements that modify boundary flux compuations.
        # Problem is how to handle the last substrate layer when it is almost 
        # gone, especially since all funtions to this function have computed 
        # rates, not total volumes.  Given the rate at which the boundary is 
        # moving down (i.e. the vertical movement of the bed into the 
        # substrate), it should be possible to compute here when during the 
        # timestep the last substrate layer is destroyed.  At this point, the 
        # alluvial cover on the bed is excatly equal to 100%.  Then, it is 
        # straightforward to modify the vertical boundary flux terms.
    
        # The simplest approach would be to compute the length of time 
        # necessary for full removal of the substrate at the computed incision
        # rate.  Once this is known, the geometry can be updated using the 
        # reduced time step length, and the fractional cover on the partly 
        # alluvial bed can then be computed using the net sediment efflux 
        # (already known) during the remainder of the timestep.  To go the 
        # other way, from a partly alluvial case to a fully alluvial one, it 
        # will be nessary to implement this in the load computation, which 
        # should compute the overall outflux from a node within the load as a
        # product of the fraction of alluvial cover and the computed tranport 
        # rate assuming 100% cover.  The, in the Exner equation, the rate of 
        # change of fractional alluvial cover will be tracked.  If this is 
        # increasing at a rate sufficient to fully cover the bed by the end of
        # the timestep, then all computations should be done for a duration 
        # equal to this length.  Then, a new substrate layer would be spawned,
        # with a volume and size distribution equal to the sediment that would
        # have been associated with the widening of the active zone over this
        # time. It is probably OK to ignore vertical exchange fluxes during the
        # short timestep since the transition back to fully alluvial should 
        # occur only infrequently.  However, for full mass conservation, it 
        # would be possible to recompute all fluxes, re-estimate vertical 
        # fluxes, and then complete the timestep.  Doing this wouldn't take 
        # too long.
    
        # An even simpler method that does not perfectly conserve mass is to
        # simply do mass conservation as normal, and to turn a node into a 
        # partly-alluvial node when the incision rate is high enough to work
        # through the remaining substrate in a given timestep.  Bed elevation
        # would be specified as the top of substrate plus the active layer 
        # thickness at the end of such a timestep, resulting in a small 
        # inconsistency in mass conservation for such a step.  From then on,
        # rate of fraction cover change would be computed in Exner rather than
        # bed elevation change.  In the case when fraction cover becomes equal
        # to one, the node would simply be switched back to fully alluvial 
        # (with no associated aggradation until the next timestep).
