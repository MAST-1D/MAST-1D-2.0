
from clsSubstratePairClass import clsSubstratePairClass
from clsReservoir import clsReservoir
from clsDurationCurve import clsDurationCurve
from clsLoad import clsLoad
from copy import deepcopy
import pdb
import numpy as np

class clsNode(object):
    """ 
    Attributes:
        ActiveLayer -- clsReservoir
        Substrate -- [clsSubstratePairClass]
        Floodplain -- clsReservoir
        Load -- clsLoad (This represents sediment moving out of the node)
        Slope -- float (Down channel slope)
        H -- float
        Bf -- float (floodplain width)
        Bc -- float (channel width)
        cbank -- float (Bank migration rate normal to down-channel direction)
        DeltaEtaB -- float (Rate of bed elevation change)
        InitialBedElev -- float (Initial bed elevation)
        Hpb -- float (point bar thickness)
        Bcrate -- float (widening rate)
        etabav -- float (bed elevation)
        ChSin -- float (channel sinusity)
        xc -- float (down channel coordinate)
        x -- float (down valley coordinate)
        Dfjk -- [float] (deposition rate on the floodplain in each bin of the
                flow duration curve j and for each size k--could be a
                property of load)
        SLatSourcejk -- [float] (lateral sources of sediment in m3/s in each 
                        bin of the flow duration curve j in size k)
        SLatSourceAv -- [float] (mean annual lateral source of sediment in 
                        m3/s in size k)
        Dfav -- [float] (duration-averaged deposition rate on the floodplain 
                in size k)
        FixedElev -- bool (flag to determine if the node's elevation is fixed 
        (e.g. bedrock or boundary condition))
        
        ****Note that pTLatSource could be adjusted as a function of time to 
        represent individual erosion events.  This would require changing the 
        duration curve at various times in the computation as well.
        TLatSourceAv -- [float] (Duration-averaged concentration of cosmogenic
                        nuclides in the sediment coming from lateral sources 
                        in size k and tracer l)
        DC -- clsDurationCurve (Duration curve object used to store flow and
              sediment flux in each bin of flow duration curve.)
        Initialized -- bool (flag to determine if arrays have been 
                       redimensioned appropriately)
        lambdap -- float (porosity)
        sigma -- float (subsidence rate)
        dxc -- float (Down channel spacing to next node represented by this 
               node)
        Cfc -- float (Friction coefficient for channel)
        Cff -- float (Friction coefficient for floodplain)
        FSandSusp -- float (Fraction sand load that is suspended)
        Kbar -- float (Parameter controlling fraction washload in point bar 
                deposits)
        AlphaBar -- float (Parameter controlling similarity between bed 
                    material load and bar deposition.)
        AlphaPartlyAlluvial -- float (Parameter controlling similarity between
        bed material load and deposition in the active layer of a partly
        alluvial node)
        nc -- float (Manning's n for channel)
        ncAddons -- float (Form drag addition for manning's n)
        ncMultiplier -- float (Sinuosity multiplier for manning's n)
        nf -- float (Manning's n for floodplain)
        PointBarSubsurfaceGSD -- clsGSD
        Flmud -- float (Floodplain number for mud)
        Flbed -- float (Floodplain number for bed material sizes)
        Data-structure related parameters
        NSizes -- int (number of sizes)
        NTracers -- int (number of tracers)
        NFlows -- int (number of flows in FDC)
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
            Here, Kbar defines the weighting of overall washload vs. bed 
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
                print('Bed material load is zero.  Washload fraction in ' + \
                    'point bar set equal to 1.')
            
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
        Return: float
        """
        # Katie make it so that mud is not in roughness 
        #calculation--for partly alluvial problem
        mud = self.ActiveLayer.GSD.F[0] 
        self.ActiveLayer.GSD.F[0] = 0.
        
        nGrain =  (self.ActiveLayer.GSD.D65) ** (1. / 6.) * 0.0146
        self.ActiveLayer.GSD.F[0] = mud
        #self.ActiveLayer.GSD.UpdateStatistics()
        return nGrain
    
    def nc(self): # Grain roughness
        """
        Return: float
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
        return self.ActiveLayer.Volume / self.ActiveLayer.L / self.Bc / \
            self.dxc
    
    FractionAlluvial = property(getFractionAlluvial)

    def SplitOrCombineSubstrate(self, LMinAfterSplit, LMinBeforeRemove, \
            Spacing):
        """
        Arguments:
            LMinAfterSplit -- float
            LMinBeforeRemove -- float
            Spacing -- float
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
            self.Substrate[-1] = deepcopy(self.Substrate[-2])
            self.Substrate[-2].C.L = Spacing
            self.Substrate[-2].F.L = Spacing
            self.Substrate[-2].C.Volume = Spacing * self.Bc * self.dxc
            self.Substrate[-2].F.Volume = Spacing * self.Bf * self.dxc / \
                self.ChSin
            self.Substrate[-1].C.L = NewL
            self.Substrate[-1].F.L = NewL
            self.Substrate[-1].C.Volume = NewL * self.Bc * self.dxc
            self.Substrate[-1].F.Volume = NewL * self.Bf * self.dxc / self.ChSin # Katie fix from 'Me.ChSin'
            print('SubstrateSplit')
        
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
    
    def __init__(self, NLayers, NTracers, BinBdySizes, NFlows):
        """
        Arguments:
            NLayers -- int
            NTracers -- int
            BinBdySizes -- [float]
            NFlows -- int
        """

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
            self.CumulativeNarrowing = 0. # Katie add
            self.CumulativeWidening = 0. # Katie add
            self.CumulativeTotFeed = 0. # Katie add
            self.CumulativeWideningVolume = 0. # Katie add
            self.CumulativeChannelVolumeChange = 0. # Katie add
            self.PartlyAlluvial = False # Katie add
            self.Canyon = False # Katie add
            self.ControlSubstrate = clsSubstratePairClass() # Katie add to hold original substrate (see Add/Remove Substrate function)

            self.Initialized = True
        else:
            raise RuntimeError('Tried to initiate clsNode twice.')
    
    def EquilibriumMudFloodplainNumber(self, DeltaEta, cbank):
        """
        This function can be used to compute the floodplain number for mud at 
        perfect equilibrium, when the bed is not changing and there is no net 
        divergence in load.  It must be called after the hydraulics have been 
        updated and mud feed has been specified for the node.
        The variable DeltaEta should represent the mud fraction of the 
        thickness of the overbank sediment deposit.
        
        Arguments:
            DeltaEta -- float
            cbank -- float
        
        Return: float
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
        This subroutine computes the normal water depth as a function of 
        channel width, floodplain thickness, floodplain width, friction 
        coefficients in channel and on floodplain, and slope, and partitions
        flow into channel and floodplain zones.

        Can perform calculation either using Manning's equation (for Manning
         = True) or Chezy equation (For Manning = False)
        
        Arguments:
            Qw -- float
            threshold -- float
            Manning -- bool
            
        Return: [float]
        """

        g = 9.81
        Czc = 1 / self.Cfc ** 0.5
        Czf = 1/ self.Cff ** 0.5


        if Manning:
            HcGuess = (self.nc() * Qw / (self.Bc * self.Slope ** 0.5)) ** \
                (3. / 5.)

        else:
            HcGuess = (Qw / (self.Bc * Czc * (g * self.Slope) ** 0.5)) ** \
                (1 / 1.5)
        if HcGuess < (self.Floodplain.L - self.ActiveLayer.L):
            Hc = HcGuess
            Qc = Qw

        else:
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
        Arguments:
            H -- float
            Bc -- float
            Czc -- float
            Sc -- float
            Bf -- float
            Czf -- float
            Tf -- float
            Sf -- float
        
        Return: float
        """
        g = 9.81
        if H > Tf:
            return Bc * Czc * H * (g * H * Sc) ** 0.5 + Bf * Czf * (H - Tf) * \
                (g * (H - Tf) * Sf) ** 0.5
        else:
            return Bc * Czc * H * (g * H * Sc) ** 0.5
    
    def ManningDischarge(self, H, Bc, nc, Sc, Bf, nf, Tf, Sf):
        """
        Arguments:
            H -- float
            Bc -- float
            nc -- float
            Sc -- float
            Bf -- float
            nf -- float
            Tf -- float
            Sf -- float
        
        Return: float
        """

        if H > Tf:
            return 1. / nc * Bc * H * H ** (2. / 3.) * Sc ** (1. / 2.) + 1. \
                / nf * Bf * (H - Tf) * (H - Tf) ** (2. / 3.) * Sf ** (1. / 2.)
        else:
            return 1. / nc * Bc * H * H ** (2. / 3.) * Sc ** (1. / 2.)
    
    def UpdateDepthAndDischargeAtAllFlows(self, Manning):
        """
        Arguments:
            Manning -- bool
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
    
    def UpdateLateralSedFluxes(self):
        """
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
                if self.Bcrate < 0.:
                    self.Substrate[m].C.ExSed.OutWidthChange[k] = \
                        -(1. - self.lambdap) * self.Substrate[m].C.GSD.F[k] * \
                        self.Substrate[m].C.L * self.Bcrate * self.dxc # Katie:  good:  this produces a positive number
                else:
                    self.Substrate[m].C.ExSed.OutWidthChange[k] = 0.
                # For floodplain zone # Katie comment out
                self.Substrate[m].F.ExSed.OutMigration[k] = \
                    (1. - self.lambdap) * self.Substrate[m].F.GSD.F[k] * \
                    self.Substrate[m].F.L * self.cbank * self.dxc # Katie unindent beginning of this line bloc
                
                # The following ensures that if channel is widening, the
                # floodplain zone substrate is exporting sediment to channel
                # zone
                if self.Bcrate > 0.:
                    self.Substrate[m].F.ExSed.OutWidthChange[k] = \
                        (1. - self.lambdap) * self.Substrate[m].C.GSD.F[k] \
                        * self.Substrate[m].F.L * self.Bcrate * self.dxc # Katie get rid of negative sign at the beginning of this code bloc--number needs to be positive for UpdateFandT function.
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
            else:
                self.ActiveLayer.ExSed.OutMigration[k] = \
                    (1. - self.lambdap) * (self.ActiveLayer.GSD.F[k] * \
                    self.ActiveLayer.L + \
                    self.FkPointBarDepositAnnualAverage[k] * self.Hpb) * \
                    self.cbank * self.dxc
        
            Floodplain = self.Floodplain
            Floodplain.ExSed.OutMigration[k] = (1. - self.lambdap) * \
                Floodplain.GSD.F[k] * Floodplain.L * self.cbank * self.dxc
        
        # Set up outgoing fluxes for floodplain, active layer, and water column
        # (i.e. load) due to width change of channel zone
        
        for k in range(self.NSizes + 1):
            if self.Bcrate < 0.:
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
                    self.ActiveLayer.ExSed.OutWidthChange[k] = 0.
                    self.Load.ExSed.OutWidthChange[k] = (1. - self.lambdap) * \
                        self.FkPointBarDepositAnnualAverage[k] * self.Hpb * \
                        self.Bcrate * self.dxc
                else:
                    self.ActiveLayer.ExSed.OutWidthChange[k] = \
                        -(1 - self.lambdap) * (self.ActiveLayer.GSD.F[k] * \
                        self.ActiveLayer.L + \
                        self.FkPointBarDepositAnnualAverage[k] * self.Hpb) * \
                        self.Bcrate * self.dxc # Katie add negative sign to make the whole value positive for the UpdateFandT function
            
            else: 
                self.ActiveLayer.ExSed.OutWidthChange[k] = 0. #if widening,
                    # there is not flux from active layer to floodplain due to
                    # narrowing
    
            Floodplain = self.Floodplain
            if self.Bcrate > 0.:
                # Assumes narrowing only occurs by point bar deposition. 
                # Flux also accounts for movement of vertical boundary 
                # between active layer and floodplain.
                Floodplain.ExSed.OutWidthChange[k] = (1. - self.lambdap) \
                    * Floodplain.GSD.F[k] * Floodplain.L * self.Bcrate * \
                    self.dxc # Katie:  this is already positive and doesn't need the negtive sign (see line 546)
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
        This subroutine is necessary because the subroutine that updates tracer
        concentrations based on mass conservation at each node operates as a
        method of clsReservoir and thus does not have access to the tracer 
        concentrations in the adjacent reservoirs. 
        Consequently, these must be set up at node level.
    
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
                if self.Bcrate < 0.:
                    if k == 0:
                        self.Load.ExTracer[L].OutWidthChange[k] = \
                            self.Load.TTemporalAvMudFeed[L]
                    else:
                        self.ActiveLayer.ExTracer[L].OutWidthChange[k] = \
                            self.ActiveLayer.T[k, L]
                Floodplain = self.Floodplain
                if self.Bcrate > 0.:
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
        This subroutine considers only fluxes associated with boundary movement
        of bottom of active layer/floodplain. Net fluxes associated with 
        overbank deposition are handled in UpdateSedimentSourcesAndSinks
    
        This routine results in a source of mud to the water column from the 
        bed if the channel is degrading into a substrate that contains mud. 
        Only the duration-averaged supply rate of mud to the water column is 
        specified.
        
        Arguments:
            BedAggredationRate -- float 
            alpha -- float
        """
        
        for k in range(self.NSizes + 1):
            
            # Outgoing flux from active layer
            
            if k == 0.:
                # no mud infiltrates for now, but this should be updated if 
                # possible
                self.ActiveLayer.ExSed.OutVerticalChange[0] = 0.
            else:
                if BedAggradationRate > 0.:
                    self.ActiveLayer.ExSed.OutVerticalChange[k] = (alpha * \
                        self.ActiveLayer.GSD.F[k] + (1. - alpha) * \
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
                self.Substrate[-1].C.ExSed.OutVerticalChange[k] = \
                    -self.Substrate[-1].C.GSD.F[k] * \
                    BedAggradationRate * self.Bc * self.dxc * \
                    (1. - self.lambdap)
            else:
                self.Substrate[-1].C.ExSed.OutVerticalChange[k] = 0.
            
            
            # Setup incoming fluxes
            
            self.Substrate[-1].F.ExSed.InVerticalChange[k] = \
                self.Floodplain.ExSed.OutVerticalChange[k]
            self.Substrate[-1].C.ExSed.InVerticalChange[k] = \
                self.ActiveLayer.ExSed.OutVerticalChange[k]
            self.Floodplain.ExSed.InVerticalChange[k] = \
                self.Substrate[-1].F.ExSed.OutVerticalChange[k]
            self.ActiveLayer.ExSed.InVerticalChange[k] = \
                self.Substrate[-1].C.ExSed.OutVerticalChange[k]

        # set up mud fluxes to water column due to incision
        self.Load.ExSed.InVerticalChange[0] = \
            self.Substrate[-1].C.ExSed.OutVerticalChange[0]
        # this assumes mud does NOT infiltrate into bed if bed is aggrading

    def UpdateVerticalExchangeTracerConcentrations(self):
        """
        This subroutine is necessary because the subroutine that updates tracer
        concentrations based on mass conservation at each node operates as a 
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
    
    def UpdateNetSedimentSourcesAndSinks(self):
        """
        Floodplain deposition is handled as a sink here
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
                Outflux += self.Load.QsAvkLoad[k]

        # **********************End computation of local net flux *************
        
        # Net source/sink fluxes to all reservoirs except water column are 
        # computed based on averages across the duration curve
        for k in range(1, self.NSizes + 1):
            self.ActiveLayer.SourceLatSed[k] = self.SLatSourceAv[k]
            # net shaving included in the lateral exchange fluxes
            self.ActiveLayer.SinkLatSed[k] = self.Dfav[k] * self.dxc 
                # Additional sinks such as dredging could be added here.
            self.ActiveLayer.SourceFeedSed[k] = self.Load.QsAvkFeed[k]
    
            # **********MASS CONSERVATION IF BED ELEVATION IS FIXED************
        
            if self.FixedElev:
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
                if NetInfluxTotal > Outflux:
                    # if the system is gaining sediment, make sure the sediment
                    # stored is similar to the size distribution of the active 
                    # layer
                    # The sediment passed downstream in the load is then equal 
                    # to the total influx in size k minus the sediment in size
                    # k that is stored.
                    self.ActiveLayer.SinkLoadSed[k] = NetInflux[k] - \
                        (NetInfluxTotal - Outflux) * \
                        (self.ActiveLayer.GSD.F[k] * \
                        self.AlphaPartlyAlluvial + \
                        self.Load.GSDBedloadAv.F[k] * \
                        (1. - self.AlphaPartlyAlluvial))
                else:
                    # If the system is sedient starved, assume the only outflux
                    # is that computed from the fraction of the bed that is 
                    # alluvial.
                    self.ActiveLayer.SinkLoadSed[k] = self.Load.QsAvkLoad[k]
            
                if self.ActiveLayer.SinkLoadSed[k] < 0.:
                    self.ActiveLayer.SinkLoadSed[k] = 0.
#                    print('negative size specific load in partly alluvial \ # Katie comment out to decrease computation time
#                        computation in UpdateSedimentSourcesAndSinks of Node \
#                        was converted to zero')
            else:
                self.ActiveLayer.SinkLoadSed[k] = self.Load.QsAvkLoad[k]
    
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
        This subroutine is necessary because the mass conservation at each node
        operates as a method of clsReservoir and thus does not have access to 
        the node's lateral sources and sinks.
        
        This routine also computes the duration-averaged tracer concentration 
        in the load.
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
        This routine computes overbank deposition rates (m/s), outgoing 
        suspended mud flux (m3/s), and fraction mud in point bars (no units) 
        for a single set of steady discharge and sediment supply terms at a 
        single node.  All terms are assumed constant over the interval of the 
        computation.
        If the computed outgoing suspended sediment flux is less than zero, the
        deposition rate is reduced so that the outgoing sediment flux is zero 
        and mass is conserved.
        
        Arguments:
            C -- float (Average suspended sediment concentration above 
                floodplain level)
            Qwf -- float (Water discharge across floodplain (m3/s))
            Fl -- float (Floodplain number (can be size specific))
            Bf -- float (Floodplain width (m))
        
        Return: float
        """

        return Fl * C * Qwf / Bf
        # Note that D is a volume rate of sediment particles deposited on 
        # floodplain per unit channel length, so has units of L2/T.  It is not
        # the vertical rate of change floodplain elevation, and thus does not 
        #include a porosity term.
    
    def RouseFractionSuspendedLoadAboveFloodplainLevel(self, SettlingVel, \
        ustar, Hc, Lf, Intervals):
        """
        Units of settling velocity and ustar must be consistent
        Units of Hc and Hf must be consistent
        
        Arguments:
            SusSedLoad -- float (Volumetric suspended sediment load in size 
                class of interest)
            Qc -- float (Water discharge in channel)
            SettlingVel -- float (Settling velocity for suspended sediment)
            Ustar -- float (Shear velocity in channel)
            Hc -- float (Flow depth in channel)
            Lf -- float (height of floodplain with respect to channel bed)
            Intervals -- float (number of intervals over which integral is 
                numerically integrated.)
        
        Return: float
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
                print('floodplain elevation is below near bed elevation so \
                    Rouse Profile can not be computed.  Concentration set to \
                    load/discharge.')
            
            # Integrate for lower part of water column up to floodplain level.
            # Assume near bed concentration is evaluated at zeta = 0.05 and
            # that no suspended transport occurs below this level.
            dzita = (Lf / Hc - 0.05) / Intervals
            sum1 = 0. # Katie change variable name from sum to sum1 to avoid using the name of a built-in function.
            for n in range(Intervals):
                ZitaDummy = 0.5 + n * dzita - dzita / 2.
                sum1 = sum1 + dzita * self.RouseIntegrand(ZitaDummy, 0.05, \
                    SettlingVel, ustar)
            IntegralBelowFloodplainLevel = sum1
            
            return 0.95 / (1. - (Lf / Hc)) * IntegralAboveFloodplainLevel / \
                (IntegralAboveFloodplainLevel + IntegralBelowFloodplainLevel)

    def RouseIntegrand(self, Zita, Zitab, Vs, ustar):
        """
        Arguments:
            Zita -- float
            ZItab -- float
            Vs -- float
            ustar -- float
            
        Return: float
        """
        return ((1. - Zita) * Zitab / ((1. - Zitab) * Zita)) ** (Vs / (0.4 * \
            ustar))
    
    def UpdateDMudj(self):
        """
        Updates mud deposition rate for each bin and computes flow-duration 
        averaged value of overbank mud deposition rate (k = 0).
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
        Updates sand deposition rate for each bin and computes flow-duration 
        averaged value of overbank sand deposition rate (k = represntative of
        all sand sizes).
        """
        g = 9.81
            
        for k in range(1, self.ActiveLayer.GSD.NBedSizes + 1):
            # Only used for sand size classes, but could be used for gravel too
            # if shear velocity was high enough. 
            # Actually, this routing could work for mud too as long as a 
            # reasonable settling velocity was available.
            # It is written as a separate subroutine since the floodplain 
            # numbers for sand and mud can be different.

            if self.ActiveLayer.GSD.D[k] < 2.:
                self.Dfav[k] = 0.
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
                        CRatio = self.\
                            RouseFractionSuspendedLoadAboveFloodplainLevel(\
                            self.ActiveLayer.GSD.Vs[k], ustar, self.DC.Hc[j], \
                            self.Floodplain.L - self.ActiveLayer.L, 20.)
    
    def WidthChange(self, Manning, thresholdQ, W, ErodeT, Bp, dt):
        """
        Katie added this function.  It is an optional add-on to replace the constant
        migration rate (Node.cbank).  Bank erosion and floodplain sequestration are
        treated as separate, independent processes.
        
        Arguments:
        
        Manning -- bool
        thresholdQ -- float
        W -- float
        ErodeT -- float
        Bp -- float
        dt -- float
        """
        
        self.cbank = 0. # Turns off constant migration rate so that ExchangeType values can be filled with this function        
        
        # Bank erosion function.  This is based on the work from Sarah Davidson's
        # PhD thesis but modified to fulfill our mass conservation requirements
        # and smaller timestep.  The channel is widened incrementally until the
        # D84 is no longer fully mobile (full mobility is defined as twice the 
        # critical shear stress for the D84--see Sarah's thesis for details.)
        
        # In order to calculate the shear stress on the widened bed, we need to
        # recalculate hydraulics.  For simplicity, normal flow is assumed within
        # each node.
        
        # Right now, this also only works with a hydrograph (or single flow in the
        # duration curve).  Would have to modify it to use it with a flow duration
        # curve.
        
        HybridD84 = (1./3.)*self.ActiveLayer.GSD.D84 + (2./3.)*self.Floodplain.GSD.D84 # Takes average of floodplain and active layer D84s    .        
        
        Tau84 = self.Load.findWC_ReferenceShear(self.Load.TauRm, HybridD84, self.ActiveLayer.GSD.D50)       
        TauFull = Tau84 # Shear stress where the D84 reaches full mobility
        
        oldBc = deepcopy(self.Bc) # Must save to use the old channel width for the conservation of mass later on.        
        
        while TauFull < self.Load.TauPrime:

            self.Bc = self.Bc + .001*oldBc # Widens channel by .1% each iteration. Percentage was arbitrary.            	
            self.UpdateDepthAndDischargeAtAllFlows(Manning)
            self.Load.TauPrime = self.Load.findWC_TauPrime(9.81, self.DC.Sf[0], self.DC.Uc[0], self.ActiveLayer.GSD.D65)
        
		#  Will cause bank erosion if floodplain height gets above certain threshold.
        #L = self.Floodplain.L
        #while L - self.ActiveLayer.L > 3.:

        #    self.Bc = self.Bc + .001*oldBc # Widens channel by .1% each iteration. Percentage was arbitrary.            	
        #    eroded = .001 * self.Hpb
        #    dz = eroded / (1. - self.lambdap)
        #    L = L - dz
			
        widthchangeE = (self.Bc - oldBc)/ErodeT # Reduces widening rate to assume that widening occurs over the course of 12 hours. Unit is width change/second
        if oldBc + widthchangeE*dt >= self.Bf: # Makes sure that width does not exceed valley width
            widthchangeE = (self.Bf-oldBc)/dt # Adjusts Bcrate to max widening
            
        # Bank depostion function.  This is a somewhat arbitrary equation that is
        # designed to be calibrated to observed narrowing rates from air photos.
        # The idea is that the channel narrows through enroachment of vegetation
        # onto point bars.  The rate of enroachment is a function of the width
        # of the point bar (which is defined as an arbitrary percentage of the 
        # channel width) and the relative particle size of the bar (right now as 
        # a ratio between the active layer and substrate D50--should later change
        # to point bar D50).  Narrowing occurs at a user-defined rate and is faster
        # when 1) the rate is higher, 2) the grainsize distribution of the active
        # layer is smaller (so that it is easier for vegetation to take hold), and/
        # or 3) the width of available bar is greater.  Narrowing is only assumed
        # to occur at low flows, which the user defines based on a flow threshold.
                
        #if self.DC.Qw[0] < thresholdQ and self.DC.Qwf[0] <= 0: # Narrows when discharges drops below a threshold but is not flooded by a reservoir
        widthchangeN = 0.
        #Calculate sand fraction to normalize narrowing rate by sand fraction.		

        if self.DC.Hc[0] < .25*self.Floodplain.L and oldBc > 40: #.465 # Make narrowing a function of channel depth, not discharge.  Narrows when flow depth is quarter bankfull.
#            Fs = 0. # Adjust narrowing rate by sand content--not currently in use
#            i = 0            
#            while self.ActiveLayer.GSD.D_upper[i] <= 2.: # Katie change from .002 (which is the lower boundary)
#                Fs = self.ActiveLayer.GSD.FfinerThanD_upper[i] # Katie change from parentheses to brackets
#                i = i + 1		
            widthchangeN = -(self.Bc-40.)*W/(365.25*24*60*60) #Katie multiply by percent sand # Change per second--simplest solution; only one calibration term which is percent widening per year.

        self.Bcrate = widthchangeE + widthchangeN
        self.Bc = oldBc

        
    def Avulsion(self, threshold, spacing):
        
        """
        Katie added this function.
        This function simulates avulsion on a node-by-node basis.  When the channel
        elevation approaches that of the bank (how close is a user-defined threshold),
        the uppermost substrate layer will be incorporated back into the floodplain.
        The active layer will also join the Floodplain reservoir.  The bed elevation
        (etabav) will be adjusted so that it is lowered by an increment of 1 substrate
        spacing.
        
        If this occurs repeatedly in a given node, it should simulate valley filling,
        like that of a delta.
        
        Attributes:
        
            -threshold (float)--difference between the bank (Floodplain.L-ActiveLayer.L)
                the channel elevation after which avulsion will occur
        """
        
        if self.Floodplain.L-self.ActiveLayer.L < threshold:
            print 'Avulsion!'
            print self.Floodplain.L
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
            # the old active layer, and both substrate channel and floodplain reservoirs.
            # An active-layer size chunk is taken out to form the new active layer.
            NewFpF = np.zeros(self.NSizes+1)
            NewFpT = np.zeros((self.NSizes+1, self.NTracers))
            
            for k in range(self.NSizes+1):
                NewFpF[k] = (self.Floodplain.GSD.F[k]*self.Floodplain.Volume + \
                    self.Substrate[-1].C.GSD.F[k]*self.Substrate[-1].C.Volume +\
                    self.Substrate[-1].F.GSD.F[k]*self.Substrate[-1].F.Volume +\
                    self.ActiveLayer.GSD.F[k]*self.Bc*self.dxc*self.ActiveLayer.L -\
                    self.Floodplain.GSD.F[k]*self.Bc*self.dxc*self.ActiveLayer.L)/\
                    NewFpVolume
                    
                # Do the same thing for tracers                
                for m in range(self.NTracers):
                    NewFpT[k, m] = (self.Floodplain.T[k, m]*self.Floodplain.Volume + \
                    self.Substrate[-1].C.T[k, m]*self.Substrate[-1].C.Volume +\
                    self.Substrate[-1].F.T[k, m]*self.Substrate[-1].F.Volume +\
                    self.ActiveLayer.T[k, m]*self.Bc*self.dxc*self.ActiveLayer.L -\
                    self.Floodplain.T[k, m]*self.Bc*self.dxc*self.ActiveLayer.L)/\
                    NewFpVolume
            
            # The active layer is given the GSD and tracer content of the (old) floodplain
            # Might think about changing this later--would be technically more
            # correct to have GSD a combination of Floodplain and Substrate,
            # but the difference is likely to be miniscule.
            self.ActiveLayer.GSD.F = deepcopy(self.Floodplain.GSD.F)
            self.ActiveLayer.T = deepcopy(self.Floodplain.T)
            self.ActiveLayer.GSD.UpdateStatistics()
            
            # Update floodplain reservoir 
            self.Floodplain.Volume = NewFpVolume
            self.Floodplain.GSD.F = NewFpF
            self.Floodplain.T = NewFpT
            self.Floodplain.GSD.UpdateStatistics()
            
            # Lower channel elevation to base of former substrate level 
            self.etabav = self.etabav-self.Substrate[-1].C.L
            self.Floodplain.L = self.Floodplain.Volume / (self.Bf * self.dxc / \
                self.ChSin)
            
            # Get rid of old uppermost substrate
            del self.Substrate[-1]

    def ExnerBed(self):

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

#        NetQsIn = np.sum(self.ActiveLayer.ExSed.InMigration[1:self.NSizes + 1]\
#            - self.ActiveLayer.ExSed.OutMigration[1:self.NSizes + 1] + \
#            self.ActiveLayer.ExSed.InWidthChange[1:self.NSizes + 1] + \
#            self.ActiveLayer.SourceLatSed[1:self.NSizes + 1] - \
#            self.ActiveLayer.SinkLatSed[1:self.NSizes + 1] + \
#            self.ActiveLayer.SourceFeedSed[1:self.NSizes + 1] - \
#            self.ActiveLayer.SinkLoadSed[1:self.NSizes + 1])
             
        # This is where a formulation accounts for the possibility that the bed
        # is partially alluvial:
        if self.FixedElev:
            self.DeltaEtaB = 0.
        else:
            self.DeltaEtaB = NetQsIn / self.dxc / self.Bc / (1. - self.lambdap)
            

#  Katie: This didn't work.
#        # Katie add:  if bed is aggrading, some will go into channel narrowing.  
#        # To be used when channel narrowing is not in the Exner equation.
#
#        ExtraNarrowing = 0. # Katie add
#        if self.DeltaEtaB > 0.:
#            ExtraNarrowing = (.5*self.DeltaEtaB*self.Bc*self.dxc)/(self.Hpb*self.dxc)
#            self.Bcrate = self.Bcrate - ExtraNarrowing
#        
#        #  Add extra material back to floodplain
#        for k in range(self.NSizes + 1):
#            if ExtraNarrowing > 0.:
#                self.ActiveLayer.ExSed.OutWidthChange[k] = \
#                    self.ActiveLayer.ExSed.OutWidthChange[k] + ExtraNarrowing*(1 * \
#                    self.ActiveLayer.GSD.F[k] + (1. - 1) * \
#                    self.Load.GSDBedloadAv.F[k])
                    
        
        # Katie add variables for cumulative channel volume change--as an alternative
        # to CumulativeBedChange when the Avulsion fuction is working.
        self.CumulativeChannelVolumeChange = self.CumulativeChannelVolumeChange + NetQsIn


    def UpdateVolumeSizeFractionsandTracersInAllReservoirs(self, dt, \
        TracerProperties):
        """
        Arguments:
            dt -- float
            TracerProperties -- clsTracerProperties
        """
        self.ActiveLayer.UpdateFandT(dt, self.lambdap, TracerProperties)
        self.Floodplain.UpdateFandT(dt, self.lambdap, TracerProperties)
        for m in range(self.NLayers()):
            self.Substrate[m].C.UpdateFandT(dt, self.lambdap, TracerProperties)
            self.Substrate[m].F.UpdateFandT(dt, self.lambdap, TracerProperties)
        # Should include tracer production and decay either here or in a
        # separate sub
    
    def UpdateGeometricParameters(self, dt):
        """
        Arguments:
            dt -- float
        """
        
        self.etabav = self.etabav + self.DeltaEtaB * dt
            
        #self.etabav = round(self.etabav,4) # Katie add
        #self.Bc = self.Bc + self.Bcrate * dt # Katie:  take this and move it to WidthChange function to put it in correct place before Exner equation is done.        
        #self.Bf = self.Bf - self.Bcrate * self.ChSin * dt # Katie:  take this and move it to WidthChange function to put it in correct place before Exner equation is done.
        
        self.Bc = self.Bc + self.Bcrate * dt
        self.Bf = self.Bf - self.Bcrate * self.ChSin * dt
        
                    
        # Katie add to keep track of how much narrowing and bank erosion occurs
        if self.Bcrate < 0:
            self.CumulativeNarrowing = self.CumulativeNarrowing + self.Bcrate*dt
        if self.Bcrate > 0:
            self.CumulativeWidening = self.CumulativeWidening + self.Bcrate*dt
            self.CumulativeWideningVolume = self.CumulativeWideningVolume + (self.Bcrate*dt*self.Floodplain.L*self.dxc)
            
       # Katie add--keep track of total feed
        self.CumulativeTotFeed = self.CumulativeTotFeed + (np.sum(self.Load.QsAvkFeed)*dt)
#        if self.Bcrate < 0:
#            self.Bc = self.Bc + self.Bcrate * dt # Note that this and the following line used to be in the UpdateGeometricParameters function but were moved because the width needs to change before the Exner equation is performed.        
#            if self.Bc < 35.:
#                self.Bc = 35.
#            self.Bf = self.Bf - self.Bcrate * self.ChSin * dt 
        
        self.Floodplain.L = self.Floodplain.Volume / (self.Bf * self.dxc / \
            self.ChSin)
        self.H = self.Floodplain.L
        
#        # Katie add to try to prevent crashing by delta in Elwha run
#        if self.Floodplain.L-self.CumulativeBedChange < .1:
#            self.etabav = self.etabav - 1.
        
        self.Substrate[-1].C.L = self.Substrate[-1].C.L + self.DeltaEtaB * dt
        self.Substrate[-1].F.L = self.Substrate[-1].F.L + self.DeltaEtaB * dt
        self.Substrate[-1].C.Volume = self.Substrate[-1].C.L * self.Bc * \
            self.dxc
        self.Substrate[-1].F.Volume = self.Substrate[-1].F.L * self.Bf * \
            self.dxc / self.ChSin
        if abs(self.Substrate[-1].F.L - self.Substrate[-1].C.L) > 0.001:
            print('Problem with substrate thickness')
        
        # This method should probably be modified for a partly alluvial case.
        # If dt is so long that channel incises through all substrate, then 
        # need to reduce alluvial cover fraction.  Conversely, for partly 
        # alluvial case, if bed cover fraction would increase to greater than
        # 1, then need to allow bed elevation to change and also to compute 
        # appropriate vertical boundary exchange fluxes.  Should be easy to 
        # write if-then statements that modify boundary flux compuations.
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
