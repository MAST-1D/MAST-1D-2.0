
from clsExchangeTypes import clsExchangeTypes
from clsGSD import clsGSD 
import numpy as np
from copy import deepcopy


class clsReservoir:
    """
    Attributes:
        L -- float (Thickness)
        T -- [float] (Tracer concentration in size k)
        ExSed -- clsExchangeTypes (Size specific sediment fluxes into and out
                 of reservoir (units of volume/time))
        ExTracer -- [clsExchangeTypes] (Size and tracer specific tracer 
                    concentrations in sediment moving into and out of 
                    reservoir)
        SourceLatSed -- [float] (Lateral Sediment Source (e.g. bluffs, etc.))
        SourceLatTracer -- [float] (Tracer concentration in net lateral tracer
                           source)
        SinkLatSed -- [float] (Sediment flux to lateral sink (e.g. floodplain
                      deposition, dredging, etc.))
        SinkLatTracer -- [float] (Tracer concentraiton in lateral sediment sink
                         --may not be necessary since outgoing T is just T)
        SourceFeedSed -- [float] (Upstream sediment feed)
        SourceFeedTracer -- [float] (Tracer concentration in upstream sediment
                            feed)
        SinkLoadSed -- [float] (Sediment flux out of node in load)
        SinkLoadTracer -- [float] (Tracer concentration in load)
        NTracers -- int (number of tracers)
        NSizes -- int (number of bed material grain sizes.  note that k = 0 
                  represents washload.)
        GSD -- clsGSD (grain size distribution for reservoir)
        Volume -- float (total reservoir volume including voids, for all sizes)
        
        Initialized -- bool
    """
    Initialized = False
    
    def getNTracers(self):
        return self.T.shape[1]
    
    NTracers = property(getNTracers)
    
    def getNSizes(self):
        return self.ExSed.NSizes
    
    NSizes = property(getNSizes)
    
    def __init__(self, BinBdySizes, NumTracers):
        """
        Arguments:
            BinBdySizes -- [float]
            NTracers -- int
        """
        if not self.Initialized:
            #self.NSizes = len(BinBdySizes) - 2
            sizes = len(BinBdySizes)-2
            #self.NTracers = NTracers
            
            #self.ExSed = clsExchangeTypes(self.NSizes)
            self.ExSed = clsExchangeTypes(sizes)
            self.ExTracer = [clsExchangeTypes(sizes) for i in range(NumTracers)]
            self.GSD = clsGSD(BinBdySizes)
            self.T = np.zeros((sizes + 1, NumTracers))
            self.Initialized = True
            self.DeltaS = 0. # Katie add to keep track of mass conservation
            self.L = 0.
            self.Volume = 0.
            
            self.SourceLatSed = np.zeros(sizes + 1)
            self.SourceLatTracer = np.zeros((sizes + 1, NumTracers))
            self.SinkLatSed = np.zeros(sizes + 1)
            self.SinkLatTracer = np.zeros((sizes + 1, NumTracers))
            self.SourceFeedSed = np.zeros(sizes + 1)
            self.SourceFeedTracer = np.zeros((sizes + 1, NumTracers))
            self.SinkLoadSed = np.zeros(sizes + 1)
            self.SinkLoadTracer = np.zeros((sizes + 1, NumTracers))
        else:
            raise RuntimeError('Tried to initiate clsReservoir twice.')
    
    def UpdateFandT(self, dt, lambdap, TracerProperties):
        """
        Arguments:
            dt -- float
            lambdap -- float
            TracerProperties -- clsTracerProperties
        """
        DeltaSedVolume = np.zeros(self.NSizes + 1)
        DeltaTVolume = np.zeros((self.NSizes + 1, self.NTracers))
        OldSedVolume = np.zeros(self.NSizes + 1)
        OldTVolume = np.zeros((self.NSizes + 1, self.NTracers))
        NewSedVolume = np.zeros(self.NSizes + 1)
        NewTVolume = np.zeros((self.NSizes + 1, self.NTracers))
        #print self.ExSed.OutVerticalChange
        NewSedVolumeTotal = deepcopy(self.Volume)
        for k in range(self.NSizes + 1):
            # Determine change in sediment volume in each size class
            DeltaSedVolume[k] = (self.ExSed.InMigration[k] + \
                self.ExSed.InWidthChange[k] + self.ExSed.InVerticalChange[k] \
                - self.ExSed.OutMigration[k] - self.ExSed.OutWidthChange[k] - \
                self.ExSed.OutVerticalChange[k] + self.SourceLatSed[k] - \
                self.SinkLatSed[k] + self.SourceFeedSed[k] - \
                self.SinkLoadSed[k]) * dt / (1. - lambdap)
            NewSedVolumeTotal = NewSedVolumeTotal + DeltaSedVolume[k]      

            # Determine change in volume of tracer sediment in each size class
            # for each tracer
            if self.NTracers > 0.:
                for L in range(self.NTracers):
                    DeltaTVolume[k, L] =(self.ExSed.InMigration[k] \
                        * self.ExTracer[L].InMigration[k] \
                        + self.ExSed.InWidthChange[k] \
                        * self.ExTracer[L].InWidthChange[k] \
                        + self.ExSed.InVerticalChange[k] \
                        * self.ExTracer[L].InVerticalChange[k] \
                        - self.ExSed.OutMigration[k] \
                        * self.ExTracer[L].OutMigration[k] \
                        - self.ExSed.OutWidthChange[k] \
                        * self.ExTracer[L].OutWidthChange[k] \
                        - self.ExSed.OutVerticalChange[k] \
                        * self.ExTracer[L].OutVerticalChange[k] \
                        + self.SourceLatSed[k] * self.SourceLatTracer[k, L] \
                        - self.SinkLatSed[k] * self.SinkLatTracer[k, L] \
                        + self.SourceFeedSed[k] * self.SourceFeedTracer[k, L] \
                        - self.SinkLoadSed[k] * self.SinkLoadTracer[k, L]) \
                        * dt / (1. - lambdap)

        if self.NTracers > 0.:
            for L in range(self.NTracers):
                for k in range(self.NSizes + 1):
                    OldTVolume[k, L] = self.Volume * self.GSD.F[k] \
                        * self.T[k, L]
                    NewTVolume[k, L] = OldTVolume[k, L] + DeltaTVolume[k, L]
        # Set Grain size fractions, Tracer Concentrations, and Total Volume
                    
        for k in range(self.NSizes + 1):
            OldSedVolume[k] = self.Volume * self.GSD.F[k]
            NewSedVolume[k] = OldSedVolume[k] + DeltaSedVolume[k]
            self.GSD.F[k] = NewSedVolume[k] / NewSedVolumeTotal
            if self.GSD.F[k] < 0.:
#                print k
                self.GSD.F[k] = 0.
            if self.NTracers > 0.:
                for L in range(self.NTracers):
                    # Tracer concentration is defined as tracer volume divided
                    # by sediment volume in a given size class
                    if NewSedVolume[k] > 0.:
                        self.T[k, L] = NewTVolume[k, L] / NewSedVolume[k]
                    else:
                        self.T[k, L] = 0.        
        self.DeltaS = NewSedVolumeTotal - self.Volume
        self.Volume = NewSedVolumeTotal
        self.GSD.UpdateStatistics()
        
        # Account for tracer production and decay

        # Note that production functions are not yet implemented.  For 
        # cosmogenic production, this will require knowledge of the total mass
        # of cover above the reservoir. It is not clear how to implement this.  
        # Perhaps this could be computed at node level and then passed to the 
        # reservoir as part of the call for the method.  It may also save time 
        # to implement this as a new subroutine since then it would be possible 
        # to allow production only in the reservoirs where it is relevant.
        
        for L in range(self.NTracers):
            for k in range(self.NSizes + 1):
                self.T[k, L] -= TracerProperties[L].DecayConst * dt / \
                    (365.25 * 24. * 60. * 60.)