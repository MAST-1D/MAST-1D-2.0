from clsExchangeTypes import clsExchangeTypes
from clsGSD import clsGSD 
import numpy as np
from copy import deepcopy

class clsReservoir: 
    """
    A class defining a sediment storage reservoir.
    
    Defines the volume and nature of sediment stored within a given 
    sediment storage reservoir and includes the method for performing size-specific 
    mass balance computations for sediment and tracer mass.
    
    Parameters
    ----------
        BinBdySizes : array_like(float)
            Sediment grain size at each bin boundary. Note that the number
            of sediment sizes NSizes = BinBdySizes - 2.
        NTracers : int
            Number of tracers.
    
    Attributes
    ----------
        NTracers : int 
            Number of tracers.
        NSizes : int 
            Number of bed material grain sizes. Note that k = 0 represents washload.)
        GSD : :obj:`MAST_1D.clsGSD`
            Grain size distribution of sediment for reservoir.
        Volume : float 
            Total reservoir volume including voids, for all sizes.
        L : float 
            Vertical thickness of reservoir
        T : array_like (float, length = NTracers)
            Tracer concentration in size k.
        ExSed : :obj:`MAST_1D.clsExchangeTypes`
            Size specific sediment fluxes representing processes considered in the 
            clsExchangeTypes class moving into and out of the reservoir. 
            (units of volume/time).
        ExTracer : array_like(:obj:`MAST_1D.clsExchangeTypes`, length = NTracers) 
            Size and tracer specific tracer concentrations in sediment
            moving into and out of reservoir. Defined for each tracer type.
        SourceLatSed : array_like(float, length=Nsizes)
            Size-specific lateral sediment source (e.g. bluffs, etc.).
        SourceLatTracer : array_like (float, ndim = 2)
            Tracer concentration in net lateral tracer source for each tracer type
            and sediment size.  Size of array is NSizes x NTracers.
        SinkLatSed : array_like(float, length=Nsizes)
            Size specific sediment flux to lateral sink (e.g. floodplain deposition,
            dredging, etc.).
        SinkLatTracer : array_like (float, ndim = 2)
            Tracer concentration in lateral sediment sink for each size and tracer time.
            Size of array is NSizes x NTracers.
            May not be necessary since for cosmogenic tracer, loss is proportional to T.
        SourceFeedSed : array_like(float, length=Nsizes)
            Upstream sediment feed for each sediment size.
        SourceFeedTracer : array_like (float, ndim = 2)
            Tracer concentration in upstream sediment feed, for each 
            sediment size and tracer type. Size of array is NSizes x NTracers.
        SinkLoadSed : array_like(float, length=Nsizes)
            Sediment flux out of node in load for each sediment size.
        SinkLoadTracer : array_like (float, ndim = 2)
            Tracer concentration in load for each 
            sediment size and tracer type. Size of array is NSizes x NTracers
        Initialized : bool
            Flag to determine if reservoir has been initialized
            
    """
    Initialized = False
    
    def getNTracers(self):
        return self.T.shape[1]
    
    NTracers = property(getNTracers)
    
    def getNSizes(self):
        return self.ExSed.NSizes
    
    NSizes = property(getNSizes)
    
    def __init__(self, BinBdySizes, NumTracers):
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
        Apply mass conservation to the reservoir to estimate sediment size fractions and
        tracer concentrations in future timestep.
        
        Arguments
        ---------
            dt : float
                Timestep (seconds).
            lambdap : float
                Porosity of reservoir deposit.
            TracerProperties : array_like(:obj:clsTracerProperties, length = NTracers)
                Production and decay properties for each tracer type.
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