import numpy as np

class clsDurationCurve:
    """
    Defines discharge distribution.

    Stores total discharge and various hydraulic parameters associated with
    various discharge bins.

    Parameters
    ----------
    NFlows : int
        Number of discharge bins.
        
    Attributes
    ----------
    Qw : array_like(float, length = NFlows) 
        Water Discharge in bin.
    p : array_like(float, length = NFlows)
        Fraction of time bin occurs.
    Uc : array_like(float, length = NFlows) 
        Water velocity in channel in bin.
    Uf : array_like(float, length = NFlows)
        Water velocity on floodplain in bin.
    Hc : array_like(float, length = NFlows)
        Water depth in channel in bin.
    Hf : array_like(float, length = NFlows)
        Water depth on floodplain in bin.
    Qwc : array_like(float, length = NFlows)
        Discharge of water in the channel in bin.
    Qwf : array_like(float, length = NFlows)
        Discharge of water on the floodplain in bin.
    Qs : array_like(float, length = NFlows)
        Sediment discharge in bin. POTENTIAL CODE
        UPDATE: probably should remove
        this since all load is handled in clsLoad
    WSE : array_like(float, length = NFlows)
        Water surface elevation.
    Sf : array_like(float, length = NFlows)
        Friction Slope.
    Dfloodplain : array_like(float, length = NFlows)
        Floodplain deposition in bin--should check 
        whether this is stored elsewhere.
    MigrationFactor : array_like(float, length = NFlows)
        Factor that determines how much of the
        annual average migration occurs for each of duration curve.
    Initialized : bool
        Flag to determine if object was initialized.
    """
    Initialized = False

    def __init__(self, NFlows):
        if not self.Initialized:
            self.Qw = np.zeros(NFlows)
            self.p = np.zeros(NFlows)
            self.Uc = np.zeros(NFlows)
            self.Uf = np.zeros(NFlows)
            self.Hc = np.zeros(NFlows)
            self.Hf = np.zeros(NFlows)
            self.Qs = np.zeros(NFlows)
            self.Qwc = np.zeros(NFlows)
            self.Qwf = np.zeros(NFlows)
            self.Dfloodplain = np.zeros(NFlows)
            self.WSE = np.zeros(NFlows)
            self.Sf = np.zeros(NFlows)
            self.Migrationfactor = np.zeros(NFlows)
            self.Initialized = True
        else:
            raise RuntimeError('Tried to initiate clsDurationCurve twice.')

    def NFlows(self):
        return len(self.p)

    def FroudChannel(self, j):
        """
        Froude number in channel in jth discharge bin.
        
        Parameters
        ----------
        j : int
            Index of discharge bin.
            
        Returns
        -------
        float   
            Froude number.
        """
        g = 9.81
        return self.Uc[j] / (g * self.Hc[j]) ** 0.5

    def WeightByFDC(self, x):
        """
        Weights a variable defined at all flows in the discharge distribution
        by the fractions of time each discharge in the distribution occurs.
        
        Parameters
        ----------
        x : array_like(length = NFlows)
            Variable to be weighted.
            
        Returns
        -------
        float   
            Weighted variable.
        """
        WeigthByFDC = 0.
        for i in range(len(p)):
            WeightByFDC = x[j] * self.p[j] + WeightByFDC
        return WeightByFDC