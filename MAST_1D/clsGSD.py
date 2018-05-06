import numpy as np
from math import log10, log

class clsGSD(object):
    """
    Attributes:
        D_lower -- [float] (size (mm) at bin lower bin boundary)
        D_upper -- [float] (size (mm) at upper bin boundary)
        FfinerThanD_upper -- [float] (cumulative fraction finer than upper 
            limit of bin boundary)
        D -- [float] (geometric mean size in bin (mm))
        F -- [float] (fraction in bin)
        Vs -- [float] (Settling velocity of geometric mean size in bin (m/s))
        D84 -- float (Size (mm) for which 84% is finer)
        D64 -- float (Size (mm) for which 64% is finer)
        D50 -- float (median diameter (mm))
        Dg -- float (geometric mean diameter (mm))
        SigmaG -- float (geometric standard deviation (no units))
        SigmaSG -- float (standard deviation on psi scale (no units))
        Initialized -- bool
        NBedSizes -- int
        UpToDate -- bool (flag that determines whether statistics need to be 
            recomputed.  Set to false if F or D are changed.)
    """
    Initialized = False
    
    def getF(self):
        self.UpToDate = False # Any time F changes, set flag that indicates 
        # statistics need to be recomputed
        return self._F

    def setF(self, value):
        self.UpToDate = False # Any time F changes, set flag that indicates 
        # statistics need to be recomputed
        self._F = value
    
    F = property(getF, setF)
    
    def getFfinerThanD_upper(self):
        """
        This is a read only property to ensure consistency with rest of GSD
        """
        if not self.UpToDate: self.UpdateStatistics() # recompute cumulative 
        # totals and statics if necessary
        return self._FfinerThanD_upper
    
    FfinerThanD_upper = property(getFfinerThanD_upper)
    
    def getNBedSizes(self):
        return len(self.D) - 1
    
    NBedSizes = property(getNBedSizes)
    
    def getD95(self):
        """
        This is a read only property to ensure consistency with rest of GSD
        """
        if not self.UpToDate: self.UpdateStatistics() # recompute cumulative 
        # totals and statics if necessary
        return self._D95
        
    D95 = property(getD95) 
    
    def getD90(self):
        """
        This is a read only property to ensure consistency with rest of GSD
        """
        if not self.UpToDate: self.UpdateStatistics() # recompute cumulative 
        # totals and statics if necessary
        return self._D90
    
    D90 = property(getD90)    
    
    def getD84(self):
        """
        This is a read only property to ensure consistency with rest of GSD
        """
        if not self.UpToDate: self.UpdateStatistics() # recompute cumulative 
        # totals and statics if necessary
        return self._D84
    
    D84 = property(getD84)
    
    def getD65(self):
        """
        This is a read only property to ensure consistency with rest of GSD
        """
        if not self.UpToDate: self.UpdateStatistics() # recompute cumulative 
        # totals and statics if necessary
        return self._D65
    
    D65 = property(getD65)
    
    def getD50(self):
        """
        This is a read only property to ensure consistency with rest of GSD
        """
        if not self.UpToDate: self.UpdateStatistics() # recompute cumulative 
        # totals and statics if necessary
        return self._D50
    
    D50 = property(getD50)
    
    def getDg(self):
        """
        This is a read only property to ensure consistency with rest of GSD
        """
        if not self.UpToDate: self.UpdateStatistics() # recompute cumulative 
        # totals and statics if necessary
        return self._Dg
    
    Dg = property(getDg)
    
    def getSigmag(self):
        """
        Geometric standard deviation, transformed to length scale (mm)
        """
        if not self.UpToDate: self.UpdateStatistics() # recompute cumulative 
        # totals and statics if necessary
        return self._Sigmag
    
    Sigmag = property(getSigmag)
    
    def getSigmaSG(self):
        """
        Standard deviation of sediment on psi scale (no units)
        """
        if not self.UpToDate: self.UpdateStatistics() # recompute cumulative 
        # totals and statics if necessary
        return self._SigmaSG
    
    SigmaSG = property(getSigmaSG)
    
    def D_x(self, x):
        """
        Geometrically interpolates the kth percentile of a grain size 
        distribution Uses a grain size distribution and the kth percentila 
        Provides as output the diameter for which k percent of the distribution
        is finer, as interpolated in log space
        
        Arguments:
            x -- float
        """
        if x < 0. or x > 100.:
            raise RuntimeError('Error:  Index of size distribution out of \
                allowable range in clsGSD.D_x')
        
        if not self.UpToDate:
            self.UpdateFfiner()
        
        n = len(self.D)
        psi_low = 0.
        psi_high = 0.
        F_high = 0.
        for k in range(n): 
            if k == 0: 
                F_low = 0.
            else:
                F_low = self._FfinerThanD_upper[k - 1]
            if x / 100. > F_low and x / 100. <= self._FfinerThanD_upper[k]:
                psi_low = log(self.D_lower[k]) / log(2.)
                psi_high = log(self.D_upper[k]) / log(2.)
                F_high = self._FfinerThanD_upper[k]
                break
        
        psi_x = (psi_high - psi_low) * (x / 100. - F_low) / (F_high - F_low) \
            + psi_low
        return 2. ** psi_x

    def UpdateDg(self):
        """
        computes and stores geometric mean grain size
        """
        psi_geometric_mean = np.sum(np.log(self.D) / log(2.) * self._F)
        
        self._Dg = 2. ** psi_geometric_mean
        
    def UpdateSigmaSG(self):
        """
        computes and stores geometric standard deviation of grain size on psi 
        scale -- and transforms back to length scale
        """
        var_psi = np.sum((np.log(self.D / self._Dg) / log(2.)) ** 2 * self._F) 
        
        self._SigmaSG = var_psi ** 0.5
        self._Sigmag = 2. ** self.SigmaSG

    def UpdateFfiner(self):
        """
        Computes cumulative fraction finer given the grain size fractions in 
        the GSD.  Renormalizes if necessary to ensure that the sum of fractions
        across all bins is 1.
        """
        Cumulative_sum = np.sum(self._F) 
        self._FfinerThanD_upper = np.cumsum(self._F) 
        if Cumulative_sum != 1.:
            self._F /= Cumulative_sum
            self._FfinerThanD_upper /= Cumulative_sum

    def __init__(self, BinBdySizes):
        """
        This creates an array starting at index 0 representing all sediment
        sizes one size class finer than on the bed (index 0).
        Array BinBdySizes should start at index 0 to characterize boundaries
        of washload bin
        
        Arguments:
            BinBdySizes -- [float]
        """
        if not self.Initialized:
            #self.UpToDate = True
            NBedSizes = len(BinBdySizes) - 2
            self._NBedSizes = NBedSizes
            self._F = np.zeros(NBedSizes + 1)
            self.D = np.zeros(NBedSizes + 1)
            self.Vs = np.zeros(NBedSizes + 1)
            self.D_lower = np.zeros(NBedSizes + 1)
            self.D_upper = np.zeros(NBedSizes + 1)
            for k in range(NBedSizes + 1):
                self.D_lower[k] = BinBdySizes[k]
                self.D_upper[k] = BinBdySizes[k + 1]
                self.D[k] = (self.D_upper[k] * self.D_lower[k]) ** 0.5
                self.Vs[k] = self.DietrichSettlingVelocity(self.D[k], \
                    1. / 1300000., 9.81, 1000., 2650.)
    
            self._FfinerThanD_upper = np.zeros(NBedSizes + 1)
            self.Initialized = True
            
        else:
            raise RuntimeError('Tried to initiate clsGSD twice.')

    def UpdateStatistics(self):
        """
        Populates basic statistics of the size distribution.
        """
        self.UpdateFfiner()
        self.UpdateDg()
        self._D50 = self.D_x(50.)
        self._D65 = self.D_x(65.)
        self._D84 = self.D_x(84.)
        self._D90 = self.D_x(90.)
        self._D95 = self.D_x(95.)
        self.UpToDate = True
        self.UpdateSigmaSG()
    
    def DietrichSettlingVelocity(self, D, nu, g, rho, rhos):
        """
        Returns settling velocity in m/s
        
        Arguments:
            D -- float (particle diameter in mm)
            nu -- float (kinematic viscosity (m2/s))
            g -- float (gravity (m/s2))
            rho -- float (water density (kg/m3))
            rhos -- float (sediment density (kg/m3))
        """
        Dstar = (rhos - rho) * g * (D / 1000) ** 3 / (rho * nu ** 2)
        if Dstar ** 2 > 0.05:
            a1 = 1.92944 * log10(Dstar)
            a2 = -0.09815 * log10(Dstar) ** 2
            a3 = -0.00575 * log10(Dstar) ** 3
            a4 = 0.00056 * log10(Dstar) ** 4
            Wstar = 10 ** (-3.76715 + a1 + a2 + a3 + a4)
        else:
            Wstar = Dstar ** 2 / 5832
            
        return (Wstar * (rhos - rho) * g * nu / rho) ** (1 / 3)