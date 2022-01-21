import numpy as np
from math import log10, log

class clsGSD(object):
    """
    Defines a grain size distribution.
    
    Stores grain size boundaries, representative sizes, and grain size
    fractions in each bin of the distribution.  Provides methods for
    computing various parameters and moments.
    
    Parameters
    ----------
    BinBdySizes : array_like(float)
        Sediment grain size at each bin boundary.
    
    Attributes
    ----------
    NBedSizes : int
        Number of bed material size classes. Note that the distribution
        includes the possibility of storing sediment in a size finer
        than the finest bed material size. 
    F : array_like(float, length = NBedSizes + 1)
        Fraction of distribution in bin.
    D : array_like(float, length = NBedSizes + 1)
        Representative diameter (mm) for bin.  Uses geometric mean
        of bin boundary sizes.
    Vs : array_like(float, length = NBedSizes + 1)
        Settling velocity for representative size (m/s).
    D_lower : array_like(float, length = NBedSizes + 1)
        Diameter (mm) at lower bin boundary.
    D_upper : array_like(float, length = NBedSizes + 1)
        Diameter (mm) at upper bin boundary.
    FfinerThanD_upper : array_like(float, length = NBedSizes + 1), read only
        Fraction of distribution finer than the upper boundary
        for the bin.
    D95 : float, read only
        Diameter (mm) for which 95% is finer.
    D90 : float, read only
        Diameter (mm) for which 90% is finer.
    D84 : float, read only
        Diameter (mm) for which 84% is finer.
    D65 : float, read only
        Diameter (mm) for which 65% is finer.
    D50 : float, read only
        Diameter (mm) for which 50% is finer.
    Dg : float, read only
        Geometric mean diameter (mm).
    Sigmag : float, read only
        Geometric standard deviation (mm).
    SigmaSG : float, read only
        Standard deviation on psi scale (no units).
    Initialized : bool
        Flag to indicate if class is initialized.
    UpToDate : bool
        Flag that determines whether statistics are updated and thus
        do not need to be recomputed.  Set to false if F or D are changed.
    
    Notes
    -----
    POTENTIAL CODE UPDATE: In code profiling, it appears as though
    a significant fraction of MAST-1D's computational time 
    is spent in the clsGSD class. There may be opportunities to
    optimize the time used by the code.  The present approach assumes
    that certain statistics from a distribution will be needed many times,
    so it would save time to simply compute the statistic once and stores
    in memory, recalling the previosuly computed value when needed, 
    and not re-computing each time it is needed.  This is implemented
    by setting a status flag whenever the bin 
    probablities are changed. The flag indicates that statistics are 
    not up to date. Then, if a statistic is needed, all common statistics are
    updated, and the flag tracking whether the statistics are up to date
    is set to True.  Then, if another statistic is called, it does not 
    need to be recomputed until the size distribution changes again.  
    The hope was that this would ultimately save time, but it is possible
    that it would be faster simply to compute statistics on an as-needed basis, 
    for example by using the D(x) method. It may also be possible
    to include some C code snippets to speed some of the computations.
    Optimizing the way size distribution statistics are computed in
    this class may reprsent the largest opportunity for reducing
    execution time for MAST-1D.
    
    """
    Initialized = False
    
    def getF(self):
        """
        Fraction in bin.
        """
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
        Cumulative fraction finer than upper bin boundary.
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
        D95 (mm).
        """
        if not self.UpToDate: self.UpdateStatistics() # recompute cumulative 
        # totals and statics if necessary
        return self._D95
        
    D95 = property(getD95) 
    
    def getD90(self):
        """
        D90 (mm).
        """
        if not self.UpToDate: self.UpdateStatistics() # recompute cumulative 
        # totals and statics if necessary
        return self._D90
    
    D90 = property(getD90)    
    
    def getD84(self):
        """
        D84 (mm).
        """
        if not self.UpToDate: self.UpdateStatistics() # recompute cumulative 
        # totals and statics if necessary
        return self._D84
    
    D84 = property(getD84)
    
    def getD65(self):
        """
        D65 (mm).
        """
        if not self.UpToDate: self.UpdateStatistics() # recompute cumulative 
        # totals and statics if necessary
        return self._D65
    
    D65 = property(getD65)
    
    def getD50(self):
        """
        Median diameter (mm).
        """
        if not self.UpToDate: self.UpdateStatistics() # recompute cumulative 
        # totals and statics if necessary
        return self._D50
    
    D50 = property(getD50)
    
    def getDg(self):
        """
        Geometric mean diameter (mm).
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
        if not self.UpToDate: self.UpdateStatistics() # recompute cumulative 
        # totals and statics if necessary
        return self._SigmaSG
    
    SigmaSG = property(getSigmaSG)
    
    def D_x(self, x):
        """
        Compute sediment diameter finer than a given percentile x.
        
        Geometrically interpolates a percentile of a grain size 
        distribution. Performs interpolation in log space (so in psi space).
        
        Parameters
        ----------
        x : float
            Percentile to compute.
        
        Returns
        -------
        float
            Percentile
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
        Computes and stores geometric mean grain size.
        """
        psi_geometric_mean = np.sum(np.log(self.D) / log(2.) * self._F)
        
        self._Dg = 2. ** psi_geometric_mean
        
    def UpdateSigmaSG(self):
        """
        Computes and stores geometric standard deviation of grain size on psi 
        scale -- and transforms back to length scale and stores that.
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
        
        The method updates cumulative fractions finer, geometric
        mean size, geometric standard deviation, and
        D50, D65, D84, D90, and D95.  Once complete, sets a flag
        indicating that statistics are OK and can be used without
        recomputation.
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
        Computes particle settling velocity.
        
        Parameters
        ----------
        D : float
            Particle diameter (mm).
        nu : float
            Kinematic viscosity (m2/s).
        g : float
            Gravitational constant (m/s2)
        rho : float
            Water density (kg/m3).
        rhos : float
            Sediment density (kg/m3).
        
        Returns
        -------
        Settling Velocity (m/s).
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