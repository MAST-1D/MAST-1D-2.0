
from clsExchangeTypes import clsExchangeTypes
from clsGSD import clsGSD
import numpy as np
from math import exp
#from decimal import* #Katie add
#getcontext().prec = 6
#from numba.decorators import jit, autojit
class clsLoad(object):
    """
    This class defines the sediment load and feed for a given node in all sizes
    and provides methods for computing bed material load.  It also provides 
    methods conserving mass and tracer in suspended washload sediment for a 
    given flow (bin of duration curve) and for averaging washload tracer 
    concentration across the flow duration curve. Right now, only Wilcock and 
    Crowe is implemented for sediment transport, but new features can be added 
    in the future as separate methods that modify the basic data elements in 
    this class, the Qsjk(j,k) and TMudj(j,k) arrays.

    This class does not provide rubust methods for tracking bed material tracer 
    concentration in the load.  Instead, in class node, it is assumed that bed 
    material tracer concentration in the load is the same as in the bed.  A 
    variable is provided in load to track duration-averaged tracer 
    concentration.

    In principle, it would be possible to track tracer concentration in the 
    load in each size class and to allow this to differ from tracer 
    concentration in the active layer.  This would require the specification of 
    an exchange rate between the active layer and the load, which could be 
    determined for a given bin of the duration curve as a function of near bed 
    concentration, and settling velocity (or entrainment rate could be solved 
    as a function of particle step length--a function of hydraulics--and 
    computed load).  For now, it is assumed that all bed material classes 
    undergo step lengths that are short relative the spacing between nodes, so
    that a large amount of mixing occurs between the active layer and load over 
    the distance dxc and very little bypassing occurs.  For this reason, tracer 
    concentration in the load is only tracked for mud. Elsewhere in the code, 
    tracer concentration in the load is set to that of the active layer.

    Not sure that ExSED and ExTracer are needed since mass conservation due to 
    morphodynamic movement occurs at node level.
    
    Attributes:
        TMudj -- [float] (Tracer concentration in flow j for size k=0 for 
                 tracer type l)
        TMudFeedj -- [float] (Tracer concentration of tracer type l in size k=0 
                     in feed for bin j of duration curve)
        TBedFeedk -- [float] (Tracer concentration of tracer type l in size k.  
                     Should be set equal to T of upstream active layer or feed 
                     for reach.)
        ExSed -- clsExchangeTypes (Size specific sediment fluxes into and out 
                 of reservoir (units of volume/time))
        ExTracer -- [clsExchangeTypes] (Size specific tracer concentrations for
                    sediment moving into and out of reservoir)
        LatSed -- [float] (Volumetric lateral sediment source for reservoir 
                  (units of volume/time))
        LatTracer -- [float] (Concentration of tracers in lateral sediment 
                     source)
        SinkSed -- [float] (Volumetric sediment sink term.)
        NTracers -- int (number of tracers)
        NSizes -- int (number of bed material grain sizes.  note that k = 0 
                  represents washload.)
        NFlows -- int (number of flows in flow duration curve)
        GSDBedloadj -- [clsGSD] (grain size distribution of load for bin j of
                      FDC)
        GSDBedloadAv -- clsGSD (grain size distribution for annual average 
                        load)
        Qsjk -- [float] (Sediment load in bin j and size k)
        Qsjkfeed -- [float] (Sediment load in bin j and size k)
        QsjTot -- [float] (Sediment load in bin j of flow duration curve, 
                  summed across all grain sizesL)
        QsAvkLoad -- [float] (Sediment load capacity averaged across all bins 
                     of FDC in size k.)
        QsAvkFeed -- [float] (Sediment feed averaged across all bins of FDC in 
                     size k)
        TVolumeAvMudFeed -- [float] (Sediment feed tracer concentration in size
                            k=0 and tracer l averaged across all volume 
                            (including weighting by duration curve))
        TTemporalAvMudFeed -- [float] (Sediment feed tracer concentration in 
                              size k=0 and tracer l averaged across all time in
                              duration curve)
        TVolumeAvMudLoad -- [float] (Sediment load tracer concentration in size
                            k=0 and tracer l averaged across all volume 
                            (including weighting by duration curve))
        TTemporalAvMudLoad -- [float] (Sediment load tracer concentration in 
                              size k=0 and tracer l averaged across all time in 
                              duration curve)
        TDepositionAvMudFeed -- [float] (Sediment feed tracer concentration in 
                                size k = 0 and tracer l averaged across all 
                                overbank mud deposition and duration curve)
        TDepositionAvMudLoad -- [float] (Sediment load tracer concentration in 
                                size k = 0 and tracer l averaged across all 
                                overbank mud deposition and duration curve)
        TauPrime -- float (channel-averaged shear stress)
        TauRm -- float (critical shear stress for the median particle size)
        Initialized -- bool
    """
    Initialized = False
    
    def __getattribute__(self, name):
        if name == 'NTracers':
            return len(self.ExTracer)
        elif name == 'NSizes':
            return self.ExSed.NSizes
        elif name == 'QsavBedTot':
            return np.sum(self.QsAvkLoad[1:])
        elif name == 'QsavBedTotFeed':
            return np.sum(self.QsAvkFeed[1:])
        elif name == 'QsavTotAllFeed':
            return np.sum(self.QsAvkFeed[0:])
        elif name == 'QsavTot':
            return np.sum(self.QsAvkLoad)
        else:
            return object.__getattribute__(self, name)
    
    def __init__(self, NFlows, BinBdySizes, NTracers):
        """
        Arguments:
            NFlows -- float
            BinBdySizes -- [float]
            NTracers -- int
        """
        if not self.Initialized:
            NSizes = len(BinBdySizes) - 2
            self.NSizes = NSizes
            self.NTracers = NTracers
            self.NFlows = NFlows
            
            self.ExSed = clsExchangeTypes(NSizes)
            self.ExTracer = [clsExchangeTypes(NSizes) for i in \
                            range(NTracers)]
            self.LatSed = np.zeros((NFlows, NSizes + 1))
            self.LatTracer = np.zeros((NFlows, NSizes + 1, NTracers))
            self.SinkSed = np.zeros((NFlows, NSizes + 1))
            
            self.GSDBedloadAv = clsGSD(BinBdySizes)
            self.GSDBedloadj = [clsGSD(BinBdySizes) for i in range(NFlows)]
            
            self.TMudj = np.zeros((NFlows, NTracers))
            self.TMudFeedj = np.zeros((NFlows, NTracers))
            self.TBedFeedk = np.zeros((NSizes + 1, NTracers))
            
            self.Qsjk = np.zeros((NFlows, NSizes + 1))
            self.Qsjkfeed = np.zeros((NFlows, NSizes + 1))
            self.QsjTot = np.zeros(NFlows)
            self.QsAvkLoad = np.zeros(NSizes + 1)
            self.QsAvkFeed = np.zeros(NSizes + 1)
            
            self.TVolumeAvMudFeed = np.zeros(NTracers)
            self.TTemporalAvMudFeed = np.zeros(NTracers)
            self.TVolumeAvMudLoad = np.zeros(NTracers)
            self.TTemporalAvMudLoad = np.zeros(NTracers)
            self.TDepositionAvMudFeed = np.zeros(NTracers)
            self.TDepostiionAvMudLoad = np.zeros(NTracers)
            
            # Katie add shear stresses for use in width function.  They would make
            # more sense in the Node class but cannot be modified from the WC
            # function because it returns an array instead of modifying the Node
            # class directly.
            self.TauPrime = np.zeros(NFlows)
            self.TauRm = 0.
            self.Eshear = np.zeros(NSizes + 1)#Katie add--size-specific excess shear stress
            self.MudEnt = 0. # Katie add            
            self.Initialized = True
        else:
            raise RuntimeError('Tried to initiate clsLoad twice.')

    def WilcockCroweLoadBySize(self, U, Sf, GSD, rho_w, SG, Bc, TrinityFit, \
        CalibrationFactor, IndD = 0): # Katie add DRange so single sizes can be calculated for width func.
        """
        This function returns a 1d array of size specific loads of length
        GSD.nsizes (i.e. the bed material sizes in GSD).  Load is in units of
        m3/s.
        It uses the reference stress form of Wilcock and Crowe 2003 unless the
        flag TrinityFit is true, in which case it uses the reference stress of
        Guiman 2009.
        This should probably be updated to allow any calibrated reference 
        stress to be used.  For now, it allows for a calibration factor that 
        multiplies the reference stress from either of the above equations by a
        constant.
        
        Arguments:
            U -- float
            Sf -- float
            GSD -- clsGSD
            rho_w -- float
            SG -- float
            Bc -- float
            TrinityFit -- boolean
            CalibrationFactor -- float
        
        Return: [float]
        """      
        #GSD.UpToDate = True
        g = 9.81
        ResultArray = np.zeros(GSD.NBedSizes + 1) 
        #x = GSD.D65
        #y = 'filler'
        
        TauPrime = self.findWC_TauPrime(g, Sf, U, GSD.D65) # Katie change; moved TauPrime to its own function
        
        ustar = (TauPrime / rho_w) ** 0.5
        # ************find sand fraction *****************
        
        
        # ***************end find sand fraction************
        
        # This is where the form of the reference stress equation is specified
        if TrinityFit == False:
            i = 1
            Fs = 0.
            while GSD.D_upper[i] <= 2.: # Katie change from .002 (which is the lower boundary)
                Fs = GSD.FfinerThanD_upper[i] # Katie change from parentheses to brackets
                i = i + 1
            TauStarRm = 0.021 + 0.015 * exp(-20. * Fs) # This is the Wilcock
            # and Crowe 2003 form
        else:
            TauStarRm = 0.03 + (0.052 - 0.03) / (1. + exp(7.1 * (GSD.SigmaSG \
                - 1.66))) # This is the Geuman 2009 form for Trinity River--Katie:  SG represents geometric SD; SigmaG does not appear to be used.
        
        # ****This is where calibration is  performed on reference shear ******
        TauStarRm = TauStarRm * CalibrationFactor

        # ********************end calibration******************************
        
        TauRm = TauStarRm * (SG - 1.) * rho_w * g * GSD.Dg / 1000.
        self.TauRm = TauRm
        
        DRange = range(1, GSD.NBedSizes + 1)
        if IndD != 0:
            DRange =  DRange[IndD:]
            for i in DRange:
                ResultArray[i], filler = self.calc_wStari(TrinityFit, TauRm, GSD.D[i], GSD.F[i], GSD, TauPrime, ustar, Bc, SG, g)
        else:            
            for i in DRange:
                ResultArray[i], self.Eshear[i] = self.calc_wStari(TrinityFit, TauRm, GSD.D[i], GSD.F[i], GSD, TauPrime, ustar, Bc, SG, g)

        return ResultArray
        
    def findWC_TauPrime(self, g, Sf, U, D65):
        
        """
        Katie add.  Calculates the channel-averaged shear stress.  Was moved from
        the main Wilcock-Crowe function.
        
        Arguments:
            Sf -- float
            g -- float
            D65
            U -- float
        """
        TauPrime = 1000. * g * (0.013) ** 1.5 * (Sf * 2 * D65) ** 0.25 * \
            U ** 1.5 # eq 2.14 in BAGS primer: n = 0.0146D65 for D65 in mm
            
        return TauPrime
        
    def calc_wStari(self, TrinityFit, TauRm, D, F, GSD, TauPrime, ustar, Bc, SG, g):
        TauRi = self.findWC_ReferenceShear(TrinityFit, TauRm, D, GSD.Dg, GSD.D50) # Katie add; moved TauRi to different function
        Eshear = self.return_excess_shear(TauPrime, TauRi) # Katie add
        phi = TauPrime / TauRi
        if phi < 1.35:
            wStari = 0.002 * phi ** 7.5
        else:
            wStari = 14. * (1. - 0.894 / phi ** 0.5) ** 4.5
            
        # Reference transport rate
#            if wStari < 0.002 and i == GSD.NBedSizes: # Katie add reference transport rate from WC
#                ResultArray[i] = 0.
#            else:
#                ResultArray[i] = wStari * GSD.F[i] * ustar ** 3 * Bc / (SG - 1.) \
#                / g
        qs = wStari * F * ustar ** 3 * Bc / (SG - 1.) \
            / g
            
        return qs, Eshear
    
    def return_excess_shear(self, TauPrime, TauRi): # Katie add
        """
        Calculates and returns excess shear stress for a specific size class for
        plotting purposes
        """
        #print TauPrime/TauRi
        return TauPrime/TauRi 
        
    def findWC_ReferenceShear(self, TrinityFit, TauRm, D, Dg, D50):
        """
        Katie add.  This function calculates the 
        the reference shear stress (TauRi,equations 3 and 4 in Wilcock and Crowe,
        2003) given inputs of the size in question, the D50, and the critical shear
        stress for the mean bed size.  It was made independent from the Wilcock 
        and Crowe function to make it easier to calculate the critical shear stress 
        for the D84 for the width change function.
        
        Arguments:
            TauRm -- float
            D -- float
            D50 -- float
        """
        b = 0.
        if TrinityFit == True:    
            b = 0.7 / (1. + exp(1.9 - D / (3*Dg)))
        else: # Katie add--before, only WC form of the hiding function was used
            b = 0.67 / (1. + exp(1.5 - D / Dg))
        
        TauRi = TauRm * (D / D50) ** b
        
        return TauRi

#############################################################################
    def WrightParkerLoadBySize(self, Q, H, U, Sf, GSD, rho_w, SG, Bc):
        
        """
        This function returns a 1d array of size specific loads of length
        GSD.nsizes (i.e. the bed material sizes in GSD).  Load is in units of
        m3/s.
        
        Load is calculated using the relation of Wright and Parker (2004).
        Currently, no calibration factor is included.
        
        Arguments:
            H -- float
            U -- float
            Sf -- float
            GSD -- clsGSD
            rho_w -- float
            SG -- float
            Bc -- float
        """

        ResultArray = np.zeros(GSD.NBedSizes + 1)    
        CArray = np.zeros(GSD.NBedSizes + 1)
        
        q = Q/Bc                
        g = 9.81
        qStar = q/((g*(GSD.D50/1000.))**0.5*(GSD.D50/1000.))
        kappa = 0.4 # von Karman constant      
        hsk = H
        uStarsk = 0
        WPalpha = 0
        vwater = 1.307*10**(-6) # Kinematic viscosity of water, 10 deg C
        kc = 0        
        error = 1
        Count = 0
        

        # ************ Guess depth related to skin friction *****************
        
        while abs(error) >= 1e-03 and Count < 100:        
        
            D50 = GSD.D50/1000.
            D84 = GSD.D84/1000.
            R = (hsk*Bc)/(2*hsk+Bc) # Hydraulic radius
            TauStar = (R*Sf)/((SG-1)*D50) # Shields number for mixture (does not assume wide channel)
            Fr = U/(hsk*g)**0.5 # Froude number
            uStarsk = (g*hsk*Sf)**0.5 # Shear velocity due to skin friction
            
            TauStarSkin = 0.05 + 0.7*(TauStar*Fr**0.7)**0.8 # Wright and Parker eq. 17
            kc = (3*D84)*(TauStar/TauStarSkin)**4   # Wright and Parker eq. 22. Sand grain roughness height is 3*D84     
            
            C5t = 0
            
            #  Grainsize specific entrainment
            
            for i in range(1, GSD.NBedSizes + 1):
                
                D = GSD.D[i]/1000. 
                F = GSD.F[i]

                vsi = (g*(SG-1)*(D)**2)/(18*vwater) # Stokes Law UNITS
                Rpi = ((SG-1)*g*D)**0.5*D/vwater # particle Reynolds number
                Xi = (uStarsk*(Rpi**0.6)/vsi)*(Sf**0.08)*(D/D50)**0.2 # Wright and Parker eq. 19
                Esi = ((7.8*10**-7)*((1-0.28*GSD.SigmaSG)*Xi)**5)\
                    /(1+(7.8*10**-7*(((1-0.28*GSD.SigmaSG)*Xi)**5)/0.3)) # Wright and Parker entrainment function, eq. 18

                C5i = (Esi*F) # Divided by in paper, but I think it's wrong...
                C5t = C5t + C5i # Sediment concentration at 5% of flow depth from bed
                CArray[i] = C5i
                
            # Determination of new depth 
            if C5t/Sf <= 10: # Wright and Parker eq. 9
                WPalpha = 1.-0.06*(C5t/Sf)**0.77
            else:
                WPalpha = 0.67-0.0025*(C5t/Sf)
            
            newhsk = D50*((WPalpha*qStar*(kc/D50)**(1./6.))/(8.32*Sf**0.5))**(3./5.) # Wright and Parker eq. 21            
            error = (hsk-newhsk)
            Count = Count + 1
            hsk = hsk - (error/3.)
    
        if Count > 99:
            print ('New skin-depth did not converge in Wright and Parker calculation')

        # ************ Calculate depth-integrated transport based on concentration *****************  
        
        # Size-specific transport
        
        for i in range(1, GSD.NBedSizes + 1):
            D = GSD.D[i]/1000.            
            F = GSD.F[i]
            C5i = CArray[i]
            
        # Estimation of depth-integral
            
            # Settling velocity--from Jimenez and Madsen, 2003 (Equations 2-48a and b in Garcia (ed) Sedimentation Engineering)
            Sstar = (D/0.9)*(g*(SG-1)*(D/0.9))**0.5/(4*vwater)
            vsi = (g*(SG-1)*(D/0.9))**0.5*(0.954+(5.12/Sstar))**-1
            
            ZRi = vsi/(WPalpha*kappa*uStarsk) # Stratification-adjusted Rouse number
            
            if ZRi <= 1: #  This is Wright and Parker's approximation, eq. 26; negates need to discretize over depth
                I = 0.679*exp(-2.23*ZRi)
            else:
                I = 0.073*ZRi**(-1.44)
        
        # Final size-specific load calculation
        
            qsi = (9.70*uStarsk*hsk*C5i*I*(hsk/kc)**(1./6.)/WPalpha)*Bc # Transport for size i, Wright and Parker eq. 24, multiplied by channel width
            ResultArray[i] = qsi
            
        return ResultArray

#############################################################################

    def UpdateSedimentLoadByDurationAndSize(self, DC, SurfaceGSD, rhow, \
        SG, Bc, FractionAlluvial, TransFunc, TrinityFit, CalibrationFactor):
        """
        Populates sediment load array in load object with wilcock crowe loads
        (m3/s) for each bin of FDC and for each size in bed GSD.
        Requires that hydraulics of node have already been populated
        
        Arguments:
            DC -- clsDurationCurve
            SurfaceGSD -- clsGSD
            rhow -- float
            SG -- float
            Bc -- float
            FractionAlluvial -- float
            TrinityFit -- bool
            CalibrationFactor -- float
        """
        for j in range(DC.NFlows()):
            
            Qk =   np.zeros(SurfaceGSD.NBedSizes + 1)          
            
            #  Calls selected transport function
            if TransFunc == 'WilcockCrowe':
                Qk = self.WilcockCroweLoadBySize(DC.Uc[j], DC.Sf[j], SurfaceGSD, \
                    rhow, SG, Bc, TrinityFit, CalibrationFactor)

            elif TransFunc == 'WrightParker':
                Qk = self.WrightParkerLoadBySize(DC.Qwc[j], DC.Hc[j], DC.Uc[j], DC.Sf[j], \
                    SurfaceGSD, rhow, SG, Bc)
            else:
                print('Transport function is not supported')
            
            for k in range(1, self.NSizes + 1):
                if Qk[k] < 0.:
                    print('Load < 0')
                # Adjust for fraction of bed that is alluvial
                Qk[k] *= FractionAlluvial

            self.Qsjk[j, 1:SurfaceGSD.NBedSizes + 1] = \
                Qk[1:SurfaceGSD.NBedSizes + 1]
            self.QsjTot[j] = np.sum(Qk)
        
        # This loop populates the size distribution for the bedload in each
        # bin of the FDC
        for j in range(DC.NFlows()):
            for k in range(1, SurfaceGSD.NBedSizes + 1):
                self.GSDBedloadj[j].F[k] = self.Qsjk[j, k] / self.QsjTot[j]
        
        # This loop populates total annual bed material load and average load
        # for each size class,

        for k in range(1, SurfaceGSD.NBedSizes + 1): # start at index k = 1 
            # since k = 0 represent washload
            self.QsAvkLoad[k] = 0.
            for j in range(DC.NFlows()):
                self.QsAvkLoad[k] = self.QsAvkLoad[k] + self.Qsjk[j, k] * \
                    DC.p[j]
                if self.QsAvkLoad[k] < 0.:
                    print('load < 0 for k = ' + str(k)) 
            self.QsAvkLoad[k] = float(self.QsAvkLoad[k])

        # This loop updates GSD for annual bed material load
        # start at index k = 1 since k = 0 represents washload
        self.GSDBedloadAv.F[1:SurfaceGSD.NBedSizes + 1] = \
            self.QsAvkLoad[1:SurfaceGSD.NBedSizes + 1] / self.QsavBedTot
    
    def CurranLoad(self, Q, Section='Middle'):
        
        """
        Katie add!
        
        Calculates suspended sediment load for feed based on Curran's rating curve 
        calibrated to the reservoir sediment presented in Konrad's pre-removal
        modeling paper.  Includes the discharge correction that needs to be applied
        when using the gauge at McDonald Bridge.
        
        Arguments:
            Q -- float
            Section -- str (determines whether or not to normalize the discharge because feed
                is determined by the Upper River discharge)
        """
        NormQ = 0.
        if Section == 'Middle':
            NormQ = .78*Q
        elif Section == 'Upper':
            NormQ = Q
        else:
            print('Unsupported discharge for Curran suspended load')

        Qssc = (0.1*10**(-4)*(NormQ)**2.5)

        return Qssc/2700. # Converts from kg/s to m^3/s
        
    def UpdateFeedDurationCurve(self, ControlNode, Multiplier):
        
        """
        Standard method for updating feed--designed for the flow duration curve.  
        Feed is a user-defined fraction of sediment transport capacity 
        
        Attributes:
        
        -ControlNode--clsNode (Node object that holds sediment feed)
        -Multiplier--float (Multiplier for feed)
        """

        KFeed = np.zeros(self.NSizes + 1)
        JKFeed = np.zeros((ControlNode.DC.NFlows(), self.NSizes + 1))
        for k in range(0, self.NSizes + 1):               
            KFeed[k] = (ControlNode.Load.QsAvkFeed[k])*Multiplier
            for j in range(ControlNode.DC.NFlows()):
                JKFeed[j, k] = Multiplier * \
                    ControlNode.Load.Qsjkfeed[j, k]            
        return KFeed, JKFeed
        
    def UpdateFeedRatingCurve(self, ControlNode, Multiplier, MudFraction, Q,\
            vfunc, TrinityFit, CalibrationFactor, ControlGSD):
                
        """
        Rating curve method--the standard for hydrograph runs.
        Feed is a user-defined fraction of sediment transport capacity
        for a given flow.  Mud feed is a set proportion of feed for the
        next finest size class.  Returns size specific, and flow/size specific
        feed, although the function can only be used when there is only one flow
        in the duration curve (i.e. when a hydrograph is used).
        
        Attributes:
        
        -ControlNode--clsNode (Node object that holds sediment feed)
        -Multiplier--float (Multiplier for feed)
        -MudFraction--float (ratio of mud to next finest size class)
        -Q--float (discharge)
        -vfunc--bool (True for Manning, False for Chezy)
        -TrinityFit--bool (True for Trinity River sed transport calc, False for 
            Wilcock and Crowe)
        -CalibrationFactor-float (calibration for sediment transport equation)
        """
        
        KFeed = np.zeros(self.NSizes + 1)
        JKFeed = np.zeros((1, self.NSizes + 1))
        
        # Calculate hydraulics for Q on control node assuming normal flow
        NormFlow = ControlNode.NormalChannelDepthAndDischarge(Q,.01, vfunc)                            
        U = NormFlow[1]/(NormFlow[0]*ControlNode.Bc)
        
        # Calculate sediment transport capacity for Q
        FeedArray = ControlNode.Load.WilcockCroweLoadBySize(U, ControlNode.Slope,\
            ControlGSD, 1000., 2.7, ControlNode.Bc, TrinityFit,\
            CalibrationFactor)
        
        for k in range(1, ControlNode.NSizes + 1):                   
            KFeed[k] = FeedArray[k]*Multiplier
            JKFeed[0, k] = Multiplier * \
            FeedArray[k] # Only works when there is one flow bin.
                    
        # Calculate mud load as a set constant of next finest size class
        KFeed[0] = KFeed[1]*MudFraction
        JKFeed[0, 0] = KFeed[0] # Only works when there is one flow bin.

        return KFeed, JKFeed
        
    def UpdateElwhaFeedRatingCurve(self, ControlNode, Multiplier, MudFraction, Q,\
            vfunc, TrinityFit, CalibrationFactor, ControlGSD, Section = 'Middle', \
            Removal = False, PercSus = 0, FixedCapacity = False, FirstPulse = False):
                
        """
        Uses the UpdateFeedRatingCurve in this object to calculate the bed material
        feed, then uses Curran's emperical rating curve for the mudload.
        
        Section determines whether the discharge needs to be normalized for the Curran
        relation.  See UpdateFeedRatingCurve for other attributes.
        """
                
        KFeed, JKFeed = self.UpdateFeedRatingCurve(ControlNode, Multiplier, MudFraction, Q,\
            vfunc, TrinityFit, CalibrationFactor, ControlGSD)

        if FixedCapacity == True:
            KFeed = ControlNode.Load.QsAvkFeed*Multiplier
            JKFeed[0,:] = KFeed # Only works with one flow bin

        Qssc = 0.
        if Removal == False:
            Qssc = self.CurranLoad(Q, Section)*Multiplier
        elif FirstPulse == True: # Only sand/fine gravel in load
            TotFeed = sum(KFeed[1:7])/(1-PercSus)
            #KFeed[1]=0.
            #JKFeed[0,1]=0.
            for k in range(7,len(KFeed)):
                KFeed[k] = 0.
                JKFeed[0,k] =0.
            Qssc = TotFeed*PercSus
        else:
            TotFeed = sum(KFeed[1:])/(1-PercSus) # Reservoir 'suspended' ratio
            Qssc = TotFeed*PercSus # PercSus is the percentage of total load that is suspended.
            #Qssc = self.CurranLoad(Q, Section)*Multiplier
            
        KFeed[0] = Qssc
        JKFeed[0, 0] = Qssc # Only works when there is one flow bin.

        return KFeed, JKFeed
		
    def UpdateElwhaFeedFineRatingCurve(self, ControlNode, Multiplier, MudFraction, Q,\
            vfunc, TrinityFit, CalibrationFactor, ControlGSD, FineGSD, Section = 'Middle', \
            Removal = False, PercSus = 0, FixedCapacity = False,  FirstPulse = True):

        """
        Like UpdateElwhaRatingCurve but for alternate removal scenario.  All sizes
		are supplied at capacity and an additional pulse consists of just fine sediment
        """
        KFeed, JKFeed = self.UpdateFeedRatingCurve(ControlNode, 1, MudFraction, Q,\
            vfunc, TrinityFit, CalibrationFactor, ControlGSD)
			
        KFeedFine, JKFeedFine = self.UpdateFeedRatingCurve(ControlNode, Multiplier-1, MudFraction, Q,\
            vfunc, TrinityFit, CalibrationFactor, FineGSD)
			
        for k in range(7,len(KFeedFine)):
            KFeedFine[k] = 0.
            JKFeedFine[0,k] =0.
			
        KFeed = KFeed + KFeedFine
        JKFeed = JKFeed + JKFeedFine
		
        TotFeed = sum(KFeedFine[1:])/(1-PercSus) # Reservoir 'suspended' ratio
        Qssc2 = TotFeed*PercSus # PercSus is the percentage of total load that is suspended.
        Qssc1 = self.CurranLoad(Q, Section)*Multiplier
        Qssc = Qssc1+Qssc2

        KFeed[0] = Qssc
        JKFeed[0, 0] = Qssc # Only works when there is one flow bin.

        return KFeed, JKFeed		
        
        
    def UpdateElwhaDamRemovalDurationCurve(self, ControlNode, Multiplier):
        
        """
        An exponential decay function to estimate sediment exported from Lake
        Mills.
        
        Attributes:
        -ControlNode-clsNode (Node with original feed characteristics)
        -C--float (Decay multiplier)
        -t--float (Time since dam removal)
        -tau--float (Exponent that determines time to recovery--same unit as t)

        Currently obsolete function.
        """
        
        KFeed = np.zeros(self.NSizes + 1)
        KFeedFine = np.zeros(self.NSizes + 1)
        JKFeed = np.zeros((ControlNode.DC.NFlows(), self.NSizes + 1))
        JKFeedFine = np.zeros((ControlNode.DC.NFlows(), self.NSizes + 1))
        
        for k in range(0, ControlNode.NSizes + 1):                   
            KFeed[k] = ControlNode.Load.QsAvkFeed[k]
            KFeedFine[k] = ControlNode.Load.QsAvkFeed[k]*Multiplier
            for j in range(ControlNode.DC.NFlows()):
                JKFeed[j, k] = ControlNode.Load.Qsjkfeed[j, k]
                JKFeedFine[j, k] = ControlNode.Load.Qsjkfeed[j, k]*Multiplier
                
        for k in range(7,len(KFeed)):
            KFeedFine[k] = 0.
            for j in range(ControlNode.DC.NFlows()):
                JKFeedFine[j,k] =0.
                
        # Note:  this method does not use the percentage-based mudload for the removal like the rating curves
        
        KFeed = KFeed + KFeedFine
        JKFeed = JKFeed + JKFeedFine
        
        return KFeed, JKFeed
        
    def UpdateElwhaDamRemovalRatingCurve(self, ControlNode, C, t, tau,
        MudFraction, Q, vfunc, TrinityFit, CalibrationFactor):
        
        """
        Applies a multiplier to the sediment rating curve for the Elwha that decays
        exponentially.
        """        
        
        KFeed, JKFeed = self.UpdateElwhaFeedRatingCurve(ControlNode, 1,\
            MudFraction, Q, vfunc, TrinityFit, CalibrationFactor)

        # Mud feed is 18.5 times sand feed
        KFeed[0] = 18.5*sum(KFeed[1:3])
        JKFeed[0,0] = KFeed[0]
        
        #  Apply the multiplier
        for k in range(0, ControlNode.NSizes + 1):                   
            KFeed[k] = (KFeed[k])*(1 + C*exp(-t/tau))
            JKFeed[0, k] = (1 + C*exp(-t/tau)) * \
                JKFeed[0, k]
        
        return KFeed, JKFeed
    
    def QWashloadAnnual(self, DC, Qsjk):
        """
        Arguments:
            DC -- clsDurationCurve
            Qsjk -- float
        
        Return: [float]
        """
        return np.sum(self.Qsjk[:, 0] * DC.p)
    
    def ApplyMassConservationToMudLoad(self, DC):
        """
        Arguments:
            DC -- clsDurationCurve
        """
        self.QsAvkLoad[0] = 0.
        for j in range(DC.NFlows()):
            # note that sinksed includes overbank deposition
            # Katie add DC.p[j] to be multiplied to QsjkFeed[j]--then removed
            self.Qsjk[j, 0] = (DC.Migrationfactor[j] * \
                (self.ExSed.InMigration[0] - self.ExSed.OutMigration[0])) / \
                DC.p[j] + self.ExSed.InVerticalChange[0] - \
                self.ExSed.OutVerticalChange[0] + self.ExSed.InWidthChange[0] \
                - self.ExSed.OutWidthChange[0] + self.Qsjkfeed[j, 0] + \
                self.LatSed[j, 0] - self.SinkSed[j, 0]
            
            if self.Qsjk[j, 0] < 0.:
                error = - self.Qsjk[j, 0]
                self.Qsjk[j, 0] = 0.
#                print('Outgoing mud flux was computed to be ' + str(error) + \
#                    'm3/s less than zero in conversion equation for water \
#                    column for bin j = ' + str(j)) # Katie change 'conversation' to 'conversion' and comment out
            self.QsAvkLoad[0] = self.QsAvkLoad[0] + DC.p[j] * self.Qsjk[j, 0]#+self.MudEnt # Katie add--entrainment of mud from active layer

    def ApplyTracerConservationToMudLoad(self, DC):
        """
        WL May 5, 2013: Note that it should also be possible to apply tracer 
        concentration to the long-term average bedload as long as an exchange 
        rate is known between the active layer and the load.  This can be 
        specified as a function of suspended sediment transport rate and 
        settling velocity (which allows computation of near-bed concentration
        and thus entrainment rate) or as a function of particle step length 
        and the computed load.  In the latter case, the issue would be then 
        making particle step length a function of the hydraulics--I need to 
        review the literature.  As is, the active layer mixes with the load 
        only when net duration-averaged morphodynamic change occurs, so the 
        model is not yet a gravel tracer model (except to the extent that 
        lateral exchange is greater than vertical exchange, which is certainly
        possible given the right conditions). However, even the, the model 
        pulls the gravel that is being put into storage from both the active
        layer and load at rather arbitrary rates, so it is not clear whether
        it is pulling the correct amount of traver gravel from the active 
        layer. In other words, if there is minimal bed change, it is possible
        for tracer particles to move too far downstream since there is no 
        exchange between load and active layer.
    
        This routine sets Tjk in the load for k = 0
    
        This assumes that the sediment tracer concentration in the water column
        is always in steady state for a given flow
        This should be reasonable since the volume of suspended sediment in the
        water column at any given time is small relative to the exchange rates.
        
        Arguments:
            DC -- clsDuratinCurve
        """
        for L in range(self.NTracers):
            for j in range(DC.NFlows()):
                # note that both sediment outflux of washload sediment and 
                # tracer is computed using feed, not load.
            
                if self.Qsjk[j, 0] > 0.:

                    # At steady state, can use either influx or outflux for 
                    # this computation
                
                    TracerInflux = (DC.Migrationfactor[j] * \
                        self.ExSed.InMigration[0] * \
                        self.ExTracer[L].InMigration[0] + \
                        self.ExSed.InVerticalChange[0] * \
                        self.ExTracer[L].InVerticalChange[0] + \
                        self.ExSed.InWidthChange[0] * \
                        self.ExTracer[L].InWidthChange[0]) / DC.p[j] + \
                        self.LatSed[j, 0] * self.LatTracer[j, 0, L] + \
                        self.Qsjkfeed[j, 0] * self.TMudFeedj[j, L]
                
                    # Tracer deposition assumes storage occurs with 
                    # concentration equal to feed concentration
                    TracerStorage = (DC.Migrationfactor[j] * \
                        self.ExSed.OutMigration[0] * \
                        self.ExTracer[L].OutMigration[0] + \
                        self.ExSed.OutWidthChange[0] * \
                        self.ExTracer[L].OutWidthChange[0] + \
                        self.ExSed.OutVerticalChange[0] * \
                        self.ExTracer[L].OutVerticalChange[0]) / DC.p[j] + \
                        self.SinkSed[j, 0] * self.TMudFeedj[j, L]
                
                    self.TMudj[j, L] = (TracerInflux - TracerStorage) / \
                        self.Qsjk[j, 0]
                else: # set tracer concentration to zero if outgoing sediment 
                #flux is zero
                    self.TMudj[j, 0] = 0.
    
    def UpdateMudTracerAverages(self, Dj, DC):
        """
        for completeness, this computes Temporal, Volumetric, and \
        Deposition-weighted Averages for Tracer Concentrations in mud Load \
        and Feed for all sizes.
        
        Arguments:
            Dj -- [float]
            DC -- clsDurationCurve
        """
            
        self.TVolumeAvMudFeed = np.sum(self.Qsjkfeed[:, 0] * \
            np.transpose(self.TMudFeedj) * DC.p / self.QsAvkFeed[0], \
            axis = 1)
        self.TTemporalAvMudFeed = np.sum(np.transpose(self.TMudFeedj) * \
            DC.p, axis = 1)
        self.TTemporalAvMudLoad = np.sum(np.transpose(self.TMudj) * DC.p, \
            axis = 1)
    
