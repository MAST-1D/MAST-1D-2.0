import numpy as np

class clsTracerProperties(object):
    """
	Properties of radioisotopic tracer that may move with sediment
    
    Represents the properties of a given radioactive tracer that may
    be associated with sediment particles in a given size class.  
    Parameters are available for representing tracers produced in-situ
    in quartz by cosmogenic production or through atmospheric fallout. 
    Not yet fully implemented.
    
    Attributes
    ----------
        ProductionRate :  float
            Production rate of cosmogenic radionuclides in atoms/gram of SiO2/yr
        Lcj : array_like(float, length = 3)
            E-folding lengths (m) for the production of cosmogenic nuclides 
            in m for production process n (spallation, fast muons, slow muons)
        coj : array_like(float, length = 3)
            Contribution of the mechanism in the production of cosmogenic 
            nuclides for production proce(ss n (spallation, fast muons, slow muons)
        DecayConst : float 
            Radioactive decay constant for the cosmogenic nuclides (1/yr)
        FalloutRate : float 
            Fallout rate for tracer
        Name : str
            Name of tracer
            
    """
    Lcj = np.zeros(3)
    coj = np.zeros(3)