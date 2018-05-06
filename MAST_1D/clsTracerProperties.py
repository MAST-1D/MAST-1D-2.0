import numpy as np

class clsTracerProperties(object):
	"""
	Attributes
		ProductionRate -- float (production rate of cosmogenic nuclides in 
						  atoms/gram of SiO2/yr)
		Lcj[3] -- [float] (e-folding lengths (m) for the production of 
						  cosmogenic nuclides in m for production process n)
		coj[3] -- [float] (contribution of the mechanism in the production 
						  of cosmogenic nuclides for production process n 
						  (spallation, fast muons, slow muons))
		DecayConst -- float (radioactive decay constant for the cosmogenic
					  nuclides in 1/yr)
		FalloutRate -- float (Fallout rate of cosmogenic nuclides)
		Name -- str
	"""
	Lcj = np.zeros(3)
	coj = np.zeros(3)