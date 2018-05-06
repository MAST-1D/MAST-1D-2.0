import numpy as np

class clsDurationCurve:
	"""
	Attributes:
		Qw -- [float] Water Discharge in bin
		p -- [float] Fraction of time bin occurs
		Uc -- [float] Water velocity in channel in bin
		Uf -- [float] Water velocity on floodplain in bin
		Hc -- [float] Water depth in channel in bin
		Hf -- [float] Water depth on floodplain in bin
		Qwc -- [float] Discharge of water in the channel in bin
		Qwf -- [float] Discharge of water on the floodplain in bin
		Qs -- [float] Sediment discharge in bin -- probably should remove
			this since all load is handled in clsLoad
		WSE -- [float] Water surface elevation
		Sf -- [float] Friction Slope
		Dfloodplain -- [float] Floodplain deposition in bin--should check 
			whether this is stored elsewhere
		MigrationFactor -- [float] Factor that determines how much of the
			annual average migration occurs for each of duration curve
		Initialized -- bool
	"""
	Initialized = False
	
	def __init__(self, NFlows):
		"""
		Arguments:
			NFlows -- float
		"""
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
		g = 9.81
		return self.Uc[j] / (g * self.Hc[j]) ** 0.5
	
	def WeightByFDC(self, x):
		WeigthByFDC = 0.
		for i in range(len(p)):
			WeightByFDC = x[j] * self.p[j] + WeightByFDC
		return WeightByFDC