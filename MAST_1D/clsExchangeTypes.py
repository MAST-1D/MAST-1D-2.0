import numpy as np

class clsExchangeTypes:
	"""
	This defines the main types of even boundary-movement exchanges.
	Non-boundary movement exchanges are handled elsewhere.
	This can be used for storing either sediment volume fluxes or sediment
	tracer concentrations.
	
	Attributes:
		OutMigration -- [float]
		OutWidthChange -- [float]
		OutVerticalChange -- [float]
		InMigration -- [float]
		InWidthChange -- [float]
		InVerticalChange -- [float]
		NSizes -- int
		Initialized -- bool
	"""
	Initialized = False
	
	def __init__(self, NSizes):
		"""
		Arguments:
			NSizes -- int
		"""
		if not self.Initialized:
			self.NSizes = NSizes
			self.OutMigration = np.zeros(NSizes + 1)
			self.OutWidthChange = np.zeros(NSizes + 1)
			self.OutVerticalChange = np.zeros(NSizes + 1)
			self.InMigration = np.zeros(NSizes + 1)
			self.InWidthChange = np.zeros(NSizes + 1)
			self.InVerticalChange = np.zeros(NSizes + 1)
			self.Initialized = True
		else:
			raise RuntimeError('Tried to initiate clsExchangeType twice.')