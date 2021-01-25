import numpy as np

class clsExchangeTypes:
    """
    Defines the main types of boundary-movement exchanges.
    
    Stores size-specific fluxes of sediment volume or tracer associated with
    lateral exchange in a reservoir.  Non-boundary movement exchanges
    are handled elsewhere. Can be used for storing either sediment
    volume fluxes or sediment tracer concentrations.
    
    Parameters
    ----------
    NSizes : int
        Number of bed material sediment size classes. CONSIDER 
        RENAMING AS NBedSizes FOR CONSISTENCY WITH GSD.
    
    Attributes
    ----------
    OutMigration : array_like(float, length = NSizes + 1)
        Movement out of reservoir associated with lateral channel
        movement that does not lead to net width change. Size
        with index k = 0 represents washload.
    OutWidthChange : array_like(float, length = NSizes + 1)
        Movement out of reservoir associated with with net widening
        or narrowing. Size with index k = 0 represents washload.
    OutVerticalChange : array_like(float, length = NSizes + 1)
        Movement out of reservoir associated with vertical
        movement of one of the reservoir boundaries. Size
        with index k = 0 represents washload.
    InMigration : array_like(float, length = NSizes + 1)
        Movement into reservoir associated with lateral channel
        movement that does not lead to net width change. Size
        with index k = 0 represents washload.
    InWidthChange : array_like(float, length = NSizes + 1)
        Movement into reservoir associated with lateral channel
        movement that does not lead to net width change. Size
        with index k = 0 represents washload.    
    InVerticalChange : array_like(float, length = NSizes + 1)
        Movement into reservoir associated with vertical
        movement of one of the reservoir boundaries. Size
        with index k = 0 represents washload.
    NSizes : int
        Number of bed material sediment size classes. CONSIDER 
        RENAMING AS NBedSizes FOR CONSISTENCY WITH GSD.
    Initialized : bool.
        Flag indicating if object has been initialized.
        
    """
    
    Initialized = False
    
    def __init__(self, NSizes):
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