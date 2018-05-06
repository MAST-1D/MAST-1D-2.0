# -*- coding: utf-8 -*-
"""
Created on Fri May 13 11:18:12 2016

@author: geography
"""
import os
import pickle
from clsOutputSpecs import clsOutputSpecs
            
class clsExport_object(object):
    
    """
    Unpickles the output object in the same directory as the object to avoid pickle problems
    """
    
    def __init__(self,folder):
        path = os.path.join(folder,'save.OutputSpecs')
        Output = pickle.load(open(path, "rb"))
        
        return Output