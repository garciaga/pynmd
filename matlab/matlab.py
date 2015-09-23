"""
Tools for managing matfiles or matlab suff

Gabriel Garcia Medina
ggarcia@coas.oregonstate.edu
May 2015

Dependencies:
-------------
datetime, numpy

Internal dependencies:
----------------------
none

"""

# Import modules
import datetime
import numpy as np

# ==============================================================================
# Tools to convert from datenum (matlab) to datetime (python)
# ==============================================================================
def datenum_to_datetime(ma_datenum):
    '''
    Tool to convert from Matlab's datenum to python's datetime
    
    PARAMETERS
    ----------
    ma_datenum    : matlab datetime numpy array
    
    RETURNS
    -------
    py_datetime   : numpy array with ma_datenum in python's datetime format
    
    '''
    
   
    py_datetime = np.array([datetime.datetime.fromordinal(int(ma_datenum[aa])) +
                   datetime.timedelta(days=ma_datenum[aa]%1) -
                   datetime.timedelta(days=366)
                   for aa in range(ma_datenum.shape[0])])
                   
    return py_datetime
    
