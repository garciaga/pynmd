# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 15:50:55 2015

@author: ggarcia

Tools for managing Arugs information and imagery
"""

import numpy as np
import datetime

#==============================================================================
# Tool to convert from epoch to datetime
#==============================================================================
def epoch_to_datetime(epoch):
    '''
    Tool to convert from UNIX Epoch time as used inn the Argus system into
    python's datetime.
    
    PARAMETERS
    ----------
    ma_datenum    : array with epoch times
    
    RETURNS
    -------
    py_datetime   : numpy array with ma_datenum in python's datetime format
    
    '''
    
    # Convert to double
    epoch = np.double(epoch)
    
    py_datetime = np.array([datetime.timedelta(seconds=epoch[aa]) +                             
                            datetime.datetime(1970,1,1,0,0,0)
                            for aa in range(epoch.shape[0])])
                   
    return py_datetime

#==============================================================================
# Tool to convert from datetime to epoch
#==============================================================================
def datetime_to_epoch(py_datetime):
    '''
    Tool to convert from UNIX Epoch time as used inn the Argus system into
    python's datetime.
    
    PARAMETERS
    ----------
    py_datetime   : numpy array with ma_datenum in python's datetime format   
    
    RETURNS
    -------
    epoch         : numpy array with epoch times
    
    '''
           
    # Convert to epoch time
    # Check if it is just one value
    if np.size(py_datetime) == 1:
        epoch = (py_datetime - 
                 datetime.datetime(1970,1,1,0,0,0)).total_seconds()
        epoch = np.array([epoch])
                         
    else:
        epoch = np.array([(aa - datetime.datetime(1970,1,1,0,0,0)).total_seconds()
                          for aa in py_datetime])

    return epoch
