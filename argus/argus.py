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
    ma_datenum    : matlab datetime numpy array
    
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