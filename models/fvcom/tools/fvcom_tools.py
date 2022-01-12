"""
Collection of functions to postprocess fvcom data
"""

import numpy as _np
import datetime as _datetime

# ==============================================================================
# Get the different fvcom times
# ==============================================================================
def datetimeToMJD(timeVec):
    """
    Return Modified Julian Dates based on a vector of datetime objects

    PARAMETERS:
    -----------
    timeVec: Vector of datetime objects

    OUTPUT:
    -------
    Dictionary containing
    time : float
        Fractional days since epoch
    Itime : int
        Integer number of days since epoch
    Itime2 : int
        Milliseconds after Itime
    Times : string with length 26
        YYYY-mm-ddTHH:MM:SS.ffffff        

    """
    
    # Assert array
    timeVec = _np.asarray(timeVec)

    # Preallocate output dictionary
    out = {'time':_np.zeros((timeVec.shape[0]),dtype=_np.float),
           'Itime':_np.zeros((timeVec.shape[0]),dtype=_np.int),
           'Itime2':_np.zeros((timeVec.shape[0]),dtype=_np.int),
           'Times':_np.zeros((timeVec.shape[0],26),dtype=_np.str)
           }

    # Set the MJD epoch
    epoch = _datetime.datetime(1858,11,17,0,0,0)

    # Transformations
    for ii in range(timeVec.shape[0]):

        # MJD
        dt = timeVec[ii] - epoch
        out['time'][ii] = dt.days + (dt.seconds / 86400.0)

        # Integer time
        out['Itime'][ii] = _np.int(_np.floor(out['time'][ii]))

        # Milliseconds after integer time
        out['Itime2'][ii] = ((out['time'][ii] - out['Itime'][ii]) * 
                              24 * 3600 * 1000)
        
        # String
        strTime = timeVec[ii].strftime('%Y-%m-%dT%H:%M:%S.%f')
        for jj in range(len(strTime)):
            out['Times'][ii,jj] = strTime[jj]

    return out
