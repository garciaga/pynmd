# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 14:10:34 2016

@author: ggarcia

Tools for time management

Dependencies:
-------------
datetime

Internal Dependencies:
----------------------


"""

from __future__ import division,print_function
import datetime as _datetime
import numpy as _np

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
    epoch = _np.double(epoch)
    
    if epoch.size == 1:
        py_datetime = (_datetime.timedelta(seconds=epoch) + 
                       _datetime.datetime(1970,1,1,0,0,0))
    else:
        py_datetime = _np.array([_datetime.timedelta(seconds=epoch[aa]) +                             
                                 _datetime.datetime(1970,1,1,0,0,0)
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
    if _np.size(py_datetime) == 1:
        epoch = (py_datetime - 
                 _datetime.datetime(1970,1,1,0,0,0)).total_seconds()
        epoch = _np.array([epoch])
                         
    else:
        epoch = _np.array([(aa - _datetime.datetime(1970,1,1,0,0,0)).total_seconds()
                           for aa in py_datetime])

    return epoch


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

    # Convert to double
    ma_datenum = _np.double(ma_datenum)
       
    if ma_datenum.size == 1:
        py_datetime = (_datetime.datetime.fromordinal(int(ma_datenum)) + 
                       _datetime.timedelta(days=ma_datenum%1) - 
                       _datetime.timedelta(days=366))
    else:   
        py_datetime = _np.array([_datetime.datetime.fromordinal(int(ma_datenum[aa])) +
                                 _datetime.timedelta(days=ma_datenum[aa]%1) -
                                 _datetime.timedelta(days=366)
                                 for aa in range(ma_datenum.shape[0])])

    return py_datetime
    

# ==============================================================================
# Convert from datetime (python) to datenum (matlab)
# ==============================================================================
def datetime_to_datenum(py_datetime):
    '''
    Tool to convert from Matlab's datenum to python's datetime
    
    PARAMETERS
    ----------
    py_datetime   : numpy array with ma_datenum in python's datetime format
    
    RETURNS
    -------
    ma_datenum    : matlab datetime numpy array
    
    NOTES:
    -------
    Matlab datenum is a fractional number that represents dates as days from
    0 January 0000.
    
    '''
   
    # Define constants
    datetimeBase = _datetime.datetime(1,1,1)     # Lowest possible datetime
    matlabOffset = _datetime.timedelta(days=367) # Matlab starts at year 0
    
    if _np.size(py_datetime) == 1:
        ma_datetime = (py_datetime - datetimeBase + 
                       matlabOffset).total_seconds()/(3600.0 * 24.0)
    else:
        ma_datetime = _np.array([(py_datetime[aa] - datetimeBase + 
                                  matlabOffset).total_seconds()/(3600.0 * 24.0)
                                for aa in range(py_datetime.shape[0])])
                       
    return ma_datetime
    

# Wrappers to link pynmd.argus and pynmd.matlab time tools ---------------------
def datenum_to_epoch(ma_datenum):
    """
    Tool to convert from Matlab's datenum to UNIX Epoch time
    
    PARAMETERS
    ----------
    ma_datenum    : matlab datetime numpy array
    
    RETURNS
    -------
    epoch         : numpy array with epoch times
    
    """
    
    # Convert to python time 
    py_datetime = datenum_to_datetime(ma_datenum)
    
    # Convert to epoch 
    epoch = datetime_to_epoch(py_datetime)
    
    return epoch


# Wrappers to link pynmd.argus and pynmd.matlab time tools ---------------------
def epoch_to_datenum(epoch):
    """
    Tool to convert from UNIX Epoch time to Matlab's datenum
    
    PARAMETERS
    ----------
    epoch         : numpy array with epoch times
    
    RETURNS
    -------
    ma_datenum    : matlab datetime numpy array
    
    """

    # Convert to python time 
    py_datetime = epoch_to_datetime(epoch)
    
    # Convert to epoch 
    ma_datenum = datetime_to_datenum(py_datetime)
    
    return ma_datenum
    
# Day of year ------------------------------------------------------------------
def dayOfYear(year,month,day):
    """
    Returns the day of the year
    
    PARAMETERS
    ----------
    year      : year
    month     : month
    day       : day
    
    RETURNS:
    --------
    dayOfYear : day of the year
    
    """    
    
    a = _datetime.datetime(year,month,day) - _datetime.datetime(year,1,1)
    
    return a.days + 1


def roundTime(dt=None, dateDelta=_datetime.timedelta(minutes=1)):
    """Round a datetime object to a multiple of a timedelta
    dt : datetime.datetime object, default now.
    dateDelta : timedelta object, we round to a multiple of this, default 1 minute.
    Author: Thierry Husson 2012 - Use it as you want but don't blame me.
            Stijn Nevens 2014 - Changed to use only datetime objects as variables
    """
    roundTo = dateDelta.total_seconds()

    if dt == None : dt = _datetime.datetime.now()
    seconds = (dt - dt.min).seconds
    # // is a floor division, not a comment on following line:
    rounding = (seconds+roundTo/2) // roundTo * roundTo
    return dt + _datetime.timedelta(0,rounding-seconds,-dt.microsecond)
#===============================================================================
# Monthly vector
#===============================================================================
def monthVector(yr1,yr2):
    """
    Create a monthly vector in datetime format between both given years
    
    PARAMETERS:
    -----------
    yr1 : Start Year
    yr2 : End Year (inclusive)
    
    RETURNS:
    --------
    dateVec : Vector of monthly entries including Jan yr1 and Dec yr2
    
    NOTES:
    ------
    yr1 and yr2 are integers
    """
    
    # Compute a safety range in case this function goes crazy
    time1 = _datetime.datetime(yr1,1,1)
    time2 = _datetime.datetime(yr2,12,31)   
    dateRange = time2 - time1
    
    # Preallocate variables
    timeVec = []
    cnt = 1 # Global counter and safety variable
    monthCnt = 0
    tmpYear = yr1
    
    while cnt < dateRange.days:
        
        cnt += 1       # global counter
        monthCnt += 1  # Next month
        
        # Month 13 is january of next year
        if monthCnt%13 == 0:
            monthCnt = 1
            tmpYear += 1
    
        # Current month as datetime object
        tmpTime = _datetime.datetime(tmpYear,monthCnt,1)
                
        # Check if the final month has been reached
        if tmpYear > yr2:
            break
        
        # If everything looks good allocate
        timeVec.append(tmpTime)
        
    # Convert to array
    timeVec = _np.array(timeVec)
    
    return timeVec
