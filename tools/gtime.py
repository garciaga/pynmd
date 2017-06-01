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

# Import internal modules ------------------------------------------------------

# Datetime and epoch tools
import pynmd.argus as _gargus
datetime_to_epoch = _gargus.datetime_to_epoch
epoch_to_datetime = _gargus.epoch_to_datetime


# Matlab to datetime tools
import pynmd.matlab as _gmat
datenum_to_datetime = _gmat.datenum_to_datetime
datetime_to_datenum = _gmat.datetime_to_datenum

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
