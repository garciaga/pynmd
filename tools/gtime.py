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


def roundTime(dt=None, dateDelta=datetime.timedelta(minutes=1)):
    """Round a datetime object to a multiple of a timedelta
    dt : datetime.datetime object, default now.
    dateDelta : timedelta object, we round to a multiple of this, default 1 minute.
    Author: Thierry Husson 2012 - Use it as you want but don't blame me.
            Stijn Nevens 2014 - Changed to use only datetime objects as variables
    """
    roundTo = dateDelta.total_seconds()

    if dt == None : dt = datetime.datetime.now()
    seconds = (dt - dt.min).seconds
    # // is a floor division, not a comment on following line:
    rounding = (seconds+roundTo/2) // roundTo * roundTo
    return dt + datetime.timedelta(0,rounding-seconds,-dt.microsecond)
