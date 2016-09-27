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
    