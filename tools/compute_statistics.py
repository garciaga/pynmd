#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Compute statistics

"""
__author__ = "Saeed Moghimi"
__copyright__ = "Copyright 2017, UCAR/NOAA"
__license__ = "GPL"
__version__ = "1.0"
__email__ = "moghimis@gmail.com"


import numpy as np
import netCDF4 as n4


def find_nearest2d(xvec,yvec,xp,yp):
    """
    In: xvec, yvec of the grid and xp,yp of the point of interst
    Retun: i,j,proximity of the nearset grid point
    """

    dist = np.sqrt((xvec-xp)**2+(yvec-yp)**2)
    i,j  = np.where(dist==dist.min())
    return i[0],j[0],dist.min()
    
def find_nearest1d(xvec,yvec,xp,yp):
    """
    In: xvec, yvec of the grid and xp,yp of the point of interst
    Retun: i,j,proximity of the nearset grid point
    """

    dist = np.sqrt((xvec-xp)**2+(yvec-yp)**2)
    i = np.where(dist==dist.min())
    return i[0],dist.min()

def find_nearest_time(dates,date):
    """
    Return nearest time indexs
    Inp: dates: timedate vector
         date: timedate pint
    Out: index
    
    """
    dtsec = []
    for tmp in dates:
        dtsec.append( np.abs( (tmp-date).total_seconds()  ) )  
    return np.argmin(dtsec)    

    
        
def statatistics(data,model):
    """ 
        Calculate statistics 
    """

    delta     = model- data
    r2        = 1.-((data-model)**2).sum()/((data-data.mean())**2).sum()
    bias      = model.mean()-data.mean()
    rbias     = bias / data.mean()
    rmse      = np.sqrt((delta**2).mean())
    nrmse1    = rmse / (data.max() - data.min())
    nrmse2    = rmse/ abs(data.mean() )
    nrmse3    = np.sqrt( (delta**2).sum() / (data**2).sum() )
    big_error = (np.abs(delta)).max()
    mae       = np.abs(delta).mean()
    cor      = (((data-data.mean())*(model-model.mean())).mean()/data.std()/model.std())
    #cor       = np.corrcoef([model,data])
    skill     = 1 -np.sum((data-model)**2) / np.sum((np.abs(model-np.mean(data)) + np.abs(data-np.mean(data)))**2);
    peak      = model.max() - data.max()
    # index of agreement
    ia = 1 -  (np.sum((data-model)**2))/(np.sum((np.abs(model-np.mean(data))+np.abs(data-np.mean(data)))**2))

    return dict(bias=bias,rmse=rmse,r2=r2,skill=skill,peak=peak,mae=mae,cor=cor,rbias=rbias,ia=ia)

    
def model_on_data(data_dates, model_dates, model_val):
    print '  >>>>>>>>>>>>>>>   '
    units     = 'seconds since 2012-04-01 00:00:00'
    data_sec  = n4.date2num(data_dates , units)
    model_sec = n4.date2num(model_dates, units)
    return np.interp(data_sec, model_sec, model_val)