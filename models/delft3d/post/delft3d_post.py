__author__ = "Saeed Moghimi"
__copyright__ = "Copyright 2015, Portland State University"
__license__ = "GPL"
__version__ = "1.0"
__email__ = "moghimis@gmail.com"

"""
#slice_plot moved to pynmd.plotting.plot_routines 

"""



import netCDF4
import netcdftime

import matplotlib.pyplot as plt
import numpy as np
import os,sys

from pynmd.plotting.vars_param import *
import pynmd.plotting.plot_settings as ps
import datetime

#import seaborn as sns
#plt.style.use('seaborn-white')





### Funcs
def read_d3d_trim(fname, all = True):
    """
    To read D3D trim*.nc file
    Reads: x,y,alfa(local rotation) of cells
    also returns nc.varibles for reading other vars in case
    
    return xcor,ycor,z_cen,alfa,maskxy,ncv
    """
    
    nc  = netCDF4.Dataset(fname)
    ncv = nc.variables
    
    if all:
        xcor  = ncv['XCOR'][:]
        ycor  = ncv['YCOR'][:]
        dep   = ncv['DP0'][:]
        dep = np.ma.masked_where(dep<-100,dep)
        maskxy = dep.mask
        zeta   = ncv['S1'][:]
        nt     = zeta.shape[0]
        alfa  = np.deg2rad (ncv['ALFAS'][:]) 
        #
        sig_layer_face   = ncv['SIG_INTF'][:]
        n_sig_face = len(sig_layer_face)
        #z_face = np.zeros((n_sig_face,xcor.shape[0],xcor.shape[1]))
        #for iz in range(n_sig_face):
        #    z_face[iz] =  sig_layer_face[iz] * dep     
        #
        #sigl = sig_layer_face[:-1] -sig_layer_face[1:]
        
        sig_layer_center = ncv['SIG_LYR'][:]
    #     sigh  =   sig_layer_face  [:-1]   - sig_layer_face[1:]
    #     
    #     h     = np.zeros((sigh.shape[0],xcor.shape[0],xcor.shape[1]))
    #     z_cen = np.zeros((nt,sigh.shape[0],xcor.shape[0],xcor.shape[1]))
    #     for it in range(nt):
    #         for iz in range(sigh.shape[0]):
    #             h[iz] = sigh[iz] * (dep +   zeta[it])
    #         
    #         for iz in range(sigh.shape[0]):
    #             z_cen[it,iz,:] =  np.cumsum(h,0)[iz] - (dep + zeta[it])      
      
        
        
        n_sig_cen = len(sig_layer_center)
        #
        z_cen = np.zeros((nt,n_sig_cen,xcor.shape[0],xcor.shape[1]))
        for it in range(nt):
            for iz in range(n_sig_cen):
                z_cen[it,iz] =  (1+sig_layer_center[iz])  * (dep  +   zeta[it]) - dep  
     
        z_cen = np.ma.masked_where(z_cen < -10, z_cen)
        #z_cen = 10 - z_cen
        
        return xcor,ycor,z_cen,alfa,maskxy,ncv
    else:
        return ncv

def read_d3d_time(netcdf_vars):
    """
    Read time from netcdf file and return Datetime vector
    """
    utime = netcdftime.utime(netcdf_vars['time'].units)
    dates = utime.num2date(netcdf_vars['time'][:])
    return dates

def rotate_3d_vel(u1,v1,w1=None,ang=0):
    import okean.calc as calc

    """
    horizontal rotation of D3D 3D flow field. Fro Curvilinear grids 
    """
    if w1 is not None:
        u1,v1,w1   = calc.rot3d (u1,v1,w1,ang,0.,inverse=False)    
    else:
        u1,v1,w1   = calc.rot3d (u1,v1,v1*0,ang,0.,inverse=False)    
    return u1,v1,w1




    
##### end of Funcs
