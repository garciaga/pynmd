# -*- coding: utf-8 -*-
"""
Tools to manage xbeach output

External dependencies:
 

Internal dependencies:

"""

from __future__ import division,print_function

import numpy as np

# Compute runup
def nc_runup(nc,r_depth=0.01):
    """
    
    Function to compute runup from netcdf file.
    
    Parameters:
    -----------
    nc           : NetCDF file handle
    r_depth      : Runup depth [m] (defaults to 0.01m)
    
    Output:
    -------
    runup        : Water surface elevation time series relative to SWL given
                   a contour depth [m]
    x_runup      : Across-shore location of runup time series [m]
                       
    Notes:
    ------
    
    """
    
    # Load variables
    h   = nc.variables['zb'][:]
    eta = nc.variables['zs'][:]
    d   = eta - h
    x   = nc.variables['globalx'][:]    
    
    # Need to manage 2d files
    
    # Preallocate runup variable
    z_runup = np.zeros((d.shape[0],d.shape[1]))
    x_runup = np.zeros_like(z_runup)
    
    # Loop over time, which is assumed to be on the first dimension.
    for aa in range(z_runup.shape[0]):
        
        # Loop over space
        for bb in range(z_runup.shape[1]):
            
            # Find the runup contour (search from left to right) 
            wdepth_ind = np.argmin(d[aa,bb,:]>r_depth)-1
            
            # Store the water surface elevation in matrix
            z_runup[aa,bb]= eta[aa,bb,wdepth_ind]
            
            # Store runup position
            x_runup[aa,bb] = x[bb,wdepth_ind]
    
    # Done
    return z_runup,x_runup
    

