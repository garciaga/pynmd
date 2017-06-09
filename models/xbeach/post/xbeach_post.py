# -*- coding: utf-8 -*-
"""
Tools to manage xbeach output

External dependencies:
 

Internal dependencies:

"""

from __future__ import division,print_function

import numpy as np

# Compute runup
def nc_runup(nc,r_depth=0.01,yind=None):
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
    1. Test for 1D simulations.
    """
    
    # See if this is a 1d or 3d simulation
    simShape = nc.variables['zb'].shape
    # 2DH simulation
    if len(simShape) == 3:
        
        if yind is None:
            yind = np.round(simShape[1]/2)
            
        h   = nc.variables['zb'][:,yind,:].squeeze()
        eta = nc.variables['zs'][:,yind,:].squeeze()
        d   = eta - h
        x   = nc.variables['globalx'][yind,:].squeeze()    
    
    # 1D simulation   
    else:        
        
        # Load variables
        h   = nc.variables['zb'][:]
        eta = nc.variables['zs'][:]
        d   = eta - h
        x   = nc.variables['globalx'][:]    
       
    # Preallocate runup variable
    z_runup = np.zeros((d.shape[0],))
    x_runup = np.zeros_like(z_runup)
    
    # Loop over time, which is assumed to be on the first dimension.
    for aa in range(z_runup.shape[0]):
                   
        # Find the runup contour (search from left to right) 
        wdepth_ind = np.argmin(d[aa,:]>r_depth)-1
        
        # Store the water surface elevation in matrix
        z_runup[aa]= h[aa,wdepth_ind]
        
        # Store runup position
        x_runup[aa] = x[wdepth_ind]
    
    # Done
    return z_runup,x_runup
    

