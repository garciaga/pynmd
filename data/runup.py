# -*- coding: utf-8 -*-
"""
A series of tools to manage runup

Author:
-------
Gabriel Garcia Medina
    ggarcia@coas.oregonstate.edu
    
Dependencies:
-------------
numpy

Internal dependencies:
----------------------
   
"""

# Import modules
import numpy as np

# ==============================================================================
# Runup maxima
# ==============================================================================
def runup_maxima(x,ot,sten):
    """
    Function to identify local minima from a water surface elevation time series
    
    USAGE:
    ------
    max_ind = runup_maxima(x,ot,sten)
    
    PARAMETERS:
    -----------
    x       : Runup time series [m]
    ot      : Time vector, must have the same length as x and time in seconds. 
    sten    : Minimum runup period to consider [s]
    
    RETURNS:
    --------
    max_ind : Vector of indices corresponding to the local maxima of the input
              time series.    
    
    NOTES:
    ------
    All inputs should be numpy arrays
    
    """

    # Local minima analysis to identify the waves -----------------------------

    # Compute the first derivative of the data
    dx           = np.zeros_like(x)
    dx[1:]       = x[1:] - x[0:-1]

    # Take the second derivative
    dx2          = np.zeros_like(dx) * np.NAN
    dx2[1:-1]    = (x[2:] - 2*x[1:-1] + x[0:-2])

    # Find local extrema by finding where the derivative of the data is zero
    # Numerically it is best to find where the derivative changes sign and 
    # record the position where this happens. 
    dx_sign      = np.zeros_like(dx)
    dx_sign[1:]  = dx[1:] * dx[0:-1]
    ind_ext      = np.where(dx_sign<0)[0] - 1

    # Identify the local minima
    # ind_min_tmp  = dx2[ind_ext]>0
    # ind_min      = ind_ext[ind_min_tmp]

    # Identify the local maxima
    ind_max_tmp  = dx2[ind_ext]<0
    ind_max      = ind_ext[ind_max_tmp]

    # Small wave filter
    wave_per     = np.ones_like(ind_max) * (sten + 1)
    wave_per[1:] = ot[ind_max][1:] - ot[ind_max][0:-1]
    ind_noise    = np.where(wave_per < sten)[0]

    # Filter small waves by looking at the identified minima and the previous 
    # one. The point closer to a local maxima should be removed.
    if sum(ind_noise)>0:
        
        # Flag for filtering more than once
        small_runup = True
        
        # Counter variable for maximum number of iterations
        kk = 0
        
        while small_runup:
            
            # Counter variable to avoid infinite loops
            kk += 1
            if kk > 4:
                print('Maximum number of iterations reached. Quitting...')
                break
            
            # Local maxima to remove
            rmv_ind = []
            
            # Loop over short waves
            for aa in range(ind_noise.shape[0]):
                
                # Time of the consecutive local minima
                ot_maxs = ot[ind_max][ind_noise[aa]-1:ind_noise[aa]+1]
                
                # Keep the maximum runup only 
                tmp_min_ind = np.argmin(x[ind_max][ind_noise[aa]-1:
                                                   ind_noise[aa]+1])
                
                if tmp_min_ind == 0:
                    rmv_ind.append(ind_noise[aa]-1)
                else:
                    rmv_ind.append(ind_noise[aa]) 
                
            # Remove short runup events 
            ind_max = np.delete(ind_max,rmv_ind)
            
            # Small wave filter
            wave_per     = np.ones_like(ind_max) * (sten + 1)
            wave_per[1:] = ot[ind_max][1:] - ot[ind_max][0:-1]
            ind_noise    = np.where(wave_per < sten)[0]
            
            if sum(ind_noise) > 0:
                small_runup = True
            else:
                small_runup = False
                        

    # Return minimum indices
    return ind_max


