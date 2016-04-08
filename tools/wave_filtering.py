# -*- coding: utf-8 -*-
"""
Functions for wave filtering

Dependencies:
-------------
numpy

Internal Dependencies:
----------------------


Created on Mon Apr  4 21:49:54 2016

@author: ggarcia
"""

from __future__ import division,print_function

# Import modules
import numpy as np

# ==============================================================================
# Infragravity wave filtering based on Guza et al 1984
# ==============================================================================
def long_wave_separate_guza84(eta,u,h):
    '''
    Separate long waves into seaward and shoreward propagating components 
    based on collocated measurements of water surface elevation and velocity.
    
    USAGE:
    ------
    etaInc,etaRef = long_wave_separate_guza84(eta,u,h)
    
    PARAMETERS:
    -----------
    eta    : Water surface elevation [m]
    u      : Depth-averaged velocity in the direction of wave propagation [m]
    h      : Water depth at each point [m]
    
    RETURNS:
    --------
    etaInc : Incident wave component [m]
    etaRef : Reflected (seaward propagating) wave component [m]

    FORMULAE:
    ---------
    etaInc = 0.5 * (eta + u * (h/g)**0.5)
    etaRef = 0.5 * (eta - u * (h/g)**0.5)
    
    NOTES:
    ------
    - Only works for longwaves in the signal. You may need to low pass filter
      your signal.
    - Has been tested for one dimensional wave propagation
    - Will operate on the last dimension of eta and u
    
    
    REFERENCES:
    -----------
    Guza, R. T., E. B. Thornton, R. A. Holman, 1984: Swash on steep and shallow
        beaches. ICCE 19th, 708 - 723.
     
    '''
    
    # Convert variables to arrays just in case
    h = np.asarray(h)
    eta = np.asarray(eta)
    u = np.asarray(u)
    
    if eta.ndim == 1:
        
        etaInc = (eta + ((h/9.81)**0.5) * u)/2.0
        etaRef = (eta - ((h/9.81)**0.5) * u)/2.0
        
    elif eta.ndim == 2:
        
        # Preallocate variables
        etaInc = np.zeros_like(eta) * np.NaN
        etaRef = np.zeros_like(eta) * np.NaN
    
        for aa in range(eta.shape[1]):
            etaInc[:,aa] = (eta[:,aa] + ((h[aa]/9.81)**0.5) * u[:,aa])/2.0
            etaRef[:,aa] = (eta[:,aa] - ((h[aa]/9.81)**0.5) * u[:,aa])/2.0
    
    else:        
        print('Error: eta.ndim > 2')
        etaInc = np.array([])
        etaRef = np.array([])
        
    return etaInc,etaRef


    