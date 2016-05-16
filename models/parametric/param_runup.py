# -*- coding: utf-8 -*-
"""
Parametric runup tools

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

from __future__ import division,print_function

import numpy as np

#===============================================================================
# Stockdon 2006
#===============================================================================
def stockdon2006(H,L,B):
    '''
    Estimate runup based on the parametric relation by Stockdon et al 2006.
    
    USAGE:
    ------
    R2 = stockdon2006(H,L,B)
    
    PARAMETERS:
    -----------
    H  : Offshore wave height [m]
    L  : Offshore wave length [L]
    B  : Beach slope
    
    RETURNS:
    --------
    R2       : 2% Runup exceedence [m]
    setup    : Wave-induced setup [m]
    incSwash : Incident wave swash magnitude [m]
    igSwash  : Infragravity wave swash magnitude [m]
    
    NOTES:
    ------
    All values are multiplied by 1.1 to account for the skewness of the runup.
    
    REFERENCES:
    -----------
    Stockdon, H. F, R. A. Holman, P. A. Howd, A. H Sallenger, 2006: Empirical
        parameterization of setup, swash, and runup. Coastal Engineering, 53,
        573-588.
        
    '''
    
    # Make sure parameters are double
    H = np.double(H)
    L = np.double(L)
    B = np.double(B)
    
    # Compute incident swash (equation 11) 
    incSwash = 1.1 / 2 * 0.75 * B * (H*L)**0.5
    
    # Infragravity swash (Equation 12)
    igSwash = 1.1 / 2 * 0.06 * (H*L)**0.5
    
    # Compute R2% (Equation 19)
    setup = 1.1 * 0.35 * B * ((H * L)**0.5)
    swash = 1.1 / 2.0 * (H*L * (0.563 * B**2 + 0.004))**0.5    
    r2 = setup + swash
    
    return r2,setup,incSwash,igSwash
