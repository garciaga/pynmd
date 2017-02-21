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


def stockdon2006Dissip(H,L):
    '''
    Estimate runup based on the parametric relation by Stockdon et al 2006 for
    dissipative beaches.
    
    The deep water surf similarity parameter should be less than 0.3 to use this
    formula.
    
    USAGE:
    ------
    R2 = stockdon2006Dissip(H,L)
    
    PARAMETERS:
    -----------
    H  : Offshore wave height [m]
    L  : Offshore wave length [L]
    
    RETURNS:
    --------
    R2 : 2% Runup exceedence [m]
    
    NOTES:
    ------
    R2 = 0.043*(HL)**0.5 for ssp < 0.3
    
    REFERENCES:
    -----------
    Stockdon, H. F, R. A. Holman, P. A. Howd, A. H Sallenger, 2006: Empirical
        parameterization of setup, swash, and runup. Coastal Engineering, 53,
        573-588.
        
    '''
    
    # Make sure parameters are double precision
    H = np.double(H)
    L = np.double(L)
    
    # Compute R2% (Equation 18)   
    r2 = 0.043 * ((H * L)**0.5)
    
    return r2

def ruggiero2001(H,L,B):
    '''
    Estimate runup based on the parametric relation by Ruggiero et al 2001.

    
    USAGE:
    ------
    R2 = ruggiero2001(H,L,B)
    
    PARAMETERS:
    -----------
    H  : Offshore wave height [m]
    L  : Offshore wave length [L]
    B  : Beach slope
    
    RETURNS:
    --------
    R2 : 2% Runup exceedence [m]
    
    NOTES:
    ------
    R2 = 0.27*(BHL)**0.5
    
    REFERENCES:
    -----------
    Ruggiero, P, P. D Komar, W. G. McDougal, J. J. Marra, and R. A. Beach, 2001:
        Wave Runup, Extreme Water Levels and the Erosion of Properties Backing
        Beaches. Journal of Coastal Research, 17 (2), 407-419.
        
    '''
    
    # Make sure parameters are double precision
    H = np.double(H)
    L = np.double(L)
    B = np.double(B)
    
    # Compute R2% (Equation 5)   
    r2 = 0.27 * ((B *H * L)**0.5)
    
    return r2

#===============================================================================
# Guza and Feddersen
#===============================================================================
def guza2012(H,L,Fw,Ds=1.0):
    '''
    Estimate the expected infragravity runup based on the parametric relation 
    by Guza and Feddersen 2012.

    
    USAGE:
    ------
    R2 = guza2012(H,L,B,Fw,Ds)
    
    PARAMETERS:
    -----------
    H  : Offshore wave height [m]
    L  : Offshore wave length [L]
    Fw : Frequency width [Hz]
    Ds : Deep water directional spread [radians] 
    
    RETURNS:
    --------
    R2 : 2% Runup exceedence [m]
    
    NOTES:
    ------
    R2 = (-0.013 ln(Fp*Ds/Fw) + 0.058)*(HL)**0.5
    
    This equation was derived for iribarren number < 0.4.
    
    
    REFERENCES:
    -----------
    Guza, R. T. and F. Feddersen, 2012: Effect of wave frequency and directional
        spread on shoreline runup. Geophysical Research Letters, 39, L11607, 
        doi:10.1029/2012GL051959.
        
    '''
    
    # Make sure parameters are double precision
    H  = np.double(H)
    L  = np.double(L)
    B  = np.double(B)
    Fw = np.double(Fw)
    Ds = np.double(Ds)
    
    Fp = (2.0 * np.pi * L / 9.81)**-0.5
    
    # Compute R2% IG (Equation 4)   
    r2 =  (-0.013 * np.log(Fp/Fw*Ds) + 0.058) * ((H * L)**0.5)
    
    return r2
