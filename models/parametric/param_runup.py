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
physics.waves
   
"""

from __future__ import division,print_function

import numpy as np
import pynmd.physics.waves as _gwaves

#===============================================================================
# Stockdon 2006
#===============================================================================
def stockdon2006(H,L,B):
    """
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
        
    """
    
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

#===============================================================================
# Peter Ruggiero's Equation
#===============================================================================
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
    """
    Estimate the expected infragravity swash based on the parametric relation 
    by Guza and Feddersen 2012.

    
    USAGE:
    ------
    Sig = guza2012(H,L,Fw,Ds)
    
    PARAMETERS:
    -----------
    H  : Offshore wave height [m]
    L  : Offshore wave length [L]
    Fw : Frequency width [Hz]
    Ds : Deep water directional spread [radians] 
    
    RETURNS:
    --------
    Sig : Infragravity swash [m]
    
    NOTES:
    ------
    R2% = A(setup + S/2) = A(setup + 0.5*(Sig**2 + Sinc**2)**0.5)

    Sig = (-0.013 ln(Fp*Ds/Fw) + 0.058)*(HL)**0.5 (Eq 4)
    
    This equation was derived for iribarren number < 0.4.
    
    REFERENCES:
    -----------
    Guza, R. T. and F. Feddersen, 2012: Effect of wave frequency and directional
        spread on shoreline runup. Geophysical Research Letters, 39, L11607, 
        doi:10.1029/2012GL051959.
        
    """
    
    # Make sure parameters are double precision
    H  = np.double(H)
    L  = np.double(L)
    Fw = np.double(Fw)
    Ds = np.double(Ds)
    
    Fp = (2.0 * np.pi * L / 9.81)**-0.5
    
    # Compute R2% IG (Equation 4)   
    r2 =  (-0.013 * np.log(Fp/Fw*Ds) + 0.058) * ((H * L)**0.5)
    
    return r2

#===============================================================================
# Hajime Mase's Equation
#===============================================================================
def mase1989(H,L,B):
    '''
    Estimate runup based on the parametric relation by Mase 1989.

    
    USAGE:
    ------
    R = mase1989(H,L,B)
    
    PARAMETERS:
    -----------
    H  : Offshore wave height [m]
    L  : Offshore wave length [L]
    B  : Beach slope
    
    RETURNS:
    --------
    Dictionary containing
    Rmax    : Maximum runup [m]
    R2      : 2% Runup exceedence [m]
    R10     : 10% Runup exceedence [m]
    R33     : Significant Runup [m] 
    Rmean   : Mean runup [m]
    boreRed : Bore reduction factor
    Tmean   : Mean runup period [s]
    
    NOTES:
    ------
    Rmax  = H*2.32*(ssp)**0.77
    R2    = H*1.86*(ssp)**0.71
    R10   = H*1.70*(ssp)**0.71
    R33   = H*1.38*(ssp)**0.70
    Rmean = H*0.69*(ssp)**0.69
    ssp   = Surf similarity parameter
    Tmean = [(2*pi*L/g)**0.5]/boreRed
    
    REFERENCES:
    -----------
    Mase, H., 1989: Random Wave Runup Height on Gentle Slope. Journal of 
        Waterway, Port, Coastal, and Ocean Engineering, 115 (5), 649-661.
        
    '''
    
    # Make sure parameters are double precision
    H = np.double(H)
    L = np.double(L)
    B = np.double(B)
    
    # Surf similarity parameter
    ssp = _gwaves.iribarren(B, H, L, verbose=False)
    
    # Compute Runup (Equation 6)   
    rmax  = H * 2.32 * (ssp**0.77)
    r2    = H * 1.86 * (ssp**0.71)
    r10   = H * 1.70 * (ssp**0.71)
    r33   = H * 1.38 * (ssp**0.70)
    rmean = H * 0.69 * (ssp**0.69)
    
    # Compute bore reduction (Equation 7)
    if ssp > 3.57:
        boreRed = 1.0
    elif ssp > 0.91 and ssp <= 3.57:
        boreRed = 0.70 * (ssp**0.28)
    else:
        boreRed = 0.72 * (ssp**0.58)
    
    # Mean runup crest period (Equation 8)
    Tmean = ((2.0*np.pi*L/9.81)**0.5)/boreRed
    
    return {'Rmax':rmax,'R2':r2,'R10':r10,'R33':r33,'Rmean':rmean,
            'boreRed':boreRed,'Tmean':Tmean}

def blenkinsopp2022(H,dtoe_swl,Bberm,Bsand=None,gamma=0.87,gammas_sw=0.59):
    """
    Estimate runup for composite beaches based on the parametric relation 
    by Blenkinsopp et al 2022.

    
    USAGE:
    ------
    R = blenkinsopp2022(H,dtoe,Bberm, Bsand)
    
    PARAMETERS:
    -----------
    H        : Offshore wave height [m]
    dtoe_swl : still water depth at cobble berm toe [m]
    Bsand    : Cobble berm slope
    Bbeach   : (optional) Beach slope
    gamma    : (optional) Saturated wave breaking paramter
               default=0.87 per Section 4.2 of paper
    
    RETURNS:
    --------
    R2_17 : Equation 17 2% Exceedance runup [m]
    R2_21 : Equation 21 2% Exceedance runup [m]
    R2_23 : Equation 23 2% Exceedance runup [m]
    eta   : wave induced setup [m]

    
    NOTES:
    ------
    Conlin et al. (2025) found this formula performs well only when the berm
    is on inundation regime (swash zone entirely on the berm).
    
    REFERENCES:
    -----------
    Blenkinsopp, C. E., Bayle, P. M., Mariins, K., Foss, O. W., Almeida, L. P.,
        Kaminsky, G. M., Schimmels, S., Matsumoto, H., 2022: Wave runup on 
        composite beaches and dynamic cobble berm revetments. Coastal
        Engineering, 176, https://doi.org/10.1016/j.coastaleng.2022.104148
        
    """

    # Equation 17
    R2_17 = 4.59 * Bberm * dtoe_swl + 0.75

    # Setup (Equation 17)
    etabar = 0.17 * H

    # Equations 21 and 23
    if Bsand is not None:

        # Equation 13
        lsz = (5/3*H - dtoe_swl) / np.tan(Bsand) + dtoe_swl / np.tan(Bberm)

        # Eta toe (Equation 22)
        eta_toe = 3.33 * 10**-4 * lsz + 0.12

        # Water depth at toe of structure
        dtoe = dtoe_swl + eta_toe
    
        # Wave height at the toe of the structure
        # Hm0 following the paper nomenclature
        if dtoe_swl > 0:
            Hm0 = np.min([H,dtoe*gamma])
        else:
            Hm0 = 0.0

        # Equation 21
        R2_21 = 0.19 * H + 3.11 * Hm0 * np.tan(Bberm) + 0.26
      
        # Equation 9 (typo in preprint refers to Equation 8)
        Ssw = 0.23 + 6.79 * Hm0 * np.tan(Bberm)

        # Equation 15 infragravity swash
        Sig = 0.0030*lsz + 0.20

        # Runup (Equation 23)
        R2_23 = 1.1 * (etabar + 0.5 * (Ssw**2 + Sig**2)**0.5)
        
    else:
        R2_21 = np.NAN
        R2_23 = np.NAN
    
    return R2_17,R2_21,R2_23,etabar

def conlin2025(H,L,Bberm,Bbeach):
    """
    Estimate runup for composite beaches based on the parametric relation 
    by Conlin et al 2025.

    
    USAGE:
    ------
    R = conlin2025(H,L,Bberm, Bbeach)
    
    PARAMETERS:
    -----------
    H        : Offshore wave height [m]
    L        : Offshore wave length [L]
    Bberm    : Cobble berm slope
    Bbeach   : Beach slope
    
    RETURNS:
    --------
    R2     : 2% Exceedance runup [m] (Equation 4a)
    etabar : Wave induced setup [m]
    S      : Swash [m]
    
    NOTES:
    ------

    
    REFERENCES:
    -----------
    Conlin, M. P., G. Wilson, H. Bond, D. Ardag, and C. Arnowil, 2025:
        Predicting wave runup on composite beaches. Coastal Engineering, 
        199, 104743, https://doi.org/10.1016/j.coastaleng.2025.104743
    """

    # Average Slope
    Bavg = 0.5 * (Bberm + Bbeach)

    # Swash
    S = 2.99 * Bavg * H + 1.28

    # Setup
    etabar = 0.92 * Bbeach * H *((H/L)**-0.3)

    # Runup
    R2 = 1.3 * (etabar + S/2)

    return R2,etabar,S

def poate2016(H,T,B):
    """
    Estimate runup for gravel beaches based on the parametric relation 
    by Poate et al 2025.

    
    USAGE:
    ------
    R = poate2016(H,T,B)
    
    PARAMETERS:
    -----------
    H        : Offshore wave height [m]
    T        : Offshore peak wave period [s]
    B        : (gravel,cobble) Beach slope
    
    RETURNS:
    --------
    R2     : 2% Exceedance runup [m] (Equation 4a)
    
    NOTES:
    ------

    
    REFERENCES:
    -----------
    Poate, T. G., R. T. McCall, G. Masselink, 2016: A new parameterisation
        for runup on gravel beaches. Coastal Engineering, 117, 176-190,
        https://doi.org/10.1016/j.coastaleng.2016.08.003
    """

    # Solves Equation 12 for peak wave period C=0.33

    return 0.33 * np.tan(B)**0.5 * H * T

    
