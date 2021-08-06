# -*- coding: utf-8 -*-
"""
Functions related to CDIP wave data

Authors:
-------
Gabriel Garcia Medina
    PNNL Marine and Coastal Research Laboratory

Log of edits:
-------------
August 2021 - Created module
    Gabriel Garcia Medina
    
"""
import numpy as np

def mem_est(a1,a2,b1,b2,ener_dens):
    """
    Subroutine to make MEM estimate of the directional spectrum.
    This version reads in normalized four. coeff., so no division by a0

    PARAMETERS
    ----------
    a1, a2, b1, b2 : array
        Frequency dependent arrays of normalized directional fourier coefficient
    
    ener_dens : array
        Frequency dependent energy spectrum

    RETURNS
    -------
    mem_out : 2d array
         Directional spectrum in 5 degree directional bins
         (in "compass heading" coordinate frame).
    
    dirs : 1d array
        72 directional bands. 1=0 degress....72=355 degrees
        where the direction is the compass heading from which
        the wave energy is arriving. e.g. 270=from west, 180=
        from south etc.

    REFERENCE
    ---------
    This script translates mem.f90 into python. The original script was 
    provided by Randolph Bucciarelli (CDIP) on 21 January 2021.

    """

    # Parameters
    numDir = 72 # Hardcoded to 5 degrees
    numFreq = a1.shape[-1]

    mem_out = np.zeros((numFreq,numDir))    

    # Loop thru freq bands
    # moments are first four fourier coeff. of the directional
    # distribution, MEM uses normalized values

    for ii in range(numFreq):

        # Get MEM estimate (d) with 1 deg resolution, d is returned in
        # compass coordinates
        s,chk = mem(a1[ii],a2[ii],b1[ii],b2[ii])

        # merge into 5 deg directional bins, multiply by 0.2 to get
        # units of m^2/Hz-deg
        for jj in np.arange(5,356,5):

            # Direction index
            dInd = np.int(jj/5) 
            
            # Average over neighboring directions
            mem_out[ii,dInd] = 0.2 * ener_dens[ii] * np.sum(s[jj-2:jj+3])

        # Average first direction
        mem_out[ii,0] = 0.2 * ener_dens[ii] * (s[-2] + s[-1] + s[0] + s[1] + s[2])

    # Direction vector
    dirs = np.arange(0,356,5)

    return mem_out,dirs


def mem(a1,a2,b1,b2,begin=0,ndeg=360,res=1.0):
    """
    The function to compute the MEM method. Not to be used in stand alone mode.

    PARAMETERS
    ----------
    a1,a2,b1,b2 : array
        Normalized fourier coefficients
        (normalized ala Long [1980]) - TRIG COORDINATES

    RETURNS
    -------
    s : array
        MEM estimate of directional spectrum, 1 deg.
        resolution - COMPASS COORDINATES

    chk : float
        check factor: MEM likes to make narrow spikes for
        directional spectra if it can (it fits the a1,a2,
        b1 abd b2's exactly).  If it does this, then the
        .1 deg stepsize for making the 1 deg resolution
        estimate is too coarse.  The check factor should
        be close to 1.


    REFERENCE
    ---------
    This script translates mem.f90 into python. The original script was 
    provided by Randolph Bucciarelli (CDIP) on 21 January 2021. 

    CDIP DOCUMENTATION
    ------------------
    Maximum entropy method (MEM) for the estimation of a directional
    distribution from pitch-and-roll data. By default, covers the
    full 360 degrees with 1.0-degree resolution. (Call MEMD directly for
    settings other than the default.)
 
    (Lygre and Krogstad, JPO v16 1986: NOTE - there is a typo in the
    paper...BOTH exponetials in eq. 13 should have negative signs.
    This is fixed in Krogstad, 1991- THH has the book)
 
    """

    # Preallocate spectrum variable
    dirs = np.arange(begin,begin + ndeg)
    s = np.zeros(dirs.shape[0])

    # Degree to radian conversion factor
    dr = 0.0174533

    # Lygre & Korgstad notation (used in CDIP script but not adopted here)
    # d1 = a1
    # d2 = b1
    # d3 = a2
    # d4 = b2

    # Complex notation
    c1 = (1+0j) * a1 + (0+1j) * b1
    c2 = (1+0j) * a2 + (0+1j) * b2

    p1 = (c1 - c2 * np.conj(c1)) / (1 - np.abs(c1)**2)
    p2 = c2 - (c1 * p1)

    x = 1 - p1 * np.conj(c1) - p2 * np.conj(c2)

    # Sum over ndeg in steps, get distribution with res degree resolution
    tot = 0.0
    offset = 0.5 * (1.0 - res)

    for rn in np.arange(begin - offset, begin + ndeg + offset, res):

        a = rn * dr  # Angle in radians
        e1 = (1+0j) * np.cos(a) - (0+1j) * np.sin(a)
        e2 = (1+0j) * np.cos(2*a) - (0+1j) * np.sin(2*a)
        y = np.abs((1+0j) - p1 * e1 - p2 * e2)**2

        # put in proper 1 deg directional band
        ndir = np.int(rn)

        # Wrap to 360 degrees
        ndir = np.mod(ndir,360)

        # Normalize by 360/(step size) if not running full 360 degrees
        if ndeg != 360:
            s[ndir] = s[ndir] + (np.abs(x/y) / (360. / res))
        else:
            s[ndir] = s[ndir] + np.abs(x/y)

        tot = tot + np.abs(x/y)

    # Normalize spectrum for full 360 degree run
    if ndeg == 360:
        for ii in np.arange(0,360):
            s[ii] = s[ii] / tot        
        chk = 1
    else:
        # tot should = 360.  If directional peak is extremely narrow then
        # 1 deg resolution may be insufficient and tot .ne. 360
        chk = tot / (360./res)
    
    return s,chk
