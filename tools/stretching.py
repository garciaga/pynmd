#!/usr/bin/python
"""

 STRETCHING:  Compute ROMS vertical coordinate stretching function

 [s,C]=stretching(Vstretching, theta_s, theta_b, hc, N, kgrid, report)

 Given vertical terrain-following vertical stretching parameters, this
 routine computes the vertical stretching function used in ROMS vertical
 coordinate transformation. Check the following link for details:

    https://www.myroms.org/wiki/index.php/Vertical_S-coordinate

 On Input:

    Vstretching   Vertical stretching function:
                    Vstretching = 1,  original (Song and Haidvogel, 1994)
                    Vstretching = 2,  A. Shchepetkin (UCLA-ROMS, 2005)
                    Vstretching = 3,  R. Geyer BBL refinement
                    Vstretching = 4,  A. Shchepetkin (UCLA-ROMS, 2010)
    theta_s       S-coordinate surface control parameter (scalar)
    theta_b       S-coordinate bottom control parameter (scalar)
    hc            Width (m) of surface or bottom boundary layer in which
                    higher vertical resolution is required during
                    stretching (scalar)
    N             Number of vertical levels (scalar)
    kgrid         Depth grid type logical switch:
                    kgrid = 0,        function at vertical RHO-points
                    kgrid = 1,        function at vertical W-points
 On Output:

    s             S-coordinate independent variable, [-1 <= s <= 0] at
                    vertical RHO- or W-points (vector)
    C             Nondimensional, monotonic, vertical stretching function,
                    C(s), 1D array, [-1 <= C(s) <= 0]

#--------------------------------------------------------------------#
  NOTES:
          - Adapted from Hernan Arango's Matlab script.
#--------------------------------------------------------------------#
"""
__author__ = "Cigdem Akan"
__email__ = "akanc@ucla.edu"
__group__ = "UCLA ROMS"

import numpy as np
import pylab as pl

def stretching(Vstr, thts, thtb, hc, N, kgrid):
    s=[]
    C=[]

    Np=N+1
    #-----------------------------------------------------------------
    # Compute ROMS S-coordinates vertical stretching function
    #-----------------------------------------------------------------
    
    # Original vertical stretching function (Song and Haidvogel, 1994).
    if (Vstr == 1):
        ds = 1.0/N

        if (kgrid == 1):
            Nlev = Np
            lev  = np.linspace(0.0,Nlev,Nlev,endpoint=True)
            s    = (lev-Nlev)*ds
        else:
            Nlev = N
            lev  = np.linspace(1.0,Nlev,Nlev,endpoint=True)-0.5
            s    = (lev-Nlev)*ds
        
        if (thts > 0):
            Ptheta = np.sinh(thts*s)/np.sinh(thts)
            Rtheta = np.tanh(thts*(s+0.5))/(2.0*np.tanh(0.5*thts))-0.5
            C      = (1.0-thtb)*Ptheta+thtb*Rtheta
        else:
            C=s

    # A. Shchepetkin (UCLA-ROMS, 2005) vertical stretching function.
    if (Vstr==2):        
        alfa = 1.0
        beta = 1.0
        ds   = 1.0/N

        if (kgrid == 1):
            Nlev = Np
            lev  = np.linspace(0.0,Nlev,Nlev,endpoint=True)
            s    = (lev-Nlev)*ds
        else:
            Nlev = N
            lev  = np.linspace(1.0,Nlev,Nlev,endpoint=True)-0.5
            s    = (lev-Nlev)*ds
        
        if (thts > 0):
            Csur = (1.0-np.cosh(thts*s))/(np.cosh(thts)-1.0)
            if (thtb > 0):
                Cbot   = -1.0+np.sinh(thtb*(s+1.0))/np.sinh(thtb)
                weigth = (s+1.0)**alfa*(1.0+(alfa/beta)*(1.0-(s+1.0)**beta))
                C      = weigth*Csur+(1.0-weigth)*Cbot
            else:
                C=Csur
        else:
            C=s

    # R. Geyer BBL vertical stretching function.
    if (Vstr==3):
        ds   = 1.0/N
       
        if (kgrid == 1):
            Nlev = Np
            lev  = np.linspace(0.0,Nlev,Nlev,endpoint=True)
            s    = (lev-Nlev)*ds
        else:
            Nlev = N
            lev  = np.linspace(1.0,Nlev,Nlev,endpoint=True)-0.5
            s    = (lev-Nlev)*ds

        if (thts > 0):
            exp_s = thts   # surface stretching exponent
            exp_b = thtb   # bottom  stretching exponent
            alpha = 3      # scale factor for all hyperbolic functions
            Cbot  = np.log(np.cosh(alpha*(s+1.0)**exp_b))/np.log(np.cosh(alpha))-1.0
            Csur  = -np.log(cosh(alpha*abs(s)**exp_s))/log(cosh(alpha))
            weight= (1-np.tanh( alpha*(s+0.5)))/2.0
            C     = weight*Cbot+(1.0-weight)*Csur
        else:
            C=s

    # A. Shchepetkin (UCLA-ROMS, 2010) double vertical stretching function
    # with bottom refinement
    if (Vstr == 4):
        ds   = 1.0/N

        if (kgrid == 1):
            Nlev = Np
            lev  = np.linspace(0.0,Nlev,Nlev,endpoint=True)
            s    = (lev-Nlev)*ds
        else:
            Nlev = N
            lev  = np.linspace(1.0,Nlev,Nlev,endpoint=True)-0.5
            s    = (lev-Nlev)*ds

        if (thts > 0):
            Csur = (1.0-np.cosh(thts*s))/(np.cosh(thts)-1.0)
        else:
            Csur = -s**2

        if (thtb > 0):
            Cbot = (np.exp(thtb*Csur)-1.0)/(1.0-np.exp(-thtb))
            C    = Cbot
        else:
            C    = Csur

    return (s,C)

