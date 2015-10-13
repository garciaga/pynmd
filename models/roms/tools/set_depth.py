#!/usr/bin/python
"""
 SET_DEPTH:  Compute ROMS grid depth from vertical stretched variables

 Given a batymetry (h), free-surface (zeta) and terrain-following
 parameters, this function computes the 3D depths for the requested
 C-grid location. If the free-surface is not provided, a zero value
 is assumed resulting in unperturb depths.  This function can be
 used when generating initial conditions or climatology data for
 an application. Check the following link for details:

    https://www.myroms.org/wiki/index.php/Vertical_S-coordinate

 On Input:

    Vtransform    Vertical transformation equation:

                    Vtransform = 1,   original transformation

                    z(x,y,s,t)=Zo(x,y,s)+zeta(x,y,t)*[1+Zo(x,y,s)/h(x,y)]

                    Zo(x,y,s)=hc*s+[h(x,y)-hc]*C(s)

                    Vtransform = 2,   new transformation

                    z(x,y,s,t)=zeta(x,y,t)+[zeta(x,y,t)+h(x,y)]*Zo(x,y,s)

                    Zo(x,y,s)=[hc*s(k)+h(x,y)*C(k)]/[hc+h(x,y)]

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

    igrid         Staggered grid C-type (integer):
                    igrid=1  => density points
                    igrid=2  => streamfunction points
                    igrid=3  => u-velocity points
                    igrid=4  => v-velocity points
                    igrid=5  => w-velocity points

    h             Bottom depth, 2D array at RHO-points (m, positive),
                    h(1:Lp+1,1:Mp+1)

    zeta          Free-surface, 2D array at RHO-points (m), OPTIONAL,
                    zeta(1:Lp+1,1:Mp+1)

 On Output:

    z             Depths (m, negative), 3D array
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
from stretching import stretching

def set_depth( Vtr, Vstr, thts, thtb, hc, N, igrid, h, zeta ):
    Np      = N+1
    Lp,Mp   = np.shape(h)
    L       = Lp-1
    M       = Mp-1
    if (igrid==5):
        z   = np.empty((Lp,Mp,Np))
    else:
        z   = np.empty((Lp,Mp,N))

    hmin    = np.min(h)
    hmax    = np.max(h)

    if (igrid == 5):
        kgrid=1
    else:
        kgrid=0
    
    s,C = stretching(Vstr, thts, thtb, hc, N, kgrid)
        
    #-----------------------------------------------------------------------
    #  Average bathymetry and free-surface at requested C-grid type.
    #-----------------------------------------------------------------------

    if (igrid==1):
        hr    = h
        zetar = zeta
    elif (igrid==2):
        hp    = 0.25*(h[0:L,0:M]+h[1:Lp,0:M]+h[0:L,1:Mp]+h[1:Lp,1:Mp])
        zetap = 0.25*(zeta[0:L,0:M]+zeta[1:Lp,0:M]+zeta[0:L,1:Mp]+zeta[1:Lp,1:Mp]) 
    elif (igrid==3):
        hu    = 0.5*(h[0:L,0:Mp]+h[1:Lp,0:Mp])
        zetau = 0.5*(zeta[0:L,0:Mp]+zeta[1:Lp,0:Mp])
    elif (igrid==4):
        hv    = 0.5*(h[0:Lp,0:M]+h[0:Lp,1:Mp])
        zetav = 0.5*(zeta[0:Lp,0:M]+zeta[0:Lp,1:Mp])
    elif (igrid==5):
        hr    = h
        zetar = zeta

    #----------------------------------------------------------------------
    # Compute depths (m) at requested C-grid location.
    #----------------------------------------------------------------------
    if (Vtr == 1):
        if (igrid==1):
            for k in range (0,N):
                z0 = (s[k]-C[k])*hc + C[k]*hr
                z[:,:,k] = z0 + zetar*(1.0 + z0/hr)
        elif (igrid==2):
            for k in range (0,N):
                z0 = (s[k]-C[k])*hc + C[k]*hp
                z[:,:,k] = z0 + zetap*(1.0 + z0/hp)
        elif (igrid==3):
            for k in range (0,N):
                z0 = (s[k]-C[k])*hc + C[k]*hu
                z[:,:,k] = z0 + zetau*(1.0 + z0/hu)
        elif (igrid==4):
            for k in range (0,N):
                z0 = (s[k]-C[k])*hc + C[k]*hv
                z[:,:,k] = z0 + zetav*(1.0 + z0/hv)                      
        elif (igrid==5):
            z[:,:,0] = -hr
            for k in range (0,Np):
                z0 = (s[k]-C[k])*hc + C[k]*hr
                z[:,:,k] = z0 + zetar*(1.0 + z0/hr)
    elif (Vtr==2):
        if (igrid==1):
            for k in range (0,N):
                z0 = (hc*s[k]+C[k]*hr)/(hc+hr)
                z[:,:,k] = zetar+(zeta+hr)*z0
        elif (igrid==2):
            for k in range (0,N):
                z0 = (hc*s[k]+C[k]*hp)/(hc+hp)
                z[:,:,k] = zetap+(zetap+hp)*z0
        elif (igrid==3):
            for k in range (0,N):
                z0 = (hc*s[k]+C[k]*hu)/(hc+hu)
                z[:,:,k] = zetau+(zetau+hu)*z0
        elif (igrid==4):
            for k in range (0,N):
                z0 = (hc*s[k]+C[k]*hv)/(hc+hv)
                z[:,:,k] = zetav+(zetav+hv)*z-1
        elif (igrid==5):
            for k in range (0,Np):
                z0 = (hc*s[k]+C[k]*hr)/(hc+hr)
                z[:,:,k] = zetar+(zetar+hr)*z0

    return z
