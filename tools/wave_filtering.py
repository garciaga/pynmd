# -*- coding: utf-8 -*-
"""
Functions for wave filtering

Dependencies:
-------------
numpy

Internal Dependencies:
----------------------
physics.waves

Created on Mon Apr  4 21:49:54 2016

@author: ggarcia
"""

from __future__ import division,print_function

# Import modules
import numpy as np
import pynmd.physics.waves as _gwaves

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


#===============================================================================
# Separation of incident and reflected waves
#===============================================================================
def sep_battjes(eta,ot,h,x,sten,vd=True,verbose=False):
    """    
    PARAMETERS:
    -----------
    eta     : long wave water surface elevation matrix (time,points)
    ot      : time vector [s]
    h       : water depth [m]
    x       : point location [m]
    sten    : Point averaging stencil
    vd      : Use the extension proposed by van Dongeren et al. (2007). 
              Defaults to True.
    verbose : (optional) True for some printed info
    
    RETURNS:
    --------
    etaInc  : incident wave time series [m]
    etaRef  : reflected wave time series [m]
    
    NOTES:
    ------
    1. eta has to be sampled at the same rate for all points
    2. numpy arrays are assumed
    3. x axis should point landward (i.e. eta_inc travel with increasing x)
    4. van Dongeren et al. (2007) extended the separation method of 
       Battjes et al. (2004) by considering wave shoaling and phase speed
       effects.
    
    REFERENCES:
    -----------
    Battjes, J. A., H. J. Bakkenes, T. T. Janssen, and A. R. van Dongeren, 2004:
        Shoaling of subharmonic gravity waves. Journal of Geophysical Research,
        109, C02009, doi:10.1029/2003JC001863.
    van Dongeren, A., J. Battjes, T. Janssen, J. van Noorloos, K. Steenhauer, 
        G. Steenbergen, and A. Reniers, 2007: Shoaling and shoreline 
        dissipation of low-frequency waves. Journal of Geophysical Research,
        112, C02011, doi:10.1029/2006JC003701.
    """
    
    # Compute general Fourier parameters
    N  = ot.shape[0]
    dt = np.mean(ot[1:] - ot[:-1])
    fN = 1.0/(2*dt)
        
    # Compute group velocities for incoming waves and shallow water propagation
    # velocities for reflected waves. 
    freq = np.fft.fftfreq(N,dt)
    k    = np.zeros_like(eta)
    c    = np.zeros_like(eta)
    cg   = np.zeros_like(eta)
    clin = np.zeros_like(eta)
    
    # Loop over frequencies (may take advantage of the symmetry properties)
    if verbose:
        print('Long wave incident and reflected wave separation')
        if vd:
            print('  Using the van Dongeren et al. (2007) methodology')
        else:
            print('  Using the Battjes et al. (2004) methodology')
        print('Computing phase and group velocities based on linear theory')
        print('  Be patient ...')
        
    # Time loop        
    for aa in range(k.shape[0]):
            
        # Loop over points in the array
        for bb in range(k.shape[1]):
            
            # Compute wave number
            if np.isclose(freq[aa],0):
                k[aa,bb]    = 0.0
                c[aa,bb]    = 0.0
                cg[aa,bb]   = 0.0
                clin[aa,bb] = 0.0
            else:
                k[aa,bb] = _gwaves.dispersion(np.abs(1.0/freq[aa]),h[bb])
                # Compute the wave celerity and group velocity
                c[aa,bb],n,cg[aa,bb] = _gwaves.celerity(np.abs(1.0/freq[aa]),
                                                       h[bb])
                # Compute celerity based on shallow water equations
                clin[aa,bb] = np.sqrt(9.81*h[bb])
    
    
    # Equation numbers based on van Dongeren et al 2007
    # First working on equation A1 ---------------------------------------------
    
    # Compute the Fourier transform and arrange the frequencies
    z_mp = np.fft.fft(eta,axis=0)
    freq = np.fft.fftfreq(N,dt)
    
    # The reflected waves will be found by looking at finite stencils. 
    eta_i_b = np.zeros_like(eta) * np.NAN
    eta_r_b = np.zeros_like(eta) * np.NAN
    
    # Compute where to start the counter based on the stencil and the dx
    half_sten = int((sten - 1)/2)
    ind_min = half_sten
    ind_max = eta.shape[1] - ind_min
    
    # Loop over array
    for ii in range(ind_min,ind_max):
    
        # Progress update
        if verbose:
            print('  Point ' + np.str(ii-ind_min+1) + ' of ' + 
                  np.str(ind_max-ind_min))
    
        # Work on a subset of eta now
        subset_ind = np.arange(ii-half_sten,ii+half_sten+1,1)
        cg_work    = cg[:,subset_ind].copy()
        clin_work  = clin[:,subset_ind].copy()
        h_work     = h[subset_ind].copy()
        x_work     = x[subset_ind].copy()
        z_mp_work  = z_mp[:,subset_ind].copy()
        
        # Equation A2 ----------------------------------------------------------
        # z_mp = q_i_mpn z_i_mRn + Q_r_mpn Z_r_mRn + e_mpn
        # The Fourier transform is considered as the sum of an incoming wave
        # component and a reflected wave component with an error term. 
        # The next step is to assume that the amplitude and phase factors are 
        # close to the values produced by linear wave theory. 
            
        # Estimate the factors Q (Equations A4,A5,A6)
        freq_mat = np.repeat(np.expand_dims(freq,axis=-1),
                             h_work.shape[0],axis=-1)
        q_i_mp   = np.zeros_like(freq_mat)*1j
        q_r_mp   = np.zeros_like(freq_mat)*1j
        for aa in range(x_work.shape[0]):
            q_i_mp[:,aa] = np.exp(-1j*np.trapz(2.0*np.pi*freq_mat[:,0:aa+1]/
                                               cg_work[:,0:aa+1],
                                               x_work[0:aa+1],axis=-1))
            q_r_mp[:,aa] = np.exp(1j*np.trapz(2.0*np.pi*freq_mat[:,0:aa+1]/
                                              clin_work[:,0:aa+1],
                                              x_work[0:aa+1],axis=-1))
            
    
        # Solve equation A3 to estimate z_i_mR and z_r_mR by using the least
        # squares method since the system is overdetermined --------------------
        z_i_mR = np.zeros((ot.shape[0],)) * 1j
        z_r_mR = np.zeros_like(z_i_mR)
        
        for aa in range(ot.shape[0]):
            # Put system of equations in the form Ax = B
            q_mat = np.vstack((q_i_mp[aa,:],q_r_mp[aa,:])).T
            if np.sum(np.isnan(q_mat)) > 0:
                continue
                z_i_mR[aa] = 0.0
                z_r_mR[aa] = 0.0
            else:
                # Use least squares method to find solution
                z_i_mR[aa],z_r_mR[aa] = np.linalg.lstsq(q_mat,
                                                        z_mp_work[aa,:])[0]
        
        
        # So far the methodology applied is the one described in Battjes (2004). 
        # If the user selected true then the modifications of 
        # van Dongeren et al. (2007) will be implemented. 
        if vd:
            
            # Correction for incident waves ------------------------------------
            # Phase correction equation (A7) 
            q_i_mp *= z_mp_work / np.repeat(np.expand_dims(np.abs(z_i_mR),
                                                           axis=-1),
                                            z_mp_work.shape[1],axis=-1)

            # Solve equation (A3)
            z_i_mR,z_r_mR = _vanDongeren2007A3(q_i_mp,q_r_mp,z_mp_work)
            
            # Correction for amplitude (A8)
            q_i_mp *= (np.abs(z_mp_work) / 
                       np.repeat(np.expand_dims(np.abs(z_i_mR),axis=-1),
                                 z_mp_work.shape[1],axis=-1))
                
            # Solve equation (A3)
            z_i_mR,z_r_mR = _vanDongeren2007A3(q_i_mp,q_r_mp,z_mp_work)                     
    
            # Correction for reflected waves -----------------------------------
            # Phase correction equation (A7) 
            q_r_mp *= (z_mp_work / 
                       np.repeat(np.expand_dims(np.abs(z_r_mR),axis=-1),
                                 z_mp_work.shape[1],axis=-1))
                                                    
            # Solve equation (A3)
            z_i_mR,z_r_mR = _vanDongeren2007A3(q_i_mp,q_r_mp,z_mp_work)
            
            # Correction for amplitude (A8)
            q_r_mp *= (np.abs(z_mp_work) / 
                       np.repeat(np.expand_dims(np.abs(z_r_mR),axis=-1),
                                 z_mp_work.shape[1],axis=-1)) 
                                    
            # Solve equation (A3)
            z_i_mR,z_r_mR = _vanDongeren2007A3(q_i_mp,q_r_mp,z_mp_work)
                            
        
        # Compute the inverse Fourier transform to get the time series of the 
        # reflected and incident waves.
        eta_i_b[:,ii] = np.fft.ifft(z_i_mR).real
        eta_r_b[:,ii] = np.fft.ifft(z_r_mR).real
        
    # End of function
    return eta_i_b,eta_r_b


def _vanDongeren2007A3(q_i_mp,q_r_mp,z_mp_work):
    """
    Equation A3 from van Dongeren et al. 2007. See sep_battjes for more info.
    """

    # Equation (A3)
    z_i_mR = np.zeros((q_i_mp.shape[0],)) * 1j
    z_r_mR = np.zeros_like(z_i_mR)
    
    for aa in range(z_i_mR.shape[0]):
        # Put system of equations in the form Ax = B
        q_mat = np.vstack((q_i_mp[aa,:],q_r_mp[aa,:])).T
        if np.sum(np.isnan(q_mat)) > 0:
            continue
            z_i_mR[aa] = 0.0
            z_r_mR[aa] = 0.0
        else:
            # Use least squares method to find solution
            z_i_mR[aa],z_r_mR[aa] = np.linalg.lstsq(q_mat,z_mp_work[aa,:])[0]
                                    
    return z_i_mR,z_r_mR
    