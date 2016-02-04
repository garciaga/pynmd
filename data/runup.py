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



#===============================================================================
# Compute mean setup
#===============================================================================
def runup_params(runup,ot):
    """
    
    Parameters:
    ----------
    runup        : Water surface elevation time series relative to SWL [m]
    ot           : Time stamp vector [s]
    
    Output:
    -------
    Dictionary containing    
    setup        : Mean water surface elevation [m]
    r2_combined  : 2% runup exceedence value [m]
    r2_cdf       : 2% runup exceedence value computed from runup CDF [m]
    r1_cdf       : 1% runup exceedence value computed from runup CDF [m]
    ig           : Significant infragravity swash elevation [m]
    in           : Significant incident swash elevation [m]
    r_max        : Maximum runup [m]
    r_var        : Runup variance [m2]          
               
    Notes:
    ------
    r2_combined = 1.1*(setup + 0.5*((ig**2 + in**2)**0.5))
    
    See also:
    ---------
    runup
    
    TODO:
    -----
    r2_max_cdf   : 2% runup exceedence value computed from runup maxima CDF [m]
    r1_max_cdf   : 1% runup exceedence value computed from runup maxima CDF [m]   

    """
    
    # Compute the setup
    setup = np.mean(runup)
    
    # Compute swash time series 
    swash = runup - setup
    
    # Compute spectrum
    freq = np.fft.fftshift(np.fft.fftfreq(ot.shape[0],ot[1]-ot[0]))
    Nt = freq.shape[0]
    if np.mod(Nt,2) == 1:
        zero_ind = np.int((Nt - 1.0)/2.0)
    else:
        zero_ind = np.int(Nt/2.0)  
    freq_amp = freq[zero_ind:]
    
    # Compute spectrum
    ff = np.fft.fftshift(np.fft.fft(swash))
    #sf = (2.0/Nt*(ff[zero_ind:].real**2 + ff[zero_ind:].imag**2)**0.5)
    sf = (ff[zero_ind:].real**2 + ff[zero_ind:].imag**2)/Nt*(ot[1]-ot[0])
    
    # Compute significant infragravity swash
    freq_ig = freq_amp<0.05
    freq_in = freq_amp>=0.05
    swash_ig = 4.0*(np.trapz(sf[freq_ig],freq_amp[freq_ig])**0.5)
    swash_in = 4.0*(np.trapz(sf[freq_in],freq_amp[freq_in])**0.5)

    # Compute R2% real and from formula
    r2_combined = (setup + 0.5*((swash_ig**2 + swash_in**2)**0.5))*1.1
    
    # Compute R2% from the cumulative distribution function
    r2_ind = np.int(np.floor(0.98*runup.shape[0]))
    r2_cdf = np.sort(runup)[r2_ind]
    
    # Compute R1% from the cumulative distribution function
    r1_ind = np.int(np.floor(0.99*runup.shape[0]))
    r1_cdf = np.sort(runup)[r1_ind]    
    
    # Compute the R2% from the runup maxima cumulative distribution function
    
    # Generate output
    return {'setup':setup, 'ig':swash_ig,'in':swash_in,
            'r2_combined':r2_combined,'r2_cdf':r2_cdf,
            'r1_cdf':r1_cdf,'r_max':runup.max(),'r_var':np.var(runup)}

