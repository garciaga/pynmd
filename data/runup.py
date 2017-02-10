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
scipy
itertools

Internal dependencies:
----------------------
gsignal
   
"""

# Import modules
import numpy as np
import scipy as _spi
import itertools as _itertools

# Internal modules
import signal as _gsignal

# ==============================================================================
# Runup maxima
# ==============================================================================
def runup_maxima(x,ot,sten,upcross=False):
    """
    Function to identify local minima from a water surface elevation time series
    
    USAGE:
    ------
    max_ind = runup_maxima(x,ot,sten,upcross=False)
    
    PARAMETERS:
    -----------
    x       : Runup time series [m]
    ot      : Time vector, must have the same length as x and time in seconds. 
    sten    : Minimum runup period to consider [s]
    uprcross: Use zero-upcrossing to identify runup events (optional)
    
    RETURNS:
    --------
    max_ind : Vector of indices corresponding to the local maxima of the input
              time series.    
    
    NOTES:
    ------
    - All inputs should be numpy arrays
    - If upcross = True then ot and sten will not be used. x must cross zero so
      remove the setup.
    
    TODO:
    -----
    Add a function to get runup extrema with a filter for consecutive minima,
      maxima. 
    
    """
    
    if upcross:
        
        # Find the upcrossing locations
        indCross = _gsignal.zero_crossing(x)
        
        # Find the index of the maximum runup between upcrossings
        ind_max = [np.argmax(x[indCross[aa]:indCross[aa+1]]) + indCross[aa] 
                   for aa in range(indCross.shape[0]-1)]
        
        return np.array(ind_max)
    
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

# ==============================================================================
# Runup maxima
# ==============================================================================
def runupUprushSpeed(x,ot):
    """
    This function finds the runup maxima and the uprush speed from the previous
    minima and from the setup line. Zero upcrossing is used to find runup
    maxima and minima.
    
    USAGE:
    ------
    max_ind,min_ind,speedMinima,speedSetup = runup_maxima(x,ot)
    
    PARAMETERS:
    -----------
    x       : Runup time series [m]
    ot      : Time vector, must have the same length as x and time in seconds. 
    
    RETURNS:
    --------
    max_ind : Vector of indices corresponding to the local maxima of the input
              time series.
    min_ind : Vector of indices corresponding to the local maxima of the input
              time series.              
    speedMinima : uprush speed with respect to previous minima
    speedSetup  : uprush speed with respect to setup
    
    NOTES:
    ------
    - All inputs should be numpy arrays
    - x must cross zero. You can for instance remove the setup from the time
      series.
    - The first uprush speed with respect to the minima is always NAN.
    """
        
    # Find the upcrossing locations
    indCross = _gsignal.zero_crossing(x)
    
    # Find the index of the maximum runup between upcrossings
    ind_max = [np.argmax(x[indCross[aa]:indCross[aa+1]]) + indCross[aa] 
               for aa in range(indCross.shape[0]-1)]
    ind_min = [np.argmin(x[indCross[aa]:indCross[aa+1]]) + indCross[aa] 
               for aa in range(indCross.shape[0]-1)]
    
    # Runup uprush speed with respect to previous minima (check me!)
    speedMinima  = np.zeros((len(ind_max))) * np.NAN
    for aa in range(1,len(ind_max)):
        tmpIndMax = ind_max[aa]
        tmpIndMin = ind_min[aa-1]
        speedMinima[aa] = ((x[tmpIndMax] - x[tmpIndMin]) / 
                           (ot[tmpIndMax] - ot[tmpIndMin]))
    
    # Runup uprush speed with respect to setup
    speedSetup  = np.zeros((len(ind_max))) * np.NAN    
    for aa in range(len(ind_max)):
        tmpIndMax   = ind_max[aa]
        tmpIndCross = indCross[aa]
        speedSetup[aa] = x[tmpIndMax] / (ot[tmpIndMax]-ot[tmpIndCross])
    
    return np.array(ind_max),np.array(ind_min),speedMinima,speedSetup
    

#===============================================================================
# Compute mean setup
#===============================================================================
def runup_params(runup,ot,detrend=True,window=True,igcut=0.05):
    """
    
    Parameters:
    ----------
    runup        : Water surface elevation time series relative to SWL [m]
    ot           : Time stamp vector [s]
    detrend      : Linearly detrend the data
    window       : Apply a hanning window to the runup time series
    igcut        : Infragravity wave frequency cutoff
    
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
    
    """
    
    # Sampling time rate
    dt = np.mean(ot[1:] - ot[:-1])
    
    # Compute the setup
    setup = np.mean(runup)
    
    # Compute swash time series 
    if detrend:
        swash = _spi.signal.detrend(runup)
    else:
        swash = runup - setup
    
    # Compute spectrum
    if window:
        
        hwind = _spi.signal.hanning(ot.shape[0])
        ff,sf = _gsignal.psdraw(swash*hwind,dt,False)
        sf /= np.mean(hwind**2)

    else:
        
        ff,sf = _gsignal.psdraw(swash,dt,False)
        
    # Compute significant infragravity swash
    freq_ig = ff<igcut
    freq_in = ff>=igcut
    swash_ig = 4.0*(np.trapz(sf[freq_ig],ff[freq_ig])**0.5)
    swash_in = 4.0*(np.trapz(sf[freq_in],ff[freq_in])**0.5)

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


#===============================================================================
# Find unexpected runup events
#===============================================================================
def unexpected_event(x,alpha,n,m,xMean=False):
    """
    Decides if the given event is unexpected by comparing each x value to the
    previous n occurrences excluding the last m. If a given x value is greater
    than all alpha*x[n...m] returns a true value.
    
    USAGE:
    ------
    unexpected = unexpected_event(x,alpha,n,m)
    
    PARAMETERS:
    -----------
    x         : Time series of runup maxima [m]
    alpha     : Unexpected runup scale
    n         : The last n values to consider
    m         : Exclude the last m values
    xMean     : (optional) value to use as datum (e.g. mean runup)
    
    RETURNS:
    --------
    unexpected  : Boolean array for values that satisfy the input conditions
    
    NOTES:
    ------
    - n must be larger than m
    - If xMean is not given then the minimum value of x will be used as datum
    
    REFERENCES:
    -----------
    Gemmrich, J. and C. Garrett, 2008: Unexpected Waves. Journal of Physical
        Oceanography, 38, 2330 - 2336.
        
    """
    
    # Check the input data
    if x.shape[0] <= n:
        print("x must be larger than n")
        return np.array([])
    elif m >= n:
        print('n must be larger than m')
        return np.array([])
    elif n <= 0:
        print('n must be a positive integer')
        return np.array([])
        
    # Change to the correct variable type
    n     = np.int(n)                   # Should be integer  
    m     = np.int(m)                   # Should be integer
    alpha = np.double(alpha)            # Should be double
    
    # Preallocate unexpected event
    unex = [False]*x.shape[0]
    
    # Find unexpected events
    for aa in range(n,x.shape[0]):
        
        # Extract stencil to evaluate
        x_tmp = x[aa-n:aa+1]
        
        # Normalize with respect to minimum runup value
        if xMean == False:
            x_tmp = x_tmp - np.min(x_tmp)
        else:
            x_tmp - x_tmp - xMean
        
        unex[aa] = np.all(x_tmp[-1]>x_tmp[0:-1-m]*alpha)
        
    # Return unexpected event matrix
    return np.array(unex)

#==============================================================================
# ROUS Wrapper        
#==============================================================================
def rousWrapper(runupMatrix,alpha,n,m,runNumberMatrix,
                speedMinima=False,velocityMatrix=None,
                periodMinima=False,timeMatrix=None,
                randomize=False,xMean=False):
    """
    This function provides more options for computing unexpected events.
    
    PARAMETERS:
    -----------
    runupMatrix     : list of arrays of runup maxima
    alpha           : Array (see unexpected_runup)
    n               : idem
    m               : idem
    runNumberMatrix : Matrix to ID the simulations 
    speedMinima     : Array of minimum uprush speed to consider
    velocityMatrix  : Matrix like runupMatrix with uprush velocity
    periodMinima    : Minimum period to consider
    timeMatrix      : Crest-crest time matrix
    randomize       : Randomize the sequence (Boolean)
    xMean           : Legacy stuff, should be false now.
    
    NOTES:
    ------
    - I have not tested this for xMean=True
    - All inputs should be numpy arrays because n.shape[0] type of identifier
      is used.
      
    """
    
    # Check input --------------------------------------------------------------
    if periodMinima is not False and timeMatrix is None:
        print('Need timeMatrix if periodMinima is requested ')
        return
    if speedMinima is not False and velocityMatrix is None:
        print('Need velocityMatrix if speedMinima is requested ')
        return    
    
    # Preallocate matrices -----------------------------------------------------
    rous = np.zeros((len(runupMatrix),alpha.shape[0],n.shape[0],m.shape[0],
                     speedMinima.shape[0],periodMinima.shape[0]))
    numRunup   = np.zeros_like(rous)
    rousRunNum = np.zeros_like(rous,dtype=np.object)
    rousDist   = np.zeros_like(rous,dtype=np.object)

    # Loop over the stacks -----------------------------------------------------
    
    # Master      
    iterList = (np.arange(rous.shape[0]),np.arange(rous.shape[1]),
                np.arange(rous.shape[2]),np.arange(rous.shape[3]))
    
    for aa,bb,cc,dd in itertools.product(*iterList):
        
        print('  Stack: ' + np.str(aa+1) + ' of ' + np.str(rous.shape[0]))
        if randomize:        
            randInd = np.random.permutation(runupMatrix[aa].shape[0])
            runupMatrix[aa] = runupMatrix[aa][randInd]
            runNumberMatrix[aa] = runNumberMatrix[aa][randInd]
            velocityMatrix[aa] = velocityMatrix[aa][randInd]
            periodMinima = False

        try:
            # Find rous
            tmpR = unexpected_event(runupMatrix[aa],alpha[bb],
                                    n[cc],m[dd],xMean)
        except:
            for ii in range(rousRunNum.shape[4]):
                for jj in range(rousRunNum.shape[5]):
                    rousRunNum[aa,bb,cc,dd,ii,jj] = np.array([])
                    rousDist[aa,bb,cc,dd,ii,jj]   = np.array([])
            continue
        
        # No ROUS Found
        if np.size(tmpR) < 1:
            for ii in range(rousRunNum.shape[4]):
                for jj in range(rousRunNum.shape[5]):
                    rousRunNum[aa,bb,cc,dd,ii,jj] = np.array([])
                    rousDist[aa,bb,cc,dd,ii,jj]   = np.array([])
            continue
        
        # Filters on the data --------------------------------------------------
        # This should be another function at some point, to many
        # loops here.
        tmpROrig = np.copy(tmpR)
        for ii in range(speedMinima.shape[0]):
            
            for jj in range(periodMinima.shape[0]):
                
                # Load the default again
                tmpR = np.copy(tmpROrig) 
                
                if speedMinima is not False:
                    tmpR[velocityMatrix[aa]<speedMinima[ii]] = False
        
                # Time filter
                if periodMinima is not False:
                    tTM = np.zeros_like(tmpR,dtype=np.float64)
                    tTM[n[cc]:] = (timeMatrix[aa][n[cc]:] - 
                                   timeMatrix[aa][:-n[cc]])
                    tmpR[tTM<periodMinima[jj]] = False
            
                # Allocate information -----------------------------
                rous[aa,bb,cc,dd,ii,jj] = np.sum(tmpR)
                rousRunNum[aa,bb,cc,dd,ii,jj] = runNumberMatrix[aa][tmpR]
                
                # Number of runups evaluated
                numRunup[aa,bb,cc,dd,ii,jj] = (len(runupMatrix[aa]) - n[cc])
                        
                # Excursion of that runup event
                tmpRD = [runupMatrix[aa][ee] - 
                         np.max(runupMatrix[aa][ee-n[cc]:ee-m[dd]])
                         for ee in np.where(tmpR)[0]]
                rousDist[aa,bb,cc,dd,ii,jj] = tmpRD

    return rous,rousRunNum,rousDist,numRunup
