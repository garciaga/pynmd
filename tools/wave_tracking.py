"""
Functions for wave tracking

Dependencies:
-------------
numpy, scipy

Internal dependencies:
----------------------
pynmd.data.signal
pynmd.physics.waves

"""

from __future__ import division,print_function

# Import modules
import numpy as np
import scipy as spi

# Internal dependencies
import pynmd.data.signal as gsignal
import pynmd.physics.waves as _gwaves


"""
List of edits:
--------------
3 September 2015 - Bug fixes in wave_tracks. Looping over all the waves found
                   in the offshore-most array. Does not crash if the number of
                   minima is larger than the number of waves. This difference
                   should be 1 at most.
"""

#===============================================================================
# Find local minima
#===============================================================================
def local_extrema(x,ot,sten,clean=True):
    """
    Function to identify local minima from a water surface elevation time series

    USAGE:
    ------
    min_ind,max_ind = local_extrema(x,ot,sten,clean)

    PARAMETERS:
    -----------
    x       : Water surface elevation time series [m]
    ot      : Time vector, must have the same length as x and time in seconds.
    sten    : Minimum wave period to consider [s]
    clean   : (Boolean: Defaults to True) Turning this flag on will ensure that
              the local minima and maxima differ by at most one value and that
              the first identified extrema is a local minima.

    RETURNS:
    --------
    min_ind : Vector of indices corresponding to the local minima.
    max_ind : Vector of indices corresponding to the local maxima.

    NOTES:
    ------
    All inputs should be numpy arrays

    TODO:
    -----
    Add consecutive local minima (maxima) filter.

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
    # record the position of where this happens.
    dx_sign      = np.zeros_like(dx)
    dx_sign[1:]  = dx[1:] * dx[0:-1]
    ind_ext      = np.where(dx_sign<0)[0] - 1

    # Identify the local minima
    ind_min_tmp  = dx2[ind_ext]>0
    ind_min      = ind_ext[ind_min_tmp]

    # Identify the local maxima
    ind_max_tmp  = dx2[ind_ext]<0
    ind_max      = ind_ext[ind_max_tmp]

    # Small wave filter for local minima ---------------------------------------
    wave_per       = np.ones_like(ind_min) * (sten + 1)
    wave_per[1:]   = ot[ind_min][1:] - ot[ind_min][0:-1]
    ind_noise_min  = np.where(wave_per < sten)[0]

    # Filter small waves by looking at the identified minima and the previous
    # one. The point closer to a local maxima should be removed.
    if sum(ind_noise_min)>0:

        # Local minima to remove
        rmv_ind_min = []

        # Loop over short waves
        for aa in range(ind_noise_min.shape[0]):

            # Time of the consecutive local minima
            ot_mins = ot[ind_min][ind_noise_min[aa]-1:ind_noise_min[aa]+1]

            # Find the index of local maxima enclosed by the previous local
            # minima
            tmp_ind = np.where(np.logical_and(ot[ind_min][ind_noise_min[aa]-1]<
                                              ot[ind_max],
                                              ot[ind_min][ind_noise_min[aa]]>
                                              ot[ind_max]))[0]

            # Find the local maxima that encloses the previous local minima
            if tmp_ind.size == 0:
                # We have an empty array, remove both minima
                rmv_ind_min.append(ind_noise_min[aa]-1)
                rmv_ind_min.append(ind_noise_min[aa])
            elif tmp_ind[0] == 0:
                ot_maxs = ot[ind_max][tmp_ind:tmp_ind+2]
                d1 = np.abs(ot_maxs[0] - ot_mins[0])
                d2 = np.min(np.abs(ot_maxs[0:] - ot_mins[1]))
            else:
                ot_maxs = ot[ind_max][tmp_ind[0]-1:tmp_ind[0]+2]

                # Identify the local minima to remove
                d1 = np.min(np.abs(ot_maxs[0:2] - ot_mins[0]))
                d2 = np.min(np.abs(ot_maxs[1:] - ot_mins[1]))


            # Identify and flag the shortest wave
            if tmp_ind.size == 0:
                continue
            elif d1 < d2:
                rmv_ind_min.append(ind_noise_min[aa]-1)
            else:
                rmv_ind_min.append(ind_noise_min[aa])


    # Small wave filter for local maxima ---------------------------------------
    wave_per       = np.ones_like(ind_max) * (sten + 1)
    wave_per[1:]   = ot[ind_max][1:] - ot[ind_max][0:-1]
    ind_noise_max  = np.where(wave_per < sten)[0]

    # Filter small waves by looking at the identified maxima and the previous
    # one. The point closer to a local minima should be removed.
    if sum(ind_noise_max)>0:

        # Local maxima to remove
        rmv_ind_max = []

        # Loop over short waves
        for aa in range(ind_noise_max.shape[0]):

            # Time of the consecutive local maxima
            ot_maxs = ot[ind_max][ind_noise_max[aa]-1:ind_noise_max[aa]+1]

            # Find the index of local minima enclosed by the previous local
            # maxima
            tmp_ind = np.where(np.logical_and(ot[ind_max][ind_noise_max[aa]-1]<
                                              ot[ind_min],
                                              ot[ind_max][ind_noise_max[aa]]>
                                              ot[ind_min]))[0]

            # Find the local maxima that encloses the previous local minima
            if tmp_ind.size == 0:
                # We have an empty array, remove both minima
                rmv_ind_max.append(ind_noise_max[aa]-1)
                rmv_ind_max.append(ind_noise_max[aa])
            elif tmp_ind[0] == 0:
                ot_mins = ot[ind_min][tmp_ind:tmp_ind+2]
                d1 = np.abs(ot_mins[0] - ot_mins[0])
                d2 = np.min(np.abs(ot_mins[0:] - ot_mins[1]))
            else:
                ot_mins = ot[ind_min][tmp_ind[0]-1:tmp_ind[0]+2]

                # Identify the local minima to remove
                d1 = np.min(np.abs(ot_mins[0:2] - ot_maxs[0]))
                d2 = np.min(np.abs(ot_mins[1:] - ot_maxs[1]))


            # Identify and flag the shortest wave
            if tmp_ind.size == 0:
                continue
            elif d1 < d2:
                rmv_ind_max.append(ind_noise_max[aa]-1)
            else:
                rmv_ind_max.append(ind_noise_max[aa])


    # Remove short waves ------------------------------------------------------
    if sum(ind_noise_min)>0:
        ind_min = np.delete(ind_min,rmv_ind_min)

    if sum(ind_noise_max)>0:
        ind_max = np.delete(ind_max,rmv_ind_max)

    # If no waves where found then exit ----------------------------------------
    if np.size(ind_max) < 1 or np.size(ind_min) < 1:
        return ind_min,ind_max

    # Index clean up (optional) ------------------------------------------------
    if clean:

        # First extrema should be a local minima
        if ind_max[0] <= ind_min[0]:
            ind_max = ind_max[1:]

        for bb in range(10):

            # Consecutive local minima filter
            for aa in range(ind_min.shape[0]-1):

                # Need to exit prematurely if points had to be removed
                if aa >= ind_min.shape[0]-1:
                    break

                # Look for consecutive local minima
                if np.sum(np.logical_and(ind_max>ind_min[aa],
                                         ind_max<ind_min[aa+1]))>0:
                    continue
                else:
                    # Consecutive local minima without a local maxima in between
                    # identified. Keep the minimum of the two.
                    if x[ind_min[aa]]>x[ind_min[aa+1]]:
                        ind_min = np.delete(ind_min,aa)
                    else:
                        ind_min = np.delete(ind_min,aa+1)


            # Consecutive local maxima filter
            for aa in range(ind_max.shape[0]-1):

                # Need to exit prematurely if points had to be removed
                if aa >= ind_max.shape[0]-1:
                    break

                # Look for consecutive local maxima
                if np.sum(np.logical_and(ind_min>ind_max[aa],
                                         ind_min<ind_max[aa+1]))>0:
                    # A local minima was found between consecutive local maxima
                    continue
                else:
                    # Consecutive local maxima without a local minima in between
                    # identified. Keep the maximum of the two.
                    if x[ind_max[aa]]>x[ind_max[aa+1]]:
                        ind_max = np.delete(ind_max,aa+1)
                    else:
                        ind_max = np.delete(ind_max,aa)

             # First extrema should be a local minima
            if ind_max[0] <= ind_min[0]:
                ind_max = ind_max[1:]

            # Check the size of the vectors
            if np.abs(ind_max.shape[0]-ind_min.shape[0])<2:
                break

        # Message about number of iterations
        if bb == 9:
            print('maximum number of iterations reached')

    # Return minimum indices
    return ind_min,ind_max


#===============================================================================
# Compute wave heights
#===============================================================================
def wave_height(x,ot,min_ind,max_ind):
    """
    Function to compute wave height from a time series based on previously
    identified local maxima and minima.

    USAGE:
    ------
    H = wave_height(x,ot,min_ind,max_ind)

    PARAMETERS:
    -----------
    x       : Water surface elevation time series [m]
    ot      : Time vector [s]
    min_ind : Vector of indices corresponding to the local minima.
    max_ind : Vector of indices corresponding to the local maxima.

    RETURNS:
    --------
    H       : vector with wave heights [m]

    SEE ALSO:
    ---------
    local_extrema

    """

    # Make sure that the first value in the array is a local minima
    if min_ind[0]>max_ind[0]:
        max_ind = max_ind[1:]

    # The length local minima and maxima vectors must differ by at most one
    if np.abs(max_ind.shape[0]-min_ind.shape[0])>1:
        print('The length of local minima and maxima vectors must differ')
        print('  by at most one. Quitting.')
        return

    # Compute the wave height
    if min_ind.shape[0]>max_ind.shape[0]:
        H = x[max_ind] - x[min_ind[0:-1]]
    elif max_ind.shape[0]>min_ind.shape[0]:
        H = np.zeros_like(x[max_ind]) * np.NAN
    else:
        H = x[max_ind] - x[min_ind]

    return H


#===============================================================================
# Compute average time lag
#===============================================================================
def time_lag(eta,ot,lags=None):
    """
    Function to compute average time lag between the wave staffs

    USAGE:
    ------
    ot_lag = time_lag(eta,ot,lags)

    PARAMETERS:
    -----------
    eta    : Numpy array of water surface elevation and time
             eta.shape = (time,points)
    ot     : Time vector (numpy array)
    lags   : Number of lags to compute

    RETURNS:
    --------
    ot_lag : Numpy array of the same dimensions as eta with the time lagged
             arrays.

    DEPENDENCIES:
    -------------
    gsignal.cross_corr
    
    NOTES:
    ------
    eta    : Must me float not integer.

    """

    # Verify the requested lags
    if not lags:
        lags = np.floor(ot.shape[0]/2)

    # Cumulative lag time
    cum_lag_time = np.zeros((eta.shape[1],))

    # Time interval
    dt = np.mean(ot[1:] - ot[:-1])

    # Loop over points
    for aa in range(1,cum_lag_time.shape[0]):

        # Find the time lagged cross-correlation to adjust the time series
        rho,stats = gsignal.cross_corr(eta[:,aa-1],eta[:,aa],lags)

        # Identify the maximum auto correlation
        if np.max(rho) < 0.8:
            print('Warning: Correlation is less than 0.8')
            print('  aa = ' + np.str(aa))
            print('  r = ' + np.str(np.max(rho)))

        # Compute cumulative lag time
        cum_lag_time[aa] = cum_lag_time[aa-1] + stats[np.argmax(rho),0] * dt

    # Create output array based on lag time
    ot_lag = np.zeros_like(eta)
    for aa in range(cum_lag_time.shape[0]):
        ot_lag[:,aa] = ot - cum_lag_time[aa]

    # Exit function
    return ot_lag


#===============================================================================
# Compute wave tracks
#===============================================================================
def wave_tracks(local_extrema,ot_lag,twind,wh=None):
    """
    Code to track the wave troughs or crests through time

    USAGE:
    ------
    wave_tracks = wave_tracks(local_extrema,ot_lag,twind)

    PARAMETERS:
    -----------
    local_extrema : Array of local extrema indices (see local_extrema)
    ot_lag        : Array of time lags (see time_lag)
    twind         : Time window to track the waves
    wh            : (Optional) Generate a wave height vector. Use with local
                    minima only.

    RETURNS:
    --------
    wave_tracks  : Best approximation of the wave position across shore.
    wave_heights : Wave heights are given if wh passed.

    """

    # Track the waves starting from the offshore most point

    # Initialize output list
    wave_tracks = []
    if wh:
        wave_heights = []

    # The number of waves to track depends on the smallest of the local minima
    # arrays.
    #num_waves = [x.shape[0] for x in local_extrema]
    num_waves = local_extrema[0].shape[0]

    # Loop over waves in the offshore most sensor
    for aa in range(np.min(num_waves)):

        # Track the waves across the other wave gauges
        tmp_ind = np.ones((ot_lag.shape[1],)).astype(int) * -999999
        tmp_ind[0] = local_extrema[0][aa]

        if wh is not None and aa < wh[0].shape[0]:
            tmp_wh = np.ones((ot_lag.shape[1],)) * np.NAN
            tmp_wh[0] = wh[0][aa]

        # Loop over sensors
        for bb in range(1,ot_lag.shape[1]):

            # Find another index within the window we are looking at
            tmp_ot = ot_lag[:,bb-1][tmp_ind[bb-1]]

            # Flag to transect loop
            break_bb = False

            # Find time difference between the current wave and the evaluated
            # Ones
            cmax = 3
            for cc in range(cmax):
				
                # Last instrument reached, get out
                if len(local_extrema) <= (bb + cc):
                    break_bb = True
                    break
                
				# If no extrema is found try next index
                if local_extrema[bb+cc].size < 1:
                    break_bb = True
                    break
                
                tmp_dt = np.abs(tmp_ot - ot_lag[:,bb+cc][local_extrema[bb+cc]])

                # If no waves are detected break
                if tmp_dt.size < 1:
                    break_bb = True
                    break

                # Find the minimum value
                tmp_min_ind = np.argmin(tmp_dt)

                # If the time corresponding to the minimum is within the window
                # of interest then allocate the value.
                # Otherwise:
                #   1. Try one more step because sometimes there are numerical
                #      issues with the local maxima filter.
                #   2. If the maximum number of iterations is reached break the
                #      two inner loops and track the next wave.
                #      Returning NANs for the rest of the values as we have
                #      lost the ability to tack that particular wave.

                if np.min(tmp_dt) > twind*(cc+1.0):
                    break_bb = True
                    continue
                else:
                    # If cc = 0, then this is straightforward.
                    # If cc > 0, then linearly interpolate and round
                    if cc == 0:
                        tmp_ind[bb] = local_extrema[bb][tmp_min_ind]
                        if wh is not None and tmp_min_ind < wh[bb].shape[0]:
                            tmp_wh[bb] = wh[bb][tmp_min_ind]
                    else:
                        # Interpolate the index
                        raw_ind = np.interp(0,[-1,cc],[tmp_ind[bb-1],
                                            local_extrema[bb+cc][tmp_min_ind]])
                        tmp_ind[bb] = np.int(np.round(raw_ind))
                        if wh is not None and tmp_min_ind < wh[bb+cc].shape[0]:
                            raw_ind = np.interp(0,[-1,cc],[tmp_ind[bb-1],
                                                   wh[bb+cc][tmp_min_ind]])
                            if np.isfinite(raw_ind):
                                tmp_wh[bb] = np.int(np.round(raw_ind))

                    # If we are able to find a point within the stencil then
                    # break the cc loop without breaking the bb loop
                    break_bb = False
                    break

            # If the maximum spatial window has been met then exit the inner
            # loop
            if break_bb:
                break

        # Allocate in output array
        wave_tracks.append(tmp_ind)
        if wh is not None and aa < wh[0].shape[0]:
            wave_heights.append(tmp_wh)

    # Exit function
    if wh is not None:
        return wave_tracks,wave_heights
    else:
        return wave_tracks


#===============================================================================
# Function to identify bore bore capture events
#===============================================================================
def bore_bore_capture(wave_tracks,ot,twind):
    '''
    Function to flag bore-bore capture events.
    
    USAGE:
    ------
    bbc_ind,bbc_flag = wave_tracking(wave_tracks,ot,twind)
    
    INPUT:
    ------
    wave_tracks: see wave_tracks
    ot         : Time vector [s]
    twind      : Maximum spacing for bores to be considered intependent [s]  
    
    RETURNS:
    --------
    bbc_ind    : Matrix that contains the spatial index at which the current
                 wave captures the previous one
    bbc_flag   : True/False matrix that identifies whether the bore was captured    

    '''

    # Matrix will contain the spatial index at which the next wave captures
    # the current one
    bbc_ind  = np.ones(len(wave_tracks),).astype('int') * -999999
    # True if the bore was captured
    bbc_flag = np.zeros(len(wave_tracks),).astype('bool')

    for aa in range(len(wave_tracks)-1):

        # Find bore time tracks
        wtind_0 = np.argmin(wave_tracks[aa]>0) - 1
        btime_0 = ot[wave_tracks[aa][:wtind_0]]

        # Find next bore time tracks
        wtind_1 = np.argmin(wave_tracks[aa+1]>0) - 1
        btime_1 = ot[wave_tracks[aa+1][:wtind_1]]

        # If the next bore is not tracked as far as the current bore then no
        # capture occurs
        if btime_0.shape[0] > btime_1.shape[0]:
            continue

        # Find the time delta between both bores at the same locations
        dbtime = btime_1[:wtind_0] - btime_0

        # Find if there was bore-bore capture by comparing the time difference
        # within a given time window.
        dtind = dbtime < twind

        # Save the index of the first location where the two bores are closer
        # than the specfied time window
        if np.sum(dtind) > 0:
            tmpind        = np.argmax(dtind)
            bbc_flag[aa]  = True
            bbc_ind[aa+1] =  tmpind


    return bbc_ind,bbc_flag


#===============================================================================
# Compute wave tracks
#===============================================================================
def wave_tracks_predictor(local_maxima,wave_height,ot,xInst,x,h,wp=None):
    """
    Code to track the wave crests througout a linear instrument array using
    bathymetric data and linear wave theory as predictors

    USAGE:
    ------
    wave_tracks = wave_tracks(local_extrema,ot_lag,twind)

    PARAMETERS:
    -----------
    local_maxima  : Array of local maxima indices (see local_extrema)
    wave_height   : Array of wave heights (see wave_height)
    ot            : Time vector [s]
    xInst         : Instrument easting
    x             : Easting of topography/bathymetry [m]
                    Must increase landward
    h             : Water depth [m]
    wp            : (optional) Representative wave period [s]

    RETURNS:
    --------
    wave_tracks  : Best approximation of the wave position across shore.
    wave_ind     : Incides of the wave tracks for slicing local_maxima and 
                   wave_height
    
    NOTES:
    ------
    1. Do not include NANs in the x and h variables. 
    2. xInst should be within the limits of x
    3. If waves were not tracked an index of -999999 will be used.
    4. If wp is not provided, the shallow water wave celerity will be used 
       as predictor.
    """
    
    # Compute the wave celerity based on linear wave theory if the wave period
    # is provided, otherwise the shallow water celerity is computed
    if wp:
        cel = np.zeros_like(h)
        for aa in range(cel.shape[0]):
            tmpCel = _gwaves.celerity(wp,h[aa])
            cel[aa] = tmpCel[0]
            del tmpCel
        
    else:
        cel = (9.81*h)**0.5
        
    # Compute time of travel
    timeTravel = np.zeros_like(x)
    timeTravel[1:] = np.cumsum(2.0/(cel[1:]+cel[:-1])*(x[1:] - x[:-1]))
    
    # Interpolate the location of instruments
    timeInt    = spi.interpolate.interp1d(x,timeTravel)
    timeOffset = timeInt(xInst)
    
    # Normalize the time travel vector
    timeOffset -= timeOffset[0]

    # Find the water depth at the instruments 
    depthInt = spi.interpolate.interp1d(x,h)
    depth    = depthInt(xInst) 
    
    # Expected average velocity between instrument locations
    # Forward differences of course.
    meanCel     = np.zeros_like(xInst) * np.NAN
    meanCel[1:] = np.abs(np.diff(xInst)) / np.diff(timeOffset) 
    
    # Loop over local maxima and preallocate the variables
    wave_tracks = np.ones((len(local_maxima[0]),xInst.shape[0]),
                          dtype=np.int64) * -999999
    wave_ind = np.copy(wave_tracks)                      
    
    for aa in range(wave_tracks.shape[0]):
     
        # Preallocate temporary variables
        # tmp_ind will be allocated into wave_tracks 
        # tmp_wave_ind will be allocated into wave_ind
        tmp_ind = np.ones((wave_tracks.shape[1],)).astype(int) * -999999
        tmp_wave_ind = np.copy(tmp_ind)
        
        tmp_ind[0] = local_maxima[0][aa]
        tmp_wave_ind[0] = aa
        
        # Loop over sensors
        for bb in range(1,tmp_ind.shape[0]):
            
            # Previous wave time
            tmp_ot = ot[tmp_ind[bb-1]]
    
            # Find the next wave based on a celerity estimate
            tmp_dt = ot[local_maxima[bb]] - tmp_ot
            tmpCel = np.abs(xInst[bb-1] - xInst[bb]) / tmp_dt
            tmp_min_ind = np.argmin(np.abs(tmpCel - meanCel[bb]))
            
            # Check if we are tracking well into the past based on half of the
            # time series length
            futureFlag = ((ot[local_maxima[bb][tmp_min_ind]] - tmp_ot) < 
                          (ot[0] - ot[-1])/2.0)
            if futureFlag:
                break
            
            # Wave breaking flag based on saturated breaking (H=kh)
            dWH = ((wave_height[bb-1][tmp_wave_ind[bb-1]] - wave_height[bb]) /
                   wave_height[bb-1][tmp_wave_ind[bb-1]])
            maxdWH = (depth[bb-1] - depth[bb])/depth[bb-1]
            
            # Three checks here to consider the previous wave as the correct one
            # All must be true
            # 1. If the wave moves slower than the shallow water celerity
            # 2. The new trajectory does not exceed amplitude dispersion
            # 3. The wave is not much smaller than the one it maps to  
            
            # Find amplitude dispersion effect
            ampDisp = (9.81*wave_height[bb-1][tmp_wave_ind[bb-1]])**0.5   
            if (tmpCel[tmp_min_ind] < meanCel[bb] and
                tmpCel[tmp_min_ind-1] < (meanCel[bb]+ampDisp) and
                dWH[tmp_min_ind-1]<=dWH[tmp_min_ind]):
                tmp_min_ind -= 1
            
            # Check if we are tracking into the past      
            futureFlag = (ot[local_maxima[bb][tmp_min_ind]] - tmp_ot) < 0
            if futureFlag:
                tmp_min_ind += 1
                
            # Stand alone wave breaking flag
            if (dWH[tmp_min_ind] > maxdWH and 
                dWH[tmp_min_ind-1] < dWH[tmp_min_ind]):

                # Pick the previous wave if the speed makes sense
                dCel = (meanCel[bb] - tmpCel[tmp_min_ind-1]) / meanCel[bb]
                
                # Check if we are tracking into the past            
                futureFlag = ot[local_maxima[bb][tmp_min_ind-1]] > tmp_ot 
                
                # Negative means faster (need to do something objective here)
                if dCel > -0.5 and futureFlag:
                    tmp_min_ind -= 1       
            
            # Finally check for wave crossing (brute force bore capture)
            if aa > 0:
                if local_maxima[bb][tmp_min_ind] < wave_tracks[aa-1,bb]:
                    tmp_min_ind = wave_ind[aa-1,bb]                
                    
            # Allocate in arrays
            tmp_wave_ind[bb] = tmp_min_ind
            tmp_ind[bb] = local_maxima[bb][tmp_min_ind]
            
        # Allocate wave data        
        wave_tracks[aa,:] = np.copy(tmp_ind)
        wave_ind[aa,:] = np.copy(tmp_wave_ind)

    # Get out
    return wave_tracks,wave_ind
