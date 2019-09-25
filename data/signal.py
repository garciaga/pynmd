# -*- coding: utf-8 -*-
"""
A series of tools for data analysis

Authors:
-------
Gabriel Garcia Medina
    Nearshore Modeling Group
    ggarcia@coas.oregonstate.edu
Saeed Moghimi

Log of edits:
-------------
April 2014 - Created module
    Gabriel Garcia Medina
17 September 2015
    Gabriel Garcia Medina

Dependencies:
-------------
    numpy
    scipy
    sys

Internal dependencies:
----------------------
    none

"""

from __future__ import division, print_function

__author__ = 'Gabriel Garcia Medina'
__email__ = 'ggarcia@coas.oregonstate.edu'
__group__ = 'Nearshore Modeling Group (Oregon State University)'

# Import modules
import numpy as np
import sys
import scipy as spi
import scipy.stats
import scipy.signal

# =============================================================================
# Cross correlation Function
# =============================================================================
def cross_corr(x,y,lags,norma=1.0):
    """

    CROSS_CORR

    This function computes the cross correlation between to input time series.

    Input
    -----
       - x, y  : are time series of the same length. Missing values should be
                 flagged as np.nan.
       - lags  : number of lags desired
       - norma : normalization, takes 0 or 1. If not provided the function
                 defaults to 1 (based estimate, N).

    Results
    -------
       - rho   : cross(auto)-correlation coefficient (aka R)
       - stats : lag, cross(auto)-covariance, x standard deviation,
                 y standard deviation, number of observations

    Notes
    -----
       - For auto-correlation use x = y
       - Positive lags indicate y follows x (i.e. x leads y)
       - Use at your own risk

    """

    # Quick data check
    lags = np.int(lags)
    norma = np.float(norma)

    if (norma != 0.0) and (norma != 1.0):
        print("A normalization value of " + np.str(norma) +
              " is not supported. Will compute with norma = 1.0\n")
        norma = 1.0

    # Preallocate variables
    rho = np.zeros((2*lags+1,))
    stats = np.zeros((2*lags+1,5))

    # loop over lags
    for aa in range(-lags,lags+1,1):

        # Get data subset
        if aa > 0:
            # Positive lags
            ytmp = np.array(y[aa:])
            xtmp = np.array(x[0:-aa])
        elif aa < 0:
            # Negative lags
            ytmp = np.array(y[0:aa])
            xtmp = np.array(x[-aa:])
        else:
            xtmp = np.array(x)
            ytmp = np.array(y)


        # Make sure NaNs are consistent in both subsets
        xtmp[np.isnan(ytmp)] = np.nan
        ytmp[np.isnan(xtmp)] = np.nan

        # Compute mean values
        nn = np.sum(np.isfinite(xtmp) - 1.0 + norma)
        xtmpmean = np.nansum(xtmp)/nn
        ytmpmean = np.nansum(ytmp)/nn

        # Cross-covariance
        crosscov = np.nansum((xtmp - xtmpmean) * (ytmp - ytmpmean))/nn
        rho[aa+lags] = crosscov / np.nanstd(xtmp) / np.nanstd(ytmp)

        # Compute and save statistics
        stats[aa + lags,0] = aa
        stats[aa + lags,1] = crosscov
        stats[aa + lags,2] = np.nanstd(xtmp)
        stats[aa + lags,3] = np.nanstd(ytmp)
        stats[aa + lags,4] = np.sum(np.isfinite(xtmp))

    return rho,stats


#===============================================================================
# Variance spectrum
#===============================================================================
def psdraw(ts,dt=1,demean=False):
    """

    [freq,Sf] = psdraw(ts,dt,demean)

    Compute the variance spectrum

    PARAMETERS:
    -----------
    ts    : time series
    dt    : sampling rate (in time domain)
    demean: Remove mean

    RETURNS:
    --------
    freq  : Spectral frequencies (Positive Fourier frequencies)
    Sf    : Variance spectrum

    NOTES:
    ------
    var(ts) = integrate(freq,Sf)
        >>> np.var(ts)
        >>> spi.integrate.trapz(Sf,freq)

    Nyquist frequency = (2*dt)**-1
    
    Frequency resolution = (N*dt)**-1
    """

    # Remove mean
    if demean:
        ts -= ts.mean()

    # Compute record length
    N = np.int(ts.shape[0])
    N2 = np.int(N/2)
    N12 = np.int((N+1)/2)

    # Another way that is equivalent for reference
    """
    # Compute frequency information
    freq = np.fft.fftshift(np.fft.fftfreq(N,dt))

    if np.mod(N,2) == 1:
        zero_ind = np.int((N - 1.0)/2.0)
    else:
        zero_ind = np.int(N/2.0)
    freq_amp = freq[zero_ind:]

    # Compute spectrum
    ff = np.fft.fftshift(np.fft.fft(ts))
    sf = 2.0/N*(ff[zero_ind:].real**2 + ff[zero_ind:].imag**2)


    """

    # Compute fourier frequencies
    fj = np.fft.fftfreq(N,dt)

    # Compute power spectral density (Cooley-Tukey Method)
    yf = np.fft.fft(ts)/N
    psd = N*dt*yf*np.conjugate(yf)

    # One sided psd from dft
    if np.mod(N,2) == 0:
       sf = np.concatenate((np.array([psd[0]]),2.0*psd[1:N2],
                            np.array([psd[N2]])))
       freq_amp = np.abs(np.concatenate((np.array([fj[0]]),fj[1:N2],
                                         np.array([fj[N2]]))))
    else:
       sf = np.concatenate((np.array([psd[0]]),2.0*psd[1:N12]))
       freq_amp = np.abs(np.concatenate((np.array([fj[0]]),fj[1:N12])))


    # End of function
    return freq_amp,sf.real


# =============================================================================
# Confidence levels
# =============================================================================
def psd_ci(spec,cl,dof):
    """
    Compute condidence intervals for Power Spectral Densities

    USAGE:
    ------
    ci = psd_ci(spec,cl,dof)

    PARAMETERS:
    -----------
    spec   : PSD
    cl     : decimal confidence level (e.g cl=0.95 for 95% confidence)
    dof    : degrees of freedom

    RETURNS:
    --------
    ci     : Confidence intervals for spec

    NOTES:
    ------
    Increase in Degrees of Freedom can be achieved by band averaging or
      ensemble averaging. For example averaging over 5 bands results in
      10 degrees of freedom.
    EXAMPLE:
    --------
    >>> freq,spec = psdraw(yt,dt)
    >>> ci = psd_ci(spec,0.95,2)
    >>> pl.errorbar(freq,spec,yerr=ci,fmt='o',color='k')

    """

    # Compute confidence level
    alpha = 1.0 - cl

    p95_upper = spi.stats.chi2.ppf(alpha/2.0,dof)
    p95_lower = spi.stats.chi2.ppf(1.0-alpha/2.0,dof)

    conf_lims_lower = spec*dof/p95_lower
    conf_lims_upper = spec*dof/p95_upper
    conf_lims = np.array([spec-conf_lims_lower,conf_lims_upper-spec])

    return conf_lims

#===============================================================================
# Compute psd with the pre-whitening and post-colouring technique
#===============================================================================
def psd_pw_pc(yt,dt=1,stencil=False):
    """
    Compute the power spectral density but the time series will be
    pre-whitened before computing the PSD and post-coloured after the fact.

    USAGE:
    ------
    [freq,Sf] = psd_pw_pc(ts,dt,stencil)

    PARAMETERS:
    -----------
    ts        : time series
    dt        : sampling rate (in time domain)
    stencil   : (Optional) Band averaging. If desired must pass stencil.

    RETURNS:
    -------
    freq      : Spectral frequencies (Positive Fourier frequencies)
    Sf        : Variance spectrum

    NOTES:
    ------
    This code depends on:
        psdraw
        band_averaging

    This procedure is good for time series that have a ``red'' spectrum.

    PROCEDURE:
    ----------
    1. Approximate the derivative using the first difference
       >>> xt = yt[1::] - yt[0:-1]
       >>> xt -= xt.mean()
    2. Compute the PSD of the first differenced time series
       >>> freq,Sx = psdraw(xt)
    3. Band averaging (if requested)
       >>> Sx = band_averaging(Sx,freq,stencil)
    4. Post-colour the timer series
       >>> Sy = Sx/(4*sin(pi*freq*dt)**2)
    """

    # Pre-whiten the spectrum
    # Compute first difference
    xt = yt[1::] - yt[0:-1]

    # Remove mean from the derivative
    xt -= xt.mean()

    # Compute power spectral density
    freq,Sx = psdraw(xt,dt)

    # Band averaging if requested
    if stencil:
        freq,Sx = band_averaging(Sx,freq,stencil)

    # Post-colour the spectrum
    Sf = Sx/(4.0 * np.sin(np.pi*freq*dt)**2)

    # End of function
    return freq,Sf

#===============================================================================
# Squared coherence spectrum
#===============================================================================
def squared_coherence(x,y,stencil,dt=1,cl=0.95):
    """

    freq,g2,g2crit,sxy_conj,sxy,sxx,syy,pxy,pcl = \
        squared_coherence(x,y,stencil,dt,cl)

    PARAMETERS:
    -----------
    x,y      : time series of same size
    stencil  : Stencil over which to band average
    dt       : sampling rate (in time domain)
    cl       : confidence level (Optional: Defaults to 0.95)

    RETURNS:
    --------
    freq     : Spectral frequencies (band averaged)
    g2       : Squared coherence spectrum
    g2crit   : Critical value that g2 must exceed to be significant at cl
    sxy_conj : Cross spectrum
    sxy      : Cross spectrum
    sxx      : Auto-spectrum
    syy      : Auto-spectrum
    pxy      : Phase spectrum
    pcl      : Confidence limits for the phase spectrum
    tl       : Time lags (t = pxy/(2*pi*freq))

    NOTES:
    ------
    The squared coherence is 1 at each frequency if no band averaging is
    performed.

    Positive phases indicate x leads y. Good to check in the time domain
    using cross_corr.

    DEPENDENCIES:
    -------------
    band_average

    """

    # Compute record length
    N = x.shape[0]

    # Compute fourier frequencies
    fj = np.fft.fftfreq(N,dt)

    # Compute power spectral density (Cooley-Tukey Method) for each time
    # series.
    xf = np.fft.fft(x)/N
    yf = np.fft.fft(y)/N

    sxy_conj = N*dt*xf*np.conjugate(yf)
    sxy      = N*dt*np.conjugate(xf)*yf
    sxx      = N*dt*xf*np.conjugate(xf)
    syy      = N*dt*yf*np.conjugate(yf)

    # Compute with only the positive side of the spectrum
    freq_ind = fj>=0
    fj       = fj[freq_ind]
    sxy_conj = sxy_conj[freq_ind]
    sxy      = sxy[freq_ind]
    sxx      = sxx[freq_ind]
    syy      = syy[freq_ind]

    # Compute the cross spectrum and quadrature spectrum
    #cxy   = 1.0*sxy.real
    #qxy = -1.0*sxy.imag

    # Compute the phase spectrum
    #pxy = np.arctan2(-1.0*qxy,cxy)
    #pxy = np.arctan(-1.0*qxy/cxy)

    # Band averaging
    fj_ba,sxy_conj_ba = band_averaging(sxy_conj,fj,stencil)
    _,sxy_ba = band_averaging(sxy,fj,stencil)
    _,sxx_ba = band_averaging(sxx,fj,stencil)
    _,syy_ba = band_averaging(syy,fj,stencil)
    #_,pxy_ba = band_averaging(pxy,fj,stencil)

    # Compute phase spectrum from band averaged quantities
    #pxy_ba = np.arctan(sxy_ba.imag/sxy_ba.real)
    pxy_ba = np.arctan2(sxy_ba.imag,sxy_ba.real)

    # Compute the cross spectral densities
    g2 = sxy_conj_ba * sxy_ba / sxx_ba / syy_ba

    # Confidence intervals on squared coherence
    alpha = 1.0 - cl
    g2crit = 1.0 - alpha**(2.0/(2.0*stencil - 2.0))
    #g2crit = (spi.stats.f.ppf(cl,2,2.0*stencil - 2.0)/
    #          ((2.0*stencil - 2.0)/2 +
    #           spi.stats.f.ppf(cl,2.0,2.0*stencil-2.0)))

    # Confidence intervals on phase
    pcl = np.arcsin((2.0/(2.0*stencil - 2.0) *
                    (1.0 - g2.real) / (g2.real) *
                    spi.stats.f.ppf(cl,2.0,2.0*stencil-2.0))**0.5)

    # Compute time lags
    tl = pxy_ba/(2.0 * np.pi * fj_ba)

    # End of function
    return {'freq':fj_ba,'g2':g2.real,'g2crit':g2crit,
            'sxy_conj':sxy_conj_ba,'sxy':sxy_ba,'sxx':sxx_ba,'syy':syy_ba,
            'pxy':pxy_ba,'pcl':pcl,'tl':tl}   


#===============================================================================
# Band averaging variance spectra
#===============================================================================
def band_averaging(psd,freq,stencil):
    """

    Function to band average power spectral density

    USAGE:
    ------
    freq_ba,psd_ba = band_averaging(psd,freq,stencil)

    PARAMETERS:
    -----------
    psd     : Power spectral density
    freq    : Frequencies of the PSD. It is assumed that freq[0] = 0. Thus the
              band averaging will start at freq[1].
    stencil : Number of frequencies to average.

    RETURNS:
    --------
    psd_ba  : PSD band averaged
    freq_ba : Frequency of the center points of the band averaged PSD

    NOTES:
    ------
    One sided spectra expected

    """

    # Make sure we are dealing with integers
    stencil = np.int(stencil)

    # Determine the size of the output arrays
    Nout = np.int(1 + np.ceil((psd.shape[0]-1)/stencil))

    # Preallocate variables
    if np.sum(np.iscomplex(psd)) > 0:
        psd_ba = np.zeros((Nout,))*1j
    else:
        psd_ba = np.zeros((Nout,))
    freq_ba = np.zeros((Nout,))

    # Allocate the first frequency and PSD estimate
    psd_ba[0] = psd[0]
    freq_ba[0] = freq[0]

    # Loop over the central bands to average
    for aa in range(1,Nout-1):

        # Find averaging limits
        ind_low = np.int((aa - 1.0)*stencil + 1.0)
        ind_high = np.int(ind_low + stencil)

        # Band averaging
        psd_ba[aa] = np.mean(psd[ind_low:ind_high])
        freq_ba[aa] = np.mean(freq[ind_low:ind_high])

    # Work on the last frequency bin
    psd_ba[Nout-1] = np.mean(psd[ind_high::])
    freq_ba[Nout-1] = np.mean(freq[ind_high::])

    # End of function
    return freq_ba,psd_ba



#===============================================================================
# Very simple implementation of a boxcar function for averaging purposes.
#===============================================================================
def boxcar(y,span,nanTreat=False,endTreat=True):
    """

    Usage:
    ------
    ylpf = boxcar(y,span)

    Input
    -----
       - y is the signal that will be filtered
       - span is the stencil width (should be an odd number, otherwise it will
         be forced to be so)
       - nanTreat (default = False) if true uses nansum instead of sum. Does
         not affect end treatment.

    Results
    -------
       - Returns filtered signal

    Notes
    -----
       - Boxcar running average cutoff frequency is found by
         fc = 0.6/(time_span)
       - Double running average filter:
         fc = 0.4429/(time_span)
       - Operates on first dimension of the array
       - A lot of assumptions about the data are made here, this function is by
         no means as robust as Matlab's smooth function. Only real valued
         numbers are assumed to be passed to the array and no repetition in the
         coordinate variable is assumed. Use at your own risk.

    """

    # Quick data check
    if span > y.shape[0]:
        print("Stencil of " + np.str(span) + " is larger than the " +
              "length of the array (" + np.str(y.shape[0]) + ")")
        sys.exit()


    # Span must be an odd number
    width = span - 1 + span % 2
    offset = np.int((width - 1.0)/2.0)

    # Preallocate variable
    ybox = np.zeros_like(y)

    # Find indices for averaging
    first = np.int(np.ceil(width/2.0) - 1.0)
    last = np.int(y.shape[0] - first - 1.0)

    if nanTreat:
        for aa in range(first,last+1):
            tmpW = np.sum(np.isfinite(y[aa-offset:aa+offset+1]),axis=0)
            ybox[aa] = np.nansum(y[aa-offset:aa+offset+1],axis=0)/tmpW
    else:
        for aa in range(first,last+1):        
            ybox[aa] = np.sum(y[aa-offset:aa+offset+1],axis=0)/width

    # Provide end treatment
    if endTreat:
        for aa in range(0,first):
            ybox[aa] = (np.sum(y[0:aa+offset+1],axis=0) /
                        (aa + offset + 1.0))

        for aa in range(last+1,y.shape[0]):
            ybox[aa] = (np.sum(y[aa-offset::],axis=0) /
                        (y.shape[0] - aa + offset + 0.0))
    
    else:
        ybox[:first] = y[:first]
        ybox[last:]  = y[last:]

    return ybox


# ==============================================================================
# Function to determine the frequency at which a signal that is not measured
# is aliased to.
# ==============================================================================
def aliased_frequency(dt_sampling,signal_period):
    """

    Determines the frequency into which a high frequency signal is aliased to.
    This is due to aliasing.

    PARAMETERS:
    -----------
    dt_sampling     : Sampling period [s,min,hr,day,etc]
    signal_period   : Period of high frequency signal that is aliased [s,etc]

    RETURNS:
    --------
    aliased_frequency : Frequency that will contain the information of the
                        high frequency process:
    """

    # Compute Nyquist frequency
    fN = 1.0/(2.0 * dt_sampling)

    # Frequency of the signal
    fk = 1.0/signal_period

    # Screen output
    print(" ")
    print("freq_signal = " + np.str(fk))
    print("freq_nyquist = " + np.str(fN))

    # Determine if the signal is indeed aliased or if it is resolved
    if fk <= fN:
        print("The signal is theoretically resolved given the sampling rate")
        return

    # Figure out the aliased frequency
    k = -1
    while True:
        k = k + 1
        aliased_frequency = fk - 2.0 * k * fN
        if aliased_frequency < fN:
            print('aliased_frequency = ' + np.str(aliased_frequency))
            print('aliased_period = ' + np.str(1.0/aliased_frequency))
            break
        elif k > 100000000:
            print('Maximum number of iterations reached, exiting')
            aliased_frequency = np.NaN
            break

    # Exit function
    return aliased_frequency



# =============================================================================
# Smooth1d Loeess filter
# =============================================================================
def smooth1d_loess(data,data_grid,span_x,est_grid=None):
    '''
    1-dimensional loess smoother. The smoothed value at each grid point is
    found from a weighted least-squares regression of the poins within
    +- SPAN_X of the grid point to a quadratic surface.

    PARAMETERS:
    -----------
    DATA      : Numpy vector of the data to be smoothed. FLag missing values
                as NaN
    DATA_GRID : Locations on the rid where data are located
    SPAN_X    : Filter half-power point (scalar). The larger the number, the
                more smoothing is performed. For the tricubic weighting
                function used here, the smoothing is approximately equivalent
                to using a running average of length equal to ~0.6*span_x.
                However, the spectral characteristics of this smoother are
                usually much more desirable.
    EST_GRID   : (optional) Are the locations of the grid where smoothed
                 estimates are desired. The estimate grid can be irregular
                 and non-monotonic. Any points in EST_GRID outside of the range
                 of DATA_GRID will not have SM_DATA-NaN

    RETURNS:
    --------
    SM_DATA    : Vector the same size as EST_GRID containing the smoothed
                 version of DATA at the locations specified by EST_GRID.
    FLAG       : Vector the same size as EST_GRID that is set to 1 when the
                 smoothed estimate is outside the range of the data within
                 SPAN_X of that grid point and 0 otherwise. When the smoothed
                 estimate is out of range, the estimate will be included in the
                 output SM_DATA. This will typically occur near the edges of
                 the DATA series or when SPAN_X only encompasses a small number
                 of grid points in DATA_GRID. While the smoothed estimate is
                 usually only marginally out-of-range in these cases, care
                 should be used when considering these points because the
                 smoothed estimate may not be very good at that particular
                 point. If there are many such points, consider using a larger
                 SPAN_X (smoothing over more points).

    NOTES:
    ------

    - Written in Matlab format by Larry O'Neill, May 27, 2006, based on a
      Fortran version written by Michael Schlax
    - LWO Update June 8, 2006 to include FLAG parameter
    - LWO Update June 9, 2006 to fix an error on some matlab versions with
      the specification of the weighting function and to include smoothed
      estimates in SM_DATA when FLAG==1
    - LWO Update April 17, 2007 to fix a check statement regarding
      clipping values of est_grid outside the range of data_grid. Fixes
      issues arising when trying to extrapolate outside range of DATA_GRID.
      Also made routine more efficient by rewriting some inefficient code and
      removing a few redundant computations. Now runs about 40% faster.
    - Written in python by Gabriel Garcia Medina, 19 April 2015. This code was
      provided to me in Matlab by Dr. Jonathan Nash (Oregon State University),
      as part of the OC683 course during the Spring 2015 term.

    '''

    if est_grid is None:
        est_grid = 1.0*data_grid

    # Data check (need to do some more)
    data_finite_flag = np.isfinite(data)

    # Preallocate output variables
    sm_data = np.zeros_like(est_grid)*np.NaN
    flag = np.zeros_like(sm_data)*np.NaN

    # Normalize grids by span_x
    data_grid = data_grid/span_x
    est_grid = est_grid/span_x

    # Loop through all the points in EST_GRID that are within the range of
    # DATA_GRID
    for ii in range(est_grid.shape[0]):

        # Verify that the data is within the range of the input data
        if est_grid[ii] < data_grid.min():
            continue
        elif est_grid[ii] > data_grid.max():
            return

        dx = data_grid - est_grid[ii];
        distx = np.abs(dx)
        igood = np.logical_and(distx < 1,data_finite_flag)
        ngood = np.sum(igood)


        # Need at least 3 data points for regression
        if ngood < 3:
            continue

        # Extract data to work with
        datareg = data[igood]
        dxsel = dx[igood]
        distsel = distx[igood]

        # Use tricubic weighting function for the filter weights.
        w = (1.0 - distsel**3)**3

        # Compute array of w*(1+x+x2)
        xin = np.zeros((w.shape[0],3))*np.NaN
        xin[:,0] = w
        xin[:,1] = w*dxsel
        xin[:,2] = xin[:,1]*dxsel

        # Least-squares solution to the over-determined set of equations.
        # It solves the equation wY = B*(wX), where Y is the smoothed estimate
        # of the input DATA, X is the coordinate relative to the grid point
        # of the smoothed estimate, and B are the regression coefficients.
        # Solves the QR decomposition

        B = np.linalg.lstsq(xin,w*datareg)

        # Smoothed value is just the first regression coefficient, since the
        # grid point was chosen such that it is at x=0. Also need to check
        # that the regression point is within the range of the data points
        # used in fitting the quadratic surface. It should be out of range
        # only rarely, and if it is, the smoothed estimate at that point is
        # given and the FLAG is set to 1.
        sm_data[ii] = B[0][0]

        if np.logical_and(B[0][0] > datareg.min(), B[0][0] < datareg.max()):
            flag[ii] = 1


    # Exit function
    return sm_data,flag



# ==============================================================================
# Slow Fourier Transform
# ==============================================================================
def slow_dft(yt,freq=None,dt=1):
    '''
    Slow Discrete Fourier transform code

    PARAMETERS:
    -----------
    yt          : Input time series (vector)
    freq        : Frequencies to compute the fourier transform (optional)
    dt          : Must be specified if dt is specified (yt sampling interval)

    RETURNS:
    --------
    yf          : Fourier coefficients (vector)

    NOTES:
    ------
    I have not extensively tested this function but I believe it works fine.
    It is really slow.

    '''

    # Compute length of the time series
    N = yt.shape[0]

    # Preallocate Fourier matrix (complex variable)
    if freq is None:
        yf = 1j * np.zeros_like(yt)
        N_fourier = N
    else:
        yf = 1j * np.zeros_like(freq)
        N_fourier = yf.shape[0]

    # Loop over fourier frequencies
    for aa in range(N_fourier):

        # Temporary sum variable (this will be a complex variable)
        tmpsum = np.complex(0)

        # Loop over input time series
        for bb in range(N):

            if freq is None:
                tmpsum += yt[bb] * np.exp(-2j * np.pi * aa * bb / N)
            else:
                tmpsum += yt[bb] * np.exp(-2j * np.pi * bb * freq[aa] * dt)

        # Allocate output for the given frequency
        yf[aa] = tmpsum / N

    # Exit function
    return yf



#===============================================================================
# Compute statistics from two time series
#===============================================================================
def basic_stats(x,y):
    '''
    Compute basic statistics for time series analysis.

    PARAMETERS:
    ----------
    x,y    : Time series of the same length to compare. x is thought to be the
             real (measured) data.

    RETURNS:
    --------
    A dictionary with the following parameters
    N      : Length of the time series
    rmse   : Root-mean squared error [input units squared]
    nrmse  : Normalized root-mean squared error [percentage]
    bias   : Bias [input units]
    pbias  : Percent bias
    si     : Scatter index
    r2     : Linear correlation coefficient    

    '''

    # Length of the time series
    N = np.size(x)

    # Root-mean-squared error
    rmse = (np.sum((x - y)**2)/N)**0.5

    # Normalized root mean squared error
    nrmse = 100.0*(np.sum(((x - y)/x)**2)/N)**0.5

    # Percent error
    pe = 100 / N * np.sum((y - x)/x)

    # Bias
    bias = np.sum(y - x)/N

    # Percent bias
    pbias = (np.sum(y) - np.sum(x)) / np.sum(x) * 100.0

    # Scatter index
    si = rmse/np.mean(x)

    # Linear correlation coefficient
    #r2 = np.corrcoef(x,y)
    r2,_ = cross_corr(x,y,0)
    
    # Effective length of time series using artificial skill method
    #Nstar = essize(x,y)[0]
    
    # Function output
    return {'N':N, 'rmse':rmse,'nrmse':nrmse,'pbias':pbias,
            'bias':bias,'si':si,'r2':r2[0],'pe':pe}



#===============================================================================
# zero crossing
#===============================================================================
def zero_crossing(x,d='up'):
    '''
    Find zero crossings in a signal

    PARAMETERS:
    -----------
    x  : Time series
    d  : (Optional) Upcrossing or downcrossing flag. Accepts 'up' (default)
         and 'down'.

    RETURNS:
    --------
    z  : Indices where upcrossing or downcrossings happen.

    '''

    # Upcrossing and downcrossing
    if d == 'down':
        y = np.logical_and(x[:-1]>0,x[1:]<=0)
    else:
        y = np.logical_and(x[:-1]<=0,x[1:]>0)
        if d != 'up':
            print('Could not understand input, upcrossing used')

    # Find indices where the time series changes sign
    z = np.where(y == True)[0]

    return z


#===============================================================================
# Empirical cumulative distribution function
#===============================================================================
def ecdf(x):
    '''
    Empirical cumulative distribution function

    USAGE:
    ------
    [xS,p,sInd] = ecdf(x)

    PARAMETERS:
    -----------
    x  : Time series of independent values

    RETURNS:
    --------
    xS   : Sorted x values
    p    : Probability
    sInd : Sorting index

    NOTES:
    ------
    - Adapted from code by Katy Serafin (OSU)

    '''

    # Force the array to be a numpy array
    x = np.array(x)

    # Remove nans in the data
    x = x[np.isfinite(x)]

    # Length of the finite runup values
    xLen = np.double(x.shape[0])

    # Compute probability
    p = np.arange(1.0,xLen+1.0,1.0)/xLen

    # Sort time series
    sInd = np.argsort(x)
    xS = x[sInd]
    

    # End of Function
    return xS,p,sInd


#==============================================================================
# Synthetic Time Series
#==============================================================================
def synthetic_ts(freq,spec,rseed=None):
    '''
    Generate a synthetic time series from an input spectrum using Fourier 
    techniques.

    USAGE:
    ------
    [t,syntTs] = synthetic_ts(freq,spec[,rseed])

    PARAMETERS:
    -----------
    freq   : Frequency matrix
    spec   : One-sided power spectral density function
    rseed  : (Optional) Seed for numpy's random number generator

    RETURNS:
    --------
    t      : Time vector which is a function of df and freq_max
    syntTs : Synthetic time series derived from the spectrum

    METHOD:
    ------
    1. Frequency axis is assumed to go from 0 to the Nyquist frequency
    2. Amplitude of the discrete Fourier transform (DFT)
       D = (df * S)**0.5 for freq[0] and freq[-1]
       D = (0.5 * df * S)**0.5 for freq[1:-1]
    3. Equally distributed random phases are generated from 0 and 2pi 
    4. Specify the real and imaginary parts of the spectrum
       Y = S e**(i*phases)
    5. Construct the two sided DFT (even number of entries will result)
    6. Perform an inverse DFT and generate the time series. 

    NOTES:
    -----
    1. Length of time series is determined by the frequency resolution
    2. Time resolution of the synthetic time series is determined by the
       maximum frequency. 
    3. frequency spacing must be constant from 0 to maximum frequency.

    '''

    # Find the frequency interval
    df = freq[2] - freq[1]

    # Find length of the time series.
    dt = 1.0/(2.0 * np.max(freq))
    N = 2 * (freq.shape[0] - 1)
    
    # Compute the amplitudes of the Discrete Fourier Transform
    ampDFT = np.zeros_like(spec)
    ampDFT[1:-1] = (0.5*df*spec[1:-1])**0.5
    ampDFT[0] = (df*spec[0])**0.5
    ampDFT[-1] = (df*spec[-1])**0.5

    # Obtain phases from a random number generator
    if rseed:
        print('Seeding with ' + np.str(rseed))
        np.random.seed(rseed)
        
    phases = np.random.rand(ampDFT.shape[0]) * 2.0 * np.pi

    # Incorporate the phases to the time Fourier transform
    ampDFTPhase = ampDFT * np.exp(1j*phases)

    # Find the complex conjugate for the negative part of the Fourier Transform
    conjAmpDFTPhase = np.conj(ampDFTPhase[1:-1])

    # Put together the final form of the Discrete Fourier Transform
    finalDFT = np.r_[ampDFTPhase,np.flipud(conjAmpDFTPhase)]

    # Compute the inverse Fourier Transform
    syntTs = np.fft.ifft(finalDFT*N).real
    t = np.arange(0,N*dt,dt)

    return t,syntTs



	
#==============================================================================
# Directional spread function
#==============================================================================
def dir_spread(thetad,S,ind_f):
    '''
    Calculates directional spread (as a function of frequency)

    USAGE:
    ------
    thetad, G_theta = wrapped_sf(thetad_peak,sigma_thetad[,mtheta,N,eq_bins])

    PARAMETERS:
    -----------
    thetad   : Spectral directions in degrees
    S        : Frequency-direction sprectrum, S(f,theta) [m**2/Hz/rad]
    ind_f    : Frequency indices (give peak frequency index for a single directional spread value)

    RETURNS:
    --------
    sigd    : Directional spread in degrees

    METHOD:
    -------
    1. tan[2*theta_mean(f)] = int{-pi,pi}[(sin(2theta) * S(f,theta)) * dtheta] / 
							  int{-pi,pi}[(cos(2theta) * S(f,theta)) * dtheta]
    2. sigma**2(f)= int{-pi,pi}[(sin**2[theta-theta_mean(f)] * S(f,theta)) * dtheta] / E(f)
	3. sigma(f) = (sigma**2(f)) **.5

    NOTES:
    -----

    '''
    # Convert to radians
    theta = np.deg2rad(thetad)
    
    # Calculate mean direction
    tan2thetam = np.trapz( np.squeeze(np.sin(2*theta)*S[ind_f,:]),x=theta) / np.trapz( np.squeeze(np.cos(2*theta)*S[ind_f,:]),x=theta)
    theta_m = np.arctan(tan2thetam * .5)
    
    # Calculate frequency spectrum
    ef = np.trapz(S,x=theta,axis=1)
    
    # Calculate the directional spread
    sig2 = np.trapz( np.squeeze(np.sin(theta*theta_m)**2 *S[ind_f,:]),x=theta) / ef[ind_f]
    sigd = np.rad2deg(sig2 ** .5) 
    
    return sigd



#===============================================================================
# Band pass filtering
#===============================================================================
def freq_dom_flt(y,dt,freqmin=None,freqmax=None,demean=True,window=True):
    """
    Frequency domain filtering
    
    PARAMETERS:
    -----------
    y       : Time series to filter
    dt      : Sampling interval of the time series of interest
    freqmin : (Optional) minimum frequency to keep
    freqmax : (Optional) maximum frequency to keep
    demean  : Remove mean from data before computing the Fourier Transform.
    window  : Apply a hamming window to the data
    
    RETURNS:
    --------
    yflt    : Frequency domain filtered dataset
    
    NOTES:
    ------
    - Need to provide freqmin or freqmax at least. Provide both for band pass
      filtering.
    - The mean will be added back to the time series after performing the 
      frequency domain filtering when demean=True.
    
    """
    
    # Check for frequency input
    if not freqmin and not freqmax:
        print('Need to provide one freqmin or freqmax at least')
        yflt = np.zeros_like(y) * np.NAN
        return yflt
    
    # Get fourier frequencies
    freq = np.abs(np.fft.fftfreq(y.shape[0],dt))
    
    # Find frequencies to remove
    if freqmin and freqmax:
        # Band pass filtering
        rmvInd = np.logical_or(freq<freqmin,freq>freqmax)
    elif freqmin:
        # High pass filtering
        rmvInd = freq < freqmin
    else:
        # Low pass filtering
        rmvInd = freq > freqmax
        
    # Remove mean if requested
    if demean:
        dataMean = np.mean(y)
        y = y - dataMean
    else:
        dataMean = 0.0
        
    # Window data if requested
    if window:
        hamWind = spi.signal.hamming(y.shape[0])
        y       *= hamWind
    else:
        hamWind = np.ones_like(y)
        
    # Apply frequency domain filter
    dft = np.fft.fft(y)
    dft[rmvInd] = 0.0j
    yflt = (np.fft.ifft(dft).real)/hamWind + dataMean
    
    return yflt


#===============================================================================
# Linear Regression
#===============================================================================
def linReg(x,y):
    """
    Compute linear regression using the least squares method
    
    PARAMETERS:
    -----------
    x : Coefficient matrix of size (N,M)
    y : Target Predictor of size (M,)
    
    RETURNS:
    --------
    B : Coefficents of linear regression
    p : Array of performance statistics
    
    Notes:
    ------
    y can be thought of as the measurements.
    
    Examples:
    ---------
    Remove a yearly signal from the data
    >>> y = data
    >>> months = np.arange(0,y.shape[0]) # Assuming data is sampled monthly
    >>> months[np.isnan(y)] = np.NAN
    >>> x = np.ones((2,y.shape[0]))
    >>> x[1,:] = np.sin(2.0*np.pi/12*months)
    >>> x[2,:] = np.cos(2.0*np.pi/12*months)
    >>> B,p = linReg(x,y)
    >>> yReg = np.array([B[aa]*x[aa,:] for aa in range(B.shape[0])])
    >>> yReg = np.sum(yReg,axis=0)
    >>> yClean = y - yReg
    
    """
    
    # Confidence level
    cl = 0.95
       
    # Linear regression --------------------------------------------------------
    Dmat = np.zeros((x.shape[0],x.shape[0]))
    Zmat = np.zeros((x.shape[0],))
    
    # Fill the coefficient matrix
    for aa in range(Dmat.shape[0]):
        for bb in range(Dmat.shape[1]):
            Dmat[aa,bb] = np.nanmean(x[aa,:]*x[bb,:])
        Zmat[aa] = np.nanmean(y*x[aa,:])
    B = np.linalg.lstsq(Dmat,Zmat)[0]

    # Evaluate the regression model
    rMod = np.array([B[aa]*x[aa,:] for aa in range(B.shape[0])])
    rMod = np.sum(rMod,axis=0)
    
    # Confidence interval for hindcast skill -----------------------------------
    
    # Compute the skill of the regression
    rho,stats = cross_corr(y,rMod,0)
    rhosq = rho[0]**2 # Skill
    
    # Compute the effective sample size
    Nstar = essize(y)[0]
    
    # Compute the critical skill
    M = x.shape[0] - 1
    sc = scrit(Nstar,M,cl=cl)
    
    # Confidence interval in each of the regression parameters
    try:
        Dinv = np.linalg.inv(Dmat)
        sigma_b1 = ((np.nanvar(y) * (1.0 - rhosq) * np.diagonal(Dinv))/
                    (Nstar - M - 1.0))
        db = (sigma_b1**0.5) * spi.stats.t.ppf(0.5+cl/2.0,Nstar-M-1)
    except:
        print('Singular regression matrix')
        print('  Cannot do one parameter models yet')
        db = 0.0
        
    # Prepare for ouptut -------------------------------------------------------
    p = {'skill':rhosq,'Nstar':Nstar,'scrit':sc,'db':db}
    
    return B,p
    
#===============================================================================
# Effective Sample Size
#===============================================================================
def essize(x,y=None):
    """

    This function computes the effective sample size for a given time series
    based on the 'Artificial Skill Method' and the 'Probability Density
    Function Method'
    
    PARAMETERS:
    -----------
    x : Time series, missing values should be flagged as NaNs
    y : Second time series (optional)
    
    RETURNS:
    --------
    essize : Effective sample size

    NOTES:
    ------
    essize[0]: Using artificial skill method
    essize[1]: Using the PDF method
    
    """
    
    # Check for optional parameters         
    if y is None:
        y = np.copy(x)
    
    # Estimate the effective sample size and confidence interval with the
    # artificial skill method --------------------------------------------------
    
    # Get number of samples
    Nall = np.sum(np.isfinite(x*y))
    
    # Compute the bounds for the lagged auto-correlation
    klow  = np.floor(0.6*Nall)
    khigh = np.ceil(0.8*Nall)
    
    # Cross-correlation to the appropriate bounds
    rho,stats = cross_corr(x,y,khigh)
    
    # Eliminate short lags
    keep_lags = np.abs(stats[:,0])>=klow
    rho = rho[keep_lags]
    stats = stats[keep_lags,:]
    rhosq = rho**2
    
    # Compute the 'Artificial Skill'
    neg_lag_ind = stats[:,0] < 0
    neg_lags = rhosq[neg_lag_ind] * stats[neg_lag_ind,-1]
    pos_lag_ind = stats[:,0] > 0
    pos_lags = rhosq[pos_lag_ind] * stats[pos_lag_ind,-1]
    pos_lags = pos_lags[::-1]   
    
    Ask = 1.0/(2*(khigh-klow+1)) * np.sum(pos_lags + neg_lags)
     
    # Compute nuhat
    nuhat_ask = 1.0/Ask;
     
    # Compute effective sample size
    nstar_ask = np.floor(nuhat_ask*Nall)
     
     
    # Compute N* with the 'PDF' method -----------------------------------------
    neg_lags_rho = rhosq[neg_lag_ind]
    pos_lags_rho = rhosq[pos_lag_ind][::-1]
    neg_lags_n   = stats[neg_lag_ind,-1]
    pos_lags_n   = stats[pos_lag_ind,-1][::-1]
    A1 = (1.0/(2*(khigh-klow+1)) * 
          np.sum(pos_lags_rho/(1.0-pos_lags_rho) + 
                 neg_lags_rho/(1.0-neg_lags_rho)))
    A2 = (1.0/(2*(khigh-klow+1)) * 
          np.sum(pos_lags_n*pos_lags_rho/(1.0-pos_lags_rho) +
                 neg_lags_n*neg_lags_rho/(1.0-neg_lags_rho)))
    
    # Compute nuhat
    nuhat_pdf = (1.0 + 4.0*A1)/A2
    
    # Compute effective sample size
    nstar_pdf = np.floor(nuhat_pdf*Nall)

    # Create output matrix
    return [np.int64(nstar_ask),np.int64(nstar_pdf)]


#===============================================================================
# Critical skill for significance
#===============================================================================
def scrit(N,dof,cl=0.95):
    """
    Critical skill at a given significance level
    
    PARAMETERS:
    -----------
    N    : Number of effective sample size (see signal.essize)
    dof  : Degrees of freedom (i.e. number of variables fit, do not consider
           the constant term)
    cl   : Confidence level (e.g. for 95% is 0.95)
    
    RETURNS:
    --------
    scrit : Critical skill for the model to be significant at a given 
            confidence level

    """
    
    # Inverse of Fisher's cumulative distribution function
    finv = spi.stats.f.ppf(cl,dof,N-dof-1)
    
    # Critical model skill
    scrit = (dof*finv)/((N-dof-1) + dof*finv)

    return scrit

#===============================================================================
# Running variance
#===============================================================================
def runVar(y,span,nanTreat=False):
    """

    Compute running variance
    
    Usage:
    ------
    var = runVariance(y,span)

    Input
    -----
       - y is the signal whose variance will be computed
       - span is the stencil width (should be an odd number, otherwise it will
         be forced to be so)
       - nanTreat (default = False) if true uses nansum instead of sum. Does
         not affect end treatment.

    Results
    -------
       - Returns variance along the signal

    Notes
    -----
       - Operates on first dimension of the array
       - A lot of assumptions about the data are made here, this function is by
         no means as robust as Matlab's smooth function. Only real valued
         numbers are assumed to be passed to the array and no repetition in the
         coordinate variable is assumed. Use at your own risk.

    """

    # Quick data check
    if span > y.shape[0]:
        print("Stencil of " + np.str(span) + " is larger than the " +
              "length of the array (" + np.str(y.shape[0]) + ")")
        return


    # Span must be an odd number
    width = span - 1 + span % 2
    offset = np.int((width - 1.0)/2.0)

    # Preallocate variable
    ybox = np.zeros_like(y) * np.NAN

    # Find indices for averaging
    first = np.int(np.ceil(width/2.0) - 1.0)
    last = np.int(y.shape[0] - first - 1.0)

    if nanTreat:
        for aa in range(first,last+1):
            ybox[aa] = np.nanvar(y[aa-offset:aa+offset+1],axis=0)
    else:
        for aa in range(first,last+1):
            if np.isnan(np.sum(y[aa-offset:aa+offset+1])):
                continue
            ybox[aa] = np.var(y[aa-offset:aa+offset+1],axis=0)

    return ybox

