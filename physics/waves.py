"""

Series of functions to manage wave information.

Author:
-------
Nearshore Modeling Group
Gabriel Garcia Medina
ggarcia@coas.oregonstate.edu

Log of edits:
-------------
April 2014 - Created module
  Gabriel Garcia Medina (ggarcia@coas.oregonstate.edu)

External dependencies:
  numpy
  scipy

Internal dependencies:
  gsignal

"""

from __future__ import division,print_function

# Import generic modules
import numpy as np
import scipy
import scipy.optimize
import scipy.signal

# Internal modules
import pynmd.data.signal as gsignal
import pynmd.data.angles as _gangles

#===============================================================================
# Compute radian frequency from wave number
#===============================================================================
def idispersion(k,h):
    '''
    Inverse linear dispersion relation

    USAGE:
    ------
    sigma = idispersion(k,h)

    INPUT:
    ------
    k      : Wave number [m**-1]
    h      : water depth [m]

    OUTPUT:
    -------
    sigma  : radian frequency [Hz]

    '''

    # Get radian wave frequency
    sigma = (9.81*k*np.tanh(k*h))**0.5

    # Return stokes velocities
    return sigma


#===============================================================================
# Compute wave number
#===============================================================================
def dispersion(period,h,u=None):
    '''

    Computes the linear dispersion relation.

    USAGE:
    ------
    k = dispersion(period,h,u)

    Input:
    ------
    Period   : Wave period [s]
    h        : Water depth [m]
    u        : (Optional) Ambient current velocity [m/s]

    RETURNS:
    --------
    k        : Wave number (2*pi/wave_length) [m**-1]

    '''

    if u is None:
        u = 0.0

    # Initialize wave number
#     k = np.zeros(len(period))
    sigma = (2.0*np.pi)/period
#     kinit = (((sigma**2)*h/9.81)*
#              ((1/np.tanh(((sigma**2)/9.81*h)**0.75))**(2.0/3.0)))
    kinit = 0.0001

    def f(k,sigma,h,u):
        y = (9.81*k*np.tanh(k*h))**0.5 - sigma + u*k
        return y

    if h <= 0:
        k = np.NAN
    else:
        try:
            k = scipy.optimize.newton(f,kinit,fprime=None,
                                      args=(sigma,h,u),maxiter=100)
        except RuntimeError:
            k = np.NAN

    return k


#===============================================================================
# Compute wave number using the Kirby and Dalrymple (1986) composite equation
#===============================================================================
def dispersion_kd86(period,h,wh,u=None):
    '''

    Return the wave number from the composite dispersion relation by
    Kirby and Dalrymple (1986).

    USAGE:
    ------
    k = dispersion_kd86(period,h,wh,u)

    Input:
    ------
    Period   : Wave period [s]
    h        : Water depth [m]
    wh       : Wave height [m]
    u        : (Optional) Ambient current velocity [m/s]

    RETURNS:
    --------
    k        : Wave number (2*pi/wave_length) [m**-1]

    CELERITY EQUATION:
    ------------------
    (c-u)^2 = g / k * (1 + f1 * eps^2 * D) tanh( k * h + f2 * eps)
    f1      = tanh(k*h)^5
    f2      = (k*h/sinh(k*h))^4
    eps     = k*wh/2
    D       = (8 + cosh(4*k*h) - 2 * tanh(k*h)^2)/(8*sinh(k*h)^4)

    REFERENCES:
    -----------
    Catalan, P., and M. C. Haller, 2008: Remote sensing of breaking wave phase
        speeds with application to non-linear depth inversions. Coastal
        Engineering, 55, 93 - 111.
    Kirby, J. T., and R. A. Dalrymple, 1986: An approximate model for nonlinear
        dispersion in monochromatic wave propagation models. Coastal
        Engineering, 9, 545-561.
    '''

    # Depth averaged ambient velocity
    if u is None:
        u = 0

    # Compute radian frequency
    sigma = (2.0*np.pi)/period

    # Initialize wave number
    #kinit = 0.0001                # Manually initialize
    kinit = dispersion(period,h,u) # Initialize with Airy dispersion relation

    # Kirby and Dalrymple (1986) dispersion function
    def f(k,sigma,h,wh,u):
        eps = k*wh/2
        f1  = np.tanh(k*h)**5
        f2  = (k*h/np.sinh(k*h))**4
        D   = ((8.0 + np.cosh(4.0*k*h) - 2.0*(np.tanh(k*h)**2))/
               (8.0*(np.sinh(k*h)**4)))
        y   = ((9.81*k*np.tanh(k*h+f2*eps)*(1.0+f1*(eps**2)*D))**0.5 +
               u*k - sigma)
        return y

    if h <= 0:
        k = np.NAN
    else:
        try:
            k = scipy.optimize.newton(f,kinit,fprime=None,
                                      args=(sigma,h,wh,u),maxiter=100)
        except RuntimeError:
            k = np.NAN

    return k


#===============================================================================
# Compute wave number using the Booij
#===============================================================================
def dispersion_booij(period,h,wh,u=None):
    '''

    Return the wave number from the composite dispersion relation by
    Booij (1981).

    USAGE:
    ------
    k = dispersion_booij(period,h,wh,u)

    Input:
    ------
    Period   : Wave period [s]
    h        : Water depth [m]
    wh       : Wave height [m]
    u        : (Optional) Ambient current velocity [m/s]

    RETURNS:
    --------
    k        : Wave number (2*pi/wave_length) [m**-1]

    CELERITY EQUATION:
    ------------------
    (c-u)^2 = g / k * tanh(k*(h+wh/2))

    REFERENCES:
    -----------
    Booij, N. 1981: Gravity waves on water with non-uniform depth and current.
        Tech. Rep. No. 81-1. Dept. Civil Engineering, Delft University of
        Technology.
    Catalan, P., and M. C. Haller, 2008: Remote sensing of breaking wave phase
        speeds with application to non-linear depth inversions. Coastal
        Engineering, 55, 93 - 111.
    '''

    # Depth averaged ambient velocity
    if u is None:
        u = 0.0

    # Compute radian frequency
    sigma = (2.0*np.pi)/period

    # Initialize wave number
    #kinit = 0.0001                     # Manually initialize
    kinit = dispersion(period,h+wh/2,u) # Initialize with Airy dispersion relation

    # Booij (1981) dispersion function
    def f(k,sigma,h,wh,u):
        y = (9.81*k*np.tanh(k*(h+wh/2.0)))**0.5 + u*k - sigma
        return y

    if h <= 0:
        k = np.NAN
    else:
        try:
            k = scipy.optimize.newton(f,kinit,fprime=None,
                                      args=(sigma,h,wh,u),maxiter=100)
        except RuntimeError:
            k = np.NAN

    return k


#===============================================================================
# Compute dispersion relation from Nwogu's equations
#===============================================================================
def dispersion_nwogu(period,h,alpha):
    '''
    Compute dispersion relation for the Boussinesq equations based on
    Nwogu 1993.

    PARAMETERS:
    -----------
    period    : Wave period [s]
    h         : Water depth [m]
    alpha     : Non-dimensional wave steepness parameter (see notes)

    RETURNS:
    --------
    k         : Wave number [m**-1]

    NOTES:
    ------
    alpha = 0       for celerity at the surface
    alpha = -1/3    solves the traditional depth averaged equations.
    alpha = -0.39   gives similar results to linear wave theory (Nwogu 1993)

    REFERENCES:
    -----------
    Nwogu, O., 1993: Alternative Form of Boussinesq Equations for Nearshore
      Wave Propagation. Journal of Waterway, Port, Coastal, and Ocean
      Engineering, 119, 618-638.

    '''

    # Compute randian frequency
    sigma = (2.0*np.pi)/period

    # Initialize wave number
    #kinit = dispersion(period,h) # Linear wave theory
    kinit = sigma**2 / 9.81       # Deep water equivalent
    #kinit = 0.05                 # Hard code parameter

    # Define the dispersion relation function
    def bd(k,sigma,h,alpha):
        #y = 9.81 * (k**2) * h * (1.0 - 1.0/3.0 * (k*h)**2) - sigma**2
        y = (9.81 * (k**2) * h * (1.0 - (alpha + 1.0/3.0) * ((k*h)**2)) /
            (1.0 - alpha * ((k*h)**2)) - sigma**2)
        return y

    # Use Newton-Rhapson method to find the wave number
    if h <= 0:
        k = np.NAN
    else:
        try:
            k = scipy.optimize.newton(bd,kinit,fprime=None,
                                      args=(sigma,h,alpha),maxiter=1000)
        except RuntimeError:
            k = np.NAN

    return k


#===============================================================================
# Linear wave approximations
#===============================================================================
def shallow_water_depth(period):
    '''

    h_shallow = shallow_water_depth(period)

    Find where the waves with the given period enter shallow water according to
    linear wave theory

    PARAMETERS:
    -----------
    period     : wave period [s]

    RETURNS:
    --------
    h_shallow  : Water depth where shallow water approximation is valid [m]

    '''

    # Initialize shallow water depth
    hinit = 2.0

    # Define shallow water limit
    def swd(h,period):
        y = dispersion(period,h) * h - np.pi/10.0
        return y

    # Find zeros of the swd function
    try:
        h = scipy.optimize.newton(swd,hinit,fprime=None,
                                  args=(period,),maxiter=100)
    except RuntimeError:
        h = np.NAN

    return h

#===============================================================================
# Wave length
#===============================================================================
def wave_length(period,h,verbose=True):
    '''
    Compute wave length using linear wave theory

    Parameters
    ----------
    period   : wave period [s]
    h        : water depth [m]

    Results
    -------
    wl_int   : real wave length [m]

    Screen output
    -------------
    wl_deep  : deep water wave length [m]
    wl_sha   : shallow water wave length [m]

    '''

    wl_deep = 9.81 * period**2 / 2.0 / np.pi
    wl_sha = period * np.sqrt(9.81 * h)
    k = dispersion(period,h)
    wl_int = 9.81 / 2.0 / np.pi * period**2 * np.tanh(k*h)

    if verbose:
        print(' ')
        print('---------------------------------------------------------')
        print('Wave Length deep water approx      = ' + np.str(wl_deep) + ' m')
        print('Wave Length shallow water approx   = ' + np.str(wl_sha) + ' m')
        print('Wave Length linear wave theory     = ' + np.str(wl_int) + ' m')
        print('---------------------------------------------------------')
        print(' ')

    return wl_int

#===============================================================================
# Compute wave number using the Kirby and Dalrymple (1986) composite equation
#===============================================================================
def celerity(period,h,u=None):
    '''

    Find celerity and group velocity from linear wave theory

    USAGE:
    ------
    c,n,cg = celerity(period,h,u)

    Input:
    ------
    period   : Wave period [s]
    h        : Water depth [m]
    u        : (Optional) Ambient current velocity [m/s]

    RETURNS:
    --------
    c        : celerity [m/s]
    n        : cg = c*n
    cg       : group velocity [m/s]

    '''

    # Radian frequency
    sigma = 2.0 * np.pi / period

    # Find wave number
    k = dispersion(period,h,u)

    # Celerity
    c = sigma/k

    # Group velocity
    n = 0.5 + k * h / np.sinh(2*k*h)
    cg = c * n

    return c,n,cg



#===============================================================================
# Compute stokes velocities from linear wave theory
#===============================================================================
def uvstokes(Hwave,Dwave,Lwave,WDepth,VertDisc):
    '''
    Compute stokes velocities

    Parameters
    ----------
    Hwave          : Wave height [m]
    Dwave          : Wave direction
    Lwave          : Wave length [m]
    WDepth         : Water depth [m]
    VertDisc       : Number of vertical layers

    Returns
    -------
    ustokes, vstokes

    '''

    k = 2*np.pi/Lwave                       # Compute wave number
    deno = np.sinh(2.0*k*WDepth)            # sinh(2kh)
    sigma = idispersion(k,WDepth)           # Radian frequency

    # Vertical discretization
    z = np.linspace(0,WDepth,VertDisc)
    nume = np.cosh(2.0*k*(WDepth - z))

    # Get zonal and meridional components
    kx = k*np.cos((270-Dwave)*np.pi/180)
    ky = k*np.sin((270-Dwave)*np.pi/180)

    # Stokes drift
    ustokes = (Hwave**2)/4.0*1027.0*9.81/sigma*kx*nume/deno
    vstokes = (Hwave**2)/4.0*1027.0*9.81/sigma*ky*nume/deno

    # Return stokes velocities
    return ustokes, vstokes


#===============================================================================
# TMA Spectrum
#===============================================================================
def tma(freq_peak,gamma,h,Hmo,freq_min=0.01,freq_max=1.0,freq_int=0.001,
        zeroth=True):
    '''
    Function to generate TMA spectrum

    Parameters
    ----------
    peak_freq    : Peak frequency [Hz]
    gamma        : Peak enhancement factor (3.3 is a good guess)
    h            : water depth [m]
    Hmo          : Significant wave height [m]

    Optional Parameters
    -------------------
    freq_min     : Minimum frequency to compute the spectrum [Hz]
    freq_max     : Maximum frequency to compute the spectrum [Hz]
    freq_int     : Frequency inteval [Hz]
    zeroth       : Prepend zeroth frequency (Default = True)

    Default values are 0.01, 1.0, and 0.001 Hz, respectively.

    Returns
    -------
    s_tma        : Tma spectrum as a function of frequency.
    freq         : Frequency axis [Hz]

    References
    ----------
    Bouws, E., H. Gunther, W. Rosenthal, and C. L. Vincent, 1985: Similarity
        of the wind wave spectrum in finite depth water 1. Spectral form.
        Journal of Geophysical Research, 90 (C1), 975-986.
    Hughes, S. 1985: The TMA shallow-water spectrum, description and
        applications. Technical Report CERC-84-7.

    Notes
    -----
    - The units of the spectrum still do not make sense to me, must be verified.
    - No scaling applied, alpha = 1
    - Zeroth frequency is added to the spectrum

    '''

    # For testing only
    #freq_peak = 0.1
    #gamma = 3.3
    #freq_min = 0.01
    #freq_max = 1.0
    #freq_int = 0.001
    #h = 10.0
    #Hmo = 1.0

    # Compute frequency vector
    freq = np.arange(freq_min,freq_max+freq_int,freq_int)

    # Constants for peak enhancement factor
    delta = np.ones_like(freq)*0.07
    delta[freq>freq_peak] = 0.09

    # Compute alpha parameter (equation 26,27) TMA report
    #wlen = 2.0*np.pi/dispersion(freq_peak**-1,h)
    #alpha = (2.0 * np.pi* Hmo / wlen)**2
    alpha = 1.0

    # TMA scaling factor
    omh = 2.0*np.pi*freq*(h/9.81)**0.5
    phi = 1.0 - 0.5 * (2.0 - omh)**2
    phi[omh<1.0] = (0.5*omh**2)[omh<1.0]
    phi[omh>2.0] = 1.0

    # Generate spectrum
    s_tma = (alpha * 9.81**2 * (2.0*np.pi)**-4 * (freq**-5) * phi *
             np.exp(-1.25 * (freq_peak/freq)**4) *
             gamma ** np.exp(-1.0*((freq - freq_peak)**2)
                             /(2.0*(delta**2)*freq_peak**2)))

    # Add zeroth frequency
    if zeroth:
        s_tma = np.r_[np.array([0.0]),spec_tma]
        freq = np.r_[np.array([0.0]),freq]

    # End of function
    return s_tma,freq


#===============================================================================
# JONSWAP Spectrum
#===============================================================================
def jonswap(freq_peak,Hmo,gamma=3.3,freq_min=0.01,freq_max=1.0,
            freq_int=0.001,goda=False,zeroth=True):
    '''
    Function to generate JONSWAP spectrum

    Parameters
    ----------
    peak_freq    : Peak frequency [Hz]
    Hmo          : Significant wave height [m]

    Optional Parameters
    -------------------
    gamma        : Peak enhancement factor (defaults to 3.3)
    freq_min     : Minimum frequency to compute the spectrum [Hz]
    freq_max     : Maximum frequency to compute the spectrum [Hz]
    freq_int     : Frequency inteval [Hz]
    goda         : Use Y. Goda's approximation to the Jonswap spectrum.
                   If false we force the scaling of the spectrum to give the
                   same wave height passed as argument.
    zeroth       : Prepend zeroth frequency (Default = True)

    Returns
    -------
    spec_jonswap : JONSWAP spectrum as a function of frequency.
    freq         : Frequency vector [Hz]

    References
    ----------
    Many, but a good description can be found in:
        Y. Goda, Random Seas and Design of Maritime Structures.

    Notes
    -----
    - Default frequency values are 0.01, 1.0, and 0.001 Hz.

    '''

    # For code development only -----------------------------------------------
    #freq_peak = 1.0/10.0
    #gamma = 3.3
    #freq_min = 0.01
    #freq_max = 0.5
    #freq_int = 0.001
    #Hmo = 1.0
    # -------------------------------------------------------------------------

    # Create frequency vector
    freq = np.arange(freq_min,freq_max+freq_int,freq_int)

    # Constants for peak enhancement factor
    sigma = np.ones_like(freq)*0.07
    sigma[freq>freq_peak] = 0.09

    # Goda's formulation
    if goda:
        # Beta parameter
        beta = (0.0624/(0.230 + 0.0336*gamma - 0.185*((1.9 + gamma)**-1)) *
                (1.094 - 0.01915*np.log(gamma)))

        # Generate spectrum
        spec_jonswap = (beta * (Hmo**2) * (freq_peak**4) * (freq**-5) *
                        np.exp(-1.25 * ((freq/freq_peak)**-4)) *
                        gamma ** (np.exp(-1.0 * (freq/freq_peak - 1.0)**2 /
                                         (2.0 * sigma**2))))

    else:

        # Generate spectrum
        spec_jonswap = ((freq**-5) *
                        np.exp(-1.25 * ((freq/freq_peak)**-4)) *
                        gamma ** (np.exp(-1.0 * (freq/freq_peak - 1.0)**2 /
                                         (2.0 * sigma**2))))

        # Scale parameter to match wave height energy in deep water
        # I am not sure this is the right way to proceed but I'll still do it.
        alpha = 1.0/16.0 * Hmo**2 / np.trapz(spec_jonswap,freq)
        spec_jonswap = spec_jonswap * alpha

    
    # Add zeroth frequency and expand the spectrum
    if zeroth:        
        spec_jonswap_old = np.r_[np.array([0.0]),spec_jonswap]
        freq_old = np.r_[np.array([0.0]),freq]
        
        freq = np.arange(0.0,freq_max+freq_int,freq_int)
        spec_jonswap = np.interp(freq,freq_old,spec_jonswap_old)

    # End of function
    return spec_jonswap,freq



#===============================================================================
# Add directional distribution to spectrum
#===============================================================================
def directional_spreading(spec,peak_dir,m,dirs=None):
    """
    This function computes a directionally spread spectrum from a frequency
    spectrum passed to the function.

    PARAMETERS:
    -----------
    spec         : frequency spectrum (as vector)
    peak_dir     : Peak wave direction in Nautical convention [degrees]
    m            : Directional width [cos(theta)**2m]
    dirs         : (Optional) vector of directions. If not given, the spectrum
                   will be computed every 5 degrees.

    RETURNS:
    --------
    dir_spec     : Directional spectrum in the same units given by the input
                   spectrum by degrees.
    dirs         : Vector with directions

    NOTES:
    ------
    Nautical convention refers to the direction the waves (wind) are coming
      (is blowing) from with respect to the true north measured clockwise.
    To recover the significant wave height
       4.004 * np.trapz(np.sum(dir_spec,axis=-1)*(dirs[2]-dirs[1]),freq)**0.5
    """

    # If direction vector is not passed as argument the directional distribution
    # will be computed every five degrees.
    try:
        if dirs == None:
            dirs = np.arange(0,360,5)
    except:
        pass

    # Change directions to radians for internal computations
    peak_dir_r = np.pi / 180.0 * peak_dir
    dirs_r = np.pi / 180.0 * dirs

    # Compute directional spread
    g = np.cos(0.5*(dirs_r - peak_dir_r))**(2*m)

    # Normalize the directional spread so that no energy is added
    # Scaling is done in direction space to account for the units
    # In other workds dir_spec will be of [m2/Hz-deg]
    dth = dirs[2] - dirs[1]
    g /= (np.sum(g)*dth)

    # Generate the directionally spread spectrum
    dir_spec = np.zeros((spec.shape[0],dirs.shape[0]))
    for aa in range(spec.shape[0]):
        dir_spec[aa,:] = spec[aa] * g

    # Return directional spectrum
    return dir_spec,dirs


#===============================================================================
# Compute bulk wave parameters from water surface elevation time series
#===============================================================================
def eta_bulk_params(eta,ot,band_ave=False,window=False,IGBand=[0.005,0.05]):
    """
    Compute bulk wave parameters from water surface elevation time series.

    Parameters:
    -----------
    eta        : Water surface elevation time series at a point [m]
    ot         : Time vector [s]
    band_ave   : (Optional) Bin average stencil
    Window     : (Optional) Application of a hanning window (True or False)
    IGBand     : (Optional) Infragravity wave frequency band 
                  defaults to 0.005-0.05 Hz   
    Output:
    -------
    Dictionary containing
    freq       : Spectral frequencies [Hz]
    spec       : Wave variance spectrum [m**2/Hz]
    cl         : 95% confidence levels on the spectral estimates
    Hs         : Significant wave height [m]
    H1         : Mean wave height [m]
    Tp         : Peak wave period [s]
    Tp_fit     : Peak wave period computed from second order polynomial fit
                 near Tp[s]
    Tm01       : First moment wave period [s]
    Tm02       : Second moment wave period [s]
    Te         : Energy period [s]
    Sw         : Spectral width (m0*m2/m1/m1 - 1)**2
    H1_t       : Mean wave height computed in the time domain [m]
    T1_t       : Mean wave period computed in the time domain [s]
    Tm01IG     : First moment wave period over the infragravity frequency band
    TpIG       : Peak wave period over the infragravity frequency band [s]
    skew_t     : Wave skewness from time domain
   
    Notes:
    ------
    mn are the different spectral moments

    """

    # Remove mean from the data
    etaw = eta - eta.mean()

    # Compute variance of original time series
    var_ts = np.var(etaw)
    
    # Compute in time domain
    [wh,wp,zcind] = whwpts(ot,eta)
    
    # Wave skewness
    skew_t = np.mean(etaw**3)/(np.mean(etaw**2)**1.5)

    # Data windowing
    if window:
        tmpwindow = scipy.signal.hanning(etaw.shape[0])
        etaw *= tmpwindow

    # Compute variance spectrum
    freq,spec = gsignal.psdraw(etaw,np.mean(ot[1:] - ot[:-1]))

    # If data has been windowed we must boost the variance of the spectrum
    # to match the original time series
    if window:

        # Compute variance of from the psd
        var_psd = np.sum(spec)*(freq[2]-freq[1])

        # Adjust variance from the windowed time series
        spec *= var_ts/var_psd

    # Band averaging if requested
    if band_ave:
        [freq,spec] = gsignal.band_averaging(spec,freq,band_ave)
    else:
        band_ave = 1

    # Compute confidence levels on the spectral parameters, to estimate noise
    # threshold.
    conf_lev = 0.95
    alpha = 1.0 - conf_lev
    cl_upper = scipy.stats.chi2.ppf(alpha/2,band_ave*2)
    cl_lower = scipy.stats.chi2.ppf(1.0-alpha/2.0,band_ave*2)
    cl = np.array([spec - spec * band_ave * 2.0 / cl_lower,
                   spec * band_ave * 2.0 / cl_upper - spec])

    # Compute bulk wave parameters
    bwp = fspec_bulk_params(freq,spec,IGBand)
    bwp['freq'] = freq
    bwp['spec'] = spec
    bwp['cl'] = cl

    bwp['H1_t'] = np.nanmean(wh)
    bwp['T1_t'] = np.nanmean(wp)
    bwp['skew_t'] = skew_t
    
    return bwp

#===============================================================================
# Function to compute bulk wave parameters from frequency spectrum
#===============================================================================
def fspec_bulk_params(freq,spec,IGBand=[0.005,0.05],zeroth=True):
    """
    Function to compute bulk wave parameters from frequency spectrum

    Parameters:
    -----------
    freq    : Vector of spectral frequencies [Hz]
    spec    : Frequency spectrum [m2/Hz]
    IGBand  : Infragravity wave frequency cutoff 
              defaults to 0.005-0.05 Hz
    zeroth  : If true will remove the first frequency from the analysis

    Returns:
    --------
    Dictionary containing bulk wave parameters
    Hs         : Significant wave height [m]
    H1         : Mean wave height [m]
    Tp         : Peak wave period [s]
    Tp_fit     : Peak wave period computed from second order polynomial fit
                 near Tp[s]
    Tm01       : First moment wave period [s]
    Tm02       : Second moment wave period [s]
    Te         : Energy period [s]
    Sw         : Spectral width (m0*m2/m1/m1 - 1)**2
    Tm01IG     : First moment wave period over the infragravity frequency band
    TpIG       : Peak wave period over the infragravity frequency band [s]
    HsIG       : Significant wave height in the infragravity frequency band [m]

    Notes:
    ------
    - mn are the different spectral moments

    """

    # Remove zeroth frequencies ------------------------------------------------
    if zeroth:
        spec = spec[1:]
        freq = freq[1:]

    # Compute spectral moments
    moment0 = np.trapz(spec,freq,axis=-1)
    moment1 = np.trapz(spec*freq,freq,axis=-1)
    moment2 = np.trapz(spec*(freq)**2,freq,axis=-1)
    momentn1 = np.trapz(spec*(freq)**-1,freq,axis=-1)    

    # Wave heights
    Hs = 4.004 * (moment0)**0.5
    H1 = Hs*2.0/3.0

    # Spectral width
    Sw = (moment0 * moment2 / moment1 / moment1 - 1)**0.5

    # Spectral periods -----------------------
    # Energy period
    Te = momentn1/moment0

    # Mean wave period
    Tm01 = moment0 / moment1

    # Second moment period
    Tm02 = (moment0 / moment2)**0.5

    # Peak wave period
    freq_max_ind = np.argmax(spec)
    Tp = freq[freq_max_ind]**-1

    # Peak wave period using a quadratic fit over the largest frequencies
    if freq_max_ind == 0:
        Tp_fit = np.nan
    elif freq_max_ind == freq.shape[0]-1:
        Tp_fit = np.nan
    else:
        minfreq = freq_max_ind - 1
        maxfreq = freq_max_ind + 2
        tmp_fit = np.polyfit(freq[minfreq:maxfreq],spec[minfreq:maxfreq],2)
        Tp_fit = (-1.0 * tmp_fit[1] / (2.0* tmp_fit[0]))**-1

    # Infragravity wave periods ------------------------------------------------
    igInd = np.logical_and(freq>=IGBand[0],freq<IGBand[1])
    moment0 = np.trapz(spec[igInd],freq[igInd],axis=-1)
    moment1 = np.trapz(spec[igInd]*freq[igInd],freq[igInd],axis=-1)
    Tm01IG  = moment0 / moment1
    
    HsIG = 4.004 * (moment0)**0.5     
    
    freq_max_ind = np.argmax(spec[igInd])
    TpIG = freq[igInd][freq_max_ind]**-1
    
    # Exit function
    return {'Hs':Hs,'H1':H1,'Tp':Tp,'Tp_fit':Tp_fit,'Tm01':Tm01,'Tm02':Tm02,
            'Te':Te,'Sw':Sw,'Tm01IG':Tm01IG,'TpIG':TpIG,'HsIG':HsIG}

#===============================================================================
# Bulk parameters from a directional spectrum
#===============================================================================
def spec_bulk_params(freq,dirs,spec,IGBand=[0.005,0.05],zeroth=True):
    """
    Function to compute bulk wave parameters from frequency-direction spectrum

    Parameters:
    -----------
    freq    : Vector of spectral frequencies [Hz]
    dirs    : Vector of directions (nautical convention) [Deg or Rad]
    spec    : Wave spectrum [m2/Hz/Deg or m2/Hz/Deg]
              The shape must be [freq,dirs]
    IGBand  : Infragravity wave frequency cutoff 
              defaults to 0.005-0.05 Hz
    zeroth  : If true will remove the first frequency from the analysis

    Returns:
    --------
    Dictionary containing bulk wave parameters
    Hs         : Significant wave height [m]
    H1         : Mean wave height [m]
    Tp         : Peak wave period [s]
    Tp_fit     : Peak wave period computed from second order polynomial fit
                 near Tp[s]
    Tm01       : First moment wave period [s]
    Tm02       : Second moment wave period [s]
    Te         : Energy period [s]
    Sw         : Spectral width (m0*m2/m1/m1 - 1)**2
    Tm01IG     : First moment wave period over the infragravity frequency band
    TpIG       : Peak wave period over the infragravity frequency band [s]
    HsIG       : Significant wave height in the infragravity frequency band [m]
    Dm         : Mean wave direction, second moment (Kuik et al. 1988)
    Dp         : Peak wave direction (computed from directional spectrum)
    
    Notes:
    ------
    - mn are the different spectral moments
    - First frequency will be discarded from the analysis. It is assumed to be
      the zeroth-frequency.

    REFERENCES:
    -----------
    Kuik, A.J., G.P. Van Vledder, and L.H. Holthuijsen, 1988: A method for the
        routine analysis of pitch-and-roll buoy wave data. Journal of Physical
        Oceanography, 18(7), pp.1020-1034.
    """
    
    # Get directional properties
    dirSpec = np.trapz(spec,freq,axis=0)
    
    # Find peak wave directon
    Dp = np.argmax(dirSpec)
    Dp = dirs[Dp]
    dirDict = dict(Dp=Dp)
    
    # Find mean wave direction (Kuik et al 1988)
    a1 = np.trapz(dirSpec*np.cos((270.0-dirs)*np.pi/180),dirs)
    b1 = np.trapz(dirSpec*np.sin((270.0-dirs)*np.pi/180),dirs)
    # Mean wave direction in nautical coordinates
    Dm = _gangles.wrapto360(270.0 - np.arctan2(b1,a1) * 180.0 / np.pi)
    dirDict.update(Dm=Dm)
    
    # Get the parameter from the frequency spectrum
    #freqSpec = np.trapz(spec,dirs,axis=-1)
    dth = np.abs(dirs[2] - dirs[1])
    freqSpec = np.sum(spec,axis=-1) * dth
    bp = fspec_bulk_params(freq,freqSpec,IGBand,zeroth=zeroth)
    bp.update(dirDict)
    
    return bp

#===============================================================================
# Wave Height and Period From time series
#===============================================================================
def whwpts(t,x,d='up'):
    '''
    Function to compute the wave height and wave period from a given time
    series.
    
    USAGE:
    ------
    [wh,wp,zcind] = whwpts(t,x,d)
    
    PARAMETERS:
    -----------
    t     : Time vector
    x     : Water surface elevation time series
    d     : Accepts 'up' or 'down' for zero upcrossing or downcrossing,
            respectively.
    
    RETURNS:
    --------
    wh    : Wave height time series
    wp    : Wave period time series
    zcind : zero-crossing index
    
    METHODS:
    --------
    1. Find zero crossings
    2. Find time difference between zero crossing to get wp
    3. Find the amplitude of the time series between each zero crossing to find
       wh.
    
    NOTES:
    ------
    - x and t must have the same length
    
    '''
    
    # Find the zero crossings
    zcross = gsignal.zero_crossing(x,d=d)
    
    # Find the wave period
    wp = t[zcross[1:]] - t[zcross[:-1]]
    
    # Find wave height
    wh = np.zeros_like(wp)
    for aa in range(wh.shape[0]):
        wh[aa] = (np.max(x[zcross[aa]:zcross[aa+1]]) -
                  np.min(x[zcross[aa]:zcross[aa+1]]))
        
    return wh,wp,zcross



#===============================================================================
# Iribarrren number
#===============================================================================
def iribarren(m,H,wl,verbose=False):
    """
    Function to compute the Iribarren number (aka surf similarity parameter) for 
    deep water conditions
    
    USAGE:
    ------
    ssp = iribarren(m,H,wl,verbose)
    
    PARAMETERS:
    -----------
    m       : Beach slope
    H       : Wave height [m]
    wl      : Wave length [m]
    verbose : Display result in the terminal (Optional defaults to False)
    
    RETURNS:
    --------
    ssp     : surf similarity parameter
    
    FORMULA:
    --------
    ssp = tan(m)/(H/wl)**0.5
       
    NOTES:
    ------
    - ssp > 3.3: Reflective beach (surging or collapsing breakers)
    - 0.5 < ssp < 3.3 : Intermediate beach (plunging breaker)
    - ssp < 0.5 : Dissipative beach (spilling breaker)
    
    """
    
    # Compute iribarren number
    ssp = np.tan(m) / (H/wl)**0.5
    
    # Display message
    if verbose:
        if ssp >= 3.3:
            print('ssp = ' + np.str(ssp) + ': Reflective conditions')
        elif ssp < 0.5:
            print('ssp = ' + np.str(ssp) + ': Dissipative conditions')
        else:
            print('ssp = ' + np.str(ssp) + ': Intermediate conditions')
    
    return ssp


#===============================================================================
# Battjes Parameter
#===============================================================================
def battjes04(m,h,igFreq,verbose=False):
    '''
    Function to compute the normalized bed slope proposed by Battjes et al 2004
    
    USAGE:
    ------
    nbd = battjes04(m,h,igFreq,verbose)
    
    PARAMETERS:
    -----------
    m       : Beach slope
    h       : Characteristic water depth in shoaling region [m]
    igFreq  : infragravity wave frequency [1/s]
    verbose : Display result in the terminal (Optional defaults to False)
    
    RETURNS:
    --------
    nbd     : Normalized bed slope
    
    FORMULA:
    --------
    nbd = m/(2 * pi * igFreq) * (g/h)**0.5
       
    NOTES:
    ------
    The next values are given if the breaking depth is substituted for h
    - nbd < 0.3: Mild slope regime, were breaking generation is not as important
    - nbd > 0.3: Steep slope regime
    
    References:
    -----------
    Battjes , J. A., H. J. Bakkenes, T. T. Janssen, and A. R. van Dongeren,
      2004: Shoaling of subharmonic gravity waves. Journal of Geophysical
      Research, 109, C02009, doi:10.1029/2003JC001863.
    
    '''
    
    # Compute iribarren number
    nbd = m / (2 * np.pi * igFreq) * (9.81/h)**0.5
    
    # Display message
    if verbose:
        if nbd < 0.3:
            print('nbd = ' + np.str(nbd) + ': Mild slope')
        else:
            print('nbd = ' + np.str(nbd) + ': Steep slope')
    
    return nbd


#===============================================================================
# Baldock Parameter
#===============================================================================
def baldock12(m,h,igFreq,H,wl,verbose=False):
    '''
    Function to compute the surfbeat similarity parameter by Baldock 2012
    
    USAGE:
    ------
    sbs = baldock12(m,h,igFreq,H,wl,verbose)
    
    PARAMETERS:
    -----------
    m       : Beach slope
    h       : Characteristic water depth in shoaling region [m]
    igFreq  : infragravity wave frequency [1/s]
    H       : Shortwave offshore wave height [m]
    wl      : Offshore wave length [m]
    verbose : Display result in the terminal (Optional defaults to False)
    
    RETURNS:
    --------
    sbs     : Surfbeat similarity parameter
    
    FORMULA:
    --------
    sbs = m/(2 * pi * igFreq) * (g/h)**0.5 * (H/wl)**0.5
       
    NOTES:
    ------
    Really no clear guidance is given, however the smallest the number the least
    effective breakpoint generation is.
    
    References:
    -----------
    Baldock, T. E., 2012: Dissipation of incident forced long waves in the 
      surf zone - Implications for the concept of "bound" wave realease at
      short wave breaking. Coastal Engineering, 276-285.

    
    '''
    
    # Compute iribarren number
    sbs = m / (2 * np.pi * igFreq) * (9.81/h)**0.5 * (H/wl) ** 0.5
    
    # Display message
    if verbose:
            print('sbs = ' + np.str(sbs))
    
    return sbs

#===============================================================================
# Reverse shoal waves
#===============================================================================
def deep_water_equivalent(H,h,period):
    """
    Code to find the equivalent deep water linear wave parameters based on
    intermediate to shallow water conditions
    
    PARAMETERS:
    -----------
    H      : Wave height [m]
    h      : Water depth [m]
    period : Wave period [s]
    
    RETURNS:
    --------
    Deep water equivalents
    H0     : Wave height [m]
    k0     : Wave number [m-1]
    
    NOTES:
    ------
    1. Will not work if waves have broken already
    
    """
    
    # Find the wave celerity
    c1,n1,cg1 = celerity(period,h)     # Intermediate water
    cg0 = (9.81*period) / (4 * np.pi)  # Deep water
    
    # Deep water wave length
    H0 = H * ((cg1/cg0)**0.5)
    
    # Deep water wave number
    k0 = ((2.0 * np.pi / period) ** 2) / 9.81
    
    return H0,k0


def shoal(H1,h1,period,h0):
    """
    Code to shoal waves based on linear wave theory
    
    PARAMETERS:
    -----------
    H1     : Wave height [m]
    h1     : Water depth [m]
    period : Wave period [s]
    h0     : Water depth to find equivalent conditions [m] 
    
    RETURNS:
    --------
    Equivalent parameters at h0
    H0     : Wave height [m]
    k0     : Wave number [m-1]
    
    NOTES:
    ------
    1. Will not work if waves have broken already
    
    """
    
    # Find the wave celerity
    c1,n1,cg1 = celerity(period,h1)     # Current conditions
    c0,n0,cg0 = celerity(period,h0)     # Desired conditions
    
    # Deep water wave length
    H0 = H1 * ((cg1/cg0)**0.5)
    
    # Deep water wave number
    k0 = dispersion(period,h0)
    
    return H0,k0


#===============================================================================
# Compute IEC Bulk Parameters from Spectrum
#===============================================================================
def iec_params(freq,dirs,spec,dpt,rho=1025.0):
    """
    Function to compute IEC recommended wave energy parameters from 
    frequency-direction spectrum

    Parameters:
    -----------
    freq    : Vector of spectral frequencies [Hz]
    dirs    : Vector of directions (nautical convention) [Deg]
    spec    : Wave spectrum [m2/Hz/Deg]
              The shape must be [time,npts,freq,dirs]
    dpt     : Water depth of shape [npts]
    rho     : (Optional) Density of water
    
    Returns:
    --------
    Dictionary containing bulk wave parameters
    Hs         : Significant wave height [m]
    OWP        : Omnidirectional wave power [W/m]
    Jth        : Maximum directionally resolved wave power [W/m]
    Th         : Direction of maximum directionally resolved wave power [deg]
    d          : Directionality coefficient (Jth/OWP)
    Te         : Energy period [s]
    Sw         : Spectral width (m0*m-2/m-1/m-1 - 1)**2
    
    Notes:
    ------
    - mn are the different spectral moments
    - Directional grid is assumed to be regular 
      need to implement a trapz for when this is not true.
    - Pass everything as arrays:
      >>> dpt = np.array([100.0]) # If you only have one point

    REFERENCES:
    -----------
    Marine energy - Wave, tidal and other water current converters - 
        Part 101: Wave energy resoure assessment and characterization. 
        IEC TS 62600-101
    """
    
    # Figure out the size of the array and adjust if necessary
    _,npts,nfreq,ndirs = np.shape(spec)

    # Directions must be sorted (if they exist)
    if ndirs > 1:
        sortInd = np.argsort(dirs)
        spec = spec[...,sortInd]
        dirs = dirs[sortInd]

    # Compute the group velocity for each point
    cg = np.zeros((npts,nfreq)) * np.NAN
    for aa in range(cg.shape[0]):
        for bb in range(cg.shape[1]):
            _,_,cg[aa,bb] = celerity(1/freq[bb],dpt[aa])

    # Non directional parameters -----------------------------------------------    
    if ndirs > 1:
        #freqSpec = np.trapz(spec,dirs*np.pi/180,axis=-1)
        dth = dirs[1] - dirs[0]
        freqSpec = np.sum(spec,axis=-1)*dth
    else:
        freqSpec = spec[...,0]
    
    # Compute omnidirectional wave power
    moment0 = np.trapz(freqSpec,freq,axis=-1)
    momentn1 = np.trapz(freqSpec*(freq)**-1,freq,axis=-1)
    momentn2 = np.trapz(freqSpec*(freq)**-2,freq,axis=-1)

    # Significant wave height
    bp = {}
    bp['Hs'] = 4.004 * (moment0**0.5)

    # Energy Period
    bp['Te'] = momentn1/moment0
    
    # Spectral width
    bp['Sw'] = (moment0*momentn2/momentn1/momentn1 - 1.0)**0.5

    # Omnidirectional wave power (loop extravaganza here)    
    owp = np.zeros_like(bp['Hs']) * np.NAN
    # Time loop
    for aa in range(spec.shape[0]):
        # Point loop
        for bb in range(spec.shape[1]):
            # Omnidirectional wave power
            owp[aa,bb] = np.trapz(freqSpec[aa,bb,:]*cg[bb,:],freq)

    # Get directional properties -----------------------------------------------

    if ndirs > 1:
        # Omidirectional wave power through a plane
        cg  = np.repeat(np.expand_dims(cg,axis=-1),ndirs,axis=-1)
        jth = np.zeros_like(owp)
        th  = np.zeros_like(owp) * np.NAN
        # Time loop
        for aa in range(spec.shape[0]):
            # Point loop
            for bb in range(spec.shape[1]):
                # Directional wave power
                tmpJth = 0.0
                dirSpec = np.trapz(spec[aa,bb,...]*cg[bb,...],freq,axis=-2)
                # Direction loop
                for cc in range(ndirs):                
                    fac = np.cos(np.pi/180.0 * (dirs[cc] - dirs))
                    fac[fac<0] = 0.0
                    tmpJth = np.sum(dirSpec*fac,axis=-1) * dth
                    if tmpJth > jth[aa,bb]:
                        jth[aa,bb] = tmpJth
                        th[aa,bb]  = dirs[cc]

    else:
        th  = np.zeros_like(owp) * np.NAN
        jth = np.zeros_like(owp) * np.NAN

    # Take care of units here
    owp *= 9.81 * rho
    jth *= 9.81 * rho

    # Allocate in arrays
    bp['Th']  = th
    bp['OWP'] = owp
    bp['Jth'] = jth
    bp['d']   = jth/owp      # Directionality coefficient
    
    return bp
