"""
Tools to manage funwave input

Authors:
--------
Gabriel Garcia Medina (ggarcia@coas.oregonstate.edu)
Saeed Moghimi (moghimi@coas.oregonstate.edu)

# Last edit
17 September 2015 - Gabriel Garcia Medina

External dependencies:
  numpy
  scipy
  netCDF4
  time
  sys
  getpass
  os
  collections
  imp

Internal dependencies:
  gsignal
  
"""

from __future__ import division,print_function

__author__ = "Gabriel Garcia Medina"
__email__ = "ggarcia@coas.oregonstate.edu"
__group__="Nearshore Modeling Group"


# Import Modules
import netCDF4
import time
import sys
import getpass
import os
import numpy as np
import pylab as pl
from collections import defaultdict
import imp

# Custom paths
import pynmd.physics.waves as gwaves

#===============================================================================
# Vorticity
#===============================================================================

def vorticity(x,y,u,v):
    '''
    
    Parameters:
    -----------
    x: 2d array of x locations of velocity points
    y: 2d array of y locations of velocity points
    u: 2d array of flow velocity in the x direction
    v: 2d array of flow velocity in the y direction3
    
    Output:
    -------
    q: 2d array of vertical vorticity (dv/dx - du/dy)
    
    '''
    
    # Compute gradients
    tmpgrad = np.gradient(x)
    dx = tmpgrad[1]
    tmpgrad = np.gradient(y)
    dy = tmpgrad[0]

    tmpgrad = np.gradient(u)
    du = tmpgrad[0]
    tmpgrad = np.gradient(v)
    dv = tmpgrad[1]
    
    # Compute vorticity
    return (dv/dx - du/dy)



# ==================================================================
# Create NetCDF file
# ==================================================================    
def convert_output(workfld,outfile,time_int,bathyfile=None,inpfile=None):
    '''
    
    Parameters:
    -----------
    workfld      : Path to the folder where output files reside
    outfile      : Output NetCDF file
    time_int     : Time interval between output files [s] (will be updated if
                   input file is provided)
    bathyfile    : Full path to input netcdf bathy file (optional)
    inpfile      : Funwave input file used for metadata (optional)
         
    Output:
    -------
    NetCDF File with the variables in the folder. Not all are supported so you 
    may need to edit this file.
    
    '''
    
    # For testing only
    #workfld = '/scratch/temp/ggarcia/spectralWidthFunwave/02-runs/results/'
    #outfile = '/scratch/temp/ggarcia/spectralWidthFunwave/02-runs/test.nc'
    #time_int = 0.1
    #bathyfile = '/scratch/temp/ggarcia/spectralWidthFunwave/02-runs/depth.nc'
    #inpfile = '/scratch/temp/ggarcia/spectralWidthFunwave/02-runs/input.txt'

    # Get variable information ------------------------------------------------
    archivos = os.listdir(workfld)                      # Get all files
    tmpvars = [x.split('_')[0] for x in archivos]       # All variables
    tmpvars = list(set(tmpvars))                        # Unique variables
    
    # If no variables found exit    
    if not tmpvars:
        print('No files found in ' + workfld)
        print('Quitting ...')
        return None
    
    # Make sure variables are within the supported ones
    supported_vars_time = ['eta','etamean','havg','hmax','hmin','hrms',
                           'mask','mask9','MFmax','u','umax','umean','v',
                           'vmean','VORmax']
    vars_2d = [x for x in tmpvars if x in supported_vars_time]
    
    
    # Read input file if provided ----------------------------------------------
    if inpfile:
        
        print("Reading data from input file:")
        print("  " + inpfile)
        
        # Open file
        tmpinpfile = open(inpfile,'r')
        
        # Output dictionary
        inpinfo = {}
        
        # Extract information (need to add wavemaker support)
        for tmpline in tmpinpfile:
            
            # Skip blank lines
            if len(tmpline.strip()) == 0:
                continue
            
            if "TITLE" == tmpline.split()[0]:
                inpinfo['title'] = tmpline.split()[2]
            elif "PX" == tmpline.split()[0]:
                inpinfo['px'] = tmpline.split()[2]
            elif "PY" == tmpline.split()[0]:
                inpinfo["py"] = tmpline.split()[2]
            elif "Mglob" == tmpline.split()[0]:
                inpinfo['mglob'] = tmpline.split()[2]
            elif "Nglob" == tmpline.split()[0]:
                inpinfo['nglob'] = tmpline.split()[2]
            elif "TOTAL_TIME" == tmpline.split()[0]:
                inpinfo['total_time'] = tmpline.split()[2]
            elif "PLOT_INTV" == tmpline.split()[0]:
                inpinfo["plot_intv"] = tmpline.split()[2]
                time_int = np.double(tmpline.split()[2])
            elif "DX" == tmpline.split()[0]:
                inpinfo['dx'] = tmpline.split()[2]
                dx = float(tmpline.split()[2])
            elif "DY" == tmpline.split()[0]:
                inpinfo['dy'] = tmpline.split()[2]
                dy = float(tmpline.split()[2])
            elif "WAVEMAKER" == tmpline.split()[0]:
                inpinfo['wavemaker'] = tmpline.split()[2]
            elif "PERIODIC" == tmpline.split()[0]:
                inpinfo['periodic'] = tmpline.split()[2]
            elif "SPONGE_ON" == tmpline.split()[0]:
                sponge = tmpline.split()[2]
                inpinfo['sponge'] = sponge
            elif "Sponge_west_width" == tmpline.split()[0] and sponge == 'T':
                inpinfo['sponge_west_width'] = tmpline.split()[2]
            elif "Sponge_east_width" == tmpline.split()[0] and sponge == 'T':
                inpinfo['sponge_east_width'] = tmpline.split()[2]
            elif "Sponge_wouth_width" == tmpline.split()[0] and sponge == 'T':
                inpinfo['sponge_south_width'] = tmpline.split()[2]
            elif "Sponge_worth_width" == tmpline.split()[0] and sponge == 'T':
                inpinfo['sponge_north_width'] = tmpline.split()[2]
            elif "R_sponge" == tmpline.split()[0] and sponge == 'T':
                inpinfo['r_sponge'] = tmpline.split()[2]
            elif "A_sponge" == tmpline.split()[0] and sponge == 'T':
                inpinfo['a_sponge'] = tmpline.split()[2]
            elif "DISPERSION" == tmpline.split()[0]:
                inpinfo['dispersion'] = tmpline.split()[2]
            elif "Gamma1" == tmpline.split()[0]:
                inpinfo['gamma1'] = tmpline.split()[2]
            elif "Gamma2" == tmpline.split()[0]:
                inpinfo['gamma2'] = tmpline.split()[2]
            elif "Gamma3" == tmpline.split()[0]:
                inpinfo['gamma3'] = tmpline.split()[2]
            elif "Beta_ref" == tmpline.split()[0]:
                inpinfo['beta_ref'] = tmpline.split()[2]
            elif "SWE_ETA_DEP" == tmpline.split()[0]:
                inpinfo['swe_eta_dep'] = tmpline.split()[2]
            elif "Friction_Matrix" == tmpline.split()[0]:
                inpinfo['friction_matrix'] = tmpline.split()[2]
            elif "Cd_file" == tmpline.split()[0]:
                inpinfo['cd_file'] = tmpline.split()[2]
            elif "Cd" == tmpline.split()[0]:
                inpinfo['cd'] = tmpline.split()[2]                              
            elif "Time_Scheme" == tmpline.split()[0]:
                inpinfo['time_scheme'] = tmpline.split()[2]
            elif "HIGH_ORDER" == tmpline.split()[0]:
                inpinfo['spatial_scheme'] = tmpline.split()[2]
            elif "CONSTRUCTION" == tmpline.split()[0]:
                inpinfo['construction'] = tmpline.split()[2]
            elif "CFL" == tmpline.split()[0]:
                inpinfo['cfl'] = tmpline.split()[2]                
            elif "MinDepth" == tmpline.split()[0]:
                inpinfo['min_depth'] = tmpline.split()[2]
            elif "MinDepthFrc" == tmpline.split()[0]:
                inpinfo['mindepthfrc'] = tmpline.split()[2]
            
        # Close file
        tmpinpfile.close()
            
    else:
        
        # Assume grid spacing to be 1 meter
        print("Input file not provided")
        dx = 1
        dy = 1
        inpinfo = False
    
    
    # Coordinates and depth ---------------------------------------------------
    # Bathymetry file provided
    if bathyfile: 
        
        print("Reading coordinates and depth from:")
        print("  " + bathyfile)
        
        ncfile = netCDF4.Dataset(bathyfile,'r')
        x_rho = ncfile.variables['x_rho'][:]
        if ncfile.variables.has_key('y_rho'):
            y_rho = ncfile.variables['y_rho'][:]
        else:
            y_rho = None
        h = ncfile.variables['h'][:]
        ncfile.close()
        
    # Bathymetry file not provided but have output file        
    elif os.path.isfile(workfld + '/dep.out'):
               
        # Fix this               
        h = np.loadtxt(workfld + '/dep.out')
        
        hdims = h.ndim
        if hdims == 1:
            x_rho = np.arange(0,h.shape[0],dx)
        elif hdims == 2:
            x_rho, y_rho = np.meshgrid(np.arange(0,h.shape[1],dx),
                                       np.arange(0,h.shape[0],dy))           
        else:
            print('Something is wrong with the depth file')
            return None    
    
    # No bathymetry file provided (I am not capable of reading the input file)
    else:
        print("No bathymetry file provided")
        print("You could copy your input bathymetry text file to ")
        print(workfld + '/dep.out')
        return None
        
    
    # Get dimensions of variables
    hdims = h.ndim
    
    # Create NetCDF file -------------------------------------------------------
    
    print("Creating " + outfile)
    
    # Global attributes  
    nc = netCDF4.Dataset(outfile, 'w', format='NETCDF4')
    nc.Description = 'Funwave Output'
    nc.Author = 'ggarcia@coas.oregonstate.edu'
    nc.Created = time.ctime()
    nc.Type = 'Funwave v2.1 snapshot output'
    nc.Owner = 'Nearshore Modeling Group'
    nc.Software = 'Created with Python ' + sys.version
    nc.NetCDF_Lib = str(netCDF4.getlibversion())
    nc.Source = workfld
    nc.Script = os.path.realpath(__file__)
    
    # Add more global variables to output
    if inpinfo:
        for tmpatt in inpinfo.keys():
            nc.__setattr__(tmpatt,inpinfo[tmpatt][:])
    
    # Create dimensions
    if hdims == 2:
        eta_rho, xi_rho = h.shape
        nc.createDimension('eta_rho', eta_rho)
    else:
        xi_rho = h.shape[0]
        eta_rho = 1
    nc.createDimension('xi_rho', xi_rho)            
    nc.createDimension('ocean_time',0)
    
    # Write coordinate axes ----------------------------------------------
    if hdims == 2:
        
        nc.createVariable('x_rho','f8',('eta_rho','xi_rho'))
        nc.variables['x_rho'].units = 'meter'
        nc.variables['x_rho'].longname = 'x-locations of RHO points'
        nc.variables['x_rho'][:] = x_rho
        
        nc.createVariable('y_rho','f8',('eta_rho','xi_rho'))
        nc.variables['y_rho'].units = 'meter'
        nc.variables['y_rho'].longname = 'y-locations of RHO points'
        nc.variables['y_rho'][:] = y_rho

        nc.createVariable('h','f8',('eta_rho','xi_rho'))
        nc.variables['h'].units = 'meter'
        nc.variables['h'].longname = 'bathymetry at RHO points'
        nc.variables['h'][:] = h


    else:
        
        nc.createVariable('x_rho','f8',('xi_rho'))
        nc.variables['x_rho'].units = 'meter'
        nc.variables['x_rho'].longname = 'x-locations of RHO points'
        nc.variables['x_rho'][:] = x_rho
        
        nc.createVariable('h','f8',('xi_rho'))
        nc.variables['h'].units = 'meter'
        nc.variables['h'].longname = 'bathymetry at RHO points'
        nc.variables['h'][:] = h
        
           
    # Create time vector -------------------------------------------
    tmpruns = [x for x in archivos if x.split('_')[0] == vars_2d[0]]
    tmpruns.sort()
    time_max = float(tmpruns[-1].split('_')[-1])
    time_min = float(tmpruns[0].split('_')[-1])
    
    
    twave = np.arange(time_min-1,time_int*(time_max-time_min+1)+time_min-1,
                      time_int)
    if twave.shape[0] > len(tmpruns):
        twave = twave[:-1]
        
    nc.createVariable('ocean_time','f8','ocean_time')
    nc.variables['ocean_time'].units = 'seconds since 2000-01-01 00:00:00'
    nc.variables['ocean_time'].calendar = 'julian'
    nc.variables['ocean_time'].long_name = 'beach time'
    nc.variables['ocean_time'][:] = twave
        
    # Create variables --------------------------------------------------------

    # Variable information
    varinfo = defaultdict(dict)

    varinfo['eta']['units'] = 'meter'
    varinfo['eta']['longname'] = 'water surface elevation'
    varinfo['etamean']['units'] = 'meter'
    varinfo['etamean']['longname'] = 'Mean wave induced setup'
    
    varinfo['u']['units'] = 'meter second-1'
    varinfo['u']['longname'] = 'Flow velocity in the xi direction'
    varinfo['v']['units'] = 'meter second-1'
    varinfo['v']['longname'] = 'Flow velocity in the eta direction'
    
    varinfo['umean']['units'] = 'meter second-1'
    varinfo['umean']['longname'] = 'Time-averaged flow velocity in xi direction'
    varinfo['vmean']['units'] = 'meter second-1'
    varinfo['vmean']['longname'] = 'Time-averaged flow velocity in eta direction'
    
    varinfo['umax']['units'] = 'meter second-1'
    varinfo['umax']['longname'] = 'Maximum flow velocity in xi direction'
    varinfo['vmax']['units'] = 'meter second-1'
    varinfo['vmax']['longname'] = 'Maximum flow velocity in eta direction'    
    
    varinfo['hmax']['units'] = 'meter'
    varinfo['hmax']['longname'] = 'Maximum wave height'
    varinfo['hmin']['units'] = 'meter'
    varinfo['hmin']['longname'] = 'Minimum wave height'
    varinfo['havg']['units'] = 'meter'
    varinfo['havg']['longname'] = 'Average wave height'
    varinfo['hrms']['units'] = 'meter'
    varinfo['hrms']['longname'] = 'Root mean squared wave height'
    
    varinfo['mask']['units'] = 'Boolean'
    varinfo['mask']['longname'] = 'Logical parameter for output wetting-drying'
    varinfo['mask9']['units'] = 'Boolean'
    varinfo['mask9']['longname'] = 'Logical parameter for output MASK9'
    
    varinfo['VORmax']['units'] = 'second-1'
    varinfo['VORmax']['longname'] = 'Maximum vorticity'
    varinfo['MFmax']['units'] = 'meter second-s'
    varinfo['MFmax']['longname'] = 'Maximum momentum flux'    
    
    
    # Create variables        
    if hdims == 1:
        nc_dims = ('ocean_time','xi_rho')
    else:
        nc_dims = ('ocean_time','eta_rho','xi_rho')
          
          
    print("Creating variables")          
    for aa in vars_2d:
        
        print('  ' + aa)
        
        # Create variable
        create_nc_var(nc,aa,nc_dims,varinfo[aa]['units'],
                      varinfo[aa]['longname'])
        
        try:
            tmpvar = np.loadtxt(workfld + '/' + aa + '_' + '%05.0f' % time_min)
        except ValueError:
            tmpvar = np.zeros_like(h) * np.NAN
                    
        nc.variables[aa][:] = np.expand_dims(tmpvar,axis=0)
        
        for bb in range(len(twave)):
            
            try:
                tmpvar = np.loadtxt(workfld + '/' + aa + '_' + 
                                    '%05.0f' % (bb + time_min))
            except ValueError:
                tmpvar = np.zeros_like(h) * np.NAN
                
            append_nc_var(nc,tmpvar,aa,bb-1)   
                    
    # Close NetCDF file
    print('Closing ' + outfile)
    nc.close()
    
    # End of function



#===============================================================================
# Compute mean setup
#===============================================================================
def setup(runup,ot):
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



#===============================================================================
# Compute Runup
#===============================================================================
def runup(eta,h,x,r_depth=0.01):
    """
    
    Parameters:
    ----------
    eta          : Water surface elevation time series [m]
    h            : Bathymetry [m] (positive down)
    x            : x coordinates of h [m]
    r_depth      : Runup depth [m] (defaults to 0.01m)
    
    Output:
    -------
    runup        : Water surface elevation time series relative to SWL given
                   a contour depth [m]
    x_runup      : Across-shore location of runup time series [m]
                   
    Notes:
    ------
    Really meant for 1D simulations.
                   
    """
    
    # Preallocate runup variable
    runup = np.zeros(eta.shape[0])
    x_runup = np.zeros_like(runup)
    
    # Loop over time, which is assumed to be on the first dimension.
    for aa in range(runup.shape[0]):
        
        # Water depth
        wdepth = eta[aa,:] + h
        
        # Find the runup contour (search from left to right) 
        wdepth_ind = np.argmin(wdepth>r_depth)-1
        
        # Store the water surface elevation in matrix
        runup[aa]= eta[aa,wdepth_ind]
        
        # Store runup position
        x_runup[aa] = x[wdepth_ind] 
    
    # Done
    return runup,x_runup



def nc_runup(nc,r_depth=0.01):
    """
    
    Function to compute runup from netcdf file.
    
    runup,x_runup = nc_runup(nc,r_depth)
    
    Parameters:
    -----------
    nc           : NetCDF file handle
    r_depth      : Runup depth [m] (defaults to 0.01m)
    
    Output:
    -------
    runup        : Water surface elevation time series relative to SWL given
                   a tolerance depth [m]
    x_runup      : Across-shore location of runup time series [m]
                   
    Notes:
    ------
    -  Really meant for 1D simulations.
    -  See convert_output for ASCII to NetCDF4 file conversion.
    
    """
    
    # Load variables
    h = nc.variables['h'][:]
    ot = nc.variables['ocean_time'][:]
    x = nc.variables['x_rho'][:]
    
    # Preallocate runup variable
    runup = np.zeros(ot.shape[0])
    x_runup = np.zeros_like(runup)
    
    # Loop over time, which is assumed to be on the first dimension.
    for aa in range(runup.shape[0]):
        
        # Water depth
        eta = nc.variables['eta'][aa,:]
        wdepth = eta + h
        
        # Find the runup contour (search from left to right) 
        wdepth_ind = np.argmin(wdepth>r_depth)-1
        
        # Store the water surface elevation in matrix
        runup[aa] = eta[wdepth_ind]
        
        # Store runup position
        x_runup[aa] = x[wdepth_ind] 
        
    
    # Done
    return runup,x_runup



#===============================================================================
# Compute 1D frequency spectrum at a point
#===============================================================================
def freq_spec_1d(eta,dt=1,verbose=True):
    """
    Computes the frequency spectrum from a given time series.
    
    freq,spec = freq_spec_1d(eta,dt,verbose)
    
    PARAMETERS:
    -----------
    eta      : Time series of water surface elevation [m]
    dt       : Time step [s]
    verbose  : Display computed bulk parameters to the screen
    
    RETURNS:
    --------
    freq     : Frequency vector
    spec     : Variance spectrum (Power spectrum)
    
    NOTES:
    ------
    This is really a copy of gsignal.psdraw. If results differ, trust gsignal
      this code will not be updated.
    """
    
    # Remove mean
    eta -= eta.mean()

    # Compute record length
    N = eta.shape[0]
    
    # Compute fourier frequencies
    fj = np.fft.fftfreq(N,dt)

    # Compute power spectral density (Cooley-Tukey Method)
    yf = np.fft.fft(eta)/N
    psd = N*dt*yf*np.conjugate(yf)

    # One sided psd from dft
    if np.mod(N,2) == 0:
        sf = np.concatenate((np.array([psd[0]]),2.0*psd[1:N/2],
                             np.array([psd[N/2]])))
        freq_amp = np.abs(np.concatenate((np.array([fj[0]]),fj[1:N/2],
                                          np.array([fj[N/2]]))))
    else:
        sf = np.concatenate((np.array([psd[0]]),2.0*psd[1:(N+1)/2]))
        freq_amp = np.abs(np.concatenate((np.array([fj[0]]),fj[1:(N+1)/2])))

    sf = sf.real
    
    if verbose:
        print("===============================================")
        print("Bulk Wave Parameters:")
        print("Hs = " + np.str(4.004*np.sqrt(np.trapz(sf,freq_amp))) + "m")
        print("H1 = "+np.str(4.004*np.sqrt(np.trapz(sf,freq_amp))*2.0/3.0)+"m")
        print("Spectral Parameters:")
        print("Nyquist Frequency = " + np.str(1.0/(2.0*dt)) + "Hz")
        print("Frequency interval = " + np.str(1.0/(N*dt)) + "Hz")
        print("===============================================")

    # End of function
    return freq_amp,sf



#===============================================================================
# Pyroms subroutine to write NetCDF fields
#===============================================================================
def create_nc_var(nc, name, dimensions, units=None, longname=None):
    '''
    Not for standalone use
    '''
    nc.createVariable(name, 'f8', dimensions)
    if units is not None:
        nc.variables[name].units = units
    if longname is not None:
        nc.variables[name].long_name = longname    

# Append NetCDF variable        
def append_nc_var(nc,var,name,tstep):
    '''
    Not for standalone use
    '''
    nc.variables[name][tstep,...] = var

