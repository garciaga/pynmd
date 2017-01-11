"""
Collection of codes to process funwaveC output data
"""

from __future__ import division,print_function

import sys
import os
import glob
import time
import numpy as np
import netCDF4
from collections import defaultdict

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


# ==============================================================================
# Create NetCDF file
# ==============================================================================
def convert_output(workfld,outfile,time_int=1.0,bathyfile=None,
                   inpfile=None,verbose=False):
    '''
    
    Parameters:
    -----------
    workfld      : Path to the folder where output files reside
    outfile      : Output NetCDF file
    time_int     : Time interval between output files [s] (will be updated if
                   input file is provided)
    bathyfile    : Full path to input netcdf bathy file (optional)
    inpfile      : Funwave input file used for metadata (optional)
    verbose      : Defaults to False
         
    Output:
    -------
    NetCDF File with the variables in the folder. Not all are supported so you 
    may need to edit this file.
    
    '''
    
    # For testing only
    #workfld = '/scratch/temp/ggarcia/funwaveC_ensemble/r0/'
    #outfile = '/scratch/temp/ggarcia/funwaveC_ensemble/02-runs/tmp.nc'
    #time_int = 0.5
    #bathyfile = '/scratch/temp/ggarcia/funwaveC_ensemble/r0/depth.nc'
    #inpfile = '/scratch/temp/ggarcia/funwaveC_ensemble/r0/input.init'

    # Get variable information ------------------------------------------------
    # Need to generalize
    archivos = glob.glob(workfld + '*.dat')                       # All files
    tmpvars = [x.split('/')[-1].split('.')[0] for x in archivos]  # All vars
    
    # If no variables found exit    
    if not tmpvars:
        print('No files found in ' + workfld)
        print('Quitting ...')
        return None
    
    # Make sure variables are within the supported ones
    supported_vars_time = ['eta','u']
    vars_2d = [x for x in tmpvars if x in supported_vars_time]
    
    
    # Read input file if provided ----------------------------------------------
    if inpfile:
        
        if verbose:
            print("Reading data from input file:")
            print("  " + inpfile)
        
        # Open file
        tmpinpfile = open(inpfile,'r')
        
        # Output dictionary
        inpinfo = {}
        
        # funwaveC is very structured so this is a simple way of reading
        # Dynamics line
        tmpLine = tmpinpfile.readline().rstrip()
        inpinfo['dynamics'] = tmpLine.split(' ')[-1]
        # Dimensions
        tmpLine = tmpinpfile.readline().rstrip()
        dx = np.float64(tmpLine.split(' ')[3])
        dy = np.float64(tmpLine.split(' ')[4])
        # Bottom stress
        tmpLine = tmpinpfile.readline().rstrip()
        inpinfo['bottomstress'] = tmpLine.split(' ')[2]
        # Mixing
        tmpLine = tmpinpfile.readline().rstrip()
        inpinfo['mixing'] = tmpLine.split(' ')[1] + ' ' + tmpLine.split(' ')[2]
        # Bathymetry
        tmpinpfile.readline()
        # Tide
        tmpinpfile.readline()
        # Wavemaker
        tmpLine = tmpinpfile.readline().rstrip()        
        inpinfo['eta_source'] = tmpLine[14:]
        # Wave breaking
        tmpLine = tmpinpfile.readline().rstrip()        
        inpinfo['breaking'] = tmpLine[9:]
        # Sponge layer
        tmpLine = tmpinpfile.readline().rstrip()        
        inpinfo['sponge'] = tmpLine[7:]
        # Forcing, initial conditions, tracers
        tmpinpfile.readline()
        tmpinpfile.readline()
        tmpinpfile.readline()
        tmpinpfile.readline()
        tmpinpfile.readline()
        tmpinpfile.readline()
        # Timing
        tmpLine = tmpinpfile.readline().rstrip()
        inpinfo['timing'] = tmpLine[7:]
        time_int = np.float64(tmpLine.split(' ')[5])
        
        # Close me
        tmpinpfile.close()
        
            
    else:
        
        # Assume grid spacing to be 1 meter
        if verbose:
            print("Input file not provided")
        dx = 1.0
        dy = 1.0
        inpinfo = False
    
    
    # Coordinates and depth ---------------------------------------------------
    # Bathymetry file provided
    if bathyfile: 
        
        if verbose:
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
       
        # Create u grid
        x_u       = np.zeros((x_rho.shape[0]+1,))
        x_u[0]    = x_rho[0] - dx/2
        x_u[-1]   = x_rho[-1] + dx/2
        x_u[1:-1] = (x_rho[1:] + x_rho[:-1])/2.0
        
    # Bathymetry file not provided but have output file        
    elif os.path.isfile(workfld + '/depth.txt'):
               
        # Fix this               
        h = np.loadtxt(workfld + '/depth.txt')
        
        hdims = h.ndim
        if hdims == 1:
            x_rho = np.arange(0,h.shape[0],dx)
        elif hdims == 2:
            x_rho, y_rho = np.meshgrid(np.arange(0,h.shape[1],dx),
                                       np.arange(0,h.shape[0],dy))           
        else:
            if verbose:
                print('Something is wrong with the depth file')
    
    # No bathymetry file provided
    else:
        
        if verbose:
            print("No bathymetry file provided")

    # Get dimensions of variables
    hdims = h.ndim
    
    # Create NetCDF file -------------------------------------------------------
    if verbose:
        print("Creating " + outfile)
    
    # Global attributes  
    nc = netCDF4.Dataset(outfile, 'w', format='NETCDF4')
    nc.Description = 'FunwaveC Output'
    nc.Author = 'ggarcia@coas.oregonstate.edu'
    nc.Created = time.ctime()
    nc.Type = 'FunwaveC snapshot output'
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
        xi_u   = xi_rho + 1
        eta_rho = 1
    nc.createDimension('xi_rho', xi_rho)
    nc.createDimension('xi_u',xi_u)            
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

        nc.createVariable('x_u','f8',('xi_u'))
        nc.variables['x_u'].units = 'meter'
        nc.variables['x_u'].longname = 'x-locations of U points'
        nc.variables['x_u'][:] = x_u        
        
    # Create time vector -------------------------------------------------------
    tmpvar = np.loadtxt(workfld + vars_2d[0] + '.dat')
    twave  = np.arange(time_int,time_int*tmpvar.shape[0]+time_int,time_int)
            
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
    varinfo['eta']['dims'] = ('ocean_time','xi_rho')
    
    varinfo['u']['units'] = 'meter second-1'
    varinfo['u']['longname'] = 'Flow velocity in the xi direction'
    varinfo['u']['dims'] = ('ocean_time','xi_u')
    
    varinfo['v']['units'] = 'meter second-1'
    varinfo['v']['longname'] = 'Flow velocity in the eta direction'
    varinfo['v']['dims'] = ('ocean_time','xi_v')
         
    if verbose:          
        print("Creating variables")
                  
    for aa in vars_2d:
        
        if verbose:
            print('  ' + aa)
        
                # Load variable
        try:
            tmpvar = np.loadtxt(workfld + '/' + aa + '.dat')
        except:
            if verbose:
                print("  Could not create " + aa)
                print("  Input file error")
            continue
        
        # Create variable
        create_nc_var(nc,aa,varinfo[aa]['dims'],varinfo[aa]['units'],
                      varinfo[aa]['longname'])        

        # Write variable
        nc.variables[aa][:] = tmpvar
        
                    
    # Close NetCDF file
    if verbose:
        print('Closing ' + outfile)
    nc.close()
    
    # End of function


#===============================================================================
# Compute Runup
#===============================================================================
def runup(eta,h,x,r_depth=None):
    """
    
    Parameters:
    ----------
    eta          : Water surface elevation time series [m]
    h            : Bathymetry [m] (positive down)
    x            : x coordinates of h [m]
    r_depth      : Runup depth [m] (optional)
    
    Output:
    -------
    runup        : Water surface elevation time series relative to SWL given
                   a contour depth [m]
    x_runup      : Across-shore location of runup time series [m]
                   
    Notes:
    ------
    - Really meant for 1D simulations.
    - Minimum runup depth as suggested by Salmon 2002 is estimated as
      h0 = 4 * max(dh)
                   
    """

    # Find the maximum runup depth
    if r_depth is None:
         runupInd = h < 1.0
         r_depth = 4.0*np.nanmax(np.abs(h[runupInd][1:] - h[runupInd][:-1])) 
             
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
    return runup,x_runup,r_depth



def nc_runup(nc,r_depth=None):
    """
    
    Function to compute runup from netcdf file.
    
    runup,x_runup = nc_runup(nc,r_depth)
    
    Parameters:
    -----------
    nc           : NetCDF file handle
    r_depth      : Runup depth [m] (optional)
    
    Output:
    -------
    runup        : Water surface elevation time series relative to SWL given
                   a tolerance depth [m]
    x_runup      : Across-shore location of runup time series [m]
                   
    Notes:
    ------
    -  Really meant for 1D simulations.
    -  See convert_output for ASCII to NetCDF4 file conversion.
    -  Minimum runup depth as suggested by Salmon 2002 is estimated as
       h0 = 4 * max(dh)
    
    """
    
    # Load variables
    h = nc.variables['h'][:]
    ot = nc.variables['ocean_time'][:]
    x = nc.variables['x_rho'][:]
    
    # Find the maximum runup depth
    if r_depth is None:
         runupInd = h < 1.0
         r_depth = 4.0*np.nanmax(np.abs(h[runupInd][1:] - h[runupInd][:-1])) 

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
    return runup,x_runup,r_depth

