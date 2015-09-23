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

# Custom paths
import pynmd.physics.waves as gwaves


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


# ==================================================================
# Create NetCDF file
# ==================================================================  
def write_bathy(x,y,h,path,ncsave=True):
    '''
    
    Parameters:
    ----------
    x,y         : 2D arrays of coordinates
    h           : Bathymetry
    path        : Full path where the output will be saved
    ncsave      : Save as NetCDF4 file (optional, defaults to True)
    
    Output:
    -------
    depth.nc    : NetCDF4 file with the bathymetry
    depth.txt   : Text file with the depth information for Funwave input.
    
    '''
    
    # Get bathymetry dimensions
    eta_rho, xi_rho = x.shape
    
    if ncsave:
    
        # Global attributes  
        nc = netCDF4.Dataset(path + 'depth.nc', 'w', format='NETCDF4')
        nc.Description = 'Funwave Bathymetry'
        nc.Author = getpass.getuser()
        nc.Created = time.ctime()
        nc.Owner = 'Nearshore Modeling Group (http://ozkan.oce.orst.edu/nmg)'
        nc.Software = 'Created with Python ' + sys.version
        nc.NetCDF_Lib = str(netCDF4.getlibversion())
        nc.Script = os.path.realpath(__file__)
     
        # Create dimensions
        nc.createDimension('xi_rho', xi_rho)
        nc.createDimension('eta_rho', eta_rho)
    
        # Write coordinates and depth to netcdf file
        create_nc_var(nc, 'x_rho', 
                     ('eta_rho', 'xi_rho'), 
                     'meter','x-locations of RHO-points')
        nc.variables['x_rho'][:] = x
        create_nc_var(nc, 'y_rho', 
                     ('eta_rho', 'xi_rho'), 
                     'meter','y-locations of RHO-points')
        nc.variables['y_rho'][:] = y
        create_nc_var(nc, 'h', 
                     ('eta_rho', 'xi_rho'), 
                     'meter','bathymetry at RHO-points')          
        nc.variables['h'][:] = h
                        
        # Close NetCDF file
        nc.close()

    else: 
    
        print("NetCDF file not requested")
        
        

    # Output the text file -----------------------------------------------------        
    fid = open(path + 'depth.txt','w')
    for aa in range(x.shape[0]):
        for bb in range(x.shape[1]):
            fid.write('%12.3f' % h[aa,bb])
        fid.write('\n')
    fid.close()



    #===========================================================================
    # Print input file options
    #===========================================================================
    print(' ')
    print('===================================================================')
    print('In your FUNWAVE input file:')
    print('Mglob = ' + np.str(xi_rho))
    print('Nglob = ' + np.str(eta_rho))
    print('DX = ' + np.str(np.abs(x[0,1] - x[0,0])))
    print('DY = ' + np.str(np.abs(y[1,0] - y[0,0])))
    print('Ywidth_WK > ' + str(y.max() - y.min()))
    print('DEP_WK = ' + str(h.max()))
    print('Check sponge layer widths')
    print('===================================================================')
    print(' ')
    
    # End of function
    

# ==================================================================
# Write 1D bathymetry file 
# ==================================================================  
def write_bathy_1d(x,h,path,ncsave=True):
    '''
    
    Parameters:    
    ----------
    x           : 1D array of x coordinates
    h           : Bathymetry
    path        : Full path where the output will be saved
    ncsave      : Save bathy as NetCDF file
    
    Output:
    -------
    depth.txt   : Text file with the depth information for Funwave input.
    depth.nc    : (Optional) NetCDF4 bathymetry file. 
    
    Notes:
    ------
    Variables are assumed to be on a regularly spaced grid.
    
    '''

    # Output the text file -----------------------------------------------------        
    fid = open(path + 'depth.txt','w')
    for aa in range(len(h)):
        fid.write('%12.3f' % h[aa])
    fid.write('\n')
    fid.close()

    if ncsave:
    
        # Global attributes  
        nc = netCDF4.Dataset(path + 'depth.nc', 'w', format='NETCDF4')
        nc.Description = 'Funwave Bathymetry'
        nc.Author = getpass.getuser()
        nc.Created = time.ctime()
        nc.Owner = 'Nearshore Modeling Group (http://ozkan.oce.orst.edu/nmg)'
        nc.Software = 'Created with Python ' + sys.version
        nc.NetCDF_Lib = str(netCDF4.getlibversion())
        nc.Script = os.path.realpath(__file__)
     
        # Create dimensions
        xi_rho = len(h)
        nc.createDimension('xi_rho', xi_rho)
    
        # Write coordinates and depth to netcdf file
        create_nc_var(nc, 'x_rho',('xi_rho'), 
                     'meter','x-locations of RHO-points')
        nc.variables['x_rho'][:] = x
        create_nc_var(nc,'h',('xi_rho'), 
                     'meter','bathymetry at RHO-points')
        nc.variables['h'][:] = h      
                
        # Close NetCDF file
        nc.close()

    else: 
    
        print("NetCDF file not requested")
        
        

    #===========================================================================
    # Print input file options
    #===========================================================================
    print(' ')
    print('===================================================================')
    print('In your FUNWAVE input file:')
    print('Mglob = ' + np.str(len(x)))
    print('Nglob = 1')
    print('DX = ' + np.str(np.abs(x[1] - x[0])))
    print('DY = anything larger than DX works')
    print('DEP_WK = ' + str(h.max()))
    print('Check sponge layer widths')
    print('===================================================================')
    print(' ')
    
    # End of function
    
#===============================================================================
# Sponge layer and source function information 
#===============================================================================

def source_sponge_info(Tp,h,delta=0.5):
    '''
    Gives guidance for sponge layer setup and wavemaker space
    
    PARAMETERS:
    -----------
    Tp    : Wave period [s]. For irregular wavemakers Fuwnave chooses the peak
            period in the main propagation direction to determine the source
            width. 
    h     : Depth of wave maker [m]
    delta : Sponge layer width parameter (see input.txt). Usually takes 
            values from 0.3-0.6.

    RETURNS:
    --------
    Prints information on the screen
    
    '''
    
    # Compute wave length
    wl = gwaves.wave_length(Tp,h)
    W = delta * (wl/2) 
    
    # Sponge layer widths
    slw = wl*3.0
    
    print(' ')
    print('-------------------------------------------------')
    print('Wavemaker source width      = ' + np.str(W) + 'm')
    print('Minimum sponge layer width  = ' + np.str(slw) + 'm')
    print(' ')
    print('Example 1D model west coast configuration:')
    print('Sponge_west_width           = ' + np.str(slw))
    print('Xc_WK or XWAVEMAKER         = ' + np.str(slw + W/2))
    print('-------------------------------------------------')
    print(' ')



