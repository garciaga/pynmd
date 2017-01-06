"""
Tools to manage funwaveC input

Authors:
--------
Gabriel Garcia Medina (ggarcia@coas.oregonstate.edu)

# Last edit
14 December 2016 - Gabriel Garcia Medina

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
# Write 1D bathymetry file 
# ==================================================================  
def write_bathy_1d(x,h,path,ncsave=True,dt=None):
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
        fid.write('%12.3f\n' % h[aa])
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
    print('In your funwaveC init file:')
    print('dimension ' + np.str(len(x)+1) + ' 6 ' +
          np.str(np.abs(x[1] - x[0])) + ' 1')
    print('===================================================================')
    print(' ')
    
    # Stability options
    dx = x[1] - x[0]
    stabilityCriteria(dx,h.max(),dt=dt,verbose=True)
    
    # End of function

#===============================================================================
# Model stability criteria    
#===============================================================================
def stabilityCriteria(dx,hmax,dt=None,verbose=False):
    """
    Give the user some ideas on the stability parameters
    """
    
    # CFL stability criterion (wave only) --------------------------------------
    print('Wave CFL Stability')
    cmax = (9.81*hmax)**0.5
    if not dt:
        dt = 0.2*dx/cmax
        print('  Wave CFL stability requires: dt < ' + '{:6.4f}'.format(dt))
    else:
        cfl_wave = cmax * dt / dx
        if cfl_wave > 0.2:
            print('  Warning CFL Wave = ' + '{:6.4f}'.format(cfl_wave) + 
                  ' > 0.2')
        else:
            print('  CFL Wave = ' + '{:6.4f}'.format(cfl_wave))
            print('  dt is ok')
    
    if verbose:
        print('    cmax * dt / dx < 0.2')
        print('    cmax = (9.81 * hmax)**0.5 = ' + '{:8.4f}'.format(cmax))
        
    # Sponge layer -------------------------------------------------------------
    print(' ')
    print('Sponge Layer')    
    k = 0.25 / dt
    print('  Maximum sponge layer damping (k) = ' + '{:6.4f}'.format(k))
    
    if verbose:
        print('    S_sp = k * dt < 0.25')
    
    # Biharmonic friction ------------------------------------------------------
    print(' ')
    print('Biharmonic Friction')
    gamma_bi = 0.008 * (dx**4) / dt
    print('  Maximum biharmonic friction layer damping (gamma_bi) = ' + 
          '{:6.4f}'.format(gamma_bi))
    
    if verbose:
        print('    S_bi = gamma_bi * dt / (dx**4) < 0.008')
    
    # Breaking stability -------------------------------------------------------
    if verbose:
        print(' ')
        print('Breaking Stability')
        print('    S_br = gamma_br * dt / (dx**2) < 0.1')
        
#===============================================================================
# Create input file    
#===============================================================================
def makeInput(inp,x,outfld):
    """
    Create input file to be used for funwaveC
    
    PARAMETERS:
    -----------
    inp     : Dictionary containing input file parameters
    x       : x axis of the bathymetry 
    outfld  : path to write the input file (i.e. outfld/input.init)
    
    RETURNS:
    --------
    funwaveC input file
    
    INPUT DICTIONARY:
    -----------------
    line1: Needs dynamics keyword
    
    
    NOTES:
    ------
    - Many limitations right now as I am getting started with the model
    - Bathymetry file assumed to be called depth.txt
    """
    
    # Create input file        
    fid = open(outfld + 'input.init','w')
    
    # Funwave dynamics
    fid.write('funwaveC dynamics ' + inp['dynamics'] + '\n')
    
    # Grid size and spacing
    if len(x.shape) == 1:
        dx = x[2] - x[1]
        # 1D model 
        fid.write('dimension ' + np.str(x.shape[0]+1) + ' 6 ' + 
                  np.str(dx) + ' ' + 
                  np.str(np.ceil(dx)) + '\n')
    else:
        print('Need to implement 2D model stuff')
    
    # Bottom stress
    fid.write('bottomstress const ' + np.str(inp['bottomstress']) + '\n')
    
    # Lateral friction
    fid.write('mixing ' + inp['mixing'][0] + ' ' + 
              np.str(inp['mixing'][1]) + '\n')
    
    # Bathymetry
    if len(x.shape) == 1:
        # 1D bathymetry
        fid.write('bathymetry file1d depth.txt\n')
    else:
        # 2D bathymetry assumed
        fid.write('bathymetry file2d depth.txt\n')
    
    # Tides
    fid.write('tide off\n')
    
    # Wavemaker
    if inp['eta_source'][0] == 'off':
        fid.write('eta_source off\n')
    else:
        fid.write('eta_source on ' + inp['eta_source'][1])
        for aa in range(2,len(inp['eta_source'])):
            if aa == 5 or aa == 8 or aa == 10:
                fid.write(' ' + np.str(np.int(inp['eta_source'][aa])))
            else:
                fid.write(' ' + np.str(inp['eta_source'][aa]))
        fid.write('\n')
        
    # Wave breaking
    if inp['breaking'] == 'off':
        fid.write('breaking off\n')
    else:
        fid.write('breaking on ' + inp['breaking'] + '\n')
    
    # Sponge layer
    if inp['sponge'] == 'off':
        fid.write('sponge off\n')
    else:
        fid.write('sponge on')
        for aa in range(len(inp['sponge'])):
            fid.write(' ' + np.str(inp['sponge'][aa]))
        fid.write('\n')
    
    # Forcing, initial conditions, tracers, etc not supported
    fid.write('forcing none 0.0 none 0.0\n')
    fid.write('initial_condition none 0.0 none 0.0 none 0.0\n')
    fid.write('tracer off\n')
    fid.write('stracer off\n')
    fid.write('rtracer off\n')
    fid.write('floats off\n')
    
    # Timing stuff only seconds supported
    fid.write('timing')
    for aa in range(4):
        # Hardcoded so it returns error if input is incomplete
        fid.write(' ' + np.str(inp['timing'][aa]) + ' (sec)')
    fid.write('\n')
    
    # Output stuff (generalize later)
    fid.write('level1 cross N 3 file ascii eta.dat\n')
    fid.write('level1 cross U 3 file ascii u.dat\n')
            
    # Close input file
    fid.close()

#===============================================================================
# Make a 1D bathymetry
#===============================================================================
def makeBathy1DPlanar(m,dx,hmax,hmin,flat,outFld):
    """
    Make a planar beach bathymetry
    
    PARAMETERS:
    -----------
    m    : Beach slope
    dx   : Resolution [m]
    hmax : Maximum depth [m]
    hmin : Minimum depth [m]
    flat : Length of flat part [m]
    
    RETURNS:
    --------
    x     : coordinates [m]
    h     : bathymetry [m]
    
    NOTES:
    ------
    - z is positive downwards
    - Should expand to 2D
    
    """

    # Generate bathymetry
    xmax = (hmax - hmin)/m + flat
    x = np.arange(1,xmax+dx,dx)
    h = x*m + hmin
    h[h>hmax] = hmax
    h = np.flipud(h)

    write_bathy_1d(x,h,outFld)

    return x,h
