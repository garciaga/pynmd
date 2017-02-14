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

#===============================================================================
# Backwards compatibility stuff
#===============================================================================
def makeBathy1DPlanar(m,dx,hmax,hmin,flat,outFld):
    x,h = makeBathyPlanar(m,dx,hmax,hmin,flat,outFld,dy=None,ly=None)
    return x,h

def write_bathy_1d(x,h,path,ncsave=True,dt=None):
    write_bathy(x,h,path,y=None,ncsave=ncsave,dt=dt)

# ==============================================================================
# Write 1D bathymetry file 
# ==============================================================================
def write_bathy(x,h,path,y=None,ncsave=True,dt=None):
    '''
    
    Parameters:    
    ----------
    x           : array of x coordinates
    y           : array of x coordinates (optional)
    h           : Bathymetry
    path        : Full path where the output will be saved
    ncsave      : Save bathy as NetCDF file
    dt          : (Optional) informs about the stability for a chosen dt
    
    Output:
    -------
    depth.txt   : Text file with the depth information for Funwave input.
    depth.nc    : (Optional) NetCDF4 bathymetry file. 
    
    Notes:
    ------
    1. Variables are assumed to be on a regularly spaced grid.
    2. If y is passed then 2D bathymetry is assumed. This means that 
       x,y,h have to be 2D arrays. Otherwise the code will fail.
    
    '''

    # Output the text file -----------------------------------------------------        
    fid = open(path + 'depth.txt','w')
    if y is None:
        for aa in range(len(h)):
            fid.write('%12.3f\n' % h[aa])
        fid.close()
    else:
        for aa in range(h.shape[1]):
            for bb in range(h.shape[0]):
                fid.write('%12.3f' % h[bb,aa])
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
        if y is None:
            xi_rho = len(h)
            nc.createDimension('xi_rho', xi_rho)
            varShape = ('xi_rho')
        else:
            eta_rho,xi_rho = np.shape(h)
            nc.createDimension('xi_rho', xi_rho)
            nc.createDimension('eta_rho',eta_rho)
            varShape = ('eta_rho','xi_rho')
            
            # Write y variable
            create_nc_var(nc, 'y_rho',varShape, 
                          'meter','y-locations of RHO-points')
            nc.variables['y_rho'][:] = y
                    
        # Write coordinates and depth to netcdf file
        create_nc_var(nc, 'x_rho',varShape, 
                     'meter','x-locations of RHO-points')
        nc.variables['x_rho'][:] = x
        
        create_nc_var(nc,'h',varShape, 
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
    if y is None:
        dx = np.abs(x[1] - x[0])
        print('dimension ' + np.str(len(x)+1) + ' 6 ' +
              np.str(dx) + ' 1')
    else:
        dx = np.abs(x[0,1] - x[0,0]) # For later use
        dy = np.abs(y[1,0] - y[0,0])
        print('dimension ' + np.str(x.shape[1]+1) + ' ' + 
              np.str(x.shape[0]) + ' ' +
              np.str(dx) + ' ' + 
              np.str(dy))
    print('===================================================================')
    print(' ')
    
    # Stability options    
    stabilityCriteria(dx,h.max(),dt=dt,verbose=True)
    
    # End of function


#===============================================================================
# Model stability criteria    
#===============================================================================
def stabilityCriteria(dx,hmax,dy=None,dt=None,verbose=False):
    """
    Give the user some ideas on the stability parameters
    """
    
    # CFL stability criterion (wave only) --------------------------------------
    print('Wave CFL Stability')
    cmax = (9.81*hmax)**0.5
    if not dt:
        dt = 0.2*dx/cmax # x direction
        if dy:
            dt = np.min([dt,0.2*dy/cmax])
        print('  Wave CFL stability requires: dt < ' + '{:6.4f}'.format(dt))
    else:
        cfl_wave = cmax * dt / dx # x direction
        if dy:
            cfl_wave = np.max([cfl_wave,cmax*dt/dy])
            
        if cfl_wave > 0.2:
            print('  Warning CFL Wave = ' + '{:6.4f}'.format(cfl_wave) + 
                  ' > 0.2')
        else:
            print('  CFL Wave = ' + '{:6.4f}'.format(cfl_wave))
            print('  dt is ok')
    
    if verbose:
        print('    cmax * dt / min(dx,dy) < 0.2')
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
    if dy:
        gamma_bi = np.max([gamma_bi,0.008*(dy**4)/dt])
    print('  Maximum biharmonic friction layer damping (gamma_bi) = ' + 
          '{:6.4f}'.format(gamma_bi))
    
    if verbose:
        print('    S_bi = gamma_bi * dt / min(dx**4,dy**4) < 0.008')
    
    # Breaking stability -------------------------------------------------------
    if verbose:
        print(' ')
        print('Breaking Stability')
        print('    S_br = gamma_br * dt / (dx**2) < 0.1')
        
#===============================================================================
# Create input file    
#===============================================================================
def makeInput(inp,outfld):
    """
    Create input file to be used for funwaveC
    
    PARAMETERS:
    -----------
    inp     : Dictionary containing input file parameters
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
    dx = np.str(inp['dx'])
    nx = np.str(inp['nx'])
    if 'ny' in inp.keys():
        ny = np.str(inp['ny'])
        dy = np.str(inp['dy'])
    else:
        ny = '6'
        dy = np.str(np.ceil(np.float64(dx)))

    fid.write('dimension ' + nx + ' ' + ny + ' ' + dx + ' ' + dy + '\n') 
    
    # Bottom stress
    fid.write('bottomstress const ' + np.str(inp['bottomstress']) + '\n')
    
    # Lateral friction
    fid.write('mixing ' + inp['mixing'][0] + ' ' + 
              np.str(inp['mixing'][1]) + '\n')
    
    # Bathymetry
    if np.double(ny) > 6.0:
        # 2D bathymetry assumed
        fid.write('bathymetry file2d depth.txt\n')
    else:
        # 1D bathymetry
        fid.write('bathymetry file1d depth.txt\n')
    
    # Tides
    fid.write('tide off\n')
    
    # Wavemaker
    if inp['eta_source'][0] == 'off':
        fid.write('eta_source off\n')
    elif inp['eta_source'][1] == 'random2nb':
        fid.write('eta_source on ' + inp['eta_source'][1])
        for aa in range(2,len(inp['eta_source'])):
            if aa == 5 or aa == 8 or aa == 10:
                fid.write(' ' + np.str(np.int(inp['eta_source'][aa])))
            else:
                fid.write(' ' + np.str(inp['eta_source'][aa]))
        fid.write('\n')
    elif inp['eta_source'][1] == 'random2filea':
        fid.write('eta_source on ' + inp['eta_source'][1])
        fid.write(' ' + inp['eta_source'][2]) # File name
        for aa in range(3,len(inp['eta_source'])):
            fid.write(' ' + np.str(np.int(inp['eta_source'][aa])))
        fid.write('\n')
    else:
        print('Accepted eta_source parameters:')
        print('  [random2nb,random2filea]')
        print('  Turning wavemaker off')
        fid.write('eta_source off\n')
        
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
def makeBathyPlanar(m,dx,hmax,hmin,flat,outFld,dy=None,ly=None):
    """
    Make a planar beach bathymetry
    
    PARAMETERS:
    -----------
    m    : Beach slope
    dx   : Resolution [m]
    hmax : Maximum depth [m]
    hmin : Minimum depth [m]
    flat : Length of flat part [m]
    
    2D PARAMETRS (optional):
    ------------------------
    dy   : Grid resolution in y direction [m]
    ly   : Grid alongshore length [m]
        
    RETURNS:
    --------
    x     : cross-shore coordinates [m]
    y     : alongshore coordinates [m] (if dy and ly are provided)
    h     : bathymetry [m]
    
    NOTES:
    ------
    - z is positive downwards
    
    """

    # Generate bathymetry
    xmax = (hmax - hmin)/m + flat
    x = np.arange(1,xmax+dx,dx)
    h = x*m + hmin
    h[h>hmax] = hmax
    h = np.flipud(h)

    if dy is None:
        write_bathy(x,h,outFld,y=None)        
        return x,h
    else:
        y = np.arange(1.0,ly+dy,dy)
        x,y = np.meshgrid(x,y)
        h = np.repeat(np.expand_dims(h,axis=0),y.shape[0],axis=0)
        write_bathy(x,h,outFld,y=y)
        return x,y,h

#===============================================================================
# Write spectra for input
#===============================================================================
def writeSpec1D(freq,spec,theta,spread,outFld):
    """
    Write a 1D spectra for funwaveC
    
    PARAMETERS:
    -----------
    freq   : frequency array [Hz]
    spec   : spectra array [m**2/Hz]
    theta  : Wave direction per frequency [deg]
    spread : Directional spread per frequency [deg]
    outFld : Output folder
    
    RETURNS:
    --------
    spec.txt with columns (freq fourierAmplitude theta spread)
    
    NOTES:
    ------
    1. The spectrum should not have a zeroth frequency. FunwaveC will not run.
    2. Frequency should be equally spaced
    3. Numpy arrays are expected
    """
    
    # Scale the spectrum
    # FunwaveC needs fourier amplitudes not spectra
    df = freq[2] - freq[1]
    fouramp = (df * spec * 2.0)**0.5 
    
    # Write the spectra file
    fid = open(outFld + 'spec.txt','w')
    for aa in range(fouramp.shape[0]):
        fid.write('{:12.8f}'.format(freq[aa]) + ' ' +
                  '{:12.8f}'.format(fouramp[aa]) + ' ' +
                  '{:12.3f}'.format(theta[aa]) + ' ' +
                  '{:12.3f}'.format(spread[aa]) + '\n')        
    fid.close()
    
    