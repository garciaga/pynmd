# -*- coding: utf-8 -*-
"""
Tools to manage xbeach input

External dependencies:
  netCDF4, getpass, time, os, sys, numpy

Internal dependencies:
  swan.pre

"""

from __future__ import division,print_function
__author__ = "Gabriel Garcia Medina"

# Import modules
import netCDF4
import getpass
import time
import os
import sys
import numpy as np

# Internal dependencies
from pynmd.models.swan import pre as gswan

#===============================================================================
# Write Grid
#===============================================================================
def write_bathy(x,z,outfld,y=None,ncsave=True):
    """
    
    Parameters:
    -----------
    x      : Vector or gridded x coordinates
    y      : Vector or gridded y coordinates (optional)
    outfld : Full path to where the files will be written
    z      : Vector or gridded depth
    ncsave : Save as netcdf file as well
    
    Output:
    -------
    Writes two or three files depending on the configuration
      x.grd, y.grd, z.dep
    
    Notes:
    ------
    Right now has been tested in 1d mode only
    
    """
    

    # Output the text file -----------------------------------------------------
    
    if y is not None:
        
        print("Writing 2D grids")
        fidx = open(outfld + 'x.grd','w')
        fidy = open(outfld + 'y.grd','w')
        fidz = open(outfld + 'z.dep','w')
                
        for aa in range(z.shape[0]):
            for bb in range(z.shape[1]):
                fidx.write('%16.4f' % x[aa,bb])
                fidy.write('%16.4f' % y[aa,bb])
                fidz.write('%16.4f' % z[aa,bb])
        
            fidx.write('\n')
            fidy.write('\n')
            fidz.write('\n')
        
        fidx.close()
        fidy.close()
        fidz.close()
        
    else:
        
        print("Writing 1D grids")
        
        # Water depth        
        fid = open(outfld + 'z.dep','w')
        for aa in range(len(z)):
            fid.write('%16.4f' % z[aa])
        fid.write('\n')
        fid.close()
        
        # Coordinates
        fid = open(outfld + 'x.grd','w')
        for aa in range(len(x)):
            fid.write('%16.4f' % x[aa])
        fid.write('\n')
        fid.close()        

    # NetCDF File --------------------------------------------------------------
    if ncsave:
    
        # Global attributes  
        nc = netCDF4.Dataset(outfld + 'depth.nc', 'w', format='NETCDF4')
        nc.Description = 'Xbeach Bathymetry'
        nc.Author = getpass.getuser()
        nc.Created = time.ctime()
        nc.Owner = 'Nearshore Modeling Group (http://ozkan.oce.orst.edu/nmg)'
        nc.Software = 'Created with Python ' + sys.version
        nc.NetCDF_Lib = str(netCDF4.getlibversion())
        nc.Script = os.path.realpath(__file__)
     
        # Create dimensions
        if y is not None:
            tmpDims = ('eta_rho','xi_rho')
            xi_rho  = z.shape[1]
            eta_rho = z.shape[0]
            nc.createDimension('eta_rho',eta_rho)
        else:
            tmpDims = ('xi_rho')
            xi_rho = z.shape[-1]
        
        nc.createDimension('xi_rho', xi_rho)
    
        # Write coordinates and depth to netcdf file
        nc.createVariable('x_rho','f8',tmpDims)
        nc.variables['x_rho'].units = 'meter'
        nc.variables['x_rho'].long_name = 'x-locations of RHO-points'
        nc.variables['x_rho'][:] = x
        
        nc.createVariable('h','f8',tmpDims)
        nc.variables['h'].units = 'meter'
        nc.variables['h'].long_name = 'bathymetry at RHO-points'
        nc.variables['h'][:] = z
        
        if y is not None:
            nc.createVariable('y_rho','f8',tmpDims)
            nc.variables['y_rho'].units = 'meter'
            nc.variables['y_rho'].long_name = 'y-locations of RHO-points'
            nc.variables['y_rho'][:] = y
                        
        # Close NetCDF file
        nc.close()
        

    #===========================================================================
    # Print input file options
    #===========================================================================
    print(' ')
    print('===================================================================')
    print('In your params.txt:')
    print('  nx = ' + np.str(xi_rho - 1))
    if y is not None:
        print('  ny = ' + np.str(eta_rho - 1))
    else:
        print('  ny = 0')
    print('===================================================================')
    
    # End of function ----------------------------------------------------------
    

#===============================================================================
# Write boundary spectrum
#===============================================================================
def write_boundary_spec(freq,spec,locations,outfile,waveTime=None,
                        sdir=None):
    """
    Calls pynmd.models.swan.pre.write_boundary_spec
    """
    
    gswan.write_boundary_spec(freq,spec,locations,outfile,waveTime=None,
                              sdir=None)
    