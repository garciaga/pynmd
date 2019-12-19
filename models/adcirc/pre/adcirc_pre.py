# -*- coding: utf-8 -*-
"""
A series of tools to process and prepare adcirc input

Authors:
-------
Fadia Ticona Rollano
    PNNL Marine Sciences Laboratory

Log of edits:
-------------
December 2019 - Created module
    Fadia Ticona Rollano

"""

from __future__ import division,print_function

#import glob,os
import sys,time
#sys.path.append(os.getcwd())
import getpass
import numpy as np
import netCDF4

# Custom paths
import pynmd.models.tools.unstructured as ustr

# ==============================================================================
# Read Fort 63 ASCII Files
# ==============================================================================
def fort14_to_nc(fort14,**kwargs):
    """ 
    Script to read fort.63-type files and store in a netcdf4 file

    PARAMETERS:
    -----------
    fort14: Path to fort14 file

    RETURNS:
    --------
    Netcdf containing
    var   : variable recorded at the nodes of an unstructured grid. Defaults as 
            'zeta'. Size: [time,nodes]

    """
    grid = ustr.read_fort14(fort14)
    
    # Create the file and add global attributes
    if 'savename' in kwargs:
        ncfile = kwargs['savename']
    else:
        ncfile = fort14 + '.nc'
    nc = netCDF4.Dataset(ncfile, 'w', format='NETCDF4')
    
    # Global attributes 
    nc.Author = getpass.getuser()
    nc.Created = time.ctime()  
    nc.Software = 'Created with Python ' + sys.version
    nc.NetCDF_Lib = str(netCDF4.getlibversion())
    
    # Get record lengths
    node = grid['x'].size
    nele = grid['triang'].size
    nvertex = 3
#    nope = grid['nope']
    neta = grid['neta']    
#    nbou = grid['nbou']
    nvel = grid['nvel']
    
    # Create dimensions
    nc.createDimension('node',node)       # Number of nodes
    nc.createDimension('nele',nele)       # Number of elements
    nc.createDimension('nvertex',nvertex) # Number of vertex
    nc.createDimension('neta',neta)       # Number of elevation boundary nodes
    nc.createDimension('nvel',nvel)       # Number of normal flow boundary nodes
    
    # Create variables
    nc.createVariable('x','f8','node')
    nc.variables['x'].long_name = 'longitude'
    nc.variables['x'].units = 'degrees_east'
    nc.variables['x'].positive = 'east'
    nc.variables['x'][:] = np.float64(grid['x'])
    
    nc.createVariable('y','f8','node')
    nc.variables['y'].long_name = 'latitude'
    nc.variables['y'].units = 'degrees_north'
    nc.variables['y'].positive = 'north'
    nc.variables['y'][:] = np.float64(grid['y'])
    
    nc.createVariable('element','i4',('nele','nvertex'))
    nc.variables['element'].long_name = 'element'
    nc.variables['element'].units = 'nondimensional'
    nc.variables['element'][:,:] = np.int32(grid['triang'])
    
    nc.createVariable('nbdv','i4','neta')
    nc.variables['nbdv'].long_name = 'node numbers on each elevation specified boundary segment'
    nc.variables['nbdv'].units = 'nondimensional'
    nc.variables['nbdv'][:] = np.int32(np.hstack(grid['nbdv']))
    
    nc.createVariable('nbvv','i4','nvel')
    nc.variables['nbvv'].long_name = 'node numbers on normal flow boundary segment'
    nc.variables['nbvv'].units = 'nondimensional'
    nc.variables['nbvv'][:] = np.int32(np.hstack(grid['nbvv']))
    
    nc.createVariable('depth','f8','node')
    nc.variables['depth'].long_name = 'distance below geoid'
    nc.variables['depth'].units = 'm'
    nc.variables['depth'][:] = np.float64(grid['z'])
    
    
    
    # All done here
    nc.close()