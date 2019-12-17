# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 16:13:02 2019

@author: tico327
"""

from __future__ import division,print_function

#import glob,os
import sys,time
#sys.path.append(os.getcwd())
import getpass
import numpy as np
import netCDF4

# ==============================================================================
# Read Fort 63 ASCII Files
# ==============================================================================
def read_fort63(fort63,varname='zeta',
                longvarname='water surface elevation above geoid',varunits='m'):
    """ 
    Script to read fort14 files

    PARAMETERS:
    -----------
    fort14: Path to fort14 file

    RETURNS:
    --------
    Dictionary containing
    x     : x coordinates of nodes (ordered)
    y     : y coordinates of nodes (ordered)
    z     : Elevation at nodes
    triang: Triangulation 
    nbdv  : Indices for open boundary nodes
    neta  : number of open boundary nodes
    nvel  : Indices for land boundaries
    nvel  : number of land boundary nodes

    NOTES:
    ------
    1. Zero counting convention used.
    2. Nodes are returned in order so that you can combine with triangulation

    """
    
    fobj = open(fort63,'r')
    
    # Create the file and add global attributes
    nc = netCDF4.Dataset(fort63 + '.nc', 'w', format='NETCDF4')
    
    # Global attributes 
    nc.Author = getpass.getuser()
    nc.Created = time.ctime()
    tmpline = fobj.readline()
    nc.description = tmpline[:32]
    nc.rundes = tmpline[:32]
    nc.runid = tmpline[32:56]
    nc.model = 'ADCIRC'    
    nc.Software = 'Created with Python ' + sys.version
    nc.NetCDF_Lib = str(netCDF4.getlibversion())
    
    tmpline = fobj.readline().split()
    ntsteps = np.int(tmpline[0])
    nodes = np.int(tmpline[1])
    
    # Create dimensions
    nc.createDimension('node',nodes)       # Number of nodes
    nc.createDimension('time',0)           # The unlimited dimension
    
    # Create time vector
    nc.createVariable('time','f8',('time'))
    nc.variables['time'].long_name = 'model time'
    nc.variables['time'].standard_name = 'time'
    nc.variables['time'].units = 'seconds since beginning of run'
    
    # Create the rest of the variables
    nc.createVariable(varname,'f8',('time','node'))
    nc.variables[varname].long_name = longvarname
    nc.variables[varname].units = varunits
        
    for tt in range(ntsteps):
        # Store variables
        nc.variables['time'][tt] = np.float64(fobj.readline().split()[0])

        for aa in range(nodes):
            nc.variables[varname][tt,aa] = np.float64(fobj.readline().split()[1])

    # All done here
    fobj.close()
