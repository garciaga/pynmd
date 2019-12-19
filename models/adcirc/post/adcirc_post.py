# -*- coding: utf-8 -*-
"""
A series of tools to process and analyze adcirc output

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

# ==============================================================================
# Read Fort 63 ASCII files and save as nc file
# ==============================================================================
def fort63_to_nc(fort63,varname='zeta',longname='water surface elevation above geoid',
                varunits='m',**kwargs):
    """ 
    Script to read fort.63-type files and store in a netcdf4 file

    PARAMETERS:
    -----------
    fort63: Path to fort63 file

    RETURNS:
    --------
    Netcdf containing
    time  : seconds since beginning of run 
    var   : temporal variable (called 'varname') recorded at the nodes of
            an unstructured grid. Size: [time,nodes]

    """
    
    fobj = open(fort63,'r')
    
    # Create the file and add global attributes
    if 'savename' in kwargs:
        ncfile = kwargs['savename']
    else:
        ncfile = fort63 + '.nc'
    nc = netCDF4.Dataset(ncfile, 'w', format='NETCDF4')
    
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
    
    # Record number of time steps and nodes
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
    nc.variables[varname].long_name = longname
    nc.variables[varname].units = varunits
        
    for tt in range(ntsteps):
        # Store variables
        nc.variables['time'][tt] = np.float64(fobj.readline().split()[0])

        for aa in range(nodes):
            nc.variables[varname][tt,aa] = np.float64(fobj.readline().split()[1])

    # All done here
    fobj.close()
    nc.close()
    