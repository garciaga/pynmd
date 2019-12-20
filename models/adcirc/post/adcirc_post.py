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

import os#,glob
import sys,time
#sys.path.append(os.getcwd())
import getpass
import numpy as np
import netCDF4

# Custom paths
import pynmd.models.adcirc.pre as adcpre

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
    time     : seconds since beginning of run 
    variable : temporal variable (called 'varname') recorded at the nodes of
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


# ==============================================================================
# Read Fort 64 ASCII files and save as nc file
# ==============================================================================
def fort64_to_nc(fort64,varname_xy=['u-vel','v-vel'],
                 longname='water column vertically averaged',
                 varunits='m s-1',**kwargs):
    """ 
    Script to read fort.64-type files and store in a netcdf4 file

    PARAMETERS:
    -----------
    fort64: Path to fort64 file

    RETURNS:
    --------
    Netcdf containing
    time         : seconds since beginning of run 
    x,y variable : temporal variables (name provided in 'xy_varname') recorded
                   at the nodes of an unstructured grid. Size: [time,nodes]
    """
    
    fobj = open(fort64,'r')
    
    # Create the file and add global attributes
    if 'savename' in kwargs:
        ncfile = kwargs['savename']
    else:
        ncfile = fort64 + '.nc'
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
    nc.createVariable(varname_xy[0],'f8',('time','node'))
    nc.variables[varname_xy[0]].long_name = longname + ' e/w velocity'
    nc.variables[varname_xy[0]].units = varunits
    
    nc.createVariable(varname_xy[1],'f8',('time','node'))
    nc.variables[varname_xy[1]].long_name = longname + ' n/s velocity'
    nc.variables[varname_xy[1]].units = varunits
        
    for tt in range(ntsteps):
        # Store variables
        nc.variables['time'][tt] = np.float64(fobj.readline().split()[0])

        for aa in range(nodes):
            tmpline = fobj.readline().split()
            nc.variables[varname_xy[0]][tt,aa] = np.float64(tmpline[1])
            nc.variables[varname_xy[1]][tt,aa] = np.float64(tmpline[2])

    # All done here
    fobj.close()
    nc.close()


# ==============================================================================
# Read all ASCII files and save as nc file
# ==============================================================================
def all_ascii2nc(runFld,ncFld):
    """ 
    Reads (most*) ADCIRC ASCII files and stores them in netcdf4 files

    Input:
    ------
    runFld: folder containing all ascii files
    ncFld:  folder where to save all nc files
    
    Notes:
    ------        
    *Still working on converting some adcirc files
    """
    # Unstructured grid + bathy ------------------------------------------------
    if os.path.exists(runFld+'fort.14'):
        if not os.path.exists(ncFld+'fort.14.nc'):
            print('Creating fort.14.nc')
            adcpre.fort14_to_nc(runFld + 'fort.14',savename=ncFld+'fort.14.nc')
    
    # fort.63-type files (scalars) ---------------------------------------------
    # Water surface elevation
    if os.path.exists(runFld+'fort.63'):
        if not os.path.exists(ncFld+'fort.63.nc'):
            print('Creating fort.63.nc')
            fort63_to_nc(runFld + 'fort.63',savename=ncFld+'fort.63.nc')
    # Atmospherec pressure
    if os.path.exists(runFld+'fort.73'):
        if not os.path.exists(ncFld+'fort.73.nc'):
            print('Creating fort.73.nc')
            fort63_to_nc(runFld + 'fort.73',varname='pressure',
                         longname='air pressure at sea level',
                         varunits='meters of waver',savename=ncFld+'fort.73.nc')
    # Significant wave height
    if os.path.exists(runFld+'swan_HS.63'):
        if not os.path.exists(ncFld+'swan_HS.63.nc'):
            print('Creating swan_HS.63.nc')
            fort63_to_nc(runFld + 'swan_HS.63',varname='swan_HS',
                         longname='significant wave height',
                         varunits='m',savename=ncFld+'swan_HS.63.nc')
    # Mean wave direction
    if os.path.exists(runFld+'swan_DIR.63'):
        if not os.path.exists(ncFld+'swan_DIR.63.nc'):
            print('Creating swan_DIR.63.nc')
            fort63_to_nc(runFld + 'swan_DIR.63',varname='swan_DIR',
                         longname='mean wave direction',
                         varunits='degrees',savename=ncFld+'swan_DIR.63.nc')
    # Mean absolute wave period
    if os.path.exists(runFld+'swan_TMM10.63'):
        if not os.path.exists(ncFld+'swan_TMM10.63.nc'):
            print('Creating swan_TMM10.63.nc')
            fort63_to_nc(runFld + 'swan_TMM10.63',varname='swan_TMM10',
                         longname='mean absolute wave period',
                         varunits='s',savename=ncFld+'swan_TMM10.63.nc')
    # Smoothed peak period
    if os.path.exists(runFld+'swan_TPS.63'):
        if not os.path.exists(ncFld+'swan_TPS.63.nc'):
            print('Creating swan_TPS.63.nc')
            fort63_to_nc(runFld + 'swan_TPS.63',varname='swan_TPS',
                         longname='smoothed peak period',
                         varunits='s',savename=ncFld+'swan_TPS.63.nc')
    
    # fort.64-type files (vectors) ---------------------------------------------
    # Depth average velocity 
    if os.path.exists(runFld+'fort.64'):
        if not os.path.exists(ncFld+'fort.64.nc'):
            print('Creating fort.64.nc')
            fort64_to_nc(runFld + 'fort.64',savename=ncFld+'fort.64.nc')
    # Wind stress or velocity 
    if os.path.exists(runFld+'fort.74'):
        if not os.path.exists(ncFld+'fort.74.nc'):
            print('Creating fort.74.nc')
            fort64_to_nc(runFld + 'fort.74',varname_xy=['windx','windy'],
                         longname='wind',savename=ncFld+'fort.74.nc')
    