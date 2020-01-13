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
# Read fort.61-type ASCII files and save as nc file
# ==============================================================================
def fort61_to_nc(fort61,staname,x,y,varname='zeta',
                 longname='water surface elevation above geoid',
                 varunits='m',ncdate='0000-00-00 00:00:00 UTC',**kwargs):
    """ 
    Script to read fort.61-type (station scalar) files and store in a netcdf4 file

    PARAMETERS:
    -----------
    fort61 : Path to fort61-type file
    x,y    : Station coordinates
    ncdate : cold start date/time in CF standard: yyyy-MM-dd hh:mm:ss tz

    RETURNS:
    --------
    Netcdf containing
    time     : seconds since beginning of run 
    variable : temporal variable (called 'varname') recorded at stations. 
               Size: [time,station]
    """
        
    fobj = open(fort61,'r')
    
    # Create the file and add global attributes
    if 'savename' in kwargs:
        ncfile = kwargs['savename']
    else:
        ncfile = fort61 + '.nc'
    nc = netCDF4.Dataset(ncfile, 'w', format='NETCDF4')
    
    # Global attributes 
    nc.Author = getpass.getuser()
    nc.Created = time.ctime()
    tmpline = fobj.readline()
    nc.description = tmpline[2:34]
    nc.rundes = tmpline[2:34]
    nc.runid = tmpline[36:60]
    nc.model = 'ADCIRC'    
    nc.Software = 'Created with Python ' + sys.version
    nc.NetCDF_Lib = str(netCDF4.getlibversion())
    
    # Record number of time steps and station
    tmpline = fobj.readline().split()
    ntsteps = np.int(tmpline[0])
    sta = np.int(tmpline[1])
    
    # Create dimensions
    nc.createDimension('time',0)           # The unlimited dimension
    nc.createDimension('station',sta)      # Number of stations
    nc.createDimension('namelen',50)       # Length of station names
    
    # Create time vector
    nc.createVariable('time','f8',('time'))
    nc.variables['time'].long_name = 'model time'
    nc.variables['time'].standard_name = 'time'
    nc.variables['time'].units = 'seconds since ' + ncdate
    
    # Create and store spatial variables
    nc.createVariable('station_name','S1',('station','namelen'))
    nc.variables['station_name'].long_name = 'station name'
    nc.variables['station_name'][:] = netCDF4.stringtochar(np.array(staname,dtype='S50'))
    nc.createVariable('x','f8','station')
    nc.variables['x'].long_name = 'longitude'
    nc.variables['x'].units = 'degrees east'
    nc.variables['x'].positive = 'east'
    nc.variables['x'][:] = x
    nc.createVariable('y','f8','station')
    nc.variables['y'].long_name = 'latitude'
    nc.variables['y'].units = 'degrees north'
    nc.variables['y'].positive = 'north'
    nc.variables['y'][:] = y
    
    # Create station variable
    nc.createVariable(varname,'f8',('time','station'))
    nc.variables[varname].long_name = longname
    nc.variables[varname].units = varunits
        
    # Store time-series variables
    for tt in range(ntsteps):
        nc.variables['time'][tt] = np.float64(fobj.readline().split()[0])
        for aa in range(sta):
            nc.variables[varname][tt,aa] = np.float64(fobj.readline().split()[1])

    # All done here
    fobj.close()
    nc.close()


# ==============================================================================
# Read fort.63-type ASCII files and save as nc file
# ==============================================================================
def fort63_to_nc(fort63,varname='zeta',
                 longname='water surface elevation above geoid',
                 varunits='m',ncdate='0000-00-00 00:00:00 UTC',**kwargs):
    """ 
    Script to read fort.63-type (scalar) files and store in a netcdf4 file

    PARAMETERS:
    -----------
    fort63: Path to fort63-type file
    ncdate : cold start date/time in CF standard: yyyy-MM-dd hh:mm:ss tz

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
    nc.description = tmpline[2:34]
    nc.rundes = tmpline[2:34]
    nc.runid = tmpline[36:60]
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
    nc.variables['time'].units = 'seconds since ' + ncdate
    
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
# Read fort.64-type ASCII files and save as nc file
# ==============================================================================
def fort64_to_nc(fort64,varname_xy=['u-vel','v-vel'],
                 longname='water column vertically averaged',
                 varunits='m s-1',ncdate='0000-00-00 00:00:00 UTC',**kwargs):
    """ 
    Script to read fort.64-type (vector) files and store in a netcdf4 file

    PARAMETERS:
    -----------
    fort64: Path to fort64-type file
    ncdate : cold start date/time in CF standard: yyyy-MM-dd hh:mm:ss tz

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
    nc.description = tmpline[2:34]
    nc.rundes = tmpline[2:34]
    nc.runid = tmpline[36:60]
    nc.model = 'ADCIRC'    
    nc.Software = 'Created with Python ' + sys.version
    nc.NetCDF_Lib = str(netCDF4.getlibversion())
    
    # Record number of time steps and nodes
    tmpline = fobj.readline().split()
    ntsteps = np.int(tmpline[0])
    nodes = np.int(tmpline[1])
    
    # Create dimensions
    nc.createDimension('time',0)           # The unlimited dimension
    nc.createDimension('node',nodes)       # Number of nodes
    
    # Create time vector
    nc.createVariable('time','f8',('time'))
    nc.variables['time'].long_name = 'model time'
    nc.variables['time'].standard_name = 'time'
    nc.variables['time'].units = 'seconds since ' + ncdate
    
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
# Read Global Maximum and Minimum ASCII files and save as nc file
# ==============================================================================
def max63_to_nc(max63,varname='zeta',
                 longname='water surface elevation above geoid',
                 varunits='m',ncdate='0000-00-00 00:00:00 UTC',**kwargs):
    """ 
    Script to read max63-type (max-min) files and store in a netcdf4 file

    PARAMETERS:
    -----------
    max63: Path to max63-type file
    ncdate : cold start date/time in CF standard: yyyy-MM-dd hh:mm:ss tz

    RETURNS:
    --------
    Netcdf containing
    time         : seconds since beginning of run 
    x,y variable : temporal variables (name provided in 'xy_varname') recorded
                   at the nodes of an unstructured grid. Size: [time,nodes]
    """
    
    fobj = open(max63,'r')
    
    # Create the file and add global attributes
    if 'savename' in kwargs:
        ncfile = kwargs['savename']
    else:
        ncfile = max63 + '.nc'
    nc = netCDF4.Dataset(ncfile, 'w', format='NETCDF4')
    
    # Global attributes 
    nc.Author = getpass.getuser()
    nc.Created = time.ctime()
    tmpline = fobj.readline()
    nc.description = tmpline[2:34]
    nc.rundes = tmpline[2:34]
    nc.runid = tmpline[36:60]
    nc.model = 'ADCIRC'    
    nc.Software = 'Created with Python ' + sys.version
    nc.NetCDF_Lib = str(netCDF4.getlibversion())
    
    # Record number of time steps and nodes
    nodes = np.int(fobj.readline().split()[1])
    
    # Create dimensions
    nc.createDimension('time',0)           # The unlimited dimension
    nc.createDimension('node',nodes)       # Number of nodes
    
    # Create time vector
    nc.createVariable('time','f8',('time'))
    nc.variables['time'].long_name = 'model time'
    nc.variables['time'].standard_name = 'time'
    nc.variables['time'].units = 'seconds since ' + ncdate
    
    # Create the rest of the variables
    nc.createVariable(varname+'_max','f8','node')
    nc.variables[varname+'_max'].long_name = 'maximum ' + longname
    nc.variables[varname+'_max'].units = varunits
    
    nc.createVariable('time_of_'+varname+'_max','f8','node')
    nc.variables['time_of_'+varname+'_max'].long_name = 'time of maximum ' + longname
    nc.variables['time_of_'+varname+'_max'].units = 'sec'
    
    # Store last time-step in run for reference
    nc.variables['time'][0] = np.float64(fobj.readline().split()[0])

    # Store max variable observations 
    for aa in range(nodes):
        nc.variables[varname+'_max'][aa] = np.float64(fobj.readline().split()[1])
    
    # Skip repeated time 
    fobj.readline()
    
    for aa in range(nodes):
        nc.variables['time_of_'+varname+'_max'][aa] = np.float64(fobj.readline().split()[1])
    
    # All done here
    fobj.close()
    nc.close()
    
    
# ==============================================================================
# Read all ASCII files and save as nc file
# ==============================================================================
def all_ascii2nc(runFld,ncFld,ncdate='0000-00-00 00:00:00 UTC'):
    """ 
    Reads (most*) ADCIRC ASCII files and stores them in netcdf4 files

    Input:
    ------
    runFld: folder containing all ascii files
    ncFld:  folder where to save all nc files
    
    Notes:
    ------        
    *Still working on converting some adcirc files/formats
    """
    # Read input file to retrieve station information
    fort15 = adcpre.readsta_fort15(runFld +'fort.15')
    
    # Unstructured grid + bathy ------------------------------------------------
    if os.path.exists(runFld+'fort.14'):
        if not os.path.exists(ncFld+'fort.14.nc'):
            print('Creating fort.14.nc')
            adcpre.fort14_to_nc(runFld + 'fort.14',savename=ncFld+'fort.14.nc')
    
    # fort.61-type files (scalars) ---------------------------------------------
    # Water surface elevation
    if os.path.exists(runFld+'fort.61'):
        if not os.path.exists(ncFld+'fort.61.nc'):
            print('Creating fort.61.nc')
            fort61_to_nc(runFld + 'fort.61',fort15['nameel'],fort15['xel'],
                         fort15['yel'],ncdate=ncdate,
                         savename=ncFld+'fort.61.nc')
            
    # fort.63-type files (scalars) ---------------------------------------------
    # Water surface elevation
    if os.path.exists(runFld+'fort.63'):
        if not os.path.exists(ncFld+'fort.63.nc'):
            print('Creating fort.63.nc')
            fort63_to_nc(runFld + 'fort.63',ncdate=ncdate,
                         savename=ncFld+'fort.63.nc')
    # Atmospherec pressure
    if os.path.exists(runFld+'fort.73'):
        if not os.path.exists(ncFld+'fort.73.nc'):
            print('Creating fort.73.nc')
            fort63_to_nc(runFld + 'fort.73',varname='pressure',
                         longname='air pressure at sea level',
                         varunits='meters of waver',ncdate=ncdate,
                         savename=ncFld+'fort.73.nc')
    # Significant wave height
    if os.path.exists(runFld+'swan_HS.63'):
        if not os.path.exists(ncFld+'swan_HS.63.nc'):
            print('Creating swan_HS.63.nc')
            fort63_to_nc(runFld + 'swan_HS.63',varname='swan_HS',
                         longname='significant wave height',
                         varunits='m',ncdate=ncdate,
                         savename=ncFld+'swan_HS.63.nc')
    if os.path.exists(runFld+'swan_HS_max.63'):
        if not os.path.exists(ncFld+'swan_HS_max.63.nc'):
            print('Creating swan_HS_max.63.nc')
            fort63_to_nc(runFld + 'swan_HS_max.63',varname='swan_HS_max',
                        longname='maximum significant wave height',
                        varunits='m',ncdate=ncdate,
                        savename=ncFld+'swan_HS_max.63.nc')
    # Mean wave direction
    if os.path.exists(runFld+'swan_DIR.63'):
        if not os.path.exists(ncFld+'swan_DIR.63.nc'):
            print('Creating swan_DIR.63.nc')
            fort63_to_nc(runFld + 'swan_DIR.63',varname='swan_DIR',
                         longname='mean wave direction',
                         varunits='degrees',ncdate=ncdate,
                         savename=ncFld+'swan_DIR.63.nc')
    if os.path.exists(runFld+'swan_DIR_max.63'):
        if not os.path.exists(ncFld+'swan_DIR_max.63.nc'):
            print('Creating swan_DIR_max.63.nc')
            fort63_to_nc(runFld + 'swan_DIR_max.63',varname='swan_DIR_max',
                         longname='maximum mean wave direction',
                         varunits='degrees',ncdate=ncdate,
                         savename=ncFld+'swan_DIR_max.63.nc')
    # Mean absolute wave period
    if os.path.exists(runFld+'swan_TMM10.63'):
        if not os.path.exists(ncFld+'swan_TMM10.63.nc'):
            print('Creating swan_TMM10.63.nc')
            fort63_to_nc(runFld + 'swan_TMM10.63',varname='swan_TMM10',
                         longname='mean absolute wave period',
                         varunits='s',ncdate=ncdate,
                         savename=ncFld+'swan_TMM10.63.nc')
    if os.path.exists(runFld+'swan_TMM10_max.63'):
        if not os.path.exists(ncFld+'swan_TMM10_max.63.nc'):
            print('Creating swan_TMM10_max.63.nc')
            fort63_to_nc(runFld + 'swan_TMM10_max.63',varname='swan_TMM10_max',
                         longname='maximum TMM10 mean wave period',
                         varunits='s',ncdate=ncdate,
                         savename=ncFld+'swan_TMM10_max.63.nc')
    # Smoothed peak period
    if os.path.exists(runFld+'swan_TPS.63'):
        if not os.path.exists(ncFld+'swan_TPS.63.nc'):
            print('Creating swan_TPS.63.nc')
            fort63_to_nc(runFld + 'swan_TPS.63',varname='swan_TPS',
                         longname='smoothed peak period',
                         varunits='s',ncdate=ncdate,
                         savename=ncFld+'swan_TPS.63.nc')
    if os.path.exists(runFld+'swan_TPS_max.63'):
        if not os.path.exists(ncFld+'swan_TPS_max.63.nc'):
            print('Creating swan_TPS_max.63.nc')
            fort63_to_nc(runFld + 'swan_TPS_max.63',varname='swan_TPS_max',
                         longname='maximum smoothed peak period',
                         varunits='s',ncdate=ncdate,
                         savename=ncFld+'swan_TPS_max.63.nc')
    
    # fort.64-type files (vectors) ---------------------------------------------
    # Depth average velocity 
    if os.path.exists(runFld+'fort.64'):
        if not os.path.exists(ncFld+'fort.64.nc'):
            print('Creating fort.64.nc')
            fort64_to_nc(runFld + 'fort.64',ncdate=ncdate,savename=ncFld+'fort.64.nc')
    # Wind stress or velocity 
    if os.path.exists(runFld+'fort.74'):
        if not os.path.exists(ncFld+'fort.74.nc'):
            print('Creating fort.74.nc')
            fort64_to_nc(runFld + 'fort.74',varname_xy=['windx','windy'],
                         longname='wind',ncdate=ncdate,savename=ncFld+'fort.74.nc')
    
    # max.63-type files (max/min values) ---------------------------------------
    # Water surface elevation
    if os.path.exists(runFld+'maxele.63'):
        if not os.path.exists(ncFld+'maxele.63.nc'):
            print('Creating maxele.63.nc')
            max63_to_nc(runFld + 'maxele.63',ncdate=ncdate,savename=ncFld+'maxele.63.nc')
    if os.path.exists(runFld+'maxvel.63'):
        if not os.path.exists(ncFld+'maxvel.63.nc'):
            print('Creating maxvel.63.nc')
            max63_to_nc(runFld + 'maxvel.63',varname='vel',
                        longname='water velocity',
                        varunits='m s-1',ncdate=ncdate,savename=ncFld+'maxvel.63.nc')
    if os.path.exists(runFld+'maxwvel.63'):
        if not os.path.exists(ncFld+'maxwvel.63.nc'):
            print('Creating maxwvel.63.nc')
            max63_to_nc(runFld + 'maxwvel.63',varname='wind',
                        longname='wind velocity',
                        varunits='m s-1',ncdate=ncdate,savename=ncFld+'maxwvel.63.nc')