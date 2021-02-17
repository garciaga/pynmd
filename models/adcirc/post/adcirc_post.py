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
import pynmd.data.angles as gangles

# ==============================================================================
# Read fort.42-type ASCII files and save as nc file
# ==============================================================================
def fort42_to_nc(fort42,varnames=['u','v','w'],units=['m s-1','m s-1','m s-1'],
                 longnames=['u velocity component','v velocity component',
                            'vertical velocity component'],
                 ncdate='0001-01-01 00:00:00 UTC',**kwargs):
    """
    Script to read fort.42-type (3d station vector) files and store in a
    netcdf4 file

    PARAMETERS:
    -----------
    fort42 : Path to unit 41 file
    x,y    : Station coordinates
    ncdate : cold start date/time in CF standard: yyyy-MM-dd hh:mm:ss tz

    RETURNS:
    --------
    Netcdf containing
    time     : seconds since beginning of run
    variable : temporal variable (called 'varname') recorded at stations.
               Size: [time,station]
    """

    fobj = open(fort42,'r')

    # Create the file
    if 'savepath' in kwargs:
        ncfile = kwargs['savepath']
    else:
        ncfile = fort42 + '.nc'
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
#    ndsets = np.int(tmpline[0])
    nsta = np.int(tmpline[1])
    nfen = np.int(tmpline[4])

    # Create dimensions
    nc.createDimension('time',0)           # The unlimited dimension
    nc.createDimension('nsta',nsta)        # Number of station nodes
    nc.createDimension('vnode',nfen)       # Number of vertical nodes/layers

    # Create time vector
    nc.createVariable('time','f8',('time'))
    nc.variables['time'].long_name = 'model time'
    nc.variables['time'].standard_name = 'time'
    nc.variables['time'].units = 'seconds since ' + ncdate
    nc.variables['time'].base_date = ncdate

    # Create and store spatial variables
    nc.createVariable('sigma','f8','vnode')
    nc.variables['sigma'].long_name = 'level of the vertical grid from -1 (bottom) to +1 (surface)'
    nc.variables['sigma'].units = 'dimensionless'

#    nc.createVariable('x','f8','nsta')
#    nc.variables['x'].long_name = 'longitude'
#    nc.variables['x'].units = 'degrees east'
#    nc.variables['x'].positive = 'east'
#    nc.variables['x'][:] = x
#
#    nc.createVariable('y','f8','nsta')
#    nc.variables['y'].long_name = 'latitude'
#    nc.variables['y'].units = 'degrees north'
#    nc.variables['y'].positive = 'north'
#    nc.variables['y'][:] = y

    # Create station variables
    nc.createVariable(varnames[0],'f8',('time','nsta','vnode'))
    nc.variables[varnames[0]].long_name = longnames[0]
    nc.variables[varnames[0]].units = units[0]

    nc.createVariable(varnames[1],'f8',('time','nsta','vnode'))
    nc.variables[varnames[1]].long_name = longnames[1]
    nc.variables[varnames[1]].units = units[1]

    nc.createVariable(varnames[2],'f8',('time','nsta','vnode'))
    nc.variables[varnames[2]].long_name = longnames[2]
    nc.variables[varnames[2]].units = units[2]

    # Store time-series variables
#    for tt in range(ntsteps):    #this would work if there are no errors with the simulation
#        nc.variables['time'][tt] = np.float64(fobj.readline().split()[0])
#        for aa in range(nsta):
#            nc.variables[varname][tt,aa] = np.float64(fobj.readline().split()[1])
    tt = -1
    line = fobj.readline()
    nc.variables['sigma'][:] = np.asarray([float(b) for b in line.split()[2:nfen*3:3]])
    while line:
        tt += 1
        nc.variables['time'][tt] = np.float64(line.split()[0])
        line = fobj.readline()
        for aa in range(nsta):
            if line:
                nc.variables[varnames[0]][tt,aa,:] = np.float64(line.split()[1:nfen*3:3])
                nc.variables[varnames[1]][tt,aa,:] = np.float64(line.split()[2:nfen*3:3])
                nc.variables[varnames[2]][tt,aa,:] = np.float64(line.split()[3:nfen*3+1:3])
                line = fobj.readline()

    # All done here
    fobj.close()
    nc.close()


# ==============================================================================
# Read fort.61-type ASCII files and save as nc file
# ==============================================================================
def fort61_to_nc(fort61,staname,x,y,varname='zeta',
                 longname='water surface elevation above geoid',
                 varunits='m',ncdate='0001-01-01 00:00:00 UTC',**kwargs):
    """
    Script to read fort.61-type (station scalar) files and store in a netcdf4 file

    PARAMETERS:
    -----------
    fort61 : Path to unit 61 file
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

    # Create the file
    if 'savepath' in kwargs:
        ncfile = kwargs['savepath']
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
#    ntsteps = np.int(tmpline[0])
    nsta = np.int(tmpline[1])

    # Create dimensions
    nc.createDimension('time',0)           # The unlimited dimension
    nc.createDimension('nsta',nsta)        # Number of stations
    nc.createDimension('namelen',50)       # Length of station names

    # Create time vector
    nc.createVariable('time','f8',('time'))
    nc.variables['time'].long_name = 'model time'
    nc.variables['time'].standard_name = 'time'
    nc.variables['time'].units = 'seconds since ' + ncdate
    nc.variables['time'].base_date = ncdate

    # Create and store spatial variables
    nc.createVariable('station_name','S1',('nsta','namelen'))
    nc.variables['station_name'].long_name = 'station name'
    nc.variables['station_name'][:] = netCDF4.stringtochar(np.array(staname,dtype='S50'))
    nc.createVariable('x','f8','nsta')
    nc.variables['x'].long_name = 'longitude'
    nc.variables['x'].units = 'degrees east'
    nc.variables['x'].positive = 'east'
    nc.variables['x'][:] = x
    nc.createVariable('y','f8','nsta')
    nc.variables['y'].long_name = 'latitude'
    nc.variables['y'].units = 'degrees north'
    nc.variables['y'].positive = 'north'
    nc.variables['y'][:] = y

    # Create station variable
    nc.createVariable(varname,'f8',('time','nsta'))
    nc.variables[varname].long_name = longname
    nc.variables[varname].units = varunits

    # Store time-series variables
#    for tt in range(ntsteps):    #this would work if there are no errors with the simulation
#        nc.variables['time'][tt] = np.float64(fobj.readline().split()[0])
#        for aa in range(nsta):
#            nc.variables[varname][tt,aa] = np.float64(fobj.readline().split()[1])
    tt = -1
    line = fobj.readline()
    while line:
        tt += 1
        nc.variables['time'][tt] = np.float64(line.split()[0])
        line = fobj.readline()
        for aa in range(nsta):
            if line:
                nc.variables[varname][tt,aa] = np.float64(line.split()[1])
                line = fobj.readline()

    # All done here
    fobj.close()
    nc.close()


# ==============================================================================
# Read fort.63-type ASCII files and save as nc file
# ==============================================================================
def fort63_to_nc(fort63,varname='zeta',
                 longname='water surface elevation above geoid',
                 varunits='m',ncdate='0001-01-01 00:00:00 UTC',**kwargs):
    """
    Script to read fort.63-type (scalar) files and store in a netcdf4 file

    PARAMETERS:
    -----------
    fort63: Path to fort.63-type file
    ncdate : cold start date/time in CF standard: yyyy-MM-dd hh:mm:ss tz

    RETURNS:
    --------
    Netcdf containing
    time     : seconds since beginning of run
    variable : temporal variable (called 'varname') recorded at the nodes of
               an unstructured grid. Size: [time,nodes]
    """

    fobj = open(fort63,'r')

    # Create the file
    if 'savepath' in kwargs:
        ncfile = kwargs['savepath']
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
#    ntsteps = np.int(tmpline[0])
    nodes = np.int(tmpline[1])

    # Create dimensions
    nc.createDimension('node',nodes)       # Number of nodes
    nc.createDimension('time',0)           # The unlimited dimension

    # Create time vector
    nc.createVariable('time','f8',('time'))
    nc.variables['time'].long_name = 'model time'
    nc.variables['time'].standard_name = 'time'
    nc.variables['time'].units = 'seconds since ' + ncdate
    nc.variables['time'].base_date = ncdate

    # Create the rest of the variables
    nc.createVariable(varname,'f8',('time','node'))
    nc.variables[varname].long_name = longname
    nc.variables[varname].units = varunits

    # Store time-series variables at nodes
#    for tt in range(ntsteps):    #this would work if there are no errors with the simulation
#        nc.variables['time'][tt] = np.float64(fobj.readline().split()[0])
#        for aa in range(nodes):
#            nc.variables[varname][tt,aa] = np.float64(fobj.readline().split()[1])
    tt = -1
    line = fobj.readline()
    while line:
        tt += 1
        nc.variables['time'][tt] = np.float64(line.split()[0])
        line = fobj.readline()
        for aa in range(nodes):
            if line:
                nc.variables[varname][tt,aa] = np.float64(line.split()[1])
                line = fobj.readline()


    # All done here
    fobj.close()
    nc.close()


# ==============================================================================
# Read fort.64-type ASCII files and save as nc file
# ==============================================================================
def fort64_to_nc(fort64,varname_xy=['u-vel','v-vel'],
                 longname='water column vertically averaged',
                 varunits='m s-1',ncdate='0001-01-01 00:00:00 UTC',**kwargs):
    """
    Script to read fort.64-type (vector) files and store in a netcdf4 file

    PARAMETERS:
    -----------
    fort64: Path to unit 64 file
    ncdate : cold start date/time in CF standard: yyyy-MM-dd hh:mm:ss tz

    RETURNS:
    --------
    Netcdf containing
    time         : seconds since beginning of run
    x,y variable : temporal variables (name provided in 'xy_varname') recorded
                   at the nodes of an unstructured grid. Size: [time,nodes]
    """

    fobj = open(fort64,'r')

    # Create the file
    if 'savepath' in kwargs:
        ncfile = kwargs['savepath']
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
#    ntsteps = np.int(tmpline[0])
    nodes = np.int(tmpline[1])

    # Create dimensions
    nc.createDimension('time',0)           # The unlimited dimension
    nc.createDimension('node',nodes)       # Number of nodes

    # Create time vector
    nc.createVariable('time','f8',('time'))
    nc.variables['time'].long_name = 'model time'
    nc.variables['time'].standard_name = 'time'
    nc.variables['time'].units = 'seconds since ' + ncdate
    nc.variables['time'].base_date = ncdate

    # Create the rest of the variables
    nc.createVariable(varname_xy[0],'f8',('time','node'))
    nc.variables[varname_xy[0]].long_name = longname + ' e/w velocity'
    nc.variables[varname_xy[0]].units = varunits

    nc.createVariable(varname_xy[1],'f8',('time','node'))
    nc.variables[varname_xy[1]].long_name = longname + ' n/s velocity'
    nc.variables[varname_xy[1]].units = varunits

    # Store variables
#    for tt in range(ntsteps):    #this would work if there are no errors with the simulation
#        nc.variables['time'][tt] = np.float64(fobj.readline().split()[0])
#        for aa in range(nodes):
#            tmpline = fobj.readline().split()
#            nc.variables[varname_xy[0]][tt,aa] = np.float64(tmpline[1])
#            nc.variables[varname_xy[1]][tt,aa] = np.float64(tmpline[2])
    tt = -1
    line = fobj.readline()
    while line:
        tt += 1
        nc.variables['time'][tt] = np.float64(line.split()[0])
        line = fobj.readline()
        for aa in range(nodes):
            nc.variables[varname_xy[0]][tt,aa] = np.float64(line.split()[1])
            nc.variables[varname_xy[1]][tt,aa] = np.float64(line.split()[2])
            line = fobj.readline()

    # All done here
    fobj.close()
    nc.close()


# ==============================================================================
# Read Global Maximum and Minimum ASCII files and save as nc file
# ==============================================================================
def max63_to_nc(max63,varname='zeta',
                 longname='water surface elevation above geoid',
                 varunits='m',ncdate='0001-01-01 00:00:00 UTC',**kwargs):
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

    # Create the file
    if 'savepath' in kwargs:
        ncfile = kwargs['savepath']
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
    nc.variables['time'].base_date = ncdate

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
# Read ncfile of degree scalars and return nc file of vector components
# ==============================================================================
def ncdirs_to_vec(ncfolder,varname='swan_DIR',
                  source_convention='directions to',
                  target_convention='directions to',
                  **kwargs):
    """
    Script to read a netcdf4 file with wave directions given as scalars
    (in degrees) and write a corresponding netcdf4 file with directions
    given as (u,v) vectors of magnitude 1 in their mean direction.
    Useful to show direction vectors in plots.

    PARAMETERS:
    -----------
    ncfolder: Path to folder containing direction netcdf file. Vector file
              will be stored here as well unless a different 'savepath' path is
              declared in the kwargs
    source,target convention: indicate whether directions are definde as 'to'
                              or 'from' in the source file
    RETURNS:
    --------
    Netcdf containing
    time         : seconds since beginning of run
    x,y variable : temporal variables (name provided in 'xy_varname') recorded
                   at the nodes of an unstructured grid. Size: [time,nodes]
    u,v directions
    NOTES:
    ------
    1) North is defined as a zero degree direction. If direction conventions do
       not match between source and target files, directions are inverted.
    """

    nc_dir = netCDF4.Dataset(ncfolder+varname+'.63.nc', 'r', format='NETCDF4')

    # Create the file
    if 'savepath' in kwargs:
        ncfile = kwargs['savepath']
    else:
        ncfile = ncfolder+varname+'vec.63.nc'
    nc_vec = netCDF4.Dataset(ncfile, 'w', format='NETCDF4')

    # Global attributes
    nc_vec.Author = getpass.getuser()
    nc_vec.Created = time.ctime()
    nc_vec.Software = 'Created with Python ' + sys.version
    nc_vec.NetCDF_Lib = str(netCDF4.getlibversion())

    # Copy additional global attributes from source
    nc_vec.setncatts({a:nc_dir.getncattr(a) for a in nc_dir.ncattrs() if
                      a not in ['creation_date','modification_date','host',
                                'convention','contact']})

    # Create dimensions
    nc_vec.createDimension('time',0)       # The unlimited dimension
    nc_vec.createDimension('node', len(nc_dir.dimensions['node'])) # Number of nodes


    # Copy variables
    for name, var in nc_dir.variables.items():
        if name in ['time','x','y']:
            # Create select vars
            nc_vec.createVariable(name, var.dtype, var.dimensions)
            # Copy the variable attributes
            nc_vec.variables[name].setncatts({a:var.getncattr(a) for a in var.ncattrs()})
            # Copy the variables values
            nc_vec.variables[name][:] = nc_dir.variables[name][:]

    # Create the rest of the variables
    nc_vec.createVariable(varname+'_u','f8',('time','node'))
    nc_vec.variables[varname+'_u'].long_name = 'e/w direction'
    nc_vec.variables[varname+'_u'].units = nc_dir[varname].units
    nc_vec.variables[varname+'_u'].convention = target_convention

    nc_vec.createVariable(varname+'_v','f8',('time','node'))
    nc_vec.variables[varname+'_v'].long_name = 'n/s direction'
    nc_vec.variables[varname+'_v'].units = nc_dir[varname].units
    nc_vec.variables[varname+'_v'].convention = target_convention

    for aa in range(len(nc_vec['time'])):
        if source_convention != target_convention:
            dirs = gangles.wrapto360(nc_dir[varname][aa,:].data + 180) * np.pi/180
        else:
            dirs = nc_dir[varname][aa,:].data  * np.pi/180
        nc_vec.variables[varname+'_u'][aa,:] = np.sin(dirs)
        nc_vec.variables[varname+'_v'][aa,:] = np.cos(dirs)

    # All done here
    nc_dir.close()
    nc_vec.close()


# ==============================================================================
# Read all ASCII files and save as nc file
# ==============================================================================
def all_ascii2nc(runFld,ncFld,ncdate='0001-01-01 00:00:00 UTC'):
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
            adcpre.fort14_to_nc(runFld + 'fort.14',savepath=ncFld+'fort.14.nc')

    # fort.61-type files (scalars) ---------------------------------------------
    # Water surface elevation
    if os.path.exists(runFld+'fort.61'):
        if not os.path.exists(ncFld+'fort.61.nc'):
            print('Creating fort.61.nc')
            fort61_to_nc(runFld + 'fort.61',fort15['nameel'],fort15['xel'],
                         fort15['yel'],ncdate=ncdate,
                         savepath=ncFld+'fort.61.nc')

    # fort.63-type files (scalars) ---------------------------------------------
    # Water surface elevation
    if os.path.exists(runFld+'fort.63'):
        if not os.path.exists(ncFld+'fort.63.nc'):
            print('Creating fort.63.nc')
            fort63_to_nc(runFld + 'fort.63',ncdate=ncdate,
                         savepath=ncFld+'fort.63.nc')
    # Atmospherec pressure
    if os.path.exists(runFld+'fort.73'):
        if not os.path.exists(ncFld+'fort.73.nc'):
            print('Creating fort.73.nc')
            fort63_to_nc(runFld + 'fort.73',varname='pressure',
                         longname='air pressure at sea level',
                         varunits='meters of waver',ncdate=ncdate,
                         savepath=ncFld+'fort.73.nc')
    # Significant wave height
    if os.path.exists(runFld+'swan_HS.63'):
        if not os.path.exists(ncFld+'swan_HS.63.nc'):
            print('Creating swan_HS.63.nc')
            fort63_to_nc(runFld + 'swan_HS.63',varname='swan_HS',
                         longname='significant wave height',
                         varunits='m',ncdate=ncdate,
                         savepath=ncFld+'swan_HS.63.nc')
    if os.path.exists(runFld+'swan_HS_max.63'):
        if not os.path.exists(ncFld+'swan_HS_max.63.nc'):
            print('Creating swan_HS_max.63.nc')
            fort63_to_nc(runFld + 'swan_HS_max.63',varname='swan_HS_max',
                        longname='maximum significant wave height',
                        varunits='m',ncdate=ncdate,
                        savepath=ncFld+'swan_HS_max.63.nc')
    # Mean wave direction
    if os.path.exists(runFld+'swan_DIR.63'):
        if not os.path.exists(ncFld+'swan_DIR.63.nc'):
            print('Creating swan_DIR.63.nc')
            fort63_to_nc(runFld + 'swan_DIR.63',varname='swan_DIR',
                         longname='mean wave direction',
                         varunits='degrees',ncdate=ncdate,
                         savepath=ncFld+'swan_DIR.63.nc')
    if os.path.exists(runFld+'swan_DIR_max.63'):
        if not os.path.exists(ncFld+'swan_DIR_max.63.nc'):
            print('Creating swan_DIR_max.63.nc')
            fort63_to_nc(runFld + 'swan_DIR_max.63',varname='swan_DIR_max',
                         longname='maximum mean wave direction',
                         varunits='degrees',ncdate=ncdate,
                         savepath=ncFld+'swan_DIR_max.63.nc')
    # Mean absolute wave period
    if os.path.exists(runFld+'swan_TMM10.63'):
        if not os.path.exists(ncFld+'swan_TMM10.63.nc'):
            print('Creating swan_TMM10.63.nc')
            fort63_to_nc(runFld + 'swan_TMM10.63',varname='swan_TMM10',
                         longname='mean absolute wave period',
                         varunits='s',ncdate=ncdate,
                         savepath=ncFld+'swan_TMM10.63.nc')
    if os.path.exists(runFld+'swan_TMM10_max.63'):
        if not os.path.exists(ncFld+'swan_TMM10_max.63.nc'):
            print('Creating swan_TMM10_max.63.nc')
            fort63_to_nc(runFld + 'swan_TMM10_max.63',varname='swan_TMM10_max',
                         longname='maximum TMM10 mean wave period',
                         varunits='s',ncdate=ncdate,
                         savepath=ncFld+'swan_TMM10_max.63.nc')
    # Smoothed peak period
    if os.path.exists(runFld+'swan_TPS.63'):
        if not os.path.exists(ncFld+'swan_TPS.63.nc'):
            print('Creating swan_TPS.63.nc')
            fort63_to_nc(runFld + 'swan_TPS.63',varname='swan_TPS',
                         longname='smoothed peak period',
                         varunits='s',ncdate=ncdate,
                         savepath=ncFld+'swan_TPS.63.nc')
    if os.path.exists(runFld+'swan_TPS_max.63'):
        if not os.path.exists(ncFld+'swan_TPS_max.63.nc'):
            print('Creating swan_TPS_max.63.nc')
            fort63_to_nc(runFld + 'swan_TPS_max.63',varname='swan_TPS_max',
                         longname='maximum smoothed peak period',
                         varunits='s',ncdate=ncdate,
                         savepath=ncFld+'swan_TPS_max.63.nc')

    # fort.64-type files (vectors) ---------------------------------------------
    # Depth average velocity
    if os.path.exists(runFld+'fort.64'):
        if not os.path.exists(ncFld+'fort.64.nc'):
            print('Creating fort.64.nc')
            fort64_to_nc(runFld + 'fort.64',ncdate=ncdate,savepath=ncFld+'fort.64.nc')
    # Wind stress or velocity
    if os.path.exists(runFld+'fort.74'):
        if not os.path.exists(ncFld+'fort.74.nc'):
            print('Creating fort.74.nc')
            fort64_to_nc(runFld + 'fort.74',varname_xy=['windx','windy'],
                         longname='wind',ncdate=ncdate,savepath=ncFld+'fort.74.nc')

    # max.63-type files (max/min values) ---------------------------------------
    # Water surface elevation
    if os.path.exists(runFld+'maxele.63'):
        if not os.path.exists(ncFld+'maxele.63.nc'):
            print('Creating maxele.63.nc')
            max63_to_nc(runFld + 'maxele.63',ncdate=ncdate,savepath=ncFld+'maxele.63.nc')
    if os.path.exists(runFld+'maxvel.63'):
        if not os.path.exists(ncFld+'maxvel.63.nc'):
            print('Creating maxvel.63.nc')
            max63_to_nc(runFld + 'maxvel.63',varname='vel',
                        longname='water velocity',
                        varunits='m s-1',ncdate=ncdate,savepath=ncFld+'maxvel.63.nc')
    if os.path.exists(runFld+'maxwvel.63'):
        if not os.path.exists(ncFld+'maxwvel.63.nc'):
            print('Creating maxwvel.63.nc')
            max63_to_nc(runFld + 'maxwvel.63',varname='wind',
                        longname='wind velocity',
                        varunits='m s-1',ncdate=ncdate,savepath=ncFld+'maxwvel.63.nc')

    # fort.42-type files (3D vectors) ------------------------------------------
    if os.path.exists(runFld+'fort.42'):
        if not os.path.exists(ncFld+'fort.42.nc'):
            print('Creating fort.42.nc')
            fort42_to_nc(runFld + 'fort.42',ncdate=ncdate,savepath=ncFld+'fort.42.nc')