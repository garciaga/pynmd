# -*- coding: utf-8 -*-
"""
Functions to read OWI (Oceanweather Inc.) data files.

Authors:
-------
Fadia Ticona Rollano
    PNNL Marine and Coastal Research Laboratory

Log of edits:
-------------
November 2020 - Created module
    Fadia Ticona Rollano
    
"""
import ntpath
import sys,time
import getpass
import numpy as np
import netCDF4
import datetime

def owi2nc(owiFile,basedate='0001-01-01 00:00:00 UTC',**kwargs):
    """ 
    Script to read wind or pressure grids from a OWI format file and store it
    in a netcdf file.

    PARAMETERS:
    -----------
    owiFile: Path to the input data file. If 'ext' is not given in the kwargs,
             then the extension (from which the data type is determined) is 
             assumed to be the last three letters of the owiFile.
    basedate : reference date and time in format yyyy-MM-dd hh:mm:ss tz
    Optional kwargs:
      'savepath' = path to output netcdf file, defaults to <owiFile_ext + '.nc'>
      'ext' = indicates the data type of the owiFile, either 'pre' for pressure
              or 'win' for wind.
    RETURNS:
    --------
    Netcdf containing
    time    : time vector in seconds since basedate
    lon,lat : [x,y] grid with longitudes and latitudes
    pressure OR u,v: [x,y,time] pressure or wind grid(s)
    """ 
    
    # Read file extension/determine data type
    if 'ext' in kwargs:
        ext = kwargs['ext']
    else:
        ext = owiFile[-3:]
    if ext not in ['pre','win']:
        print("Error! Must provide file extension ")
        sys.exit()
    
    # Open data file
    print('Reading',ntpath.basename(owiFile))
    fobj = open(owiFile,'r')
    
    # Read in begining/ending dates of win file
    line = fobj.readline().split()
    #date1 = datetime.datetime.strptime(line[-2],'%Y%m%d%H')
    #date2 = datetime.datetime.strptime(line[-1],'%Y%m%d%H')
    
    # Create the ncfile
    if 'savepath' in kwargs:
        ncfile = kwargs['savepath']
    else:
        ncfile = owiFile[:-4] + '_' + ext + '.nc'
    nc = netCDF4.Dataset(ncfile, 'w', format='NETCDF4')
    
    # Global attributes 
    nc.Author = getpass.getuser()
    nc.Created = time.ctime()
    nc.description = ' '.join(line[:2])
    nc.Software = 'Created with Python ' + sys.version
    nc.NetCDF_Lib = str(netCDF4.getlibversion())
    
    # Read grid specifications
    line = fobj.readline()
    iLat = int(line[5:9])
    iLon = int(line[15:19])
    dx = float(line[22:28])
    dy = float(line[31:37])
    swLat = float(line[43:51])
    swLon = float(line[57:65])
        
    # Create dimensions
    nc.createDimension('time',0)      # The unlimited dimension
    nc.createDimension('x',iLon)      # Number of longitudes
    nc.createDimension('y',iLat)      # Number of latitudes
    
    # Create time vector
    nc.createVariable('time','f8',('time'))
    nc.variables['time'].long_name = 'model time'
    nc.variables['time'].standard_name = 'time'
    nc.variables['time'].units = 'seconds since ' + basedate
    nc.variables['time'].base_date = basedate
    
    # Create and store spatial variables
    nc.createVariable('lon','f8','x')
    nc.variables['lon'].long_name = 'longitude'
    nc.variables['lon'].units = 'degrees east'
    nc.variables['lon'].positive = 'east'
    nc.variables['lon'][:] = np.asarray([swLon + (a*dx) for a in range(iLon)])
    
    nc.createVariable('lat','f8','y')
    nc.variables['lat'].long_name = 'latitude'
    nc.variables['lat'].units = 'degrees north'
    nc.variables['lat'].positive = 'north'
    nc.variables['lat'][:] = np.asarray([swLat + (a*dy) for a in range(iLat)])
    
    if ext == 'pre':
        nc.createVariable('pressure','f8',('x','y','time'))
        nc.variables['pressure'].long_name = 'atmospheric pressure'
        nc.variables['pressure'].units = 'mbar'
    elif ext == 'win':
        nc.createVariable('u','f8',('x','y','time'))
        nc.variables['u'].long_name = 'west-east velocity'
        nc.variables['u'].units = 'm s-1'
        
        nc.createVariable('v','f8',('x','y','time'))
        nc.variables['v'].long_name = 'south-north velocity'
        nc.variables['v'].units = 'm s-1'
    
    tt = -1
    bd = datetime.datetime.strptime(basedate,'%Y-%m-%d %H:%M:%S %Z')
    while line:
        tt += 1
        # Read grid date
        lCymdHM = datetime.datetime.strptime(line[68:80],'%Y%m%d%H%M')
        nc.variables['time'][tt] = (lCymdHM-bd).total_seconds()
        # Read grid
        ix = 0
        iy = 0
        if ext == 'pre':
            params = ['pressure']
        elif ext == 'win':
            params = ['u','v']
        for p in params:
            for a in range(int(np.ceil(iLat*iLon/8))):
                line = fobj.readline().split()
                for b in line:
                    if ix==iLon:
                        ix = 0
                        iy += 1
                    nc.variables[p][ix,iy,tt] = float(b)
                    ix += 1
            if p == 'u':
                ix = 0
                iy = 0
            else:
                line = fobj.readline()
            
    fobj.close()
    nc.close()