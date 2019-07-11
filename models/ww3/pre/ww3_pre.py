"""
Tools to prepare WW3 input

Authors:
-------
Gabriel Garcia Medina

Log of edits:
-------------
April 2018 - Created module
    Gabriel Garcia Medina

Dependencies:
-------------
    numpy, time, datetime, sys, os, re, netCDF4, collections, getpass   
  
Internal dependencies:
----------------------
    angles
"""

from __future__ import division,print_function

import os as _os
import sys as _sys
import time as _time
import netCDF4 as _netCDF4
import getpass as _getpass
import numpy as _np
import datetime as _datetime

import pynmd.data.angles as _gangles

# ==============================================================================
# Subroutine to write bathymetry file ------------------------------------------
# ==============================================================================
def write_bathy(outfld,lon,lat,depth,spherical=True,mapsta=None):
    '''
    Function to generate bathymetry files to be used as input for WW3 v4.18
    
    PARAMETERS:
    -----------
    outfld        : folder to write the files to
    lon           : longitudes (or x for cartesian) at every grid point 
                    (I assume they have been meshgridded)
    lat           : latitudes (or y for cartesian) at every grid point
    spherical     : True for spherical grids and False for cartesian grids
    mapsta        : Status map (optional)
                    0 = Land Point, 1 = Regular sea point
                    2 = Active boundary point, 3 = point excluded from grid
    
    OUTPUT:
    -------
    ww3_grid.bot  : ASCII file containing the grid depths
    ww3_grid.lon  : ASCII file containing the longitudes (or x locations)
    ww3_grid.lat  : ASCII file containing the latitudes (or y locations)
    ww3_grid.nc   : NetCDF file containing the grid information
    ww3_grid.inp  : Input file template to be used for the grid preprocessor. 
                    You will need to fix many things it will probably not
                    work out of the box.
    ww3_grid.mask : ASCII file containing the maps status if requested
    
    NOTES:   
    ------
    I assume that if the grid is in spherical coordiantes the latitude and 
        longitudes are given in decimal degrees (east). On the other hand if the 
        grid is cartesian I assume meters.
    
    '''
    
    # Convert longitudes to degrees east 
    lon = _gangles.wrapto360(lon)
    
    # Write ascii files --------------------------------------------------------
    # Bathymetry
    fid = open(outfld + '/ww3_grid.bot','w')
    for aa in range(depth.shape[0]):
        for bb in range(depth.shape[1]):
            fid.write('%12.4f' % depth[aa,bb])
        fid.write('\n')
    fid.close()
    
    
    # Latitude File
    fid = open(outfld + '/ww3_grid.lat','w')
    for aa in range(lat.shape[0]):
        for bb in range(lat.shape[1]):
            fid.write('%12.4f' % lat[aa,bb])
        fid.write('\n')
    fid.close()
    
    
    # Longitude File
    fid = open(outfld + '/ww3_grid.lon','w')
    for aa in range(lon.shape[0]):
        for bb in range(lon.shape[1]):
            fid.write('%12.4f' % lon[aa,bb])
        fid.write('\n')
    fid.close()

    # Map status
    if mapsta is not None:

        fid = open(outfld + '/ww3_grid.mask','w')
        for aa in range(mapsta.shape[0]):
            for bb in range(mapsta.shape[1]):
                fid.write('{:3.0f}'.format(mapsta[aa,bb]))
            fid.write('\n')
        fid.close()
    
    # Write netcdf file --------------------------------------------------------
    
    # Global attributes  
    nc = _netCDF4.Dataset(outfld + '/ww3_grid.nc', 'w', format='NETCDF4')
    nc.Description = 'Wavewatch III Bathymetry'
    nc.Author = _getpass.getuser()
    nc.Created = _time.ctime()
    nc.Owner = 'Nearshore Modeling Group (http://ozkan.oce.orst.edu/nmg)'
    nc.Software = 'Created with Python ' + _sys.version
    nc.NetCDF_Lib = str(_netCDF4.getlibversion())
    nc.Script = _os.path.realpath(__file__)           
        
    # Create dimensions
    nc.createDimension('xi_rho', lon.shape[1])
    nc.createDimension('eta_rho',lon.shape[0])
    
    # Write coordinates and depth to netcdf file    
    if spherical:
        nc.createVariable('lat_rho','f8',('eta_rho','xi_rho'))
        nc.variables['lat_rho'].units = 'degree_north'
        nc.variables['lat_rho'].long_name = 'latitude of RHO-points'
        nc.variables['lat_rho'][:] = lat
    else:
        nc.createVariable('y_rho','f8',('eta_rho','xi_rho'))
        nc.variables['y_rho'].units = 'meter'
        nc.variables['y_rho'].long_name = 'y location of RHO-points'
        nc.variables['y_rho'][:] = lat
    
    # Write longitude
    if spherical:
        nc.createVariable('lon_rho','f8',('eta_rho','xi_rho'))
        nc.variables['lon_rho'].units = 'degree_east'
        nc.variables['lon_rho'].long_name = 'latitude of RHO-points'
        nc.variables['lon_rho'][:] = lon
    else:
        nc.createVariable('x_rho','f8',('eta_rho','xi_rho'))
        nc.variables['x_rho'].units = 'meter'
        nc.variables['x_rho'].long_name = 'x location of RHO-points'
        nc.variables['x_rho'][:] = lon
        
    # Write water depth
    nc.createVariable('h','f8',('eta_rho','xi_rho'))
    nc.variables['h'].units = 'meter'
    nc.variables['h'].long_name = 'bathymetry at RHO-points'
    nc.variables['h'][:] = depth

    # Map status
    if mapsta is not None:
        nc.createVariable('MAPSTA','f8',('eta_rho','xi_rho'))
        nc.variables['MAPSTA'].units = '1'
        nc.variables['MAPSTA'].long_name = 'status map'
        nc.variables['MAPSTA'][:] = mapsta

    # Close NetCDF file
    nc.close()
    
    # Write input file for grid preprocessor ----------------------------------
    
    # Create a generic input file
    fid = open(outfld + '/ww3_grid.inp','w')
    
    fid.write('$ ----------------------------------------------------------$\n')
    fid.write('$ WAVEWATCH III Grid preprocessor input file                $\n')
    fid.write('$ ----------------------------------------------------------$\n')
    fid.write('\'Grid Name\'\n')
    fid.write('$\n')
    fid.write('1.1  0.030  40  36  0.\n')
    fid.write('$\n')
    fid.write('F T T T T T\n')
    fid.write('$\n')
    fid.write('3600. 600. 3600. 300.\n')
    fid.write('$\n')
    fid.write('END OF NAMELISTS\n')
    fid.write('$\n')
    fid.write('$ Define grid --------------------------------------------- $\n')
    fid.write('$ Five records containing :\n')
    fid.write('$  1 Type of grid, coordinate system and type of closure: GSTRG, FLAGLL,\n')
    fid.write('$    CSTRG. Grid closure can only be applied in spherical coordinates.\n')
    fid.write('$      GSTRG  : String indicating type of grid :\n')
    fid.write('$               ''RECT''  : rectilinear\n')
    fid.write('$               ''CURV''  : curvilinear\n')
    fid.write('$      FLAGLL : Flag to indicate coordinate system :\n')
    fid.write('$               T  : Spherical (lon/lat in degrees)\n')
    fid.write('$               F  : Cartesian (meters)\n')
    fid.write('$      CSTRG  : String indicating the type of grid index space closure :\n')
    fid.write('$               ''NONE''  : No closure is applied\n')
    fid.write('$               ''SMPL''  : Simple grid closure : Grid is periodic in the\n')
    fid.write('$                         : i-index and wraps at i=NX+1. In other words,\n')
    fid.write('$                         : (NX+1,J) => (1,J). A grid with simple closure\n')
    fid.write('$                         : may be rectilinear or curvilinear.\n')
    fid.write('$               ''TRPL''  : Tripole grid closure : Grid is periodic in the\n')
    fid.write('$                         : i-index and wraps at i=NX+1 and has closure at\n')
    fid.write('$                         : j=NY+1. In other words, (NX+1,J<=NY) => (1,J)\n')
    fid.write('$                         : and (I,NY+1) => (MOD(NX-I+1,NX)+1,NY). Tripole\n')
    fid.write('$                         : grid closure requires that NX be even. A grid\n')
    fid.write('$                         : with tripole closure must be curvilinear.\n')
    fid.write('$  2 NX, NY. As the outer grid lines are always defined as land\n')
    fid.write('$    points, the minimum size is 3x3.\n')
    fid.write('$  3 Unit number of file with x-coordinate.\n')
    fid.write('$    Scale factor and add offset: x <= scale_fac * x_read + add_offset.\n')
    fid.write('$    IDLA, IDFM, format for formatted read, FROM and filename.\n')
    fid.write('$  4 Unit number of file with y-coordinate.\n')
    fid.write('$    Scale factor and add offset: y <= scale_fac * y_read + add_offset.\n')
    fid.write('$    IDLA, IDFM, format for formatted read, FROM and filename.\n')
    fid.write('$  5 Limiting bottom depth (m) to discriminate between land and sea\n')
    fid.write('$    points, minimum water depth (m) as allowed in model, unit number\n')
    fid.write('$    of file with bottom depths, scale factor for bottom depths (mult.),\n')
    fid.write('$    IDLA, IDFM, format for formatted read, FROM and filename.\n')
    fid.write('$      IDLA : Layout indicator :\n')
    fid.write('$                  1   : Read line-by-line bottom to top.\n')
    fid.write('$                  2   : Like 1, single read statement.\n')
    fid.write('$                  3   : Read line-by-line top to bottom.\n')
    fid.write('$                  4   : Like 3, single read statement.\n')
    fid.write('$      IDFM : format indicator :\n')
    fid.write('$                  1   : Free format.\n')
    fid.write('$                  2   : Fixed format with above format descriptor.\n')
    fid.write('$                  3   : Unformatted.\n')
    fid.write('$      FROM : file type parameter\n')
    fid.write('$             ''UNIT'' : open file by unit number only.\n')
    fid.write('$             ''NAME'' : open file by name and assign to unit.\n')
    fid.write('$  If the Unit Numbers in above files is 10 then data is read from this file\n')
    fid.write('$\n')

    # I prefer considering only curvilinear grids.     
    if spherical:
        fid.write('\'CURV\' T \'NONE\'\n')
    else:
        fid.write('\'CURV\' F \'NONE\'\n')

    # Grid size
    fid.write(_np.str(lon.shape[1]) + '  ' + _np.str(lon.shape[0]) + '\n')
    
    # Path to grid name
    fid.write('           11 1.0 0. 1 1 \'(....)\'  \'NAME\'  \'ww3_grid.lon\'\n')
    fid.write('           12 1.0 0. 1 1 \'(....)\'  \'NAME\'  \'ww3_grid.lat\'\n')
    
    # Bottom bathymetry information
    fid.write('$ Bottom bathymetry\n')
    fid.write(' -0.5 0.05 13        1 1 \'(....)\'  \'NAME\'  \'ww3_grid.bot\'\n')
    fid.write('$ Subgrid Information\n')
    fid.write('           10        1 1 \'(....)\'  \'PART\'  \'dummy\'\n')
    
    # Input boundary points
    fid.write('$ Input boundary points\n')
    if mapsta is not None:
        fid.write('           15        1 1 \'(....)\'  \'NAME\'  \'ww3_grid.mask\'\n')
    else:        
        fid.write('  0   0   F\n')
        fid.write('$\n')
        fid.write('  0   0   F\n')
        fid.write('  0   0\n')    

    # Output boundary points
    fid.write('$Output boundary points\n')
    fid.write('  0.  0.  0.  0.  0\n')
    fid.write('$ ----------------------------------------------------------$\n')
    fid.write('$ End of input file                                         $\n')
    fid.write('$ ----------------------------------------------------------$\n')

    # Close input file
    fid.close()

# ==============================================================================
# Subroutine to write bathymetry file ------------------------------------------
# ==============================================================================
def write_wind_nc(outfld,timeVec,lon,lat,uwnd,vwnd,temp=False):
    '''
    Function to generate bathymetry files to be used as input for WW3 v4.18
    
    PARAMETERS:
    -----------
    outfld        : folder to write the files to
    timeVec       : Time vector (datetime array)
    lon           : longitudes (or x for cartesian) at every grid point 
                    (I assume they have been meshgridded)
    lat           : latitudes (or y for cartesian) at every grid point
    uwnd          : zonal component of wind [m/s]
    vwnd          : meridional component of wind [m/s]
    temp          : (optional) Air-sea temperature differences [C]
    
    OUTPUT:
    -------
    ww3_wind.nc   : NetCDF file containing the grid information
    
    NOTES:   
    ------
    I assume that if the grid is in spherical coordiantes the latitude and 
        longitudes are given in decimal degrees (east). On the other hand if the 
        grid is cartesian I assume meters.
    
    '''
    spherical = True # Placeholder for future expansion to cartesian

    # Convert longitudes to degrees east 
    lon = _gangles.wrapto360(lon)
    
    # Write netcdf file --------------------------------------------------------
    
    # Global attributes  
    nc = _netCDF4.Dataset(outfld + '/ww3_wind.nc', 'w', format='NETCDF4')
    nc.Description = 'Wavewatch III Wind Input'
    nc.Author = _getpass.getuser()
    nc.Created = _time.ctime()
    nc.Software = 'Created with Python ' + _sys.version
    nc.NetCDF_Lib = str(_netCDF4.getlibversion())
    nc.Script = _os.path.realpath(__file__)           
        
    # Create dimensions
    nc.createDimension('longitude', lon.shape[1])
    nc.createDimension('latitude',lon.shape[0])
    nc.createDimension('time',0)
    
    # Write coordinates and depth to netcdf file    
    if spherical:
        nc.createVariable('latitude','f8',('latitude','longitude'))
        nc.variables['latitude'].units = 'degree_north'
        nc.variables['latitude'].long_name = 'latitude of RHO-points'
        nc.variables['latitude'][:] = lat
    else:
        nc.createVariable('y_rho','f8',('eta_rho','xi_rho'))
        nc.variables['y_rho'].units = 'meter'
        nc.variables['y_rho'].long_name = 'y location of RHO-points'
        nc.variables['y_rho'][:] = lat
    
    # Write longitude
    if spherical:
        nc.createVariable('longitude','f8',('latitude','longitude'))
        nc.variables['longitude'].units = 'degree_east'
        nc.variables['longitude'].long_name = 'latitude of RHO-points'
        nc.variables['longitude'][:] = lon
    else:
        nc.createVariable('x_rho','f8',('eta_rho','xi_rho'))
        nc.variables['x_rho'].units = 'meter'
        nc.variables['x_rho'].long_name = 'x location of RHO-points'
        nc.variables['x_rho'][:] = lon

    # Time vector 
    nc.createVariable('time','f8',('time'))
    nc.variables['time'].units = 'seconds since 1900-01-01 00:00:00'
    nc.variables['time'].conventions = 'relative julian seconds'
    baseTime = _datetime.datetime(1900,1,1)
    timeVec = _np.array([(aa - baseTime).total_seconds() for aa in timeVec])
    nc.variables['time'][:] = timeVec
           
    # Write wind components
    nc.createVariable('UWND','f8',('time','latitude','longitude'))
    nc.variables['UWND'].units = 'meter second-1'
    nc.variables['UWND'].long_name = 'Zonal component of wind'
    nc.variables['UWND'][:] = uwnd

    nc.createVariable('VWND','f8',('time','latitude','longitude'))
    nc.variables['VWND'].units = 'meter second-1'
    nc.variables['VWND'].long_name = 'Meridional component of wind'
    nc.variables['VWND'][:] = vwnd

    if temp:
        nc.createVariable('temp','f8',('time','latitude','longitude'))
        nc.variables['temp'].units = 'degrees C'
        nc.variables['temp'].long_name = 'Air-sea temperature differences'
        nc.variables['temp'][:] = temp
        
    # Close NetCDF file
    nc.close()
    
# ==============================================================================
# Function to write spectrum in WW3 format
# ==============================================================================
def createWW3AcsiiSpec(outFile,specDict):
    """
    Document me
    """

    # Degree to radian
    deg2rad = _np.pi / 180.0
    specDict['spec'] *= (deg2rad**-1)
    dirs = _gangles.wrapto2pi(specDict['dirs']*deg2rad + _np.pi)

    # Create the output file
    fid = open(outFile,'w')
    
    # Write header
    fid.write('\'' + specDict['fileId'] + '\'' + 
              '{:7.0f}'.format(specDict['freq'].shape[0]) + 
              '{:7.0f}'.format(specDict['dirs'].shape[0]) + 
              '{:7.0f}'.format(specDict['lat'].shape[0]) + 
              ' \'' + specDict['gridName'] + '\'' + '\n')

    # Write the frequencies
    filas = _np.int(_np.floor(specDict['freq'].shape[0]/8.0))
    if _np.mod(specDict['freq'].shape[0],8.0) > 0:
        filas += 1
    cnt = -1
    for aa in range(filas):
        for bb in range(8):
            cnt += 1
            if cnt == specDict['freq'].shape[0]:
                break
            else:
                fid.write('{:9.3e}'.format(specDict['freq'][cnt]) + ' ')
        fid.write('\n')
    
    # Write the directions
    filas = _np.int(_np.floor(specDict['dirs'].shape[0]/7.0))
    if _np.mod(specDict['dirs'].shape[0],7.0) > 0:
        filas += 1
    cnt = -1
    for aa in range(filas):
        for bb in range(7):
            cnt += 1
            if cnt == specDict['dirs'].shape[0]:
                break
            else:
                fid.write('{:9.3e}'.format(dirs[cnt]))
                fid.write('  ')

        fid.write('\n')
    
    # Write the spectrum
    specPnt = specDict['dirs'].shape[0] * specDict['freq'].shape[0]
    filas = _np.int(_np.floor(specPnt/7.0))
    if _np.mod(specPnt,7.0) > 0:
        filas += 1    
    for aa in range(specDict['time'].shape[0]):
        # Write time
        fid.write(specDict['time'][aa].strftime('%Y%m%d %H%M%S') + '\n')

        # Loop over points
        for bb in range(specDict['lat'].shape[0]):
            fid.write('\'' + '{:10s}'.format(specDict['name'][bb]) + '\'  ' + 
                      '{:6.2f}'.format(specDict['lat'][bb]) + 
                      '{:7.2f}'.format(specDict['lon'][bb]) + 
                      '{:10.2f}'.format(1000.0) + 
                      '{:7.2f}'.format(0) + 
                      '{:7.2f}'.format(270) +
                      '{:7.2f}'.format(0) + 
                      '{:7.2f}'.format(270) + '\n')

            # Loop over spectrum
            tmpSpec = specDict['spec'][aa,bb,...].flatten('F')
            cnt = -1            
            for cc in range(filas):
                for dd in range(7):
                    cnt += 1
                    if cnt == specPnt:
                        break
                    else:                        
                        fid.write('{:9.3e}'.format(tmpSpec[cnt]))
                        fid.write('  ')
                fid.write('\n')

    # Close the file
    fid.close()
