"""
Tools to analyse WW3 output

Authors:
-------
Gabriel Garcia Medina
    Nearshore Modeling Group
    ggarcia@coas.oregonstate.edu
Saeed Moghimi

Log of edits:
-------------
August 2014 - Created module
    Gabriel Garcia Medina

Dependencies:
-------------
    numpy, time, datetime, sys, os, re, netCDF4, collections, getpass   
  
Internal dependencies:
----------------------
    angles
"""

from __future__ import division,print_function

__author__ = "Gabriel Garcia Medina"
__email__ = "ggarcia@coas.oregonstate.edu"
__group__ = 'Nearshore Modeling Group'

# Import Modules
import pylab as pl
import time
import datetime
import sys,os
import numpy as np
import re
import netCDF4
from collections import defaultdict
import getpass

# My functions
import pynmd.data.angles as gangles


#===============================================================================
# Read wave spectrum 
#===============================================================================
def read_spec(specfile,bulkparam=True):
    '''
    Read wave spectrum from a text file produced by WAVEWATCH and returns the
    spectral density and axes. Usage:
    
    read_spec(specfile)
    
    Parameters
    ----------
    specfile       : Full path to the spectra file to be read (string)                    
    bulkparam      : Flag to compute some bulk parameters (defaults to True)
    
    Returns
    -------
    ww3spec: dictionary
             Contains
               'name          : station name
               'latitude'     : decimal latitudes
               'longitude'    : decimal longitudes
               'dpt'          : depth time series
               'wnd'          : 10m wind at the station
               'wnddir'       : wind direction
               'cur'          : flow velocity at the station
               'curdir'       : flow direction
               'direction'    : spectral directions in radians from true north
               'frequency'    : spectral frequencies in Hz
               'spec'         : 4D array of spectral density data
               'wavetime'     : datetime array
    
    if bulkmaram is True:
               'Hsig'         : significant wave height [m]
               'sw'           : spectral width (m0*m2/m1/m1 - 1)**2
               'te'           : energy period [s] (m-1/m0)
               'tm01'         : first moment mean wave period [s] (m0/m1)
               'tp'           : peak wave period of frequency spectrum [s]
               'tp_fit'       : peak wave period computed from second order 
                                polynomial fit near tp [s]
               'mwd'          : average wave direction [rad] (see Kuik 1988)
               'tp_2d'        : peak wave period from frequency-direction 
                                spectrum [s]
               'wdir_tp_2d'   : direction of tp_2d [rad]
               
    References
    ----------
    Kuik, A. J., G. P. van Vledder, and L. H. Holthuijsen, 1988: A method for
    the routine analysis of pitch-and-roll buoy wave data. Journal of Physical
    Oceanography, 18, 1020 - 1034.

    Notes
    -----
    1. WW3 provides the directions in a modified natical convention where the 
       direction waves are going to (measured from true north) is given. 
       a. To convert WW3 to Cartesian (degrees)
          cart = 450 - WW3
       b. To convert WW3 to traditional nautical convention
          naut = 180 + WW3 
    '''
    
    # For code development only -----------------------------------------------
    # specfile = ('/home/shusin2/users/ggarcia/projects2/tillamook-climate/'+
    #             '30-gfs-fabrice/198011/out/ww3.80112100.buoy.spc.enp')
    # bulkparam = True 
    # -------------------------------------------------------------------------
    
    # Find the length of the file and get dates
    # Add open error here and exit the code
    fobj = open(specfile,'r')    
    ldates = re.findall(r'\d{8}\s\d{6}',fobj.read())
    fobj.close()
    
    # Get Datetime array
    wavetime = [datetime.datetime.strptime(x,"%Y%m%d %H%M%S") for x in ldates]
    # wavetime = np.asarray(wavetime)
    
    # Open the file
    fobj = open(specfile,'r')
    
    # Read file information
    tmpline = fobj.readline().split()
    nfreq   = int(tmpline[3])
    ndir    = int(tmpline[4])
    npts    = int(tmpline[5])
    
    # Get the frequencies    
    tmpline = fobj.readline().split()       # Get the first line
    frequency = [float(x) for x in tmpline] # Allocate the first line
    
    # Verify if all the frequencies are in the first line
    if nfreq > len(tmpline):
        for aa in range(int(np.ceil(np.double(nfreq)/len(tmpline))-1)):
            tmpline = fobj.readline().split()
            tmpfreq = [float(x) for x in tmpline]
            frequency.extend(tmpfreq)
            del tmpfreq,tmpline
    frequency = np.asarray(frequency)
    
    # Read the directions
    tmpline = fobj.readline().split()       
    direction = [float(x) for x in tmpline] 
    
    # Verify if all the frequencies are in the first line
    if ndir > len(tmpline):
        for aa in range(int(np.ceil(np.double(ndir)/len(tmpline))-1)):
            tmpline = fobj.readline().split()
            tmpdir = [float(x) for x in tmpline]
            direction.extend(tmpdir)
            del tmpdir,tmpline
    direction = np.asarray(direction)    
    
    # Sort directions
    sortind = np.argsort(direction)
    direction = direction[sortind]
    
    # Direction vector in degrees
    # dir_degree = direction * 180.0/pi
    
    # Preallocate/Initialize variables
    spec         = np.zeros((npts,len(wavetime),ndir,nfreq))
    station_name = [None]*npts
    latitude     = np.zeros(npts)
    longitude    = np.zeros(npts)
    dpt          = np.zeros((npts,len(wavetime)))
    wnd          = np.zeros((npts,len(wavetime)))
    wnddir       = np.zeros((npts,len(wavetime)))
    cur          = np.zeros((npts,len(wavetime)))
    curdir       = np.zeros((npts,len(wavetime)))
    
    # Loop over the file to get spectral data
    for aa in range(len(wavetime)):
        
        # Read date line
        tmpline = fobj.readline()
        
        # End of file reached
        if tmpline == '':
            break                 
        
        # Loop over data points
        for bb in range(npts):
            
            # Read point information
            tmpline = fobj.readline()
            
            # Store Point information
            tmpstr = tmpline.split()
            
            # Station Information
            # Note: Wavewatch III does not provide a space between the latitude
            #       and the longitude if the longitude is a negative number
            #       and has 5 digits (i.e. longitude <= -100.00)
            station_name[bb] = tmpstr[0][1:-1]
            if len(tmpstr) == 7:
                latitude[bb]    = float(tmpstr[1].rsplit('-')[0])
                longitude[bb]   = -1.0*float(tmpstr[1].rsplit('-')[1])
                dpt[bb,aa]      = float(tmpstr[2])
                wnd[bb,aa]      = float(tmpstr[3])
                wnddir[bb,aa]   = float(tmpstr[4])
                cur[bb,aa]      = float(tmpstr[5])
                curdir[bb,aa]   = float(tmpstr[6])
            else:
                latitude[bb]    = float(tmpstr[1])
                longitude[bb]   = float(tmpstr[2])
                dpt[bb,aa]      = float(tmpstr[3])
                wnd[bb,aa]      = float(tmpstr[4])
                wnddir[bb,aa]   = float(tmpstr[5])
                cur[bb,aa]      = float(tmpstr[6])
                curdir[bb,aa]   = float(tmpstr[7])
                                               
            # Read and allocate spectral data
            tmpline = fobj.readline().split() 
            tmpspec = [float(x) for x in tmpline]                      
            if nfreq*ndir > len(tmpline):                

                if len(tmpline) == nfreq:
                    tmpRange = range(ndir-1)                
                else:
                    tmpSize = np.double(nfreq)*np.double(ndir)
                    tmpRange = range(int(np.floor(tmpSize/len(tmpline))))

                for cc in tmpRange:
                    tmpline = fobj.readline().split()
                    tmpflt  = [float(x) for x in tmpline]
                    tmpspec.extend(tmpflt)
                    del tmpflt,tmpline
            tmpspec = np.asarray(tmpspec)
            tmpspec = tmpspec.reshape((ndir,nfreq))
            # Sort directions
            spec[bb,aa,:,:] = tmpspec[sortind,:]
             
            # Clean up
            del tmpstr,tmpspec
                          
                      
    # Close file
    fobj.close()
    
    
    # Compute bulk parameters
    if bulkparam:
        
        # Frequency spectra
        # freq_spec = np.trapz(spec,direction,axis=-2) # Incorrect
        dth = np.abs(direction[2] - direction[1])
        freq_spec = np.sum(spec,axis=-2) * dth
        
        # Moments of spectra
        moment0 = np.trapz(freq_spec,frequency,axis=-1)
        moment1 = np.trapz(freq_spec*frequency,frequency,axis=-1)
        moment2 = np.trapz(freq_spec*(frequency)**2,frequency,axis=-1)
        momentn1 = np.trapz(freq_spec*(frequency)**-1,frequency,axis=-1)
                       
        # Significant wave height
        Hsig = 4.004 * (moment0)**0.5
               
        # Spectral width
        sw = (moment0 * moment2 / moment1 / moment1 - 1)**0.5
        
        # Energy period
        te = momentn1/moment0
        
        # Mean wave period
        tm01 = moment0 / moment1
                
        # Peak wave period 
        # tp_fit: using a quadratic fit over the three largest freqs
        # tp maximum frequency 
        freq_max_ind = np.argmax(freq_spec,axis=-1)
        tp_fit = np.zeros_like(te)
        tp = np.zeros_like(te)
        
        for aa in range(freq_spec.shape[0]):
            for bb in range(freq_spec.shape[1]):
                
                # Find peak wave frequency from the frequency spectrum
                tp[aa,bb] = frequency[freq_max_ind[aa,bb]]**-1
                
                # Polynomial fit
                if freq_max_ind[aa,bb] == 0:
                    tp_fit[aa,bb] = np.nan
                elif freq_max_ind[aa,bb] == frequency.shape[0]-1:
                    tp_fit[aa,bb] = np.nan
                else:
                    minfreq = freq_max_ind[aa,bb] - 1
                    maxfreq = freq_max_ind[aa,bb] + 2
                    tmp_fit = np.polyfit(frequency[minfreq:maxfreq],
                                         freq_spec[aa,bb,minfreq:maxfreq],
                                         2)
                    tp_fit[aa,bb] = (-1.0 * tmp_fit[1] / (2.0* tmp_fit[0]))**-1
        
        # Direction spectra
        dir_spec = np.trapz(spec,frequency,axis=-1)
        
        # Mean wave direction (Kuik 1988)
        # WW3 outputs spectra using cartesian convention, thus angles will
        # be changed to oceanographic convention.
        a1 = np.trapz(direction,dir_spec*np.cos(np.pi/2 - direction))
        b1 = np.trapz(direction,dir_spec*np.sin(np.pi/2 -direction))
        mwd = gangles.wrapto2pi(np.pi/2 - np.arctan2(b1,a1))
        
        # Peak wave direction
        tp_2d = np.zeros_like(te)
        wdir_tp_2d = np.zeros_like(te)
        for aa in range(spec.shape[0]):
            for bb in range(spec.shape[1]):
                freq_max_ind = np.unravel_index(spec[aa,bb,:,:].argmax(),
                                                spec[aa,bb,:,:].shape)
                tp_2d[aa,bb] = frequency[freq_max_ind[1]]**-1
                wdir_tp_2d[aa,bb] = direction[freq_max_ind[0]]
        
        # Convert spectrum and directions to oceanographic convention
        wdir_tp_2d = gangles.wrapto2pi(np.pi + wdir_tp_2d)
        direction = gangles.wrapto2pi(np.pi + direction)
        sortind = np.argsort(direction)
        direction = direction[sortind]
        spec = spec[:,:,sortind,:]

        # Return values
        return {'name': station_name, 'latitude':latitude, 'longitude':longitude ,
                'dpt':dpt,'wnd':wnd, 'wnddir':wnddir,'cur':cur, 'curdir':curdir,
                'direction':direction,'frequency':frequency,'spec':spec,
                'wavetime':wavetime, 'Hsig':Hsig, 'sw':sw, 'te':te,
                'tm01':tm01, 'tp':tp, 'tp_fit':tp_fit, 'mwd':mwd, 
                'tp_2d':tp_2d, 'wdir_tp_2d':wdir_tp_2d}
        
    else: 
        
        # Convert spectrum and direction to oceanographic convention
        direction = gangles.wrapto2pi(np.pi + direction)
        sortind = np.argsort(direction)
        direction = direction[sortind]
        spec = spec[:,:,sortind,:]        
        
        # Return values
        return {'name': station_name, 'latitude':latitude, 'longitude':longitude ,
                'dpt':dpt,'wnd':wnd, 'wnddir':wnddir,'cur':cur, 'curdir':curdir,
                'direction':direction,'frequency':frequency,'spec':spec,
                'wavetime':wavetime}
        
# =============================================================================
# Read source term output
# =============================================================================
def read_src_term(specfile,version='5.16',stot=True):
    """
    Reads the source term output from Wavewatch III
       
    read_src_term(specfile)
    
    Parameters
    ----------
    specfile       : Full path to the spectra file to be read (string)                    
    
    Returns
    -------
    ww3spec: dictionary
             Contains
               'fileId'       : spectral ouptut file identifier
               'name          : station name
               'latitude'     : decimal latitudes
               'longitude'    : decimal longitudes
               'dpt'          : depth time series
               'wnd'          : 10m wind at the station
               'wnddir'       : wind direction
               'cur'          : flow velocity at the station
               'curdir'       : flow direction
               'direction'    : spectral directions in radians from true north
               'frequency'    : spectral frequencies in Hz               
               'wavetime'     : datetime array
               'sterms'       : Dictionary containing 4D array of source term
                                density data. Keys correspond to the source 
                                term.

    NOTES
    -----
    1. There are issues with the source term formatting in WW3 v5.16

    """
    
    # Find the length of the file and get dates
    # Add open error here and exit the code
    fobj = open(specfile,'r')    
    ldates = re.findall(r'\d{8}\s\d{6}',fobj.read())
    fobj.close()
        
    # Get Datetime array
    wavetime = [datetime.datetime.strptime(x,"%Y%m%d %H%M%S") for x in ldates]
    # wavetime = np.asarray(wavetime)
    
    # Open the file
    fobj = open(specfile,'r')
    
    # Read file information
    # The header includes
    # 'File ID', number of frequencies, number of directions, number of points,
    #   Flags for spectrum, Sin, Snl, Sds, Sbt, Sice, Stot
    tmpline  = fobj.readline()
    fileId   = tmpline.split('\'')[1]
    tmpline  = tmpline.split('\'')[2].split()
    nfreq    = int(tmpline[0])
    ndir     = int(tmpline[1])
    npts     = int(tmpline[2])
    
    # All source term keys
    if version == '5.16':
        stkeysAll = ['SW','Sin','Snl','Sds','Sbt','Sice','Stot']
        #stkeysAll = ['SW','Sin','Snl','Sds','Sbt','Stot'] # WW3 v4.18
        fobj.readline().split() # There is a weird line in v5.16
    else:
        stkeysAll = ['SW','Sin','Snl','Sds','Sbt','Stot'] # WW3 v4.18
    
    # Activated source term keys only
    stkeys = []
    for aa in range(3,len(tmpline)):
        if tmpline[aa] == 'T':
            stkeys.append(stkeysAll[aa-3])            
    
    # The total source term flag is not being written in v5.16.
    if version == '5.16' and stot:
        stkeys.append(stkeysAll[-1])

    # Get the frequencies
    tmpline = fobj.readline().split()       # Get the first line
    frequency = [float(x) for x in tmpline] # Allocate the first line
    
    # Verify if all the frequencies are in the first line
    if nfreq > len(tmpline):
        for aa in range(int(np.ceil(np.double(nfreq)/len(tmpline))-1)):
            tmpline = fobj.readline().split()
            tmpfreq = [float(x) for x in tmpline]
            frequency.extend(tmpfreq)
            del tmpfreq,tmpline
    frequency = np.asarray(frequency)
    
    # Read the directions (in radians and oceanographic convention)
    tmpline = fobj.readline().split()       
    direction = [float(x) for x in tmpline] 
    
    # Verify if all the frequencies are in the first line
    if ndir > len(tmpline):
        for aa in range(int(np.ceil(np.double(ndir)/len(tmpline))-1)):
            tmpline = fobj.readline().split()
            tmpdir = [float(x) for x in tmpline]
            direction.extend(tmpdir)
            del tmpdir,tmpline
    direction = np.asarray(direction)    
    
    # Sort directions
    sortind = np.argsort(direction)
    direction = direction[sortind]
    
    # Direction vector in degrees
    # dir_degree = direction * 180.0/pi
    
    # Preallocate/Initialize variables
    station_name = [None]*npts
    latitude     = np.zeros(npts)
    longitude    = np.zeros(npts)
    dpt          = np.zeros((npts,len(wavetime)))
    wnd          = np.zeros((npts,len(wavetime)))
    wnddir       = np.zeros((npts,len(wavetime)))
    cur          = np.zeros((npts,len(wavetime)))
    curdir       = np.zeros((npts,len(wavetime)))
    
    # The source term dictionary
    sterms = {}
    for aa in stkeys:
        sterms[aa] = np.zeros((npts,len(wavetime),ndir,nfreq))
    
    # Loop over the file to get source terms
    # The looping order works as follows
    # - Time -> aa
    #   - Points -> bb
    #     - Source terms -> dd
    
    for aa in range(len(wavetime)):
        
        # Read date line
        tmpline = fobj.readline()
        
        # End of file reached
        if tmpline == '':
            break                 
        
        # Loop over data points
        for bb in range(npts):
            
            # Read point information
            tmpline = fobj.readline()

            # Remove string markers
            tmpline = tmpline.translate(None,'\'')
            
            # Store Point information
            tmpstr = tmpline.split()
            
            # Station Information
            # Note: Wavewatch III does not provide a space between the latitude
            #       and the longitude if the longitude is a negative number
            #       and has 5 digits (i.e. longitude <= -100.00)
            station_name[bb] = tmpstr[0][:]            
            if len(tmpstr) == 7:
                latitude[bb]    = float(tmpstr[1].rsplit('-')[0])
                longitude[bb]   = -1.0*float(tmpstr[1].rsplit('-')[1])
                dpt[bb,aa]      = float(tmpstr[2])
                wnd[bb,aa]      = float(tmpstr[3])
                wnddir[bb,aa]   = float(tmpstr[4])
                cur[bb,aa]      = float(tmpstr[5])
                curdir[bb,aa]   = float(tmpstr[6])
            else:
                latitude[bb]    = float(tmpstr[1])
                longitude[bb]   = float(tmpstr[2])
                dpt[bb,aa]      = float(tmpstr[3])
                wnd[bb,aa]      = float(tmpstr[4])
                wnddir[bb,aa]   = float(tmpstr[5])
                cur[bb,aa]      = float(tmpstr[6])
                curdir[bb,aa]   = float(tmpstr[7])
                                               
            # Read and allocate the source term parameters
            for dd in range(len(stkeys)):
                           
                tmpline = fobj.readline().split() 
                tmpspec = [float(x) for x in tmpline]                      
                if nfreq*ndir > len(tmpline):
                    for cc in range(int(np.floor(np.double(nfreq)*np.double(ndir)
                                             /len(tmpline)))):
                        tmpline = fobj.readline().split()
                        tmpflt  = [float(x) for x in tmpline]
                        tmpspec.extend(tmpflt)
                        del tmpflt,tmpline
                tmpspec = np.asarray(tmpspec)
                tmpspec = tmpspec.reshape((ndir,nfreq))
                
                # Sort directions
                sterms[stkeys[dd]][bb,aa,...] = tmpspec[sortind,:]
                          
    # Close file
    fobj.close()    
    
    # Convert spectrum and direction to oceanographic convention 
    direction = gangles.wrapto2pi(np.pi + direction)
    sortind = np.argsort(direction)
    direction = direction[sortind]
    for aa in sterms.keys():
        sterms[aa] = sterms[aa][:,:,sortind,:]      
    
    # Return values
    return {'name': station_name, 'latitude':latitude, 'longitude':longitude ,
            'dpt':dpt,'wnd':wnd, 'wnddir':wnddir,'cur':cur, 'curdir':curdir,
            'direction':direction,'frequency':frequency,'sterms':sterms,
            'wavetime':wavetime,'fileId':fileId}
       
#===============================================================================
# Write spectral data in WW3 Format        
#===============================================================================
def write_nc_spec(latitude,longitude,spectrum,frequency,direction,time1,station,
                  fileout=None,spherical=True):
    '''
    
    Parameters
    ----------
    latitude     : Array of the latitudes of each spectra point [Degree]
                   Dimension (time,station)
    longitude    : Array of the longitudes of each spectra point [Degree east]
                   Dimension (time,station)
    spectrum     : Variance spectrum [m2/Hz/rad]
                   Dimensions (time,station,frequency,direction)
    frequency    : frequency at each bin center [Hz]
    direction    : direction of propagation of each spectral bin [rad]
    time         : vector with days from 1900-01-01 00:00:00
    station      : Array with station names (not longer than 16 characters)
    fileout      : Output netCDF file
    spherical    : Flag for metadata. For spherical coordinates True and False
                   for cartesian coordiantes. 
    
    Notes:
    ------
    Direction of where waves are traveling to.
    
    Version:
    --------
    v0.1 code created:
      Gabriel Garcia Medina, Saeed Moghimi, April 2015
                   
    '''
    
    # For testing purposes only ------------------------------------------------   
    #nctest = netCDF4.Dataset('/home/shusin2/shared/nmg/pycodes/'
    #                         'ww3.basin.201401_spec.nc','r')
    #fileout = '/home/shusin2/shared/nmg/pycodes/test2.nc'
    #time1 = nctest.variables['time'][:]
    #station = nctest.variables['station'][:1]
    #latitude = nctest.variables['latitude'][:,:1]   * 0.0  
    #longitude = nctest.variables['longitude'][:,:1] * 0.0+15.0
    #frequency = nctest.variables['frequency'][:]
    #direction = nctest.variables['direction'][:]
    #spectrum = nctest.variables['efth'][:,:1,:,:]
    # -------------------------------------------------------------------------
    
#     # Variable dictionary
#     varinfo = defaultdict(dict)
#     varinfo['frequency']['units'] = 's-1'
#     varinfo['frequency']['long_name'] = 'frequency of center band'
#     varinfo['frequency']['standard_name'] = 'sea_surface_wave_frequency'
#     varinfo['direction']['units'] = 'degree'
#     varinfo['direction']['long_name'] = 'sea surface wave to direction'
#     varinfo['direction']['standard_name'] = 'sea_surface_wave_to_direction'
#     varinfo['latitude']['units'] = 'degree_north'
#     varinfo['latitude']['long_name'] = 'latitude'
#     varinfo['latitude']['standard_name'] = 'latitude'
#     varinfo['longitude']['units'] = 'degree_east'
#     varinfo['longitude']['long_name'] = 'longitude'
#     varinfo['longitude']['standard_name'] = 'longitude'
#     varinfo['efth']['units'] = 'm2 s rad-1'
#     varinfo['efth']['long_name'] = ('sea surface wave directional variance' + 
#                                     ' spectral density')
    
    
    # Create netcdf file
    if fileout:
        nc = netCDF4.Dataset(fileout,'w')
    else:
        print("Writing NetCDF file " + os.getcwd() + "/ww3_spec_inp.nc")
        nc = netCDF4.Dataset('./ww3_spec_inp.nc','w')
        
    # Create general attributes
    nc.Description = 'Wavewatch III 4.18 input spectra file'
    nc.Author = getpass.getuser()
    nc.Created = time.ctime()
    nc.Owner = 'Nearshore Modeling Group (http://ozkan.oce.orst.edu/nmg)'
    nc.Software = 'Created with Python ' + sys.version
    nc.NetCDF_Lib = str(netCDF4.getlibversion())
    nc.Script = os.path.realpath(__file__)            
        
    # Create dimensions
    nc.createDimension('time',0)
    nc.createDimension('frequency',frequency.shape[0])
    nc.createDimension('direction',direction.shape[0])
    nc.createDimension('string16',16)
    nc.createDimension('station',len(station))
    
    # Write time
    nc.createVariable('time','f8',('time'))
    nc.variables['time'].long_name = 'julian day (UT)'
    nc.variables['time'].standard_time = 'time'
    nc.variables['time'].units = 'days since 1900-01-01T00:00:00Z'
    nc.variables['time'].conventions = ('Relative julian days with decimal ' +
                                        'part (as parts of the day)')
    nc.variables['time'][:] = time1
    
    # Write station id and name
    nc.createVariable('station','i8',('station'))
    nc.variables['station'].long_name = 'station id'
    nc.variables['station'].axis = 'X'
    nc.variables['station'][:] = np.arange(0,len(station),1)
    
    # IDK why they do this in WW3 so I am trying to mimic this behaviour
    nc.createVariable('string16','i8',('string16'))
    nc.variables['string16'].long_name = 'station_name number of characters'
    nc.variables['string16'].axis = 'W'
    nc.variables['string16'][:] = np.empty((16,),dtype=int)
    
    # Write station names
    nc.createVariable('station_name','S1',('station','string16'))
    nc.variables['station_name'].long_name = 'station name'
    nc.variables['station_name'].content = 'XW'
    nc.variables['station_name'].associates = 'station string16'
    for aa in range(len(station)):
        nc.variables['station_name'][aa] = station[aa]
    
    # Spatial dimensions
    nc.createVariable('latitude','f8',('time','station'))
    nc.variables['latitude'].long_name = 'latitude'
    nc.variables['latitude'].standard_name = 'latitude'
    nc.variables['latitude'].units = 'degree_north'
    nc.variables['latitude'].valid_min = -90.0
    nc.variables['latitude'].valid_max = 90.0
    nc.variables['latitude'].content = 'TX'
    nc.variables['latitude'].associates = 'time station'
    nc.variables['latitude'][:] = latitude
    
    nc.createVariable('longitude','f8',('time','station'))
    nc.variables['longitude'].long_name = 'longitude'
    nc.variables['longitude'].standard_name = 'longitude'
    nc.variables['longitude'].units = 'degree_east'
    nc.variables['longitude'].valid_min = -180.0
    nc.variables['longitude'].valid_max = 180.0
    nc.variables['longitude'].content = 'TX'
    nc.variables['longitude'].associates = 'time station'
    nc.variables['longitude'][:] = gangles.wrapto180(latitude)
    
    # Create frequency and direction vectors
    nc.createVariable('frequency','f8',('frequency'))
    nc.variables['frequency'].long_name = 'frequency of center band'
    nc.variables['frequency'].standard_name = 'sea_surface_wave_frequency'
    nc.variables['frequency'].globwave_name = 'frequency'
    nc.variables['frequency'].units = 's-1'
    nc.variables['frequency'].valid_min = 0.0
    nc.variables['frequency'].valid_max = 10.0
    nc.variables['frequency'].axis = 'Y'
    nc.variables['frequency'][:] = frequency
    
    nc.createVariable('direction','f8',('direction'))
    nc.variables['direction'].long_name = 'sea surface wave to direction'
    nc.variables['direction'].standard_name = 'sea_surface_wave_to_direction'
    nc.variables['direction'].globwave_name = 'direction'
    nc.variables['direction'].units = 'degree'
    nc.variables['direction'].valid_min = 0.0
    nc.variables['direction'].valid_max = 360.0
    nc.variables['direction'][:] = gangles.wrapto360(direction)
    
    # Write spectral data  
    nc.createVariable('efth','f8',('time', 'station', 'frequency', 'direction'))
    nc.variables['efth'].long_name = ('sea surface wave directional variance' +
                                      ' spectral density')
    nc.variables['efth'].standard_name = ('sea_surface_wave_directional_' + 
                                          'variance_spectral_density')
    nc.variables['efth'].globwave_name = 'directional_variance_spectral_density'
    nc.variables['efth'].units = 'm2 s rad-1'
    nc.variables['efth'].scale_factor = 1.0
    nc.variables['efth'].add_offset = 0.0
    nc.variables['efth'].valid_min = 0.0
    nc.variables['efth'].valid_max = 1.0e20
    nc.variables['efth'].content = 'TXYZ'
    nc.variables['efth'].associates = 'time station frequency direction'
    nc.variables['efth'][:] = spectrum

    nc.close()

# ==============================================================================
# Read Wind Files 
# ==============================================================================
def readWind(windFile,atDiff=False):
    """
    Read wind files in WW3 format
    """

    # Find the length of the file and get dates
    print('Finding time steps')
    fobj = open(windFile,'r')    
    ldates = re.findall(r'\d{8}\s\d{6}',fobj.read())
    fobj.close()
    
    # Set as wavetime
    windTime = [datetime.datetime.strptime(x,"%Y%m%d %H%M%S") for x in ldates]
    windTime = np.array(windTime)

    # Read wind data
    print('Reading wind data')
    uwnd = []      # Zonal Wind
    vwnd = []      # Meridional Wind
    astemp = []    # Air-sea temperature differences
    tmpData = []   # Tmp data container

    # Open the file
    fobj = open(windFile,'r')
    
    # Discard the first line (it contains time information)
    fobj.readline()
    dataFlag = True
    cnt = 0
    while dataFlag:
        
        # Read line
        tmpline = fobj.readline().rstrip().split()
        cnt += 1
        # Another date stamp or file ended
        if len(tmpline) <= 2:
            # Allocate the wind data
            tmpData = np.array(tmpData)
            lats = tmpData.shape[0]
            
            if atDiff:
                ind = np.int(lats/3)
                ind2=ind*2 + 1
                uwnd.append(tmpData[:ind,:])                
                vwnd.append(tmpData[ind:ind2,:])
                astemp.append(tmpData[ind2:,:])
            else:
                ind = np.int(lats/2)
                uwnd.append(tmpData[:ind,:])
                vwnd.append(tmpData[ind:,:])

            # Did we reach end of file
            if len(tmpline[0]) < 1:
                dataFlag = False
                break
            
            # Reset the container
            tmpData = []

        else:
            # Store wind in temporary array
            tmpData.append([np.float(bb) for bb in tmpline])
    
    fobj.close()

    # Generate arrays    
    ww3 = {'ot':windTime,'uwnd':np.array(uwnd),'vwnd':np.array(vwnd)}
    if atDiff:
        ww3['tDiff'] = astemp

    return ww3
