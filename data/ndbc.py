"""
Nearshore Modeling Group

Seres of functions to read NDBC text data and store as NetCDF.

Edits:
v0.1 Gabriel Garcia Medina (ggarcia@coas.oregonstate.edu) July 2014
v0.2 Gabriel Garcia Medina (ggarcia@coas.oregonstate.edu) November 2014
v0.3 Gabriel Garcia Medina (ggarcia@coas.oregonstate.edu) February 2015
     Added spec2nc
v0.4 Gabriel Garcia Medina (ggarcia@coas.oregonstate.edu) March 2016
     Fixes to spec2nc 
"""

__author__ = "Gabriel Garcia Medina"
__email__ = "ggarcia@coas.oregonstate.edu"
__group__ = 'Nearshore Modeling Group'

# Import some modules
import time
import numpy as np
import netCDF4
import pylab as pl
import sys
import datetime
import os
import getpass
import glob


#===============================================================================
# pyroms subroutine to write NetCDF fields
#===============================================================================
# def write_nc_var(var, name, dimensions, units=None, longname=None):
#     '''
#     Pyroms subroutine to write NetCDF fields
#     '''
#     nc.createVariable(name, 'f8', dimensions)
#     if units is not None:
#       nc.variables[name].units = units  
#     if longname is not None:
#       nc.variables[name].long_name = longname  
#     nc.variables[name][:] = var


#===============================================================================
# Bulk parameter code
#===============================================================================
def bulk2nc(buoyfld,buoyid,ncformat=4,verbose=True):
    '''
    Code to convert bulk parameter text files into netcdf file
    
    Usage:
    ------
    bulk2nc(buoyfld,buoyid,ncformat)
    
    Input:
    ------
    buoyfld  = Folder where the bulk paramerter text files reside.
    buoyid   = Netcdf buoy identifier (to figure out the file names)
    ncformat = set as 3 for netCDF3, set as 4 for netCDF4 (default)
    verbose  = some extra information (True is default)
    
    Notes:
    Only the bulk parameter files must be present in that directory. The code is
    not smart enough (and I do not have the time to make it so) to figure out 
    the bulk parameter files. 
    
    '''
    
    #===========================================================================
    # Read and Clean Up Data
    #===========================================================================
    
    # Get all files in folder
    archivos = os.listdir(buoyfld)
    archivos.sort()
    
    # Initialize variables
    wavetime = np.empty([1,])
    WDIR = np.empty([1,]) #AKA WD
    WSPD = np.empty([1,])
    GST = np.empty([1,])
    WVHT = np.empty([1,])
    DPD = np.empty([1,])
    APD = np.empty([1,])
    MWD = np.empty([1,])
    PRES = np.empty([1,]) #AKA BAR
    ATMP = np.empty([1,])
    WTMP = np.empty([1,])
    
    
    # Loop over files
    for tmpFile in archivos:
        
        if tmpFile.endswith('.txt'):
            
            if verbose:
                print(tmpFile)

            # Read header lines to determine the location of variables
            with open(buoyfld + '/' + tmpFile,'r') as f:
                header1 = f.readline()
                header2 = f.readline()
                
                # Determine the number of header lines (NDBC uses one or two)
                if header1.startswith('Y') or header1.startswith('#'):
                    headcnt = 1
                if header2.startswith('Y') or header2.startswith('#'):
                    headcnt = 2
                    
                
            # Collapse multiple spaces into one
            header3 = ' '.join(header1.split())
            header4 = header3.split()
            
            # Load buoy data
            tmpdata = pl.loadtxt(buoyfld + '/' + tmpFile,skiprows=headcnt)
            
            
            # ======================  Allocate variables  ==================== #
            # Wind direction
            if any(tt == 'WD' for tt in header4):
                tind = header4.index('WD')
                WDIR = np.concatenate((WDIR,tmpdata[:,tind]))
            elif any(tt == 'WDIR' for tt in header4):
                tind = header4.index('WDIR')
                WDIR = np.concatenate((WDIR,tmpdata[:,tind]))
            else:
                tmparray = np.ones([tmpdata.shape[0],1])
                tmparray[:] = np.NaN
                WDIR = np.concatenate((WDIR,tmparray))
                del tmparray 
            
            # Wind Speed
            if any(tt == 'WSPD' for tt in header4):
                tind = header4.index('WSPD')
                WSPD = np.concatenate((WSPD,tmpdata[:,tind]))
            else:
                tmparray = np.ones([tmpdata.shape[0],1])
                tmparray[:] = np.NaN
                WSPD = np.concatenate((WSPD,tmparray))
                del tmparray 
                
            # Wind Gust
            if any(tt == 'GST' for tt in header4):
                tind = header4.index('GST')
                GST = np.concatenate((GST,tmpdata[:,tind]))
            else:
                tmparray = np.ones([tmpdata.shape[0],1])
                tmparray[:] = np.NaN
                GST = np.concatenate((GST,tmparray))
                del tmparray 
            
            # Wave Height
            if any(tt == 'WVHT' for tt in header4):
                tind = header4.index('WVHT')
                WVHT = np.concatenate((WVHT,tmpdata[:,tind]))
            else:
                tmparray = np.ones([tmpdata.shape[0],1])
                tmparray[:] = np.NaN
                WVHT = np.concatenate((WVHT,tmparray))
                del tmparray 
                
            # Dominant Wave Period
            if any(tt == 'DPD' for tt in header4):
                tind = header4.index('DPD')
                DPD = np.concatenate((DPD,tmpdata[:,tind]))
            else:
                tmparray = np.ones([tmpdata.shape[0],1])
                tmparray[:] = np.NaN
                DPD = np.concatenate((DPD,tmparray))
                del tmparray
                
            # Average Wave Period
            if any(tt == 'APD' for tt in header4):
                tind = header4.index('APD')
                APD = np.concatenate((APD,tmpdata[:,tind]))
            else:
                tmparray = np.ones([tmpdata.shape[0],1])
                tmparray[:] = np.NaN
                APD = np.concatenate((APD,tmparray))
                del tmparray 
                            
            # Mean Wave Direction
            if any(tt == 'MWD' for tt in header4):
                tind = header4.index('MWD')
                MWD = np.concatenate((MWD,tmpdata[:,tind]))
            else:
                tmparray = np.ones([tmpdata.shape[0],1])
                tmparray[:] = np.NaN
                MWD = np.concatenate((MWD,tmparray))
                del tmparray 
    
            # Sea Level Pressure
            if any(tt == 'PRES' for tt in header4):
                tind = header4.index('PRES')
                PRES = np.concatenate((PRES,tmpdata[:,tind]))
            elif any(tt == 'BAR' for tt in header4):
                tind = header4.index('BAR')
                PRES = np.concatenate((PRES,tmpdata[:,tind]))
            else:
                tmparray = np.ones([tmpdata.shape[0],1])
                tmparray[:] = np.NaN
                PRES = np.concatenate((PRES,tmparray))
                del tmparray 
    
            # Air temperature
            if any(tt == 'ATMP' for tt in header4):
                tind = header4.index('ATMP')
                ATMP = np.concatenate((ATMP,tmpdata[:,tind]))
            else:
                tmparray = np.ones([tmpdata.shape[0],1])
                tmparray[:] = np.NaN
                ATMP = np.concatenate((ATMP,tmparray))
                del tmparray  
    
            # Sea surface temperature
            if any(tt == 'WTMP' for tt in header4):
                tind = header4.index('WTMP')
                WTMP = np.concatenate((WTMP,tmpdata[:,tind]))
            else:
                tmparray = np.ones([tmpdata.shape[0],1])
                tmparray[:] = np.NaN
                WTMP = np.concatenate((WTMP,tmparray))
                del tmparray                   
      
      
            # ================== Time management =================
            # Years
            if any(tt == 'YY' for tt in header4):
                tind = header4.index('YY')
                years = tmpdata[:,tind] + 1900   
            elif any(tt == '#YY' for tt in header4):
                tind = header4.index('#YY')
                years = tmpdata[:,tind]             
            else:
                tind = header4.index('YYYY')
                years = tmpdata[:,tind]
    
            # Minutes
            if any(tt == 'mm' for tt in header4):
                tind = header4.index('mm')
                mm = tmpdata[:,tind]
            else:
                mm = np.zeros([tmpdata.shape[0],1])
    
    
            # Create time vector
            # Seconds from 1900-01-01 will be used
            for aa in range(tmpdata.shape[0]):
                tmptime = datetime.datetime(int(years[aa]),
                                            int(tmpdata[aa,1]),
                                            int(tmpdata[aa,2]),
                                            int(tmpdata[aa,3]),
                                            int(mm[aa]),
                                            0)
                sectime = tmptime - datetime.datetime(1900,1,1,0,0,0)
                secsecs = sectime.total_seconds()
                wavetime = np.concatenate([wavetime,np.array([secsecs])])            
                del tmptime,sectime,secsecs
            
            # ================ Clean up ================
            del header1,header2,header3,header4,years,mm
    
    
    # Remove first term of each array
    wavetime = np.delete(wavetime,0)
    WDIR = np.delete(WDIR,0)
    WSPD = np.delete(WSPD,0)
    GST = np.delete(GST,0)
    WVHT = np.delete(WVHT,0)
    DPD = np.delete(DPD,0)
    APD = np.delete(APD,0)
    MWD = np.delete(MWD,0)
    PRES = np.delete(PRES,0)
    ATMP = np.delete(ATMP,0)
    WTMP = np.delete(WTMP,0)
    
    # Clean up variables
    WDIR[WDIR==999] = np.NaN
    WSPD[WSPD==99] = np.NaN
    GST[GST==99] = np.NaN
    WVHT[WVHT==99] = np.NaN
    DPD[DPD==99] = np.NaN
    APD[APD==99] = np.NaN
    MWD[MWD==999] = np.NaN
    PRES[PRES==9999] = np.NaN
    ATMP[ATMP==99] = np.NaN
    WTMP[WTMP==99] = np.NaN
    
    # Order chronologically
    sorted_index = np.argsort(wavetime)
    wavetime = [wavetime[i] for i in sorted_index]
    WDIR = [WDIR[i] for i in sorted_index]
    WSPD = [WSPD[i] for i in sorted_index]
    GST = [GST[i] for i in sorted_index]
    WVHT = [WVHT[i] for i in sorted_index]
    DPD = [DPD[i] for i in sorted_index]
    APD = [APD[i] for i in sorted_index]
    MWD = [MWD[i] for i in sorted_index]
    PRES = [PRES[i] for i in sorted_index]
    ATMP = [ATMP[i] for i in sorted_index]
    WTMP = [WTMP[i] for i in sorted_index]
    
    
    #===========================================================================
    # Save as NetCDF 
    #===========================================================================
      
      
    # Global attributes  
    if ncformat == 4:
        print("Saving the buoy data with NetCDF4 format")
        nc = netCDF4.Dataset(buoyfld + '/' + buoyid + '.nc', 'w',
                             format='NETCDF4')
    else:
        print("Saving the buoy data with NetCDF3 format")
        nc = netCDF4.Dataset(buoyfld + '/' + buoyid + '.nc', 'w', 
                             format='NETCDF3_CLASSIC')
    nc.Description = buoyid + ' NDBC Bulk Parameter Data'
    nc.rawdata = 'National Data Buoy Center \nwww.ndbc.noaa.gov'
    nc.Author = 'ggarcia@coas.oregonstate.edu \nNearshore Modeling Group'
    nc.Created = time.ctime()
    nc.Software = 'Created with Python ' + sys.version
    nc.NetCDF_Lib = str(netCDF4.getlibversion())
    nc.Script = os.path.realpath(__file__)
    
    # Create dimensions  
    nc.createDimension('wave_time',None)
    
    
    # pyroms subroutine to write NetCDF fields
    def write_nc_var(var, name, dimensions, units=None, longname=None):
        nc.createVariable(name, 'f8', dimensions)
        if units is not None:
            nc.variables[name].units = units  
        if longname is not None:
            nc.variables[name].long_name = longname  
        nc.variables[name][:] = var
          
    
    # Write Variables To NetCDF file
    write_nc_var(wavetime,'wave_time','wave_time',
                 'seconds since 1900-01-01 00:00:00','measurement time UTC')
    write_nc_var(WDIR, 'WDIR', 'wave_time', 'degrees',
                 'Wind direction (direction the wind is coming from in ' + 
                 'degrees clockwise from true North')
    write_nc_var(WSPD, 'WSPD', 'wave_time', 'meter second-1',
                 'Wind speed averaged over an eight-minute period')
    write_nc_var(GST, 'GST','wave_time','meter second-1','Peak gust speed')
    write_nc_var(WVHT,'WVHT','wave_time','meter',
                 'Significant wave height during the 20 minute sampling period')
    write_nc_var(DPD,'DPD','wave_time','second',
                 'Dominant wave period (period with the maximum wave energy)')
    write_nc_var(APD,'APD','wave_time','second',
                 'Average wave period of all waves during the 20 minute' + 
                 ' sampling period')
    write_nc_var(MWD,'MWD','wave_time','degrees',
                 'The direction from which the waves at the dominant period ' + 
                 'are coming in degrees from true North, increasing clockwise')
    write_nc_var(PRES,'PRES','wave_time','hPa','Sea level pressure')
    write_nc_var(ATMP,'ATMP','wave_time','Celsius','Air temperature')
    write_nc_var(WTMP,'WTMP','wave_time','Celsius','Sea surface temperature')
      
    # Close NetCDF File 
    nc.close()
    


#===============================================================================
# Spectral parameters code
#===============================================================================

def spec2nc(buoyfld,dtheta=5):
    '''
    Code to convert NDBC spectral data files to netCDF format. 
    
    Usage:
    ------
    spec2nc(buoyfld,dtheta)
    
    Input:
    ------
    buoyfld  : Folder where the text files reside. Those should be the only
               files in the folder.
    dtheta   : Directional resolution for the reconstruction of the frequency-
               direction spectrum. Defaults to 5 degrees. 
    
    Notes:
      1. NetCDF4 file will be generated
      2. Code is not optimized since this is not something you will want to be
         running often. Beware of slow performance for large datasets.
    
    References:
    Kuik, A.J., G.Ph. van Vledder, and L.H. Holthuijsen, 1998: "Method for
      the Routine Analysis of Pitch-and-Roll Buoy Wave Data", Journal of
      Physical Oceanography, 18, 1020-1034.
      
    TODO:
    - Only works with newer formats (YY MM DD hh mm)
    - Fix negative energy (the fix is in the code just incorporate)
    '''
    
    # For testing only ---------------------------------------------------------
    #buoyfld = '/home/shusin2/users/ggarcia/data/wave/b46029/spec/'
    #dtheta = 20
    # --------------------------------------------------------------------------
    
    # Construct directional angle
    angles = np.arange(0.0,360.0,dtheta)
    
    # Time reference
    basetime = datetime.datetime(1900,1,1,0,0,0)

    #===========================================================================
    # Read file information
    #===========================================================================
    
    # Get all files in folder
    archivos = glob.glob(buoyfld + '/*.txt')
    archivos = [x.split('/')[-1] for x in archivos]
    
    # Year information
    years = [x.split('.')[0][-4:] for x in archivos]    # Get all year stamps
    years = list(set(years))                            # Find unique years
    years.sort()                                        # Sort years
    
    # Get buoy ID information
    buoyid = [x[0:5] for x in archivos]                 # Find buoy ids
    buoyid = list(set(buoyid))                          # Find unique ids    
    if len(buoyid)>1:
        print('This code does not support conversion for multiple buoys')
        buoyid = buoyid[0]
        print('  ' +  buoyid + ' will be processed')
    else:
        buoyid = buoyid[0]
        
    # Info
    print('Found ' + np.str(len(years)) + ' files')
    print(buoyfld)    
    
    # Create output netcdf file ------------------------------------------------
    # Global attributes  
    nc = netCDF4.Dataset(buoyfld + '/' + buoyid + '_spec.nc',
                         'w',format='NETCDF4')
    nc.Description = buoyid + ' NDBC Spectral Data'
    nc.Rawdata = 'National Data Buoy Center \nwww.ndbc.noaa.gov'
    nc.Author = getpass.getuser()
    nc.Created = time.ctime()
    nc.Software = 'Created with Python ' + sys.version
    nc.NetCDF_Lib = str(netCDF4.getlibversion())
    nc.Script = os.path.realpath(__file__)
    nc.Notes = 'Nautical convention used for directions'    
    
    # Reconstruct the spectrum -------------------------------------------------        

    # counter variable to create variables in the netcdf file
    masterCnt = 0
    cnt_freq = 0 
    cnt_dir = 0
    tstep_freq = 0
    tstep_dir = 0
    
    # This variable will change if the format changes
    formIdCnt = 0
    formId = '0'
            
    # Frequency array to find if the reported frequencies changed
    freqArray = []
    
    # Loop over years
    for aa in years:
        
        # Info
        print('  Working on year ' + aa)
        
        # Load spectral density files
        # Check if file exists
        tmpfile = buoyfld + buoyid + 'w' + aa + '.txt'
        if os.path.isfile(tmpfile) == False:
            # No spectral density found for the given year, go to next one
            continue
        
        # Info
        print('    Spectral density data')
        
        # Increase master counter variable
        masterCnt += 1
        
        # Read spectral density data (frequency spectra) and identify the time
        # information given
        f_w = open(tmpfile,'r')
        freq = f_w.readline().split()
               
        # Allocate the frequency array and determine if the format changed ----
        freqArray.append(freq)
        if masterCnt > 1:
            if freq != freqArray[masterCnt-2]:
                
                # Reset counter variables                
                cnt_freq = 0 
                cnt_dir = 0
                tstep_freq = 0
                tstep_dir = 0
                
                # Update form Id counter
                formIdCnt += 1
                
                if formIdCnt > 9:
                    print('\n10 Different Formats Found')
                    print('Check your data, quitting ...')
                    nc.close()
                    sys.exit()
                    
                # Update form Id text
                formId = '%01.0f' % formIdCnt
                
                # Message to user
                print('    Different format found')
                print('      New variables introduced')
                
               
        # Find if minutes are given
        if freq[4] == 'mm':
            freqInd0 = 5
        else:
            freqInd0 = 4
               
        # Read frequencies
        freq = np.array(freq[freqInd0:],dtype=float)
        f_w.close()
        
        # Load spectral density
        freq_spec = np.loadtxt(tmpfile,skiprows=1)
        
        # Allocate time and spectral density data    
        freq_time = np.zeros((freq_spec.shape[0]))  
        for bb in range(freq_time.shape[0]):
            tmpYear  = np.int(freq_spec[bb,0])
            tmpMonth = np.int(freq_spec[bb,1])
            tmpDay   = np.int(freq_spec[bb,2])
            tmpHour  = np.int(freq_spec[bb,3])
            if freqInd0 == 4:
                tmpMin = np.int(0)
            else:
                tmpMin = np.int(freq_spec[bb,4])
            
            if tmpYear < 100:
                tmpYear = tmpYear + 1900
            
            freq_time[bb] = (datetime.datetime(tmpYear,tmpMonth,tmpDay,
                                               tmpHour,tmpMin) -
                             basetime).total_seconds()

        freq_spec = freq_spec[:,freqInd0:]
        
        # No Data Filter (NDBC uses 999.00 when there is no data)
        goodDataInd = freq_spec[:,1] < 990.00
        freq_time = freq_time[goodDataInd]
        freq_spec = freq_spec[goodDataInd,:]
        
        # Create frequency spectra variables
        cnt_freq += 1
        if cnt_freq == 1:
            
            # Create dimensions  (NetCDF4 supports multiple unlimited dimensions)
            nc.createDimension('wave_time'+formId,None)
        
            # Create bulk parameter variables
            nc.createVariable('Hsig'+formId,'f8','wave_time'+formId)
            nc.variables['Hsig'+formId].units = 'meter'
            nc.variables['Hsig'+formId].long_name = 'Significant wave height'
            
            # Create frequency dimension
            nc.createDimension('freq'+formId,freq.shape[0])
            
            nc.createVariable('wave_time'+formId,'f8','wave_time'+formId)
            nc.variables['wave_time'+formId].units = \
            "seconds since 1900-01-01 00:00:00"
            nc.variables['wave_time'+formId].calendar = "julian"
            
            nc.createVariable('freq_spec'+formId,'f8',
                              ('wave_time'+formId,'freq'+formId))
            nc.variables['freq_spec'+formId].units = 'meter2 second'
            nc.variables['freq_spec'+formId].long_name = 'Frequency variance spectrum'            
            
            nc.createVariable('frequency'+formId,'f8',('freq'+formId))
            nc.variables['frequency'+formId].units = 'Hz'
            nc.variables['frequency'+formId].long_name = 'Spectral frequency'
            nc.variables['frequency'+formId][:] = freq
        
    
        # Information
        print('    Computing Bulk Parameters')
        
        # Compute bulk parameters
        moment0 = np.trapz(freq_spec.T,freq,axis=0)
        Hsig = 4.004*(moment0)**0.5
        
        # Write to NetCDF file
        if cnt_freq == 1:
            nc.variables['Hsig'+formId][:] = Hsig
            nc.variables['freq_spec'+formId][:] = freq_spec
            nc.variables['wave_time'+formId][:] = freq_time
        else:
            nc.variables['Hsig'+formId][tstep_freq:] = Hsig
            nc.variables['freq_spec'+formId][tstep_freq:,:] = freq_spec
            nc.variables['wave_time'+formId][tstep_freq:] = freq_time

                
        # Check if directional data exists -------------------------------------
        tmp_alpha_1 = buoyfld  + buoyid + 'd' + aa + '.txt'
        tmp_alpha_2 = buoyfld  + buoyid + 'i' + aa + '.txt'
        tmp_r_1 = buoyfld  + buoyid + 'j' + aa + '.txt'
        tmp_r_2 = buoyfld  + buoyid + 'k' + aa + '.txt'
    
    
        if (os.path.isfile(tmp_alpha_1) and os.path.isfile(tmp_alpha_2) and
            os.path.isfile(tmp_r_1) and os.path.isfile(tmp_r_2)):
            
            # Information
            print('    Directional Data')
            
            # Read frequency of the directional spectra (not always agree with
            # the spectral densities)
            f_w2 = open(tmp_alpha_1,'r')
            freqDirSpec = f_w2.readline().split()
            freqDirSpec = np.array(freqDirSpec[freqInd0:],dtype=float)
            f_w2.close()  
            
            # Create directional spectra variables
            cnt_dir += 1
            if cnt_dir == 1:
                nc.createDimension('dir_time'+formId,None)
                nc.createDimension('dir'+formId,angles.shape[0])
                
                # Create frequency dimension
                nc.createDimension('freqDir'+formId,freqDirSpec.shape[0])
            
                nc.createVariable('dir_time'+formId,'f8','dir_time'+formId)
                nc.variables['dir_time'+formId].units = \
                "seconds since 1900-01-01 00:00:00"
                nc.variables['dir_time'+formId].calendar = "julian"
            
                nc.createVariable('dir_spec'+formId,'f8',
                                  ('dir_time'+formId,'freqDir'+formId,
                                   'dir'+formId))
                nc.variables['dir_spec'+formId].units = 'meter2 second degree-1'
                nc.variables['dir_spec'+formId].long_name = \
                    'Frequency-Direction variance spectrum'  
                    
                nc.createVariable('direction'+formId,'f8',('dir'+formId))
                nc.variables['direction'+formId].units = 'degree'
                nc.variables['direction'+formId].long_name = \
                    'Degrees from true north in oceanographic convention'
                nc.variables['direction'+formId][:] = angles
                
                nc.createVariable('frequencyDir'+formId,'f8',('freqDir'+formId))
                nc.variables['frequencyDir'+formId].units = 'Hz'
                nc.variables['frequencyDir'+formId].long_name = 'Spectral frequency for dir_spec'
                nc.variables['frequencyDir'+formId][:] = freqDirSpec
            
            
            # Read spectral data
            alpha_1 = np.loadtxt(tmp_alpha_1,skiprows=1)
            alpha_2 = np.loadtxt(tmp_alpha_1,skiprows=1)
            r_1 = np.loadtxt(tmp_alpha_1,skiprows=1) * 0.01
            r_2 = np.loadtxt(tmp_alpha_1,skiprows=1) * 0.01
    
            # Allocate date
            dir_time = np.zeros((alpha_1.shape[0]))              
            for bb in range(dir_time.shape[0]):
                tmpYear  = np.int(alpha_1[bb,0])
                tmpMonth = np.int(alpha_1[bb,1])
                tmpDay   = np.int(alpha_1[bb,2])
                tmpHour  = np.int(alpha_1[bb,3])
                if freqInd0 == 4:
                    tmpMin = np.int(0)
                else:
                    tmpMin = np.int(alpha_1[bb,4])
                                
                if tmpYear < 100:
                    tmpYear = tmpYear + 1900                                
                
                dir_time[bb] = (datetime.datetime(tmpYear,tmpMonth,tmpDay,
                                                  tmpHour,tmpMin) - 
                                basetime).total_seconds()

                             
            # Read data
            alpha_1 = alpha_1[:,freqInd0:]
            alpha_2 = alpha_2[:,freqInd0:]
            r_1 = r_1[:,freqInd0:]
            r_2 = r_2[:,freqInd0:]
            
            # No Data Filter (NDBC uses 999.00 when there is no data)
            goodDataInd = np.logical_and(alpha_1[:,1] != 999.00,
                                         alpha_1[:,2] != 999.00)
            alpha_1  = alpha_1[goodDataInd,:]
            alpha_2  = alpha_2[goodDataInd,:]
            r_1      = r_1[goodDataInd,:]
            r_2      = r_2[goodDataInd,:]
            dir_time = dir_time[goodDataInd]
            
            
            # Find where dir_time and freq_time match and compute those values
            # only
            repInd    = np.in1d(dir_time,freq_time)
            alpha_1   = alpha_1[repInd]
            alpha_2   = alpha_2[repInd]
            r_1       = r_1[repInd]
            r_2       = r_2[repInd]
            dir_time  = dir_time[repInd]
            
            repInd    = np.in1d(freq_time,dir_time)
            freq_spec = freq_spec[repInd] 

            # Interpolate density spectrum into directional bins            
            if not np.array_equal(freq,freqDirSpec):
                freqSpecAll = np.copy(freq_spec)
                freq_spec   = np.zeros_like((r_1)) * np.NAN
                for bb in range(freq_spec.shape[0]):
                    freq_spec[bb,:] = np.interp(freqDirSpec,freq,
                                                freqSpecAll[bb,:])
            
            # Construct 2D spectra
            # See http://www.ndbc.noaa.gov/measdes.shtml
            wspec = np.NaN * np.zeros((alpha_1.shape[0],
                                       alpha_1.shape[1],angles.shape[0]))
                        
            # Time loop
            for bb in range(wspec.shape[0]):
                # Frequency loop  
                for cc in range(wspec.shape[1]):
                    # Direction loop
                    for dd in range(wspec.shape[2]):
                        wspec[bb,cc,dd] = (freq_spec[bb,cc] * np.pi/180.0 *
                                           (1.0/np.pi) * 
                                           (0.5 + r_1[bb,cc] * 
                                            np.cos((angles[dd]-alpha_1[bb,cc])*
                                                   np.pi/180.0) +
                                            r_2[bb,cc] * 
                                            np.cos(2 * np.pi / 180.0 * 
                                                   (angles[dd]-alpha_2[bb,cc])))
                                           )
#             CLEAN ME UP---------------
#             # Get the positive energy
#             tmpSpec = spec.copy()
#             tmpSpec[tmpSpec<0] = 0.0
#             tmpFreqSpec = np.sum(tmpSpec,axis=-1)*(dirs[2]-dirs[1])
#             posEnergy = np.trapz(np.abs(tmpFreqSpec),freq,axis=-1)
#             
#             # Get the total energy
#             tmpFreqSpec = np.sum(spec,axis=-1)*(dirs[2]-dirs[1])
#             totalEnergy = np.trapz(np.abs(tmpFreqSpec),freq,axis=-1)
#             
#             # Scale the spectrum
#             spec = np.array([spec[aa,...] * totalEnergy[aa]/posEnergy[aa] 
#                              for aa in range(spec.shape[0])])
#             spec[spec<0] = 0.0

            # Write to file
            if cnt_dir == 1:
                nc.variables['dir_spec'+formId][:] = wspec
                nc.variables['dir_time'+formId][:] = dir_time
            else:
                nc.variables['dir_spec'+formId][tstep_dir:,:,:] = wspec
                nc.variables['dir_time'+formId][tstep_dir:] = dir_time
            
            tstep_dir += dir_time.shape[0]
                    
        # Update frequency time step and go to next year        
        tstep_freq += freq_time.shape[0]
        
       
    # Wrap up ------------------------------------------------------------------
    # Information
    print('Data stored as:')
    print('  ' + buoyfld + '/' + buoyid + '.nc')
    
    # Close NetCDF File     
    nc.close()
        

    