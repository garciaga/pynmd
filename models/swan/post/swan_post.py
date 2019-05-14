"""
Tools to manage SWAN output

Authors:
-------
Gabriel Garcia Medina
    Nearshore Modeling Group
    ggarcia@coas.oregonstate.edu
Saeed Moghimi

Log of edits:
-------------
February 2015 - Created module
    Gabriel Garcia Medina

Dependencies:
-------------
    numpy
  
Internal dependencies:
----------------------
    none

"""

from __future__ import division,print_function

__author__ = "Gabriel Garcia Medina"
__email__ = "ggarcia@coas.oregonstate.edu"
__group__ = 'Nearshore Modeling Group'

# Import Modules
import netCDF4 as _netCDF4
import time as _time
import sys as _sys
import getpass as _getpass
import pynmd.tools.gtime as _gtime

import numpy as np
import datetime

#===============================================================================
# Read wave spectrum 
#===============================================================================
def read_spec(specfile):
    '''
    Read wave spectrum from a text file produced by SWAN and returns the
    spectral density and axes. Usage:
    
    read_spec(specfile)
    
    Parameters
    ----------
    specfile: string
              Full path to the spectra file to be read.
    
    Returns
    -------
    dictionary with contents
               'coords'       : x,y coordiate of spectra points
               'dirs'         : spectral directions in radians from true north
               'freq'         : spectral frequencies in Hz
               'spec'         : 4D array of spectral density data
                                (time,locations,frequencies,directions)
               'ot'           : Array of time stamps
               'info'         : General information
               
    Notes:
        Tested on SWAN 41.01A only
        
    '''

    # Find the length of the file and get dates
    fobj = open(specfile,'r')    

    # Read header
    fobj.readline()
    info = {}    
    tmpline = fobj.readline().split()
    info['version'] = tmpline[-1]
    tmpline = fobj.readline().split()
    info['Project'] = tmpline[2]   

    
    # Find if time dependent data is present
    tmpline = fobj.readline().split(' ')[0]
    if tmpline == 'TIME':
        timeDep = True
        # Do not need the extra two lines
        fobj.readline() # Time coding option (may need later)
        fobj.readline() # Locations header
    else:
        timeDep = False

    # Get number of locations
    tmpline = fobj.readline().split()
    num_loc = int(tmpline[0])
    
    # Get coordinates
    coords = np.zeros((num_loc,2))
    for aa in range(int(num_loc)):
        tmpline = fobj.readline().split()
        tmpcoord = [float(x) for x in tmpline]
        coords[aa,:] = tmpcoord[:]

        
    # Frequencies
    tmpline = fobj.readline().split()
    info['freq_type'] = tmpline[0]
    tmpline = fobj.readline().split()        
    num_freq = int(tmpline[0])
    freq = np.zeros((num_freq,))
    for aa in range(int(num_freq)):
        freq[aa] = float(fobj.readline())
        
    # Directions
    tmpline = fobj.readline().split()  
    tmpstr = ' '
    info['angle_convention'] = tmpstr.join(tmpline[1:])
       
    tmpline = fobj.readline().split()
    num_dir = int(tmpline[0])
    dirs = []
    for aa in range(int(num_dir)):        
        dirs.extend([float(fobj.readline())])
    dirs = np.asarray(dirs)    
    
    # Sort directions
    negind = dirs<0
    dirs[negind] = dirs[negind] + 360.
    sortind = np.argsort(dirs)
    dirs = dirs[sortind]
    
    # Number of quantities
    fobj.readline()
    fobj.readline()
    
    # Energy information
    fobj.readline()
    tmpline = fobj.readline().split()
    info['spec_units'] = tmpline[0]
    
    # Scale the spectrum depending on the output units
    #if info['spec_units'] == 'J/m2/Hz/degr':
    #    facun = 1.0/1025.0/9.81;
    #else:
    #    # No scaling necessary if units are m2/Hz/degr
    #    facun = 1.0;

    
    # Exception value
    exception = float(fobj.readline().split()[0])
    
    # Preallocate variables    
    ot = []
    allSpec = []
    dataFlag = True
    while dataFlag:
        
        # Time dependent lines
        if timeDep:
            
            # Read line            
            tmpline = fobj.readline().split(' ')[0]
            if len(tmpline) < 1:
                dataFlag = False
                break
    
            # Get date and time information
            ot.append(datetime.datetime.strptime(tmpline,'%Y%m%d.%H%M%S'))

        # Loop over points
        spec = np.zeros((num_loc,num_freq,num_dir))    
        
        for aa in range(int(num_loc)):
            
            # Read and allocate spectral data            
            # Factor or NODATA
            tmpline = fobj.readline().rstrip()
            if tmpline == 'NODATA':
                spec[aa,...] *= np.NAN
                continue
            elif tmpline == 'ZERO':
                spec[aa,...] *= 0.0
                continue
                
            # Get the scale factor (multiplier)
            tmpfactor = float(fobj.readline())
            
            # Preallocate the spectrum
            tmpspec = np.zeros((num_freq,num_dir))        
            
            for bb in range(int(num_freq)):
                tmpline = fobj.readline().split()
                tmpline = [float(x) for x in tmpline]
                tmpspec[bb,:] = np.asarray(tmpline)
    
            # Remove missing data
            tmpspec[tmpspec==exception] = 0.0
            
            # Scale and allocate the spectrum
            spec[aa,:,:] = tmpfactor * tmpspec[:,sortind] #* facun
                 
        # Allocate in array
        allSpec.append(spec)

        if not timeDep:
            dataFlag = False              
    
    # Close the text file
    fobj.close() 
    
    # Manage vectors
    ot = np.asarray(ot)
    allSpec = np.asarray(allSpec)

    # Return values
    return {'spec':allSpec, 'freq':freq, 'dirs':dirs,
            'coords':coords,'info':info,'ot':ot}
    

# ==============================================================================
# Convert UnSWAN stations to netCDF4
# ==============================================================================
def unSWANCompToNetcdf(swanFile,swanGrid):
    """

    Convert UnSWAN results over all the computational grid to NetCDF4

    PARAMETERS:
    -----------
    swanFile - Ascii file with bulk parameter output over all computational
               grid
    swanGrid - The fort.14 file (for now)
    
    RETURNS:
    --------
    NetCDF file with the same name as swanFile with .nc extension    

    NOTES:
    ------
    1. Only works with fort.14 files
    2. It is assumed that the first three variables in swanFile are:
       Time Xp Yp ...

    TODO:
    -----
    1. Implement conversion from stationary conditions
    2. Implement ability to read from Triangle
    
    """
   
    # Read grid file -----------------------------------------------------------
    fobj = open(swanGrid,'r')
    fobj.readline()
    tmpline = fobj.readline().split()
    nele = np.int(tmpline[0])
    npt = np.int(tmpline[1])

    xp = np.zeros((npt,))
    yp = np.zeros_like(xp)

    # Read the grid (just in case)
    for aa in range(npt):
        tmpline = fobj.readline().split()
        xp[aa] = np.float(tmpline[1])
        yp[aa] = np.float(tmpline[2])

    # Read triangles
    triang = np.zeros((nele,3),dtype=np.int)
    for aa in range(nele):
        tmpline = fobj.readline().split()
        triang[aa,0] = np.int(tmpline[2])
        triang[aa,1] = np.int(tmpline[3])
        triang[aa,2] = np.int(tmpline[4])

    # All done here
    fobj.close()

    # Remember zero counting in python
    triang -= 1

    # Read the computational grid output ---------------------------------------
    fobj = open(swanFile,'r')

    # Discard the first two lines
    fobj.readline()
    fobj.readline()

    # Get the model information
    modInfo = fobj.readline().split(':')
    modInfo = {'Run':modInfo[1].split(' ')[0],
            'Table':modInfo[2].split(' ')[0],
            'SWAN':modInfo[3][:-1]}

    # Empty line
    fobj.readline()

    # Variables
    varkeys = fobj.readline()[1:].split()
    # varkeys = fobj.readline()[1:-1]
    # varkeys = [item for item in filter(None,varkeys.split(' '))]

    # Find if the grid information is in the file
    gridFlag = 'Xp' in varkeys and 'Yp' in varkeys

    # Units
    units = fobj.readline()[1:].split()[1:]
    units[0] = '[UTC]' # Manually overwrite

    # Empty line
    fobj.readline()

    # Find all nodes in file
    tmpLine = fobj.readline().split()
    date0 = tmpLine[0]
    if gridFlag:
        xp = []
        yp = []
        xp.append(np.float64(tmpLine[1]))
        yp.append(np.float64(tmpLine[2]))
    outpnt = 1
    for line in fobj:
        if line.split()[0] == date0:
            if gridFlag:
                xp.append(np.float64(line.split()[1]))
                yp.append(np.float64(line.split()[2]))      
            outpnt += 1
        else:
            break

    if gridFlag:
        xp = np.asarray(xp)
        yp = np.asarray(yp)

    # Close the file
    fobj.close()

    # Sanity check here
    #if outpnt != npt:
    #    print("Dimensions of grid and output are different")
    #    return


    # Read File and Prepare NetCDF File ----------------------------------------
    print("Creating netCDF file")

    # Create the file and add global attributes
    outFile = swanFile + '.nc'
    nc = _netCDF4.Dataset(outFile, 'w', format='NETCDF4')
    nc.Description = 'UnSWAN Output'
    nc.Author = _getpass.getuser()
    nc.Created = _time.ctime()
    nc.Software = 'Created with Python ' + _sys.version
    nc.NetCDF_Lib = str(_netCDF4.getlibversion())
    nc.Source = swanFile
    nc.Run = modInfo['Run']
    nc.Version = modInfo['SWAN']
    nc.Table = modInfo['Table']

    # Create dimensions
    nc.createDimension('three',3)
    nc.createDimension('nele',nele) # Number of elements
    nc.createDimension('npnt',npt) # Number of nodes
    nc.createDimension('time',0) # The unlimited dimension

    # Create time vector
    nc.createVariable('ocean_time','f8',('time'))
    nc.variables['ocean_time'].units = 'seconds since 1900-01-01 00:00:00'
    nc.variables['ocean_time'].calendar = 'julian'

    nc.createVariable('matlab_time','f8',('time'))
    nc.variables['matlab_time'].units = 'Days since 0000-01-01 00:00:00'

    # Write coordinates
    nc.createVariable(varkeys[1],'f8',('npnt'))
    nc.variables[varkeys[1]].units = units[1][1:-1]
    nc.variables[varkeys[1]][:] = xp

    nc.createVariable(varkeys[2],'f8',('npnt'))
    nc.variables[varkeys[2]].units = units[2][1:-1]
    nc.variables[varkeys[2]][:] = yp

    # Write the elements
    nc.createVariable('elements','f8',('nele','three'))
    nc.variables['elements'].long_name = 'Triangulation'
    nc.variables['elements'][:] = triang

    # Create the rest of the variables
    for aa in range(3,len(varkeys)):
        nc.createVariable(varkeys[aa],'f8',('time','npnt'))
        nc.variables[varkeys[aa]].units = units[aa][1:-1]

    # ==============================================================================
    # Read and write variables to netCDF file
    # ==============================================================================

    print("  Reading and writing the file contents ...")

    # First three variables are already considered
    varkeys = varkeys[3:]

    # Load the swan file again
    fobj = open(swanFile)

    # Ignore the first seven lines
    for aa in range(7):
        fobj.readline()

    # Preallocate variable container
    #   If the grid is not crazy large this should be ok
    tmpvars = np.zeros((len(varkeys),npt))
    cnt = 0
    tstep = -1 # Time step counter (location in netcdf file)

    # Loop until the end of the file
    for line in fobj:
        
        # Increase counter variable
        cnt += 1
        
        # Temporarily allocate variables
        aa = np.asarray([np.float64(aa) for aa in line.split()])
        tmpvars[:,cnt-1] = aa[3:]

        # Write variables if all nodes have been read for the current time
        if cnt == npt:

            # Quick update
            print('  Storing ' + line.split()[0])

            # Increase time step variable
            tstep += 1

            # Time management
            timePython = datetime.datetime.strptime(line.split()[0],
                                                    "%Y%m%d.%H%M%S")
            nc.variables['matlab_time'][tstep] = _gtime.datetime_to_datenum(timePython)
            ocean_time = timePython - datetime.datetime(1900,1,1)
            nc.variables['ocean_time'][tstep] = ocean_time.total_seconds()

            # Write variables to netcdf file
            for aa in range(len(varkeys)):
                nc.variables[varkeys[aa]][tstep,...] = tmpvars[aa,:]

            # Reset counter variable
            cnt = 0

            # Reset variable container
            tmpvars *= np.NAN


    # All done here
    fobj.close()
    nc.close()


# ==============================================================================
# Convert UnSWAN stations to netCDF4
# ==============================================================================
def readTable(tableFile):
    """
    Read SWAN point output in the form of table

    PARAMETERS:
    -----------
    tableFile - Ascii file with bulk parameter output over requested points
    
    RETURNS:
    --------
    Dictionary with the variables in the file

    """
    
    # Read the computational grid output ---------------------------------------
    print('Reading and processing header lines')
    fobj = open(tableFile,'r')

    # Discard the first two lines
    fobj.readline()
    fobj.readline()

    # Get the model information
    modInfo = fobj.readline().split(':')
    swanOut = {'Run':modInfo[1].split(' ')[0],
               'Table':modInfo[2].split(' ')[0],
               'SWAN':modInfo[3][:-1]}

    # Empty line
    fobj.readline()

    # Variables
    varkeys = fobj.readline()[1:].split()

    # Find if the grid information is in the file
    gridFlag = 'Xp' in varkeys and 'Yp' in varkeys

    # Units
    units = fobj.readline()[1:].split()[1:]
    units[0] = '[UTC]' # Manually overwrite

    # Take care of empty units
    aa = -1
    while aa < len(units):
        # Counter variable increase
        aa += 1
        if len(units[aa]) == 1:
            units[aa] = '[]'
            aa += 1
            units.pop(aa)
    
    tmpOut = {}
    for aa in range(len(units)):
        tmpOut[varkeys[aa]] = units[aa]
    swanOut['units'] = tmpOut

    # Empty line
    fobj.readline()

    # Find all stations in file and store the coordinates if they are present    
    tmpLine = fobj.readline().split()
    date0 = tmpLine[0]
    if gridFlag:
        xp = []
        yp = []
        xInd = varkeys.index('Xp')
        yInd = varkeys.index('Yp')
        xp.append(np.float64(tmpLine[xInd]))
        yp.append(np.float64(tmpLine[yInd]))
    outpnt = 1
    for line in fobj:
        if line.split()[0] == date0:
            if gridFlag:
                xp.append(np.float64(line.split()[xInd]))
                yp.append(np.float64(line.split()[yInd]))      
            outpnt += 1
        else:
            break

    if gridFlag:
        xp = np.asarray(xp)
        yp = np.asarray(yp)

    # Close the file
    fobj.close()


    # ==============================================================================
    # Read and write variables to netCDF file
    # ==============================================================================

    print("  Reading the file contents ...")
    
    # Read the file
    swanVar = np.genfromtxt(tableFile,comments='%')

    # Reshape the file
    numTimes = np.int(swanVar.shape[0]/outpnt)
    swanVar = np.reshape(swanVar,(numTimes,outpnt,swanVar.shape[1]))

    # time variables
    ot = []
    for aa in range(swanVar.shape[0]):
        line = '{:15.6f}'.format(swanVar[aa,0,0])
        ot.append(datetime.datetime.strptime(line,"%Y%m%d.%H%M%S"))
    swanOut['Time'] = np.array(ot)

    # Store the coordinates
    if gridFlag:
        swanOut['Xp'] = xp
        swanOut['Yp'] = yp

    # Allocate the other variables
    for aa in varkeys:
        if aa == 'Time' or aa == 'Xp' or aa == 'Yp':
            continue
        ind = varkeys.index(aa)
        swanOut[aa] = swanVar[:,:,ind]

    return swanOut


# ==============================================================================
# Read UNSWAN computational grid
# ==============================================================================
def readUnSWANComp(swanFile,verbose=False):
    """

    Read UnSWAN computational grids

    WARNING: read only what you can afford

    PARAMETERS:
    -----------
    swanFile - Ascii file with bulk parameter output over all computational
               grid
    
    RETURNS:
    --------
    dictionary with variables in the computational grid
    
    NOTES:
    ------
    This script assumes that the first variables in the computational grid
    output file are [time,xp,yp], in that order.
    """
   
    # Read the computational grid output ---------------------------------------
    fobj = open(swanFile,'r')

    # Discard the first two lines
    fobj.readline()
    fobj.readline()

    # Get the model information
    modInfo = fobj.readline().split(':')
    modInfo = {'Run':modInfo[1].split(' ')[0],
               'Table':modInfo[2].split(' ')[0],
               'SWAN':modInfo[3][:-1]}

    # Empty line
    fobj.readline()

    # Variables
    varkeys = fobj.readline()[1:].split()

    # Find if the grid information is in the file
    gridFlag = 'Xp' in varkeys and 'Yp' in varkeys

    # Units
    units = fobj.readline()[1:].split()[1:]
    units[0] = '[UTC]' # Manually overwrite

    # Empty line
    fobj.readline()

    # Find all nodes in file
    tmpLine = fobj.readline().split()
    date0 = tmpLine[0]
    if gridFlag:
        xp = []
        yp = []
        xp.append(np.float64(tmpLine[1]))
        yp.append(np.float64(tmpLine[2]))
    npt = 1
    for line in fobj:
        if line.split()[0] == date0:
            if gridFlag:
                xp.append(np.float64(line.split()[1]))
                yp.append(np.float64(line.split()[2]))      
            npt += 1
        else:
            break

    if gridFlag:
        xp = np.asarray(xp)
        yp = np.asarray(yp)

    # Close the file
    fobj.close()

    # ==========================================================================
    # Read and write variables to netCDF file
    # ==========================================================================

    # First three variables are already considered
    if gridFlag:
        varkeys = varkeys[3:]
    else:
        # Still assuming the first entry is time
        varkeys = varkeys[1:]

    # Initialize dictionary
    swan = dict.fromkeys(varkeys)
    for aa in varkeys:
        swan[aa] = []
    swan['ocean_time'] = []
    if gridFlag:
        swan['Xp'] = xp
        swan['Yp'] = yp

    # Load the swan file again
    fobj = open(swanFile,'r')

    # Ignore the first seven lines
    for aa in range(7):
        fobj.readline()

    # Preallocate variable container
    #   If the grid is not crazy large this should be ok    
    tmpvars = np.zeros((len(varkeys),npt))
    cnt = 0
    
    # Loop until the end of the file
    for line in fobj:
        
        # Increase counter variable
        cnt += 1
        
        # Temporarily allocate variables
        aa = np.asarray([np.float64(aa) for aa in line.split()])
        tmpvars[:,cnt-1] = aa[3:]

        # Write variables if all nodes have been read for the current time
        if cnt == npt:
            
            # Quick update
            if verbose:
                print('  Reading ' + line.split()[0])

            # Time management
            timePython = datetime.datetime.strptime(line.split()[0],
                                                    "%Y%m%d.%H%M%S")
            swan['ocean_time'].append(timePython)

            # Store variables in array
            for aa in range(len(varkeys)):
                swan[varkeys[aa]].append(tmpvars[aa,:])

            # Reset counter variable
            cnt = 0

            # Reset variable container
            tmpvars = np.zeros((len(varkeys),npt))
            
    # File read
    fobj.close()

    # Make arrays
    for aa in swan.keys():
        swan[aa] = np.array(swan[aa])

    return swan
