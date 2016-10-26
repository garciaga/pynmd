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
    num_loc = float(tmpline[0])
    
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
    num_freq = float(tmpline[0])
    freq = np.zeros((num_freq,))
    for aa in range(int(num_freq)):
        freq[aa] = float(fobj.readline())
        
    # Directions
    tmpline = fobj.readline().split()  
    tmpstr = ' '
    info['angle_convention'] = tmpstr.join(tmpline[1:])
       
    tmpline = fobj.readline().split()
    num_dir = float(tmpline[0])
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
    



