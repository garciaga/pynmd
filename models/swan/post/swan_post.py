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
    ww3spec: dictionary
             Contains
               'coords'       : x,y coordiate of spectra points
               'dirs'         : spectral directions in radians from true north
               'freq'         : spectral frequencies in Hz
               'spec'         : 4D array of spectral density data
               'info'         : General information
               
    Notes:
        Tested on SWAN 41.01A only.        
        Only works for stationary output.
        
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

    
    # Change for non stationary file
    fobj.readline()

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
    
    # Exception value
    exception = float(fobj.readline().split()[0])
    


    # Loop over points
    spec = np.zeros((num_loc,num_freq,num_dir))
    
    
    for aa in range(int(num_loc)):
        
        # Read and allocate spectral data
        fobj.readline()
        tmpfactor = float(fobj.readline())        
        tmpspec = np.zeros((num_freq,num_dir))        
        
        for bb in range(int(num_freq)):
            tmpline = fobj.readline().split()
            tmpline = [float(x) for x in tmpline]
            tmpspec[bb,:] = np.asarray(tmpline)

        # Remove missing data
        tmpspec[tmpspec==exception] = 0.0
        
        # Scale and allocate the spectrum
        spec[aa,:,:] = tmpfactor * tmpspec[:,sortind]
                
    
    # Close the text file
    fobj.close() 
    

    # Return values
    return {'spec':spec, 'freq':freq, 'dirs':dirs,
            'coords':coords,'info':info}
    



