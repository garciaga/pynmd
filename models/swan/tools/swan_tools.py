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

__author__ = "Gabriel Garcia Medina"
__email__ = "ggarcia@coas.oregonstate.edu"
__group__ = 'Nearshore Modeling Group'

# Import Modules
import datetime as _datetime
import numpy as _np
import re as _re

# ==============================================================================
# Read Wind Files 
# ==============================================================================
def readWind(windFile):
    """
    Read wind files in SWAN ASCII format
    """

    # Find the length of the file and get dates
    print('Finding time steps')
    fobj = open(windFile,'r')    
    ldates = _re.findall(r'\d{8}.\d{6}',fobj.read())
    fobj.close()
    
    # Set as wavetime
    windTime = [_datetime.datetime.strptime(x,"%Y%m%d.%H%M%S") for x in ldates]
    windTime = _np.array(windTime)

    # Read wind data
    print('Reading wind data')
    uwnd = []      # Zonal Wind
    vwnd = []      # Meridional Wind
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
            tmpData = _np.array(tmpData)
            lats = tmpData.shape[0]
            
            ind = _np.int(lats/2)
            uwnd.append(tmpData[:ind,:])
            vwnd.append(tmpData[ind:,:])

            # Did we reach end of file
            if len(tmpline) == 0:
                dataFlag = False
                break
            
            if len(tmpline[0]) < 1:
                dataFlag = False
                break
            
            # Reset the container
            tmpData = []

        else:
            # Store wind in temporary array
            tmpData.append([_np.float(bb) for bb in tmpline])
    
    fobj.close()

    # Generate arrays    
    ww3 = {'ot':windTime,'uwnd':_np.array(uwnd),'vwnd':_np.array(vwnd)}

    return ww3
