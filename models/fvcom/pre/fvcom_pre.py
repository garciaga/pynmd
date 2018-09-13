"""
Functions to process and create raw fvcom input files
"""

from __future__ import division,print_function

import numpy as _np

def read_grid(grdFile):
    """
    Reads Casename_grd.dat file

    PARAMETERS:
    -----------
    grdFile: Full path to grid file

    RETURNS:
    --------
    Dictionary containing
    node  : Number of nodes
    nele  : Number of elements
    x     : x coordinates
    y     : y coordinates
    h     : depth 
    triang: triangulation

    """

    # Open the grid file
    fobj = open(grdFile,'r')

    # Read the elements -----------------
    ele = []
    ele.append([int(aa) for aa in fobj.readline().split('\t')[:-1]])
    de = 1 # Delta element
    while de > 0:
        tmpLine = fobj.readline().split('\t')
        de = int(tmpLine[0]) - ele[-1][0]        
        if de > 0:
            tmpLine = [int(aa) for aa in tmpLine]
            ele.append(tmpLine[:-1])
        else:
            break

    # Convert to array
    ele = _np.array(ele)
    
    # Number of elements (useful value to have sometimes)
    nele = ele[-1,0]

    # Preallocate the coordinates
    node = _np.max(ele[:,1:])
    x = _np.zeros((node))
    y = _np.zeros_like(x)
    h = _np.zeros_like(x)

    # Allocate the first coordinates
    x[0] = _np.float(tmpLine[1])
    y[0] = _np.float(tmpLine[2])
    h[0] = _np.float(tmpLine[3])

    # Read the coordinates
    ii = 0
    while True:
        ii += 1
        tmpLine = fobj.readline().rstrip().split('\t')
        if len(tmpLine) < 2:
            break        
        x[ii] = _np.float(tmpLine[1])
        y[ii] = _np.float(tmpLine[2])
        h[ii] = _np.float(tmpLine[3])
    
    # All done here
    fobj.close()

    # Remember zero counting in python    
    ele -= 1

    return {'node':node,'x':x,'y':y,'h':h,'nele':nele,'triang':ele}
