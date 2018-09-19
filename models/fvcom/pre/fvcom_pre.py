"""
Functions to process and create raw fvcom input files
"""

from __future__ import division,print_function

import numpy as _np

# ==============================================================================
# Read grid file
# ==============================================================================
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
    nodeId = _np.zeros((node),dtype=_np.int)
    x = _np.zeros((node))
    y = _np.zeros_like(x)
    h = _np.zeros_like(x)

    # Allocate the first coordinates
    nodeId[0] = _np.int(tmpLine[0])
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
        nodeId[ii] = _np.int(tmpLine[0])
        x[ii] = _np.float(tmpLine[1])
        y[ii] = _np.float(tmpLine[2])
        h[ii] = _np.float(tmpLine[3])
    
    # All done here
    fobj.close()

    # Remember zero counting in python    
    ele -= 1

    # Output stuff
    return {'node':node,'nodeId':nodeId,'x':x,'y':y,'h':h,
            'nele':nele,'triang':ele}

# ==============================================================================
# Read river file
# ==============================================================================
def read_river(rivFile):
    """
    Reads Casename_riv.dat file

    PARAMETERS:
    -----------
    rivFile: Full path to the river file

    RETURNS:
    --------
    Dictionary containing
    nRiv    : Number of rivers
    nodeId  : Node ID of the rivers
    rivSig  : Vertical distribution of river outflow
    nTime   : Number of time steps at which river output is prescribed
    time    : Time vector of length nTime
    q       : River discharge
    temp    : River temperature
    salt    : River salinity

    """

    # Open the grid file
    fobj = open(rivFile,'r')

    # Ignore the first line for now
    fobj.readline()

    # Number of nodes with river input
    nRiv = _np.int(fobj.readline().rstrip())

    # Read the nodes where all these rivers are input
    nodeId = _np.zeros((nRiv),dtype=_np.int)
    for aa in range(nRiv):
        tmpLine = fobj.readline().rstrip().split('\t')
        nodeId[aa] = _np.int(tmpLine[0])

    # Read the vertical distribution of the river output
    rivSig = []
    for aa in range(nRiv):
        tmpLine = fobj.readline().rstrip().split('\t')
        tmpLine = [float(aa) for aa in tmpLine[1:]]
        rivSig.append(tmpLine)
    rivSig = _np.array(rivSig)

    # Number of the river discharge records
    nTime = _np.int(fobj.readline().rstrip())

    # Loop over time and rivers
    ot = _np.zeros((nTime))
    q  = _np.zeros((nTime,nRiv))
    temp = _np.zeros_like(q)
    salt = _np.zeros_like(q)

    for aa in range(nTime):
        # Read time
        ot[aa] = _np.float(fobj.readline().rstrip())
        # River flow
        tmpLine = fobj.readline().rstrip().split()
        q[aa,:] = [float(bb) for bb in tmpLine]
        # River temperature
        tmpLine = fobj.readline().rstrip().split()
        temp[aa,:] = [float(bb) for bb in tmpLine]
        # River salinity
        tmpLine = fobj.readline().rstrip().split()
        salt[aa,:] = [float(bb) for bb in tmpLine]

    # All done here
    fobj.close()

    # Output stuff
    return {'nodeId':nodeId,'nRiv':nRiv,'nTime':nTime,'vertDist':rivSig,
            'time':ot,'q':q,'temp':temp,'salt':salt}


# ==============================================================================
# Write river file
# ==============================================================================
def write_river(rivFile,rivData):
    """
    Writes Casename_riv.dat file

    PARAMETERS:
    -----------
    rivFile: Full path to the river file
    rivData: Dictionary containing 
        nRiv    : Number of rivers
        nodeId  : Node ID of the rivers
        rivSig  : Vertical distribution of river outflow
        nTime   : Number of time steps at which river output is prescribed
        time    : Time vector of length nTime
        q       : River discharge
        temp    : River temperature
        salt    : River salinity
    
    """

    # Open the grid file
    fobj = open(rivFile,'w')

    # Write the header line
    fobj.write('node  specified\n')

    # Number of nodes with river input
    try:
        nRiv = rivData['nRiv']
    except:
        nRiv = rivData['q'].shape[1]
    fobj.write(_np.str(nRiv) + '\n')

    # Write the nodes where all the rivers discharge
    for aa in range(nRiv):
        fobj.write(_np.str(rivData['nodeId'][aa]) + '\n')
    
    # Write the vertical distribution of the river output
    for aa in range(rivData['vertDist'].shape[0]):
        fobj.write(_np.str(rivData['nodeId'][aa]) + ' ')
        for bb in range(rivData['vertDist'].shape[1]):
            fobj.write('{:8.4f}'.format(rivData['vertDist'][aa,bb]))
        fobj.write('\n')

    # Number of the river discharge records
    try:
        nTime = rivData['nTime']
    except:
        nTime = rivData['q'].shape[0]
    fobj.write(_np.str(nTime) + '\n')

    # Loop over time and rivers
    for aa in range(nTime):
        
        # Write current time
        fobj.write('{:20.8f}'.format(rivData['time'][aa]) + '\n')
        
        # Discharge
        for bb in range(nRiv):
            fobj.write('{:12.4f}'.format(rivData['q'][aa,bb]))
        fobj.write('\n')
        
        # Temperature
        for bb in range(nRiv):
            fobj.write('{:12.4f}'.format(rivData['temp'][aa,bb]))
        fobj.write('\n')
        
        # Salinity
        for bb in range(nRiv):
            fobj.write('{:12.4f}'.format(rivData['salt'][aa,bb]))
        fobj.write('\n')

    # All done here
    fobj.close()

