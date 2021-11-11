"""
Functions to process and create raw fvcom input files
"""

import numpy as _np
import os as _os

# ==============================================================================
# Read grid file
# ==============================================================================
def read_grid_2(grdFile):
    """
    Reads Casename_grd.dat file compatible with FVCOM 2

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


def read_grid(grdFile):
    """
    Reads Casename_grd.dat file compatible with FVCOM 3+

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
    triang: triangulation

    """

    # Open the grid file
    fobj = open(grdFile,'r')

    # Read the elements -----------------
    nnode = _np.int(fobj.readline().rstrip().split('=')[-1])
    nele = _np.int(fobj.readline().rstrip().split('=')[-1])

    ele = _np.zeros((nele,3),dtype=_np.int)
    for ii in range(ele.shape[0]):
        tmpLine = fobj.readline().rstrip().split()
        ele[ii,0] = _np.int(tmpLine[1])
        ele[ii,1] = _np.int(tmpLine[2])
        ele[ii,2] = _np.int(tmpLine[3])

    # Zero counting
    ele -= 1

    # Read the node coordinates
    nodeId = _np.zeros((nnode),dtype=_np.int)
    x = _np.zeros((nnode))
    y = _np.zeros((nnode))

    for ii in range(x.shape[0]):
        tmpLine = fobj.readline().rstrip().split()
        nodeId[ii] = _np.int(tmpLine[0])
        x[ii] = _np.float(tmpLine[1])
        y[ii] = _np.float(tmpLine[2])

    # Zero counting
    nodeId -= 1

    # All done here
    fobj.close()

    # Output stuff
    return {'nodeId':nodeId,'x':x,'y':y,
            'nele':nele,'triang':ele}

def read_dep(depFile):
    """
    Reads Casename_dep.dat file compatible with FVCOM 2

    PARAMETERS:
    -----------
    grdFile: Full path to grid file

    RETURNS:
    --------
    Dictionary containing
    x     : x coordinates
    y     : y coordinates
    h     : depth 

    """

    # Open the grid file
    fobj = open(depFile,'r')

    # Number of nodes
    nnode = _np.int(fobj.readline().rstrip().split('=')[-1])

    # Read the node coordinates
    x = _np.zeros((nnode))
    y = _np.zeros((nnode))
    h = _np.zeros((nnode))

    for ii in range(x.shape[0]):
        tmpLine = fobj.readline().rstrip().split()
        x[ii] = _np.float(tmpLine[0])
        y[ii] = _np.float(tmpLine[1])
        h[ii] = _np.float(tmpLine[2])
   
    # All done here
    fobj.close()

    # Output stuff
    return {'x':x,'y':y,'h':h}


# ==============================================================================
# Write grid file
# ==============================================================================
def write_grid(grd,outFld,casename):
    """
    Writes Casename_grd.dat file

    PARAMETERS:
    -----------
    grd: Dictionary containing
         x     : x or longitude of nodes
         y     : y or latitude of nodes
         triang: Elements (use zero counting)
    outFld : Path to directory where file will be written
    casename: String with name of file

    NOTES:
    ------
    1. Zero counting convention used.

    """

    # Create the output file
    outFile = outFld + '/' + casename + '_grd.dat'
    print('Creating: ' + outFile)
    
    fid = open(outFile,'w')

    # Write header lines
    fid.write('Node Number = {:.0f}\n'.format(grd['x'].shape[0]))
    fid.write('Cell Number = {:.0f}\n'.format(grd['triang'].shape[0]))

    # Write triangulation
    for ii in range(grd['triang'].shape[0]):
        fid.write('{:.0f} {:.0f} {:.0f} {:.0f} 1\n'.format(ii+1,
                        grd['triang'][ii,0] + 1,
                        grd['triang'][ii,1] + 1,
                        grd['triang'][ii,2] + 1))
    
    # Write nodes
    for ii in range(grd['x'].shape[0]):
        fid.write('{:.0f} {:.8f} {:.8f}\n'.format(ii + 1,
                                                  grd['x'][ii],
                                                  grd['y'][ii]))
    
    # Close now
    fid.close()

# ==============================================================================
# Write depth file
# ==============================================================================
def write_dep(grd,outFld,casename):
    """
    Writes Casename_dep.dat file

    PARAMETERS:
    -----------
    grd: Dictionary containing
         x     : x or longitude of nodes
         y     : y or latitude of nodes
         d     : depth [meter]
    outFld : Path to directory where file will be written
    casename: String with name of file

    NOTES:
    ------
    1. Depth in FVCOM is specified at the node point

    """

    # Create the output file
    outFile = outFld + '/' + casename + '_dep.dat'
    print('Creating: ' + outFile)
    
    fid = open(outFile,'w')

    # Write header lines
    fid.write('Node Number = {:.0f}\n'.format(grd['x'].shape[0]))

    # Write nodes
    for ii in range(grd['x'].shape[0]):
        fid.write('{:.8f} {:.8f} {:.3f}\n'.format(grd['x'][ii],
                                                  grd['y'][ii],
                                                  grd['z'][ii]))
    
    # Close now
    fid.close()

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
    nodeId  : Node ID of the rivers (you'll need to subtract 1 for python)
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
        # Check if there are tabs in the file
        tmpLine = fobj.readline().rstrip().split()
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

# ==============================================================================
# Read meteorological forcing
# ==============================================================================
def read_mc(meteoFile):
    """
    Reads Casename_mc.dat file

    PARAMETERS:
    -----------
    meteoFile: Full path to the file

    RETURNS:
    --------
    Dictionary containing
    time   = Time in hours measured relative to start time of model
    qprec  = Amount of precipitation (m/year)
    qevap  = Amount of evaopration (m/year)
    wds    = Wind speed (m/s)
    wdd    = Wind direction from which wind blows clockwise from north (deg)
    hflux  = Net surface heat flux (W/m2)
    hshort = Shortwave radiation flux at the surface (W/m2)

    """

    # Find out how long the file is
    nLine = 0
    fobj = open(meteoFile,'r')
    line = fobj.readline()
    while line:
        nLine += 1
        line = fobj.readline()        
    fobj.close()

    # Number of entries in the file
    nEnt = _np.int((nLine-1)/2)

    # Open the grid file
    fobj = open(meteoFile,'r')

    # Ignore the first line for now
    fobj.readline()

    # Preallocate variables
    ot     = _np.zeros((nEnt,)) # Time in hours measured relative to start time
    qprec  = _np.zeros_like(ot) # Amount of precipitation m/year
    qevap  = _np.zeros_like(ot) # Amount of evaopration m/year
    wds    = _np.zeros_like(ot) # Wind speed (m/s)
    wdd    = _np.zeros_like(ot) # Wind directio in degrees from which wind blows
                                # Clockwise from north
    hflux  = _np.zeros_like(ot) # Net surface heat flux (W/m2)
    hshort = _np.zeros_like(ot) # Shortwave radiation flux at the surface (W/m2)

    for aa in range(nEnt):
        # Read time
        tmpLine = fobj.readline()
        if tmpLine == '':
            break
        ot[aa] = _np.float(tmpLine.rstrip())

        # Read variables
        tmpLine = fobj.readline()
        if tmpLine == '':
            break
        tmpLine    = tmpLine.rstrip().split()
        qprec[aa]  = _np.float(tmpLine[0])
        qevap[aa]  = _np.float(tmpLine[1])
        wds[aa]    = _np.float(tmpLine[2])
        wdd[aa]    = _np.float(tmpLine[3])
        hflux[aa]  = _np.float(tmpLine[4])
        hshort[aa] = _np.float(tmpLine[5])

    # Close the file
    fobj.close()

    # Give back
    return {'time':ot,'qprec':qprec,'qevap':qevap,'wds':wds,'wdd':wdd,
            'hflux':hflux,'hshort':hshort}

# ==============================================================================
# Read meteorological forcing
# ==============================================================================
def read_obc(obcFile):
    """
    Reads open boundary node (Casename_obc.dat) file

    PARAMETERS:
    -----------
    meteoFile: Full path to the file

    RETURNS:
    --------
    nx2 array containing node id (zero based), node type

    """

    xy = _np.genfromtxt(obcFile,skip_header=1)[:,1:]
    
    # Convert to zero counting
    xy[:,0] -= 1
    
    # Create integers
    xy = _np.asarray(xy,dtype=_np.int)

    return xy
