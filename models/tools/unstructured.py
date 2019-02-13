"""
Basic tools for dealing with unstructured grids
"""

from __future__ import division,print_function

import numpy as _np

def areaInt(xy,triang,z=None):
    """
    Compute the total area of the grid. If z is given it will give an 
    integrated value. 
    
    """

    # Compute the area of the elements
    area = _np.zeros((triang.shape[0],))
    for ii,aa in enumerate(triang):
        
        # Get the triangle points
        pnt = xy[aa]

        # Get two sides of the triangle
        side = pnt[1:] - pnt[0]

        # Vertex method to compute the area of a polygon
        # area = | (x1 * y2 - y1 * x2) + ... + (xn * y1 - yn * x1)| /2
        area[ii] = _np.abs(side[0,0] * side[1,1] - side[0,1] * side[1,0]) / 2.0
    
    # Check if we need to compute an area integral
    if z is not None:
        vol = 0.0
        for ii,aa in enumerate(triang):

            # Compute element mean and multiply by area
            vol += (area[ii] * _np.mean(z[aa]))

        return area,area.sum(),vol
    else:
        return area,area.sum()
    
# ==============================================================================
# Read Fort 14 Files
# ==============================================================================
def read_fort14(fort14):
    """ 
    Script to read fort14 files

    PARAMETERS:
    -----------
    fort14: Path to fort14 file

    RETURNS:
    --------
    Dictionary containing
    x     : x coordinates of nodes (ordered)
    y     : y coordinates of nodes (ordered)
    z     : Elevation at nodes
    triang: Triangulation 
    nbdv  : Indices for open boundary nodes
    neta  : number of open boundary nodes
    nvel  : Indices for land boundaries
    nvel  : number of land boundary nodes

    NOTES:
    ------
    1. Zero counting convention used.
    2. Nodes are returned in order so that you can combine with triangulation

    """

    fobj = open(fort14,'r')
    fobj.readline() # Header line
    tmpline = fobj.readline().split()
    nele = _np.int(tmpline[0])
    npt = _np.int(tmpline[1])

    x = _np.zeros((npt,))
    y = _np.zeros_like(x)
    z  = _np.zeros_like(x)

    # Read the grid (just in case)
    for aa in range(npt):
        tmpline = fobj.readline().split()
        x[aa] = _np.float(tmpline[1])
        y[aa] = _np.float(tmpline[2])
        z[aa]  = _np.float(tmpline[3])

    # Read triangles
    triang = _np.zeros((nele,3),dtype=_np.int)
    for aa in range(nele):
        tmpline = fobj.readline().split()
        triang[aa,0] = _np.int(tmpline[2])
        triang[aa,1] = _np.int(tmpline[3])
        triang[aa,2] = _np.int(tmpline[4])

    # Remember zero counting in python
    triang -= 1

    # Open boundary conditions -------------------------------------------------
    # Number of open boundaries
    tmpline = fobj.readline().split()
    nope = _np.int(tmpline[0])
    # Number of open boundary nodes
    tmpline = fobj.readline().split()
    neta = _np.int(tmpline[0])
    # Read open boundary nodes - loop over different open boundaries
    nbdv = []
    for aa in range(nope):
        tmpline = fobj.readline().split()
        tmpnodes = _np.int(tmpline[0])
        tmpbc = _np.zeros((tmpnodes,),dtype=_np.int)
        for bb in range(tmpnodes):
            tmpline = fobj.readline().split()
            # Allocate and automatically conver to python counting convention
            tmpbc[bb] = _np.int(tmpline[0]) - 1
        nbdv.append(tmpbc)

    # Land boundaries ----------------------------------------------------------
    # Number of land boundaries
    tmpline = fobj.readline().split()
    nbou = _np.int(tmpline[0])
    # Number of land boundary nodes
    tmpline = fobj.readline().split()
    nvel = _np.int(tmpline[0])
    # Read open boundary nodes - loop over different open boundaries
    nbvv = []
    for aa in range(nbou):
        tmpline = fobj.readline().split()
        tmpnodes = _np.int(tmpline[0])
        tmpbc = _np.zeros((tmpnodes,),dtype=_np.int)
        for bb in range(tmpnodes):
            tmpline = fobj.readline().split()
            # Allocate and automatically conver to python counting convention
            tmpbc[bb] = _np.int(tmpline[0]) - 1
        nbvv.append(tmpbc)    

    # All done here
    fobj.close()

    # Return the results
    return {'x':x,'y':y,'z':z,'triang':triang,
            'nbdv':nbdv,'neta':neta,'nope':nope,
            'nbvv':nbvv,'nvel':nvel}
