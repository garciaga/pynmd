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
    