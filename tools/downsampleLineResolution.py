"""
Function to find points along a line at a desired alongline resolution
"""

from __future__ import division,print_function

import numpy as _np

# ==============================================================================
# Downsample to model resolution
# ==============================================================================
def downsampleLineResolution(x,y,res):

    # Allocate the first point
    xds = []
    yds = []
    xds.append(x[0])
    yds.append(y[0])

    # Initialize distance counter
    dist = 0.0

    # Loop over points
    for aa in range(1,x.shape[0]):
        
        # Get distance between points
        dx = x[aa] - x[aa-1]
        dy = y[aa] - y[aa-1]
        tmpDist = (dx**2 + dy**2)**0.5
        
        # Check if distance along line is less than model resolution
        if (dist + tmpDist) < res:
            dist += tmpDist
            continue          

        # If we are here it means that we need to interpolate
        # Find the angle
        theta = _np.arctan2(dy,dx)

        # Find how far we need to interpolate
        intDist = res - dist

        # Compute the new coordinates
        xds.append(_np.cos(theta)*intDist + x[aa-1])
        yds.append(_np.sin(theta)*intDist + y[aa-1])
        
        # Reset distance counter based on the interpolated distance
        dist = tmpDist - intDist

        # Original data is too coarse. I need to reset the distance counter and 
        # start from the freshly interpolated coordinates.
        if dist > res:
            # Figure out how many more points are needed
            numPnt = _np.int(_np.floor(dist/res))
            # Interpolate a few more points
            for bb in range(numPnt):
                xds.append(_np.cos(theta)*res + xds[-1])
                yds.append(_np.sin(theta)*res + yds[-1])
            # Calculate the final distance
            dx = x[aa] - xds[-1]
            dy = y[aa] - yds[-1]
            dist = (dx**2 + dy**2)**0.5
            
    # Allocate in array
    xds = _np.array(xds)
    yds = _np.array(yds)

    return xds,yds
