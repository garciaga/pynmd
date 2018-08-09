"""
General plotting tools meant for interacting with matplotlib
"""

from __future__ import division,print_function

import numpy as _np
import matplotlib.pyplot as _plt
from matplotlib.dates import date2num as _date2num


def align_yticks(ax1,ax2,color1=True,color2=True):
    """
    Make ax2 have the same ticks as ax1
    
    PARAMETERS:
    -----------
    ax1    : Handle of the primary axis
    ax2    : Handle of the secondary axis (the one that will be manipulated)
    color# : If true makes ytick and ylabel take the color of the last plotted
             line.
             
    RETURNS:
    --------
    ax1,ax2
    
    TODO:
    -----
    align_xticks
    
    """
    
    # Get the data from the axes
    # Primary axis
    lims1  = ax1.get_ylim()
    ticks1 = ax1.get_yticks()
    
    # Secondary axis
    lims2  = ax2.get_ylim()
    ticks2 = ax2.get_yticks()
    
    # Non-dimensionalize the primary ticks
    oldTicks = (ticks1 - lims1[0])/(lims1[1] - lims1[0])
    
    # Set the new ticks
    newTicks = oldTicks * (lims2[1] - lims2[0]) + lims2[0]
    ax2.set_yticks(newTicks)
    
    # Set the limits again
    # This is relevant if there are ticks outside of the plotting window
    # I.E. newTicks < 0 or newTicks > 1
    ax2.set_ylim(lims2)
    
    # The new grid is now redundant
    ax2.grid('off')

    # Color the ticks and label based on the last plot color
    if color1:

        # get the last color plotted
        tmpHandle = ax1.get_lines()[-1]
        ax1.tick_params('y',colors=tmpHandle.get_color())
        ax1.set_ylabel(ax1.get_ylabel(),color=tmpHandle.get_color())
    
    # Color the ticks and label based on the last plot color
    if color2:
        
        # get the last color plotted
        tmpHandle = ax2.get_lines()[-1]
        ax2.tick_params('y',colors=tmpHandle.get_color())
        ax2.set_ylabel(ax2.get_ylabel(),color=tmpHandle.get_color())
        
    return ax1,ax2
    

def stick_plot(time, u, v, **kw):
    """
    Create a stick plot

    PARAMETERS:
    -----------
    time  - Datetime vector
    u     - zonal component of wind or whatever you are plotting as numpy array
    v     - same as u but for meridional component
    
    Keyword arguments:
    ------------------
    'width'          - vector width (Default = 0.002)
    'headwidth'      - (Default = 0)
    'headlength'     - (Default = 0)
    'headaxislength' - (Default = 0)
    'ax'             - Axis to plot otherwise a new one will be created
    'ref'            - Reference vector length for the legend (Default = 1)
    'units'          - Label text (Default = r"$m s^{-1}$")

    RETURNS:
    --------
    q     - quiver handle
    qk    - quiver label handle
    ax    - axis handle
    
    Notes:
    ------
    1. The original script was developed by Filipe Fernandes and can be found:
       https://ocefpaf.github.io/python4oceanographers/blog/2014/09/15/stick_plot/
    
    """

    # Read keyword arguments
    width = kw.pop('width', 0.002)
    headwidth = kw.pop('headwidth', 0)
    headlength = kw.pop('headlength', 0)
    headaxislength = kw.pop('headaxislength', 0)
    angles = kw.pop('angles', 'uv')
    ax = kw.pop('ax', None)
    ref = kw.pop('ref',1)
    units = kw.pop('units',r"$m s^{-1}$")
    
    if angles != 'uv':
        raise AssertionError("Stickplot angles must be 'uv' so that"
                             "if *U*==*V* the angle of the arrow on"
                             "the plot is 45 degrees CCW from the *x*-axis.")

    time, u, v = map(_np.asanyarray, (time, u, v))
    if not ax:
        fig, ax = _plt.subplots()
    
    q = ax.quiver(_date2num(time), [[0]*len(time)], u, v,
                  angles='uv', width=width, headwidth=headwidth,
                  headlength=headlength, headaxislength=headaxislength,
                  **kw)

    ax.axes.get_yaxis().set_visible(False)
    ax.xaxis_date()

    qk = ax.quiverkey(q, 0.1, 0.85, ref,
                      _np.str(ref) + ' ' + units,
                    labelpos='N', coordinates='axes')    
    
    return q,qk,ax
