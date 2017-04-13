"""
General plotting tools meant for interacting with matplotlib
"""

from __future__ import division,print_function

# import numpy as _np

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
    