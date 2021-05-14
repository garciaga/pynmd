"""
General plotting tools meant for interacting with matplotlib
"""

from __future__ import division,print_function

import numpy as _np
import matplotlib.pyplot as _plt
from matplotlib.dates import date2num as _date2num
import cartopy.crs as _ccrs

import pynmd.data.angles as _gangles

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


def align_xticks(ax1,ax2,color1=True,color2=True):
    """
    Make ax2 have the x axis ticks at the same location as ax1
    
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
    lims1  = ax1.get_xlim()
    ticks1 = ax1.get_xticks()
    
    # Secondary axis
    lims2  = ax2.get_xlim()
    ticks2 = ax2.get_xticks()
    
    # Non-dimensionalize the primary ticks
    oldTicks = (ticks1 - lims1[0])/(lims1[1] - lims1[0])
    
    # Set the new ticks
    newTicks = oldTicks * (lims2[1] - lims2[0]) + lims2[0]
    ax2.set_xticks(newTicks)
    
    # Set the limits again
    # This is relevant if there are ticks outside of the plotting window
    # I.E. newTicks < 0 or newTicks > 1
    ax2.set_xlim(lims2)
    
    # The new grid is now redundant
    ax2.grid('off')

    # Color the ticks and label based on the last plot color
    if color1:

        # get the last color plotted
        tmpHandle = ax1.get_lines()[-1]
        ax1.tick_params('x',colors=tmpHandle.get_color())
        ax1.set_xlabel(ax1.get_xlabel(),color=tmpHandle.get_color())
    
    # Color the ticks and label based on the last plot color
    if color2:
        
        # get the last color plotted
        tmpHandle = ax2.get_lines()[-1]
        ax2.tick_params('x',colors=tmpHandle.get_color())
        ax2.set_xlabel(ax2.get_xlabel(),color=tmpHandle.get_color())
        
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


def matchSubplotWidth(ax1,ax2):
    """
    Function to match the width of two subplots

    PARAMETERS:
    -----------
    ax1: Main axis
    ax2: Axis to be scaled
    """

    pos1 = ax1.get_position()
    pos2 = ax2.get_position()
    ax2.set_position([pos1.x0,pos2.y0,pos1.width,pos2.height])

def matchSubplotHeight(ax1,ax2):
    """
    Function to match the height of two subplots

    PARAMETERS:
    -----------
    ax1: Main axis
    ax2: Axis to be scaled
    """

    pos1 = ax1.get_position()
    pos2 = ax2.get_position()
    ax2.set_position([pos2.x0,pos1.y0,pos2.width,pos1.height])
    

def ticksNonRectangularProj(ax,parallels,meridians,proj,minLon=120,maxLon=240):
    """
    Function to draw ticklines and grid for non-rectangular projections 
    generated from Cartopy

    This function is under development, there are some things that need to be
    automated

    PARAMETERS:
    -----------
    ax        : axes handle
    parallels : Parallels to plot [list]
    meridians : Meridians to plot [list]
    proj      : Plot projection

    RETURNS:
    --------
    ax       : axes handle
    """
    
    dataProj = _ccrs.PlateCarree()

    projYticks = []
    for ii,aa in enumerate(parallels):
        tmpX = _np.arange(minLon,maxLon,0.001) # Need to find a way to automate
        tmpY = _np.ones_like(tmpX) * aa
        parTrans = proj.transform_points(dataProj,tmpX,tmpY)
        ax.plot(parTrans[:,0],parTrans[:,1],color='gray',linewidth=0.25)

        # Find if the line intersects the bounding box
        [x1,x2] = ax.get_xlim()
        [y1,y2] = ax.get_ylim()
        ind = _np.argmin(_np.abs(parTrans[:,0] - x1))
        if _np.logical_and(parTrans[ind,1]>y1,parTrans[ind,1]<y2):
            # Interpolate to axes position
            yInt = _np.interp(x1,parTrans[:,0],parTrans[:,1])
            projYticks.append([aa,yInt])

    if len(projYticks) > 0:
        projYticks = _np.array(projYticks)
        ax.set_yticks(projYticks[:,1])
        labelStr = []
        for aa in projYticks[:,0]:
            if _np.isclose(aa,0):
                hem = ''
            elif aa > 0:
                hem = 'N'
            else:
                hem = 'S'
            labelStr.append('{:2.0f}'.format(_np.abs(aa)) + r'$^{\circ}$' + hem)
        ax.set_yticklabels(labelStr)

    # Deal with meridians    
    projXticks = []
    for ii,aa in enumerate(meridians):
        tmpY = _np.arange(-90,90,0.001)
        tmpX = _np.ones_like(tmpY) * aa
        parTrans = proj.transform_points(dataProj,tmpX,tmpY)
        ax.plot(parTrans[:,0],parTrans[:,1],color='gray',linewidth=0.25)

        # Find if the line intersects the bounding box
        [x1,x2] = ax.get_xlim()
        [y1,y2] = ax.get_ylim()
        ind = _np.argmin(_np.abs(parTrans[:,1] - y1))
        if _np.logical_and(parTrans[ind,0]>x1,parTrans[ind,0]<x2):
            # Interpolate to axes position
            xInt = _np.interp(y1,parTrans[:,1],parTrans[:,0])
            projXticks.append([aa,xInt])

    if len(projXticks) > 0:
        projXticks = _np.array(projXticks)
        ax.set_xticks(projXticks[:,1])
        labelStr = []
        for aa in projXticks[:,0]:
            if _np.isclose(aa,180.0):
                hem = ''
            elif aa < 180.0:
                hem = 'E'
            else:
                hem = 'W'
            aa = _np.abs(_gangles.wrapto180(aa))
            labelStr.append('{:3.0f}'.format(aa) + r'$^{\circ}$' + hem)
        ax.set_xticklabels(labelStr)

    return ax
