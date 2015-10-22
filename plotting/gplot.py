# Import modules
from __future__ import division,print_function
import pylab as pl
import numpy as np

# Jet colorbar without green
def jetWoGn(reverse=False):
    """
    jetWoGn(reverse=False)
       - returning a colormap similar to cm.jet, but without green.
         if reverse=True, the map starts with red instead of blue.
         
    Notes 
    -----
    Courtesy of Dr. Saeed Moghimi
    
    """
    m=18 # magic number, which works fine
    m0=pl.floor(m*0.0)
    m1=pl.floor(m*0.2)
    m2=pl.floor(m*0.2)
    m3=pl.floor(m/2)-m2-m1

    b_ = np.hstack( (0.4*np.arange(m1)/(m1-1.)+0.6, np.ones((m2+m3,)) ) )
    g_ = np.hstack( (np.zeros((m1,)),np.arange(m2)/(m2-1.),np.ones((m3,))) )
    r_ = np.hstack( (np.zeros((m1,)),np.zeros((m2,)),np.arange(m3)/(m3-1.)))

    r = np.hstack((r_,pl.flipud(b_)))
    g = np.hstack((g_,pl.flipud(g_)))
    b = np.hstack((b_,pl.flipud(r_)))

    if reverse:
        r = pl.flipud(r)
        g = pl.flipud(g)
        b = pl.flipud(b)

    ra = pl.linspace(0.0,1.0,m)

    cdict = {'red': zip(ra,r,r),
            'green': zip(ra,g,g),
            'blue': zip(ra,b,b)}

    return pl.matplotlib.colors.LinearSegmentedColormap('new_RdBl',cdict,256)


#===============================================================================
# Get colormap colors
#===============================================================================
def get_colormap_colors(N,cmapname='jet'):
    """
    
    PARAMETERS:
    -----------
    N          : Number of colors to return
    cmapname   : Name of colormap (defaults to jet)
    
    RETURNS:
    --------
    colors     : Nx4 list of colors [(R,G,B,1.0)]
    
    """
    
    # Get colormap
    cl_map = pl.cm.get_cmap(cmapname)
     
    # Extract all the colours
    cl_col = np.array([cl_map(ii) for ii in range(cl_map.N)])
     
    # Linear interpolation
    r_N = np.interp(np.arange(N)/(N-1.0),
                    np.arange(cl_map.N)/(cl_map.N-1.0),cl_col[:,0])
    g_N = np.interp(np.arange(N)/(N-1.0),
                    np.arange(cl_map.N)/(cl_map.N-1.0),cl_col[:,1])
    b_N = np.interp(np.arange(N)/(N-1.0),
                    np.arange(cl_map.N)/(cl_map.N-1.0),cl_col[:,2])          
     
    # Allocate in array
    colors = []
    for aa in range(r_N.shape[0]):
        colors.append((r_N[aa],g_N[aa],b_N[aa],1.0))
         
    return colors

         
         
     
