__author__ = "Saeed Moghimi"
__copyright__ = "Copyright 2015, Portland State University"
__license__ = "GPL"
__version__ = "1.0"
__email__ = "moghimis@gmail.com"



import netCDF4
import netcdftime

import matplotlib.pyplot as plt
import numpy as np
import os,sys
from   pynmd.plotting.vars_param import *
import pynmd.plotting.plot_settings as ps
import datetime
import matplotlib.tri as Tri


#import seaborn as sns
#plt.style.use('seaborn-white')


### Funcs
def ReadDates(DirName):
    fname =  DirName + 'fort.63.nc'  
    nc  = netCDF4.Dataset(fname)
    ncv = nc.variables
    t_var   = nc.variables['time'] 
    dates = netCDF4.num2date(t_var[:],t_var.units)
    return dates


def ReadElev(DirName, tind):
    """
    fname: fort.63.nc file
    tind: time index
    """
    fname =  DirName + 'fort.63.nc'  
    nc  = netCDF4.Dataset(fname)
    ncv  = nc.variables
    t_var= nc.variables['time'] 
    date = netCDF4.num2date(t_var[tind],t_var.units)
    elev = nc.variables['zeta'][tind].squeeze()
    nc.close()
    return {'ncv':ncv,'date':date,'elev':elev}    

def ReadUV(DirName, tind):
    """
    fname: fort.64.nc file
    tind: time index
    """
    fname =  DirName + 'fort.64.nc'  
    nc  = netCDF4.Dataset(fname)
    ncv = nc.variables
    t_var= nc.variables['time'] 
    date = netCDF4.num2date(t_var[tind],t_var.units)
    u    = nc.variables['u-vel'][tind].squeeze()
    v    = nc.variables['v-vel'][tind].squeeze()
    uv   = np.sqrt(u*u + v*v)
    nc.close()

    return {'ncv':ncv,'date':date,'u':u,'v':v,'uv':uv}   

def ReadVar(fname='',varname='',time_name=None , tind = None):
    """
    fname: ncfile name
    varname:
    """
    out = {}
    nc  = netCDF4.Dataset(fname)
    ncv = nc.variables
    try:
        if time_name is not None:
            tname = time_name
        else:
            tname =   'time' 
        t_var= nc.variables[tname] 
        if tind is None:
           date = netCDF4.num2date(t_var[:],t_var.units)
        else:
           date = netCDF4.num2date(t_var[tind],t_var.units)
        
        out.update({'date':date})
    
    except:
        pass

    if tind is None:
        var = nc.variables[varname][:].squeeze()
    else:
        var = nc.variables[varname][tind].squeeze()

    x   = nc.variables['longitude'][:]
    y   = nc.variables['latitude'][:]
    # read connectivity array
    el  = nc.variables['tri'][:] - 1
    # create a triangulation object, specifying the triangle connectivity array
    tri = Tri.Triangulation(x,y, triangles=el)
    nc.close()
    out.update({'ncv':ncv,varname:var, 'tri':tri}) 
    return out   

def ReadTri(DirName):

    """
    fname: one of fort.*.nc file
    tind: time index
    """ 
    fname =  DirName + 'fort.63.nc'  
    nc  = netCDF4.Dataset(fname)
    ncv = nc.variables
    x   = nc.variables['x'][:]
    y   = nc.variables['y'][:]
    # read connectivity array
    el  = nc.variables['element'][:] - 1
    # create a triangulation object, specifying the triangle connectivity array
    tri = Tri.Triangulation(x,y, triangles=el)
    nc.close()
    return x,y,tri

def time_series_plot(ax,data):

    """
    plot timeseries as point
    ax      : handel to axis 
    data    : a dict contain required data (x,y,u,v,val,dens,var,var_def,limits)
    """
    #plot options
    width      = 0.015
    linewidth  = 0.1
    pcolor_plot = True
    
    # read main inputs
    xx     = data['xx'][:]
    val    = data['val'][:]
    #
    var   = data['var']
    lim   = data['lim']
    #
    if 'xlab' in data.keys():
        xlab  = data['xlab']
    else:
        ax.set_xlabel('X[m]')
    
    if 'title' in data.keys():
        title = data['title']
        print title
    else:
        title = ' '
    #
    ax.plot(xx,val ,'k-o',label=var['label'])
    #plt.plot(x2[0],data,'k-o',label='data' )
    leg=ax.legend(loc='best')
    try:
        frame=leg.get_frame()
        frame.set_edgecolor('None')
        frame.set_facecolor('None')
    except:
        pass
    ax.set_title(title)

    # the linewidth of the rectangular axis frame
    fr_linewidth=0.4
    [i.set_linewidth(fr_linewidth) for i in ax.spines.itervalues()]
   
    #ax.xaxis.set_ticks(ticks=range(len(point_list))) 
    #ax.xaxis.set_ticklabels(ticklabels=point_list)  #,fontsize=18)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_rotation(30)

    ax.set_ylim(var['vmin'],var['vmax'])
    
    #if not ax.is_last_row():
    #   plt.setp( ax, 'xticklabels', [] )
    #   ax.set_xlabel('')
    #if not plt.gca().is_first_col():
    #   plt.setp( plt.gca(), 'yticklabels', [] )
    #   plt.ylabel('') 
    #print '   > plot line'
    return 


def MapPlot(fig,axp,data,args):
    """
    plot slice 2D plot
    fig     : handel to figure
    axp     : handel to axis 
    data    : a dict contain required data (x,y,u,v,val,dens,var,var_def,limits)
    args    : a dict contain optional info for plotting
        plot_cb : plot colorbar
        vec     : plot vectors
        skipx , skipy 
    """
    #plot options
    width      = 0.015
    linewidth  = 0.1
    
    # read main inputs
    xx    = data['tri'].x
    yy    = data['tri'].y
    tri   = data['tri']
    #
    uu    = data['uu'][:]
    vv    = data['vv'][:]
    val   = data['val'][:]
    #
    var   = data['var']
    lim   = data['lim']
    #
    
    ### Check for args
    if 'subsamp' in args.keys():
        subsamp = args['subsamp']
    else:
        subsamp = 1
    
    if 'scale' in args.keys():
        scale   = args['scale']
    else:
        scale = 1
    
    if 'vec' in args.keys():
        vec = args['vec']
    else:
        vec = False
    
    if 'pcolor_plot' in args.keys():
        pcolor_plot = args['pcolor_plot']
    else:
        pcolor_plot = True        
        
    
    if 'xlab' in args.keys():
        xlab  = args['xlab']
    else:
        axp.set_xlabel('Lon [deg]')
    
    if 'ylab' in args.keys():
        ylab  = args['ylab']
    else:
        axp.set_ylabel('Lat [deg]')
    
    #tri.y is Lat.
    axp.set_aspect=1.0/np.cos(tri.y.mean() * np.pi / 180.0)
    
    if 'title' in args.keys():
        title = args['title']
    else:
        title = ' '
    
    if 'u_vec_scale' in args.keys():
        u_vec_scale = args['u_vec_scale']
    else:
        u_vec_scale = 0.5
    
    if 'plot_cb' in args.keys():
        plot_cb = args['plot_cb']
    else:
        plot_cb = False
    
    if 'v_vec_scale' in args.keys():
        v_vec_scale = args['v_vec_scale']
    else:
        v_vec_scale = 0.5
    
    if 'cb_dx'  in args.keys():
        cb_dx = args['cb_dx']
    else:
        cb_dx = 1.02
    
    if 'cb_dy'  in args.keys():
        cb_dy = args['cb_dy']
    else:
        cb_dy = 0.25
    
    if 'cb_size'  in args.keys():
        cb_size = args['cb_size']
    else:
        cb_size = 0.5
     
    
    #set plot limits based on var information        
    xmin = lim['xmin']
    xmax = lim['xmax']
    ymin = lim['ymin']
    ymax = lim['ymax']
    
    cmin   = var['vmin']
    cmax   = var['vmax']
    cbtext = var['label']
    cmap   = var['cmap']
    
    # set scale for vertical vector plots
    pos_ax   = np.array (axp.get_position ())
    aheight  = pos_ax[1][1] -pos_ax[0][1]
    awidth   = pos_ax[1][0] -pos_ax[0][0]
    #only for vector plot
    fwidth  = fig.get_size_inches()[0]
    fheight = fig.get_size_inches()[1]
    wscale = aheight*fheight/(awidth*fwidth) * (xmax-xmin)/(ymax-ymin)
    #
    xsub1 = pos_ax[0][0]
    xsub2 = pos_ax[1][0]
    ysub1 = pos_ax[0][1]            
    ysub2 = pos_ax[1][1]
    
        
    #plot    
    if pcolor_plot:
        pc1   = axp.tripcolor(tri, val,cmap=cmap)
        pc1.set_clim(cmin , cmax)
    else:
        levels = np.linspace(cmin,cmax,21)
        pc1    = axp.tricontourf(tri, val,shading='faceted',levels=levels,cmap=cmap,extend="both")
        pc2    = axp.tricontour (tri, val,levels=levels[::1],colors='k', linewidths=0.3)
    
    #plot colorbar
    if plot_cb:
        cbax    = fig.add_axes([xsub1+ cb_dx * (xsub2-xsub1), ysub1 + cb_dy * (ysub2-ysub1), 0.01, cb_size]) 
        cbticks = np.linspace(cmin, cmax, 3, endpoint=True)
        cb      = plt.colorbar(pc1,cax=cbax,ticks=cbticks,format='%1.4g',orientation='vertical') #
        [i.set_linewidth(0.01) for i in cb.ax.spines.itervalues()]
        cb.set_label(cbtext)
    
    #plot bottom line
    axp.plot(xx[-1],yy[-1],'k',lw=2)
    
    #plot surface line
    axp.plot(xx[0],yy[0],'k',lw=0.5)
    
    #plot density contours
    if 'cont' in data.keys():
        cont = data['cont'][:]
        dmin = data['cvar']['vmin']
        dmax = data['cvar']['vmax']
        if True:
            levels = [1,2,10,20,30]
        else:
            levels = np.linspace(dmin,dmax,7)
       
        Cabc   = axp.contour (xx,yy,cont,levels,colors='k',linewidths=20*linewidth)
        #plt.clabel(Cabc,fontsize=6,fmt='%.4g')
    
    #plot vectors
    if vec: 
        ind = np.argwhere((xx>xmin) & (xx<xmax) & (yy>ymin) & (yy<ymax))
        np.random.shuffle(ind)
        Nvec = int(len(ind) / subsamp)
        idv = ind[:Nvec]
    
        Q = axp.quiver(xx[idv],yy[idv],uu[idv], wscale * vv[idv],
                   pivot='middle',units='inches',scale=scale,width=width)  #,color=acol)
        # vel reference
        xvec1 = xmin + 0.1  * (xmax - xmin)
        yvec1 = ymin + 0.8  * (ymax - ymin)
    
        if False:
            axp.quiver([0,xvec1,],[0,yvec1],[0,u_vec_scale,],[0,0.0,],pivot='tail',
                       units='inches',scale=scale,width=width)
            axp.text  ( 1.0 * xvec1,   yvec1- 0.05 * (ymax - ymin) ,str(u_vec_scale)+'[m/s]',size=9)
        
            axp.quiver([0,xvec1,],[0,yvec1],[0,0.0,],[0,v_vec_scale* wscale,],
                       pivot='tail',units='inches',scale=scale,
                       width=width)
            axp.text  ( xvec1+ 0.001*(xmax-xmin),   yvec1 + 0.05 * (ymax - ymin) ,
                        str(v_vec_scale)+'[m/s]',size=9)
        else:
            maxstr='%3.1f m/s' % u_vec_scale
            qk = axp.quiverkey(Q,0.2,0.9,u_vec_scale,maxstr,labelpos='W')
    
    
    for tick in axp.xaxis.get_major_ticks():
        tick.label.set_rotation(20) 
    
    if True:
        #text22=ps.ordinal[nn] +' '+ title
        text22 = title
        print text22
        axp.set_title(text22,loc='left')
    
    axp.set_xlim(xmin,xmax)
    axp.set_ylim(ymin,ymax)
    
    if not axp.is_last_row():
       plt.setp( axp, 'xticklabels', [] )
       axp.set_xlabel('')
    if not axp.is_first_col():
       plt.setp( axp, 'yticklabels', [] )
       axp.set_ylabel('') 
    
    return 


  
##### end of Funcs
