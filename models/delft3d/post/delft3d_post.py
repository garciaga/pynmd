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

from pynmd.plotting.vars_param import *
import pynmd.plotting.plot_settings as ps
import datetime

#import seaborn as sns
#plt.style.use('seaborn-white')


### Funcs
def read_d3d_trim(fname, all = True):
    """
    To read D3D trim*.nc file
    Reads: x,y,alfa(local rotation) of cells
    also returns nc.varibles for reading other vars in case
    
    return xcor,ycor,z_cen,alfa,maskxy,ncv
    """
    
    nc  = netCDF4.Dataset(fname)
    ncv = nc.variables
    
    if all:
        xcor  = ncv['XCOR'][:]
        ycor  = ncv['YCOR'][:]
        dep   = ncv['DP0'][:]
        dep = np.ma.masked_where(dep<-100,dep)
        maskxy = dep.mask
        zeta   = ncv['S1'][:]
        nt     = zeta.shape[0]
        alfa  = np.deg2rad (ncv['ALFAS'][:]) 
        #
        sig_layer_face   = ncv['SIG_INTF'][:]
        n_sig_face = len(sig_layer_face)
        #z_face = np.zeros((n_sig_face,xcor.shape[0],xcor.shape[1]))
        #for iz in range(n_sig_face):
        #    z_face[iz] =  sig_layer_face[iz] * dep     
        #
        #sigl = sig_layer_face[:-1] -sig_layer_face[1:]
        
        sig_layer_center = ncv['SIG_LYR'][:]
    #     sigh  =   sig_layer_face  [:-1]   - sig_layer_face[1:]
    #     
    #     h     = np.zeros((sigh.shape[0],xcor.shape[0],xcor.shape[1]))
    #     z_cen = np.zeros((nt,sigh.shape[0],xcor.shape[0],xcor.shape[1]))
    #     for it in range(nt):
    #         for iz in range(sigh.shape[0]):
    #             h[iz] = sigh[iz] * (dep +   zeta[it])
    #         
    #         for iz in range(sigh.shape[0]):
    #             z_cen[it,iz,:] =  np.cumsum(h,0)[iz] - (dep + zeta[it])      
      
        
        
        n_sig_cen = len(sig_layer_center)
        #
        z_cen = np.zeros((nt,n_sig_cen,xcor.shape[0],xcor.shape[1]))
        for it in range(nt):
            for iz in range(n_sig_cen):
                z_cen[it,iz] =  (1+sig_layer_center[iz])  * (dep  +   zeta[it]) - dep  
     
        z_cen = np.ma.masked_where(z_cen < -10, z_cen)
        #z_cen = 10 - z_cen
        
        return xcor,ycor,z_cen,alfa,maskxy,ncv
    else:
        return ncv

def read_d3d_time(netcdf_vars):
    """
    Read time from netcdf file and return Datetime vector
    """
    utime = netcdftime.utime(netcdf_vars['time'].units)
    dates = utime.num2date(netcdf_vars['time'][:])
    return dates

def rotate_3d_vel(u1,v1,w1=None,ang=0):
    import okean.calc as calc

    """
    horizontal rotation of D3D 3D flow field. Fro Curvilinear grids 
    """
    if w1 is not None:
        u1,v1,w1   = calc.rot3d (u1,v1,w1,ang,0.,inverse=False)    
    else:
        u1,v1,w1   = calc.rot3d (u1,v1,v1*0,ang,0.,inverse=False)    
    return u1,v1,w1

def slice_plot(fig,axp,data,args):
    """
    plot slice 2D plot
    fig     : handel to figure
    axp     : handel to axis 
    data    : a dict contan required data (x,y,u,v,val,dens,var,var_def,limits)
    args    : a dict containing optional info for plotting
        plot_cb : plot colorbar
        vec     : plot vectors
        skipx , skipy 
    """
    #plot options
    width      = 0.015
    linewidth  = 0.1
    
    # read main inputs
    xx    = data['xx'][:]
    yy    = data['yy'][:]
    #
    uu    = data['uu'][:]
    vv    = data['vv'][:]
    val   = data['val'][:]
    #
    var   = data['var']
    lim   = data['lim']
    slice = data['slice']
    #
    
    ### Check for args
    if 'skipx' in args.keys():
        skipx = args['skipx']
    else:
        skipx = 1
    
    if 'skipy' in args.keys():
        skipy = args['skipy']
    else:
        skipy = 1
    
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
        
    if 'river_mile' in args.keys():
        river_mile = args['river_mile']
    else:
        river_mile = False         
            
    
    #    
    if 'dens' in data.keys():
        dens =  data['dens'][:]
        dmin,dmax = defs['dens']['vmin'],defs['dens']['vmax']

    
    #
    if   slice == 'k':
        u_vec_scale = 0.5
        v_vec_scale = 0.05
        #        
        axp.set_xlabel ('X [m]')
        axp.set_ylabel ('Y [m]')
        
    elif slice == 'j':     
        axp.set_xlabel ('X [m]')
        axp.set_ylabel ('Z [m]')
        #
        u_vec_scale = 0.5
        v_vec_scale = 1e-4
    else:
        axp.set_xlabel ('Y [m]')
        axp.set_ylabel ('Z [m]')
        #
        u_vec_scale = 0.2
        v_vec_scale = 1e-3
    
    if 'u_vec_scale' in args.keys():
        u_vec_scale = args['u_vec_scale']

    if 'v_vec_scale' in args.keys():
        v_vec_scale = args['v_vec_scale']
   
    
    #color bar related
    if 'plot_cb' in args.keys():
        plot_cb = args['plot_cb']
    else:
        plot_cb = False
    
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
            
    
    #
    #set plot limits based on var information        
    xmin = lim['xmin']
    xmax = lim['xmax']
    ymin = lim['zmin']
    ymax = lim['zmax']

    if river_mile and slice in ['k','j']:
        axp.set_xlabel ('RM [mile]')
        xx   = xx   / 1609.34
        xmin = xmin / 1609.34
        xmax = xmax / 1609.34

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
    xsub1=pos_ax[0][0]
    xsub2=pos_ax[1][0]
    ysub1=pos_ax[0][1]            
    ysub2=pos_ax[1][1]
    
    #plot    
    if pcolor_plot:
        pc1  = axp.pcolor(xx,yy,val,cmap=cmap)
        pc1.set_clim(cmin , cmax)
    else:
        levels = np.linspace(cmin,cmax,21)
        pc1    = axp.contourf(xx,yy,val,levels=levels,cmap=cmap,extend="both")
        pc2    = axp.contour (xx,yy,val,levels=levels[::1],colors='k', linewidths=0.3)
    
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
        axp.quiver(xx[  ::skipy,::skipx],         yy[  ::skipy,::skipx],
                   uu[  ::skipy,::skipx],wscale * vv[  ::skipy,::skipx],
                   pivot='middle',units='inches',scale=scale,width=width)  #,color=acol)
        # vel reference
        xvec1 = xmin + 0.75  * (xmax - xmin)
        yvec1 = ymin + 0.1   * (ymax - ymin)

        axp.quiver([0,xvec1,],[0,yvec1],[0,u_vec_scale,],[0,0.0,],pivot='tail',
                   units='inches',scale=scale,width=width)
        axp.text  ( 1.0 * xvec1,   yvec1- 0.05 * (ymax - ymin) ,str(u_vec_scale)+'[m/s]',size=9)
    

        axp.quiver([0,xvec1,],[0,yvec1],[0,0.0,],[0,v_vec_scale* wscale,],
                   pivot='tail',units='inches',scale=scale,
                   width=width)
        axp.text  ( xvec1+ 0.001*(xmax-xmin),   yvec1 + 0.05 * (ymax - ymin) ,
                    str(v_vec_scale)+'[m/s]',size=9)
    
    for tick in axp.xaxis.get_major_ticks():
        tick.label.set_rotation(20) 
    
    if False:
        #set plot title
        text22=ps.ordinal[nn] +' '+ title
        xx11=xmin + (xmax-xmin) * 0.00
        yy11=ymax + (ymax-ymin) * 0.05
        axp.set_title(text22)#,horizontalalignment='left')
    
    axp.set_xlim(xmin,xmax)
    axp.set_ylim(ymin,ymax)
    
    #if not axp.is_last_row():
    #   plt.setp( axp, 'xticklabels', [] )
    #   axp.set_xlabel('')
    #if not axp.is_first_col():
       #plt.setp( axp, 'yticklabels', [] )
       #axp.set_ylabel('') 
    
    return 

def time_series_plot(ax,data):

    """
    plot timeseries as point
    ax     : handel to axis 
    data    : a dict contan required data (x,y,u,v,val,dens,var,var_def,limits)
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
    
    if 'xlab' in data.keys():
        title = data['title']
    else:
        title = ' '
    #
    #if not isinstance(xx[0], datetime.datetime.date):
    #    xmin = lim['xmin']
    #    xmax = lim['xmax']
    #    ax.set_xlim(xmin,xmax)
    
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
    
    ax.set_ylabel(var['label'])     
    
    ax.grid('off')
    #if not ax.is_last_row():
    #   plt.setp( ax, 'xticklabels', [] )
    #   ax.set_xlabel('')
    #if not plt.gca().is_first_col():
    #   plt.setp( plt.gca(), 'yticklabels', [] )
    #   plt.ylabel('') 
    #print '   > plot line'
    return 


    
##### end of Funcs
