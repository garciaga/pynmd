# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import pynmd.plotting.plot_settings as ps
from   pynmd.tools.compute_statistics import statatistics
import datetime
import matplotlib.tri as Tri
import matplotlib.dates as mdates


import cartopy.crs as ccrs
from cartopy.io import shapereader
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

from math import floor
from matplotlib import patheffects

"""
Simple routines for plooting time_series, pcolors for structured and trianglur data (from adcirc and delf3d moved here) 
"""
__author__ = "Saeed Moghimi"
__copyright__ = "Copyright 2017, UCAR/NOAA"
__license__ = "GPL"
__version__ = "1.0"
__email__ = "moghimis@gmail.com"



def slice_plot(fig,axp,data,args):
    """
    plot slice 2D plot
    fig     : handel to figure
    axp     : handel to axis 
    data      : a dict contan required data (x,y,u,v,val,dens,var,var_def,limits)
        xx    : x axis
        yy    : y axis
        uu    : u velocity
        vv    : v velocity
        val   : data to plot
        var   : var_def dict
        lim   : limits dict
        slice : 'i' , 'j' and 'k' slice indicator 
        dens  : density field (mainly for vertical cross sections)
        
    args    : a dict containing optional info for plotting
        plot_cb : plot colorbar
        vec     : plot vectors
        skipx , skipy : skip points in x and y direction
        
        vec     : True or False  for plotting vectors
        scale   : vector scale 0.01 ~ 4
        u_vec_scale: size of refrence u vector to plot 
        v_vec_scale: size of refrence v vector to plot 
        
        pcolor_plot: True and False for countor plot
        river_mile: To add river miles    
        
        plot_cb: True/False: to add colorbar 
        cb_dx  : x location (0~1)
        cb_dy  : y location (0~1) 
        cb_size: size (0.1 ~ 1)
        
         
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


def TriMapPlot(fig,axp,data,args):
    """
    plot 2D maps for trangular data
    
    fig     : handel to figure
    axp     : handel to axis 
    data      : a dict contan required data (x,y,u,v,val,dens,var,var_def,limits)
        xx    : x axis
        yy    : y axis
        uu    : u velocity
        vv    : v velocity
        val   : data to plot
        var   : var_def dict
        lim   : limits dict
        slice : 'i' , 'j' and 'k' slice indicator 
        dens  : density field (mainly for vertical cross sections)
        
    args    : a dict containing optional info for plotting
        plot_cb : plot colorbar
        vec     : plot vectors
        skipx , skipy : skip points in x and y direction
        
        vec     : True or False  for plotting vectors
        scale   : vector scale 0.01 ~ 4
        u_vec_scale: size of refrence u vector to plot 
        v_vec_scale: size of refrence v vector to plot 
        
        pcolor_plot: True and False for countor plot
        river_mile: To add river miles    
        
        plot_cb: True/False: to add colorbar 
        cb_dx  : x location (0~1)
        cb_dy  : y location (0~1) 
        cb_size: size (0.1 ~ 1)
        
         

    """
    #plot options
    width      = 0.015
    linewidth  = 0.1
    
    # read main inputs
    xx    = data['tri'].x
    yy    = data['tri'].y
    tri   = data['tri']
    #
    if 'uu' in data.keys():
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

def TimeSeriesPlot(ax,data,args):

    """
    plot timeseries as point
    ax     : handel to axis 
    data    : a dict contan required data (x,y,u,v,val,dens,var,var_def,limits)
        data['xx']  : x axis could be datetime
        data['val'] : y axis
              var   : var_def dict
              lim   : limits dict    
    
    
    args    : plot setting stuff
       'xlab'     : x axis label. (if not assumed to be DATETIME) 
       'title'    : panel title
       'obs'      : is it Observation or model (True/False)
       'label'    : Legend label
       'color'    : Line color
       'panel_num': panle number for fetchin ordinal a,b,c,...
       'leg_loc'  : legend location number (best if not given)

    #TODO: check if the xaxis is not date omit the title! 
    """
    
    
    #plot options
    width      = 0.015
    linewidth  = 0.1
    
    var   = data['var']
    lim   = data['lim']
    #
    if 'xlab' in args.keys():
        xlab  = args['xlab']
    else:
        ax.set_xlabel('DateTime')
    
    if 'title' in args.keys():
        title = args['title']
    else:
        title = ' '
    
    if 'obs' in args.keys() and args['obs']:
        marker0    = 'o'
        linestyle0 = 'None'
        lw = 0.5
        #ms = 2  #orig
        #ms = 1  #reviewer 2
        label ='OBS'
    else:
        marker0    = None
        linestyle0 = '-'
        #ms = 2  #orig
        lw = 1.5    

    if 'ms' in args.keys():
        ms = args['ms']
    else:
        ms = 2        
    
    if 'alpha' in args.keys():
        alpha = args['alpha']
    else:
        alpha = 0.75    
    
    
    if 'label' in args.keys():
        label = args['label']
        
    if 'color' in args.keys():
        color = args['color']
    else:
        color = 'k'
        
    if 'panel_num' in args.keys():
        panel_num = args['panel_num']
    else:
        panel_num = None 
    
    xmin,xmax  = lim['xmin'],lim['xmax']        

    if type(data['xx'][0]) is datetime.date:
        if type(xmin) is not datetime.date:
            print (' This is a plot_date call make sure xlims are datetime objs ... ')
            sys.exit('ERROR !')
    
    ax.plot(data['xx'],data['val'],color=color,label=label,
            linestyle=linestyle0,marker=marker0,ms=ms,lw = lw,alpha=alpha)
    

    ax.set_xlim(xmin,xmax)
    ax.set_ylim(var['vmin'],var['vmax'])
    
    if 'plot_leg' in args.keys() and args['plot_leg']:
        if 'leg_loc' in  args.keys():
            loc = args['leg_loc']
        else:
            loc = 'best'    

        leg=ax.legend(loc=loc,ncol=3)
        try:
            frame=leg.get_frame()
            frame.set_edgecolor('None')
            frame.set_facecolor('None')
        except:
            pass
    
    ax.set_title(title)

    # the linewidth of the rectangular axis frame
    #fr_linewidth=0.4
    #[i.set_linewidth(fr_linewidth) for i in ax.spines.itervalues()]
   
    #ax.xaxis.set_ticks(ticks=range(len(point_list))) 
    #ax.xaxis.set_ticklabels(ticklabels=point_list)  #,fontsize=18)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_rotation(30)

    ax.set_ylim(var['vmin'],var['vmax'])
    
    ax.set_ylabel(var['label'])     
 
    xfmt = mdates.DateFormatter('%Y-%m-%d %H:%M')
    ax.xaxis.set_major_formatter(xfmt)
    #to add a) , b) and ..  to top of the panels
    if panel_num is not None:
        ytext = var['vmax'] + 0.05 * (var['vmax']-var['vmin'])         
        if type(xmax) is datetime.date or type(xmax) is datetime.datetime:
            ddt =  (xmax - xmin).total_seconds()
            xtext = xmax - datetime.timedelta(0.05 * ddt/86400.0)
        else:
            xtext = xmax - 0.05 * (xmax-xmin)
        
        ax.text (xtext, ytext, '('+ps.ordinal[panel_num])

    print (var['label'])
    ax.set_ylabel(var['label'])
    ax.set_xlabel('Time')
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,3))
    
    ax.grid('off')
    if not ax.is_last_row():
    #   plt.setp( ax, 'xticklabels', [] )
       ax.set_xlabel('')
    if not plt.gca().is_first_col():
    #   plt.setp( plt.gca(), 'yticklabels', [] )
       plt.ylabel('') 
    #print '   > plot line'
    return 



#### below this line not well tested yet

def read_salt_temp(filename,ncvar,wher,dates,data):
    """
    Find model point based on min rms error 
    
    """

    nc    = netCDF4.Dataset(filename)
    ncv   = nc.variables
    tvar  ='ocean_time'
    utd   = netcdftime.utime(ncv[tvar].units)
    dates_r = utd.num2date(ncv[tvar][:] )
    ind_sim   =  np.where (   (dates_r>=date_first) & (dates_r<=date_last )   )
    option = 'find_min_rmse'
    #option = 'find_min_bias'    
    val  = ncv[ncvar][:,wher,:,:].squeeze()
    nmodj,nmodi = val[0,:].shape
    if option in ['find_min_rmse' ,'find_min_bias']:
        val_rms = []
        val_bis = []
        val_all = []
        for modj in range(nmodj):
            for modi in range(nmodi):
                pre_val = val[:,modj,modi]
                pre_val_interp =  np.interp( utd.date2num(dates),
                                             utd.date2num(dates_r[ind_sim].squeeze()),
                                             pre_val[ind_sim].squeeze())
                bias,rmse,r2=statatistics(data , pre_val_interp)
                val_all.append(pre_val_interp)
                val_rms.append(rmse)
                val_bis.append(bias)
        if option in ['find_min_rmse']:
            indx    = np.argmin(np.array(val_rms))
        else:
            indx    = np.argmin(np.array(val_bis))

        print ('index>',indx) #,'   rmse> ',val_rms[indx],val_rms
        val_interp = val_all[indx]
    nc.close()

    return dates_r[ind_sim].squeeze()[:-1],val_interp  

def plot_scatter(ax,data,model,var=dict(),color='k',marker = None,nn=None, title=None):
    """
    plot scatter
    """
    
    min1 = var['vmin']
    max1 = var['vmax']
    line=[min1,max1]                
    ax.plot   (line,line,'k' ,linewidth=0.2)


    if marker is None:
        marker='.'
    else:
        marker = marker
    scat = ax.scatter(model,data,edgecolor='none', alpha=0.7,c=color,marker=marker)
    ax.set_ylabel('Data')
    ax.set_xlabel('Model')
    ax.set_xlim(min1,max1)
    ax.set_ylim(min1,max1)
    ax.set_aspect(1)    

    fit    = np.polyfit(model,data,1)
    fit_fn = np.poly1d(fit) # fit_fn is now a function which takes in x and returns an estimate for y
    ax.plot(model, fit_fn(model),c=color,lw=1 ,alpha=1.0)
    
    ddt =  max1-min1
    
    if nn is not None:
        txt = '(' + ps.ordinal[nn]
        if title is not None:
            txt = txt +' '+ title
      
        ax.text (max1 - 0.95 * ddt ,var['vmax']+ 0.05 * (var['vmax']-var['vmin']) , txt, fontsize=9)

    stat = statatistics(data,model)
    txt = 'RMSE='+"%3.3f"    % stat['rmse'] +\
          '\nbias='+"%3.3f"  % stat['bias'] +\
          '\nR2='+"%3.3f"    % stat['r2']+\
          '\nskill='+"%3.3f" % stat['skill'] +\
          '\nRB='+"%3.3f"    % stat['rbias'] +\
          '\nIA='+"%3.3f"    % stat['ia'] +\
          '\nN='+"%3i"       % len(model)

    ax.text (min1 + 0.05 * ddt ,var['vmax'] - 0.35 * (var['vmax']-var['vmin']) , txt , fontsize=7)
    print (var['label'])    

    # the linewidth of the rectangular axis frame
    #fr_linewidth=0.5
    #[i.set_linewidth(fr_linewidth) for i in ax.spines.itervalues()]
    
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_rotation(0)
        tick.label.set_fontsize(10)

    #if not ax.is_last_row():
    #   ax.set_xlabel('')
    #   plt.setp( ax, 'xticklabels', [] )
    #if not ax.is_last_col():
    #   ax.set_ylabel('')
    #   plt.setp( ax, 'yticklabels', [] )

    return scat


def ncks(param='zeta',xvar='eta_rho',yvar='xi_rho',ix=0,jy=0,filein='tmp.nc',fileout='tmp2.nc'):
      comm1='ncks -O  -v '+ param+' -d '+xvar +' '+str(jnum)+' -d '+ yvar +' '+str(inum)+' '+filein+' '+fileout
      os.system(comm1)


from matplotlib.dates import date2num
def stick_plot(time, u, v, **kw):
    width = kw.pop('width', 0.002)
    headwidth = kw.pop('headwidth', 0)
    headlength = kw.pop('headlength', 0)
    headaxislength = kw.pop('headaxislength', 0)
    angles = kw.pop('angles', 'uv')
    ax = kw.pop('ax', None)
    
    if angles != 'uv':
        raise AssertionError("Stickplot angles must be 'uv' so that"
                             "if *U*==*V* the angle of the arrow on"
                             "the plot is 45 degrees CCW from the *x*-axis.")

    time, u, v = map(np.asanyarray, (time, u, v))
    if not ax:
        fig, ax = plt.subplots()
    
    q = ax.quiver(date2num(time), [[0]*len(time)], u, v,
                  angles='uv', width=width, headwidth=headwidth,
                  headlength=headlength, headaxislength=headaxislength,
                  **kw)

    ax.axes.get_yaxis().set_visible(False)
    ax.xaxis_date()
    return q

def make_map(projection=ccrs.PlateCarree(), res='m', xylabels = False):
    
    """
    Generate fig and ax using cartopy
    input: projection
    output: fig and ax
    """


    subplot_kw = dict(projection=projection)
    fig, ax = plt.subplots(figsize=(9, 13),
                           subplot_kw=subplot_kw)
    if xylabels:
        gl = ax.gridlines(draw_labels=True)
        gl.xlabels_top = gl.ylabels_right = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
    else:
        gl = ax.gridlines(draw_labels=False)
        gl.xlines = False
        gl.ylines = False


    
    #if res is not None:
    #    if res == 'm':
    #        ax.background_img(name='BM', resolution='high')   # from local hdd you need to > import pynmd.plotting
    #    else:
    #        ax.background_img(name='BMH', resolution='high')   # from local hdd you need to > import pynmd.plotting

    return fig, ax



def utm_from_lon(lon):
    """
    utm_from_lon - UTM zone for a longitude

    Not right for some polar regions (Norway, Svalbard, Antartica)

    :param float lon: longitude
    :return: UTM zone number
    :rtype: int
    """
    return floor( ( lon + 180 ) / 6) + 1

def scale_bar(ax, proj, length, location=(0.5, 0.05), linewidth=3,
              units='km', m_per_unit=1000):
    """

    http://stackoverflow.com/a/35705477/1072212
    ax is the axes to draw the scalebar on.
    proj is the projection the axes are in
    location is center of the scalebar in axis coordinates ie. 0.5 is the middle of the plot
    length is the length of the scalebar in km.
    linewidth is the thickness of the scalebar.
    units is the name of the unit
    m_per_unit is the number of meters in a unit
    """
    # find lat/lon center to find best UTM zone
    x0, x1, y0, y1 = ax.get_extent(proj.as_geodetic())
    # Projection in metres
    utm = ccrs.UTM(utm_from_lon((x0+x1)/2))
    # Get the extent of the plotted area in coordinates in metres
    x0, x1, y0, y1 = ax.get_extent(utm)
    # Turn the specified scalebar location into coordinates in metres
    sbcx, sbcy = x0 + (x1 - x0) * location[0], y0 + (y1 - y0) * location[1]
    # Generate the x coordinate for the ends of the scalebar
    bar_xs = [sbcx - length * m_per_unit/2, sbcx + length * m_per_unit/2]
    # buffer for scalebar
    buffer = [patheffects.withStroke(linewidth=5, foreground="w")]
    # Plot the scalebar with buffer
    ax.plot(bar_xs, [sbcy, sbcy], transform=utm, color='k',
        linewidth=linewidth, path_effects=buffer)
    # buffer for text
    buffer = [patheffects.withStroke(linewidth=3, foreground="w")]
    # Plot the scalebar label
    t0 = ax.text(sbcx, sbcy, str(length) + ' ' + units, transform=utm,
        horizontalalignment='center', verticalalignment='bottom',
        path_effects=buffer, zorder=2)
    left = x0+(x1-x0)*0.05
    # Plot the N arrow
    t1 = ax.text(left, sbcy, u'\u25B2\nN', transform=utm,
        horizontalalignment='center', verticalalignment='bottom',
        path_effects=buffer, zorder=2)
    # Plot the scalebar without buffer, in case covered by text buffer
    ax.plot(bar_xs, [sbcy, sbcy], transform=utm, color='k',
        linewidth=linewidth, zorder=3)

def readTrack ( atcfFile ):
    """
    Reads ATCF-formatted file
    Args:
        'atcfFile': (str) - full path to the ATCF file
    Returns:
        dict: 'lat', 'lon', 'vmax', 'mslp','dates'
    """
    lines = open(atcfFile).readlines()
        
    myOcn  = []
    myCy   = []
    myDate = []
    myLat  = []
    myLon  = []
    myVmax = []
    myMSLP = []
    for line in lines:
        r = line.rstrip().split(',')
        myOcn.append  (r[0])
        myCy.append   (int(r[1]))
        myDate.append (datetime.strptime(r[2].strip(),'%Y%m%d%H'))
        latSign = -1.0
        if 'N' in r[6]:
            latSign = 1.0     
        myLat.append  (latSign*0.1*float(r[6][:-1]))
        lonSign = -1.0
        if 'E' in r[7]:
            lonSign = 1.0
        myLon.append  (lonSign*0.1*float(r[7][:-1]))
        myVmax.append (float(r[8]))
        myMSLP.append (float(r[9]))
    
    return { 
            'basin' : myOcn,    'cy' : myCy, 'dates' : myDate, 
            'lat'   : myLat,   'lon' : myLon,
            'vmax'  : myVmax, 'mslp' : myMSLP }


def read_track(ax,path,date):
    ike_track_file = '/scratch4/COASTAL/coastal/save/Saeed.Moghimi/models/NEMS/NEMS_inps/data/tracks/ike_bal092008.dat'
    track = readTrack(ike_track_file)
    keys = ['dates', 'lon', 'vmax', 'lat']
    for key in keys:
        tmp   = pd.DataFrame(track[key],columns=[key])

        #dfh   = df
        if 'trc' not in globals():
            trc = tmp
        else:
            trc  = pd.concat([trc,tmp],axis=1,join_axes=[trc.index])    
    
    
    
    trc = trc.drop_duplicates(subset='dates',keep='first')
    trc = trc.set_index (trc.dates)
    trc = trc.resample('H').interpolate()
    trc.drop('dates',axis=1,inplace=True)
    
    dates = datetime64todatetime(trc.index)
    
    return dates,trc.lon.values, trc.lat.values
    

def plot_track(ax,track,date=None):
    
    if date is not None:
        dates = np.array(track['dates'])
        ind = np.array(np.where((dates==date))).squeeze().item()
        ax.scatter(lon[ind],lat[ind],s=50,c='r',alpha=50)
    ax.plot(track['lon'],track['lat'],lw=2,c='r')


#def imscatter(x, y, image=None, ax=None, zoom=0.05):
#    """
#    Plot image as marker at x and y like scatter plots
#    """
#    from matplotlib.offsetbox import OffsetImage, AnnotationBbox
   
#    if image is None:
#        image = '/scratch4/COASTAL/coastal/save/Saeed.Moghimi/models/NEMS/NEMS_inps/01_data/cliparts/hurricane-1085673_960_720_red.png'
    
    
#    if ax is None:
#        ax = plt.gca()
#    try:
#        image = plt.imread(image)
#    except TypeError:
#        # Likely already an array...
#        pass
#    im = OffsetImage(image, zoom=zoom)
#    x, y = np.atleast_1d(x, y)
#    artists = []
#    for x0, y0 in zip(x, y):
#        ab = AnnotationBbox(im, (x0, y0), xycoords='data', frameon=False)
#        artists.append(ax.add_artist(ab))
#    ax.update_datalim(np.column_stack([x, y]))
#    ax.autoscale()
#    return artists

if __name__ == '__main__':

    ax = plt.axes(projection=ccrs.Mercator())
    plt.title('Cyprus')
    ax.set_extent([31, 35.5, 34, 36], ccrs.Geodetic())
    ax.stock_img()
    ax.coastlines(resolution='10m')

    scale_bar(ax, ccrs.Mercator(), 100)  # 100 km scale bar
    # or to use m instead of km
    # scale_bar(ax, ccrs.Mercator(), 100000, m_per_unit=1, units='m')
    # or to use miles instead of km
    # scale_bar(ax, ccrs.Mercator(), 60, m_per_unit=1609.34, units='miles')
