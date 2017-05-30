# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import pynmd.plotting.plot_settings as ps
import datetime
import matplotlib.tri as Tri


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
        ms = 2  #orig
        #ms = 1  #reviewer 2
        label ='OBS'
    else:
        marker0    = None
        linestyle0 = '-'
        ms = 2  #orig
        lw = 1.5    

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
    
    ax.plot(data['xx'],data['val'] ,color=color,label=label,
            linestyle=linestyle0,marker=marker0,ms=ms,lw = lw)
    
    xmin,xmax  = lim['xmin'],lim['xmax']
    date_first = xmin
    date_last  = xmax
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
    fr_linewidth=0.4
    [i.set_linewidth(fr_linewidth) for i in ax.spines.itervalues()]
   
    #ax.xaxis.set_ticks(ticks=range(len(point_list))) 
    #ax.xaxis.set_ticklabels(ticklabels=point_list)  #,fontsize=18)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_rotation(30)

    ax.set_ylim(var['vmin'],var['vmax'])
    
    ax.set_ylabel(var['label'])     
 
    #to add a) , b) and ..  to top of the panels
    if panel_num is not None:
        ytext = var['vmax']+ 0.05 * (var['vmax']-var['vmin'])         
        if type(date_last) is datetime.date:
            ddt =  (date_last - date_first).total_seconds()
            xtext = date_last - datetime.timedelta(0.05 * ddt/86400.0)
        else:
            xtext = xmax - 0.05 (xmax-xmin)
        
        ax.text (xtext, ytext, ps.ordinal[panel_num])

    print var['label']
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




