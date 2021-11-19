# -*- coding: utf-8 -*-

from __future__ import division,print_function

__author__ = "Saeed Moghimi"
__copyright__ = "Copyright 2015, Portland State University"
__license__ = "GPL"
__version__ = "1.0"
__email__ = "moghimis@gmail.com"


"""
A series of tools to process and analyze adcirc output

Authors:
-------
Fadia Ticona Rollano
    PNNL Marine Sciences Laboratory

Log of edits:
-------------
December 2019 - Created module
    Fadia Ticona Rollano

"""


import time
#sys.path.append(os.getcwd())
import getpass
import numpy as np
import netCDF4

# Custom paths
import pynmd.models.adcirc.pre as adcpre
import pynmd.data.angles as gangles


import matplotlib.pyplot as plt
import os,sys
from   pynmd.plotting.vars_param import *
import pynmd.plotting.plot_settings as ps
import datetime
import matplotlib.tri as Tri

import cartopy.crs as ccrs
from cartopy.io import shapereader
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

from math import floor
from matplotlib import patheffects

#import seaborn as sns
#plt.style.use('seaborn-white')


### Funcs
def ReadDates(DirName):
    fname =  os.path.abspath(DirName + 'fort.63.nc')  
    #print fname
    nc  = netCDF4.Dataset(fname)
    ncv = nc.variables
    t_var   = nc.variables['time'] 
    dates = netCDF4.num2date(t_var[:],t_var.units)
    nc.close()
    return dates


def ReadElev(DirName, tind):
    """
    fname: fort.63.nc file
    tind: time index
    """
    fname =  os.path.abspath(DirName + 'fort.63.nc') 
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
    fname =  os.path.abspath(DirName + 'fort.64.nc' ) 
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
    fn  = os.path.abspath(fname)
    #print fn
    nc  = netCDF4.Dataset(fn)
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

    # x   = nc.variables['longitude'][:]
    # y   = nc.variables['latitude'][:]
    # read connectivity array
    # el  = nc.variables['tri'][:] - 1
    # create a triangulation object, specifying the triangle connectivity array
    # tri = Tri.Triangulation(x,y, triangles=el)
    nc.close()
    out.update({'ncv':ncv,varname:var}) 
    return out   

def maskTri(Tri,zmask):
    """
    Inputs: 
    tri object
    mask array of vertex
    
    Returned: maksed tri object
    """
    
    print ('[info:] Mask Tri ... ')
    mask = np.ones(len(Tri.triangles), dtype=bool)
    count = 0
    for t in Tri.triangles:
        count+=1
        ind = t
        if np.any(zmask[ind-1]):
            mask[count-1] = False    
    Tri.set_mask = mask
    return Tri

def maskTri_v2 (Tri,mask):
    m = np.all(mask[Tri.triangles],axis=1) 
    Tri.set_mask = m
    return Tri

def ReadTri_v1(DirName):

    """
    fname: one of fort.*.nc file
    tind: time index
    """ 
    fname =  os.path.abspath(DirName + 'fort.63.nc'  )
    nc  = netCDF4.Dataset(fname)
    ncv = nc.variables
    x   = nc.variables['x'][:]
    y   = nc.variables['y'][:]
    # read connectivity array
    el  = nc.variables['element'][:] - 1
    # create a triangulation object, specifying the triangle connectivity array
    print ('[info:] Generate Mask ...')
    tri  = Tri.Triangulation(x,y, triangles=el)
    try:
        zeta = nc.variables['zeta'][0].squeeze()
        #zeta = np.ma.masked_where(np.isnan(zeta),zeta)
        tri = tri_mask(tri,zeta.mask)
    except:
        print (' Tri mask did not applied !')
        pass
    nc.close()
    return x,y,tri


def ReadTri(DirName):

    """
    fname: one of fort.*.nc file
    tind: time index
    """ 
    fname =  os.path.abspath(DirName + '/maxele.63.nc'  )
    nc  = netCDF4.Dataset(fname)
    ncv = nc.variables
    x   = nc.variables['x'][:]
    y   = nc.variables['y'][:]
    # read connectivity array
    el  = nc.variables['element'][:] - 1
    # create a triangulation object, specifying the triangle connectivity array
    print ('[info:] Generate Tri ...')
    tri  = Tri.Triangulation(x,y, triangles=el)
    if False:
        try:
            zeta = nc.variables['zeta_max'][:].squeeze()
            zeta = np.ma.masked_where(np.isnan(zeta),zeta)
            tri = tri_mask(tri,zeta.mask)
            print ('[info:] Generate Tri.mask ...')

        except:
            print (' Tri mask did not applied !')
            pass
    nc.close()
    return x,y,tri


def ReadFort80(dir):
    """
    Read fort.80 file for domain decomposition information
    
    return a dictionary including:
    IMAP_EL_LG
    IMAP_NOD_LG
    IMAP_NOD_GL (negative values for not owned elements)
   
    """
    fdata = open( dir + '/fort.80' ,  'r')
    while True:
    #for  line in fdata.readlines():
        line = fdata.readline()
        if 'Number of processors'     in line: nproc = int(line.split()[0]) 
        if 'Total # elements & nodes' in line: 
            nelem = int(line.split()[0])
            nnode = int(line.split()[1])          
        if 'NWLON, NWLAT'             in line:
            print (line)
            break

    #allocate
    IMAP_NOD_LG  = []
    IMAP_NOD_GL   = np.zeros((nnode,3),dtype=np.int)
    IMAP_EL_LG   = []
    #IMAP_STAE_LG = []
    #IMAP_STAV_LG = []
    #IMAP_STAC_LG = []
    #IMAP_STAM_LG = []
    pe_all        = []
    print ('[info:] read nodes local2global')
    for inp in  range( nproc ):
        line1       = fdata.readline()
        #print line1
        pe          = int(line1.split()[0])
        nnodp       = int(line1.split()[1])                
        nod_res_tot = int(line1.split()[2])                
        pe_all.append(pe)
        tmpa = np.array([])
        proc_read = True
        while proc_read:
           line1 = fdata.readline()
           tmp = np.array([int(v) for v in line1.split()])
           tmpa = np.r_[tmpa,tmp]
           if len(tmpa) == nnodp:
               IMAP_NOD_LG.append(tmpa)
               proc_read = False
      
    print ('[info:] read nodes local2global')
    line1       = fdata.readline()
    #print line1
    for il in range(nnode):
        line1       = fdata.readline()
        node_globa = int(line1.split()[0])
        pe         = int(line1.split()[1])                
        node_local = int(line1.split()[2]) 
        IMAP_NOD_GL[il,:] = node_globa , pe ,  node_local

    print ('[info:] read elements local2global')
    for inp in  range( nproc ):
        line1       = fdata.readline()
        #print line1
        pe          = int(line1.split()[0])   # pe number
        nelmp       = int(line1.split()[1])   # element on pe             

        pe_all.append(pe)
        tmpa = np.array([])
        proc_read = True
        while proc_read:
           line1 = fdata.readline()
           tmp = np.array([int(v) for v in line1.split()])
           tmpa = np.r_[tmpa,tmp]
           if len(tmpa) == nelmp:
               IMAP_EL_LG.append(tmpa)
               proc_read = False
           
    fdata.close() 
    return (dict (nproc = nproc, nelem = nelem, nnode = nnode ,
                  IMAP_EL_LG  = np.array(IMAP_EL_LG) ,
                  IMAP_NOD_LG = np.array(IMAP_NOD_LG),
                  IMAP_NOD_GL = np.array(IMAP_NOD_GL) )) 


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


    
    if res is not None:
        if res == 'm':
            ax.background_img(name='BM', resolution='high')   # from local hdd you need to > import pynmd.plotting
        else:
            ax.background_img(name='BMH', resolution='high')   # from local hdd you need to > import pynmd.plotting

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


def maskDryElements(grid):
    dry_masked = np.ma.masked_where(grid['depth']<=0., grid['depth'])
    return np.all(dry_masked.mask[grid['Elements']-1],axis=1)

def maskTri(tri,mask):
    return np.all(mask[tri.triangles],axis=1)

def maskTolExceed(OldGrid,NewGrid,tol=0.1):
    diff = OldGrid['depth'] - NewGrid['depth']
    diffm = np.ma.masked_where(np.abs(diff) <= tol*np.abs(OldGrid['depth']), diff)
    return np.all(diffm.mask[NewGrid['Elements']-1],axis=1)


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



def find_hwm(xgrd,ygrd,maxe,xhwm,yhwm,elev_hwm,convert2msl=None,bias_cor=None ,flag='pos'):
    from pynmd.tools.compute_statistics import find_nearest1d

    
    """
    In: xgrd,ygrd,maxele: model infos
        xhwm,yhwm,elev_hwm:data infos
        flag: how to treat data model comparison
        flag = all :    find nearset grid point
             = valid:   find nearset grid point with non-nan value
             = pos:     find nearset grid point with positive value
             = neg:     find nearset grid point with negative value
        
        
    Retun: model and data vector
    
    
    #  Delta or convert2msl is always for going from vertical datum to msl by an addition to that datum
    # MSL = Vert_datam + convert2msl 

    
    """
    if   flag == 'valid':
        maxe = np.ma.masked_where(maxe==elev_max.fill_value, maxe)
        mask =  maxe.mask   
    elif flag == 'pos':
        mask = [maxe  < 0.0]
    elif flag == 'neg':
        mask = [maxe  > 0.0]
    elif flag == 'all':
        mask = np.isnan(xgrd) 

        #mask = [maxe < -900.0]

    else:
        print ('Choose a valid flag > ')
        print ('flag = all :    find nearset grid point ')
        print ('     = valid:   find nearset grid point with non-nan value')
        print ('     = pos:     find nearset grid point with positive value')
        print ('     = neg:     find nearset grid point with negative valueChoose a valid flag > ')
        sys.exit('ERROR') 
    
    mask  = np.array(mask).squeeze()

    xgrd = xgrd[~mask]
    ygrd = ygrd[~mask]
    maxe = maxe[~mask]
    #
    if convert2msl is not None:
        convert2msl = convert2msl[~mask]
    else:
        convert2msl = np.zeros_like(xgrd)
    #
    if bias_cor is not None:
        bias_cor = bias_cor[~mask]
    else:
        bias_cor = np.zeros_like(xgrd)  

    data  = []
    model = [] 
    prox  = []
    xmodel = []
    ymodel = []
    
    for ip in range(len(xhwm)):
        i,pr  = find_nearest1d(xvec = xgrd,yvec = ygrd,xp = xhwm[ip],yp = yhwm[ip])
        data.append (elev_hwm [ip] + convert2msl[i])
        model.append(maxe[i]+bias_cor[i])
        xmodel.append(xgrd[i].item())
        ymodel.append(ygrd[i].item())
        
        prox.append(pr)
    
    data   = np.array(data ).squeeze()
    model  = np.array(model).squeeze()
    prox   = np.array(prox ).squeeze()
    xmodel = np.array(xmodel).squeeze()
    ymodel = np.array(ymodel).squeeze()
    
    
    #
    #maskf = [model < 0.0]
    #maskf  = np.array(maskf).squeeze()
    #return data[~maskf],model[~maskf],prox[~maskf],xhwm[~maskf],yhwm[~maskf]
    return data,xhwm,yhwm,model,xmodel,ymodel,prox




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

  
##### end of Funcs


# ==============================================================================
# Read fort.42-type ASCII files and save as nc file
# ==============================================================================
def fort42_to_nc(fort42,varnames=['u','v','w'],units=['m s-1','m s-1','m s-1'],
                 longnames=['u velocity component','v velocity component',
                            'vertical velocity component'],
                 ncdate='0001-01-01 00:00:00 UTC',**kwargs):
    """
    Script to read fort.42-type (3d station vector) files and store in a
    netcdf4 file

    PARAMETERS:
    -----------
    fort42 : Path to unit 41 file
    x,y    : Station coordinates
    ncdate : cold start date/time in CF standard: yyyy-MM-dd hh:mm:ss tz

    RETURNS:
    --------
    Netcdf containing
    time     : seconds since beginning of run
    variable : temporal variable (called 'varname') recorded at stations.
               Size: [time,station]
    """

    fobj = open(fort42,'r')

    # Create the file
    if 'savepath' in kwargs:
        ncfile = kwargs['savepath']
    else:
        ncfile = fort42 + '.nc'
    nc = netCDF4.Dataset(ncfile, 'w', format='NETCDF4')

    # Global attributes
    nc.Author = getpass.getuser()
    nc.Created = time.ctime()
    tmpline = fobj.readline()
    nc.description = tmpline[2:34]
    nc.rundes = tmpline[2:34]
    nc.runid = tmpline[36:60]
    nc.model = 'ADCIRC'
    nc.Software = 'Created with Python ' + sys.version
    nc.NetCDF_Lib = str(netCDF4.getlibversion())

    # Record number of time steps and nodes
    tmpline = fobj.readline().split()
#    ndsets = np.int(tmpline[0])
    nsta = np.int(tmpline[1])
    nfen = np.int(tmpline[4])

    # Create dimensions
    nc.createDimension('time',0)           # The unlimited dimension
    nc.createDimension('nsta',nsta)        # Number of station nodes
    nc.createDimension('vnode',nfen)       # Number of vertical nodes/layers

    # Create time vector
    nc.createVariable('time','f8',('time'))
    nc.variables['time'].long_name = 'model time'
    nc.variables['time'].standard_name = 'time'
    nc.variables['time'].units = 'seconds since ' + ncdate
    nc.variables['time'].base_date = ncdate

    # Create and store spatial variables
    nc.createVariable('sigma','f8','vnode')
    nc.variables['sigma'].long_name = 'level of the vertical grid from -1 (bottom) to +1 (surface)'
    nc.variables['sigma'].units = 'dimensionless'

#    nc.createVariable('x','f8','nsta')
#    nc.variables['x'].long_name = 'longitude'
#    nc.variables['x'].units = 'degrees east'
#    nc.variables['x'].positive = 'east'
#    nc.variables['x'][:] = x
#
#    nc.createVariable('y','f8','nsta')
#    nc.variables['y'].long_name = 'latitude'
#    nc.variables['y'].units = 'degrees north'
#    nc.variables['y'].positive = 'north'
#    nc.variables['y'][:] = y

    # Create station variables
    nc.createVariable(varnames[0],'f8',('time','nsta','vnode'))
    nc.variables[varnames[0]].long_name = longnames[0]
    nc.variables[varnames[0]].units = units[0]

    nc.createVariable(varnames[1],'f8',('time','nsta','vnode'))
    nc.variables[varnames[1]].long_name = longnames[1]
    nc.variables[varnames[1]].units = units[1]

    nc.createVariable(varnames[2],'f8',('time','nsta','vnode'))
    nc.variables[varnames[2]].long_name = longnames[2]
    nc.variables[varnames[2]].units = units[2]

    # Store time-series variables
#    for tt in range(ntsteps):    #this would work if there are no errors with the simulation
#        nc.variables['time'][tt] = np.float64(fobj.readline().split()[0])
#        for aa in range(nsta):
#            nc.variables[varname][tt,aa] = np.float64(fobj.readline().split()[1])
    tt = -1
    line = fobj.readline()
    nc.variables['sigma'][:] = np.asarray([float(b) for b in line.split()[2:nfen*3:3]])
    while line:
        tt += 1
        nc.variables['time'][tt] = np.float64(line.split()[0])
        line = fobj.readline()
        for aa in range(nsta):
            if line:
                nc.variables[varnames[0]][tt,aa,:] = np.float64(line.split()[1:nfen*3:3])
                nc.variables[varnames[1]][tt,aa,:] = np.float64(line.split()[2:nfen*3:3])
                nc.variables[varnames[2]][tt,aa,:] = np.float64(line.split()[3:nfen*3+1:3])
                line = fobj.readline()

    # All done here
    fobj.close()
    nc.close()


# ==============================================================================
# Read fort.61-type ASCII files and save as nc file
# ==============================================================================
def fort61_to_nc(fort61,staname,x,y,varname='zeta',
                 longname='water surface elevation above geoid',
                 varunits='m',ncdate='0001-01-01 00:00:00 UTC',**kwargs):
    """
    Script to read fort.61-type (station scalar) files and store in a netcdf4 file

    PARAMETERS:
    -----------
    fort61 : Path to unit 61 file
    x,y    : Station coordinates
    ncdate : cold start date/time in CF standard: yyyy-MM-dd hh:mm:ss tz

    RETURNS:
    --------
    Netcdf containing
    time     : seconds since beginning of run
    variable : temporal variable (called 'varname') recorded at stations.
               Size: [time,station]
    """

    fobj = open(fort61,'r')

    # Create the file
    if 'savepath' in kwargs:
        ncfile = kwargs['savepath']
    else:
        ncfile = fort61 + '.nc'
    nc = netCDF4.Dataset(ncfile, 'w', format='NETCDF4')

    # Global attributes
    nc.Author = getpass.getuser()
    nc.Created = time.ctime()
    tmpline = fobj.readline()
    nc.description = tmpline[2:34]
    nc.rundes = tmpline[2:34]
    nc.runid = tmpline[36:60]
    nc.model = 'ADCIRC'
    nc.Software = 'Created with Python ' + sys.version
    nc.NetCDF_Lib = str(netCDF4.getlibversion())

    # Record number of time steps and station
    tmpline = fobj.readline().split()
#    ntsteps = np.int(tmpline[0])
    nsta = np.int(tmpline[1])

    # Create dimensions
    nc.createDimension('time',0)           # The unlimited dimension
    nc.createDimension('nsta',nsta)        # Number of stations
    nc.createDimension('namelen',50)       # Length of station names

    # Create time vector
    nc.createVariable('time','f8',('time'))
    nc.variables['time'].long_name = 'model time'
    nc.variables['time'].standard_name = 'time'
    nc.variables['time'].units = 'seconds since ' + ncdate
    nc.variables['time'].base_date = ncdate

    # Create and store spatial variables
    nc.createVariable('station_name','S1',('nsta','namelen'))
    nc.variables['station_name'].long_name = 'station name'
    nc.variables['station_name'][:] = netCDF4.stringtochar(np.array(staname,dtype='S50'))
    nc.createVariable('x','f8','nsta')
    nc.variables['x'].long_name = 'longitude'
    nc.variables['x'].units = 'degrees east'
    nc.variables['x'].positive = 'east'
    nc.variables['x'][:] = x
    nc.createVariable('y','f8','nsta')
    nc.variables['y'].long_name = 'latitude'
    nc.variables['y'].units = 'degrees north'
    nc.variables['y'].positive = 'north'
    nc.variables['y'][:] = y

    # Create station variable
    nc.createVariable(varname,'f8',('time','nsta'))
    nc.variables[varname].long_name = longname
    nc.variables[varname].units = varunits

    # Store time-series variables
#    for tt in range(ntsteps):    #this would work if there are no errors with the simulation
#        nc.variables['time'][tt] = np.float64(fobj.readline().split()[0])
#        for aa in range(nsta):
#            nc.variables[varname][tt,aa] = np.float64(fobj.readline().split()[1])
    tt = -1
    line = fobj.readline()
    while line:
        tt += 1
        nc.variables['time'][tt] = np.float64(line.split()[0])
        line = fobj.readline()
        for aa in range(nsta):
            if line:
                nc.variables[varname][tt,aa] = np.float64(line.split()[1])
                line = fobj.readline()

    # All done here
    fobj.close()
    nc.close()


# ==============================================================================
# Read fort.63-type ASCII files and save as nc file
# ==============================================================================
def fort63_to_nc(fort63,varname='zeta',
                 longname='water surface elevation above geoid',
                 varunits='m',ncdate='0001-01-01 00:00:00 UTC',**kwargs):
    """
    Script to read fort.63-type (scalar) files and store in a netcdf4 file

    PARAMETERS:
    -----------
    fort63: Path to fort.63-type file
    ncdate : cold start date/time in CF standard: yyyy-MM-dd hh:mm:ss tz

    RETURNS:
    --------
    Netcdf containing
    time     : seconds since beginning of run
    variable : temporal variable (called 'varname') recorded at the nodes of
               an unstructured grid. Size: [time,nodes]
    """

    fobj = open(fort63,'r')

    # Create the file
    if 'savepath' in kwargs:
        ncfile = kwargs['savepath']
    else:
        ncfile = fort63 + '.nc'
    nc = netCDF4.Dataset(ncfile, 'w', format='NETCDF4')

    # Global attributes
    nc.Author = getpass.getuser()
    nc.Created = time.ctime()
    tmpline = fobj.readline()
    nc.description = tmpline[2:34]
    nc.rundes = tmpline[2:34]
    nc.runid = tmpline[36:60]
    nc.model = 'ADCIRC'
    nc.Software = 'Created with Python ' + sys.version
    nc.NetCDF_Lib = str(netCDF4.getlibversion())

    # Record number of time steps and nodes
    tmpline = fobj.readline().split()
#    ntsteps = np.int(tmpline[0])
    nodes = np.int(tmpline[1])

    # Create dimensions
    nc.createDimension('node',nodes)       # Number of nodes
    nc.createDimension('time',0)           # The unlimited dimension

    # Create time vector
    nc.createVariable('time','f8',('time'))
    nc.variables['time'].long_name = 'model time'
    nc.variables['time'].standard_name = 'time'
    nc.variables['time'].units = 'seconds since ' + ncdate
    nc.variables['time'].base_date = ncdate

    # Create the rest of the variables
    nc.createVariable(varname,'f8',('time','node'))
    nc.variables[varname].long_name = longname
    nc.variables[varname].units = varunits

    # Store time-series variables at nodes
#    for tt in range(ntsteps):    #this would work if there are no errors with the simulation
#        nc.variables['time'][tt] = np.float64(fobj.readline().split()[0])
#        for aa in range(nodes):
#            nc.variables[varname][tt,aa] = np.float64(fobj.readline().split()[1])
    tt = -1
    line = fobj.readline()
    while line:
        tt += 1
        nc.variables['time'][tt] = np.float64(line.split()[0])
        line = fobj.readline()
        for aa in range(nodes):
            if line:
                nc.variables[varname][tt,aa] = np.float64(line.split()[1])
                line = fobj.readline()


    # All done here
    fobj.close()
    nc.close()


# ==============================================================================
# Read fort.64-type ASCII files and save as nc file
# ==============================================================================
def fort64_to_nc(fort64,varname_xy=['u-vel','v-vel'],
                 longname='water column vertically averaged',
                 varunits='m s-1',ncdate='0001-01-01 00:00:00 UTC',**kwargs):
    """
    Script to read fort.64-type (vector) files and store in a netcdf4 file

    PARAMETERS:
    -----------
    fort64: Path to unit 64 file
    ncdate : cold start date/time in CF standard: yyyy-MM-dd hh:mm:ss tz

    RETURNS:
    --------
    Netcdf containing
    time         : seconds since beginning of run
    x,y variable : temporal variables (name provided in 'xy_varname') recorded
                   at the nodes of an unstructured grid. Size: [time,nodes]
    """

    fobj = open(fort64,'r')

    # Create the file
    if 'savepath' in kwargs:
        ncfile = kwargs['savepath']
    else:
        ncfile = fort64 + '.nc'
    nc = netCDF4.Dataset(ncfile, 'w', format='NETCDF4')

    # Global attributes
    nc.Author = getpass.getuser()
    nc.Created = time.ctime()
    tmpline = fobj.readline()
    nc.description = tmpline[2:34]
    nc.rundes = tmpline[2:34]
    nc.runid = tmpline[36:60]
    nc.model = 'ADCIRC'
    nc.Software = 'Created with Python ' + sys.version
    nc.NetCDF_Lib = str(netCDF4.getlibversion())

    # Record number of time steps and nodes
    tmpline = fobj.readline().split()
#    ntsteps = np.int(tmpline[0])
    nodes = np.int(tmpline[1])

    # Create dimensions
    nc.createDimension('time',0)           # The unlimited dimension
    nc.createDimension('node',nodes)       # Number of nodes

    # Create time vector
    nc.createVariable('time','f8',('time'))
    nc.variables['time'].long_name = 'model time'
    nc.variables['time'].standard_name = 'time'
    nc.variables['time'].units = 'seconds since ' + ncdate
    nc.variables['time'].base_date = ncdate

    # Create the rest of the variables
    nc.createVariable(varname_xy[0],'f8',('time','node'))
    nc.variables[varname_xy[0]].long_name = longname + ' e/w velocity'
    nc.variables[varname_xy[0]].units = varunits

    nc.createVariable(varname_xy[1],'f8',('time','node'))
    nc.variables[varname_xy[1]].long_name = longname + ' n/s velocity'
    nc.variables[varname_xy[1]].units = varunits

    # Store variables
#    for tt in range(ntsteps):    #this would work if there are no errors with the simulation
#        nc.variables['time'][tt] = np.float64(fobj.readline().split()[0])
#        for aa in range(nodes):
#            tmpline = fobj.readline().split()
#            nc.variables[varname_xy[0]][tt,aa] = np.float64(tmpline[1])
#            nc.variables[varname_xy[1]][tt,aa] = np.float64(tmpline[2])
    tt = -1
    line = fobj.readline()
    while line:
        tt += 1
        nc.variables['time'][tt] = np.float64(line.split()[0])
        line = fobj.readline()
        for aa in range(nodes):
            nc.variables[varname_xy[0]][tt,aa] = np.float64(line.split()[1])
            nc.variables[varname_xy[1]][tt,aa] = np.float64(line.split()[2])
            line = fobj.readline()

    # All done here
    fobj.close()
    nc.close()


# ==============================================================================
# Read Global Maximum and Minimum ASCII files and save as nc file
# ==============================================================================
def max63_to_nc(max63,varname='zeta',
                 longname='water surface elevation above geoid',
                 varunits='m',ncdate='0001-01-01 00:00:00 UTC',**kwargs):
    """
    Script to read max63-type (max-min) files and store in a netcdf4 file

    PARAMETERS:
    -----------
    max63: Path to max63-type file
    ncdate : cold start date/time in CF standard: yyyy-MM-dd hh:mm:ss tz

    RETURNS:
    --------
    Netcdf containing
    time         : seconds since beginning of run
    x,y variable : temporal variables (name provided in 'xy_varname') recorded
                   at the nodes of an unstructured grid. Size: [time,nodes]
    """

    fobj = open(max63,'r')

    # Create the file
    if 'savepath' in kwargs:
        ncfile = kwargs['savepath']
    else:
        ncfile = max63 + '.nc'
    nc = netCDF4.Dataset(ncfile, 'w', format='NETCDF4')

    # Global attributes
    nc.Author = getpass.getuser()
    nc.Created = time.ctime()
    tmpline = fobj.readline()
    nc.description = tmpline[2:34]
    nc.rundes = tmpline[2:34]
    nc.runid = tmpline[36:60]
    nc.model = 'ADCIRC'
    nc.Software = 'Created with Python ' + sys.version
    nc.NetCDF_Lib = str(netCDF4.getlibversion())

    # Record number of time steps and nodes
    nodes = np.int(fobj.readline().split()[1])

    # Create dimensions
    nc.createDimension('time',0)           # The unlimited dimension
    nc.createDimension('node',nodes)       # Number of nodes

    # Create time vector
    nc.createVariable('time','f8',('time'))
    nc.variables['time'].long_name = 'model time'
    nc.variables['time'].standard_name = 'time'
    nc.variables['time'].units = 'seconds since ' + ncdate
    nc.variables['time'].base_date = ncdate

    # Create the rest of the variables
    nc.createVariable(varname+'_max','f8','node')
    nc.variables[varname+'_max'].long_name = 'maximum ' + longname
    nc.variables[varname+'_max'].units = varunits

    nc.createVariable('time_of_'+varname+'_max','f8','node')
    nc.variables['time_of_'+varname+'_max'].long_name = 'time of maximum ' + longname
    nc.variables['time_of_'+varname+'_max'].units = 'sec'

    # Store last time-step in run for reference
    nc.variables['time'][0] = np.float64(fobj.readline().split()[0])

    # Store max variable observations
    for aa in range(nodes):
        nc.variables[varname+'_max'][aa] = np.float64(fobj.readline().split()[1])

    # Skip repeated time
    fobj.readline()

    for aa in range(nodes):
        nc.variables['time_of_'+varname+'_max'][aa] = np.float64(fobj.readline().split()[1])

    # All done here
    fobj.close()
    nc.close()


# ==============================================================================
# Read ncfile of degree scalars and return nc file of vector components
# ==============================================================================
def ncdirs_to_vec(ncfolder,varname='swan_DIR',
                  source_convention='directions to',
                  target_convention='directions to',
                  **kwargs):
    """
    Script to read a netcdf4 file with wave directions given as scalars
    (in degrees) and write a corresponding netcdf4 file with directions
    given as (u,v) vectors of magnitude 1 in their mean direction.
    Useful to show direction vectors in plots.

    PARAMETERS:
    -----------
    ncfolder: Path to folder containing direction netcdf file. Vector file
              will be stored here as well unless a different 'savepath' path is
              declared in the kwargs
    source,target convention: indicate whether directions are definde as 'to'
                              or 'from' in the source file
    RETURNS:
    --------
    Netcdf containing
    time         : seconds since beginning of run
    x,y variable : temporal variables (name provided in 'xy_varname') recorded
                   at the nodes of an unstructured grid. Size: [time,nodes]
    u,v directions
    NOTES:
    ------
    1) North is defined as a zero degree direction. If direction conventions do
       not match between source and target files, directions are inverted.
    """

    nc_dir = netCDF4.Dataset(ncfolder+varname+'.63.nc', 'r', format='NETCDF4')

    # Create the file
    if 'savepath' in kwargs:
        ncfile = kwargs['savepath']
    else:
        ncfile = ncfolder+varname+'vec.63.nc'
    nc_vec = netCDF4.Dataset(ncfile, 'w', format='NETCDF4')

    # Global attributes
    nc_vec.Author = getpass.getuser()
    nc_vec.Created = time.ctime()
    nc_vec.Software = 'Created with Python ' + sys.version
    nc_vec.NetCDF_Lib = str(netCDF4.getlibversion())

    # Copy additional global attributes from source
    nc_vec.setncatts({a:nc_dir.getncattr(a) for a in nc_dir.ncattrs() if
                      a not in ['creation_date','modification_date','host',
                                'convention','contact']})

    # Create dimensions
    nc_vec.createDimension('time',0)       # The unlimited dimension
    nc_vec.createDimension('node', len(nc_dir.dimensions['node'])) # Number of nodes


    # Copy variables
    for name, var in nc_dir.variables.items():
        if name in ['time','x','y']:
            # Create select vars
            nc_vec.createVariable(name, var.dtype, var.dimensions)
            # Copy the variable attributes
            nc_vec.variables[name].setncatts({a:var.getncattr(a) for a in var.ncattrs()})
            # Copy the variables values
            nc_vec.variables[name][:] = nc_dir.variables[name][:]

    # Create the rest of the variables
    nc_vec.createVariable(varname+'_u','f8',('time','node'))
    nc_vec.variables[varname+'_u'].long_name = 'e/w direction'
    nc_vec.variables[varname+'_u'].units = nc_dir[varname].units
    nc_vec.variables[varname+'_u'].convention = target_convention

    nc_vec.createVariable(varname+'_v','f8',('time','node'))
    nc_vec.variables[varname+'_v'].long_name = 'n/s direction'
    nc_vec.variables[varname+'_v'].units = nc_dir[varname].units
    nc_vec.variables[varname+'_v'].convention = target_convention

    for aa in range(len(nc_vec['time'])):
        if source_convention != target_convention:
            dirs = gangles.wrapto360(nc_dir[varname][aa,:].data + 180) * np.pi/180
        else:
            dirs = nc_dir[varname][aa,:].data  * np.pi/180
        nc_vec.variables[varname+'_u'][aa,:] = np.sin(dirs)
        nc_vec.variables[varname+'_v'][aa,:] = np.cos(dirs)

    # All done here
    nc_dir.close()
    nc_vec.close()


# ==============================================================================
# Read all ASCII files and save as nc file
# ==============================================================================
def all_ascii2nc(runFld,ncFld,ncdate='0001-01-01 00:00:00 UTC'):
    """
    Reads (most*) ADCIRC ASCII files and stores them in netcdf4 files

    Input:
    ------
    runFld: folder containing all ascii files
    ncFld:  folder where to save all nc files

    Notes:
    ------
    *Still working on converting some adcirc files/formats
    """
    # Read input file to retrieve station information
    fort15 = adcpre.readsta_fort15(runFld +'fort.15')

    # Unstructured grid + bathy ------------------------------------------------
    if os.path.exists(runFld+'fort.14'):
        if not os.path.exists(ncFld+'fort.14.nc'):
            print('Creating fort.14.nc')
            adcpre.fort14_to_nc(runFld + 'fort.14',savepath=ncFld+'fort.14.nc')

    # fort.61-type files (scalars) ---------------------------------------------
    # Water surface elevation
    if os.path.exists(runFld+'fort.61'):
        if not os.path.exists(ncFld+'fort.61.nc'):
            print('Creating fort.61.nc')
            fort61_to_nc(runFld + 'fort.61',fort15['nameel'],fort15['xel'],
                         fort15['yel'],ncdate=ncdate,
                         savepath=ncFld+'fort.61.nc')

    # fort.63-type files (scalars) ---------------------------------------------
    # Water surface elevation
    if os.path.exists(runFld+'fort.63'):
        if not os.path.exists(ncFld+'fort.63.nc'):
            print('Creating fort.63.nc')
            fort63_to_nc(runFld + 'fort.63',ncdate=ncdate,
                         savepath=ncFld+'fort.63.nc')
    # Atmospherec pressure
    if os.path.exists(runFld+'fort.73'):
        if not os.path.exists(ncFld+'fort.73.nc'):
            print('Creating fort.73.nc')
            fort63_to_nc(runFld + 'fort.73',varname='pressure',
                         longname='air pressure at sea level',
                         varunits='meters of waver',ncdate=ncdate,
                         savepath=ncFld+'fort.73.nc')
    # Significant wave height
    if os.path.exists(runFld+'swan_HS.63'):
        if not os.path.exists(ncFld+'swan_HS.63.nc'):
            print('Creating swan_HS.63.nc')
            fort63_to_nc(runFld + 'swan_HS.63',varname='swan_HS',
                         longname='significant wave height',
                         varunits='m',ncdate=ncdate,
                         savepath=ncFld+'swan_HS.63.nc')
    if os.path.exists(runFld+'swan_HS_max.63'):
        if not os.path.exists(ncFld+'swan_HS_max.63.nc'):
            print('Creating swan_HS_max.63.nc')
            fort63_to_nc(runFld + 'swan_HS_max.63',varname='swan_HS_max',
                        longname='maximum significant wave height',
                        varunits='m',ncdate=ncdate,
                        savepath=ncFld+'swan_HS_max.63.nc')
    # Mean wave direction
    if os.path.exists(runFld+'swan_DIR.63'):
        if not os.path.exists(ncFld+'swan_DIR.63.nc'):
            print('Creating swan_DIR.63.nc')
            fort63_to_nc(runFld + 'swan_DIR.63',varname='swan_DIR',
                         longname='mean wave direction',
                         varunits='degrees',ncdate=ncdate,
                         savepath=ncFld+'swan_DIR.63.nc')
    if os.path.exists(runFld+'swan_DIR_max.63'):
        if not os.path.exists(ncFld+'swan_DIR_max.63.nc'):
            print('Creating swan_DIR_max.63.nc')
            fort63_to_nc(runFld + 'swan_DIR_max.63',varname='swan_DIR_max',
                         longname='maximum mean wave direction',
                         varunits='degrees',ncdate=ncdate,
                         savepath=ncFld+'swan_DIR_max.63.nc')
    # Mean absolute wave period
    if os.path.exists(runFld+'swan_TMM10.63'):
        if not os.path.exists(ncFld+'swan_TMM10.63.nc'):
            print('Creating swan_TMM10.63.nc')
            fort63_to_nc(runFld + 'swan_TMM10.63',varname='swan_TMM10',
                         longname='mean absolute wave period',
                         varunits='s',ncdate=ncdate,
                         savepath=ncFld+'swan_TMM10.63.nc')
    if os.path.exists(runFld+'swan_TMM10_max.63'):
        if not os.path.exists(ncFld+'swan_TMM10_max.63.nc'):
            print('Creating swan_TMM10_max.63.nc')
            fort63_to_nc(runFld + 'swan_TMM10_max.63',varname='swan_TMM10_max',
                         longname='maximum TMM10 mean wave period',
                         varunits='s',ncdate=ncdate,
                         savepath=ncFld+'swan_TMM10_max.63.nc')
    # Smoothed peak period
    if os.path.exists(runFld+'swan_TPS.63'):
        if not os.path.exists(ncFld+'swan_TPS.63.nc'):
            print('Creating swan_TPS.63.nc')
            fort63_to_nc(runFld + 'swan_TPS.63',varname='swan_TPS',
                         longname='smoothed peak period',
                         varunits='s',ncdate=ncdate,
                         savepath=ncFld+'swan_TPS.63.nc')
    if os.path.exists(runFld+'swan_TPS_max.63'):
        if not os.path.exists(ncFld+'swan_TPS_max.63.nc'):
            print('Creating swan_TPS_max.63.nc')
            fort63_to_nc(runFld + 'swan_TPS_max.63',varname='swan_TPS_max',
                         longname='maximum smoothed peak period',
                         varunits='s',ncdate=ncdate,
                         savepath=ncFld+'swan_TPS_max.63.nc')

    # fort.64-type files (vectors) ---------------------------------------------
    # Depth average velocity
    if os.path.exists(runFld+'fort.64'):
        if not os.path.exists(ncFld+'fort.64.nc'):
            print('Creating fort.64.nc')
            fort64_to_nc(runFld + 'fort.64',ncdate=ncdate,savepath=ncFld+'fort.64.nc')
    # Wind stress or velocity
    if os.path.exists(runFld+'fort.74'):
        if not os.path.exists(ncFld+'fort.74.nc'):
            print('Creating fort.74.nc')
            fort64_to_nc(runFld + 'fort.74',varname_xy=['windx','windy'],
                         longname='wind',ncdate=ncdate,savepath=ncFld+'fort.74.nc')

    # max.63-type files (max/min values) ---------------------------------------
    # Water surface elevation
    if os.path.exists(runFld+'maxele.63'):
        if not os.path.exists(ncFld+'maxele.63.nc'):
            print('Creating maxele.63.nc')
            max63_to_nc(runFld + 'maxele.63',ncdate=ncdate,savepath=ncFld+'maxele.63.nc')
    if os.path.exists(runFld+'maxvel.63'):
        if not os.path.exists(ncFld+'maxvel.63.nc'):
            print('Creating maxvel.63.nc')
            max63_to_nc(runFld + 'maxvel.63',varname='vel',
                        longname='water velocity',
                        varunits='m s-1',ncdate=ncdate,savepath=ncFld+'maxvel.63.nc')
    if os.path.exists(runFld+'maxwvel.63'):
        if not os.path.exists(ncFld+'maxwvel.63.nc'):
            print('Creating maxwvel.63.nc')
            max63_to_nc(runFld + 'maxwvel.63',varname='wind',
                        longname='wind velocity',
                        varunits='m s-1',ncdate=ncdate,savepath=ncFld+'maxwvel.63.nc')

    # fort.42-type files (3D vectors) ------------------------------------------
    if os.path.exists(runFld+'fort.42'):
        if not os.path.exists(ncFld+'fort.42.nc'):
            print('Creating fort.42.nc')
            fort42_to_nc(runFld + 'fort.42',ncdate=ncdate,savepath=ncFld+'fort.42.nc')
