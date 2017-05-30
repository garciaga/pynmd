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
    fname =  os.path.abspath(DirName + 'fort.63.nc')  
    #print fname
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
    fname =  os.path.abspath(DirName + 'fort.63.nc'  )
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


  
##### end of Funcs
