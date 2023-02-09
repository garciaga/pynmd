# -*- coding: utf-8 -*-
"""
A series of tools to process and prepare adcirc input

Authors:
-------
Fadia Ticona Rollano
    PNNL Marine Sciences Laboratory

Log of edits:
-------------
December 2019 - Created module
    Fadia Ticona Rollano

"""
from __future__ import division,print_function

__author__ = "Saeed Moghimi"
__copyright__ = "Copyright 2015, Portland State University"
__license__ = "GPL"
__version__ = "1.0"
__email__ = "moghimis@gmail.com"



import os#,glob
import sys,time
#sys.path.append(os.getcwd())
import getpass
import numpy as np
import netCDF4

# Custom paths
import pynmd.models.tools.unstructured as ustr



def fort14_update_bathy(data, path=None, file_name='fort.14'):
    """
    Write out bathymetry in data to ``fort.14`` formated file, by updating
    path/file_name accordingly

    :type data: dict(x,y,bathy) in vector format arranged based on node numbers
    :param data: python object to save the ``fort.14`` data to
    :type bathymetry: array or None
    :param bathymetry: if None then use the bathymetry in ``data``
    :type path: string or None
    :param path: path to the``fort.14`` fortmatted file
    :param string  file_name: file name

    """
    if path is None:
        path = os.getcwd()

    file_name = os.path.join(path, file_name)
    tmp = os.path.join(path, 'temp.14')

    # this currrently uses the present fort.14 as a template for formatting
    # purposes
    bathymetry = data['depth']
    x          = data['x']
    y          = data['y']
    nodes      = 1 + range(len(x))  
    # pylint: disable=C0103
    with open(file_name, 'r') as f, open(tmp, 'w') as fw:
        fw.write(f.readline())
        fw.write(f.readline())
        for i in range(len(x)):
            # pylint: disable=C0103
            fw.write('{:<7d} {:9.8E} {:9.8E} {:7.2f}\n'.format(nodes[i], x[i], y[i],
                                                               bathymetry[i]))
            f.readline()
        for line in f:
            fw.write(line)
    #os.rename(tmp, file_name)




# ==============================================================================
# Read fort.14 ASCII file and save as nc file (not necessary to run adcirc)
# ==============================================================================
def fort14_to_nc(fort14,**kwargs):
    """ 
    Script to read an unstructured grid in fort.14 format and store it in a 
    netcdf4 file

    PARAMETERS:
    -----------
    fort14: Path to fort14 file
    savepath(optional)

    RETURNS:
    --------
    Netcdf containing
    x,y:     node coordinates
    element: triangular elements
    nbdv:    node numbers on each elevation specified boundary segment
    nbvv:    node numbers on normal flow boundary segment
    depth:   depth at each node
    
    """
    grid = ustr.read_fort14(fort14)
    
    # Create the file and add global attributes
    if 'savepath' in kwargs:
        ncfile = kwargs['savepath']
    else:
        ncfile = fort14 + '.nc'
    nc = netCDF4.Dataset(ncfile, 'w', format='NETCDF4')
    
    # Global attributes 
    nc.Author = getpass.getuser()
    nc.Created = time.ctime()  
    nc.Software = 'Created with Python ' + sys.version
    nc.NetCDF_Lib = str(netCDF4.getlibversion())
    
    # Get record lengths
    node = grid['x'].size
    nele = grid['triang'].shape[0]
    nvertex = grid['triang'].shape[1]
#    nope = grid['nope']
    neta = grid['neta']    
#    nbou = grid['nbou']
    nvel = grid['nvel']
    
    # Create dimensions
    nc.createDimension('node',node)       # Number of nodes
    nc.createDimension('nele',nele)       # Number of elements
    nc.createDimension('nvertex',nvertex) # Number of vertex
    nc.createDimension('neta',neta)       # Number of elevation boundary nodes
    nc.createDimension('nvel',nvel)       # Number of normal flow boundary nodes
    
    # Create variables
    nc.createVariable('x','f8','node')
    nc.variables['x'].long_name = 'longitude'
    nc.variables['x'].units = 'degrees_east'
    nc.variables['x'].positive = 'east'
    nc.variables['x'][:] = np.float64(grid['x'])
    
    nc.createVariable('y','f8','node')
    nc.variables['y'].long_name = 'latitude'
    nc.variables['y'].units = 'degrees_north'
    nc.variables['y'].positive = 'north'
    nc.variables['y'][:] = np.float64(grid['y'])
    
    nc.createVariable('element','i4',('nele','nvertex'))
    nc.variables['element'].long_name = 'element'
    nc.variables['element'].units = 'nondimensional'
    # Add 1 back to node indexing to match standard adcirc files convention
    nc.variables['element'][:,:] = np.int32(grid['triang']+1)
    
    nc.createVariable('nbdv','i4','neta')
    nc.variables['nbdv'].long_name = 'node numbers on each elevation specified boundary segment'
    nc.variables['nbdv'].units = 'nondimensional'
    nc.variables['nbdv'][:] = np.int32(np.hstack(grid['nbdv']))
    
    nc.createVariable('nbvv','i4','nvel')
    nc.variables['nbvv'].long_name = 'node numbers on normal flow boundary segment'
    nc.variables['nbvv'].units = 'nondimensional'
    nc.variables['nbvv'][:] = np.int32(np.hstack(grid['nbvv']))
    
    nc.createVariable('depth','f8','node')
    nc.variables['depth'].long_name = 'distance below geoid'
    nc.variables['depth'].units = 'm'
    nc.variables['depth'][:] = np.float64(grid['z'])
        
    # All done here
    nc.close()

# ==============================================================================
# Read station information from fort.15 ASCII files
# ==============================================================================
def readsta_fort15(fort15):
    """ 
    Script to read fort.15 input file and store station information

    PARAMETERS:
    -----------
    fort15: Path to fort15 file

    RETURNS:
    --------
    Dictionary containing the number, coordinates, and name of stations 
    as follows:
        
    nstae,xel,yel,nameel:  water surface elevation 
    nstav,xev,yev,nameev:  velocity 
    nstac,xec,yec,nameec:  concentrations
    nstam,xem,yem,nameem:  meterology
    
    """
    
    # Get "neta" (total number of elevation boundary nodes) from grid file
    fort14 = os.path.dirname(fort15)+'/fort.14'
    if not os.path.isfile(fort14):
        print ('fort.14 is not loacated in the same directory as fort.15')
        fort14 = input("Please enter path to fort.14 file:\n")
    try:
        neta = ustr.read_fort14(fort14)['neta']
    except IOError:
        print("Grid file (fort.14) not accessible")
            
    # Open Model Parameter and Periodic Boundary Condition File
    fobj = open(fort15,'r')
    
    # Skip lines (checks are needed because the number of lines varies)
    for ii in range(7):
        fobj.readline() #RUENDES -- ICS
    im = np.int(fobj.readline().split()[0])
    if im==21:
        fobj.readline() #IDEN
    for ii in range(4):
        fobj.readline() #NOLIBF -- NOLICAT
    tmpline = np.int(fobj.readline().split()[0]) #NWP
    for ii in range(tmpline+2):
        fobj.readline() #AttrName -- NTIP
    nws = np.int(fobj.readline().split()[0])
    fobj.readline() #NRAMP
    fobj.readline() #G
    tmpline = np.float64(fobj.readline().split()[0]) #TAU0
    if tmpline==-5.:
        fobj.readline() #Tau0FullDomainMin, Tau0FullDomainMax
    for ii in range(3):
        fobj.readline() #DTDP -- REFTIME
    if nws not in [0,1,9,11]:
        fobj.readline() #WTIMINC
    for ii in range(6):
        fobj.readline() #RNDAY -- TAU/CD/CF...FGMMA
    if im==0 or im==1 or im==2 or im==10:
        fobj.readline() #ESLM/ESLM,ESLC
    fobj.readline() #CORI
    ntif = np.int(fobj.readline().split()[0])
    for ii in range(ntif*2):
        fobj.readline() #-TIPOTAG -- ...,FACET(k)
    nbfr = np.int(fobj.readline().split()[0])
    for ii in range(nbfr*3 + nbfr*neta):
        fobj.readline() #BOUNTAG -- ...,EFA(k,j)
    fobj.readline() #ANGINN
    # Our runs don't consider river or ocean inflow (IBTYPE 2,12,22,32,52)
    # so I skip trying to read those lines (NFFR--...,ENPH(k,j)) here
    
    # Elevation recording stations (fort.61 output)
    fobj.readline() #NOUTE, TOUTSE, TOUTFE, NSPOOLE
    nstae = np.int(fobj.readline().split()[0])
    xel = np.zeros(nstae,)
    yel = np.zeros(nstae,)
    nameel = ['' for x in range(nstae)]
    for ii in range(nstae):
        tmpline = fobj.readline().split()
        xel[ii] = np.float64(tmpline[0])
        yel[ii] = np.float64(tmpline[1])
        n = ''
        for s in tmpline[2:]:
            n = n+s
        nameel[ii] = n[:50].ljust(50)
    dix = {'nstae':nstae,'xel':xel,'yel':yel,'nameel':nameel}
    
    # Velocity recording stations (fort.62 output)
    fobj.readline() #NOUTV, TOUTSV,TOUTFV, NSPOOLV
    nstav = np.int(fobj.readline().split()[0])
    xev = np.zeros(nstav,)
    yev = np.zeros(nstav,)
    nameev = ['' for x in range(nstav)]
    for ii in range(nstav):
        tmpline = fobj.readline().split()
        xev[ii] = np.float64(tmpline[0])
        yev[ii] = np.float64(tmpline[1])
        n = ''
        for s in tmpline[2:]:
            n = n+s
        nameev[ii] = n[:50].ljust(50)
    dix.update({'nstav':nstav,'xev':xev,'yev':yev,'nameev':nameev})
    
    # Concentration recording stations (fort.91 output)
    if im==10:
        fobj.readline() #NOUTC, TOUTSC, TOUTFC, NSPOOLC
        nstac = np.int(fobj.readline().split()[0])
        xec = np.zeros(nstac,)
        yec = np.zeros(nstac,)
        nameec = ['' for x in range(nstac)]
        for ii in range(nstac):
            tmpline = fobj.readline().split()
            xec[ii] = np.float64(tmpline[0])
            yec[ii] = np.float64(tmpline[1])
            n = ''
            for s in tmpline[2:]:
                n = n+s
            nameec[ii] = n[:50].ljust(50)
        dix.update({'nstac':nstac,'xev':xec,'yev':yec,'nameec':nameec})
        
    # Meterology recording stations (units 71&72 output)
    if nws != 0:
        fobj.readline() #NOUTM, TOUTSM, TOUTFM, NSPOOLM
        nstam = np.int(fobj.readline().split()[0])
        xem = np.zeros(nstam,)
        yem = np.zeros(nstam,)
        nameem = ['' for x in range(nstam)]
        for ii in range(nstam):
            tmpline = fobj.readline().split()
            xem[ii] = np.float64(tmpline[0])
            yem[ii] = np.float64(tmpline[1])
            n = ''
            for s in tmpline[2:]:
                n = n+s
            nameem[ii] = n[:50].ljust(50)
        dix.update({'nstam':nstam,'xem':xem,'yem':yem,'nameem':nameem})
    
    # Close ASCII file
    fobj.close()
    
    # Generate output
    return dix

# ==============================================================================
# Read fort.22 ASCII file and save as nc file (not necessary to run adcirc)
# ==============================================================================
def fort22nws5_to_nc(fort22,wtiminc,ncdate='0001-01-01 00:00:00 UTC',**kwargs):
    """ 
    Script to read fort.22 input file corresponding to a simulation of type
    NWS=5 (wind velocity, every node, every WTIMINC) and store it in a 
    netcdf4 file

    PARAMETERS:
    -----------
    fort22:   path to fort22 file
    wtiminc:  time increment between meteorological forcing data sets (in sec)
    ncdate:   cold start date/time in CF standard: yyyy-MM-dd hh:mm:ss tz
    nodes :   total number of nodes (optional). If not provided, code will 
              attempt to read fort14 to obtain the node count
    fort14:   path to fort14 file (optional, checked only if 'nodes' is not
              provided, assumed to be in the same directory as fort22)
    savepath: path to output nc file
    
    RETURNS:
    --------
    Netcdf containing values at each corresponding node
    wvx,wvy: applied horizontal wind velocity in the x,y directions. An 
             oceanographic convention is used where velocity is positive when
             it is blowing toward positive coordinate directions. 
    prn:     applied atmospheric pressure at the free surface
    
    """
    
    # Get total number of elevation boundary nodes from grid file
    if 'nodes' in kwargs:
        nodes = kwargs['nodes']
    else:
        if 'fort14' in kwargs:
            fort14 = kwargs['fort14']
        else:
            fort14 = os.path.dirname(fort22)+'/fort.14'
        if not os.path.isfile(fort14):
            print ('Grid file (fort.14) was not specified and is not loacated '+ 
                   'in the same directory as fort.22')
            fort14 = input('Please enter path to fort.14 file (or quit and '+ 
                           'provide total number of nodes):\n')
        try:
            with open(fort14,'r') as fobj:
                fobj.readline() # Header line
                nodes = np.int(fobj.readline().split()[1])
        except IOError:
            print("Grid file not accessible")
        
    # Create the file and add global attributes
    if 'savepath' in kwargs:
        ncfile = kwargs['savepath']
    else:
        ncfile = fort22 + '.nc'
    nc = netCDF4.Dataset(ncfile, 'w', format='NETCDF4')
    
    # Global attributes 
    nc.Author = getpass.getuser()
    nc.Created = time.ctime()  
    nc.Software = 'Created with Python ' + sys.version
    nc.NetCDF_Lib = str(netCDF4.getlibversion())
    
    # Create dimensions
    nc.createDimension('time',0)          # The unlimited dimension
    nc.createDimension('node',nodes)      # Number of nodes
    
    # Create time vector
    nc.createVariable('time','f8',('time'))
    nc.variables['time'].long_name = 'model time'
    nc.variables['time'].standard_name = 'time'
    nc.variables['time'].units = 'seconds since ' + ncdate
    nc.variables['time'].base_date = ncdate
    
    # Create variables
    nc.createVariable('wvx','f8',('time','node'))
    nc.variables['wvx'].long_name = 'applied horizontal wind velocity in the x direction'
    nc.variables['wvx'].units = 'meters of water'
    
    nc.createVariable('wvy','f8',('time','node'))
    nc.variables['wvy'].long_name = 'applied horizontal wind velocity in the y direction'
    nc.variables['wvy'].units = 'meters of water'
    
    nc.createVariable('prn','f8',('time','node'))
    nc.variables['prn'].long_name = 'applied atmospheric pressure at the free surface'
    nc.variables['prn'].units = 'meters of water'
    
    # Open ascii file and store data in nc file
    fobj = open(fort22,'r')
    tt = 0
    secs = 0
    line = fobj.readline()
    print('Working on time step:')
    while line:
        nc.variables['time'][tt] = secs
        print('    '+str(tt))
        for nn in range(0,nodes):
            tmpline = [np.float64(a) for a in line.split(',')]
            nc.variables['wvx'][tt,nn] = np.float64(tmpline[1])
            nc.variables['wvy'][tt,nn] = np.float64(tmpline[2])
            nc.variables['prn'][tt,nn] = np.float64(tmpline[3])
            
            line = fobj.readline()
        tt += 1
        secs += wtiminc
    fobj.close()
    
    # All done here
    nc.close()
    fobj.close()
    
