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

import os#,glob
import sys,time
#sys.path.append(os.getcwd())
import getpass
import numpy as np
import netCDF4

# Custom paths
import pynmd.models.tools.unstructured as ustr

# ==============================================================================
# Read fort.14 ASCII file and save as nc file
# ==============================================================================
def fort14_to_nc(fort14,**kwargs):
    """ 
    Script to read an unstructured grid in fort.14 format and store it in a 
    netcdf4 file

    PARAMETERS:
    -----------
    fort14: Path to fort14 file
    savename(optional)

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
    if 'savename' in kwargs:
        ncfile = kwargs['savename']
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
    
    grid = ustr.read_fort14(os.path.dirname(fort15)+'/fort.14')    
    neta = grid['neta']
    
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
    dic_el = {'nstae':nstae,'xel':xel,'yel':yel,'nameel':nameel}
    
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
    dic_ev = {'nstav':nstav,'xev':xev,'yev':yev,'nameev':nameev}
    
    # Concentration recording stations (fort.91 output)
    dic_ec ={} #dummy
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
        dic_ec = {'nstac':nstac,'xev':xec,'yev':yec,'nameec':nameec}
        
    # Meterology recording stations (units 71&72 output)
    dic_em ={} #dummy
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
        dic_em = {'nstam':nstam,'xem':xem,'yem':yem,'nameem':nameem}
    
    # Close ASCII file
    fobj.close()
    
    # Generate output
    return {**dic_el,**dic_ev,**dic_ec,**dic_em}