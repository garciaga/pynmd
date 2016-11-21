from __future__ import division,print_function,absolute_import

import sys

import netCDF4
import numpy as np
import time
import getpass
from collections import defaultdict


def gridVars():
    '''
    Returns grid variable information, list is incomplete
    See
    Data/ROMS/CDL/
    '''
    gridInfo = defaultdict(dict)
    
    gridInfo['x_rho']['long_name'] = 'x-locations of RHO-points'
    gridInfo['x_rho']['units'] = 'meter'
    gridInfo['x_rho']['field'] = 'x_rho, scalar'
    gridInfo['x_rho']['dimensions'] = 'r2dvar'
    
    gridInfo['y_rho']['long_name'] = 'y-locations of RHO-points'
    gridInfo['y_rho']['units'] = 'meter'
    gridInfo['y_rho']['field'] = 'x_rho, scalar'
    gridInfo['y_rho']['dimensions'] = 'r2dvar'

    gridInfo['x_u']['long_name'] = 'x-locations of U-points'
    gridInfo['x_u']['units'] = 'meter'
    gridInfo['x_u']['field'] = 'x_u, scalar'
    gridInfo['x_u']['dimensions'] = 'u2dvar'
    
    gridInfo['y_u']['long_name'] = 'y-locations of U-points'
    gridInfo['y_u']['units'] = 'meter'
    gridInfo['y_u']['field'] = 'x_u, scalar'
    gridInfo['y_u']['dimensions'] = 'u2dvar'
    
    gridInfo['x_v']['long_name'] = 'x-locations of V-points'
    gridInfo['x_v']['units'] = 'meter'
    gridInfo['x_v']['field'] = 'x_rho, scalar'
    gridInfo['x_v']['dimensions'] = 'v2dvar'
    
    gridInfo['y_v']['long_name'] = 'y-locations of V-points'
    gridInfo['y_v']['units'] = 'meter'
    gridInfo['y_v']['field'] = 'x_rho, scalar'
    gridInfo['y_v']['dimensions'] = 'v2dvar'
    
    gridInfo['x_psi']['long_name'] = 'x-locations of PSI-points'
    gridInfo['x_psi']['units'] = 'meter'
    gridInfo['x_psi']['field'] = 'x_rho, scalar'
    gridInfo['x_psi']['dimensions'] = 'p2dvar'
    
    gridInfo['y_psi']['long_name'] = 'y-locations of PSI-points'
    gridInfo['y_psi']['units'] = 'meter'
    gridInfo['y_psi']['field'] = 'x_rho, scalar'
    gridInfo['y_psi']['dimensions'] = 'p2dvar'
    
    gridInfo['h']['long_name'] = 'bathymetry at RHO-points'
    gridInfo['h']['units'] = 'meter'    
    gridInfo['h']['field'] = 'bath, scalar'
    gridInfo['h']['dimensions'] = 'r2dvar'
    
    gridInfo['f']['long_name'] = 'Coriolis parameter at RHO-points'
    gridInfo['f']['units'] = 'second-1'
    gridInfo['f']['field'] = 'coriolis, scalar'
    gridInfo['f']['dimensions'] = 'r2dvar'
    
    return gridInfo


def gridDims():
    '''
    Returns grid dimensions, incomplete list many nulvars to add
    '''
    g = defaultdict(dict)
    
    # 2D Vars
    g['r2dvar'] = ('eta_rho','xi_rho')
    g['p2dvar'] = ('eta_psi','xi_psi')
    g['u2dvar'] = ('eta_u','xi_u')
    g['v2dvar'] = ('eta_v','xi_v')
    
    # 3D Vars
    g['r3dvar'] = ('s_rho','eta_rho','xi_rho')
    g['p3dvar'] = ('s_rho','eta_psi','xi_psi')    
    g['u3dvar'] = ('s_rho','eta_u','xi_u')
    g['v3dvar'] = ('s_rho','eta_v','xi_v')
    g['w3dvar'] = ('s_w','eta_rho','xi_rho')
    
    # Nulvars
    g['nulvar']['zeta_west']  = ('eta_rho',)
    g['nulvar']['zeta_east']  = ('eta_rho',)
    g['nulvar']['zeta_north'] = ('xi_rho',)
    g['nulvar']['zeta_south'] = ('xi_rho',)

    g['nulvar']['ubar_west']  = ('eta_u',)
    g['nulvar']['ubar_east']  = ('eta_u',)
    g['nulvar']['ubar_north'] = ('xi_u',)
    g['nulvar']['ubar_south'] = ('xi_u',)
    
    g['nulvar']['vbar_west']  = ('eta_v',)
    g['nulvar']['vbar_east']  = ('eta_v',)
    g['nulvar']['vbar_north'] = ('xi_v',)
    g['nulvar']['vbar_south'] = ('xi_v',)
       
    g['nulvar']['u_west']  = ('s_rho','eta_u')
    g['nulvar']['u_east']  = ('s_rho','eta_u')
    g['nulvar']['u_north'] = ('s_rho','xi_u')
    g['nulvar']['u_south'] = ('s_rho','xi_u')
    
    g['nulvar']['v_west']  = ('s_rho','eta_v')
    g['nulvar']['v_east']  = ('s_rho','eta_v')
    g['nulvar']['v_north'] = ('s_rho','xi_v')
    g['nulvar']['v_south'] = ('s_rho','xi_v')

    g['nulvar']['temp_west']  = ('s_rho','eta_rho')
    g['nulvar']['temp_east']  = ('s_rho','eta_rho')
    g['nulvar']['temp_north'] = ('s_rho','xi_rho')
    g['nulvar']['temp_south'] = ('s_rho','xi_rho')

    g['nulvar']['salt_west']  = ('s_rho','eta_rho')
    g['nulvar']['salt_east']  = ('s_rho','eta_rho')
    g['nulvar']['salt_north'] = ('s_rho','xi_rho')
    g['nulvar']['salt_south'] = ('s_rho','xi_rho')
       
    # Done    
    return g    
    
    
def cfl(maxdep,dt_3d,num_dt_2d,dxx,dyy,S3D=True):
    """
    -------------------------------------------------------------------------
     Courant number calculation for ROMS:
     Cu = c * dt * SQRT (1/dx^2 + 1/dy^2)
     
     where c=SQRT(g*h) is phase speed for barotropic mode, and dx, dy
     are grid spacing in each direction.
    
    -------------------------------------------------------------------------
     Written by: Cigdem Akan (cakan@coas.oregonstate.edu)
     2012, Last Revision: 19-March-2013
    
    -------------------------------------------------------------------------#
    Call like:
    cfl(maxdep,dt_3d,num_dt_2d,dxx,dyy,S3D=True)
    
    params:
    maxdep      maximum depth in the domain (m)
    dt_3d       if S3D=True => baroclinic time step , else => barotropic time step if 2D
    num_dt_2d   number of 2d in between
    dxx         minimum dx in meter
    dyy         minimum dy in meter
    S3D=True    baroclinic (3D) switch
    example:
    cfl(10,2,4,10,10,S3D=True)

    """
    
    #S3D = 1 # 1 for 3D simulation, 0 for 2D simulation   

    g   = 9.81       # m/s2
    h   = maxdep * 1.0        # maximum depth in the domain (m)
    dt1 = dt_3d  * 1.0        # baroclinic time step if 3D, barotropic time step if 2D
    dt2 = num_dt_2d  * 1.0

    if S3D:
        print(' Evaluate a simple CFL check for baroclinic run')
        dt  = 1.0*dt1/dt2  
    else:
        print(' Evaluate a simple CFL check for 2D run')
        dt = dt1

    dx  = dxx * 1.0 # meters
    dy  = dyy * 1.0 # meters

    c  = np.sqrt(g*h) # m/s
    Cu = c*dt*np.sqrt((1./dx**2)+(1./dy**2))

    # I am not sure what gamma is exactly, I could not find a very good 
    # explanation for it but, apparently it is better if it is greater than 0.
    gamma = 0.284*(1.-10./dt1)

    print('    > The Courant number = ',Cu)
    print('    > Gamma              = ',gamma)
    print('    > nHis   (30 min)    = ',30*60/ dt1)
    print('    > nHis   (15 min)    = ',15*60/ dt1)
    print('    > ninfo  ( 1 min)    = ',1*60 / dt1)
    print('    > ndefHis(daily)     = ',24*60*60/dt1)


def writeROMSGrid(outFile,variables,varinfo,
                  timeinfo=None,verbose=False):
    """
    Code to write roms grids based on input variables
    
    INPUT:
    ------
    outFile    : Full path to output netcdf file
    variables  : Dictionary containing variables. (see notes)
    varinfo    : Path to ROMS varinfo.dat. Usually under:
                 /path/to/roms/trunk/ROMS/External/varinfo.dat
    timeinfo   : Time vector in nested dictionary form (see example)
    verbose    : some info displayed on the terminal                 
    
    RETURNS:
    --------
    
    NOTES:
    ------
    -  variables need rho variables (x_rho,y_rho)
    
    TODO:
    -----
    -  Include a varinfo.dat with the pynmd package
    
    EXAMPLE:
    ---------
    >>> import numpy as np
    >>> from collections import defaultdict
    >>> import pynmd.models.roms.pre as groms    
    >>>
    >>> x_rho,y_rho = np.meshgrid(np.arange(0,10,1),np.arange(0,10,1)
    >>> variables = {}
    >>> variables['x_rho'] = x_rho
    >>> variables['y_rho'] = y_rho
    >>>
    >>> variables['type'] = 'ROMS FORCING FILE'
    >>>
    >>> timeinfo = defaultdict(time)
    >>> timeinfo['sms_time']['units'] = 'seconds since 1970-01-01 00:00:00 UTC'
    >>> timeinfo['sms_time']['data'] = np.arange(0,1000,10)
    >>> timeinfo['sms_time']['long_name'] = 'surface momentum stress time'
    >>>
    >>> groms.writeROMSGrid(outFile,...)
    
    """

    # Testing only
    #outFile = ('/home/shusin3/users/ggarcia/projects/other/' + 
    #           'columbiaRiverPlume/11-romsIdealized/03-extendedNR2/' + 
    #           'orIdealizedExtStr.nc')
    #varinfo = ('/home/shusin3/users/ggarcia/projects/other/' + 
    #           'columbiaRiverPlume/94-roms/trunk/ROMS/External/varinfo.dat')
    #verbose = True
    
   
    # Read varinfo -------------------------------------------------------------
    if verbose:
        print('Reading Varinfo')
        print('  ' + varinfo)
    
    # Create dictionary to get information from
    varinfodict = defaultdict(dict)
    
    # Structure of the variables 
    # hardcoded because I couldn't find a generalized description
    varStruct = ['long_name','units','field','time','id','dimensions']
       
    # Read the varinfo
    with open(varinfo) as f:
        
        # Read into a list of lines
        lines = f.readlines()
        
        # Loop through lines
        for line in lines:
            
            #Skip lines starting with !
            if line.startswith("!"):
                continue
            
            # Skip blank lines
            if line.isspace():
                continue
            
            # Skip other non variable lines
            if line.startswith("'$"):
                continue
            
            # Read variables and allocate properties
            if line.startswith("'"):
                tmpVar = line.split(' ')[0][1:-1]
                cnt = -1
                continue
            
            # Add counter variable and allocate only the necessary fields
            cnt += 1            
            if cnt > (len(varStruct) - 1):
                continue
            
            varinfodict[tmpVar][varStruct[cnt]] =\
              line.split('  ')[1].split('\n')[0][1:-1]
            
    f.close()    
    
    # Import grid variables and update the dictionary
    gridVariables = gridVars()
    varinfodict.update(gridVariables)
    
    # Create netCDF file -------------------------------------------------------
    if verbose:
        print('Creating: ' + outFile)
        
    nc = netCDF4.Dataset(outFile, 'w', format='NETCDF4')    
    nc.Author = getpass.getuser()
    nc.Created = time.ctime()
    #nc.Owner = 'Nearshore Modeling Group (http://ozkan.oce.orst.edu/nmg)'
    nc.Software = 'Created with Python ' + sys.version
    nc.NetCDF_Lib = str(netCDF4.getlibversion())    
    
    if 'type' in variables.keys():
        nc.type = variables['type']
    
    # Create dimensions --------------------------------------------------------
    
    if verbose:
        print("  Creating Dimensions ...")
    # Horizontal dimensions
    nc.createDimension('xi_psi',np.size(variables['x_rho'],1)-1)
    nc.createDimension('xi_rho',np.size(variables['x_rho'],1))
    nc.createDimension('xi_u',np.size(variables['x_rho'],1)-1)
    nc.createDimension('xi_v',np.size(variables['x_rho'],1))
    nc.createDimension('eta_psi',np.size(variables['x_rho'],0)-1)
    nc.createDimension('eta_rho',np.size(variables['x_rho'],0))
    nc.createDimension('eta_u',np.size(variables['x_rho'],0))
    nc.createDimension('eta_v',np.size(variables['x_rho'],0)-1)
    
    # Vertical dimensions
    if 's_rho' in variables.keys():
        nc.createDimension('s_rho',variables['s_rho'].size)
        nc.createDimension('s_w',variables['s_w'].size)
    
    # Time dimension and variable (unlimited,netCDF4 supports multiple)
    if timeinfo:
        
        if verbose:
            print('  Writing Time Variables:')
            
        for aa in timeinfo.keys():
            
            if verbose:
                print('    ' + aa)
                
            # Create dimension
            nc.createDimension(aa,0)
            nc.createVariable(aa,'f8',(aa))
            
            # Write variable information
            for bb in timeinfo[aa].keys():
                
                if bb == 'data':
                    nc.variables[aa][:] = timeinfo[aa]['data'] 
                else:
                    nc.variables[aa].__setattr__(bb,timeinfo[aa][bb])
        
    # Write other variables
    if verbose:
        print('  Writing variables:')

    # Loading dimensions
    dimInfo = gridDims()

    for aa in variables.keys():
        
        # Determine dimensions
        try:
            tmpdims = dimInfo[varinfodict[aa]['dimensions']]
            if varinfodict[aa]['dimensions'] == 'nulvar':
                tmpdims = tmpdims[aa]
        except:
            continue
        
        if verbose:
            print('    ' + aa)
            
        if 'time' in varinfodict[aa].keys():
            tmpdims = (varinfodict[aa]['time'],) + tmpdims
        
        # Create Variable
        nc.createVariable(aa,'f8',tmpdims)
        
        # Write attributes
        for bb in varinfodict[aa].keys():
            if bb == 'dimensions' or bb == 'id':
                continue
            nc.variables[aa].__setattr__(bb,varinfodict[aa][bb])
        
        # Write data
        nc.variables[aa][:] = variables[aa]
                
        
    # Close netcdf file
    nc.close()
    
    if verbose:
        print('Done')
        
    