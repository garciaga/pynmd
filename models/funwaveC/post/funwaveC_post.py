"""
Collection of codes to process funwaveC output data
"""

from __future__ import division,print_function

import numpy as np
import netCDF4

# ==================================================================
# Create NetCDF file
# ==================================================================    
def convert_output(workfld,outfile,time_int,bathyfile=None,inpfile=None):
    '''
    
    Parameters:
    -----------
    workfld      : Path to the folder where output files reside
    outfile      : Output NetCDF file
    time_int     : Time interval between output files [s] (will be updated if
                   input file is provided)
    bathyfile    : Full path to input netcdf bathy file (optional)
    inpfile      : Funwave input file used for metadata (optional)
         
    Output:
    -------
    NetCDF File with the variables in the folder. Not all are supported so you 
    may need to edit this file.
    
    '''
    
    # For testing only
    #workfld = '/scratch/temp/ggarcia/spectralWidthFunwave/02-runs/results/'
    #outfile = '/scratch/temp/ggarcia/spectralWidthFunwave/02-runs/test.nc'
    #time_int = 0.1
    #bathyfile = '/scratch/temp/ggarcia/spectralWidthFunwave/02-runs/depth.nc'
    #inpfile = '/scratch/temp/ggarcia/spectralWidthFunwave/02-runs/input.txt'

    # Get variable information ------------------------------------------------
    archivos = os.listdir(workfld)                      # Get all files
    tmpvars = [x.split('_')[0] for x in archivos]       # All variables
    tmpvars = list(set(tmpvars))                        # Unique variables
    
    # If no variables found exit    
    if not tmpvars:
        print('No files found in ' + workfld)
        print('Quitting ...')
        return None
    
    # Make sure variables are within the supported ones
    supported_vars_time = ['eta','etamean','havg','hmax','hmin','hrms',
                           'mask','mask9','MFmax','u','umax','umean','v',
                           'vmean','VORmax']
    vars_2d = [x for x in tmpvars if x in supported_vars_time]
    
    
    # Read input file if provided ----------------------------------------------
    if inpfile:
        
        print("Reading data from input file:")
        print("  " + inpfile)
        
        # Open file
        tmpinpfile = open(inpfile,'r')
        
        # Output dictionary
        inpinfo = {}
        
        # Extract information (need to add wavemaker support)
        for tmpline in tmpinpfile:
            
            # Skip blank lines
            if len(tmpline.strip()) == 0:
                continue
            
            if "TITLE" == tmpline.split()[0]:
                inpinfo['title'] = tmpline.split()[2]
            elif "PX" == tmpline.split()[0]:
                inpinfo['px'] = tmpline.split()[2]
            elif "PY" == tmpline.split()[0]:
                inpinfo["py"] = tmpline.split()[2]
            elif "Mglob" == tmpline.split()[0]:
                inpinfo['mglob'] = tmpline.split()[2]
            elif "Nglob" == tmpline.split()[0]:
                inpinfo['nglob'] = tmpline.split()[2]
            elif "TOTAL_TIME" == tmpline.split()[0]:
                inpinfo['total_time'] = tmpline.split()[2]
            elif "PLOT_INTV" == tmpline.split()[0]:
                inpinfo["plot_intv"] = tmpline.split()[2]
                time_int = np.double(tmpline.split()[2])
            elif "DX" == tmpline.split()[0]:
                inpinfo['dx'] = tmpline.split()[2]
                dx = float(tmpline.split()[2])
            elif "DY" == tmpline.split()[0]:
                inpinfo['dy'] = tmpline.split()[2]
                dy = float(tmpline.split()[2])
            elif "WAVEMAKER" == tmpline.split()[0]:
                inpinfo['wavemaker'] = tmpline.split()[2]
            elif "PERIODIC" == tmpline.split()[0]:
                inpinfo['periodic'] = tmpline.split()[2]
            elif "SPONGE_ON" == tmpline.split()[0]:
                sponge = tmpline.split()[2]
                inpinfo['sponge'] = sponge
            elif "Sponge_west_width" == tmpline.split()[0] and sponge == 'T':
                inpinfo['sponge_west_width'] = tmpline.split()[2]
            elif "Sponge_east_width" == tmpline.split()[0] and sponge == 'T':
                inpinfo['sponge_east_width'] = tmpline.split()[2]
            elif "Sponge_wouth_width" == tmpline.split()[0] and sponge == 'T':
                inpinfo['sponge_south_width'] = tmpline.split()[2]
            elif "Sponge_worth_width" == tmpline.split()[0] and sponge == 'T':
                inpinfo['sponge_north_width'] = tmpline.split()[2]
            elif "R_sponge" == tmpline.split()[0] and sponge == 'T':
                inpinfo['r_sponge'] = tmpline.split()[2]
            elif "A_sponge" == tmpline.split()[0] and sponge == 'T':
                inpinfo['a_sponge'] = tmpline.split()[2]
            elif "DISPERSION" == tmpline.split()[0]:
                inpinfo['dispersion'] = tmpline.split()[2]
            elif "Gamma1" == tmpline.split()[0]:
                inpinfo['gamma1'] = tmpline.split()[2]
            elif "Gamma2" == tmpline.split()[0]:
                inpinfo['gamma2'] = tmpline.split()[2]
            elif "Gamma3" == tmpline.split()[0]:
                inpinfo['gamma3'] = tmpline.split()[2]
            elif "Beta_ref" == tmpline.split()[0]:
                inpinfo['beta_ref'] = tmpline.split()[2]
            elif "SWE_ETA_DEP" == tmpline.split()[0]:
                inpinfo['swe_eta_dep'] = tmpline.split()[2]
            elif "Friction_Matrix" == tmpline.split()[0]:
                inpinfo['friction_matrix'] = tmpline.split()[2]
            elif "Cd_file" == tmpline.split()[0]:
                inpinfo['cd_file'] = tmpline.split()[2]
            elif "Cd" == tmpline.split()[0]:
                inpinfo['cd'] = tmpline.split()[2]                              
            elif "Time_Scheme" == tmpline.split()[0]:
                inpinfo['time_scheme'] = tmpline.split()[2]
            elif "HIGH_ORDER" == tmpline.split()[0]:
                inpinfo['spatial_scheme'] = tmpline.split()[2]
            elif "CONSTRUCTION" == tmpline.split()[0]:
                inpinfo['construction'] = tmpline.split()[2]
            elif "CFL" == tmpline.split()[0]:
                inpinfo['cfl'] = tmpline.split()[2]                
            elif "MinDepth" == tmpline.split()[0]:
                inpinfo['min_depth'] = tmpline.split()[2]
            elif "MinDepthFrc" == tmpline.split()[0]:
                inpinfo['mindepthfrc'] = tmpline.split()[2]
            
        # Close file
        tmpinpfile.close()
            
    else:
        
        # Assume grid spacing to be 1 meter
        print("Input file not provided")
        dx = 1
        dy = 1
        inpinfo = False
    
    
    # Coordinates and depth ---------------------------------------------------
    # Bathymetry file provided
    if bathyfile: 
        
        print("Reading coordinates and depth from:")
        print("  " + bathyfile)
        
        ncfile = netCDF4.Dataset(bathyfile,'r')
        x_rho = ncfile.variables['x_rho'][:]
        if ncfile.variables.has_key('y_rho'):
            y_rho = ncfile.variables['y_rho'][:]
        else:
            y_rho = None
        h = ncfile.variables['h'][:]
        ncfile.close()
        
    # Bathymetry file not provided but have output file        
    elif os.path.isfile(workfld + '/dep.out'):
               
        # Fix this               
        h = np.loadtxt(workfld + '/dep.out')
        
        hdims = h.ndim
        if hdims == 1:
            x_rho = np.arange(0,h.shape[0],dx)
        elif hdims == 2:
            x_rho, y_rho = np.meshgrid(np.arange(0,h.shape[1],dx),
                                       np.arange(0,h.shape[0],dy))           
        else:
            print('Something is wrong with the depth file')
            return None    
    
    # No bathymetry file provided (I am not capable of reading the input file)
    else:
        print("No bathymetry file provided")
        print("You could copy your input bathymetry text file to ")
        print(workfld + '/dep.out')
        return None
        
    
    # Get dimensions of variables
    hdims = h.ndim
    
    # Create NetCDF file -------------------------------------------------------
    
    print("Creating " + outfile)
    
    # Global attributes  
    nc = netCDF4.Dataset(outfile, 'w', format='NETCDF4')
    nc.Description = 'Funwave Output'
    nc.Author = 'ggarcia@coas.oregonstate.edu'
    nc.Created = time.ctime()
    nc.Type = 'Funwave v2.1 snapshot output'
    nc.Owner = 'Nearshore Modeling Group'
    nc.Software = 'Created with Python ' + sys.version
    nc.NetCDF_Lib = str(netCDF4.getlibversion())
    nc.Source = workfld
    nc.Script = os.path.realpath(__file__)
    
    # Add more global variables to output
    if inpinfo:
        for tmpatt in inpinfo.keys():
            nc.__setattr__(tmpatt,inpinfo[tmpatt][:])
    
    # Create dimensions
    if hdims == 2:
        eta_rho, xi_rho = h.shape
        nc.createDimension('eta_rho', eta_rho)
    else:
        xi_rho = h.shape[0]
        eta_rho = 1
    nc.createDimension('xi_rho', xi_rho)            
    nc.createDimension('ocean_time',0)
    
    # Write coordinate axes ----------------------------------------------
    if hdims == 2:
        
        nc.createVariable('x_rho','f8',('eta_rho','xi_rho'))
        nc.variables['x_rho'].units = 'meter'
        nc.variables['x_rho'].longname = 'x-locations of RHO points'
        nc.variables['x_rho'][:] = x_rho
        
        nc.createVariable('y_rho','f8',('eta_rho','xi_rho'))
        nc.variables['y_rho'].units = 'meter'
        nc.variables['y_rho'].longname = 'y-locations of RHO points'
        nc.variables['y_rho'][:] = y_rho

        nc.createVariable('h','f8',('eta_rho','xi_rho'))
        nc.variables['h'].units = 'meter'
        nc.variables['h'].longname = 'bathymetry at RHO points'
        nc.variables['h'][:] = h


    else:
        
        nc.createVariable('x_rho','f8',('xi_rho'))
        nc.variables['x_rho'].units = 'meter'
        nc.variables['x_rho'].longname = 'x-locations of RHO points'
        nc.variables['x_rho'][:] = x_rho
        
        nc.createVariable('h','f8',('xi_rho'))
        nc.variables['h'].units = 'meter'
        nc.variables['h'].longname = 'bathymetry at RHO points'
        nc.variables['h'][:] = h
        
           
    # Create time vector -------------------------------------------
    tmpruns = [x for x in archivos if x.split('_')[0] == vars_2d[0]]
    tmpruns.sort()
    time_max = float(tmpruns[-1].split('_')[-1])
    time_min = float(tmpruns[0].split('_')[-1])
    
    
    twave = np.arange(time_min-1,time_int*(time_max-time_min+1)+time_min-1,
                      time_int)
    if twave.shape[0] > len(tmpruns):
        twave = twave[:-1]
        
    nc.createVariable('ocean_time','f8','ocean_time')
    nc.variables['ocean_time'].units = 'seconds since 2000-01-01 00:00:00'
    nc.variables['ocean_time'].calendar = 'julian'
    nc.variables['ocean_time'].long_name = 'beach time'
    nc.variables['ocean_time'][:] = twave
        
    # Create variables --------------------------------------------------------

    # Variable information
    varinfo = defaultdict(dict)

    varinfo['eta']['units'] = 'meter'
    varinfo['eta']['longname'] = 'water surface elevation'
    varinfo['etamean']['units'] = 'meter'
    varinfo['etamean']['longname'] = 'Mean wave induced setup'
    
    varinfo['u']['units'] = 'meter second-1'
    varinfo['u']['longname'] = 'Flow velocity in the xi direction'
    varinfo['v']['units'] = 'meter second-1'
    varinfo['v']['longname'] = 'Flow velocity in the eta direction'
    
    varinfo['umean']['units'] = 'meter second-1'
    varinfo['umean']['longname'] = 'Time-averaged flow velocity in xi direction'
    varinfo['vmean']['units'] = 'meter second-1'
    varinfo['vmean']['longname'] = 'Time-averaged flow velocity in eta direction'
    
    varinfo['umax']['units'] = 'meter second-1'
    varinfo['umax']['longname'] = 'Maximum flow velocity in xi direction'
    varinfo['vmax']['units'] = 'meter second-1'
    varinfo['vmax']['longname'] = 'Maximum flow velocity in eta direction'    
    
    varinfo['hmax']['units'] = 'meter'
    varinfo['hmax']['longname'] = 'Maximum wave height'
    varinfo['hmin']['units'] = 'meter'
    varinfo['hmin']['longname'] = 'Minimum wave height'
    varinfo['havg']['units'] = 'meter'
    varinfo['havg']['longname'] = 'Average wave height'
    varinfo['hrms']['units'] = 'meter'
    varinfo['hrms']['longname'] = 'Root mean squared wave height'
    
    varinfo['mask']['units'] = 'Boolean'
    varinfo['mask']['longname'] = 'Logical parameter for output wetting-drying'
    varinfo['mask9']['units'] = 'Boolean'
    varinfo['mask9']['longname'] = 'Logical parameter for output MASK9'
    
    varinfo['VORmax']['units'] = 'second-1'
    varinfo['VORmax']['longname'] = 'Maximum vorticity'
    varinfo['MFmax']['units'] = 'meter second-s'
    varinfo['MFmax']['longname'] = 'Maximum momentum flux'    
    
    
    # Create variables        
    if hdims == 1:
        nc_dims = ('ocean_time','xi_rho')
    else:
        nc_dims = ('ocean_time','eta_rho','xi_rho')
          
          
    print("Creating variables")          
    for aa in vars_2d:
        
        print('  ' + aa)
        
        # Create variable
        create_nc_var(nc,aa,nc_dims,varinfo[aa]['units'],
                      varinfo[aa]['longname'])
        
        try:
            tmpvar = np.loadtxt(workfld + '/' + aa + '_' + '%05.0f' % time_min)
        except ValueError:
            tmpvar = np.zeros_like(h) * np.NAN
                    
        nc.variables[aa][:] = np.expand_dims(tmpvar,axis=0)
        
        for bb in range(len(twave)):
            
            try:
                tmpvar = np.loadtxt(workfld + '/' + aa + '_' + 
                                    '%05.0f' % (bb + time_min))
            except ValueError:
                tmpvar = np.zeros_like(h) * np.NAN
                
            append_nc_var(nc,tmpvar,aa,bb-1)   
                    
    # Close NetCDF file
    print('Closing ' + outfile)
    nc.close()
    
    # End of function
