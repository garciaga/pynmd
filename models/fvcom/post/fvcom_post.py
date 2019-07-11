"""
Collection of functions to postprocess fvcom data
"""

from __future__ import division,print_function

import numpy as _np
import netCDF4 as _netCDF4

# ==============================================================================
# Create restart files for temperature and salinity
# ==============================================================================
def restart_its(resFile,ncFile,outFld):
    """
    Create the initial condition file for temperature and salinity for a 
    restart. The time step is taken from resFile and the temperature and 
    salinity are taken from a netcdf output file. 

    PARAMETERS:
    -----------
    resFile: Full path to restart file. The filename should be in the 
             Casename_restart.dat format. 
    ncFile : netcdf file with model output. Salinity and Temperature fields
             must be present in the file. 
    outFld : Path to folder where Casename_its.dat will be written

    OUTPUT:
    -------
    Casename_its.dat : File will be created from the closest timestep 
                       available in ncFile.

    NOTES:
    ------
    1. To use the initial temperature and salinity file for the restart
         RESTART = Hot_cold_s

    """

    # Read the time step from the restart file
    fid = open(resFile,'r')
    tInd = _np.int(fid.readline().rstrip())
    fid.close()

    # Read temperature and salinity from the netCDF file    
    nc = _netCDF4.Dataset(ncFile,'r')
    # Find the time to load
    ot = _np.array(nc.variables['iint'][:])
    tLoad = _np.argmin(_np.abs(ot-tInd))
    # Warning if the restart time is not the same as output time
    dt = _np.abs(ot[tLoad] - tInd)
    if dt >=1:
        print('Restart time and model output differ')
        print('  Restart Time Step: ' + _np.str(tInd))
        print('   Output Time Step: ' + _np.str(ot[tLoad]))

    # Load data
    temp = _np.array(nc.variables['temp'][tLoad,...])
    salt = _np.array(nc.variables['salinity'][tLoad,...])    

    # Get the run ID 
    runId = ncFile.split('/')[-1].split('_')[0]

    # Write the Casename_its.dat file
    fid = open(outFld + '/' + runId + '_its.dat','w')
    fid.write('TEMPERATURE AND SALINITY FOR STEP: ' + _np.str(tInd) + '\n')
    fid.write('observed\n')

    # Loop over node points
    for aa in range(temp.shape[1]):
        # Write temperature for each layer
        for bb in range(temp.shape[0]):
            fid.write('{:8.2f}'.format(temp[bb,aa]))
        fid.write('\n')
        # Write Salinity for each layer
        for bb in range(temp.shape[0]):
            fid.write('{:8.2f}'.format(salt[bb,aa]))
        fid.write('\n')

    # Close the file
    fid.close()

def computeVorticity(u,v,a1u,a2u,ele):
    """
    Compute vorticity
    
    Parameters:
    -----------
    u:   Eastward directed velocity (nele,)
    v:   Northward directed velocity (nele,)
    a1u, a2u: Shape parameters from FVCOM
    ele: Elements

    Returns:
    --------
    q = dvdx - dudy

    Notes:
    ------
    Right now this is meant only for a snapshot of u and v. It should be trivial
    to extend to time and depth dependent arrays. 

    Under development! Wrong answers given at the moment
    """

    # Get neighboring elements
    nbe = np.zeros_like(ele)


    # Preallocate vertical vorticity variable
    dvdx = (a1u[0,:] * v + 
            a1u[1,:] * v[ele[:,0]] + 
            a1u[2,:] * v[ele[:,1]] + 
            a1u[3,:] * v[ele[:,2]])
    dudy = (a2u[0,:] * u + 
            a2u[1,:] * u[ele[:,0]] + 
            a2u[2,:] * u[ele[:,1]] + 
            a2u[3,:] * u[ele[:,2]])
    return dvdx - dudy
