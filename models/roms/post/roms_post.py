"""

Tools to manage ROMS output.

Authors:
-------
Gabriel Garcia Medina
    Nearshore Modeling Group
    ggarcia@coas.oregonstate.edu
Saeed Moghimi


Dependencies:
-------------
numpy,netCDF4

Internal dependencies:
----------------------
none

"""

__author__ = "Gabriel Garcia Medina"
__email__ = "ggarcia@coas.oregonstate.edu"
__group__ = "Nearshore Modeling Group"


# Modules used
import numpy as np
import netCDF4


# Curl of the quantities will be computed at psi points. 
def curl(x,y,u,v):
    '''
    Compute curl of parameters in terms of ROMS variables. Output is given at 
    PSI points. 
    
    q = curl(x_v,y_u,u,v)
    '''
    return ( np.diff(v, axis=-1)/np.diff(x, axis=-1) -
             np.diff(u, axis=-2)/np.diff(y, axis=-2) )


# Compute divergence
def divergence(x,y,u,v):
    '''
    Compute the divergence of the parameters in terms of ROMS variables
    
    div = divergence(x_u,y_v,u,v)
    '''
    return ( (np.diff(u, axis=-1)/np.diff(x, axis=-1))[...,1:-1,:] + 
             (np.diff(v, axis=-2)/np.diff(y, axis=-2))[...,:,1:-1] )

def curlxy_old(x,y,z,u,v,w):
    """
    Compute horizontal vorticity components in terms of ROMS variables.
    
    qx,qy = curlxy_old(x_rho,y_rho,z(rho),u,v,w)
    
    Parameters:
    -----------
    x_rho, y_rho : x and y locations of rho points
    z            : layer elevation matrix (can be time varying) at rho points
    u,v,w        : flow component matrices on their corresponding grids
    
    Notes:
    ------
    qx = dw/dy - dv/dz
    qy = du/dz - dw/dx
    
    qx is computed at inner w points in the vertical and v poins in eta
    qy is computed at inner w points in the vertical and u points in xi
    
    See also:
    --------- 
    groms.curl
    octant.roms.nc_depths
    
    """
    
    # Compute water surface elevation at y points
    ztmp = 0.5*(z[...,0:-1,:] + z[...,1:,:])        
    qx = ((np.diff(w,axis=-2)/np.diff(y,axis=-2))[...,1:-1,:,:] - 
          np.diff(v,axis=-3)/np.diff(ztmp,axis=-3))
    # Compute water surface elevation at x pointsj
    ztmp = 0.5*(z[...,0:-1] + z[...,1:])
    qy = (-1.0*(np.diff(w,axis=-1)/np.diff(x,axis=-1))[...,1:-1,:,:] + 
          np.diff(u,axis=-3)/np.diff(ztmp,axis=-3))
    
    return qx,qy
    


def curlxy(x,y,z,u,v,w):
    """
    Compute horizontal vorticity components in terms of ROMS variables.
    
    qx,qy = curlxy(x_rho,y_rho,z(rho),u,v,w)
    
    Parameters:
    -----------
    x_rho, y_rho : x and y locations of rho points
    z            : layer elevation matrix (can be time varying) at rho points
    u,v,w        : flow component matrices on their corresponding grids
    
    Notes:
    ------
    qx = dw/dy - dv/dz
    qy = du/dz - dw/dx
    
    qx is computed at inner rho points in the vertical and v poins in eta
    qy is computed at inner rho points in the vertical and u points in xi
    
    similar to curlxy_old but computed using gradients and not differences.
    
    See also:
    --------- 
    groms.curl
    octant.roms.nc_depths
    
    
    """
    
    # Compute water surface elevation at y points
    ztmp = 0.5*(z[...,0:-1,:] + z[...,1:,:])        

    # qx
    dwdy = np.gradient(w)[-2]/np.gradient(y)[-2]
    dwdy = 0.5*(dwdy[...,1:,:,:] + dwdy[...,0:-1,:,:])
    dwdy = 0.5*(dwdy[...,1:,:] + dwdy[...,0:-1,:])
    dvdz = np.gradient(v)[-3]/np.gradient(ztmp)[-3]
    qx = dwdy - dvdz
    
    # Compute water surface elevation at x pointsj
    ztmp = 0.5*(z[...,0:-1] + z[...,1:])
    
    # qy
    dwdx = np.gradient(w)[-1]/np.gradient(x)[-1]
    dwdx = 0.5*(dwdx[...,1:,:,:] + dwdx[...,0:-1,:,:])
    dwdx = 0.5*(dwdx[...,1:] + dwdx[...,0:-1])
    dudz = np.gradient(u)[-3]/np.gradient(ztmp)[-3] 
    qy = dudz - dwdx
    
    return qx,qy
        
        

#===============================================================================
# Write SWAN bathymetry file from either ROMS history file or bathymetry input 
#===============================================================================
def roms_to_swan_bathy_rect(hisfile,outfld,sstacks):
  ''' 
  Generate a SWAN bathymetry file from either a ROMS history or bathymetry input
  file. 
  
  roms_to_swan_bathy(hisfile,outfld,sstacks)
  
  Parameters
  ----------
  hisfile  : ROMS history or bathymetry input netCDF file
  outfld   : Folder to save the output file
  sstacks  : Number of times to clone the bathymetry (must be an odd number)
  
  Returns
  -------
  A text file (swan_bathy.bot) that contains a regular bathymetry file for SWAN.
  
  Notes
  -----
  If the file is detected to be in curvilinear coordinates it will interpolate
  to rectangular coordinates and generate a grid with a spacing equal to the 
  minimum one found in the file. 
    
  '''  
  
  # Extend the grid
  if sstacks % 2 == 0:
    print "sstacks input cannot be an even number, quitting..."
    return
  
  # Load variables of interest from the ocean_his.nc file
  ncfile = netCDF4.Dataset(hisfile,'r')  
  h = ncfile.variables['h'][:]
  x_rho = ncfile.variables['x_rho'][:]
  y_rho = ncfile.variables['y_rho'][:]  
  ncfile.close()
  
  # Verify is the grid is curvilinear
  dif1 = np.diff(x_rho,axis=0)
  dif2 = np.diff(x_rho,axis=1)
  dif3 = np.diff(y_rho,axis=0)
  dif4 = np.diff(y_rho,axis=1)
  
  if (len(np.unique(dif1))!=1 or len(np.unique(dif2))!=1 or 
      len(np.unique(dif3))!=1 or len(np.unique(dif4))!=1):
      
      # Generate new x_rho axis with minimum grid spacing
      dx_int = np.sort(np.array((dif1.min(),dif2.min())))
      if dx_int[0] == 0:
          dx_int = dx_int[1]
      else:
          dx_int = dx_int[0]
      x_int = np.linspace(x_rho.min(),x_rho.max(),
                          (x_rho.max() - x_rho.min())/dx_int + 1)
      

      # Generate new y_rho axis with minimum grid spacing
      dy_int = np.sort(np.array((dif3.min(),dif4.min())))
      if dy_int[0] == 0:
          dy_int = dy_int[1]
      else:
          dy_int = dy_int[0]
      y_int = np.linspace(y_rho.min(),y_rho.max(),
                          (y_rho.max() - y_rho.min())/dy_int + 1)
      
      # Interpolate bathymetry to new rho points (simple bilinear interpolation)      
      ff = spi.RectBivariateSpline(y_rho[:,0],x_rho[0,:],h)
      h = ff(y_int,x_int)
      [x_rho,y_rho]=np.meshgrid(x_int,y_int)
      
  
  # Grid extension
  if sstacks > 1:
    x_rho = np.repeat(x_rho,sstacks,axis=0)    
    y_dx = y_rho[1,0] - y_rho[0,0]
    y_ll = ((sstacks - 1)/2) * (y_rho[0,0] - y_rho[-1,0]) - y_dx
    y_ur = ((sstacks - 1)/2) * (y_rho[-1,0] - y_rho[0,0]) + y_rho[-1,0] + y_dx
    y_tmp = np.linspace(y_ll,y_ur,(y_ur-y_ll)/y_dx+1)
    [cols,rows] = y_rho.shape
    y_rho = np.array([y_tmp,]*rows).transpose()
   
  # Print text file with extended and interpolated bathymetry
  fid = open(outfld+'/swan_bathy.bot', 'w')
  for cc in range(sstacks):
    for aa in range(h.shape[0]):
      for bb in range(h.shape[1]):
        fid.write('%10.3f' % h[aa,bb])
      fid.write('\n')
  fid.close()
  
  #---------------------------------------------------------- Output for swan.in
  print ' '
  print "======================================================================"
  print "Created swan_coord.grd and swan_bathy.bot"
  print ('CGRID REGULAR ' + np.str(x_rho.min()) + ' ' + 
         np.str(y_rho.min()) + ' 0 ' + np.str(x_rho.max() - x_rho.min()) + ' ' +
         np.str(y_rho.max() - y_rho.min()) + ' ' + 
         np.str(x_rho.shape[1]-1) + ' ' + np.str(x_rho.shape[0]-1) + 
         ' CIRCLE ...')
  print ('INPGRID BOTTTOM REGULAR ' + np.str(x_rho.min()) + ' ' +
           np.str(y_rho.min()) + ' 0 ' + np.str(x_rho.shape[1]-1) + ' ' +
           np.str(x_rho.shape[0]-1) + ' ' + 
           np.str(np.abs(x_rho[0,1] - x_rho[0,0])) + ' ' +
           np.str(np.abs(y_rho[1,0] - y_rho[0,0])) + ' EXC ...')
  print "======================================================================"    
  
  
  
  
#===============================================================================
# Write SWAN bathymetry file from either ROMS history file or bathymetry input 
#===============================================================================
def roms_to_swan_bathy_curv(hisfile,outfld):
  ''' 
  Generate a SWAN bathymetry file from either a ROMS history or bathymetry input
  file. 
  
  roms_to_swan_bathy_curv(hisfile,outfld)
  
  Parameters
  ----------
  hisfile  : ROMS history or bathymetry input netCDF file
  outfld   : Folder to save the output files
  
  Returns
  -------
  Two text files (swan_bathy.bot, swan_coord.grd) that contain the bathymetry 
  and coordinates of the grid for SWAN input. 
  
  Notes
  -----
    
    
  '''  
   
  # Load variables of interest from the ocean_his.nc file
  ncfile = netCDF4.Dataset(hisfile,'r')  
  h = ncfile.variables['h'][:]
  x_rho = ncfile.variables['x_rho'][:]
  y_rho = ncfile.variables['y_rho'][:]  
  ncfile.close()
  
   
  # Print text file with extended and interpolated bathymetry
  fid = open(outfld+'/swan_bathy.bot', 'w')  
  for aa in range(h.shape[0]):
      for bb in range(h.shape[1]):
          fid.write('%12.4f' % h[aa,bb])
      fid.write('\n')
  fid.close()
  
  # Print text file with extended and interpolated bathymetry
  fid = open(outfld+'/swan_coord.bot', 'w')  
  for aa in range(x_rho.shape[0]):
      for bb in range(x_rho.shape[1]):
          fid.write('%12.6f' % x_rho[aa,bb])
          fid.write('\n')
  for aa in range(y_rho.shape[0]):
      for bb in range(y_rho.shape[1]):
          fid.write('%12.6f' % y_rho[aa,bb])
          fid.write('\n')          
  fid.close()  
  
  #---------------------------------------------------------- Output for swan.in
  print ' '
  print "========================================================"
  print "Created swan_coord.grd and swan_bathy.bot"
  print ('CGRID CURVILINEAR ' + np.str(h.shape[1]-1) + ' ' + 
         np.str(h.shape[0]-1) + ' CIRCLE ...')
  print ('INPGRID BOTTTOM CURVILINEAR 0 0 ' + np.str(h.shape[1]-1) + ' ' + 
         np.str(h.shape[0]-1) + ' EXC ...')
  print "========================================================"  


                
