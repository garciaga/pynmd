__author__ = "Saeed Moghimi"
__copyright__ = "Copyright 2015, Portland State University"
__license__ = "GPL"
__version__ = "1.0"
__email__ = "moghimis@gmail.com"

import numpy as np
import os,sys

def write_init(fname,data):
    """
    Write Delft3D init txt file
    - data    : a dict contan required data (elev,u,v,temp,salt)
        u   (nz,ny,nx)
        v   (nz,ny,nx)
        temp(nz,ny,nx)
        salt(nz,ny,nx)

    - filename  
    
    
    #We need all this fields in one text file
    #1) water level
    #2) u-velocity (for each (k)-layer)
    #3) v-velocity (for each (k)-layer)
    #4) salinity (for each (k)-layer, if selected)
    #5) temperature (for each (k)-layer)
    #6) constituent 1(for each (k)-layer)
    #7) constituent 2 (for each (k)-layer)
    #8) constituent n (for each (k)-layer)
    """
    
    os.system('rm '+ fname)
    # read main inputs
    elev = data['elev'][:]
    uu   = data['uu']  [:]
    vv   = data['vv']  [:]
    salt = data['salt'][:]
    temp = data['temp'][:]
    
    if 'sed' in data.keys():
       sed  = data['sed'][:]
    
    #
    nlayer =   temp.shape[0]       
    ## Write
    fmt = '%.4g'
    f   = open(fname,'a')
    
    #1)water level
    np.savetxt(f,elev.T , fmt=fmt)
    
    #2) u-velocity (for each (k)-layer)
    for il in range (nlayer):
        np.savetxt(f,uu[il].T, fmt=fmt)
    
    #3) v-velocity (for each (k)-layer)
    for il in range (nlayer):
        np.savetxt(f,vv[il].T, fmt=fmt)
    
    #4) salinity (for each (k)-layer, if selected)
    for il in range (nlayer):
        np.savetxt(f,salt[il].T, fmt=fmt)
    
    #5) Temp (for each (k)-layer, if selected)
    for il in range (nlayer):
        np.savetxt(f,temp[il].T, fmt=fmt)
   
    if 'sed' in data.keys():
    #6) one sed component  (for each (k)-layer, if selected)
        if len(sed.shape) > 3 : 
            if len(sed.shape) > 3 : 
                for ised in range (len(sed)): #sys.exit('ERR : only one sed comp is implemented')
                    for il in range (nlayer):
                        np.savetxt(f,sed[ised,il].T, fmt=fmt)
        else:        
            for il in range (nlayer):
                np.savetxt(f,sed[il].T, fmt=fmt)
    
    f.close()
    
    

def write_rgh(fname,data):
    """
    Write Delft3D roughness txt file
    - data    : a dict contan required data (elev,u,v,temp,salt)
        rfgx (ny,nx)
        rfgy

    - filename  
    
    """
    os.system('rm '+ fname)
    # read main inputs
    rghx = data['rghx'][:]
    rghy = data['rghy'][:]
    #
    ## Write
    fmt = '%.4g'
    f   = open(fname,'a')
    
    #1)roughness for ux and uy
    np.savetxt(f,rghx.T , fmt=fmt)
    np.savetxt(f,rghy.T , fmt=fmt)
   
   
    f.close()
    
def write_dep(fname,data):
    """
    Write Delft3D dep txt file
    - data    : a dict contan  dep (ny,nx)

    - filename  
    
    """
    os.system('rm '+ fname)
    # read main inputs
    dep = data['dep'][:]
    #
    ## Write
    fmt = '%.4g'
    f   = open(fname,'a')
    
    #1)roughness for ux and uy
    np.savetxt(f,dep.T , fmt=fmt)
  
    f.close()        

