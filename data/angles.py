"""

Functions to manipulate angles
   
# Dependencies:
  numpy
  
"""

__author__ = "Gabriel Garcia Medina"
__email__ = "ggarcia@coas.oregonstate.edu"
__group__ = "Nearshore Modeling Group"

# Import modules
import numpy as np

def wrapto360(x):
    '''
    Wraps angles to 360 degrees.
    '''
    
    # Make sure you are working with floating point numbers
    x = np.double(x)
    
    # Identify angles larger than 360
    x = np.mod(x,360.0)
    
    # Return wrapped dimensions  
    return(x)
    
    
def wrapto2pi(x):
    '''
    Wraps angles from 0 to 2*pi
    '''
    
    # Make sure you are working with floating point numbers
    x = np.double(x)
    
    # Identify angles larger than 2*pi
    x = np.mod(x,2*np.pi)
    
    return(x)
    

def wrapto180(x):
    '''
    Wraps angles from -180 to 180
    '''
    
    # Make sure you are working with floating point numbers
    x = np.double(x) + 180.0
    
    # Identify angles larger than 2*pi
    x = np.mod(x,360.0) - 180.0
    
    return(x)


def wraptopi(x):
    '''
    Wraps angles from -pi to pi
    '''
    
    # Make sure you are working with floating point numbers
    x = np.double(x) + np.pi
    
    # Identify angles larger than 2*pi
    x = np.mod(x,2*np.pi) - np.pi
    
    return(x)
            
def cartToNautDeg(x):
    """
    Converts from cartesian to nautical convention and wraps to 360

    If the direction that waves are travelling to is given it will return
    the nautical equivalent of that. In other words the output would be
    direction waves are traveling to in nautical coordinates.
    """

    x = 90.0 - np.double(x)

    return wrapto360(x)
    
def decimalDegreesToDMS(x,nDec=0):
    """
    Converts decimal degrees to Degree, Minute, Seconds

    PARAMETERS:
    -----------
    x:    Decimal degrees
    nDEC: Number of decimals in the seconds

    RETURNS:
    --------
    d: Degrees (integer)
    m: Minutes (integer)
    s: Seconds (Float rounded at nDec)
    """
    
    d = np.int(x)
    m = (x - d)*60
    s = np.round((m - np.int(m))*60,decimals=nDec)
    m = np.int(m)

    return d,m,s
    