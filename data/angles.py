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
            
         
# End of module
