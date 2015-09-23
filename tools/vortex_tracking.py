"""

Script to track surface vortices.

References:

McWilliams, J. C. 1990: The vortices of two-dimensional turbulence.
    Journal of Fluid Mechanics, 219, 361-385.
Long and Ozkan-Haller 2009: Low frequency characteristics of wave group-forced
    vortices. Journal of Geophysical Research, 114 (C08004).

    
# Dependencies:
  numpy
  
# TO DO:
  - Under certain conditions a vortex is counted twice (or more). Need to
    implement a filter for that.   
  - When searching for vortices (flagind) need to implement a boundary
    criterion to stop the search.  


# Version information:
v0.1   : First code development  
         Gabriel Garcia Medina (ggarcia@coas.oregonstate.edu), November 2014
v0.1.1 : Added centroid output and minor changes to documentation
         Gabriel Garcia Medina (ggarcia@coas.oregonstate.edu), January 2015    
"""

__author__ = "Gabriel Garcia Medina"
__email__ = "ggarcia@coas.oregonstate.edu"
__group__ = "Nearshore Modeling Group"
__version__ = '0.1.1'

# Import modules
import numpy as np


#===============================================================================
# Main vortex tracking algorithm
#===============================================================================
def vortex_tracking_main(x,y,q,ot,vthresh,maxpnt,stencil,boundary,minarea,
                         rcirc,verbose):
    '''
    
    Usage:
    coords = vortex_tracking(x,y,q,ot,vthresh,stencil,boundary,minarea,
                             rcirc,verbose)
    
        Parameters
        ----------
        x,y      : 2D arrays of coordinates
        q        : Array of vorticity computed at x,y points
        ot       : Time vector
        vthresh  : Tracking threshold 
        maxpnt   : Maximum number of points in a given vortex
        stencil  : Initial area filtering stencil 
        boundary : Width of the boundary (indices)
        minarea  : Minimum area to accept a vortex
        rcirc    : Circularity criterion {C/[2*(pi*area)**0.5]}
        verbose  : Print messages to the screen
    
        Output:
        ------
        coords   : List of dictionaries with vortex boundaries
        
        Fields in coords:
        -----------------
        ocean_time : time stamp [s]
        xind       : x indices of vortex boundary
        yind       : y indices of vortex boundary
        x          : x location of vortex boundary [m]
        y          : y location of vortex boundary [m]
        area       : vortex area [m2]
        rvort      : Vortex circularity {Circumference/[2*(pi*area)**0.5]}
        qe         : Vorticity extrema [s-1]
        xeind      : Vorticity extrema index in x direction
        yeind      : Vorticity extrema index in y direction
        xe         : x location of vorticity extrema [m]
        ye         : y location of vorticity extrema [m]
        cx         : x location of vorticity centroid [m]
        cy         : y location of vorticity centroid [m]
        
        
    '''
    
    # Initialize variables
    coords = []
    
    # Loop over time
    for aa in range(len(ot)):
        
        # Identify vorticity extrema
        qi = q[aa,:,:]
        # Get 2d indices of values that exceed threshold
        yind,xind = np.where(qi>vthresh)
        xe = x[yind,xind]
        ye = y[yind,xind]
        qe = qi[yind,xind]        
        
        # Sort from strongest to weakest
        sort_ind = np.flipud(np.argsort(qe))
        xe = xe[sort_ind]
        ye = ye[sort_ind]
        qe = qe[sort_ind]
        xind = xind[sort_ind]
        yind = yind[sort_ind]       
        
        # Area filtering (remove adjacent points that exceed threshold based on
        # the stencil size)
        indtmp = np.ones_like(xind,dtype=bool)
        for bb in range(len(indtmp)):
            if indtmp[bb]:
                xtmp = np.abs(xind - xind[bb])>stencil
                ytmp = np.abs(yind - yind[bb])>stencil
                xytmp = np.logical_and(xtmp,ytmp)
                indtmp[(bb+1):] = xytmp[(bb+1):] * indtmp[(bb+1):]
        # Remove adjacent points
        xe = xe[indtmp]
        ye = ye[indtmp]
        qe = qe[indtmp]
        xind = xind[indtmp]
        yind = yind[indtmp]
        
        
        # Remove points adjacent to the boundaries
        indtmp = np.logical_and(xind>=boundary,xind<(qi.shape[1]-boundary))
        xe = xe[indtmp]
        ye = ye[indtmp]
        qe = qe[indtmp]
        xind = xind[indtmp]
        yind = yind[indtmp]
        indtmp = np.logical_and(yind>=boundary,yind<(qi.shape[0]-boundary))
        xe = xe[indtmp]
        ye = ye[indtmp]
        qe = qe[indtmp]
        xind = xind[indtmp]
        yind = yind[indtmp]
        
        # If there is at least a point over threshold
        if len(qe) > 1:
            
            # Find the boundaries of the vortex
            # From largest to smallest intensity
            # Counter clockwise search will be employed
            
            # Looping over vortex points
            for bb in range(len(qe)):
                
                # Check if the evaluated vortex lies within the boundaries of a
                # previously computed vortex
                # TO DO SOON!!!
                
                # Get boundary of vortex
                # Returns False if the vortex is not found or does not meet the
                # selected criteria.
                bndtmp = get_boundary(x,y,qi,xind[bb],yind[bb],vthresh,
                                      maxpnt,minarea,rcirc,verbose)
                
                # Prepare output list
                if bndtmp is not False:
                                                           
                    # Add fields to dictionary
                    bndtmp['ocean_time'] = ot[aa]
                    bndtmp['qe'] = qe[bb]
                    bndtmp['xeind'] = xind[bb]
                    bndtmp['yeind'] = yind[bb]
                    bndtmp['xe'] = xe[bb]
                    bndtmp['ye'] = ye[bb]
                                                       
                    # Add to output matrix
                    coords.append(bndtmp)
                    
        else:
            if verbose:
                print("  No vorticity values exceed the selected threshold " + 
                      "of " + np.str(vthresh) + "s^{-1}")
                print("    after the area and boundary filters are applied.")
                print("    Field " + np.str(aa+1) + " of " + np.str(len(ot)))


    # End of function
    return coords
   

#===============================================================================
# Get boundary of vortex given a center point and a threshold
#===============================================================================
def get_boundary(x,y,q,xind,yind,vthresh,maxpnt,minarea,rcirc,verbose):
    '''
    
    Usage:
    get_boundary(x,y,q,xind,yind,vthresh,maxpnt,minarea,rcirc,verbose)
    
        Parameters:
        ----------
        x,y        : 2D coordinate arrays of vorticity locations (psi points)
        q          : Array of vorticity computed at x,y points
        xind,yind  : x and y indices of the vorticity extrema to track
        vthresh    : vorticity intensity to track
        maxpnt     : Maximum number of points in the circumference of vortex 
                     Mainly used to avoid infinite loops 
        minarea    : Minimum area to accept the vortex 
        rcirc      : Circularity criterion {C/[2*(pi*area)**0.5]}
        verbose    : Display output messages
        
        Output:
        -------
        bound      : Dictionary of closed boundary points surrounding a vortex.
                     Returns None if the boundary does not meet specific 
                     criteria (see codes for details).
    '''
    
    # Initialize variables
    out1 = [xind]               # Temporary output parameter 1 (x index)
    out2 = [yind]               # Temporary output parameter 2 (y index)
    out3 = [x[yind,xind]]       # Temporary output parameter 3 (x coordinate)
    out4 = [y[yind,xind]]       # Temporary output parameter 4 (y coordinate)
    exitflag = False            # Flag for loop around contour
    bound = True                # Output variable initialization
    counter = 0                 # Counter variable
    searchord = [0,3,2,1]       # Search order (right,down,left,up)

    # Go to the right boundary closest to the velocity extrema point and update
    # output arrays
    tmpind = np.argmin(q[yind,(xind+1):]>=vthresh)
    xind = xind + tmpind
    out1.append(xind)
    out2.append(yind)
    out3.append(x[yind,xind])
    out4.append(y[yind,xind])
        
    # Loop around the edge
    while exitflag is False:
        
        # Counter variable increase
        counter += 1
        brkflag = False
        
        # Start searching for vorticity values counter-clockwise from the index
        # points of the vorticity extrema
        for aa in range(5):

            if aa == searchord[0]:
                # Search to the right
                # Returns the first false index value 
                flagind = q[yind,xind+1]>=vthresh
                if flagind:
                    xind = xind + 1
                    brkflag = True                                  
            elif aa == searchord[1]:
                # Search down
                flagind = q[yind-1,xind]>=vthresh
                if flagind:
                    yind = yind - 1
                    brkflag = True
                    #searchord = [0,1,2,3]   # right, down, left, up                   
            elif aa == searchord[2]:
                # Search to the left
                flagind = q[yind,xind-1]>=vthresh
                if flagind:
                    xind = xind - 1
                    brkflag = True                    
            elif aa == searchord[3]:
                # Search upwards
                flagind = q[yind+1,xind]>=vthresh
                if flagind:
                    yind = yind + 1
                    brkflag = True                        
            else:
                # No suitable vorticity contour is found around the point being
                # considered.
                if verbose:
                    print("Isolated point found")
                exitflag = True
                bound = False
            
            # Debugging messages
            if verbose:    
                print('counter ' + np.str(counter))
                print('xind ' + np.str(xind))
                print('yind ' + np.str(yind))
                print('aa ' + np.str(aa))
                print('flagind ' + np.str(flagind))
                print('q[yind,xind] ' + np.str(q[yind,xind]))
            
            if brkflag:
                # Change search order
                if aa == 0:
                    # Shift 90 degrees right
                    searchord = [searchord[3],searchord[0],
                                 searchord[1],searchord[2]]
                elif aa == 2:
                    # Shift 90 degrees left
                    searchord = [searchord[1],searchord[2],
                                 searchord[3],searchord[0]]
                elif aa == 3:
                    # Shift 180 degrees
                    searchord = [searchord[2],searchord[3],
                                 searchord[0],searchord[1]]
                
                # Break loop
                break
        
        
        # Update the index dictionary 
        if exitflag is False:
            
            # Update output list
            out1.append(xind)
            out2.append(yind)
            out3.append(x[yind,xind])
            out4.append(y[yind,xind])
                
            # Check for maximum number of iterations
            if counter > maxpnt:
                if verbose:
                    print('Maximum number of boundary points (' + 
                          np.str(maxpnt) + ') have been exceeded')
                exitflag = True
                bound = False                              # Bug fixed v0.1.1
                
            # Verify if the new indices correspond to the first boundary points
            elif xind == out1[1] and yind == out2[1]:
                if verbose:
                    print('Back to the origin')
                # Remove initial (center) point
                out1 = out1[1:]
                out2 = out2[1:]
                out3 = out3[1:]
                out4 = out4[1:]
                exitflag = True
                
                
        
    # Employ area and shape checks ---------------------------------------------
    if bound:
        
        # Compute area of the vortex
        out1 = np.array(out1)
        out2 = np.array(out2)
        out3 = np.array(out3)
        out4 = np.array(out4)
        # Shoelace formula
        area = 0.5*(np.sum(out3[0:-1]*out4[1:]) - 
                    np.sum(out3[1:]*out4[0:-1]))
        if area < minarea:
            if verbose:
                print('Vortex rejected on area criterion')
                print('  Vortex area = ' + np.str(area) + 
                      'm2 and minimum area = ' + np.str(minarea) + 'm2\n')
            bound = False
        
    # Circumference test        
    if bound:
        perimeter = np.sum(((out3[1:]-out3[0:-1])**2 + 
                            (out4[1:]-out4[0:-1])**2)**0.5)
        rvort = perimeter/(2*np.sqrt(np.pi*area))
        if rvort > rcirc:
            if verbose:
                print('Vortex rejected based on circularity')
                print('  Rvortex = ' + np.str(rvort) + 
                      ' and Rmax = ' + np.str(rcirc) + '\n')
            bound = False
        
    # Compute centroid of the polygon (v0.1.1)
    if bound:
        cx = 1.0/(6.0*area)*np.sum((out3[0:-1] + out3[1:])*
                                   (out3[0:-1]*out4[1:] - 
                                    out3[1:]*out4[0:-1]))
        cy = 1.0/(6.0*area)*np.sum((out4[0:-1] + out4[1:])*
                                   (out3[0:-1]*out4[1:] - 
                                    out3[1:]*out4[0:-1]))
             
    # Wrap-up ------------------------------------------------------------------
    if bound:
        bound = {'xind':out1,'yind':out2,'x':out3,'y':out4,
                 'area':area,'rvort':rvort,'cx':cx,'cy':cy}
    return bound
       



#===============================================================================
# Get boundary of vortex given a center point and a threshold
#===============================================================================
def uniqueness(coords,maxspeed,maxdt,verbose=False,maxiter=10000):
    '''
    
    Usage:
    indices = uniqueness(coords,maxspeed,maxdt,verbose,maxiter)
    
        Parameters:
        ----------
        coords     : List of dictionaries containing the vortex boundary and
                     center information
        maxspeed   : Maximum speed allowed for vortices [m/s]
        maxdt      : Maximum time between sucessive vortices [s]
        verbose    : Display output messages (defaults to False)
        maxiter    : Maximum number of iterations (defaults to 10000)
        
        Output:
        -------
        indices    : List of lists that contain temporal indices for unique 
                     vortices
                     
    '''
       
    # Extract time and coordinate values from the vorticity extrema point 
    xe = []
    ye = []
    ot = []
    for aa in range(len(coords)):
        xe.append(coords[aa]['xe'])
        ye.append(coords[aa]['ye'])
        ot.append(coords[aa]['ocean_time'])
    
    # Loop forward in time tracking each vortex in the order found
    counter = 0                                 # Vortex counter variable
    iter = 0                                    # Emergency maximum iteration
    vortexflag = False                          # Flag for considered vortices 
    indices = []                                # List of lists
    indexedvort = np.zeros_like(xe,dtype=bool)  # Array of values to manipulate
                                                # counter variable
    tmpindices = [0]                            # Indices for current vortex
    indexedvort[0] = True                       # First one is indexed            
                                                    
    if verbose:
        print('\n')
        print('='*60)
        print('New vortex found')
                                                        
    # Loop forward from the second vortex center
    while vortexflag is not True:
        
        # Parent loop counter
        iter += 1
        samevort = True
        if verbose:
            print('  Iteration number ' + np.str(iter))            
        
        # Find vortex life and distance with respect to current evaluated
        # point. I should to filtering of the already considered vortex values.
        tmplife = ot[(counter+1):] - ot[counter]
        tmpspeed = (((((xe[(counter+1):] - xe[counter])**2) + 
                      (ye[(counter+1):] - ye[counter])**2)**0.5)/
                    tmplife)
        
        # Loop forward in time and space to find a vortex that satisfies the
        # life and speed limits. If this fails the counter variable will
        # increase and a new list will be generated to find the next vortex.
        for aa in range(len(tmplife)):

            # Verify vortex age and distance
            if (tmplife[aa] < maxdt) and (tmpspeed[aa] < maxspeed):
                # If true then this is the same vortex and the index will
                # be added to the temporary list                
                counter += aa + 1
                tmpindices.append(counter)
                indexedvort[counter] = True     
                if verbose:
                    print('  counter = ' + np.str(counter)) 
                break
            elif aa == (len(tmplife) - 1):
                # This is a new vortex
                samevort = False
        
        # Check if current vortex reached end of file
        # Search for new vortex (current one has already been allocated)
        if counter == (len(xe)-1):
            if verbose:
                print('End of file reached on current vortex')               
            samevort = False
                
        # If the evaluated vortex is a new vortex then add the previous 
        # vortex to the output list and reset the temporary list
        if samevort is not True:
            indices.append(list(tmpindices))
            counter = np.argmin(indexedvort)
            tmpindices = [counter]
            indexedvort[counter] = True
            
            if verbose:
                print('\n')
                print('='*60)
                print('New vortex found')
                print('  counter = ' + np.str(counter))
            
              
        # Check if all the vortices have been considered and allocate current
        # array
        if False not in indexedvort:
            if verbose:
                print('\n')
                print('='*60)
                print('All vortices have been indexed')           
            vortexflag = True        
            indices.append(list(tmpindices))
        # Check for iteration limits         
        elif iter > maxiter:
            if verbose:
                print('\n')
                print('='*60)
                print('Warning: Maximum number of iterations reached')
            vortexflag = True       
                    
                    
    
    # Print summary to screen --------------------------------------------------
    if verbose:
        print("A total of " + np.str(len(indices)) + 
              " different vortices found")
        
   
    # End of uniqueness code
    return indices

# End of module