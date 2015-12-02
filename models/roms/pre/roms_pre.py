
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
    
    import numpy as np
    #S3D = 1 # 1 for 3D simulation, 0 for 2D simulation   

    g   = 9.81       # m/s2
    h   = maxdep * 1.0        # maximum depth in the domain (m)
    dt1 = dt_3d  * 1.0        # baroclinic time step if 3D, barotropic time step if 2D
    dt2 = num_dt_2d  * 1.0

    if S3D:
        print ' Evaluate a simple CFL check for baroclinic run'
        dt  = 1.0*dt1/dt2  
    else:
        print ' Evaluate a simple CFL check for 2D run'
        dt = dt1

    dx  = dxx * 1.0 # meters
    dy  = dyy * 1.0 # meters

    c  = np.sqrt(g*h) # m/s
    Cu = c*dt*np.sqrt((1./dx**2)+(1./dy**2))

    # I am not sure what gamma is exactly, I could not find a very good 
    # explanation for it but, apparently it is better if it is greater than 0.
    gamma = 0.284*(1.-10./dt1)

    print '    > The Courant number = ',Cu
    print '    > Gamma              = ',gamma
    print '    > nHis   (30 min)    = ',30*60/ dt1
    print '    > nHis   (15 min)    = ',15*60/ dt1
    print '    > ninfo  ( 1 min)    = ',1*60 / dt1
    print '    > ndefHis(daily)     = ',24*60*60/dt1


