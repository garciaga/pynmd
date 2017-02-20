# -*- coding: utf-8 -*-
"""
A series of tools related to retirement investments

Author:
-------
Gabriel Garcia Medina

Dependencies:
-------------
numpy
locale
time

Internal dependencies:
----------------------
ivestments
"""

from __future__ import division,print_function

# Import modules
import numpy as _np
import locale
import time as _time

# Internal modules
import investments as _invest

#===============================================================================
# Estimates for retirement accounts
#===============================================================================
def lifecycle(N,c,p=0.0,r=0.07,ir=0.025,ar=0.0,rR=0.0,
              currYear=None,cIAdj=False):
    """
    Estimate the growth and lifecycle of an account
    
    PARAMETERS:
    -----------
    N  : Number of times to compound interest
    p  : Principal (starting amount of money)
    r  : Interest rate
    c  : Contribution amount at each unit of time
         Can be float or array of size (N,)
    ir : Inflation rate (expected annual inflation rate)
    ar : Annual withdrawals from the account (today's dollars)
    rR : Interest rate during retirement. This accounts for choosing a more 
         conservative investment portfolio during retirement
    cIAdj : Adjust contributions to inflation
    
    RETURNS:
    --------
    tC: Time vector
    cT: Total contributions
    aP: Annual value of principal
    aW: Annual withrawings
    
    NOTES:
    ------
    1. My treatment of inflation is raw, be careful.
    2. A reasonable guess for USD inflation rate is 2.5%.
        
    """
    
    # Manage input
    tC = _np.arange(1,N+1) # Time vector
    if _np.size(c) == 1:
        c = _np.ones_like(tC) * c
    elif _np.size(c) != N:
        print('Number of entries in c is different than N')
        print('  Assuming yearly contributions of ' + _np.str(c[0]))        
        c = _np.ones_like(tC) * c[0]

    # Adjust contributions to inflation
    if cIAdj:
        tmpIAdj = (1.0+ir)**(_np.arange(N))
        c *= tmpIAdj
    
    # Preallocate variables
    rM = []                     # Money retired in today's dollars
    aP = []                     # Principal each year
    
    # First estimate the account growth while making contributions
    for aa in range(tC.shape[0]):
        
        if aa == 0:
            aP.append(_invest.compoundInterest(1.0,p=p,r=r,c=c[aa]))
        else:
            aP.append(_invest.compoundInterest(1.0,p=aP[aa-1],r=r,c=c[aa]))

    # Allocate contributions matrix
    c = _np.cumsum(c)
    
    # Value of today's dollars in the year of your last contribution
    tmpDV = (1.0-ir)**(N-1)
        
    print('Value of the account at retirement:')
    print('  Total Contributions: ' + 
          locale.format("%d",c[-1],monetary=True,grouping=True))
    print('  Raw:                 ' + 
          locale.format("%d",aP[-1],monetary=True,grouping=True))
    print('  Inflation Adjusted:  ' + 
          locale.format("%d",aP[-1]*tmpDV,monetary=True,grouping=True))
    
    # Start retiring money until it runs out -------------
    tmpAP = _np.copy(aP[-1])
    
    tmpW = [] # Allocate withdrawings
    aa = 0
    richFlag = False # True if money will not run out
    while tmpAP > 0.0:
        
        # Increase counter variable
        aa += 1
        
        # Adjust the price of today's dollar for inflation
        tmpDV *= (1.0-ir) # Need to figure this out later
        
        # Find how much money you are actually withdrawing
        tmpA = ar / tmpDV
        tmpW.append(tmpA)
        
        # Withdraw the money
        tmpAP = aP[-1] - tmpA        
        
        if tmpAP < 0:
            aP.append(0.0)
            break
        else:
            # Principal keeps getting interest
            aP.append(_invest.compoundInterest(1.0,p=tmpAP,r=rR,c=0.0))
        
        if aP[-1] > aP[-2]:
            print('Money supply will not run out')
            return tC,c,_np.array(aP),tmpA

    # Shift the time vectors        
    if not currYear:
        currYear = _time.localtime()[0]

    # Allocate retirement year matrix
    tR = _np.arange(N+1,N+aa+1)
    aW = _np.r_[_np.zeros_like(tC),_np.array(tmpW)]
    tC = _np.r_[tC,tR] + (currYear-1)
    c  = _np.r_[c,_np.zeros_like(tR)]
    aP = _np.array(aP)
    
    # Summary of investments
    print('Money supply will last for ' + _np.str(aa) + ' years after last ' + 
          'contribution')
    print('Withdrawals:')
    print('  ' + '{:4d}'.format(tC[N]) + ' ' + 
          locale.format("%d",aW[N],monetary=True,grouping=True))
    print('  ' + '{:4d}'.format(tC[-2]) + ' ' + 
          locale.format("%d",aW[-2],monetary=True,grouping=True))
    print('  ' + '{:4d}'.format(tC[-1]) + ' ' + 
          locale.format("%d",aW[-1],monetary=True,grouping=True))

    # Exit the function
    return tC,c,aP,aW
    