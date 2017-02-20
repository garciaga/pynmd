# -*- coding: utf-8 -*-
"""
A series of tools related to investments

Author:
-------
Gabriel Garcia Medina

Notes:
------
I am a novice economist thus do not take these formulae seriously.

Dependencies:
-------------
numpy
pandas

Internal dependencies:
----------------------

"""

from __future__ import division,print_function

# Import modules
import numpy as _np

#===============================================================================
# Compound interest
#===============================================================================
def compoundInterest(N,p=0.0,r=0.07,c=0.0):
    """
    Simple Formula to compute compound intrest with the option adding
    equal contributions during N
    
    PARAMETERS
    ----------
    N: Number of times to compound interest
    p: Principal (starting amount of money)
    r: Interest rate
    c: Contribution amount at each unit of time
    
    RETURNS:
    --------
    B: Ending balance
    
    Notes:
    ------
    
    """
    
    # Ending balance is principal growth rate plus contributions
    balance = (p * (1.0+r)**N + 
               c * (((1+r)**(N+1) - (1+r)) / r))
    
    return balance
    