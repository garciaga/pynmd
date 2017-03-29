#!/usr/bin/python
"""
Convert physical quantities

"""
__author__ = "Saeed Moghimi"
__copyright__ = "Copyright 2017, NOAA/UCAR"
__license__ = "GPL"
__version__ = "1.0"
__email__ = "moghimis@gmail.com"



def Chezy_to_z0(C, H):
    """
    Convert roughness height Chezy "C" to z_0 which is a
    function of the water depth

    Inputs:
        C = Chezy "C" (non-dimensional)
        H = water depth (meters) (can be vector)
    Outputs:
        z0 = roughness height (meters)

    """
    k_s = 12 * H / (10 ** (C/18.0))
    z0 = k_s / 30.0
    return z0



