# -*- coding: utf-8 -*-
"""
A series of tools for data clustering

Authors:
-------
Gabriel Garcia Medina
    Nearshore Modeling Group
    ggarcia@coas.oregonstate.edu


Dependencies:
-------------
    numpy
    itertools
    
Internal dependencies:
----------------------
    signal

"""

from __future__ import division,print_function

# Python modules
import numpy as np
import itertools

# Internal dependencies
import signal as gsignal

#===============================================================================
# Simple Maximum Dissimilarity Algorithm
#===============================================================================
def mda_simp(x,numClust,dirvar,seed=None):
    """
    Simple maximum dissimilarity algorithm implementation
    
    PARAMETERS:
    -----------
    x        : array-like data of shape [numSamples,numCharacteristics] 
               to be clustered
    numClust : Number of clusters desired
    dirvar   : Boolean array of shape [numCharacteristics] which is true for
               directional variables only.
    seed     : (optional) index of the data to be used to seed the algorithm.
               If it is not passed a randomly selected value within x will be 
               used.
               
    RETURNS:
    --------
    clust    : array-like of shape [numClust,numCharacteristics] of the points
               identified by the algorithm.
    
    EXAMPLE
    -------
    >> import numpy as np
    >> import pynmd.data.clustering as gclust
    >> x = np.c_[waveHeight,wavePeriod,waveDirection]
    >> numClust = 10
    >> dirvar = np.array([0,0,1],dtype=np.bool)
    >> seed = np.argmax(waveHeight)
    >> clusters = gclust.clustering(x,numClust,dirvar,seed)
    
    """
    
    # Preallocate variables
    clust = np.zeros((numClust,x.shape[1])) * np.NAN
    
    # Seed 
    if seed:
        clustInd = seed
        clust[0,:] = x[seed,:]
    else:
        clustInd = np.int64(np.round(np.random.rand() * x.shape[0]))
        clust[0,:] = x[clustInd,:]
    
    # Remove seed from input data
    allInd = np.arange(0,x.shape[0],1) != clustInd
    x = x[allInd,:]
    
    # Find number of characteristics
    numChar = x.shape[1]
    
    # Find range of data for normalizing
    normVals = np.zeros((numChar,)) * np.NAN
    for aa in range(numChar):
        # Directional variables
        if dirvar[aa]:
            normVals[aa] = 180.0
        # Scalar variables        
        else:
            normVals[aa] = np.max(x[:,aa]) - np.min(x[:,aa])
    
    # Preallocate difference matrix
    distance = np.ones((x.shape[0],)) * 10.0**20 
    
    # Loop over desired clusters
    for aa in range(1,numClust):
        
        # Find the normalized distance between the other elements in the array -----
        dif = np.zeros_like(x) * np.NAN
        
        # Loop over characteristics
        for bb in range(numChar):
            
            if dirvar[bb]:
                # Directional varialbles
                tmp1 = np.abs(x[:,bb] - clust[aa-1,bb])
                tmp2 = 360.0 - tmp1
                dif[:,bb] = np.min(np.c_[tmp1,tmp2],axis=1)/normVals[bb]
                del tmp1,tmp2
            else:
                # Scalar variables
                dif[:,bb] = (x[:,bb] - clust[aa-1,bb])/normVals[bb] 
    
        # Square the differences and sum over them
        distance = np.min(np.c_[distance,np.sum(dif**2,axis=1)],axis=1)
        
        # Find the farthest most point ---------------------------------------------
        clustInd = np.argmax(distance)
        
        # Allocate the new cluster
        clust[aa,:] = x[clustInd,:]
        
        # Remove current cluster from array
        allInd = np.arange(0,x.shape[0],1) != clustInd
        x = x[allInd,:]
        
        # Remove the seed from the distance matrix
        distance = distance[allInd]
    
    return clust


#===============================================================================
# Equal probability clustering
#===============================================================================

def eqProbClust(x,numBin):
    """
    Clustering data by keeping the probability of occurrence constant
    
    PARAMETERS:
    -----------
    x       : Array of data to cluster of size [numSamples,numCharacteristics]
    numBin  : Array of number of bins per characteristic
    
    RETURNS:
    --------
    groupNo   : Group number of each data point  
    edges     : Edges of bins for all groups
    centroids : Centroids of all the identified groups
    
    
    NOTES:
    ------
    - The data will be clustered in the order of the characteristics
    
    """
        
    # Preallocate variables ----------------------------------------------------
    edges     = np.zeros((np.cumprod(numBin)[-1],numBin.shape[0]*2))
    groupNo   = np.zeros((x.shape[0],)) * np.NAN
    groupComb = np.array(list(itertools.product(*[range(aa) for aa in numBin])))
    centroids = np.zeros_like(groupComb) * np.NAN
    
    # Loop over combinations for binning
    for aa in range(groupComb.shape[0]):
        
        # Bin in all dimensions
        dataInd = np.ones((x.shape[0],),dtype=np.bool)    
        
        for bb in range(x.shape[1]):  
            
            # Get edges to bin data
            tmpEdges = eqProbBins(x[dataInd,bb],numBin[bb])
            
            # Store edges
            edges[aa,bb*2:bb*2+2] = tmpEdges[groupComb[aa][bb]:
                                             groupComb[aa][bb]+2]
            
            # Update data index
            if groupComb[aa,bb] == numBin[bb] - 1:
                tmpInd = np.logical_and(x[:,bb]>=edges[aa,bb*2],
                                        x[:,bb]<=edges[aa,bb*2+1])
            else:
                tmpInd = np.logical_and(x[:,bb]>=edges[aa,bb*2],
                                        x[:,bb]<edges[aa,bb*2+1])
            dataInd *= tmpInd
            
        # Allocate group number
        groupNo[dataInd] = aa
        
        # Get cluster centers
        centroids[aa,:] = np.mean(x[dataInd,:],axis=0)
    
    return groupNo,edges,centroids


#===============================================================================
# Equal Probability binning
#===============================================================================
def eqProbBins(x,numBin):
    """
    Function to get bin edges that will achieve equal probability
    
    PARAMETERS:
    -----------
    x      : 1-D array
    numBin : Number of bins desired 
    
    RETURNS:
    --------
    edges  : Bin edges that ensure equal probability of occurrence
     
    """
    
    # Compute the empirical cumulative distribution for interpolation purposes
    sortVar,probVar,sortInd = gsignal.ecdf(x)
    
    # Find the edges for equal probability clustering
    edges = np.zeros((numBin+1,))
    edges[0]  = sortVar.min() # Minimum data value is the default cluster minima
    edges[-1] = sortVar.max() # Maximum data value is the default cluster maxima
    
    # Find equal probability for intermediate bins
    prob = np.arange(1,numBin)/numBin    
    edges[1:-1] = np.interp(prob,probVar,sortVar)
    
    return edges
