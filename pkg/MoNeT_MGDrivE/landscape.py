import math
import numpy as np
import scipy.stats as stats
#import vincenty as vn


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Distances
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def euclideanDistance(a, b):
    '''
    Calculates the Euclidean distance between two-dimensional coordinates.
    '''
    dist = math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2)
    return dist


def calculateDistanceMatrix(landscape, distFun=euclideanDistance):
    '''
    Returns the distance matrix according to the provided distance function.
        Examples of these are: euclideanDistance (xy), vn.vincenty (latlong).
    '''
    coordsNum = len(landscape)
    distMatrix = np.empty((coordsNum, coordsNum))
    for (i, coordA) in enumerate(landscape):
        for (j, coordB) in enumerate(landscape):
            distMatrix[i][j] = distFun(coordA, coordB)
    return distMatrix


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Kernels
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

def inverseLinearStep(distance, params=[.75, 1]):
    '''
    This function returns a migration estimate based on the inverse of the
        distance. NOTE: This is a terrible way to do it, but it's a first
        approximation. Should be replaced with the zero-inflated exponential.
    '''
    if math.isclose(distance, 0):
        return params[0]
    else:
        return (1 / (distance * params[1]))
    return True


def zeroInflatedExponential(distance, params=[0.01848777, 1.0e-10, math.inf]):
    '''
    Calculates the zero-inflated exponential for the mosquito movement kernel
        (default parameters set to Aedes aegypti calibrations).

        params = [rate, a, b]
    '''
    if(params[1] > params[2]):
        return None

    scale = 1.0/params[0]
    gA = stats.expon.cdf(params[1], scale=scale)
    gB = stats.expon.cdf(params[2], scale=scale)
    if np.isclose(gA, gB, 0.000001):
        return None

    densNum = stats.expon.pdf(distance, scale=scale)
    densDen = gB - gA

    return round(densNum/densDen, 10)


def migrationKernel(distMat, params=[.75, 1], kernelFun=inverseLinearStep):
    '''
    Takes in the distances matrix, zero inflated value (step) and two extra
        parameters to determine the change from distances into distance-based
        migration probabilities (based on the kernel function provided).
    '''
    coordsNum = len(distMat)
    migrMat = np.empty((coordsNum, coordsNum))
    for (i, row) in enumerate(distMat):
        for (j, dst) in enumerate(row):
            migrMat[i][j] = kernelFun(dst, params=params)
        # Normalize rows to sum 1
        migrMat[i] = migrMat[i] / sum(migrMat[i])
    return migrMat
