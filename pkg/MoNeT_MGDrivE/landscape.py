import math
import numpy as np
#import vincenty as vn


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
