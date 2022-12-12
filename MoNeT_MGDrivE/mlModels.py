
import numpy as np

def adjRSquared(rSquared, samplesNum, featuresNum):
    rAdj = 1-(1-rSquared)*(samplesNum-1)/(samplesNum-featuresNum-1)
    return rAdj


def unison_shuffled_copies(a, b, size=1000):
    assert len(a) == len(b)
    p = np.random.permutation(len(a))
    return a[p][:size], b[p][:size]