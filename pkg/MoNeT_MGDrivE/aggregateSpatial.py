
import numpy as np


def filterFilesByIndex(files, ix, male=True, female=True):
    m = [files['male'][z] for z in ix] if male else []
    f = [files['female'][z] for z in ix] if male else []
    ffiles = {'male': m, 'female': f}
    return ffiles


def filterGarbageByIndex(landRepetition, indices):
    return list(map(landRepetition.__getitem__, indices))


def filterAggregateGarbageByIndex(landscapeReps, indices):
    genes = landscapeReps['genotypes']
    repsNumber = len(landscapeReps['landscapes'])
    traces = []
    for j in range(0, repsNumber):
        probe = landscapeReps['landscapes'][j]
        trace = np.sum(filterGarbageByIndex(probe, indices), axis=0)
        traces.append([trace])
    filteredLand = {'genotypes': genes, 'landscapes': traces}
    return filteredLand
