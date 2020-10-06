#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import operator as op
import MoNeT_MGDrivE.summaryStatistics as sstat


def getPopRepsRatios(base, trace, gIx):
    """
    Given a baseline population (mean) and the traces of an experiment (reps)
        this function calculates the fraction (repRto) of a given genotype
        index (gIx).
    Args:
        base (sum): Mean pop array
        trace (srp): Repetitions pop array
        gIx (int): Genotype index
    Returns:
        np.array: Fractions array.
    """
    (basePop, tracePops) = (base['population'], trace['landscapes'])
    ratioReps = [sstat.getPopRatio(trace, basePop, gIx) for trace in tracePops]
    ratioArr = np.asarray(ratioReps)
    return ratioArr


def compRatioToThreshold(repsRatios, thld, cmprOp=op.lt):
    """
    Given a population fraction, this function performs the comparison
        against a threshold and returns the bool array.
    Args:
        repsRatios (np.array): Fraction repetitions array (x: time, y: rep)
        thiS (float): Threshold to make the introgression comparison (0 to 1)
    Returns:
        np.array: bool array with info of days in which the comparison is true.
    """
    thresholdArray = np.apply_along_axis(cmprOp, 0, repsRatios, thld)
    return thresholdArray


def calcTTI(repRto, thiS, clipValue=None):
    """
    Given a population fraction, this function calculates the time it takes
        for the system to go below a given threshold.
    Args:
        repRto (np.array): Fraction repetitions array (x: time, y: rep)
        thiS (float): Threshold to make the introgression comparison (0 to 1)
    Returns:
        list: Times at which each of the repetitions goes below the given
                threshold (first time it breaks it).
    """
    if clipValue is None:
        (_, days) = repRto.shape
    else:
        days = clipValue
    thiSBool = [compRatioToThreshold(repRto, i, op.lt) for i in thiS]
    ttiS = [np.argmax(thiBool == 1, axis=1) for thiBool in thiSBool]
    clipped = [[x if x > 0 else days for x in i] for i in ttiS]
    return clipped


def calcTTO(repRto, thoS):
    """
    Given a population fraction, this function calculates the time it takes
        for the system to go above a given threshold.
    Args:
        repRto (np.array): Fraction repetitions array (x: time, y: rep)
        thoS (float): Threshold to make the outrogression comparison (0 to 1)
    Returns:
        list: Times at which each of the repetitions goes above the given
                threshold (last time it breaks it).
    """
    (reps, days) = repRto.shape
    thoSBool = [compRatioToThreshold(repRto, i, op.gt) for i in thoS]
    ttoS = [np.subtract(days, np.argmin(np.flip(thoBool), axis=1)) for thoBool in thoSBool]
    return ttoS


def calcWOP(repRto, thwS):
    """
    Given a population fraction, this function calculates number of days in
        which the repetition is below a given threshold (window of protection).
    Args:
        repRto (np.array): Fraction repetitions array (x: time, y: rep)
        thwS (float): Threshold to make the outrogression comparison (0 to 1)
    Returns:
        list: Sum of the total of days the threshold is met (not the same as
            TTO-TTI for time-varying population sizes).
    """
    thwSBool = [compRatioToThreshold(repRto, i, op.lt) for i in thwS]
    wopS = [np.sum(thwBool, axis=1) for thwBool in thwSBool]
    return wopS


def calcMinMax(repRto):
    """
    Given a population fraction, this function calculates minimum/maximum
        fractions and the days at which they happen.
    Args:
        repRto (np.array): Fraction repetitions array (x: time, y: rep)
        thwS (float): Threshold to make the outrogression comparison (0 to 1)
    Returns:
        tuple: Min/max times and values.
    """
    (mni, mxi) = (repRto.min(axis=1), repRto.max(axis=1))
    mnx = np.asarray([np.where(repRto[i] == mni[i])[0][0] for i in range(len(mni))])
    mxx = np.asarray([np.where(repRto[i] == mxi[i])[0][0] for i in range(len(mxi))])
    return ((mni, mnx), (mxi, mxx))


def getRatioAtTime(repRto, ttpS):
    """
    Given a population fraction, returns the populations at a given time.
    Args:
        repRto (np.array): Fraction repetitions array (x: time, y: rep)
        ttpS (int): Day at which we want to probe the population
    Returns:
        np.array: Population values.
    """
    return np.asarray([repRto[:, ttp] for ttp in ttpS])


def initDFsForDA(
            fPaths, header, thiS, thoS, thwS, ttpS,
            peak=['min', 'minx', 'max', 'maxx']
        ):
    fNum = len(fPaths)
    heads = [list(header)+i for i in (thiS, thoS, thwS, ttpS, peak)]
    DFEmpty = [pd.DataFrame(int(0), index=range(fNum), columns=h) for h in heads]
    return DFEmpty


def filterDFWithID(df, xpid):
    xpidz = list(zip(list(df.columns)[:7], xpid))
    filters = [df[i[0]] == i[1] for i in xpidz]
    filter = list(map(all, zip(*filters)))
    return df[filter]


def loadDFFromSummary(fName):
    df = pd.read_csv(fName)
    header = list(df.columns)
    indRan = sum([i[0] == 'i' for i in header])
    headerInd = header[:indRan]
    return (df, header, headerInd)
