import numpy as np
import operator as op


def reachedSteadtStateAtDay(
            aggData,
            safety=.01,
            finalFrame=-1
        ):
    """
    Description:
        * Calculates the point at which the aggregated alleles reach their
            steady state (defined by the final frame of the simulation).
            * aggData: Genotypes aggregated data.
    In:
        * safety: Envelope of values around the steady state that are
            considered "stable" (as a proportion of the final total allele
            composition).
        * finalFrame: Index of the day to be considered as "stable" for
            reference purposes.
    Out:
        * steadtStateReach: Day at which the dynamics of the system became
            stable.
    Notes:
        * This is a quick and simple way to calculate the steady state, but
            assumes that the system stabilizes and that the simulation was
            run for enough time for the dynamics to become stable.
    """
    finalFrame = aggData["population"][finalFrame]
    tolerance = round(sum(finalFrame) * safety)
    toleranceUp = finalFrame + tolerance
    toleranceDown = finalFrame - tolerance

    daysMax = len(aggData["population"])

    for i in range(0, daysMax):
        steadyStateReach = daysMax
        testFrame = aggData["population"][i]

        boolsUp = testFrame < toleranceUp
        boolsDown = testFrame > toleranceDown
        zeros = testFrame < 1

        if (all(boolsUp) and all(boolsDown)) or all(zeros):
            steadyStateReach = i
            break

    return steadyStateReach


def getTimeToMin(
            aggData,
            safety=.01
        ):
    """
    Description:
        * Calculates the point at which the total population reaches
            its minimum.
    In:
        * aggData: Genotypes aggregated data.
        * safety: Envelope of values around the steady state that are
            considered "stable" (as a proportion of the final total allele
            composition).
    Out:
        * time: Point in time at which the minimum is reached
        * popMin: Population size at its minimum
    """
    pop = [sum(row) for row in aggData['population']]
    for time in range(len(pop)):
        popMin = min(pop)
        if np.isclose(pop[time], popMin, atol=safety):
            break
    return (time, popMin)


def comparePopToThresholds(pop, gIx, tIx, thrs, cmprOp=op.lt):
    """Calculates if the genotypes at a desired index meet the condition
        passed as fractions of the total population. This function was created
        to calculate where the population goes below a given threshold for
        times of suppression.

    Parameters
    ----------
    pop : numpy array
        Genotypes numpy array of a signle population.
    gIx : int
        Index of the genotype of interest.
    tIx : list of integers
        List of the genotypes for the 'total population' calculation (for
            fractions to be computed).
    thrs : list of floats (0 to 1)
        List of the thresholds to compare the population fraction against.
    cmprOp : operator
        Comparsion to be perfomed between pop and thrs
            (https://docs.python.org/3/library/operator.html).

    Returns
    -------
    type
        Numpy array of bools flagging conditions being met. Each column
            represents one of the thresholds passed in tIx for the whole
            duration (time) of the sim (pop array's length).
    """
    flagsArray = np.empty((len(pop), len(thrs)), dtype=bool)
    for (i, dayData) in enumerate(pop):
        totalPop = sum(dayData[tIx])
        if (totalPop > 0):
            fraction = (dayData[gIx] / totalPop)
        else:
            fraction = dayData[gIx]
        closeFlags = [cmprOp(fraction, i) for i in thrs]
        flagsArray[i] = closeFlags
    return flagsArray
