
import numpy as np
import matplotlib.pyplot as plt
from numpy.random import uniform


def adjRSquared(rSquared, samplesNum, featuresNum):
    rAdj = 1-(1-rSquared)*(samplesNum-1)/(samplesNum-featuresNum-1)
    return rAdj


def unison_shuffled_copies(a, b, size=1000):
    assert len(a) == len(b)
    p = np.random.permutation(len(a))
    return a[p][:size], b[p][:size]


def getSamples_PDPICE(
        modelPredict, indVarIx, tracesNum=250,
        X=None, varRanges=None,
        indVarDelta=0.1, indVarStep=None
    ):
    # Auto-range variables sampling if needed ---------------------------------
    if varRanges:
        vRanges = varRanges
    elif (not varRanges) and not (X is None):
        minMax = [np.min(X, axis=0), np.max(X, axis=0)]
        vRanges = list(zip(*minMax))
    else:
        print("Error, please provide either the X vector (training inputs) or the list of variables ranges (varRanges)!")
        return {}
    # Get original sampling scheme (no X sweep yet) ---------------------------
    samples = np.zeros((len(vRanges), tracesNum))
    for (ix, ran) in enumerate(vRanges):
        samples[ix] = uniform(low=ran[0], high=ran[1], size=(tracesNum,))
    samples = samples.T
    # Get independent variable steps (x sweep) --------------------------------
    (rMin, rMax) = (vRanges[indVarIx][0], vRanges[indVarIx][1])
    stepType = type(indVarStep)
    if (stepType is list) or (stepType is tuple) or (stepType is np.ndarray):
        ivarSteps = np.array(indVarStep)
    else:
        step = indVarStep if indVarStep else (rMax-rMin)*indVarDelta
        ivarSteps = np.arange(rMin, rMax+step, step)
    # Evaluate model on steps -------------------------------------------------
    traces = np.zeros((samples.shape[0], ivarSteps.shape[0]))
    for six in range(samples.shape[0]):
        smpSubset = np.tile(samples[six], [ivarSteps.shape[0], 1])
        for (r, ivar) in enumerate(ivarSteps):
            smpSubset[r][indVarIx] = ivar
        yOut = modelPredict(smpSubset, verbose=False)
        traces[six] = yOut
    # Return dict -------------------------------------------------------------
    pdpice = {'x': ivarSteps, 'pdp': traces.T, 'ice': np.mean(traces, axis=0)}
    return pdpice


def plotPDPICE(
        pdpice, figAx=None,
        PDP=True, ICE=True, 
        YLIM=None, TITLE=None, 
        pdpKwargs={'color': '#a3cef155', 'ls': '-', 'lw': 0.15},
        iceKwargs={'color': '#ef476fff', 'ls': ':', 'lw': 3}
    ):
    if figAx:
        (fig, ax) = figAx
    else:
        (fig, ax) = plt.subplots(figsize=(5, 5))
    # Unpack variables --------------------------------------------------------
    (x, pdp, ice) = (pdpice['x'], pdpice['pdp'], pdpice['ice'])
    # Generate plots ----------------------------------------------------------
    if PDP:
        ax.plot(x, pdp, **pdpKwargs)
    if ICE:
        ax.plot(x, ice, **iceKwargs)
    # Axis and frame ----------------------------------------------------------
    ylim = YLIM if YLIM else (np.min(pdp.T), np.max(pdp.T))
    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(*YLIM)
    if TITLE:
        ax.set_title(TITLE)
    # Return figure -----------------------------------------------------------
    return (fig, ax)