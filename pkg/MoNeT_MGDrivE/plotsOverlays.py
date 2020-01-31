

def removeTicksAndLabels(ax):
    """Deletes all the ticks and labels from the plot axes passed as argument.

    Parameters
    ----------
    ax: Matplotlib axes.
        Can be obtained from a figure with: fig.get_axes()[0]

    Returns
    -------
    ax
        Returns the changed ax with the modification performed.

    """
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_xticks([])
    ax.set_yticks([])
    return ax


def setRange(ax, style):
    """Re-scales the plot to the range stated on the 'style' element.

    Parameters
    ----------
    ax : Matplotlib axes
        Can be obtained from a figure with: fig.get_axes()[0]
    style : Plot-style dictionary.
        Should contain:
            style['xRange'] = (lo, hi)
            style['yRange'] = (lo, hi)

    Returns
    -------
    ax
        Returns the changed ax with the modification performed.

    """
    ax.set_xlim(style['xRange'][0], style['xRange'][1])
    ax.set_ylim(style['yRange'][0], style['yRange'][1])
    return ax


def printVLines(ax, xCoords, width, color, alpha, lStyle):
    """Adds vertical lines to the plot in the list of x-coordinates passed.
        Created as an auxiliary function to mark days at which events happen
        (threshold crossings).

    Parameters
    ----------
    ax : Matplotlib axes
        Can be obtained from a figure with: fig.get_axes()[0]
    xCoords : type
        Description of parameter `xCoords`.
    width : type
        Description of parameter `width`.
    color : type
        Description of parameter `color`.
    alpha : type
        Description of parameter `alpha`.
    lStyle : type
        Description of parameter `lStyle`.

    Returns
    -------
    type
        Description of returned object.

    """
    for vLine in xCoords:
        ax.axvline(
                x=vLine, linewidth=width, linestyle=lStyle,
                color=color, alpha=alpha
            )
    return ax
