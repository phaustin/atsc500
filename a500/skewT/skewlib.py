"""
   construct a skewT - ln P diagram with dry adiabats
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticks

from a405.thermo.thermlib import convertTempToSkew, convertSkewToTemp, find_theta
from a405.thermo.constants import constants as c


def makeSkewDry(ax, skew=30):
    """       
      draw a skew-T lnP diagram on an axis

      Parameters
      ----------
      
      ax : matplotlib.axes
           matplotlib figure axis

      skew : float

             adjustable coefficient to make isotherms slope
             compared to adiabats

      Returns
      -------
      
      ax : matplotlib.axes
           the modified figure axis

      skew : float
             skew coefficient (K)

      """
    # majorFormatter = ticks.FormatStrFormatter('%d')
    # ax.yaxis.set_major_formatter(majorFormatter)
    yplot = range(1000, 190, -10)
    xplot = range(-300, -139)
    pvals = np.size(yplot)
    tvals = np.size(xplot)
    temp = np.zeros([pvals, tvals])
    theTheta = np.zeros([pvals, tvals])

    # lay down a reference grid that labels xplot,yplot points 
    # in the new (skewT-lnP) coordinate system .
    # Each value of the temp matrix holds the actual (data)
    # temperature label (in deg C)  of the xplot, yplot coordinate.
    # pairs. The transformation is given by W&H 3.56, p. 78.  Note
    # that there is a sign difference, because rather than
    # taking y= -log(P) like W&H, I take y= +log(P) and
    # then reverse the y axis

    for i in yplot:
        for j in xplot:
            # Note that we don't have to transform the y
            # coordinate, as it is still pressure.
            iInd = yplot.index(i)
            jInd = xplot.index(j)
            temp[iInd, jInd] = convertSkewToTemp(j, i, skew)
            Tk = c.Tc + temp[iInd, jInd]
            pressPa = i * 100.
            theTheta[iInd, jInd] = find_theta(Tk, pressPa)

            #
            # Contour the temperature matrix.
            #

            # First, make sure that all plotted lines are solid.
    mpl.rcParams['contour.negative_linestyle'] = 'solid'
    tempLabels = range(-40, 50, 10)
    ax.set_yscale('log')
    majorFormatter = ticks.FormatStrFormatter('%d')
    ax.yaxis.set_major_formatter(majorFormatter)
    ax.yaxis.set_minor_formatter(majorFormatter)
    plt.setp(ax.get_xticklabels(), weight='bold')
    plt.setp(ax.get_yticklabels(), weight='bold')
    plt.setp(ax.get_yticklabels(minor=True), weight='bold')
    tempLevs = ax.contour(xplot, yplot, temp, tempLabels, \
                          colors='k')
    #
    # Customize the plot
    #
    thetaLabels = list(range(200, 390, 10))
    thetaLevs = ax.contour(xplot, yplot, theTheta, thetaLabels, \
                      colors='b')
    # Transform the temperature,dewpoint from data coords to
    # plotting coords.
    ax.set_title('skew T - lnp chart')
    ax.set_ylabel('pressure (hPa)')
    ax.set_xlabel('temperature (deg C)')
    #
    # Crop image to a more usable size
    #
    TempTickLabels = range(-15, 40, 5)
    TempTickCoords = TempTickLabels
    skewTickCoords = convertTempToSkew(TempTickCoords, 1.e3, skew)
    ax.set_xticks(skewTickCoords)
    ax.set_xticklabels(TempTickLabels)

    skewLimits = convertTempToSkew([-15, 35], 1.e3, skew)

    ax.axis([skewLimits[0], skewLimits[1], 300, 1.e3])
    #
    # Create line labels
    #
    fntsz = 9  # Handle for 'fontsize' of the line label.
    ovrlp = True  # Handle for 'inline'. Any integer other than 0
    # creates a white space around the label.

    tempLevs.clabel(inline=ovrlp,
                    inline_spacing=0,
                    fmt='%2d',
                    fontsize=fntsz,
                    use_clabeltext=True)
    thetaLevs.clabel(inline=ovrlp,
                     inline_spacing=0,
                     fmt='%5d',
                     fontsize=fntsz,
                     use_clabeltext=True)
    ax.invert_yaxis()
    ax.yaxis.grid(True,which='minor')
    ax.figure.canvas.draw()
    return ax, skew

def plot_test():
    """
      create test plot when module is run
      """
    fig, ax = plt.subplots(1, 1)
    ax, skew = makeSkewDry(ax)
    figname = 'plot_test.png'
    fig.canvas.print_figure(figname)
    print('created {}'.format(figname))
    return ax


if __name__ == "__main__":
    ax = plot_test()
