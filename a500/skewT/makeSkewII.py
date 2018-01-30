"""
   construct a skewT - ln P diagram with rsat and thetaes lines
"""
from importlib import reload
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import a405thermo.thermlib
reload(a405thermo.thermlib)
from a405thermo.thermlib import convertTempToSkew, convertSkewToTemp, find_theta, find_thetaet
from a405thermo.constants import constants as c
from a405thermo.thermlib import find_rsat


def find_corners(temps, press=1.e3, skew=30.):
    """
        convert a pair of temperatures into skew plotting coordinates at a constant pressure level

        Parameters
        ----------

        temps : [float]
                xaxis corner temperatures (degC)

        press : float
                pressure level for skew conversion

        skew  : float
                coefficient to use in the tempToSkew conversion

        Returns
        -------

        corners : [float]
                x axis corners in plotting coordinates

      """
    corners = convertTempToSkew(temps, press, skew)
    return list(corners)


def makeSkewWet(ax, corners=[-30, 25], skew=30):
    """       
      draw a skew-T lnP diagram on an axis

      Parameters
      ----------
      
      ax : matplotlib.axes
           matplotlib figure axis

      corners : [float]
                x axis temperature limits (degC)

      skew : float

             adjustable coefficient to make isotherms slope
             compared to adiabats

      Returns
      -------
      
      ax : matplotlib.axes
           the modified figure axis

      """
    yplot = range(1000, 190, -6)  #
    xcorners = find_corners(corners, skew=skew)
    xplot = list(np.linspace(xcorners[0], xcorners[1], 45))
    pvals = np.size(yplot)
    tvals = np.size(xplot)
    temp = np.zeros([pvals, tvals])
    theTheta = np.zeros_like(temp)
    the_rsat = np.zeros_like(temp)
    theThetae = np.zeros([pvals, tvals])

    # lay down a reference grid that labels xplot,yplot points 
    # in the new (skewT-lnP) coordinate system .
    # Each value of the temp matrix holds the actual (data)
    # temperature label (in deg C)  of the xplot, yplot coordinate.
    # pairs. The transformation is given by W&H 3.56, p. 78.  Note
    # that there is a sign difference, because rather than
    # taking y= -log(P) like W&H, I take y= +log(P) and
    # then reverse the y axis

    for presshPa in yplot:  #loop over pressures
        for skewed in xplot:  #loop over skewed xcoords
            # Note that we don't have to transform the y
            # coordinate, as it is still pressure.
            iInd = yplot.index(presshPa)
            jInd = xplot.index(skewed)
            temp[iInd, jInd] = convertSkewToTemp(skewed, presshPa, skew)
            Tk = c.Tc + temp[iInd, jInd]
            pressPa = presshPa * 100.
            theTheta[iInd, jInd] = find_theta(Tk, pressPa)
            rs = find_rsat(Tk, pressPa)
            the_rsat[iInd, jInd] = rs
            theThetae[iInd, jInd] = find_thetaet(Tk,rs,Tk, pressPa)
    #
    # Contour the temperature matrix.
    #

    # First, make sure that all plotted lines are solid.
    mpl.rcParams['contour.negative_linestyle'] = 'solid'
    tempLabels = range(-40, 50, 10)
    ax.contour(xplot, yplot, temp, tempLabels, \
                          colors='k')
    #
    # contour theta
    #
    thetaLabels = list(range(200, 390, 10))
    thetaLevs = ax.contour(xplot, yplot, theTheta, thetaLabels, \
                      colors='b')
    #
    # contour rsat
    #
    rsLabels = [0.1, 0.25, 0.5, 1, 2, 3] + list(range(4, 28, 2)) #+ [26,  28]
    rsLevs = ax.contour(xplot,
                        yplot,
                        the_rsat * 1.e3,
                        levels=rsLabels,
                        colors='g',
                        linewidths=.5)

    thetaeLabels = np.arange(250, 410, 10)
    thetaeLevs = ax.contour(xplot, yplot, theThetae, thetaeLabels, \
                      colors='r')
    #
    # Customize the plot
    #
    ax.set_yscale('log')
    locs = np.array(range(100, 1100, 100))
    labels = locs
    ax.set_yticks(locs)
    ax.set_yticklabels(labels)  # Conventionally labels semilog graph.
    ax.set_ybound((200, 1000))
    plt.setp(ax.get_xticklabels(), weight='bold')
    plt.setp(ax.get_yticklabels(), weight='bold')
    ax.yaxis.grid(True)

    ax.set_title('skew T - lnp chart')
    ax.set_ylabel('pressure (hPa)')
    ax.set_xlabel('temperature (deg C)')

    #
    # Crop image to a more usable size
    #

    TempTickLabels = range(-30, 40, 5)

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

    #tempLevs.clabel(inline=ovrlp, inline_spacing=0,fmt='%2d', fontsize=fntsz,use_clabeltext=True)
    thetaLevs.clabel(inline=ovrlp,
                     inline_spacing=0,
                     fmt='%3d',
                     fontsize=fntsz,
                     use_clabeltext=True)
    rsLevs.clabel(inline=ovrlp,
                  inline_spacing=0,
                  fmt='%3.2g',
                  fontsize=fntsz,
                  use_clabeltext=True)
    thetaeLevs.clabel(thetaeLabels,
                      inline_spacing=0,
                      inline=ovrlp,
                      fmt='%5g',
                      fontsize=fntsz,
                      use_clabeltext=True)

    ax.invert_yaxis()
    #ax.figure.canvas.draw()
    xcorners = find_corners(corners, skew=skew)
    ax.set(ylim=(1000, 300), xlim=xcorners)
    return ax, skew


def plot_test():
    """
      create test plot when module is run
      """
    fig, ax = plt.subplots(1, 1)
    corners = [-25, 30]
    ax, skew = makeSkewWet(ax, corners=corners, skew=25)
    figname = 'plot_test_wet.pdf'
    fig.canvas.print_figure(figname)
    print('created {}'.format(figname))
    return ax


if __name__ == "__main__":
    plot_test()
