import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import tools.configuration_file as conf

scaling_factor = 50

def plotTube(ax, crossSection, velocity, pressure, dx, t):
    radius0 = np.sqrt(crossSection/np.pi)
    N = velocity.shape[0]
    ax.plot(np.arange(N) * dx, conf.r0 + (radius0 - conf.r0) * scaling_factor, 'k')
    ax.plot(np.arange(N) * dx, -(conf.r0 + (radius0 - conf.r0) * scaling_factor), 'k')
    iii = 0
    rects = []
    map = plt.get_cmap('RdBu')
    for x in np.arange(N)*dx:
        dy = (conf.r0 + (radius0[iii] - conf.r0) * scaling_factor)
        rect = Rectangle((x-.5*dx,-dy), dx, 2 * dy, color=map((velocity[iii] + conf.u0) / conf.ampl))
        ax.add_patch(rect)
        iii += 1
        rects.append(rect)


    #plt.quiver(np.arange(N+1)*dx,np.zeros_like(velocity),velocity,np.zeros_like(velocity))
    #plt.imshow(np.vstack((velocity,velocity,velocity,velocity)),origin="lower")
    #plt.imshow(np.vstack((velocity,velocity)),origin="upper")
    ax.set_ylim([-2,2])


def plotVar(ax, crossSection, dx, t):
    radius0 = np.sqrt(crossSection/np.pi)
    radius_mean = np.mean(np.sqrt(crossSection/np.pi))
    N = crossSection.shape[0]
    plt.plot(np.arange(N) * dx, (radius_mean-radius0)*scaling_factor)
    lim = np.max(np.abs(radius0-radius_mean))
    borders = 10**0
    ax.set_ylim([-borders,+borders])


def doPlotting(ax, crossSection0, velocity0, pressure0, dx, t):
        plotTube(ax, crossSection0, velocity0, pressure0, dx, t)
        #plotVar(ax[1], crossSection0, dx, t)
        plt.title(t)
        plt.pause(0.1)
        #ax[1].cla()
