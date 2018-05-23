from __future__ import division
import numpy as np
import configuration_file as conf

def solve_solid(pressure):
    """
    compute cross section area of the tube from pressure and elasticity module
    :param pressure:
    :return: new cross section
    """
    crossSection0 = conf.crossSection0(pressure.shape[0] - 1)
    pressure0 = conf.p0 * np.ones_like(pressure)
    crossSection = crossSection0 * ((pressure0 - 2.0 * conf.c_mk ** 2) ** 2 / (pressure - 2.0 * conf.c_mk ** 2) ** 2)
    return crossSection


def dsolve_solid(pressure):
    """
    compute derivative in time of cross section area of the tube from pressure and elasticity module
    :param pressure:
    :return: new derivative in time of cross section
    """
    crossSection0 = conf.crossSection0(pressure.shape[0] - 1)
    pressure0 = conf.p0 * np.ones_like(pressure)
    crossSection = -2 * crossSection0 * ((pressure0 - 2.0 * conf.c_mk ** 2) ** 2 / (pressure - 2.0 * conf.c_mk ** 2) ** 3)
    return crossSection
