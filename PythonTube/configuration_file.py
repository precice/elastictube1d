from __future__ import division
from enum import Enum

import numpy as np
import sympy as sp


class PlottingModes(Enum):
    OFF = 0  # no plotting over time
    VIDEO = 1  # create a video
    DEBUG = 2  # provide a debug plot over time of the simulation


# physical properties of the tube
r0 = 1/np.sqrt(np.pi)  # radius of the tube
a0 = r0**2 * np.pi  # cross sectional area
E = 10000  # elasticity module
L = 10  # length of tube/simulation domain

# physical properties of the fluid
# inflow velocity
u0 = 10  # mean velocity
ampl = 3  # amplitude of varying velocity
frequency = 10  # frequency of variation
t_shift = 0  # temporal shift of variation

p0 = 0  # pressure at outlet

# physical properties of coupled problem
c_mk = np.sqrt(E/2/r0)  # wave speed

velocity_in = lambda t: u0 + ampl * np.sin(frequency * (t+t_shift) * np.pi)  # inflow velocity

# helper function to create constant cross section
def crossSection0(N):
    return a0 * np.ones(N + 1)

# numerical properties
k_max_nonlin = 1000  # maximum number of non-linear solver iterations

# discretization
tau = 10**10  # timestep size, set it to a large value to enforce tau from precice_config.xml
n_elem = 100  # number of elements in x direction
