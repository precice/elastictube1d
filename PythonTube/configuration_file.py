from __future__ import division
from enum import Enum

import numpy as np
import sympy as sp


class PlottingModes(Enum):
    OFF = 0  # no plotting over time
    VIDEO = 1  # create a video
    DEBUG = 2  # provide a debug plot over time of the simulation


class OutputModes(Enum):
    OFF = 0  # no plotting over time
    VTK = 1  # produce VTK output
    NETCDF = 2  # produce NETCDF output


class TimeStepping(Enum):
    ImplicitEuler = 0  # implicit Euler time stepping
    TrapezoidalRule = 1  # trapezoidal rule time stepping
    TrapezoidalRuleCustom = 2  # trapezoidal rule time stepping with pseudo waveform relaxation


class InitializationProcedure(Enum):
    FromConstants = 0  # use given constants for initialization
    FromPrecomputed = 1  # use precomputed solution for initialization


# physical properties of the tube
r0 = 1/np.sqrt(np.pi)  # radius of the tube
a0 = r0**2 * np.pi  # cross sectional area
E = 100000  # elasticity module
L = .1  # length of tube/simulation domain

# physical properties of the fluid
# inflow velocity
u0 = 10  # mean velocity
ampl = 3  # amplitude of varying velocity
frequency = 1  # frequency of variation
t_shift = 0  # temporal shift of variation

p0 = 0  # pressure at outlet

initialization_procedure = InitializationProcedure.FromConstants

# physical properties of coupled problem
c_mk = np.sqrt(E/2/r0)  # wave speed

velocity_in = lambda t: u0 + ampl * np.sin(frequency * (t+t_shift) * np.pi)**2  # inflow velocity


# helper function to create constant cross section
def crossSection0(N):
    return a0 * np.ones(N + 1)


# numerical properties
k_max_nonlin = 1000  # maximum number of non-linear solver iterations
tol_nonlin = 10e-12  # tolerance for non-linear solver

# discretization
tau = 10**10  # timestep size, set it to a large value to enforce tau from precice_config.xml
n_elem = 10  # number of elements in x direction
time_stepping_scheme = TimeStepping.ImplicitEuler  # time stepping scheme used
