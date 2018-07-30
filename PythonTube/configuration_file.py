from __future__ import division
from enum import Enum

import numpy as np


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


class InitializationProcedure(Enum):
    FromConstants = 0  # use given constants for initialization
    FromPrecomputed = 1  # use precomputed solution for initialization


class CouplingAlgorithm(Enum):
    PartitionedPythonImplicit = 0
    PartitionedPythonExplicit = 1
    PartitionedPythonCustomized = 2
    PartitionedPreCICE = 3
    PartitionedPreCICECustomized = 4
    Monolitic = 5


class CrossSectionGeometry(Enum):
    Straight = 1
    LinearlyGrowing = 2
    Sinus = 3


class VelocityInType(Enum):
    Constant = 1
    Growing = 2
    Sinus = 3
    SinusSquared = 4
    RampedSinusSquared = 5


## STRUCTURE SOLVER PARAMETERS
# physical properties of the tube
r0 = 1/np.sqrt(np.pi)  # radius of the tube
a0 = r0**2 * np.pi  # cross sectional area
E = 1000  # elasticity module
L = 1  # length of tube/simulation domain

# helper function to create constant cross section
def crossSection0(N, params={'a0':a0}, geometry=CrossSectionGeometry.Straight):
    if geometry is CrossSectionGeometry.Straight:
        return params['a0'] * np.ones(N + 1)
    elif geometry is CrossSectionGeometry.LinearlyGrowing:
        return params['a0'] * np.ones(N + 1) + np.linspace(0, params['a1']-params['a0'], N+1)
    elif geometry is CrossSectionGeometry.Sinus:
        return params['a0'] * np.ones(N + 1) + params['a0'] * np.sin(np.linspace(0,np.pi,N+1))
    else:
        print("unknown geometry type!")
        quit()

## FLUID SOLVER PARAMETERS
# physical properties of the fluid

# inflow velocity
u0 = 10  # mean velocity
ampl = 1  # amplitude of varying velocity
frequency = 1  # frequency of variation
t_shift = 0  # temporal shift of variation

default_velocity_type = VelocityInType.SinusSquared

def velocity_in(t, type=default_velocity_type):
    if type is VelocityInType.Constant:
        return u0 # constant inflow velocity
    elif type is VelocityInType.Growing:
        return u0 + t * u0
    elif type is VelocityInType.Sinus:
        return u0 + ampl * np.sin(frequency * (t+t_shift) * np.pi)  # inflow velocity
    elif type is VelocityInType.SinusSquared:
        return u0 + ampl * np.sin(frequency * (t+t_shift) * np.pi)**2  # inflow velocity
    elif type is VelocityInType.RampedSinusSquared:
        if t < 1.0/frequency:
            return u0
        else:
            return velocity_in(t, type=VelocityInType.SinusSquared)
    else:
        print("unknown type!")
        quit()


def dvelocity_in(t, type=default_velocity_type):
    if type is VelocityInType.Constant:
        return 0 # constant inflow velocity
    elif type is VelocityInType.Growing:
        return u0
    elif type is VelocityInType.Sinus:
        return ampl * np.cos(frequency * (t+t_shift) * np.pi) * frequency * np.pi
    elif type is VelocityInType.SinusSquared:
        return u0 + 2 * ampl * np.sin(frequency * (t+t_shift) * np.pi) * np.cos(frequency * (t+t_shift) * np.pi) * frequency * np.pi  # inflow velocity
    elif type is VelocityInType.RampedSinusSquared:
        if t < 1.0/frequency:
            return 0
        else:
            return dvelocity_in(t, type=VelocityInType.SinusSquared)
    else:
        print("unknown type!")
        quit()

p0 = 0  # pressure at outlet

initialization_procedure = InitializationProcedure.FromConstants
precomputed_filename = "ENTER FILENAME HERE!"

# numerical properties for non-linear fluid solver
k_max_nonlin = 1000  # maximum number of non-linear solver iterations
tol_nonlin = 10e-6  # tolerance for non-linear solver

## COUPLING and COMMON QUANTITIES
# physical properties of coupled problem
c_mk = np.sqrt(E/2/r0)  # wave speed

# coupling numerical properties (not used in preCICE)
underrelaxation_factor=.0000001  # underrelaxation factor
relax_newton = .1
k_max_coupling = 1000  # maximum number of coupling iterations per timestep
e_coupling = 10**-10  # error tolerance in coupling
coupling_mode = CouplingAlgorithm.PartitionedPythonImplicit

# discretization
T_max = 1  # total simulation time
tau0 = .5
n_elem = 10  # number of elements in x direction
time_stepping_scheme = TimeStepping.ImplicitEuler  # time stepping scheme used

# experimental setup
n_tau = 4
