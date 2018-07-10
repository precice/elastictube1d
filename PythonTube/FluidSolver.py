from __future__ import division
import os
import sys
import argparse
import configuration_file as config
from thetaScheme import perform_partitioned_implicit_trapezoidal_rule_step, perform_partitioned_implicit_euler_step
import numpy as np
import datetime

import matplotlib.pyplot as plt
import matplotlib.animation as manimation

from output import create_output, create_video
from input import load_precomputed

def linear_interpolation(x0, x1, t0, dt):
    return lambda t: x1 * (t-t0)/dt + x0 * (t0+dt-t)/dt

# check if PRECICE_ROOT is defined
if not os.getenv('PRECICE_ROOT'):
   print "ERROR: PRECICE_ROOT not defined!"
   exit(1)

precice_root = os.getenv('PRECICE_ROOT')
precice_python_adapter_root = precice_root+"/src/precice/adapters/python"
sys.path.insert(0, precice_python_adapter_root)

from mpi4py import MPI
import PySolverInterface
from PySolverInterface import *

parser = argparse.ArgumentParser()
parser.add_argument("configurationFileName", help="Name of the xml precice configuration file.", type=str)

try:
    args = parser.parse_args()
except SystemExit:
    print("")
    print("Did you forget adding the precice configuration file as an argument?")
    print("Try $python FluidSolver.py precice-config.xml")
    quit()

print "Starting Fluid Solver..."

configFileName = args.configurationFileName

n_elem = config.n_elem
dx = config.L / n_elem  # element length

print "N: " + str(n_elem)

solverName = "FLUID"

print "Configure preCICE..."
interface = PySolverInterface(solverName, 0, 1)
interface.configure(configFileName)
print "preCICE configured..."

dimensions = interface.getDimensions()

if config.initialization_procedure is config.InitializationProcedure.FromConstants:
    crossSectionLength = config.a0 * np.ones(n_elem + 1)
    velocity = config.velocity_in(0) * crossSectionLength[0] * np.ones(n_elem + 1) / crossSectionLength  # initialize such that mass conservation is fulfilled
    pressure = config.p0 * np.ones(n_elem + 1)
elif config.initialization_procedure is config.InitializationProcedure.FromPrecomputed:
    velocity_of_t, pressure_of_t, crossSectionLength_of_t = load_precomputed(config.precomputed_filename)
    velocity = velocity_of_t(0)
    pressure = pressure_of_t(0)
    crossSectionLength = crossSectionLength_of_t(0)
else:
    print "invalid initialization procedure!"
    quit()

plotting_mode = config.PlottingModes.OFF
output_mode = config.OutputModes.NETCDF

if plotting_mode == config.PlottingModes.VIDEO:
    fig, ax = plt.subplots(1)
    FFMpegWriter = manimation.writers['imagemagick']
    writer = FFMpegWriter(fps=15)
    writer.setup(fig, "writer_test.mp4", 100)

meshID = interface.getMeshID("Fluid_Nodes")
crossSectionLengthID = interface.getDataID("CrossSectionLength", meshID)
pressureID = interface.getDataID("Pressure", meshID)

vertexIDs = np.zeros(n_elem + 1)
grid = np.zeros([dimensions, n_elem + 1])

grid[0,:] = np.linspace(0, config.L, n_elem + 1)  # x component
grid[1,:] = 0  # y component, leave blank

interface.setMeshVertices(meshID, n_elem + 1, grid.flatten('F'), vertexIDs)

t = 0

print "Fluid: init precice..."
precice_tau = interface.initialize()

if interface.isActionRequired(PyActionWriteInitialData()):
    interface.writeBlockScalarData(pressureID, n_elem + 1, vertexIDs, pressure)
    interface.fulfilledAction(PyActionWriteInitialData())

interface.initializeData()

if interface.isReadDataAvailable():
    interface.readBlockScalarData(crossSectionLengthID, n_elem + 1, vertexIDs, crossSectionLength)

velocity_n = np.copy(velocity)
pressure_n = np.copy(pressure)
crossSectionLength_n = np.copy(crossSectionLength)

sim_start_time = datetime.datetime.now()

while interface.isCouplingOngoing():
    # When an implicit coupling scheme is used, checkpointing is required
    if interface.isActionRequired(PyActionWriteIterationCheckpoint()):
        interface.fulfilledAction(PyActionWriteIterationCheckpoint())

    if config.time_stepping_scheme is config.TimeStepping.ImplicitEuler:
        if config.coupling_mode is config.CouplingAlgorithm.PartitionedPreCICE:
            velocity, pressure, success = perform_partitioned_implicit_euler_step(velocity_n, pressure_n, crossSectionLength_n, crossSectionLength, dx, precice_tau, config.velocity_in(t + precice_tau))
        elif config.coupling_mode is config.CouplingAlgorithm.PartitionedPreCICECustomized:
            crossSectionLength_of_t = linear_interpolation(crossSectionLength_n, crossSectionLength, t, precice_tau)
            velocity, pressure, success = perform_partitioned_implicit_euler_step(velocity_n, pressure_n, crossSectionLength_of_t(t), crossSectionLength_of_t(t+precice_tau), dx, precice_tau, config.velocity_in(t + precice_tau))
        else:
            raise Exception("invalid combination of time stepping scheme [%s] and coupling mode [%s]!" % (config.time_stepping_scheme.name, config.coupling_mode.name))
    elif config.time_stepping_scheme is config.TimeStepping.TrapezoidalRule:
        if config.coupling_mode is config.CouplingAlgorithm.PartitionedPreCICE:
            velocity, pressure, success = perform_partitioned_implicit_trapezoidal_rule_step(velocity_n, pressure_n, crossSectionLength_n, crossSectionLength, dx, precice_tau, config.velocity_in(t + precice_tau), custom_coupling=False)
        elif config.coupling_mode is config.CouplingAlgorithm.PartitionedPreCICECustomized:
            crossSectionLength_of_t = linear_interpolation(crossSectionLength_n, crossSectionLength, t, precice_tau)
            velocity, pressure, success = perform_partitioned_implicit_trapezoidal_rule_step(velocity_n, pressure_n, crossSectionLength_of_t(t), crossSectionLength_of_t(t+precice_tau), dx, precice_tau, config.velocity_in(t + precice_tau), custom_coupling=True)
        else:
            raise Exception("invalid combination of time stepping scheme [%s] and coupling mode [%s]!" % (config.time_stepping_scheme.name, config.coupling_mode.name))
    else:
        raise Exception("invalid time stepping scheme [%s]!" % (config.time_stepping_scheme.name))

    interface.writeBlockScalarData(pressureID, n_elem + 1, vertexIDs, pressure)
    interface.advance(precice_tau)
    interface.readBlockScalarData(crossSectionLengthID, n_elem + 1, vertexIDs, crossSectionLength)

    if interface.isActionRequired(PyActionReadIterationCheckpoint()): # i.e. not yet converged
        interface.fulfilledAction(PyActionReadIterationCheckpoint())
    else: # converged, timestep complete
        t += precice_tau
        velocity_n = np.copy(velocity)
        pressure_n = np.copy(pressure)
        crossSectionLength_n = np.copy(crossSectionLength)

        metadata = {
            'created_on': str(sim_start_time),
            'tau': precice_tau,
            'dx': dx,
            'timestepping': config.time_stepping_scheme.name,
            'coupling': config.coupling_mode.name,
            'elasticity_module': config.E,
            'length': config.L,
            'n_elem': config.n_elem,
            'inflow_frequency': config.frequency,
            'inflow_amplitude': config.ampl,
            'inflow_mean': config.u0
        }

        create_output(t, velocity, pressure, crossSectionLength, metadata, output_mode)

        if plotting_mode is config.PlottingModes.VIDEO:
            create_video(t, velocity, pressure, crossSectionLength, metadata, writer)


print "Exiting FluidSolver"

if plotting_mode is config.PlottingModes.VIDEO:
    writer.finish()

interface.finalize()
