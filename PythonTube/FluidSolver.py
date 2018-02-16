from __future__ import division
import os
import sys
import argparse
import configuration_file as config
from thetaScheme import perform_partitioned_implicit_trapezoidal_rule_step, perform_partitioned_implicit_euler_step
import numpy as np
import tubePlotting

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as manimation

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

print "Starting Fluid Solver..."

parser = argparse.ArgumentParser()
parser.add_argument("configurationFileName", help="Name of the xml config file.", type=str)
args = parser.parse_args()

configFileName = args.configurationFileName
N = config.n_elem
dx = config.L / N  # element length

print "N: " + str(N)

solverName = "FLUID"

print "Configure preCICE..."
interface = PySolverInterface(solverName, 0, 1)
interface.configure(configFileName)
print "preCICE configured..."

dimensions = interface.getDimensions()

velocity = config.velocity_in(0) * np.ones(N+1)
velocity_n = config.velocity_in(0) * np.ones(N+1)
pressure = config.p0 * np.ones(N+1)
pressure_n = config.p0 * np.ones(N+1)
crossSectionLength = config.a0 * np.ones(N+1)
crossSectionLength_n = config.a0 * np.ones(N+1)

plotting_mode = config.PlottingModes.VIDEO

if plotting_mode == config.PlottingModes.VIDEO:
    fig, ax = plt.subplots(1)
    FFMpegWriter = manimation.writers['ffmpeg']
    writer = FFMpegWriter(fps=15)
    writer.setup(fig, "writer_test.mp4", 100)

meshID = interface.getMeshID("Fluid_Nodes")
crossSectionLengthID = interface.getDataID("CrossSectionLength", meshID)
pressureID = interface.getDataID("Pressure", meshID)

vertexIDs = np.zeros(N+1)
grid = np.zeros([dimensions, N+1])

grid[0,:] = np.linspace(0, config.L, N+1)  # x component
grid[1,:] = 0  # y component, leave blank

interface.setMeshVertices(meshID, N+1, grid.flatten('F'), vertexIDs)

t = 0

print "Fluid: init precice..."
precice_tau = interface.initialize()

if interface.isActionRequired(PyActionWriteInitialData()):
    interface.writeBlockScalarData(pressureID, N+1, vertexIDs, pressure)
    interface.fulfilledAction(PyActionWriteInitialData())

interface.initializeData()

if interface.isReadDataAvailable():
    interface.readBlockScalarData(crossSectionLengthID, N+1, vertexIDs, crossSectionLength)

crossSectionLength_n = np.copy(crossSectionLength)
velocity_n = config.velocity_in(0) * crossSectionLength_n[0] * np.ones(N+1) / crossSectionLength_n  # initialize such that mass conservation is fulfilled

print crossSectionLength_n

while interface.isCouplingOngoing():
    # When an implicit coupling scheme is used, checkpointing is required
    if interface.isActionRequired(PyActionWriteIterationCheckpoint()):
        interface.fulfilledAction(PyActionWriteIterationCheckpoint())

    velocity, pressure, success = perform_partitioned_implicit_euler_step(velocity_n, pressure_n, crossSectionLength_n, crossSectionLength, dx, precice_tau, config.velocity_in(t + precice_tau), custom_coupling=False)
    interface.writeBlockScalarData(pressureID, N+1, vertexIDs, pressure)
    interface.advance(precice_tau)
    interface.readBlockScalarData(crossSectionLengthID, N+1, vertexIDs, crossSectionLength)

    if interface.isActionRequired(PyActionReadIterationCheckpoint()): # i.e. not yet converged
        interface.fulfilledAction(PyActionReadIterationCheckpoint())
    else: # converged, timestep complete
        t += precice_tau
        if plotting_mode is config.PlottingModes.VIDEO:
            tubePlotting.doPlotting(ax, crossSectionLength_n, velocity_n, pressure_n, dx, t)
            writer.grab_frame()
            ax.cla()
        velocity_n = np.copy(velocity)
        pressure_n = np.copy(pressure)
        crossSectionLength_n = np.copy(crossSectionLength)

print "Exiting FluidSolver"

interface.finalize()

