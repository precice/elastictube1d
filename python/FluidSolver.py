from __future__ import division, print_function
import os
import sys
import argparse
import configuration_file as config
from thetaScheme import perform_partitioned_implicit_trapezoidal_rule_step, perform_partitioned_implicit_euler_step
import numpy as np
import tubePlotting

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as manimation

from output import writeOutputToVTK

import precice
from precice import *

parser = argparse.ArgumentParser()
parser.add_argument("configurationFileName", help="Name of the xml precice configuration file.", nargs='?', type=str, default="precice-config.xml")
parser.add_argument("--write-vtk", help="Save vtk files of each timestep in the 'VTK' folder.", action='store_true')
parser.add_argument("--enable-plot", help="Show a continuously updated plot of the tube while simulating.", action='store_true')
parser.add_argument("--write-video", help="Save a video of the simulation as 'writer_test.mp4'. \
                    NOTE: This requires 'enable_plot' to be active!", action='store_true')

try:
    args = parser.parse_args()
except SystemExit:
    print("")
    print("Did you forget adding the precice configuration file as an argument?")
    print("Try '$ python FluidSolver.py precice-config.xml'")
    quit()

output_mode = config.OutputModes.VTK if args.write_vtk else config.OutputModes.OFF
plotting_mode = config.PlottingModes.VIDEO if args.enable_plot else config.PlottingModes.OFF
if args.write_video and not args.enable_plot:
    print("")
    print("To create a video it is required to enable plotting for this run.")
    print("Please supply both the '--enable-plot' and '--write-video' flags.")
    quit()
writeVideoToFile = True if args.write_video else False

print("Starting Fluid Solver...")

configFileName = args.configurationFileName

N = config.n_elem
dx = config.L / N  # element length

print("N: " + str(N))

solverName = "FLUID"

print("Configure preCICE...")
interface = precice.Interface(solverName, configFileName, 0, 1)
print("preCICE configured...")

dimensions = interface.get_dimensions()

velocity = config.velocity_in(0) * np.ones(N+1)
velocity_n = config.velocity_in(0) * np.ones(N+1)
pressure = config.p0 * np.ones(N+1)
pressure_n = config.p0 * np.ones(N+1)
crossSectionLength = config.a0 * np.ones(N+1)
crossSectionLength_n = config.a0 * np.ones(N+1)


if plotting_mode == config.PlottingModes.VIDEO:
    fig, ax = plt.subplots(1)
    if writeVideoToFile:
        FFMpegWriter = manimation.writers['imagemagick']
        metadata = dict(title='PulseTube')
        writer = FFMpegWriter(fps=15, metadata=metadata)
        writer.setup(fig, "writer_test.mp4", 100)

meshID = interface.get_mesh_id("Fluid_Nodes")
crossSectionLengthID = interface.get_data_id("CrossSectionLength", meshID)
pressureID = interface.get_data_id("Pressure", meshID)

vertexIDs = np.zeros(N+1)
grid = np.zeros([N+1, dimensions])

grid[:, 0] = np.linspace(0, config.L, N+1)  # x component
grid[:, 1] = 0  # y component, leave blank

vertexIDs = interface.set_mesh_vertices(meshID, grid)

t = 0

print("Fluid: init precice...")
precice_tau = interface.initialize()

if interface.is_action_required(action_write_initial_data()):
    interface.write_block_scalar_data(pressureID, vertexIDs, pressure)
    interface.fulfilled_action(action_write_initial_data())

interface.initialize_data()

if interface.is_read_data_available():
    interface.read_block_scalar_data(crossSectionLengthID, vertexIDs, crossSectionLength)

crossSectionLength_n = np.copy(crossSectionLength)
velocity_n = config.velocity_in(0) * crossSectionLength_n[0] * np.ones(N+1) / crossSectionLength_n  # initialize such that mass conservation is fulfilled

print(crossSectionLength_n)

while interface.is_coupling_ongoing():
    # When an implicit coupling scheme is used, checkpointing is required
    if interface.is_action_required(action_write_iteration_checkpoint()):
        interface.fulfilled_action(action_write_iteration_checkpoint())

    velocity, pressure, success = perform_partitioned_implicit_euler_step(velocity_n, pressure_n, crossSectionLength_n, crossSectionLength, dx, precice_tau, config.velocity_in(t + precice_tau), custom_coupling=False)
    interface.write_block_scalar_data(pressureID, vertexIDs, pressure)
    precice_tau = interface.advance(precice_tau)
    interface.read_block_scalar_data(crossSectionLengthID, vertexIDs, crossSectionLength)

    if interface.is_action_required(action_read_iteration_checkpoint()): # i.e. not yet converged
        interface.fulfilled_action(action_read_iteration_checkpoint())
    else:  # converged, timestep complete
        t += precice_tau
        if plotting_mode is config.PlottingModes.VIDEO:
            tubePlotting.doPlotting(ax, crossSectionLength_n, velocity_n, pressure_n, dx, t)
            if writeVideoToFile:
                writer.grab_frame()
            ax.cla()
        velocity_n = np.copy(velocity)
        pressure_n = np.copy(pressure)
        crossSectionLength_n = np.copy(crossSectionLength)
        if output_mode is config.OutputModes.VTK:
            writeOutputToVTK(t, "fluid", dx, N+1, datanames=["velocity", "pressure", "crossSection"], data=[velocity_n, pressure_n, crossSectionLength_n])

print("Exiting FluidSolver")

if plotting_mode is config.PlottingModes.VIDEO and writeVideoToFile:
    writer.finish()

interface.finalize()
