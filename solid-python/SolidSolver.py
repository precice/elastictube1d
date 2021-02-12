from __future__ import division, print_function

import os
import sys
import argparse
import numpy as np
import configuration_file as config

import precice
from precice import *

print("Starting Solid Solver...")

parser = argparse.ArgumentParser()
parser.add_argument("configurationFileName", help="Name of the xml config file.", nargs='?', type=str,
                    default="precice-config.xml")

try:
    args = parser.parse_args()
except SystemExit:
    print("")
    print("Did you forget adding the precice configuration file as an argument?")
    print("Try '$ python SolidSolver.py precice-config.xml'")
    quit()

configFileName = args.configurationFileName
N = config.n_elem

print("N: " + str(N))

solverName = "Solid"

print("Configure preCICE...")
interface = precice.Interface(solverName, configFileName, 0, 1)
print("preCICE configured...")

dimensions = interface.get_dimensions()

pressure = config.p0 * np.ones(N + 1)
crossSectionLength = config.a0 * np.ones(N + 1)

meshID = interface.get_mesh_id("Solid-Nodes")
crossSectionLengthID = interface.get_data_id("CrossSectionLength", meshID)
pressureID = interface.get_data_id("Pressure", meshID)

vertexIDs = np.zeros(N + 1)
grid = np.zeros([N + 1, dimensions])

grid[:, 0] = np.linspace(0, config.L, N + 1)  # x component
grid[:, 1] = 0  # np.linspace(0, config.L, N+1)  # y component, leave blank

vertexIDs = interface.set_mesh_vertices(meshID, grid)

t = 0

print("Solid: init precice...")

# preCICE defines timestep size of solver via precice-config.xml
precice_dt = interface.initialize()

if interface.is_action_required(action_write_initial_data()):
    interface.write_block_scalar_data(crossSectionLengthID, vertexIDs, crossSectionLength)
    interface.mark_action_fulfilled(action_write_initial_data())

interface.initialize_data()

if interface.is_read_data_available():
    pressure = interface.read_block_scalar_data(pressureID, vertexIDs)

crossSection0 = config.crossSection0(pressure.shape[0] - 1)
pressure0 = config.p0 * np.ones_like(pressure)

while interface.is_coupling_ongoing():
    # When an implicit coupling scheme is used, checkpointing is required
    if interface.is_action_required(action_write_iteration_checkpoint()):
        interface.mark_action_fulfilled(action_write_iteration_checkpoint())

    crossSectionLength = crossSection0 * (
                (pressure0 - 2.0 * config.c_mk ** 2) ** 2 / (pressure - 2.0 * config.c_mk ** 2) ** 2)

    interface.write_block_scalar_data(crossSectionLengthID, vertexIDs, crossSectionLength)
    precice_dt = interface.advance(precice_dt)
    pressure = interface.read_block_scalar_data(pressureID, vertexIDs)

    if interface.is_action_required(action_read_iteration_checkpoint()):  # i.e. not yet converged
        interface.mark_action_fulfilled(action_read_iteration_checkpoint())
    else:
        t += precice_dt

print("Exiting SolidSolver")

interface.finalize()
