from __future__ import division, print_function

import os
import sys
import argparse
from mpi4py import MPI
import numpy as np
import configuration_file as config
 
import precice
from precice import *

print("Starting Structure Solver...")

parser = argparse.ArgumentParser()
parser.add_argument("configurationFileName", help="Name of the xml config file.", type=str)

try:
    args = parser.parse_args()
except SystemExit:
    print("")
    print("Did you forget adding the precice configuration file as an argument?")
    print("Try $python StructureSolver.py precice-config.xml")
    quit()

configFileName = args.configurationFileName
N = config.n_elem

print("N: " + str(N))

solverName = "STRUCTURE"

print("Configure preCICE...")
interface = precice.Interface(solverName, 0, 1)
interface.configure(configFileName)
print("preCICE configured...")

dimensions = interface.get_dimensions()

pressure = config.p0 * np.ones(N+1)
crossSectionLength = config.a0 * np.ones(N+1)

meshID = interface.get_mesh_id("Structure_Nodes")
crossSectionLengthID = interface.get_data_id("CrossSectionLength", meshID)
pressureID = interface.get_data_id("Pressure", meshID)

vertexIDs = np.zeros(N+1)
grid = np.zeros([dimensions, N+1])

grid[0,:] = np.linspace(0, config.L, N+1)  # x component
#grid[1,:] = np.linspace(0, config.L, N+1)  # y component, leave blank

interface.set_mesh_vertices(meshID, N+1, grid.flatten('F'), vertexIDs)

t = 0

print("Structure: init precice...")
precice_tau = interface.initialize()

if (interface.is_action_required(action_write_initial_data())):
   interface.write_block_scalar_data(crossSectionLengthID, N+1, vertexIDs, crossSectionLength)
   interface.fulfilled_action(action_write_initial_data())

interface.initialize_data()

if (interface.is_read_data_available()):
   interface.read_block_scalar_data(pressureID, N+1, vertexIDs, pressure)

crossSection0 = config.crossSection0(pressure.shape[0] - 1)
pressure0 = config.p0 * np.ones_like(pressure)

while interface.is_coupling_ongoing():
   # When an implicit coupling scheme is used, checkpointing is required
    if interface.is_action_required(action_write_iteration_checkpoint()):
        interface.fulfilled_action(action_write_iteration_checkpoint())

    crossSectionLength = crossSection0 * ((pressure0 - 2.0 * config.c_mk ** 2) ** 2 / (pressure - 2.0 * config.c_mk ** 2) ** 2)

    interface.write_block_scalar_data(crossSectionLengthID, N + 1, vertexIDs, crossSectionLength)
    interface.advance(precice_tau)
    interface.read_block_scalar_data(pressureID, N + 1, vertexIDs, pressure)

    if interface.is_action_required(action_read_iteration_checkpoint()): # i.e. not yet converged
        interface.fulfilled_action(action_read_iteration_checkpoint())
    else:
        t += precice_tau

print("Exiting StructureSolver")

interface.finalize()

