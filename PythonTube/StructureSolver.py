from __future__ import division

import os
import sys
import argparse
from mpi4py import MPI
import numpy as np
import configuration_file as config
 
import PySolverInterface
from PySolverInterface import *

print "Starting Structure Solver..."

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

print "N: " + str(N)

solverName = "STRUCTURE"

print "Configure preCICE..."
interface = PySolverInterface(solverName, 0, 1)
interface.configure(configFileName)
print "preCICE configured..."

dimensions = interface.getDimensions()

pressure = config.p0 * np.ones(N+1)
crossSectionLength = config.a0 * np.ones(N+1)

meshID = interface.getMeshID("Structure_Nodes")
crossSectionLengthID = interface.getDataID("CrossSectionLength", meshID)
pressureID = interface.getDataID("Pressure", meshID)

vertexIDs = np.zeros(N+1)
grid = np.zeros([dimensions, N+1])

grid[0,:] = np.linspace(0, config.L, N+1)  # x component
#grid[1,:] = np.linspace(0, config.L, N+1)  # y component, leave blank

interface.setMeshVertices(meshID, N+1, grid.flatten('F'), vertexIDs)

t = 0

print "Structure: init precice..."
precice_tau = interface.initialize()

if (interface.isActionRequired(PyActionWriteInitialData())):
   interface.writeBlockScalarData(crossSectionLengthID, N+1, vertexIDs, crossSectionLength)
   interface.fulfilledAction(PyActionWriteInitialData())

interface.initializeData()

if (interface.isReadDataAvailable()):
   interface.readBlockScalarData(pressureID, N+1, vertexIDs, pressure)

crossSection0 = config.crossSection0(pressure.shape[0] - 1)
pressure0 = config.p0 * np.ones_like(pressure)

while interface.isCouplingOngoing():
   # When an implicit coupling scheme is used, checkpointing is required
    if interface.isActionRequired(PyActionWriteIterationCheckpoint()):
        interface.fulfilledAction(PyActionWriteIterationCheckpoint())

    crossSectionLength = crossSection0 * ((pressure0 - 2.0 * config.c_mk ** 2) ** 2 / (pressure - 2.0 * config.c_mk ** 2) ** 2)

    interface.writeBlockScalarData(crossSectionLengthID, N + 1, vertexIDs, crossSectionLength)
    interface.advance(precice_tau)
    interface.readBlockScalarData(pressureID, N + 1, vertexIDs, pressure)

    if interface.isActionRequired(PyActionReadIterationCheckpoint()): # i.e. not yet converged
        interface.fulfilledAction(PyActionReadIterationCheckpoint())
    else:
        t += precice_tau

print "Exiting StructureSolver"

interface.finalize()

