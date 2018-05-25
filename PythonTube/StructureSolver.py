from __future__ import division

import os
import sys
import argparse
from mpi4py import MPI
import numpy as np
import configuration_file as config
 
# check if PRECICE_ROOT is defined
if not os.getenv('PRECICE_ROOT'):
   print "ERROR: PRECICE_ROOT not defined!"
   exit(1)

precice_root = os.getenv('PRECICE_ROOT')
precice_python_adapter_root = precice_root+"/src/precice/adapters/python"
sys.path.insert(0, precice_python_adapter_root)

import PySolverInterface
from PySolverInterface import *
from pythonCouplingHelpers.solid import solve_solid

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

if config.initialization_procedure is config.InitializationProcedure.FromConstants:
    pressure = config.p0 * np.ones(N+1)
    crossSectionLength = config.a0 * np.ones(N+1)
elif config.initialization_procedure is config.InitializationProcedure.FromPrecomputed:
    print "has to be implemented!"
    # todo to be implemented!
    # pressure = ...
    # crossSectionLength = ...
    quit()
else:
    print "invalid initialization procedure!"
    quit()

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

while interface.isCouplingOngoing():
   # When an implicit coupling scheme is used, checkpointing is required
    if interface.isActionRequired(PyActionWriteIterationCheckpoint()):
        interface.fulfilledAction(PyActionWriteIterationCheckpoint())

    crossSectionLength = solve_solid(pressure)

    interface.writeBlockScalarData(crossSectionLengthID, N + 1, vertexIDs, crossSectionLength)
    interface.advance(precice_tau)
    interface.readBlockScalarData(pressureID, N + 1, vertexIDs, pressure)

    if interface.isActionRequired(PyActionReadIterationCheckpoint()): # i.e. not yet converged
        interface.fulfilledAction(PyActionReadIterationCheckpoint())
    else:
        t += precice_tau

print "Exiting StructureSolver"

interface.finalize()

