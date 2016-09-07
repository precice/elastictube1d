import os
import sys
import argparse
from mpi4py import MPI
 
# check if PRECICE_ROOT is defined
if not os.getenv('PRECICE_ROOT'):
   print "ERROR: PRECICE_ROOT not defined!"
   exit(1)

precice_root = os.getenv('PRECICE_ROOT')
precice_python_adapter_root = precice_root+"/src/precice/adapters/python"
sys.path.insert(0, precice_python_adapter_root)

import PySolverInterface
from PySolverInterface import *

print "Starting Structure Solver..."

parser = argparse.ArgumentParser()
parser.add_argument("configurationFileName", help="Name of the xml config file.", type=str)
parser.add_argument("N", help="Number of mesh elements, needs to be equal for fluid and structure solver.", type=int)
args = parser.parse_args()

configFileName = args.configurationFileName
N              = args.N

print "N: " + str(N)

solverName = "STRUCTURE";

print "Configure preCICE..."
interface = PySolverInterface(solverName, 0, 1)
interface.configure(configFileName)
print "preCICE configured..."

dimensions = interface.getDimensions();

pressure               = [0.0001] * (N+1) #Pressure
crossSectionLength     = [1.0] * (N+1) #Cross-section length

meshID = interface.getMeshID("Structure_Nodes")
crossSectionLengthID = interface.getDataID("CrossSectionLength", meshID)
pressureID = interface.getDataID("Pressure", meshID)

vertexIDs = [0.0]*(N+1)
grid = [None]*(dimensions*(N+1))

for i in range(0, N+1):
   for dim in range(0, dimensions):
      grid[i*dimensions + dim]= i*(1-dim)

t = 0

interface.setMeshVertices(meshID, N+1, grid, vertexIDs)

print "Structure: init precice..."
interface.initialize()

if (interface.isActionRequired(PyActionWriteInitialData())):
   interface.writeBlockScalarData(crossSectionLengthID, N+1, vertexIDs, crossSectionLength)
   interface.fulfilledAction(PyActionWriteInitialData())

interface.initializeData();

if (interface.isReadDataAvailable()):
   interface.readBlockScalarData(pressureID, N+1, vertexIDs, pressure)

while (interface.isCouplingOngoing()):
   # When an implicit coupling scheme is used, checkpointing is required
   if (interface.isActionRequired(PyActionWriteIterationCheckpoint())):
      interface.fulfilledAction(PyActionWriteIterationCheckpoint())

   for i in range(0, N+1):
      crossSectionLength[i] = 4.0 / ((2.0 - pressure[i]) * (2.0 - pressure[i]))

   interface.writeBlockScalarData(pressureID, N+1, vertexIDs, pressure)
   interface.advance(0.01)
   interface.readBlockScalarData(crossSectionLengthID, N+1, vertexIDs, crossSectionLength)

   if (interface.isActionRequired(PyActionReadIterationCheckpoint())): # i.e. not yet converged
      interface.fulfilledAction(PyActionReadIterationCheckpoint())
   else:
      t = t + 1

print "Exiting StructureSolver"

interface.finalize()

