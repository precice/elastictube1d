#include "StructureSolver.h"
#include "precice/SolverInterface.hpp"

#include <iostream>
#include <mpi.h>

using namespace precice;
using namespace precice::constants;

int main(int argc, char** argv)
{
  std::cout << "Starting Structure Solver..." << std::endl;
  if (argc != 3) {
    std::cout << std::endl;
    std::cout << "Structure: Usage: mpiexec -np <#procs> " << argv[0] << " <configurationFileName> <N>" << std::endl;
    std::cout << std::endl;
    std::cout << "N:     Number of mesh elements, needs to be equal for fluid and structure solver." << std::endl;
    return -1;
  }

  MPI_Init(&argc, &argv);

  int domainSize, gridOffset, rank, size, chunkLength;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  domainSize = atoi(argv[2]);
  if ((domainSize + 1) % size == 0) {
    chunkLength = (domainSize + 1) / size;
    gridOffset = rank * chunkLength;
  } else if (rank < (domainSize + 1) % size) {
    chunkLength = (domainSize + 1) / size + 1;
    gridOffset = rank * chunkLength;
  } else {
    chunkLength = (domainSize + 1) / size;
    gridOffset = ((domainSize + 1) % size) * ((domainSize + 1) / size + 1) + (rank - ((domainSize + 1) % size)) * (domainSize + 1) / size;
  }

  std::vector<double> pressure(chunkLength);
  std::vector<double> crossSectionLength(chunkLength);

  std::string configFileName(argv[1]);
  std::string solverName = "STRUCTURE";

  SolverInterface interface(solverName, configFileName, rank, size);

  int meshID = interface.getMeshID("Structure_Nodes");
  int pressureID = interface.getDataID("Pressure", meshID);
  int crossSectionLengthID = interface.getDataID("CrossSectionLength", meshID);

  int dimensions = interface.getDimensions();
  std::vector<double> grid(dimensions * chunkLength);

  std::vector<int> vertexIDs(chunkLength);

  for (int i = 0; i < chunkLength; i++) {
    crossSectionLength[i] = 1.0;
    pressure[i] = 0.0;
    for (int j = 0; j < dimensions; j++) {
      grid[i * dimensions + j] = j == 0 ? gridOffset + (double)i : 0.0;
    }
  }

  interface.setMeshVertices(meshID, chunkLength, grid.data(), vertexIDs.data());

  interface.initialize();

  double t = 0;
  double dt = 0.01;

  if (interface.isActionRequired(actionWriteInitialData())) {
    interface.writeBlockScalarData(crossSectionLengthID, chunkLength, vertexIDs.data(), crossSectionLength.data());
    interface.fulfilledAction(actionWriteInitialData());
  }

  interface.initializeData();

  if (interface.isReadDataAvailable()) {
    interface.readBlockScalarData(pressureID, chunkLength, vertexIDs.data(), pressure.data());
  }

  while (interface.isCouplingOngoing()) {
    if (interface.isActionRequired(actionWriteIterationCheckpoint())) {
      interface.fulfilledAction(actionWriteIterationCheckpoint());
    }

    structureComputeSolution(rank, size, chunkLength, pressure.data(), crossSectionLength.data()); // Call Solver
                                                                                                   //structureDataDisplay(crossSectionLength, chunkLength);
                                                                                                   //structureDataDisplay(pressure, chunkLength);

    interface.writeBlockScalarData(crossSectionLengthID, chunkLength, vertexIDs.data(), crossSectionLength.data());

    interface.advance(dt);

    interface.readBlockScalarData(pressureID, chunkLength, vertexIDs.data(), pressure.data());

    if (interface.isActionRequired(actionReadIterationCheckpoint())) { // i.e. fluid not yet converged
      interface.fulfilledAction(actionReadIterationCheckpoint());
    } else {
      t += dt;
    }
  }

  interface.finalize();
  MPI_Finalize();

  return 0;
}
