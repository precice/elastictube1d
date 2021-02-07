#include "FluidSolver.h"
#include "fluid_nl.h"

#include "precice/SolverInterface.hpp"
#include <iostream>
#include <mpi.h>

using namespace precice;
using namespace precice::constants;

int main(int argc, char** argv)
{

  std::cout << "Starting Fluid Solver..." << std::endl;
  if (argc != 5 && argc != 6) {
    std::cout << std::endl;
    std::cout << "Fluid: Usage: mpiexec -np <#procs> " << argv[0] << " <configurationFileName> <N> <tau> <kappa> -parallel" << std::endl;
    std::cout << "or" << std::endl;
    std::cout << "Usage: " << argv[0] << " configurationFileName> <N> <tau> <kappa>" << std::endl;
    std::cout << std::endl;
    std::cout << "N:     Number of mesh elements, needs to be equal for fluid and Solid solver." << std::endl;
    std::cout << "tau:   Dimensionless time step size." << std::endl;
    std::cout << "kappa: Dimensionless structural stiffness." << std::endl;

    return -1;
  }

  std::string configFileName(argv[1]);
  int domainSize = atoi(argv[2]); //N
  int chunkLength = domainSize + 1; //serial run
  double tau = atof(argv[3]);
  double kappa = atof(argv[4]);

  std::string solverName = "FLUID";

  std::string outputFilePrefix = "Postproc/out_fluid"; //extra

    
  int gridOffset, rank = 0, size = 1;


  if (argc == 6){
	  MPI_Init(&argc, &argv);
  	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  	MPI_Comm_size(MPI_COMM_WORLD, &size);
    outputFilePrefix += std::to_string(rank); //extra

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
  }



  SolverInterface interface(solverName, configFileName, rank , size);
  double *velocity, *pressure,  *crossSectionLength;

  velocity             = new double[chunkLength];
  pressure             = new double[chunkLength]; 
  crossSectionLength   = new double[chunkLength];

  int *vertexIDs;
  vertexIDs = new int[chunkLength];

  int dimensions = interface.getDimensions();
  int meshID = interface.getMeshID("Fluid_Nodes");
  int pressureID = interface.getDataID("Pressure", meshID);
  int crossSectionLengthID = interface.getDataID("CrossSectionLength", meshID);

  
  double* grid;
  grid = new double[dimensions * chunkLength];
  
  for (int i = 0; i < chunkLength; i++) {
    pressure[i] = 0.0;
    crossSectionLength[i] = 1.0;
    velocity[i] = 1.0 / kappa;
  }

  if (argc==6){
    for (int i = 0; i < chunkLength; i++) {
      for (int j = 0; j < dimensions; j++) {
          grid[i * dimensions + j] = j == 0 ? gridOffset + (double)i : 0.0;
      }
    
    }
  } else {
    for (int i = 0; i < chunkLength; i++) {
      for (int j = 0; j < dimensions; j++) {
          grid[i * dimensions + j] = i * (1 - j);
      }
    }
  }


  interface.setMeshVertices(meshID, chunkLength, grid, vertexIDs);

  std::cout << "Initialize preCICE..." << std::endl;
  interface.initialize();

  double t = 0.0;
  double dt = 0.01;

  if (interface.isActionRequired(actionWriteInitialData())) {
    interface.writeBlockScalarData(pressureID, chunkLength, vertexIDs, pressure);
    interface.markActionFulfilled(actionWriteInitialData());
  }

  interface.initializeData();

  if (interface.isReadDataAvailable()) {
    interface.readBlockScalarData(crossSectionLengthID, chunkLength, vertexIDs, crossSectionLength);
  }

  int out_counter = 0;

  while (interface.isCouplingOngoing()) {
    int convergenceCounter = 0;
    if (interface.isActionRequired(actionWriteIterationCheckpoint())) {
      interface.markActionFulfilled(actionWriteIterationCheckpoint());
    }

    if (argc == 6){
	    fluidComputeSolution(rank, size, domainSize, chunkLength, kappa, tau, 0.0, t+dt,
      pressure, 
      crossSectionLength, 
      velocity);
    } else {
      fluid_nl(crossSectionLength,   
	    velocity,           
	    pressure,     
	    t, domainSize, kappa, tau); 
    }
    

    interface.writeBlockScalarData(pressureID, chunkLength, vertexIDs, pressure);

    interface.advance(dt);

    interface.readBlockScalarData(crossSectionLengthID, chunkLength, vertexIDs, crossSectionLength);

    if (interface.isActionRequired(actionReadIterationCheckpoint())) { // i.e. not yet converged
      interface.markActionFulfilled(actionReadIterationCheckpoint());
      convergenceCounter++;
    } else {
      t += dt;
      write_vtk(t, out_counter, outputFilePrefix.c_str(), chunkLength, grid, velocity, pressure, crossSectionLength);

      out_counter++;
    }
  }

  delete [] velocity;
  delete [] pressure;
  delete [] crossSectionLength;
  delete [] vertexIDs;
  delete [] grid;
  interface.finalize();
  if (argc == 6){
    MPI_Finalize();

  }

  return 0;
}