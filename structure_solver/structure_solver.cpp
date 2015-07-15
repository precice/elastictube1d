#include <iostream>
#include <stdlib.h>
#include "precice/SolverInterface.hpp"
//#include "mpi.h"

using std::cout;
using std::endl;

void printData (const std::vector<double>& data)
{
  cout << "Received data = " << data[0];
  for (size_t i=1; i < data.size(); i++){
    cout << ", " << data[i];
  }
  cout << endl;
}


int main (int argc, char **argv)
{
  cout << "Starting Structure Solver..." << endl;
  using namespace precice;
  using namespace precice::constants;

  if ( argc != 3 )
  {
    cout << endl;
    cout << "Usage: " << argv[0] << " configurationFileName N" << endl;
    cout << endl;
    cout << "N:     Number of mesh elements, needs to be equal for fluid and structure solver." << endl;
    return -1;
  }

  std::string configFileName(argv[1]);
  int N = atoi( argv[2] );

  std::cout << "N: " << N << std::endl;

  std::string dummyName = "STRUCTURE";

  SolverInterface interface(dummyName, 0, 1);
  interface.configure(configFileName);
  cout << "preCICE configured..." << endl;

  //init data
  double *displ, *sigma;
  int dimensions = interface.getDimensions();
  displ     = new double[N+1];  // Second dimension (only one cell deep) stored right after the first dimension: see SolverInterfaceImpl::setMeshVertices
  sigma     = new double[N+1];
  double *grid;
  grid = new double[dimensions*(N+1)];

  //precice stuff
  int meshID = interface.getMeshID("Structure_Nodes");
  int displID = interface.getDataID ( "Displacements", meshID );
  int sigmaID = interface.getDataID ( "Stresses", meshID );
  int *vertexIDs;
  vertexIDs = new int[N+1];

  for(int i=0; i<=N; i++)
  {
    displ[i] = 1.0;
    sigma[i] = 0.0;
    for(int dim = 0; dim < dimensions; dim++)
      grid[i*dimensions + dim] = 0;//i*(1-dim);   // Define the y-component of each grid point as zero
  }

  int t = 0;
  interface.setMeshVertices(meshID, N+1, grid, vertexIDs);
  //for (int i=0; i < (N+1); i++){
    //cout << "VertexID: " << vertexIDs[i] << " | grid: " << grid[i*dimensions + 0] << ", " << grid[i*dimensions + 1] << " | Displacements: " << displ[i*dimensions + 0] << ", " << displ[i*dimensions + 1] << " | Stresses: " << sigma[i*dimensions + 0] << ", " << sigma[i*dimensions + 1] << "\n";
  //}
  cout << "\n";
  //for(int i=0;i<=N;i++)
  //{
  //  vertexIDs[i] = interface.setMeshVertex(meshID, static_cast<const double*>(grid + i));
  //}

  cout << "Structure: init precice..." << endl;
  double dt = interface.initialize();

  if (interface.isActionRequired(actionWriteInitialData()))
  {
    interface.writeBlockScalarData(displID, N+1, vertexIDs, displ);
    //interface.initializeData();
    interface.fulfilledAction(actionWriteInitialData());
  }

  interface.initializeData();


  if (interface.isReadDataAvailable())
  {
    interface.readBlockScalarData(sigmaID , N+1, vertexIDs, sigma);
  }

  while (interface.isCouplingOngoing())
  {
    // When an implicit coupling scheme is used, checkpointing is required
    if (interface.isActionRequired(actionWriteIterationCheckpoint()))
    {
      interface.fulfilledAction(actionWriteIterationCheckpoint());
    }

    for ( int i = 0; i <= N; i++ )
    {
      displ[i]   = 4.0 / ((2.0 - sigma[i])*(2.0 - sigma[i]));
    }

    interface.writeBlockScalarData(displID, N+1, vertexIDs, displ);
    interface.advance(0.01); // Advance by dt = 0.01
    interface.readBlockScalarData(sigmaID, N+1, vertexIDs, sigma);

    if (interface.isActionRequired(actionReadIterationCheckpoint()))
    {
      cout << "Iterate" << endl;
      interface.fulfilledAction(actionReadIterationCheckpoint());
    }
    else
    {
      cout << "Advancing in time, finished timestep: " << t << endl;
      t++;
    }
  }

  interface.finalize();
  cout << "Exiting StructureSolver" << endl;

  return 0;
}
