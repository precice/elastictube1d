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
  double *crossSectionLength, *pressure;
  int dimensions = interface.getDimensions();
  crossSectionLength     = new double[N+1];  // Second dimension (only one cell deep) stored right after the first dimension: see SolverInterfaceImpl::setMeshVertices
  pressure     = new double[N+1];
  double *grid;
  grid = new double[dimensions*(N+1)];

  //precice stuff
  int meshID = interface.getMeshID("Structure_Nodes");
  int crossSectionLengthID = interface.getDataID ( "CrossSectionLength", meshID );
  int pressureID = interface.getDataID ( "Pressure", meshID );
  int *vertexIDs;
  vertexIDs = new int[N+1];

  for(int i=0; i<=N; i++)
  {
    crossSectionLength[i] = 1.0;
    pressure[i] = 0.0;
    for(int dim = 0; dim < dimensions; dim++)
      grid[i*dimensions + dim] = i*(1-dim);   // Define the y-component of each grid point as zero
  }

  int t = 0;
  interface.setMeshVertices(meshID, N+1, grid, vertexIDs);

  cout << "Structure: init precice..." << endl;
  double dt = interface.initialize();

  if (interface.isActionRequired(actionWriteInitialData()))
  {
    interface.writeBlockScalarData(crossSectionLengthID, N+1, vertexIDs, crossSectionLength);
    //interface.initializeData();
    interface.fulfilledAction(actionWriteInitialData());
  }

  interface.initializeData();


  if (interface.isReadDataAvailable())
  {
    interface.readBlockScalarData(pressureID , N+1, vertexIDs, pressure);
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
      crossSectionLength[i]   = 4.0 / ((2.0 - pressure[i])*(2.0 - pressure[i]));
    }

    interface.writeBlockScalarData(crossSectionLengthID, N+1, vertexIDs, crossSectionLength);
    interface.advance(0.01); // Advance by dt = 0.01
    interface.readBlockScalarData(pressureID, N+1, vertexIDs, pressure);

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
