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

  std::string dummyName = "STRUCTURE_1D";
  
  SolverInterface interface(dummyName, 0, 1);
  interface.configure(configFileName);

  //init data
  double *p, *a;
  p     = new double[N+1];
  a     = new double[N+1];

  //precice stuff
  int meshID = interface.getMeshID("WetSurface");
  int aID = interface.getDataID ( "CrossSectionalArea", meshID );
  int pID = interface.getDataID ( "Pressure", meshID );
  int *vertexIDs;
  vertexIDs = new int[N+1];

  for(int i=0; i<=N; i++)
  {
    vertexIDs[i]=i;
    a[i]=1.0;
    p[i]=0.0;
  }
 
  int t = 0;

  cout << "Structure: init precice..." << endl;
  double dt = interface.initialize();

  if (interface.isActionRequired(actionWriteInitialData()))
  {
    interface.writeBlockScalarData(aID, N+1, vertexIDs, a);
    //interface.initializeData();
    interface.fulfilledAction(actionWriteInitialData());
  }

  interface.initializeData();


  if (interface.isReadDataAvailable())
  {
    interface.readBlockScalarData(pID , N+1, vertexIDs, p);
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
      a[i]   = 4.0 / ((2.0 - p[i])*(2.0 - p[i]));
    }
	
    interface.writeBlockScalarData(aID, N+1, vertexIDs, a);
    interface.advance(0.01); // Advance by dt = 0.01
    interface.readBlockScalarData(pID, N+1, vertexIDs, p);
      
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
