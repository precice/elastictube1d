#include <iostream>
#include <stdlib.h>
#include "precice/SolverInterface.hpp"
#include "fluid_nl.h"

using std::cout;
using std::endl;

using namespace precice;
using namespace precice::constants;

void printData (const std::vector<double>& data)
{
  std::cout << "Received data = " << data[0];
  for (size_t i=1; i < data.size(); i++){
    std::cout << ", " << data[i];
  }
  std::cout << std::endl;
}


int main (int argc, char **argv)
{
  cout << "Starting Fluid Solver..." << endl;

  if ( argc != 5 )
  {
    cout << endl;
    cout << "Usage: " << argv[0] <<  " configurationFileName N tau kappa" << endl;
    cout << endl;
    cout << "N:     Number of mesh elements, needs to be equal for fluid and structure solver." << endl;
    cout << "tau:   Dimensionless time step size." << endl;
    cout << "kappa: Dimensionless structural stiffness." << endl;
    return -1;
  }

  std::string configFileName(argv[1]);
  int N        = atoi( argv[2] );
  double tau   = atof( argv[3] );
  double kappa = atof( argv[4] );

  std::cout << "N: " << N << " tau: " << tau << " kappa: " << kappa << std::endl;

  std::string solverName = "FLUID";

  cout << "Configure preCICE..." << endl;
  // Initialize the solver interface with our name, our process index (like rank) and the total number of processes.
  SolverInterface interface(solverName, 0, 1);
  // Provide the configuration file to precice. After configuration a usuable state of that SolverInterface is reached.
  // Reads the XML file and contacts the server, if used.
  interface.configure(configFileName);
  cout << "preCICE configured..." << endl;

  // init data
  int i;
  double *velocity, *velocity_n, *pressure, *pressure_n, *crossSectionLength, *crossSectionLength_n;
  int dimensions = interface.getDimensions();

  velocity     = new double[N+1]; // Speed
  velocity_n   = new double[N+1];
  pressure     = new double[N+1]; // Pressure
  pressure_n   = new double[N+1];
  crossSectionLength     = new double[N+1];
  crossSectionLength_n   = new double[N+1];

  //precice stuff
  int meshID = interface.getMeshID("Fluid_Nodes");
  int crossSectionLengthID = interface.getDataID("Displacements", meshID);
  int pressureID = interface.getDataID("Stresses", meshID);
  int *vertexIDs;
  vertexIDs = new int[(N+1)];
  double *grid;
  grid = new double[dimensions*(N+1)];

  for ( i = 0; i <= N; i++ )
  {
    velocity[i]   = 1.0 / ( kappa * 1.0 );
    velocity_n[i] = 1.0 / ( kappa * 1.0 );
    crossSectionLength[i]   = 1.0;
    crossSectionLength_n[i] = 1.0;
    pressure[i]   = 0.0;
    pressure_n[i] = 0.0;
    for(int dim = 0; dim < dimensions; dim++)
      grid[i*dimensions + dim]= i*(1-dim);
  }

  int t = 0; //number of timesteps

  interface.setMeshVertices(meshID, N+1, grid, vertexIDs);

  cout << "Fluid: init precice..." << endl;
  interface.initialize();


  if (interface.isActionRequired(actionWriteInitialData())) {
    interface.writeBlockScalarData(pressureID,N+1,vertexIDs,pressure);
    //interface.initializeData();
    interface.fulfilledAction(actionWriteInitialData());
  }

  interface.initializeData();

  if (interface.isReadDataAvailable()) {
    interface.readBlockScalarData(crossSectionLengthID,N+1,vertexIDs,crossSectionLength);
  }

  while (interface.isCouplingOngoing()){
    // When an implicit coupling scheme is used, checkpointing is required
    if (interface.isActionRequired(actionWriteIterationCheckpoint())){
      interface.fulfilledAction(actionWriteIterationCheckpoint());
    }

    //p_old is not used for gamma = 0.0
    fluid_nl( crossSectionLength, crossSectionLength_n, velocity, velocity_n, pressure, pressure_n, pressure, t+1, N, kappa, tau, 0.0);

    interface.writeBlockScalarData(pressureID, N+1, vertexIDs, pressure);
    interface.advance(0.01);
    interface.readBlockScalarData(crossSectionLengthID, N+1, vertexIDs, crossSectionLength);

    if (interface.isActionRequired(actionReadIterationCheckpoint())){ // i.e. not yet converged
      // cout << "Iterate" << endl;
      interface.fulfilledAction(actionReadIterationCheckpoint());
    }
    else {
      // cout << "Fluid: Advancing in time, finished timestep: " << t << endl;
      t++;

      for ( i = 0; i <=N; i++)
      {
	velocity_n[i] = velocity[i];
	pressure_n[i] = pressure[i];
	crossSectionLength_n[i] = crossSectionLength[i];
      }
    }
  }

  interface.finalize();
  cout << "Exiting FluidSolver" << endl;

  return 0;
}
