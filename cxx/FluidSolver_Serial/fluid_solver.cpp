#include "fluid_nl.h"
#include "precice/SolverInterface.hpp"
#include <iostream>
#include <stdlib.h>

using std::cout;
using std::endl;

using namespace precice;
using namespace precice::constants;

int main(int argc, char** argv)
{
  cout << "Starting Fluid Solver..." << endl;

  if (argc != 5) {
    cout << endl;
    cout << "Usage: " << argv[0] << " configurationFileName N tau kappa" << endl;
    cout << endl;
    cout << "N:     Number of mesh elements." << endl;
    cout << "tau:   Dimensionless time step size." << endl;
    cout << "kappa: Dimensionless structural stiffness." << endl;
    return -1;
  }

  std::string configFileName(argv[1]);
  int N = atoi(argv[2]);
  double tau = atof(argv[3]);
  double kappa = atof(argv[4]);

  std::cout << "N: " << N << " tau: " << tau << " kappa: " << kappa << std::endl;

  std::string solverName = "FLUID";
  
  std::string outputFilePrefix = "Postproc/out_fluid";

  cout << "Configure preCICE..." << endl;
  // Create preCICE with the solver's name, the rank, and the total number of processes.
  SolverInterface interface(solverName, configFileName, 0, 1);

  int i;
  double *velocity, *velocity_n, *pressure, *pressure_n, *crossSectionLength, *crossSectionLength_n;
  
  int dimensions = interface.getDimensions();

  velocity             = new double[N + 1];
  velocity_n           = new double[N + 1];
  pressure             = new double[N + 1]; 
  pressure_n           = new double[N + 1];
  crossSectionLength   = new double[N + 1];
  crossSectionLength_n = new double[N + 1];

  // get IDs from preCICE
  int meshID               = interface.getMeshID("Fluid_Nodes");
  int crossSectionLengthID = interface.getDataID("CrossSectionLength", meshID);
  int pressureID           = interface.getDataID("Pressure", meshID);
  
  int* vertexIDs;
  double* grid;
  vertexIDs = new int[(N + 1)];
  grid = new double[dimensions * (N + 1)];

  // init data values and mesh
  for (i = 0; i <= N; i++) {
    velocity[i]             = 1.0 / (kappa * 1.0);
    velocity_n[i]           = 1.0 / (kappa * 1.0);
    crossSectionLength[i]   = 1.0;
    crossSectionLength_n[i] = 1.0;
    pressure[i]             = 0.0;
    pressure_n[i]           = 0.0;
    
    for (int dim = 0; dim < dimensions; dim++)
      grid[i * dimensions + dim] = i * (1 - dim);
  }

  double t = 0.0; // time
  double dt = 0.01; // solver timestep size

  // tell preCICE about your coupling interface mesh
  interface.setMeshVertices(meshID, N + 1, grid, vertexIDs);

  cout << "Initialize preCICE..." << endl;
  interface.initialize();
  
  // write initial data if required
  if (interface.isActionRequired(actionWriteInitialData())) {
    interface.writeBlockScalarData(pressureID, N + 1, vertexIDs, pressure);
    interface.fulfilledAction(actionWriteInitialData());
  }

  // initial data is sent or received if necessary
  interface.initializeData();

  // read data if available
  if (interface.isReadDataAvailable()) {
    interface.readBlockScalarData(crossSectionLengthID, N + 1, vertexIDs, crossSectionLength);
  }
  int out_counter = 0;  
  
  while (interface.isCouplingOngoing()) {
    // for an implicit coupling, you can store an iteration checkpoint here (from the first iteration of a timestep)
    // this is, however, not necessary for this scenario
    if (interface.isActionRequired(actionWriteIterationCheckpoint())) {
      interface.fulfilledAction(actionWriteIterationCheckpoint());
    }
    
    fluid_nl(crossSectionLength, crossSectionLength_n,  
	     velocity, velocity_n,                      
	     pressure, pressure_n,            
	     t, N, kappa, tau); 
    
    // write pressure data to precice
    interface.writeBlockScalarData(pressureID, N + 1, vertexIDs, pressure);
    
    interface.advance(dt);
    
    // read crossSectionLength data from precice
    interface.readBlockScalarData(crossSectionLengthID, N + 1, vertexIDs, crossSectionLength);

    // set variables back to checkpoint
    if (interface.isActionRequired(actionReadIterationCheckpoint())) { 
    // i.e. not yet converged, you could restore a checkpoint here (not necessary for this scenario)      
      interface.fulfilledAction(actionReadIterationCheckpoint());
    }
    else{
      t += dt;
      for (i = 0; i <= N; i++) {
        velocity_n[i]           = velocity[i];
        pressure_n[i]           = pressure[i];
        crossSectionLength_n[i] = crossSectionLength[i];
      }      
      write_vtk(t, out_counter, outputFilePrefix.c_str(), N, grid, velocity_n, pressure_n, crossSectionLength_n);
      out_counter++;
    }
  }

  interface.finalize();

  delete [] velocity;
  delete [] velocity_n;
  delete [] pressure;
  delete [] pressure_n;
  delete [] crossSectionLength;
  delete [] crossSectionLength_n;
  delete [] vertexIDs;
  delete [] grid;

  return 0;
}
