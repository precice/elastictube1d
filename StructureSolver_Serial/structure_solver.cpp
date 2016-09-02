#include "precice/SolverInterface.hpp"
#include <iostream>
#include <stdlib.h>
//#include "mpi.h"

using std::cout;
using std::endl;

void printData(const std::vector<double>& data)
{
  cout << "Received data = " << data[0];
  for (size_t i = 1; i < data.size(); i++) {
    cout << ", " << data[i];
  }
  cout << endl;
}

int main(int argc, char** argv)
{
  cout << "Starting Structure Solver..." << endl;
  using namespace precice;
  using namespace precice::constants;

  if (argc != 3) {
    cout << endl;
    cout << "Usage: " << argv[0] << " configurationFileName N" << endl;
    cout << endl;
    cout << "N:     Number of mesh elements, needs to be equal for fluid and structure solver." << endl;
    return -1;
  }

  std::string configFileName(argv[1]);
  int N = atoi(argv[2]);

  std::cout << "N: " << N << std::endl;

  std::string dummyName = "STRUCTURE";

  SolverInterface interface(dummyName, 0, 1);
  interface.configure(configFileName);
  cout << "preCICE configured..." << endl;

  //init data
  double *crossSectionLength, *pressure;
  int dimensions = interface.getDimensions();
  crossSectionLength = new double[N + 1]; // Second dimension (only one cell deep) stored right after the first dimension: see SolverInterfaceImpl::setMeshVertices
  pressure = new double[N + 1];
  double* grid;
  grid = new double[dimensions * (N + 1)];
  
  double dt = 0.02; // solver timestep size
  double precice_dt; // maximum precice timestep size

  //precice stuff
  int meshID = interface.getMeshID("Structure_Nodes");
  int crossSectionLengthID = interface.getDataID("CrossSectionLength", meshID);
  int pressureID = interface.getDataID("Pressure", meshID);
  int* vertexIDs;
  vertexIDs = new int[N + 1];

  for (int i = 0; i <= N; i++) {
    crossSectionLength[i] = 1.0;
    pressure[i] = 0.0;
    for (int dim = 0; dim < dimensions; dim++)
      grid[i * dimensions + dim] = i * (1 - dim); // Define the y-component of each grid point as zero
  }

  int tstep_counter = 0; // number of time steps (only coupling iteration time steps)
  int t = 0;             // number of time steps (including subcycling time steps)
  int tsub = 0;          // number of current subcycling time steps
  int n_subcycles = 0;   // number of subcycles
  int t_steps_total = 0; // number of total timesteps, i.e., t_steps*n_subcycles
  
  interface.setMeshVertices(meshID, N + 1, grid, vertexIDs);

  cout << "Structure: init precice..." << endl;
  precice_dt = interface.initialize();
  
  n_subcycles = (int)(precice_dt/dt);
  t_steps_total = 100*n_subcycles;

  if (interface.isActionRequired(actionWriteInitialData())) {
    interface.writeBlockScalarData(crossSectionLengthID, N + 1, vertexIDs, crossSectionLength);
    //interface.initializeData();
    interface.fulfilledAction(actionWriteInitialData());
  }

  interface.initializeData();

  if (interface.isReadDataAvailable()) {
    interface.readBlockScalarData(pressureID, N + 1, vertexIDs, pressure);
  }

  while (interface.isCouplingOngoing()) {
    // When an implicit coupling scheme is used, checkpointing is required
    if (interface.isActionRequired(actionWriteIterationCheckpoint())) {
      
      if(tstep_counter > 0){
        cout << "Advancing in time, finished timestep: " << tstep_counter << endl;
        t += n_subcycles;        
        tsub = 0;
      }
      tstep_counter++;
      
      // write checkpoint, save state variables (not needed here, stationary solver)       
      interface.fulfilledAction(actionWriteIterationCheckpoint());
    }

    // choose smalles time step (sub-cycling if dt is smaller than precice_dt)
    dt = std::min(precice_dt, dt);
    
    // advance in time for subcycling
    tsub++;
    
    for (int i = 0; i <= N; i++) {
      crossSectionLength[i] = 4.0 / ((2.0 - pressure[i]) * (2.0 - pressure[i]));
    }

    // send crossSectionLength data to precice
    interface.writeBlockScalarData(crossSectionLengthID, N + 1, vertexIDs, crossSectionLength);
    
    // advance
    precice_dt = interface.advance(dt);
    
    // receive pressure data from precice
    interface.readBlockScalarData(pressureID, N + 1, vertexIDs, pressure);

    if (interface.isActionRequired(actionReadIterationCheckpoint())) {
      cout << "Iterate" << endl;
      tsub = 0;
      
      interface.fulfilledAction(actionReadIterationCheckpoint());
    }
  }

  interface.finalize();
  cout << "Exiting StructureSolver" << endl;

  return 0;
}
