#include "fluid_nl.h"
#include "precice/SolverInterface.hpp"
#include <iostream>
#include <stdlib.h>

using std::cout;
using std::endl;

using namespace precice;
using namespace precice::constants;

void printData(const std::vector<double>& data)
{
  std::cout << "Received data = " << data[0];
  for (size_t i = 1; i < data.size(); i++) {
    std::cout << ", " << data[i];
  }
  std::cout << std::endl;
}

int main(int argc, char** argv)
{
  cout << "Starting Fluid Solver..." << endl;

  if (argc != 5) {
    cout << endl;
    cout << "Usage: " << argv[0] << " configurationFileName N tau kappa" << endl;
    cout << endl;
    cout << "N:     Number of mesh elements, needs to be equal for fluid and structure solver." << endl;
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
  double *velocity_subcycle_n, *pressure_subcycle_n, *crossSectionLength_subcycle_n;
  
  int dimensions = interface.getDimensions();

  velocity             = new double[N + 1]; // Speed
  velocity_n           = new double[N + 1];
  velocity_subcycle_n  = new double[N + 1];
  pressure             = new double[N + 1]; // Pressure
  pressure_n           = new double[N + 1];
  pressure_subcycle_n  = new double[N + 1];
  crossSectionLength   = new double[N + 1];
  crossSectionLength_n = new double[N + 1];
  crossSectionLength_subcycle_n = new double[N + 1];
  
  double dt = 0.01; // solver timestep size
  double precice_dt; // maximum precice timestep size

  //precice stuff
  int meshID               = interface.getMeshID("Fluid_Nodes");
  int crossSectionLengthID = interface.getDataID("CrossSectionLength", meshID);
  int pressureID           = interface.getDataID("Pressure", meshID);
  
  int* vertexIDs;
  double* grid;
  vertexIDs = new int[(N + 1)];
  grid = new double[dimensions * (N + 1)];

  for (i = 0; i <= N; i++) {
    velocity[i]             = 1.0 / (kappa * 1.0);
    velocity_n[i]           = 1.0 / (kappa * 1.0);
    velocity_subcycle_n[i]  = 1.0 / (kappa * 1.0);
    crossSectionLength[i]   = 1.0;
    crossSectionLength_n[i] = 1.0;
    crossSectionLength_subcycle_n[i] = 1.0;
    pressure[i]             = 0.0;
    pressure_n[i]           = 0.0;
    pressure_subcycle_n[i]  = 0.0;
    
    for (int dim = 0; dim < dimensions; dim++)
      grid[i * dimensions + dim] = i * (1 - dim);
  }

  int tstep_counter = 0; // number of time steps (only coupling iteration time steps)
  int t = 0;             // number of time steps (including subcycling time steps)
  int tsub = 0;          // number of current subcycling time steps
  int n_subcycles = 0;   // number of subcycles
  int t_steps_total = 0; // number of total timesteps, i.e., t_steps*n_subcycles

  interface.setMeshVertices(meshID, N + 1, grid, vertexIDs);

  cout << "Fluid: init precice..." << endl;
  precice_dt = interface.initialize();
  
  n_subcycles = (int)(precice_dt/dt);
  t_steps_total = 100*n_subcycles;
  
  std::cout<<"n_subcycles: "<<n_subcycles<<" t_steps_total: "<<t_steps_total<<std::endl;

  if (interface.isActionRequired(actionWriteInitialData())) {
    interface.writeBlockScalarData(pressureID, N + 1, vertexIDs, pressure);
    //interface.initializeData();
    interface.fulfilledAction(actionWriteInitialData());
  }

  interface.initializeData();

  if (interface.isReadDataAvailable()) {
    interface.readBlockScalarData(crossSectionLengthID, N + 1, vertexIDs, crossSectionLength);
  }
  
  while (interface.isCouplingOngoing()) {
    // When an implicit coupling scheme is used, checkpointing is required
    if (interface.isActionRequired(actionWriteIterationCheckpoint())) {
      
      // save old state, save checkpoint (not needed here)
      
      interface.fulfilledAction(actionWriteIterationCheckpoint());
    }
    
    // compute adaptive time step (not needed)
    // choose smalles time step (sub-cycling if dt is smaller than precice_dt)
    dt = std::min(precice_dt, dt);
    
    // advance in time for subcycling
    tsub++;

    // p_old is not used for gamma = 0.0
    fluid_nl(crossSectionLength, crossSectionLength_subcycle_n,  // crossSectionLength
	     velocity, velocity_subcycle_n,                      // velocity
	     pressure, pressure_subcycle_n, pressure,            // pressure
	     (t + tsub)/(double)t_steps_total,                   // scaled time for inflow condition (sample sine curve)
	     N, kappa, tau, 0.0);                                // dimensionless parameters
    
    std::cout<<"scaled t: "<<(t + tsub)/(double)t_steps_total<<std::endl;
    
    // store state variables for previous subcycle
    for (i = 0; i <= N; i++) {
      velocity_subcycle_n[i] = velocity[i];
      pressure_subcycle_n[i] = pressure[i];
      crossSectionLength_subcycle_n[i] = crossSectionLength[i];
    }

    // write pressure data to precice
    interface.writeBlockScalarData(pressureID, N + 1, vertexIDs, pressure);
    
    // advance
    precice_dt = interface.advance(dt);
    
    // read crossSectionLength data from precice
    interface.readBlockScalarData(crossSectionLengthID, N + 1, vertexIDs, crossSectionLength);

    // set variables back to checkpoint
    if (interface.isActionRequired(actionReadIterationCheckpoint())) { // i.e. not yet converged      
      std::cout<<" # ITERATE # "<<std::endl;
      tsub = 0;
      for (i = 0; i <= N; i++) {
        velocity_subcycle_n[i] = velocity_n[i];
        pressure_subcycle_n[i] = pressure_n[i];
        crossSectionLength_subcycle_n[i] = crossSectionLength_n[i];
	
	// also reset current state variables, but keep crossSectionLength
	// (if no subcycling is enabled, coupling iterations are slightly worse if we reset
	//  velocity and pressure. For the case (tau, kappa) = (0.01, 10) we get 6.15 its
	//  without reset and 6.20 iterations if velocity and pressure is reset)
	velocity[i] = velocity_n[i]; 
	pressure[i] = pressure[i];
      }
      
      interface.fulfilledAction(actionReadIterationCheckpoint());
    } else {
      cout << "Fluid: Advancing in time, finished timestep: " << tstep_counter << endl;
      t += n_subcycles;
      tstep_counter++;
      tsub = 0;

      // store state variables from last time step (required in fluid_nl)
      for (i = 0; i <= N; i++) {
        velocity_n[i]           = velocity[i];
        pressure_n[i]           = pressure[i];
        crossSectionLength_n[i] = crossSectionLength[i];
	
	velocity_subcycle_n[i]           = velocity[i];
        pressure_subcycle_n[i]           = pressure[i];
        crossSectionLength_subcycle_n[i] = crossSectionLength[i];
      }
    }
  }

  interface.finalize();
  cout << "Exiting FluidSolver" << endl;

  return 0;
}
