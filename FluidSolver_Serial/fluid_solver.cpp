#include "fluid_nl.h"
#include "precice/SolverInterface.hpp"
#include <iostream>
#include <stdlib.h>
#include "../utils/AD.hpp"
#include "Eigen/Dense"
#include <assert.h>


using std::cout;
using std::endl;

using namespace AD;

using namespace precice;
using namespace precice::constants;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> matrixRM;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> vector;

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
  vector velocity                      = vector::Zero(N+1);
  vector velocity_n                    = vector::Zero(N+1);
  vector pressure                      = vector::Zero(N+1);
  vector pressure_n                    = vector::Zero(N+1);
  vector crossSectionLength            = vector::Zero(N+1);
  vector crossSectionLength_n          = vector::Zero(N+1);
  vector p_save                        = vector::Zero(N+1);
  vector cSL_save                      = vector::Zero(N+1);
  vector v_save                        = vector::Zero(N+1);
  
  vector positions                     = vector::Zero(N+1);
  
  matrix J_F = matrix::Zero(N+1, N+1);
  
  dualReal *v, *v_prev, *p, *p_prev, *cSL, *cSL_prev; 
  
  int dimensions = interface.getDimensions();
  
  v         = new dualReal[N+1];
  v_prev    = new dualReal[N+1];
  p         = new dualReal[N+1];
  p_prev    = new dualReal[N+1];
  cSL       = new dualReal[N+1];
  cSL_prev  = new dualReal[N+1];
  
  
  double dt = 0.01; // solver timestep size
  double precice_dt; // maximum precice timestep size

  //precice stuff
  int meshID               = interface.getMeshID("Fluid_Nodes");
  int crossSectionLengthID = interface.getDataID("CrossSectionLength", meshID);
  int pressureID           = interface.getDataID("Pressure", meshID);
  
  Eigen::VectorXi vertexIDs = Eigen::VectorXi::Zero(N+1);
  double *grid;
  grid = new double[dimensions * (N + 1)];

  for (i = 0; i < pressure.size(); i++) {
    velocity(i)             = 1.0 / (kappa * 1.0);
    velocity_n(i)           = 1.0 / (kappa * 1.0);
    crossSectionLength(i)   = 1.0;
    crossSectionLength_n(i) = 1.0;
    pressure(i)             = 0.0;
    pressure_n(i)           = 0.0;
    
    for (int dim = 0; dim < dimensions; dim++)
      grid[i * dimensions + dim] = i * (1 - dim);
    // init positions
    positions(i) = double(i);
  }

  int tstep_counter = 0; // number of time steps (only coupling iteration time steps)
  int t = 0;             // number of time steps (including subcycling time steps)
  int tsub = 0;          // number of current subcycling time steps
  int n_subcycles = 0;   // number of subcycles
  int t_steps_total = 0; // number of total timesteps, i.e., t_steps*n_subcycles

  interface.setMeshVertices(meshID, N + 1, grid, vertexIDs.data());

  cout << "Fluid: init precice..." << endl;
  precice_dt = interface.initialize();
  
  n_subcycles = (int)(precice_dt/dt);
  assert(n_subcycles == 1);
  
  t_steps_total = 100*n_subcycles;
 
  if (interface.isActionRequired(actionWriteInitialData())) {
    interface.writeBlockScalarData(pressureID, N + 1, vertexIDs.data(), pressure.data());
    //interface.initializeData();
    interface.fulfilledAction(actionWriteInitialData());
  }

  interface.initializeData();

  if (interface.isReadDataAvailable()) {
    interface.readBlockScalarData(crossSectionLengthID, N + 1, vertexIDs.data(), crossSectionLength.data());
  }
  
  while (interface.isCouplingOngoing()) {
    // When an implicit coupling scheme is used, checkpointing is required
    if (interface.isActionRequired(actionWriteIterationCheckpoint())) {
      
      if(tstep_counter > 0){
        cout << "Fluid: Advancing in time, finished timestep: " << tstep_counter << endl;
        t += n_subcycles;        
        tsub = 0;

        // store state variables from last time step (required in fluid_nl)
	velocity_n = velocity;
	pressure_n = pressure;
	crossSectionLength_n = crossSectionLength;        
      }
      tstep_counter++;
      
      interface.fulfilledAction(actionWriteIterationCheckpoint());
    }
    
    // compute adaptive time step (not needed)
    // choose smalles time step (sub-cycling if dt is smaller than precice_dt)
    dt = std::min(precice_dt, dt);
    
    // advance in time for subcycling
    tsub++;
    
    p_save = pressure;
    cSL_save = crossSectionLength;
    v_save = velocity;
    
    // build Jacobian J_F of Fluid solver, N+1 function evaluations
    for(int r = 0; r < crossSectionLength.size(); r++)
    {
      // get new dual datatypes
      for (i = 0; i <= N; i++) {
	double ei = (i==r) ? 1. : 0.;
       cSL[i] = dualReal(cSL_save[i], ei); 
       cSL_prev[i] = dualReal(crossSectionLength_n[i]);
       p[i] = dualReal(p_save[i]);
       p_prev[i] = dualReal(pressure_n[i]);
       v[i] = dualReal(v_save[i]);
       v_prev[i] = dualReal(velocity_n[i]);
      }
      
      // p_old is not used for gamma = 0.0
      fluid_nl(cSL, cSL_prev,                               // crossSectionLength
	       v, v_prev,                                   // velocity
	       p, p_prev, p,                                // pressure
	       (t + tsub)/(double)t_steps_total,            // scaled time for inflow condition (sample sine curve)
	       N, dualReal(kappa), dualReal(tau), 0.0);     // dimensionless parameters
      
      // write back
      for (i = 0; i <= N; i++){
       crossSectionLength[i] = cSL[i].u;
       crossSectionLength_n[i] = cSL_prev[i].u;
       pressure[i] = p[i].u;
       pressure_n[i] = p_prev[i].u;
       velocity[i] = v[i].u;
       velocity_n[i] = v_prev[i].u;
      }
      
      
      // \frac{ \partial F }{ \partial xr } is rth column of Jacobian
      for (i = 0; i <= N; i++) 
        J_F(r,i) = p[i].v;
    }
    vector ptil = pressure_n + J_F*(crossSectionLength - crossSectionLength_n);
    vector diff = pressure - ptil;
    Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", " << ", ";");
    
    //std::cout<<"J_F: "<<J_F<<std::endl;
    std::cout<<"\n norm(p-ptil): "<<diff.norm()<<"\n p := "<<pressure.format(CommaInitFmt)<<"\n ptil:= "<<ptil.format(CommaInitFmt)<<std::endl;
    
    // write pressure data to precice
    interface.writeBlockScalarData(pressureID, N + 1, vertexIDs.data(), pressure.data());
    
    // advance
    precice_dt = interface.advance(dt);
    
    // read crossSectionLength data from precice
    interface.readBlockScalarData(crossSectionLengthID, N + 1, vertexIDs.data(), crossSectionLength.data());

    // set variables back to checkpoint
    if (interface.isActionRequired(actionReadIterationCheckpoint())) { // i.e. not yet converged      
      std::cout<<" # ITERATE # "<<std::endl;
      tsub = 0;
 
      // also reset current state variables, but keep crossSectionLength
      // (if no subcycling is enabled, coupling iterations are slightly worse if we reset
      //  velocity and pressure. For the case (tau, kappa) = (0.01, 10) we get 6.15 its
      //  without reset and 6.20 iterations if velocity and pressure is reset)
      velocity = velocity_n; 
      pressure = pressure_n;
    
      interface.fulfilledAction(actionReadIterationCheckpoint());
    }
  }

  interface.finalize();
  cout << "Exiting FluidSolver" << endl;

  return 0;
}
