#include <iostream>
#include <stdlib.h>
#include "precice/SolverInterface.hpp"
#include "fluid_nl.h"
#include "NearestNeighborMapping.hpp"

using std::cout;
using std::endl;

using namespace precice;
using namespace precice::constants;

void printData(
    const std::vector<double>& data)
{
  std::cout << "Received data = " << data[0];
  for (size_t i = 1; i < data.size(); i++) {
    std::cout << ", " << data[i];
  }
  std::cout << std::endl;
}


int main(
    int argc, char **argv)
{
  cout << "Starting Fluid Solver..." << endl;

  if (argc != 6)
      {
    cout << endl;
    cout << "Usage: " << argv[0] << " configurationFileName N tau kappa" << endl;
    cout << endl;
    cout << "N:     Number of mesh elements, needs to be equal for fluid and structure solver." << endl;
    cout << "N_SM:  Number of surrogate model mesh elements, needs to be equal for fluid and structure solver." << endl;
    cout << "tau:   Dimensionless time step size." << endl;
    cout << "kappa: Dimensionless structural stiffness." << endl;
    return -1;
  }

  std::string configFileName(argv[1]);
  int N = atoi(argv[2]);
  int N_SM = atoi(argv[3]);
  double tau = atof(argv[4]);
  double kappa = atof(argv[5]);

  std::cout << "N: " << N <<" N (surrogate model): " << N_SM << " tau: " << tau << " kappa: " << kappa << std::endl;

  std::string solverName = "FLUID";

  cout << "Configure preCICE..." << endl;
  // Initialize the solver interface with our name, our process index (like rank) and the total number of processes.
  SolverInterface interface(solverName, 0, 1);
  // Provide the configuration file to precice. After configuration a usuable state of that SolverInterface is reached.
  // Reads the XML file and contacts the server, if used.
  interface.configure(configFileName);
  cout << "preCICE configured..." << endl;

  // init data (fine)
  int i;
  double *u, *u_n, *p, *p_n, *a, *a_n;
  int dimensions = interface.getDimensions();

  u = new double[N + 1]; // Speed
  u_n = new double[N + 1];
  p = new double[N + 1]; // Pressure
  p_n = new double[N + 1];
  a = new double[N + 1];
  a_n = new double[N + 1];

  // init data (coarse)
  double *u_coarse, *u_n_coarse, *p_coarse, *p_n_coarse, *a_coarse, *a_n_coarse;
  double *a_copy_coarse, *p_copy_coarse; // fine mesh

  a_copy_coarse = new double[N + 1]; // coarse displ. data on fine mesh
  p_copy_coarse = new double[N + 1]; // coarse pressure data on fine mesh

  u_coarse = new double[N_SM + 1]; // Speed
  u_n_coarse = new double[N_SM + 1];
  p_coarse = new double[N_SM + 1]; // Pressure
  p_n_coarse = new double[N_SM + 1];
  a_coarse = new double[N_SM + 1];
  a_n_coarse = new double[N_SM + 1];


  // mappings
  NearestNeighborMapping upMapping, downMapping;


  //precice stuff
  int meshID = interface.getMeshID("Fluid_Nodes");
  int aID = interface.getDataID("Displacements", meshID);
  int pID = interface.getDataID("Stresses", meshID);

  int aID_coarse = interface.getDataID("Coarse_Displacements", meshID);
  int pID_coarse = interface.getDataID("Coarse_Stresses", meshID);

  int *vertexIDs, *vertexIDs_coarse;
  vertexIDs = new int[(N + 1)];
  vertexIDs_coarse = new int[(N + 1)];
  double *grid;
  grid = new double[dimensions * (N + 1)];

  // fine model
  for (i = 0; i <= N; i++)
      {
    u[i] = 1.0 / (kappa * 1.0);
    u_n[i] = 1.0 / (kappa * 1.0);
    a[i] = 1.0;
    a_copy_coarse[i] = 1.0;
    a_n[i] = 1.0;
    p[i] = 0.0;
    p_copy_coarse[i] = 0.0;
    p_n[i] = 0.0;
    for (int dim = 0; dim < dimensions; dim++)
      grid[i * dimensions + dim] = i * (1 - dim);
  }
  // surrogate model
  for (i = 0; i <= N_SM; i++)
      {
    u_coarse[i] = 1.0 / (kappa * 1.0);
    u_n_coarse[i] = 1.0 / (kappa * 1.0);
    a_coarse[i] = 1.0;
    a_n_coarse[i] = 1.0;
    p_coarse[i] = 0.0;
    p_n_coarse[i] = 0.0;
    //for (int dim = 0; dim < dimensions; dim++)
    //  grid[i * dimensions + dim] = i * (1 - dim);
  }

  int t = 0; //number of timesteps

  interface.setMeshVertices(meshID, N + 1, grid, vertexIDs);
  interface.setMeshVertices(meshID, N + 1, grid, vertexIDs_coarse); // TODO: ???

  cout << "Fluid: init precice..." << endl;
  interface.initialize();

  if (interface.isActionRequired(actionWriteInitialData())) {
    interface.writeBlockScalarData(pID, N + 1, vertexIDs, p);
    interface.writeBlockScalarData(pID_coarse, N + 1, vertexIDs_coarse, p_copy_coarse); // TODO: ???
    //interface.initializeData();
    interface.fulfilledAction(actionWriteInitialData());
  }

  interface.initializeData();

  if (interface.isReadDataAvailable()) {
    interface.readBlockScalarData(aID, N + 1, vertexIDs, a);
    interface.readBlockScalarData(aID_coarse, N + 1, vertexIDs_coarse, a_copy_coarse); // TODO: ???
  }

  while (interface.isCouplingOngoing()) {
    // When an implicit coupling scheme is used, checkpointing is required
    if (interface.isActionRequired(actionWriteIterationCheckpoint())) {
      interface.fulfilledAction(actionWriteIterationCheckpoint());
    }


    // only surrogate model evaluation for surrogate model optimization
    if(interface.hasToEvaluateSurrogateModel())
    {
      // map down:  fine --> coarse
      downMapping.map(N, N_SM, a_copy_coarse, a_coarse);
      downMapping.map(N, N_SM, p_copy_coarse, p_coarse);

      // ### surrogate model evaluation ###    p_old is not used for gamma = 0.0
      fluid_nl(a_coarse, a_n_coarse, u_coarse, u_n_coarse, p_coarse, p_n_coarse, p_coarse, t + 1, N_SM, kappa, tau, 0.0);

      // map up:  coarse --> fine
      upMapping.map(N_SM, N, a_coarse, a_copy_coarse);
      upMapping.map(N_SM, N, p_coarse, p_copy_coarse);

      // write coarse model response (on fine mesh)
      interface.writeBlockScalarData(pID_coarse, N + 1, vertexIDs_coarse, p_copy_coarse);


      // fine and coarse model evaluation (in MM iteration cycles)
    }else{

      // ### fine model evaluation ###    p_old is not used for gamma = 0.0
      fluid_nl(a, a_n, u, u_n, p, p_n, p, t + 1, N, kappa, tau, 0.0);

      // write fine model response
      interface.writeBlockScalarData(pID, N + 1, vertexIDs, p);

      // map down:  fine --> coarse
      downMapping.map(N, N_SM, a_copy_coarse, a_coarse);
      downMapping.map(N, N_SM, p_copy_coarse, p_coarse);

      // ### surrogate model evaluation ###    p_old is not used for gamma = 0.0
      fluid_nl(a_coarse, a_n_coarse, u_coarse, u_n_coarse, p_coarse, p_n_coarse, p_coarse, t + 1, N_SM, kappa, tau, 0.0);

      // map up:  coarse --> fine
      upMapping.map(N_SM, N, a_coarse, a_copy_coarse);
      upMapping.map(N_SM, N, p_coarse, p_copy_coarse);

      // write coarse model response (on fine mesh)
      interface.writeBlockScalarData(pID_coarse, N + 1, vertexIDs_coarse, p_copy_coarse);


    }

    // perform coupling using preCICE
    interface.advance(0.01);


    // read coupling data from preCICE.

    // only surrogate model evaluation for surrogate model optimization
    if(interface.hasToEvaluateSurrogateModel())
    {
      interface.readBlockScalarData(aID_coarse, N + 1, vertexIDs_coarse, a_copy_coarse);

      // fine and coarse model evaluation (in MM iteration cycles)
    }else{
      interface.readBlockScalarData(aID, N + 1, vertexIDs, a);
      interface.readBlockScalarData(aID_coarse, N + 1, vertexIDs_coarse, a_copy_coarse);
    }

    if (interface.isActionRequired(actionReadIterationCheckpoint())) { // i.e. not yet converged
      // cout << "Iterate" << endl;
      interface.fulfilledAction(actionReadIterationCheckpoint());
    }
    else {
      // cout << "Fluid: Advancing in time, finished timestep: " << t << endl;
      t++;

      for (i = 0; i <= N; i++)
      {
        u_n[i] = u[i];
        p_n[i] = p[i];
        a_n[i] = a[i];
      }
      for (i = 0; i <= N_SM; i++)
      {
        u_n_coarse[i] = u_coarse[i];
        p_n_coarse[i] = p_coarse[i];
        a_n_coarse[i] = a_coarse[i];
      }
    }
  }

  interface.finalize();
  cout << "Exiting FluidSolver" << endl;

  return 0;
}
