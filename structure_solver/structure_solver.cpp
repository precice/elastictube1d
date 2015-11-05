#include <iostream>
#include <stdlib.h>
#include "precice/SolverInterface.hpp"
#include "../mapping/NearestNeighborMapping.hpp"
#include "../mapping/LinearInterpolationMapping.hpp"
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

  int isMultilevelApproach = atoi(argv[argc-1]);
  int N = 0, N_SM = 0;

  if(!isMultilevelApproach)
  {
    if ( argc != 4 )
    {
      cout<<"argc= "<<argc<<std::endl;
      cout << endl;
      cout << "Usage: " << argv[0] << " configurationFileName N" << endl;
      cout << endl;
      cout << "N:     Number of mesh elements, needs to be equal for fluid and structure solver." << endl;
      cout << "isMultilevelApproach: 0/1" << endl;
      return -1;
    }
    N = atoi( argv[2] );
    std::cout << "N: " << N << std::endl;
  }else{

    if ( argc != 5 )
    {
      cout<<"argc= "<<argc<<std::endl;
      cout << endl;
      cout << "Usage: " << argv[0] << " configurationFileName N" << endl;
      cout << endl;
      cout << "N:     Number of mesh elements, needs to be equal for fluid and structure solver." << endl;
      cout << "N_SM:  Number of surrogate model mesh elements, needs to be equal for fluid and structure solver." << endl;
      cout << "isMultilevelApproach: 0/1" << endl;
      return -1;
    }
    N = atoi( argv[2] );
    N_SM = atoi( argv[3] );
    std::cout << "N: " << N <<" N (surrogate model): "<< N_SM << std::endl;
  }
  std::string configFileName(argv[1]);


  std::string dummyName = "STRUCTURE_1D";

  SolverInterface interface(dummyName, 0, 1);
  interface.configure(configFileName);
  cout << "preCICE configured..." << endl;

  //init data (fine model)
  double *displ, *a_n,  *sigma;
  int dimensions = interface.getDimensions();
  displ     = new double[N+1];  // Second dimension (only one cell deep) stored right after the first dimension: see SolverInterfaceImpl::setMeshVertices
  a_n       = new double[N+1];
  sigma     = new double[N+1];
  double *grid;
  grid = new double[dimensions*(N+1)];

  //init data (fine model)
  double *displ_coarse, *a_n_coarse, *sigma_coarse;
  double *displ_copy_coarse;
  double *sigma_copy_coarse;
  // mappings
  NearestNeighborMapping upMapping, downMapping;

  if(isMultilevelApproach)
  {
    displ_coarse      = new double[N_SM+1];  // Second dimension (only one cell deep) stored right after the first dimension: see SolverInterfaceImpl::setMeshVertices
    a_n_coarse        = new double[N_SM+1];
    sigma_coarse      = new double[N_SM+1];
    displ_copy_coarse = new double[N+1];
    sigma_copy_coarse = new double[N+1];
  }

  //precice stuff
  int meshID = interface.getMeshID("Structure_Nodes");
  int displID = interface.getDataID ( "Displacements", meshID );
  int sigmaID = interface.getDataID ( "Stresses", meshID );
  int *vertexIDs;

  int displID_coarse;
  int sigmaID_coarse;

  if(isMultilevelApproach)
  {
    displID_coarse = interface.getDataID ( "Coarse_Displacements", meshID );
    sigmaID_coarse = interface.getDataID ( "Coarse_Stresses", meshID );
  }

  vertexIDs = new int[N+1];


  // fine model init
  for(int i=0; i<=N; i++)
  {
    displ[i] = 0.0;
    a_n[i]   = 1.0;
    sigma[i] = 0.0;

    if(isMultilevelApproach)
      displ_copy_coarse[i] = 0.0;

    for(int dim = 0; dim < dimensions; dim++)
      grid[i*dimensions + dim] = i*(1-dim);   // Define the y-component of each grid point as zero
  }

  if(isMultilevelApproach)
  {
    // surrogate model init
    for(int i=0; i<=N_SM; i++)
    {
      displ_coarse[i] = 0.0;
      a_n_coarse[i]   = 1.0;
      sigma_coarse[i] = 0.0;
    }
  }

  int t = 0;
  interface.setMeshVertices(meshID, N+1, grid, vertexIDs);

  cout << "Structure: init precice..." << endl;
  double dt = interface.initialize();

  if (interface.isActionRequired(actionWriteInitialData()))
  {
    interface.writeBlockScalarData(displID, N+1, vertexIDs, displ);

    if(isMultilevelApproach)
      interface.writeBlockScalarData(displID_coarse, N+1, vertexIDs, displ_copy_coarse);

    //interface.initializeData();
    interface.fulfilledAction(actionWriteInitialData());
  }

  interface.initializeData();


  if (interface.isReadDataAvailable())
  {
    interface.readBlockScalarData(sigmaID , N+1, vertexIDs, sigma);

    if(isMultilevelApproach)
      interface.readBlockScalarData(sigmaID_coarse , N+1, vertexIDs, sigma_copy_coarse);
  }

  while (interface.isCouplingOngoing())
  {
    // When an implicit coupling scheme is used, checkpointing is required
    if (interface.isActionRequired(actionWriteIterationCheckpoint()))
    {
      interface.fulfilledAction(actionWriteIterationCheckpoint());
    }

    // surrogate model evaluation for surrogate model optimization or MM cycle
    if(interface.hasToEvaluateSurrogateModel())
    {
        std::cout<<"\n    ### evaluate coarse model of solid solver, t="<<t<<" ###\n"<<std::endl;

        // map down:  fine --> coarse [displ, pressure]00
        //downMapping.map(N+1, N_SM+1, displ_copy_coarse, displ_coarse);
        downMapping.map(N+1, N_SM+1, sigma_copy_coarse, sigma_coarse);

        // ### surrogate model evaluation ###
        for ( int i = 0; i <= N_SM; i++ ){
          displ_coarse[i]   = 4.0 / ((2.0 - sigma_coarse[i])*(2.0 - sigma_coarse[i]));
          displ_coarse[i]  -= a_n_coarse[i];
        }

        // map up:  coarse --> fine [displ, pressure]
        upMapping.map(N_SM+1, N+1, displ_coarse, displ_copy_coarse);
        //upMapping.map(N_SM+1, N+1, sigma_coarse, sigma_copy_coarse);

        // write coarse model response (on fine mesh)
        interface.writeBlockScalarData(displID_coarse, N+1, vertexIDs, displ_copy_coarse);
    }


    // fine model evaluation (in MM iteration cycles)
    if(interface.hasToEvaluateFineModel())
    {
      std::cout<<"\n    ### evaluate fine model of solid solver, t="<<t<<" ###\n"<<std::endl;

      // ### fine model evaluation ###
      for ( int i = 0; i <= N; i++ ){
        displ[i]   = 4.0 / ((2.0 - sigma[i])*(2.0 - sigma[i]));
        displ[i]  -= a_n[i];
      }
      // write fine model response
      interface.writeBlockScalarData(displID, N+1, vertexIDs, displ);
    }

    // perform coupling using preCICE
    interface.advance(0.01); // Advance by dt = 0.01


    // read coupling data from preCICE.

    // surrogate model evaluation for surrogate model optimization or MM cycle
    if (interface.hasToEvaluateSurrogateModel()){
      interface.readBlockScalarData(sigmaID_coarse, N + 1, vertexIDs, sigma_copy_coarse);
    }

    // fine model evaluation (in MM iteration cycles)
    if(interface.hasToEvaluateFineModel()){
      interface.readBlockScalarData(sigmaID, N+1, vertexIDs, sigma);
    }


    if (interface.isActionRequired(actionReadIterationCheckpoint()))
    {
      // i.e., not converged
      cout << "Iterate" << endl;
      interface.fulfilledAction(actionReadIterationCheckpoint());

    // i.e., converged
    }else{
      // store absolute cross sectionl area value, i.e., Lagrangian solution
      for (int i = 0; i <= N; i++)
        a_n[i] = a_n[i] + displ[i];

      if(isMultilevelApproach)
      {
        // store absolute cross sectionl area value, i.e., Lagrangian solution
        for (int i = 0; i <= N_SM; i++)
          a_n_coarse[i] = a_n_coarse[i] + displ_coarse[i];
      }


      cout << "\n\n ------------------------------------------------\n"
              " Advancing in time, Structure Solver finished timestep: "
           << t <<"\n ------------------------------------------------"<< endl;

      // advance in time
      t++;
    }
  }


  interface.finalize();
  cout << "Exiting StructureSolver" << endl;

  delete [] displ;
  delete [] sigma;
  delete [] grid;

  if (isMultilevelApproach)
  {
    delete[] displ_coarse;
    delete[] sigma_coarse;
    delete[] displ_copy_coarse;
    delete[] sigma_copy_coarse;
  }

  return 0;
}
