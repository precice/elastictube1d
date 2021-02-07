#include "fluid_nl.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>

/* 
   Function for solving the linear system 
   LAPACK is used DGESV computes the solution to a real system of linear equations
   A * x = b,
   where A is an N-by-N matrix and x and b are N-by-NRHS matrices.
*/
extern "C" {
void dgesv_(
    int* n,
    int* nrhs,
    double* A,
    int* lda,
    int* ipiv,
    double* b,
    int* ldb,
    int* info);
}

/* Function for fluid_nl i.e. non-linear */
int fluid_nl(
    double* crossSectionLength,
    double* velocity,
    double* pressure,
    double t,
    int N,
    double kappa,
    double tau)
{
  /* fluid_nl Variables */
  int i, j = 0, k, ampl;
  double alpha, dx;
  double tmp, tmp2;
  double* Res;
  double** LHS;
  double temp_sum;
  double norm_1, norm_2;
  double norm = 1.0;

  // Used as Ax = b
  // i.e. LHS*x = Res
  Res = (double*)calloc((2 * N + 2), sizeof(double));
  LHS = (double**)calloc((2 * N + 2), sizeof(double*));
  for (i = 0; i < (2 * N + 2); ++i)
    LHS[i] = (double*)calloc((2 * N + 2), sizeof(double));

  /* LAPACK Variables */
  double* A = (double*)calloc((2 * N + 2) * (2 * N + 2), sizeof(double));
  ;
  int* ipiv = (int*)calloc((2 * N + 2), sizeof(int));
  ;
  int nlhs = (2 * N + 2);
  int nrhs = 1;
  int info;

  /* Stabilization Intensity */
  alpha = (N * kappa * tau) / (N * tau + 1);
  dx = 1.0 / (N * kappa * tau);
  ampl = 100;

  k = 0;
  while (1) {
    for (i = 0; i < (2 * N + 2); i++)
      Res[i] = 0.0;

    for (i = 1; i < N; i++) {
      /* Momentum */
      Res[i] = velocity[i] * crossSectionLength[i] * dx;
      Res[i] = Res[i] - 0.25 * crossSectionLength[i + 1] * velocity[i] * velocity[i + 1] - 0.25 * crossSectionLength[i] * velocity[i] * velocity[i + 1];
      Res[i] = Res[i] - crossSectionLength[i] * dx * velocity[i] - 0.25 * crossSectionLength[i + 1] * velocity[i] * velocity[i] - 0.25 * crossSectionLength[i] * velocity[i] * velocity[i] + 0.25 * crossSectionLength[i] * velocity[i - 1] * velocity[i] + 0.25 * crossSectionLength[i - 1] * velocity[i - 1] * velocity[i];
      Res[i] = Res[i] + 0.25 * crossSectionLength[i - 1] * velocity[i - 1] * velocity[i - 1] + 0.25 * crossSectionLength[i] * velocity[i - 1] * velocity[i - 1];
      Res[i] = Res[i] + 0.25 * crossSectionLength[i - 1] * pressure[i - 1] + 0.25 * crossSectionLength[i] * pressure[i - 1] - 0.25 * crossSectionLength[i - 1] * pressure[i] + 0.25 * crossSectionLength[i + 1] * pressure[i] - 0.25 * crossSectionLength[i] * pressure[i + 1] - 0.25 * crossSectionLength[i + 1] * pressure[i + 1];

      /* Continuity */
      Res[i + N + 1] = -(crossSectionLength[i] - crossSectionLength[i]) * dx;
      Res[i + N + 1] = Res[i + N + 1] + 0.25 * crossSectionLength[i - 1] * velocity[i - 1] + 0.25 * crossSectionLength[i] * velocity[i - 1] + 0.25 * crossSectionLength[i - 1] * velocity[i] - 0.25 * crossSectionLength[i + 1] * velocity[i] - 0.25 * crossSectionLength[i] * velocity[i + 1] - 0.25 * crossSectionLength[i + 1] * velocity[i + 1];
      Res[i + N + 1] = Res[i + N + 1] + alpha * pressure[i - 1] - 2 * alpha * pressure[i] + alpha * pressure[i + 1];
    }

    i = N - 1;

    /* Boundary */

    /* Velocity Inlet is prescribed */
    tmp = sin(PI * (t+0.01)); //to not start with 0 velocity
    Res[0] = (1.0 / kappa) + (1.0 / (kappa * ampl)) * tmp * tmp - velocity[0];

    /* Pressure Inlet is lineary interpolated */
    Res[N + 1] = -pressure[0] + 2 * pressure[1] - pressure[2];

    /* Velocity Outlet is lineary interpolated */
    Res[N] = -velocity[N] + 2 * velocity[N - 1] - velocity[N - 2];

    /* Pressure Outlet is "non-reflecting" */
    tmp2 = sqrt(1 - pressure[N] / 2) - (velocity[N] - velocity[N]) / 4;
    Res[2 * N + 1] = -pressure[N] + 2 * (1 - tmp2 * tmp2);

    k += 1; // Iteration Count

    // compute norm of residual
    temp_sum = 0;
    for (i = 0; i < (2 * N + 2); i++) {
      temp_sum += Res[i] * Res[i];
    }
    norm_1 = sqrt(temp_sum);
    temp_sum = 0;
    for (i = 0; i < (N + 1); i++) {
      temp_sum += (pressure[i] * pressure[i]) + (velocity[i] * velocity[i]);
    }
    norm_2 = sqrt(temp_sum);
    norm = norm_1 / norm_2; 

    if ((norm < 1e-15 && k > 1) || k > 50) {
      printf("Nonlinear Solver break, iterations: %i, residual norm: %e\n", k, norm);
      break;
    }

    /* Initilizing the the LHS i.e. Left Hand Side */
    for (i = 0; i <= (2 * N + 1); i++)
      for (j = 0; j <= (2 * N + 1); j++)
        LHS[i][j] = 0.0;

    for (i = 1; i < N; i++) {
      // Momentum, Velocity
      LHS[i][i - 1] = LHS[i][i - 1] - 0.25 * crossSectionLength[i - 1] * velocity[i - 1] * 2 - 0.25 * crossSectionLength[i] * velocity[i - 1] * 2 - 0.25 * crossSectionLength[i] * velocity[i] - 0.25 * crossSectionLength[i - 1] * velocity[i];
      LHS[i][i] = LHS[i][i] + 0.25 * crossSectionLength[i + 1] * velocity[i + 1] + 0.25 * crossSectionLength[i] * velocity[i + 1] + crossSectionLength[i] * dx + 0.25 * crossSectionLength[i + 1] * velocity[i] * 2 + 0.25 * crossSectionLength[i] * velocity[i] * 2 - 0.25 * crossSectionLength[i] * velocity[i - 1] - 0.25 * crossSectionLength[i - 1] * velocity[i - 1];
      LHS[i][i + 1] = LHS[i][i + 1] + 0.25 * crossSectionLength[i + 1] * velocity[i] + 0.25 * crossSectionLength[i] * velocity[i];

      // Momentum, Pressure
      LHS[i][N + 1 + i - 1] = LHS[i][N + 1 + i - 1] - 0.25 * crossSectionLength[i - 1] - 0.25 * crossSectionLength[i];
      LHS[i][N + 1 + i] = LHS[i][N + 1 + i] + 0.25 * crossSectionLength[i - 1] - 0.25 * crossSectionLength[i + 1];
      LHS[i][N + 1 + i + 1] = LHS[i][N + 1 + i + 1] + 0.25 * crossSectionLength[i] + 0.25 * crossSectionLength[i + 1];

      // Continuity, Velocity
      LHS[i + N + 1][i - 1] = LHS[i + N + 1][i - 1] - 0.25 * crossSectionLength[i - 1] - 0.25 * crossSectionLength[i];
      LHS[i + N + 1][i] = LHS[i + N + 1][i] - 0.25 * crossSectionLength[i - 1] + 0.25 * crossSectionLength[i + 1];
      LHS[i + N + 1][i + 1] = LHS[i + N + 1][i + 1] + 0.25 * crossSectionLength[i] + 0.25 * crossSectionLength[i + 1];

      // Continuity, Pressure
      LHS[i + N + 1][N + 1 + i - 1] = LHS[i + N + 1][N + 1 + i - 1] - alpha;
      LHS[i + N + 1][N + 1 + i] = LHS[i + N + 1][N + 1 + i] + 2 * alpha;
      LHS[i + N + 1][N + 1 + i + 1] = LHS[i + N + 1][N + 1 + i + 1] - alpha;
    }

    /* Boundary */

    // Velocity Inlet is prescribed
    LHS[0][0] = 1;
    // Pressure Inlet is lineary interpolated
    LHS[N + 1][N + 1] = 1;
    LHS[N + 1][N + 2] = -2;
    LHS[N + 1][N + 3] = 1;
    // Velocity Outlet is lineary interpolated
    LHS[N][N] = 1;
    LHS[N][N - 1] = -2;
    LHS[N][N - 2] = 1;
    // Pressure Outlet is Non-Reflecting
    LHS[2 * N + 1][2 * N + 1] = 1;
    LHS[2 * N + 1][N] = -(sqrt(1 - pressure[N] / 2.0) - (velocity[N] - velocity[N]) / 4.0);

    /* LAPACK requires a 1D array 
       i.e. Linearizing 2D 
    */
    int counter = 0;
    for (i = 0; i <= (2 * N + 1); i++) {
      for (j = 0; j <= (2 * N + 1); j++) {
        A[counter] = LHS[j][i];
        counter++;
      }
    }

    /* LAPACK Function call to solve the linear system */
    dgesv_(&nlhs, &nrhs, A, &nlhs, ipiv, Res, &nlhs, &info);

    if (info != 0) {
      printf("Linear Solver not converged!, Info: %i\n", info);
    }

    for (i = 0; i <= N; i++) {
      velocity[i] = velocity[i] + Res[i];
      pressure[i] = pressure[i] + Res[i + N + 1];
    }

  } 
  return 0;
}

void initializeWriting(std::ofstream&filestream)
{
  filestream.setf(std::ios::showpoint);
  filestream.setf(std::ios::scientific);
  filestream << std::setprecision(16);
}

void writeHeader(std::ostream& outFile)
{
  outFile << "# vtk DataFile Version 2.0" << std::endl << std::endl
          << "ASCII" << std::endl << std::endl
          << "DATASET UNSTRUCTURED_GRID" << std::endl << std::endl;
}

void exportMesh(std::ofstream& outFile, int N_slices, double* grid)
{  
  // Plot vertices
  outFile << "POINTS " << N_slices << " float "<<std::endl << std::endl;
  
  for (int i = 0; i<N_slices; i++)
  {
	  // read x,y from grid. Set z = 0
	  // Values are stored in grid in the following way: [x_0,y_0,x_1,y_1,...x_n-1,y_n-1]	 
	  double x = grid[2 * i + 0]; 
	  double y = grid[2 * i + 1];
	  double z = 0.0;
	  outFile << x << "  " << y << "  " << z << std::endl;
  }
  outFile << std::endl;
}

void exportVectorData(std::ofstream& outFile, int N_slices, double* data, const char* dataname)
{
	outFile << "VECTORS " << dataname << " float" << std::endl;

	for(int i = 0; i < N_slices; i++)
	{ 	
		// Plot vertex data 
		// read x vector component from dataset. Set y,z = 0
		// Values are stored in dataset in the following way: [vx_0,vx_1,...vx_n-1]
		double vx = data[i]; 
		double vy = 0.0; 
		double vz = 0.0;           
		outFile << vx << "  " << vy << "  " << vz << std::endl;      
	}  
	
	outFile << std::endl;          
}

void exportScalarData(std::ofstream& outFile, int N_slices, double* data, std::string dataname)
{  
	outFile << "SCALARS " << dataname << " float" << std::endl;
	outFile << "LOOKUP_TABLE default" << std::endl;

	for(int i = 0; i < N_slices; i++)
	{ 
		// Plot vertex data            
		outFile << data[i] << std::endl;      
	}
	
	outFile << std::endl;
}

void write_vtk(double t, int iteration, const char* filename_prefix, int N_slices, double* grid, double* velocity, double* pressure, double* diameter)
{
	std::stringstream filename_stream;
	filename_stream << filename_prefix <<"_"<< iteration <<".vtk";
	std::string filename = filename_stream.str();
	printf("writing timestep at t=%f to %s\n", t, filename.c_str());
				
	std::ofstream outstream(filename);	

	initializeWriting(outstream);
	writeHeader(outstream);
	exportMesh(outstream, N_slices, grid);
	
	outstream << "POINT_DATA " << N_slices << std::endl;
	outstream << std::endl;  
  	
	exportVectorData(outstream, N_slices, velocity, "velocity");
	exportScalarData(outstream, N_slices, pressure, "pressure");
	exportScalarData(outstream, N_slices, diameter, "diameter");
	
	outstream.close();		
}
