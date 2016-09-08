#include "fluid_nl.h"
#include <math.h>
#include <cmath>
#include <cfloat>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>

#include "../utils/AD.hpp"
#include "../utils/Linsolve.hpp"
#include <assert.h>

using namespace AD;

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
    dualReal* crossSectionLength,
    dualReal* crossSectionLength_n,
    dualReal* velocity,
    dualReal* velocity_n,
    dualReal* pressure,
    dualReal* pressure_n,
    dualReal* pressure_old,
    double scaled_t,
    int N,
    dualReal kappa,
    dualReal tau,
    dualReal gamma)
{
  /* fluid_nl Variables */
  int i, j = 0, k, ampl;
  dualReal alpha, dx;
  dualReal tmp, tmp2;
  dualReal* Res, *X;
  dualReal** LHS;
  dualReal temp_sum;
  double norm_1, norm_2;
  double norm = 1.0;

  // Used as Ax = b
  // i.e. LHS*x = Res
  Res = (dualReal*)calloc((2 * N + 2), sizeof(dualReal));
  X   = (dualReal*)calloc((2 * N + 2), sizeof(dualReal)); 
  LHS = (dualReal**)calloc((2 * N + 2), sizeof(dualReal*));
  for (i = 0; i < (2 * N + 2); ++i)
    LHS[i] = (dualReal*)calloc((2 * N + 2), sizeof(dualReal));

  /* LAPACK Variables */
  
  /*
  double* A = (double*)calloc((2 * N + 2) * (2 * N + 2), sizeof(double));
  double* res = (double*)calloc((2 * N + 2), sizeof(double));
  ;
  int* ipiv = (int*)calloc((2 * N + 2), sizeof(int));
  ;
  int nlhs = (2 * N + 2);
  int nrhs = 1;
  int info;
  */

  /* Stabilization Intensity */
  alpha = (dualReal(N) * kappa * tau) / (dualReal(N) * tau + dualReal(1.0));
  dx = dualReal(1.0) / (dualReal(N) * kappa * tau);
  ampl = 100;

  k = 0;
  /* Stopping Criteria */
  // TODO: Include the residual
  while (1) {
    for (i = 0; i < (2 * N + 2); i++)
      Res[i] = dualReal(0.0);

    for (i = 1; i < N; i++) {
      /* Momentum */
      Res[i] = velocity_n[i] * crossSectionLength[i] * dx;
      Res[i] = Res[i] - dualReal(0.25) * crossSectionLength[i + 1] * velocity[i] * velocity[i + 1] - dualReal(0.25) * crossSectionLength[i] * velocity[i] * velocity[i + 1];
      Res[i] = Res[i] - crossSectionLength[i] * dx * velocity[i] - dualReal(0.25) * crossSectionLength[i + 1] * velocity[i] * velocity[i] - dualReal(0.25) * crossSectionLength[i] * velocity[i] * velocity[i] + dualReal(0.25) * crossSectionLength[i] * velocity[i - 1] * velocity[i] + dualReal(0.25) * crossSectionLength[i - 1] * velocity[i - 1] * velocity[i];
      Res[i] = Res[i] + dualReal(0.25) * crossSectionLength[i - 1] * velocity[i - 1] * velocity[i - 1] + dualReal(0.25) * crossSectionLength[i] * velocity[i - 1] * velocity[i - 1]; 
      Res[i] = Res[i] + dualReal(0.25) * crossSectionLength[i - 1] * pressure[i - 1] + dualReal(0.25) * crossSectionLength[i] * pressure[i - 1] + dualReal(0.25) * crossSectionLength[i - 1] * pressure[i] - dualReal(0.25) * crossSectionLength[i + 1] * pressure[i] - dualReal(0.25) * crossSectionLength[i] * pressure[i + 1] - dualReal(0.25) * crossSectionLength[i + 1] * pressure[i + 1];

      /* Continuity */
      Res[i + N + 1] = -(crossSectionLength[i] - crossSectionLength_n[i]) * dx + pressure_old[i] * gamma * dx;
      Res[i + N + 1] = Res[i + N + 1] + dualReal(0.25) * crossSectionLength[i - 1] * velocity[i - 1] + dualReal(0.25) * crossSectionLength[i] * velocity[i - 1] + dualReal(0.25) * crossSectionLength[i - 1] * velocity[i] - dualReal(0.25) * crossSectionLength[i + 1] * velocity[i] - dualReal(0.25) * crossSectionLength[i] * velocity[i + 1] - dualReal(0.25) * crossSectionLength[i + 1] * velocity[i + 1];
      Res[i + N + 1] = Res[i + N + 1] + alpha * pressure[i - 1] - dualReal(2) * alpha * pressure[i] - gamma * pressure[i] * dx + alpha * pressure[i + 1];
    }

    i = N - 1;

    /* Boundary */

    /* Velocity Inlet is prescribed */
    tmp = AD::sin(dualReal(PI*scaled_t));
    Res[0] = (dualReal(1.0) / kappa) + (dualReal(1.0) / (kappa * ampl)) * tmp * tmp - velocity[0];

    /* Pressure Inlet is lineary interpolated */
    Res[N + 1] = -pressure[0] + dualReal(2.) * pressure[1] - pressure[2];

    /* Velocity Outlet is lineary interpolated */
    Res[N] = -velocity[N] + dualReal(2.) * velocity[N - 1] - velocity[N - 2];

    /* Pressure Outlet is "non-reflecting" */
    tmp2 = sqrt(dualReal(1.) - pressure_n[N] / dualReal(2.)) - (velocity[N] - velocity_n[N]) / dualReal(4.);
    Res[2 * N + 1] = -pressure[N] + dualReal(2.) * (dualReal(1.) - tmp2 * tmp2);

    /* Stopping Criteria */
    k += 1; // Iteration Count

    temp_sum.u = 0;
    for (i = 0; i < (2 * N + 2); i++) {
      temp_sum = temp_sum + Res[i] * Res[i];
    }
    norm_1 = sqrt(temp_sum.u);

    temp_sum = 0;
    for (i = 0; i < (N + 1); i++) {
      temp_sum = temp_sum + (pressure[i] * pressure[i]) + (velocity[i] * velocity[i]);
    }
    norm_2 = sqrt(temp_sum.u);

    //std::cout<<"norm1: "<<norm_1<<", norm2: "<<norm_2<<std::endl;
    
    norm = norm_1 / norm_2; // Norm

    if ((norm < 1e-15 && k > 1) || k > 50) {
      printf("Nonlinear Solver break, Its: %i, norm: %e\n", k, norm);
      break;
    }

    /* Initilizing the the LHS i.e. Left Hand Side */
    for (i = 0; i <= (2 * N + 1); i++)
      for (j = 0; j <= (2 * N + 1); j++)
        LHS[i][j] = dualReal(0.0);

    for (i = 1; i < N; i++) {
      // Momentum, Velocity
      LHS[i][i - 1] = LHS[i][i - 1] - dualReal(0.25) * crossSectionLength[i - 1] * velocity[i - 1] * dualReal(2.) - dualReal(0.25) * crossSectionLength[i] * velocity[i - 1] * dualReal(2.) - dualReal(0.25) * crossSectionLength[i] * velocity[i] + dualReal(0.25) * crossSectionLength[i - 1] * velocity[i];
      LHS[i][i]     = LHS[i][i] + dualReal(0.25) * crossSectionLength[i + 1] * velocity[i + 1] + dualReal(0.25) * crossSectionLength[i] * velocity[i + 1] + crossSectionLength[i] * dx + dualReal(0.25) * crossSectionLength[i + 1] * velocity[i] * dualReal(2.) + dualReal(0.25) * crossSectionLength[i] * velocity[i] * dualReal(2.) - dualReal(0.25) * crossSectionLength[i] * velocity[i - 1] - dualReal(0.25) * crossSectionLength[i - 1] * velocity[i - 1];
      LHS[i][i + 1] = LHS[i][i + 1] + dualReal(0.25) * crossSectionLength[i + 1] * velocity[i] + dualReal(0.25) * crossSectionLength[i] * velocity[i];

      // Momentum, Pressure
      LHS[i][N + 1 + i - 1] = LHS[i][N + 1 + i - 1] - dualReal(0.25) * crossSectionLength[i - 1] - dualReal(0.25) * crossSectionLength[i];
      LHS[i][N + 1 + i]     = LHS[i][N + 1 + i] + dualReal(0.25) * crossSectionLength[i - 1] - dualReal(0.25) * crossSectionLength[i + 1];
      LHS[i][N + 1 + i + 1] = LHS[i][N + 1 + i + 1] + dualReal(0.25) * crossSectionLength[i] + dualReal(0.25) * crossSectionLength[i + 1];

      // Continuity, Velocity
      LHS[i + N + 1][i - 1] = LHS[i + N + 1][i - 1] - dualReal(0.25) * crossSectionLength[i - 1] - dualReal(0.25) * crossSectionLength[i];
      LHS[i + N + 1][i]     = LHS[i + N + 1][i] - dualReal(0.25) * crossSectionLength[i - 1] + dualReal(0.25) * crossSectionLength[i + 1];
      LHS[i + N + 1][i + 1] = LHS[i + N + 1][i + 1] + dualReal(0.25) * crossSectionLength[i] + dualReal(0.25) * crossSectionLength[i + 1];

      // Continuity, Pressure
      LHS[i + N + 1][N + 1 + i - 1] = LHS[i + N + 1][N + 1 + i - 1] - alpha;
      LHS[i + N + 1][N + 1 + i]     = LHS[i + N + 1][N + 1 + i] + dualReal(2.) * alpha + gamma * dx;
      LHS[i + N + 1][N + 1 + i + 1] = LHS[i + N + 1][N + 1 + i + 1] - alpha;
    }

    /* Boundary */

    // Velocity Inlet is prescribed
    LHS[0][0] = dualReal(1.);
    // Pressure Inlet is lineary interpolated
    LHS[N + 1][N + 1] = dualReal(1.);
    LHS[N + 1][N + 2] = dualReal(-2.);
    LHS[N + 1][N + 3] = dualReal(1.);
    // Velocity Outlet is lineary interpolated
    LHS[N][N] = dualReal(1.);
    LHS[N][N - 1] = dualReal(-2.);
    LHS[N][N - 2] = dualReal(1.);
    // Pressure Outlet is Non-Reflecting
    LHS[2 * N + 1][2 * N + 1] = dualReal(1.);
    LHS[2 * N + 1][N] = -(AD::sqrt(dualReal(1.) - pressure_n[N] / dualReal(2.)) - (velocity[N] - velocity_n[N]) / dualReal(4.));

     
    // Solve Linear System using LAPACK

    /* LAPACK requires a 1D array 
       i.e. Linearizing 2D 
    */
    
    /*
    int counter = 0;
    for (i = 0; i <= (2 * N + 1); i++) {
      for (j = 0; j <= (2 * N + 1); j++) {
        A[counter] = LHS[j][i].u;
        counter++;
      }
    }
    
    for (i = 0; i <= (2 * N + 1); i++) 
     res[i] = Res[i].u;

    // LAPACK Function call to solve the linear system 
    dgesv_(&nlhs, &nrhs, A, &nlhs, ipiv, res, &nlhs, &info);
    
    if (info != 0) {
      printf("Linear Solver not converged!!!, Info: %i\n", info);
    }
    
    std::vector<double> abl(2*N+2);
    
    for (i = 0; i <= (2 * N + 1); i++) {
     Res[i].u = res[i];
     
     abl[i] = 0.;     
     for (j = 0; j <= (2 * N + 1); j++) {
       // A*d/dx1 x 
       abl[i] = abl[i] +  Res[j].v * LHS[i][j].u;
     }
    }
    
    for (i = 0; i <= (2 * N + 1); i++)
      Res[i].v = abl[i];
    */
    
    
    LA::linsolve((2 * N + 2), (2 * N + 2), LHS, Res, X);
    
    //exit(0);
    /*
     for (i = 0; i <= (2 * N + 1); i++){
      for (j = 0; j <= (2 * N + 1); j++){
	if(std::isnan(LHS[i][j].u + LHS[i][j].v)){
          std::cout<<"nan: ("<<LHS[i][j].u<<", "<<LHS[i][j].v<<")"<<std::endl;}
      }
      //if(std::isnan(X[i].u + X[i].v))
        std::cout<<"X["<<i<<"] = ("<<X[i].u<<", "<<X[i].v<<")"<<std::endl;
    }
    exit(0);
    */

    for (i = 0; i <= N; i++) {
      velocity[i] = velocity[i] + X[i];
      pressure[i] = pressure[i] + X[i + N + 1];
    }

  } // END OF WHILE


  return 0;
}
