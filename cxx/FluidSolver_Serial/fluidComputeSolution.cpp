#include "FluidSolver.h"

#include <iostream>
#include <cmath>
#include <mpi.h>
#include <vector>
#include <iomanip>

using std::sin;
using std::sqrt;

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

void fluidComputeSolution(
    int rank,
    int size,
    int domainSize,
    int chunkLength,
    double kappa,
    double tau,
    double gamma,
    double scaled_t,
    double* pressure,
    double* crossSectionLength,
    double* velocity)
{
  /*
   * Step 1: Recieve the complete dataset in process 0.
   */
  if (rank != 0) {
    int N = domainSize;
    int i = rank;
    double chunkLength_temp;
    if ((N + 1) % size == 0) {
      chunkLength_temp = (N + 1) / size;
    } else if (i < (N + 1) % size) {
      chunkLength_temp = (N + 1) / size + 1;
    } else {
      chunkLength_temp = (N + 1) / size;
    }
    chunkLength =chunkLength_temp;
    double* velo = new double[chunkLength];
    double* press = new double[chunkLength];
    double* diam = new double[chunkLength];


    int tagStart = 7 * rank; // dont needed
    MPI_Send(pressure, chunkLength, MPI_DOUBLE, 0, tagStart + 0, MPI_COMM_WORLD);
    //MPI_Send(pressure_n, chunkLength, MPI_DOUBLE, 0, tagStart + 1, MPI_COMM_WORLD);
    //MPI_Send(pressure_old, chunkLength, MPI_DOUBLE, 0, tagStart + 2, MPI_COMM_WORLD);
    MPI_Send(crossSectionLength, chunkLength, MPI_DOUBLE, 0, tagStart + 3, MPI_COMM_WORLD);
    //MPI_Send(crossSectionLength_n, chunkLength, MPI_DOUBLE, 0, tagStart + 4, MPI_COMM_WORLD);
    MPI_Send(velocity, chunkLength, MPI_DOUBLE, 0, tagStart + 5, MPI_COMM_WORLD);
    //MPI_Send(velocity_n, chunkLength, MPI_DOUBLE, 0, tagStart + 6, MPI_COMM_WORLD);

    MPI_Status status;
    MPI_Recv(press, chunkLength, MPI_DOUBLE, 0, tagStart + 0, MPI_COMM_WORLD, &status);
    MPI_Recv(diam, chunkLength, MPI_DOUBLE, 0, tagStart + 3, MPI_COMM_WORLD, &status);
    MPI_Recv(velo, chunkLength, MPI_DOUBLE, 0, tagStart + 5, MPI_COMM_WORLD, &status);
    
    for (int k = 0; k< chunkLength; k++){
      velocity[k] = velo[k];
      pressure[k] = press[k];
      crossSectionLength[k] = diam[k];
    }
    
  } else {
    double *pressure_NLS;
    double *crossSectionLength_NLS;
    double *velocity_NLS;

    int N = domainSize;

    pressure_NLS = new double[N + 1];

    crossSectionLength_NLS = new double[N + 1];

    velocity_NLS = new double[N + 1];

    for (int i = 0; i < chunkLength; i++) {
      pressure_NLS[i] = pressure[i];
      //pressure_n_NLS[i] = pressure_n[i];
      //pressure_old_NLS[i] = pressure_old[i];
      crossSectionLength_NLS[i] = crossSectionLength[i];
      //crossSectionLength_n_NLS[i] = crossSectionLength_n[i];
      velocity_NLS[i] = velocity[i];
      //velocity_n_NLS[i] = velocity_n[i];
    }

    for (int i = 1; i < size; i++) {
      int tagStart = 7 * i;
      int chunkLength_temp;
      int gridOffset;

      if ((N + 1) % size == 0) {
        chunkLength_temp = (N + 1) / size;
        gridOffset = i * chunkLength_temp;
      } else if (i < (N + 1) % size) {
        chunkLength_temp = (N + 1) / size + 1;
        gridOffset = i * chunkLength_temp;
      } else {
        chunkLength_temp = (N + 1) / size;
        gridOffset = ((N + 1) % size) * ((N + 1) / size + 1) + (i - ((N + 1) % size)) * (N + 1) / size;
      }

      MPI_Recv(pressure_NLS + gridOffset, chunkLength_temp, MPI_DOUBLE, i, tagStart + 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      //MPI_Recv(pressure_n_NLS + gridOffset, chunkLength_temp, MPI_DOUBLE, i, tagStart + 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      //MPI_Recv(pressure_old_NLS + gridOffset, chunkLength_temp, MPI_DOUBLE, i, tagStart + 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(crossSectionLength_NLS + gridOffset, chunkLength_temp, MPI_DOUBLE, i, tagStart + 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      //MPI_Recv(crossSectionLength_n_NLS + gridOffset, chunkLength_temp, MPI_DOUBLE, i, tagStart + 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(velocity_NLS + gridOffset, chunkLength_temp, MPI_DOUBLE, i, tagStart + 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      //MPI_Recv(velocity_n_NLS + gridOffset, chunkLength_temp, MPI_DOUBLE, i, tagStart + 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }


    // LAPACK Variables here
    double *Res, **LHS, *A, alpha, dx, tmp, tmp2, temp_sum, norm_1, norm_2, norm = 1.0;
    int *ipiv, nlhs = 2 * N + 2, nrhs = 1, info, ampl;

    Res = new double[2 * N + 2];
    LHS = new double*[2 * N + 2];
    for (int i = 0; i < 2 * N + 2; i++) {
      LHS[i] = new double[2 * N + 2];
    }
    A = new double[(2 * N + 2) * (2 * N + 2)];
    ipiv = new int[2 * N + 2];

    // Stabilization intensity
    alpha = (N * kappa * tau) / (N * tau + 1);
    dx = 1.0 / (N * kappa * tau);
    ampl = 100;

    int whileLoopCounter = 0;
    while (1) { // Add stopping criterion
      for (int i = 0; i < 2 * N + 2; i++) {
        Res[i] = 0.0;
      }

      for (int i = 1; i < N; i++) { // chnaged to 1
        /* Momentum */
        Res[i] = velocity_NLS[i] * crossSectionLength_NLS[i] * dx;
        Res[i] = Res[i] - 0.25 * crossSectionLength_NLS[i + 1] * velocity_NLS[i] * velocity_NLS[i + 1] - 0.25 * crossSectionLength_NLS[i] * velocity_NLS[i] * velocity_NLS[i + 1];
        Res[i] = Res[i] - crossSectionLength_NLS[i] * dx * velocity_NLS[i] - 0.25 * crossSectionLength_NLS[i + 1] * velocity_NLS[i] * velocity_NLS[i] - 0.25 * crossSectionLength_NLS[i] * velocity_NLS[i] * velocity_NLS[i] + 0.25 * crossSectionLength_NLS[i] * velocity_NLS[i - 1] * velocity_NLS[i] + 0.25 * crossSectionLength_NLS[i - 1] * velocity_NLS[i - 1] * velocity_NLS[i];
        Res[i] = Res[i] + 0.25 * crossSectionLength_NLS[i - 1] * velocity_NLS[i - 1] * velocity_NLS[i - 1] + 0.25 * crossSectionLength_NLS[i] * velocity_NLS[i - 1] * velocity_NLS[i - 1];
        Res[i] = Res[i] + 0.25 * crossSectionLength_NLS[i - 1] * pressure_NLS[i - 1] + 0.25 * crossSectionLength_NLS[i] * pressure_NLS[i - 1] - 0.25 * crossSectionLength_NLS[i - 1] * pressure_NLS[i] + 0.25 * crossSectionLength_NLS[i + 1] * pressure_NLS[i] - 0.25 * crossSectionLength_NLS[i] * pressure_NLS[i + 1] - 0.25 * crossSectionLength_NLS[i + 1] * pressure_NLS[i + 1];

        /* Continuity */
        //Res[i + N + 1] = -(crossSectionLength_NLS[i] - crossSectionLength_n_NLS[i]) * dx + pressure_old_NLS[i] * gamma * dx;
        Res[i + N + 1] = -(crossSectionLength_NLS[i] - crossSectionLength_NLS[i]) * dx;
        Res[i + N + 1] = Res[i + N + 1] + 0.25 * crossSectionLength_NLS[i - 1] * velocity_NLS[i - 1] + 0.25 * crossSectionLength_NLS[i] * velocity_NLS[i - 1] + 0.25 * crossSectionLength_NLS[i - 1] * velocity_NLS[i] - 0.25 * crossSectionLength_NLS[i + 1] * velocity_NLS[i] - 0.25 * crossSectionLength_NLS[i] * velocity_NLS[i + 1] - 0.25 * crossSectionLength_NLS[i + 1] * velocity_NLS[i + 1];
        //Res[i + N + 1] = Res[i + N + 1] + alpha * pressure_NLS[i - 1] - 2 * alpha * pressure_NLS[i] - gamma * pressure_NLS[i] * dx + alpha * pressure_NLS[i + 1];
        Res[i + N + 1] = Res[i + N + 1] + alpha * pressure_NLS[i - 1] - 2 * alpha * pressure_NLS[i] + alpha * pressure_NLS[i + 1];

      }

      // Boundary

      // Velocity
      tmp = sin(PI * scaled_t );
      Res[0] = (1.0 / kappa) + (1.0 / (kappa * ampl)) * tmp * tmp - velocity_NLS[0];

      // Pressure Inlet is lineary interpolated
      Res[N + 1] = -pressure_NLS[0] + 2 * pressure_NLS[1] - pressure_NLS[2];

      // Velocity Outlet is lineary interpolated
      Res[N] = -velocity_NLS[N] + 2 * velocity_NLS[N - 1] - velocity_NLS[N - 2];

      // Pressure Outlet is "non-reflecting"
      tmp2 = sqrt(1 - pressure_NLS[N] / 2) - (velocity_NLS[N] - velocity_NLS[N]) / 4;
      Res[2 * N + 1] = -pressure_NLS[N] + 2 * (1 - tmp2 * tmp2);

      // Stopping Criteria
      whileLoopCounter += 1; // Iteration Count

      temp_sum = 0;
      for (int i = 0; i < (2 * N + 2); i++) {
        temp_sum += Res[i] * Res[i];
      }
      norm_1 = sqrt(temp_sum);

      temp_sum = 0;
      for (int i = 0; i < (N + 1); i++) {
        temp_sum += (pressure_NLS[i] * pressure_NLS[i]) + (velocity_NLS[i] * velocity_NLS[i]);
      }
      norm_2 = sqrt(temp_sum);

      norm = norm_1 / norm_2; // Norm

      if ((norm < 1e-15 && whileLoopCounter > 1) || whileLoopCounter > 50) {
        std::cout << "Nonlinear Solver break, Its: " << whileLoopCounter << ", norm: " << norm << std::endl;
        break;
      }

      // Initilizing the the LHS i.e. Left Hand Side
      for (int i = 0; i <= (2 * N + 1); i++) {
        for (int j = 0; j <= (2 * N + 1); j++) {
          LHS[i][j] = 0.0;
        }
      }

      for (int i = 1; i < N; i++) {
        // Momentum, Velocity
        LHS[i][i - 1] = LHS[i][i - 1] - 0.25 * crossSectionLength_NLS[i - 1] * velocity_NLS[i - 1] * 2 - 0.25 * crossSectionLength_NLS[i] * velocity_NLS[i - 1] * 2 - 0.25 * crossSectionLength_NLS[i] * velocity_NLS[i] - 0.25 * crossSectionLength_NLS[i - 1] * velocity_NLS[i];
        LHS[i][i] = LHS[i][i] + 0.25 * crossSectionLength_NLS[i + 1] * velocity_NLS[i + 1] + 0.25 * crossSectionLength_NLS[i] * velocity_NLS[i + 1] + crossSectionLength_NLS[i] * dx + 0.25 * crossSectionLength_NLS[i + 1] * velocity_NLS[i] * 2 + 0.25 * crossSectionLength_NLS[i] * velocity_NLS[i] * 2 - 0.25 * crossSectionLength_NLS[i] * velocity_NLS[i - 1] - 0.25 * crossSectionLength_NLS[i - 1] * velocity_NLS[i - 1];
        LHS[i][i + 1] = LHS[i][i + 1] + 0.25 * crossSectionLength_NLS[i + 1] * velocity_NLS[i] + 0.25 * crossSectionLength_NLS[i] * velocity_NLS[i];

        // Momentum, Pressure
        LHS[i][N + 1 + i - 1] = LHS[i][N + 1 + i - 1] - 0.25 * crossSectionLength_NLS[i - 1] - 0.25 * crossSectionLength_NLS[i];
        LHS[i][N + 1 + i] = LHS[i][N + 1 + i] + 0.25 * crossSectionLength_NLS[i - 1] - 0.25 * crossSectionLength_NLS[i + 1];
        LHS[i][N + 1 + i + 1] = LHS[i][N + 1 + i + 1] + 0.25 * crossSectionLength_NLS[i] + 0.25 * crossSectionLength_NLS[i + 1];

        // Continuity, Velocity
        LHS[i + N + 1][i - 1] = LHS[i + N + 1][i - 1] - 0.25 * crossSectionLength_NLS[i - 1] - 0.25 * crossSectionLength_NLS[i];
        LHS[i + N + 1][i] = LHS[i + N + 1][i] - 0.25 * crossSectionLength_NLS[i - 1] + 0.25 * crossSectionLength_NLS[i + 1];
        LHS[i + N + 1][i + 1] = LHS[i + N + 1][i + 1] + 0.25 * crossSectionLength_NLS[i] + 0.25 * crossSectionLength_NLS[i + 1];

        // Continuity, Pressure
        LHS[i + N + 1][N + 1 + i - 1] = LHS[i + N + 1][N + 1 + i - 1] - alpha;
        LHS[i + N + 1][N + 1 + i] = LHS[i + N + 1][N + 1 + i] + 2 * alpha ;
        //LHS[i + N + 1][N + 1 + i] = LHS[i + N + 1][N + 1 + i] + 2 * alpha + gamma * dx;

        LHS[i + N + 1][N + 1 + i + 1] = LHS[i + N + 1][N + 1 + i + 1] - alpha;
      }

      // Boundary
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
      LHS[2 * N + 1][N] = -(sqrt(1 - pressure_NLS[N] / 2.0) - (velocity_NLS[N] - velocity_NLS[N]) / 4.0);

      // Solve Linear System using LAPACK

      /* LAPACK requires a 1D array
         i.e. Linearizing 2D
      */

      int counter = 0;
      for (int i = 0; i <= (2 * N + 1); i++) {
        for (int j = 0; j <= (2 * N + 1); j++) {
          A[counter] = LHS[j][i];
          counter++;
        }
      }

      /* LAPACK Function call to solve the linear system */
      dgesv_(&nlhs, &nrhs, A, &nlhs, ipiv, Res, &nlhs, &info);

      if (info != 0) {
        std::cout << "Linear Solver not converged!, Info: " << info << std::endl;
      }

      for (int i = 0; i < N; i++) {
        velocity_NLS[i] = velocity_NLS[i] + Res[i];
        pressure_NLS[i] = pressure_NLS[i] + Res[i + N + 1];
      }
    } // End of while loop

    for (int i = 0; i < chunkLength; i++) {
      pressure[i] = pressure_NLS[i];
      //pressure_n[i] = pressure_n_NLS[i];
      //p_old[i] = p_old_NLS[i];
      crossSectionLength[i] = crossSectionLength_NLS[i];
      //crossSectionLength_n[i] = crossSectionLength_n_NLS[i];
      velocity[i] = velocity_NLS[i];
      //velocity_n[i] = velocity_n_NLS[i];
    }

    for (int i = 1; i < size; i++) {
      int tagStart = 7 * i;
      int chunkLength_temp;
      int gridOffset;

      if ((N + 1) % size == 0) {
        chunkLength_temp = (N + 1) / size;
        gridOffset = i * chunkLength_temp;
      } else if (i < (N + 1) % size) {
        chunkLength_temp = (N + 1) / size + 1;
        gridOffset = i * chunkLength_temp;
      } else {
        chunkLength_temp = (N + 1) / size;
        gridOffset = ((N + 1) % size) * ((N + 1) / size + 1) + (i - ((N + 1) % size)) * (N + 1) / size;
      }

      double* vel = new double[5];
      vel[0] = 0.0123;
      vel[1] = 1e-9+0.01;

      vel[2] = 2e-9+0.01;

      vel[3] = 3e-9+0.01;
      vel[4] = 4e-9+0.01;


      MPI_Send(pressure_NLS + gridOffset, chunkLength_temp, MPI_DOUBLE, i, tagStart + 0, MPI_COMM_WORLD);
      MPI_Send(crossSectionLength_NLS + gridOffset, chunkLength_temp, MPI_DOUBLE, i, tagStart + 3, MPI_COMM_WORLD);
      MPI_Send(velocity_NLS + gridOffset, chunkLength_temp, MPI_DOUBLE, i, tagStart + 5, MPI_COMM_WORLD);

    }


    delete [] pressure_NLS;
    delete [] crossSectionLength_NLS;
    delete [] velocity_NLS;
    delete [] Res;
    delete [] LHS;
    delete [] A;
    delete [] ipiv;
  }
}
