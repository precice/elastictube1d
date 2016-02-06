#ifndef FLUIDSOLVER_H_
#define FLUIDSOLVER_H_

#define PI 3.14159265359

#include <iostream>
#include <string>
#include <cmath>
#include <stdlib.h>
#include <mpi.h>

void fluidInit(
		int rank,
		int chunkLength,
		double *data
		);

void fluidComputeSolution(
		int rank,
		int size,
		int domainSize,
		int chunkLength,
		double kappa,
		double tau,
		double gamma,
		double t,
		double *pressure,
		double *pressure_n,
		double *pressure_old,
		double *crossSectionLength,
		double *crossSectionLength_n,
		double *velocity,
		double *velocity_n
		);

void fluidDataDisplay(
		double *data,
		int counterLength
		);

#endif
