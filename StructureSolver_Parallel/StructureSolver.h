#ifndef STRUCTURESOLVER_H_
#define STRUCTURESOLVER_H_

#include <iostream>
#include <string>
#include <stdlib.h>
#include <mpi.h>

void structureInit(int chunkLength, double *data);

void structureComputeSolution(int rank, int size, int chunkLength, double *pressure, double *crossSectionLength);

void structureDataDisplay(double *data, int length);

#endif
