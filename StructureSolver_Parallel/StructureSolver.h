#ifndef STRUCTURESOLVER_H_
#define STRUCTURESOLVER_H_

#include <iostream>
#include <mpi.h>
#include <stdlib.h>
#include <string>

void structureInit(int chunkLength, double* data);

void structureComputeSolution(int rank, int size, int chunkLength, double* pressure, double* crossSectionLength);

void structureDataDisplay(double* data, int length);

#endif
