#include "StructureSolver.h"

void structureComputeSolution(int rank, int size, int chunkLength, double *pressure, double *crossSectionLength){
  /*
   * Update displacement of membrane based on pressure data from the fluid solver
   */

  for (int i = 0; i < chunkLength; i++) {
    crossSectionLength[i] = 4.0 / ((2.0 - pressure[i]) * (2.0 - pressure[i]));;
  }
}
