#include "StructureSolver.h"

void structureDataDisplay(double* data, int length)
{
  std::cout << "\n";
  for (int i = 0; i < length; i++) {
    std::cout << data[i] << " ";
  }
  std::cout << "\n";
}
