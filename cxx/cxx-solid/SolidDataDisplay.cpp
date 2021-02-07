#include "SolidSolver.h"

#include <iostream>

void SolidDataDisplay(double* data, int length)
{
  std::cout << "\n";
  for (int i = 0; i < length; i++) {
    std::cout << data[i] << " ";
  }
  std::cout << "\n";
}
