#include "FluidSolver.h"

#include <iostream>
#include <vector>

void fluidDataDisplay(double* data, int length)
{
  std::cout << "\n";
  for (int i = 0; i < length; i++)
    std::cout << data[i] << " ";
  std::cout << "\n";
}
