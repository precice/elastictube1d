#include "FluidSolver.h"

void fluidDataDisplay(double *data, int length){
	std::cout << "\nFluid Solver Data: ";
	std::cout << "\n";
	for (int i = 0; i < length; i++)
		std::cout << data[i] << " ";
	std::cout << "\n";
}
