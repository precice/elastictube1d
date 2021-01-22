#ifndef FLUID_NL_H_
#define FLUID_NL_H_

#include <vector>

#define PI 3.14159265359

int fluid_nl(std::vector<double> crossSectionLength,
             std::vector<double> crossSectionLength_n,
             std::vector<double> velocity,
             std::vector<double> velocity_n,
             std::vector<double> pressure,
             std::vector<double> pressure_n,
             double t,
             int N,
             double kappa,
             double tau);

int linsolve(int n,
             double** A,
             double* b,
             double* x);
             
void write_vtk(double t, 
				int iteration, 
				const char* filename_prefix,
				int N_slices,
				double* grid,
				std::vector<double>  velocity, 
				std::vector<double>  pressure, 
				std::vector<double>  diameter);             

#endif
