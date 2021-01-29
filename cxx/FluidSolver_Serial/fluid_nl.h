#ifndef FLUID_NL_H_
#define FLUID_NL_H_
#include <vector>

#define PI 3.14159265359

int fluid_nl(double* crossSectionLength,
             double* crossSectionLength_n,
             double* velocity,
             double* velocity_n,
             double* pressure,
             double* pressure_n,
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
				double* velocity, 
				double* pressure, 
				double* diameter);             

#endif
