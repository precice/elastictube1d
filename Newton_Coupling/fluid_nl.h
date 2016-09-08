
#ifndef FLUID_NL_H_
#define FLUID_NL_H_

#include "../utils/AD.hpp"

#define PI 3.14159265359
using namespace AD;

int fluid_nl(dualReal* crossSectionLength,
             dualReal* crossSectionLength_n,
             dualReal* velocity,
             dualReal* velocity_n,
             dualReal* pressure,
             dualReal* pressure_n,
             dualReal* pressure_old,
             double scaled_t,
             int N,
             dualReal kappa,
             dualReal tau,
             dualReal gamma);

int linsolve(int n,
             double** A,
             double* b,
             double* x);

#endif
