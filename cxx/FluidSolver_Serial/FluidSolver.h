#pragma once

const double PI = 3.14159265359;

void fluidInit(
    int rank,
    int chunkLength,
    double* data);

void fluidComputeSolution(
    int rank,
    int size,
    int domainSize,
    int chunkLength,
    double kappa,
    double tau,
    double gamma,
    double t,
    double* pressure,
    double* pressure_n,
    double* pressure_old,
    double* crossSectionLength,
    double* crossSectionLength_n,
    double* velocity,
    double* velocity_n);

void fluidDataDisplay(
    double* data,
    int counterLength);
