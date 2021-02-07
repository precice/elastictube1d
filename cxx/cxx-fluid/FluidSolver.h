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
    double* crossSectionLength,
    double* velocity);

void fluidDataDisplay(
    double* data,
    int counterLength);
