#pragma once

void structureInit(int chunkLength, double* data);

void structureComputeSolution(int rank, int size, int chunkLength, double* pressure, double* crossSectionLength);

void structureDataDisplay(double* data, int length);
