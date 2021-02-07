#pragma once

void SolidInit(int chunkLength, double* data);

void SolidComputeSolution(int rank, int size, int chunkLength, double* pressure, double* crossSectionLength);

void SolidDataDisplay(double* data, int length);
