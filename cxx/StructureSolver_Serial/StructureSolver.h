#pragma once
#include <vector>

void structureInit(int chunkLength, std::vector<double> data);

void structureComputeSolution(int rank, int size, int chunkLength, std::vector<double> pressure, std::vector<double> crossSectionLength);

void structureDataDisplay(std::vector<double> data, int length);
