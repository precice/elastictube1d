/*
 * TestMapping.cpp
 *
 *  Created on: Oct 20, 2015
 *      Author: scheufks
 */

#include <iostream>
#include <stdlib.h>
#include "precice/SolverInterface.hpp"
#include "../mapping/NearestNeighborMapping.hpp"
#include "../mapping/LinearInterpolationMapping.hpp"
//#include "mpi.h"

using std::cout;
using std::endl;

void printData (const std::vector<double>& data)
{
  cout << "Received data = " << data[0];
  for (size_t i=1; i < data.size(); i++){
    cout << ", " << data[i];
  }
  cout << endl;
}


void testNearestNeighborMapping()
{

  // mappings
    NearestNeighborMapping upMapping, downMapping;
    int N = 10, N_SM = 5;

    double *displ_copy_coarse, *displ_coarse, *tmp;
    displ_coarse = new double[N_SM];
    tmp = new double[N_SM];
    displ_copy_coarse = new double[N];

    for (int i = 0; i < N; i++)
      displ_copy_coarse[i] = i;
    for (int i = 0; i < N_SM; i++)
          displ_coarse[i] = i;
    std::cout<<"init: ";
    for(int i = 0; i <N; i++)
      std::cout<<displ_copy_coarse[i]<<", ";
    std::cout<<"\n";

    downMapping.map(N, N_SM, displ_copy_coarse, tmp);

    std::cout<<"\n\nfine --> coarse: ";
    for(int i = 0; i <N_SM; i++)
      std::cout<<tmp[i]<<", ";

    upMapping.map(N_SM, N, tmp, displ_copy_coarse);

    std::cout<<"\n\ncoarse --> fine: ";
    for(int i = 0; i <N; i++)
      std::cout<<displ_copy_coarse[i]<<", ";
}


void testLinerInterpolationMapping()
{

  // mappings
    LinearInterpolationMapping upMapping, downMapping;
    int N = 10, N_SM = 7;

    double *displ_copy_coarse, *displ_coarse, *tmp, *tmp2;
    displ_coarse = new double[N_SM];
    tmp = new double[N_SM];
    tmp2 = new double[N];
    displ_copy_coarse = new double[N];

    for (int i = 0; i < N; i++)
      displ_copy_coarse[i] = i;
    for (int i = 0; i < N_SM; i++)
          displ_coarse[i] = i;
    std::cout<<"init: ";
    for(int i = 0; i <N; i++)
      std::cout<<displ_copy_coarse[i]<<", ";
    std::cout<<"\n";

    downMapping.map(N, N_SM, displ_copy_coarse, tmp);

    std::cout<<"\n\nfine --> coarse: ";
    for(int i = 0; i <N_SM; i++)
      std::cout<<tmp[i]<<", ";
    std::cout<<"\n";

    upMapping.map(N_SM, N, tmp, tmp2);

    std::cout<<"\n\ncoarse --> fine: ";
    for(int i = 0; i <N; i++)
      std::cout<<tmp2[i]<<", ";
    std::cout<<"\n";
}


int main (int argc, char **argv)
{
  std::coutcout << "Test Mapping ... " << std::endl;

  std::cout<<"\n .. test nearest neighbor mapping"<<std::endl;
  testNearestNeighborMapping();

  std::cout<<"\n .. test linear interpolation mapping"<<std::endl;
  testLinerInterpolationMapping();
}


