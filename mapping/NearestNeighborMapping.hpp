/*
 * NearestNeighborMapping.hpp
 *
 *  Created on: Sep 28, 2015
 *      Author: scheufks
 */

#ifndef NEARESTNEIGHBORMAPPING_HPP_
#define NEARESTNEIGHBORMAPPING_HPP_

#include "Mapping.hpp"
#include "NearestNeighborQuery.hpp"
#include <limits>
#include <vector>
#include <cmath>
#include <iostream>



class NearestNeighborMapping : public Mapping
{

public:
  NearestNeighborMapping();

  virtual ~NearestNeighborMapping(){};

  virtual void computeMapping(int inputN,
                                int outputN,
                                std::vector<double>& inputCoords,
                                std::vector<double>& outputCoords);

  virtual void map(int inputN,
           int outputN,
           double* inputData,
           double* outputData);


private:
  std::vector<int> _map;

};

#endif /* NEARESTNEIGHBORMAPPING_HPP_ */
