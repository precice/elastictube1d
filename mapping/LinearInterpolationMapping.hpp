/*
 * LinearInterpolationMapping.hpp
 *
 *  Created on: Oct 6, 2015
 *      Author: scheufks
 */

#ifndef LINEARINTERPOLATIONMAPPING_HPP_
#define LINEARINTERPOLATIONMAPPING_HPP_

#include "Mapping.hpp"
#include "NearestNeighborQuery.hpp"
#include <vector>
#include <cmath>
#include <iostream>
#include "Eigen/Dense"


class LinearInterpolationMapping : public Mapping
{

public:
  LinearInterpolationMapping();

  virtual ~LinearInterpolationMapping(){};

  virtual void computeMapping(int inputN,
                              int outputN,
                              std::vector<double>& inputCoords,
                              std::vector<double>& outputCoords);

  virtual void map(int inputN,
           int outputN,
           double* inputData,
           double* outputData);


private:
  Eigen::MatrixXd _map;
};


#endif /* LINEARINTERPOLATIONMAPPING_HPP_ */
