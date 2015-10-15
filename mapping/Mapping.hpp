/*
 * Mapping.hpp
 *
 *  Created on: Oct 6, 2015
 *      Author: scheufks
 */

#ifndef MAPPING_HPP_
#define MAPPING_HPP_

#include <limits>
#include <vector>
#include <cmath>
#include <iostream>


class Mapping
{

public:
  Mapping()
    :_isComputed(false)
 {};

  virtual ~Mapping(){};

  virtual void computeMapping(int inputN,
                                int outputN,
                                std::vector<double>& inputCoords,
                                std::vector<double>& outputCoords) = 0;

  virtual void map(int inputN,
           int outputN,
           double* inputData,
           double* outputData) = 0;

  virtual void setIsComputed(bool b)
  {
    _isComputed = false;
  }

protected:
  bool _isComputed;

};




#endif /* MAPPING_HPP_ */
