/*
 * LinearInterpolationMapping.cpp
 *
 *  Created on: Oct 6, 2015
 *      Author: scheufks
 */

#include "LinearInterpolationMapping.hpp"
#include "NearestNeighborQuery.hpp"

#include <iostream>

LinearInterpolationMapping::LinearInterpolationMapping()
  :
  _isComputed(false),
  _map()
{}

void LinearInterpolationMapping::computeMapping
(
    int inputN,
    int outputN,
    std::vector<double>& inputCoords,
    std::vector<double>& outputCoords)
{
  _map.resize(outputN, inputN);
  for (int i=0; i < outputN; i++ ){
    double outCoord = outputCoords.at(i);
    int outID = outCoord * outputN;
    NearestNeighborQuery query(outCoord);
    query.find(inputCoords);
    int closest = query.getClosestVertexID();
    int second = query.get2ndClosestVertexID();
    double lambda = (inputCoords.at(closest) - outputCoords.at(outID))/(inputCoords.at(closest) - inputCoords.at(second));
    _map(outID, closest) = (1.-lambda);
    _map(outID, second) = lambda;
    }
    _isComputed = true;
}

void LinearInterpolationMapping::map
(
    int inputN,
    int outputN,
    double* inputData,
    double* outputData)
{
  if(not _isComputed){
    std::vector<double> _inputCoords(inputN);
    std::vector<double> _outputCoords(outputN);
    for(int i=1; i<_inputCoords.size(); i++)
      _inputCoords.at(i) = (double)i/(double)inputN;
    for(int i=1; i<_outputCoords.size(); i++)
      _outputCoords.at(i) = (double)i/(double)outputN;

    computeMapping(inputN, outputN, _inputCoords, _outputCoords);
  }

  Eigen::VectorXd inData(inputN);
  Eigen::VectorXd outData(outputN);
  for (int i = 0; i < inputN; i++)
    inData(i) = inputData[i];

 // std::cout<<"Map:\n"<<_map<<std::endl;
  assert(outputN == _map.rows() && inputN == _map.cols());
  outData = _map*inData;

  for(int i=0; i<outputN; i++)
    outputData[i] = outData(i);

}




