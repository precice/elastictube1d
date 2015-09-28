/*
 * NearestNeighborMapping.cpp
 *
 *  Created on: Sep 28, 2015
 *      Author: scheufks
 */


#include "NearestNeighborMapping.hpp"

NearestNeighborMapping::NearestNeighborMapping()
  :
_isComputed(false),
_map()
{}

void NearestNeighborMapping::computeMapping
(
    std::vector<int>& inputIDs,
    std::vector<int>& outputIDs,
    std::vector<double>& inputCoords,
    std::vector<double>& outputCoords)
{
  _map.resize(outputIDs.size());
  for (int i=0; i < outputIDs.size(); i++ ){
      double outCoord = outputCoords.at(i);
      NearestNeighborQuery query(outCoord);
      query.find(inputCoords);

      _map[i] = query.getClosestVertexID();
    }
  _isComputed = true;
}

void NearestNeighborMapping::map
(
    int inputN,
    int outputN,
    double* inputData,
    double* outputData)
{
  if(_isComputed){

    for(int i=0; i<_map.size(); i++){
      outputData[i] = inputData[_map.at(i)];
    }

  }else{
    std::vector<int> inputIDs(inputN);
    std::vector<int> outputIDs(outputN);
    for(int i=0; i<inputIDs.size(); i++)
      inputIDs.at(i) = i;
    for(int i=0; i<outputIDs.size(); i++)
        outputIDs.at(i) = i;

    std::vector<double> _inputCoords(inputN);
    std::vector<double> _outputCoords(outputN);
    for(int i=0; i<_inputCoords.size(); i++)
      _inputCoords.at(i) = 1./(double)i;
    for(int i=0; i<_outputCoords.size(); i++)
        _outputCoords.at(i) = 1./(double)i;

    computeMapping(inputIDs, outputIDs, _inputCoords, _outputCoords);
  }
}

void NearestNeighborMapping::setIsComputed
(
    bool b)
{
  _isComputed = false;
}




