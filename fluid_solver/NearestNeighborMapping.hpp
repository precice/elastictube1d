/*
 * NearestNeighborMapping.hpp
 *
 *  Created on: Sep 28, 2015
 *      Author: scheufks
 */

#ifndef NEARESTNEIGHBORMAPPING_HPP_
#define NEARESTNEIGHBORMAPPING_HPP_

#include <limits>
#include <vector>
#include <cmath>

class NearestNeighborQuery
{
public:
  NearestNeighborQuery(double searchPoint)
  :
  _hasFound(false),
  _searchPoint(searchPoint),
  _closestVertexID(-1)
  {}

  ~NearestNeighborQuery(){};

  void find(std::vector<double> mesh)
  {
    double shortestDistance (std::numeric_limits<double>::max());
    for(int id=0; id<mesh.size(); id++)
    {
      double dist = std::abs((double)(_searchPoint - mesh.at(id)));
      if(dist < shortestDistance){
        _hasFound = true;
        shortestDistance = dist;
        _closestVertexID = id;
      }
    }
  }

  int getClosestVertexID()
  {
    if(_hasFound) return _closestVertexID;
    return -1;
  }

private:
  bool _hasFound;
  double _searchPoint;
  int _closestVertexID;

};


class NearestNeighborMapping
{

public:
  NearestNeighborMapping();

  ~NearestNeighborMapping(){};

  void computeMapping(std::vector<int>& inputIDs,
                      std::vector<int>& outputIDs,
                      std::vector<double>& inputCoords,
                      std::vector<double>& outputCoords);

  void map(int inputN,
           int outputN,
           double* inputData,
           double* outputData);

  void setIsComputed(bool b);

private:
  bool _isComputed;
  std::vector<int> _map;

};

#endif /* NEARESTNEIGHBORMAPPING_HPP_ */
