/*
 * NearestNeighborMapping.hpp
 *
 *  Created on: Sep 28, 2015
 *      Author: scheufks
 */

#ifndef NEARESTNEIGHBORMAPPING_HPP_
#define NEARESTNEIGHBORMAPPING_HPP_

#include "Mapping.hpp"
#include <limits>
#include <vector>
#include <cmath>
#include <iostream>

class NearestNeighborQuery
{
public:
  NearestNeighborQuery(double searchPoint)
  :
  _hasFound(false),
  _searchPoint(searchPoint),
  _closestVertexID(-1),
  _2ndClosestVertexID(-1)
  {}

  ~NearestNeighborQuery(){};

  void find(std::vector<double>& mesh)
  {
    double shortestDistance (std::numeric_limits<double>::max());
    double secondShortestDist = shortestDistance;
    for(int id=0; id<mesh.size(); id++)
    {
      double dist = std::abs((double)(_searchPoint - mesh.at(id)));
      if(dist < shortestDistance){
        _hasFound = true;
        shortestDistance = dist;
        _closestVertexID = id;
      }
    }
    for(int id=0; id<mesh.size(); id++)
    {
      double dist = std::abs((double)(_searchPoint - mesh.at(id)));
      if(dist < secondShortestDist){
        if(id != _closestVertexID)
        {
          secondShortestDist = dist;
          _2ndClosestVertexID = id;
        }
      }
    }
  }

  int getClosestVertexID()
  {
    if(_hasFound) return _closestVertexID;
    return -1;
  }

  int get2ndClosestVertexID()
  {
    if(_hasFound) return _2ndClosestVertexID;
    return -1;
  }

private:
  bool _hasFound;
  double _searchPoint;
  int _closestVertexID;
  int _2ndClosestVertexID;

};


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
