/*
 * NearestNeighborQuery.cpp
 *
 *  Created on: Oct 15, 2015
 *      Author: scheufks
 */

#include "NearestNeighborQuery.hpp"


NearestNeighborQuery::NearestNeighborQuery
(
    double searchPoint)
:
    _hasFound(false),
    _searchPoint(searchPoint),
    _closestVertexID(-1),
    _2ndClosestVertexID(-1)
{}


void NearestNeighborQuery::find
(
    std::vector<double>& mesh)
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
