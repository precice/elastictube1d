/*
 * NearestNeighborQuery.hpp
 *
 *  Created on: Oct 15, 2015
 *      Author: scheufks
 */

#ifndef NEARESTNEIGHBORQUERY_HPP_
#define NEARESTNEIGHBORQUERY_HPP_

#include <limits>
#include <vector>
#include <cmath>
#include <iostream>

class NearestNeighborQuery
{
public:
  NearestNeighborQuery(double searchPoint);

  ~NearestNeighborQuery(){};

  void find(std::vector<double>& mesh);

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



#endif /* NEARESTNEIGHBORQUERY_HPP_ */
