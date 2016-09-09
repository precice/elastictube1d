#include "AD.hpp"

#include <cmath>
#include <assert.h>

namespace AD {

  // dualReal methods ------------------------------------------
  
  dualReal operator+(const dualReal& lhs, const dualReal& rhs){
   return  dualReal(lhs.u + rhs.u, lhs.v + rhs.v);
  }
  
  dualReal operator-(const dualReal& lhs, const dualReal& rhs){
   return  dualReal(lhs.u - rhs.u, lhs.v - rhs.v);
  }
  
  dualReal operator-(const dualReal& rhs){
   return  dualReal(-rhs.u, -rhs.v);
  }
  
  dualReal operator*(const dualReal& lhs, const dualReal& rhs){
   return  dualReal(lhs.u * rhs.u, lhs.v * rhs.u + lhs.u * rhs.v);
  }
  
  dualReal operator/(const dualReal& lhs, const dualReal& rhs){
   if(lhs.u == 0 || rhs.u == 0) 
     return dualReal(0.); 
   else
     return  dualReal(lhs.u / rhs.u, (lhs.v * rhs.u - lhs.u * rhs.v)/(rhs.u * rhs.u));
  }
  
  dualReal sin(const dualReal& x){
   return dualReal(std::sin(x.u), x.v*std::cos(x.u)); 
  }
  
  dualReal cos(const dualReal& x){
   return dualReal(std::cos(x.u), -x.v*std::sin(x.u)); 
  }
  
  dualReal exp(const dualReal& x){
   return dualReal(std::exp(x.u), x.v*std::exp(x.u)); 
  }
  
  dualReal log(const dualReal& x){
   if(x.u == 0)
     return dualReal(std::log(x.u)); 
   else
     return dualReal(std::log(x.u), x.v/x.u); 
  }
  
  dualReal pow(const dualReal& x, int e){
   return dualReal(std::pow(x.u, e), e*std::pow(x.u,e-1)*x.v); 
  }
  
  dualReal abs(const dualReal& x){
   return dualReal(std::abs(x.u), x.v*(x.u/std::abs(x.u))); 
  }
  
  dualReal sqrt(const dualReal& x){
   return dualReal(std::sqrt(x.u), .5*1./std::sqrt(x.u)*x.v); 
  }
  
  // dReal_grad methods ------------------------------------------
  
  dReal_grad operator+(const dReal_grad& lhs, const dReal_grad& rhs){
    std::size_t n = std::max(lhs.v.size(), rhs.v.size());
    std::size_t m = std::min(lhs.v.size(), rhs.v.size());
    if (n == 1) 
      return  dReal_grad(lhs.u + rhs.u, lhs.v[0] + rhs.v[0]);
    std::vector<double> _v(n);
    for(std::size_t i = 0; i < n; i++){
      if(i < m)
        _v[i] = lhs.v[i] + rhs.v[i];
      else if(n == lhs.v.size())
	_v[i] = lhs.v[i]; 
      else
	_v[i] = rhs.v[i];
    }
    return  dReal_grad(lhs.u + rhs.u, _v);
  }
  
  dReal_grad operator-(const dReal_grad& lhs, const dReal_grad& rhs){
   std::size_t n = std::max(lhs.v.size(), rhs.v.size());
    std::size_t m = std::min(lhs.v.size(), rhs.v.size());
    if (n == 1) 
      return  dReal_grad(lhs.u - rhs.u, lhs.v[0] - rhs.v[0]);
    std::vector<double> _v(n);
    for(std::size_t i = 0; i < n; i++){
      if(i < m)
        _v[i] = lhs.v[i] - rhs.v[i];
      else if(n == lhs.v.size())
	_v[i] = lhs.v[i]; 
      else
	_v[i] = -rhs.v[i];
    }
    return  dReal_grad(lhs.u - rhs.u, _v);
  }
  
  void operator-(dReal_grad& rhs){
    for(std::size_t i = 0; i < rhs.v.size(); i++)
      rhs.v[i] = -rhs.v[i];
    rhs.u = -rhs.u;
  }
  
  dReal_grad operator*(const dReal_grad& lhs, const dReal_grad& rhs){
    std::size_t n = std::max(lhs.v.size(), rhs.v.size());
    std::size_t m = std::min(lhs.v.size(), rhs.v.size());
    if (n == 1) 
      return  dReal_grad(lhs.u * rhs.u, lhs.v[0] * rhs.u + lhs.u * rhs.v[0]);
    std::vector<double> _v(n);
    for(std::size_t i = 0; i < n; i++){
      if(i < m)
        _v[i] = lhs.v[i] * rhs.u + lhs.u * rhs.v[i];
      else if(n == lhs.v.size())
	_v[i] = lhs.v[i] * rhs.u; 
      else
	_v[i] = rhs.v[i] * lhs.u;
    }
    return  dReal_grad(lhs.u * rhs.u, _v);
  }
  
  dReal_grad operator/(const dReal_grad& lhs, const dReal_grad& rhs){
    std::size_t n = std::max(lhs.v.size(), rhs.v.size());
    std::size_t m = std::min(lhs.v.size(), rhs.v.size());
    if(lhs.u == 0 || rhs.u == 0){
        return dReal_grad(0.); 
    }else if (n == 1){
        return  dReal_grad(lhs.u / rhs.u, (lhs.v[0] * rhs.u - lhs.u * rhs.v[0])/(rhs.u * rhs.u));
    }else{
      std::vector<double> _v(n);
      for(std::size_t i = 0; i < n; i++){
        if(i < m)
          _v[i] = (lhs.v[i] * rhs.u - lhs.u * rhs.v[i])/(rhs.u * rhs.u);
        else if(n == lhs.v.size())
	  _v[i] = (lhs.v[i] * rhs.u)/(rhs.u * rhs.u); 
        else
	  _v[i] = (-lhs.u * rhs.v[i])/(rhs.u * rhs.u);
      }
      return  dReal_grad(lhs.u / rhs.u, _v);
    }
  }
  
  void sin(dReal_grad& x){
    for(std::size_t i = 0; i < x.v.size(); i++)
      x.v[i] = x.v[i]*std::cos(x.u);
    x.u = std::sin(x.u);
  }
  
  void cos(dReal_grad& x){
    for(std::size_t i = 0; i < x.v.size(); i++)
      x.v[i] = -x.v[i]*std::sin(x.u);
    x.u = std::cos(x.u);
  }
  
  void exp(dReal_grad& x){
    for(std::size_t i = 0; i < x.v.size(); i++)
      x.v[i] = x.v[i]*std::exp(x.u);
    x.u = std::exp(x.u);
  }
  
  void log(dReal_grad& x){
    if(x.u != 0)
      for(std::size_t i = 0; i < x.v.size(); i++)
        x.v[i] = x.v[i]/x.u;
    x.u = std::log(x.u);
  }
  
  void pow(dReal_grad& x, int e){
    for(std::size_t i = 0; i < x.v.size(); i++)
      x.v[i] = e*std::pow(x.u,e-1)*x.v[i];
    x.u = std::pow(x.u, e);
  }
  
  void abs(dReal_grad& x){
    for(std::size_t i = 0; i < x.v.size(); i++)
      x.v[i] = x.v[i]*(x.u/std::abs(x.u));
    x.u = std::abs(x.u);
  }
  
  void sqrt(dReal_grad& x){
    for(std::size_t i = 0; i < x.v.size(); i++)
      x.v[i] = .5*1./std::sqrt(x.u)*x.v[i];
    x.u = std::sqrt(x.u);
  } 
  
}
