#include "AD.hpp"

#include <cmath>
#include <assert.h>
#include <boost/concept_check.hpp>

namespace AD {

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
   if(lhs.u == 0) return dualReal(0.); 
   if(rhs.u == 0)
     return dualReal(lhs.u / rhs.u, 0);
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
  
}
