#ifndef AD_H_
#define AD_H_

namespace AD {

  struct dualReal{
    
    double u;
    double v;
    dualReal() : u(0), v(0) {}
    dualReal(double _u) : u(_u), v(0) {}
    dualReal(double _u, double _v) : u(_u), v(_v) {}
  };
  
  dualReal operator+(const dualReal& lhs, const dualReal& rhs);
  
  dualReal operator-(const dualReal& lhs, const dualReal& rhs);
  
  dualReal operator-(const dualReal& rhs);
  
  dualReal operator*(const dualReal& lhs, const dualReal& rhs);
  
  dualReal operator/(const dualReal& lhs, const dualReal& rhs);
  
  dualReal sin(const dualReal& x);
  
  dualReal cos(const dualReal& x);
  
  dualReal exp(const dualReal& x);
  
  dualReal log(const dualReal& x);
  
  dualReal pow(const dualReal& x, int e);
  
  dualReal abs(const dualReal& x);
  
  dualReal sqrt(const dualReal& x);
}
#endif
